// Packed presence-mask grouping for trait-wise NA pruning.
//
// A column of the input matrix is interpreted as a binary presence pattern
// (TRUE/1 means retained, FALSE/0 means omitted).  Patterns are represented
// as little-endian 64-bit words (row 1 is bit 0 of word 0).  Hashes are used
// only to find candidate groups; every candidate is compared word-for-word,
// so a hash collision can never merge two different masks.

#include <Rcpp.h>

#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using Word = std::uint64_t;

// A deterministic 64-bit mixer.  This is deliberately not used as the
// grouping identity: equal hashes are always resolved by exact word equality.
Word mix_word(Word x) {
  x ^= x >> 30;
  x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27;
  x *= UINT64_C(0x94d049bb133111eb);
  x ^= x >> 31;
  return x;
}

Word hash_words(const std::vector<Word>& words) {
  Word h = UINT64_C(0x243f6a8885a308d3);
  for (std::size_t i = 0; i < words.size(); ++i) {
    const Word salt = UINT64_C(0x9e3779b97f4a7c15) *
      static_cast<Word>(i + 1);
    h = mix_word(h ^ mix_word(words[i] + salt));
  }
  return mix_word(h ^ static_cast<Word>(words.size()));
}

struct MaskGroup {
  std::vector<Word> words;
  std::vector<int> columns;  // 1-based input columns
  std::vector<int> keep;     // 1-based retained rows/tips
};

bool read_presence(SEXP x, const R_xlen_t offset, int type) {
  if (type == LGLSXP) {
    const int value = LOGICAL(x)[offset];
    if (value == NA_LOGICAL) {
      Rcpp::stop("present contains NA; presence masks must be complete.");
    }
    return value != 0;
  }

  if (type == INTSXP) {
    const int value = INTEGER(x)[offset];
    if (value == NA_INTEGER) {
      Rcpp::stop("present contains NA; presence masks must be complete.");
    }
    if (value != 0 && value != 1) {
      Rcpp::stop("numeric presence masks must contain only 0 and 1.");
    }
    return value != 0;
  }

  const double value = REAL(x)[offset];
  if (R_IsNA(value) || R_IsNaN(value) || !R_FINITE(value)) {
    Rcpp::stop("present contains NA/non-finite values; presence masks must be complete.");
  }
  if (value != 0.0 && value != 1.0) {
    Rcpp::stop("numeric presence masks must contain only 0 and 1.");
  }
  return value != 0.0;
}

std::string printable_key(const std::vector<Word>& words,
                          const R_xlen_t n_rows) {
  // Include n_rows so a key remains unambiguous when printed or persisted
  // across matrices with different dimensions.  Fixed-width hexadecimal
  // words are platform independent and preserve leading zero words.
  std::ostringstream stream;
  stream << "n=" << n_rows << ";" << std::hex << std::setfill('0');
  for (std::size_t i = 0; i < words.size(); ++i) {
    if (i != 0) stream << ':';
    stream << std::setw(16) << words[i];
  }
  return stream.str();
}

}  // namespace

// Group columns of a logical (or 0/1 numeric) presence matrix.
//
// The grouping ID is assigned in first-occurrence order, starting at one.
// `columns` and `keep` are lists in that same group order and contain 1-based
// indices.  `count` is the number of columns in each group; `n_keep` is the
// number of retained rows in each group.  `key` is a stable printable debug
// label and is not used to decide equality.
//
// [[Rcpp::export]]
Rcpp::List group_na_masks_cpp(SEXP present) {
  const int type = TYPEOF(present);
  if (type != LGLSXP && type != INTSXP && type != REALSXP) {
    Rcpp::stop("present must be a logical or numeric matrix.");
  }

  SEXP dim = Rf_getAttrib(present, R_DimSymbol);
  if (Rf_isNull(dim) || Rf_length(dim) != 2) {
    Rcpp::stop("present must be a matrix with two dimensions.");
  }
  const int* dims = INTEGER(dim);
  const R_xlen_t n_rows = static_cast<R_xlen_t>(dims[0]);
  const R_xlen_t n_cols = static_cast<R_xlen_t>(dims[1]);
  if (n_rows < 0 || n_cols < 0) {
    Rcpp::stop("present matrix dimensions must be non-negative.");
  }

  const std::size_t n_word = static_cast<std::size_t>(
    (n_rows + static_cast<R_xlen_t>(63)) / static_cast<R_xlen_t>(64)
  );

  std::vector<MaskGroup> groups;
  groups.reserve(static_cast<std::size_t>(n_cols));
  std::unordered_map<Word, std::vector<std::size_t> > buckets;
  buckets.reserve(static_cast<std::size_t>(n_cols) * 2 + 1);

  Rcpp::IntegerVector group_id(n_cols);
  for (R_xlen_t column = 0; column < n_cols; ++column) {
    std::vector<Word> words(n_word, UINT64_C(0));
    for (R_xlen_t row = 0; row < n_rows; ++row) {
      const R_xlen_t offset = row + n_rows * column;  // R is column-major
      if (read_presence(present, offset, type)) {
        const std::size_t word = static_cast<std::size_t>(row >> 6);
        const unsigned int bit = static_cast<unsigned int>(row & 63);
        words[word] |= (UINT64_C(1) << bit);
      }
    }

    const Word hash = hash_words(words);
    std::size_t group = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t>& candidates = buckets[hash];
    for (std::size_t i = 0; i < candidates.size(); ++i) {
      const std::size_t candidate = candidates[i];
      if (groups[candidate].words == words) {
        group = candidate;
        break;
      }
    }
    if (group == std::numeric_limits<std::size_t>::max()) {
      group = groups.size();
      MaskGroup fresh;
      fresh.words.swap(words);
      groups.push_back(fresh);
      candidates.push_back(group);

      // Keep rows in increasing order.  This is done once per new group,
      // yielding O(n p) mask scanning plus O(n * groups) retained-index data.
      for (R_xlen_t row = 0; row < n_rows; ++row) {
        const R_xlen_t offset = row + n_rows * column;
        if (read_presence(present, offset, type)) {
          groups[group].keep.push_back(static_cast<int>(row + 1));
        }
      }
    }

    groups[group].columns.push_back(static_cast<int>(column + 1));
    group_id[column] = static_cast<int>(group + 1);
  }

  const R_xlen_t n_group = static_cast<R_xlen_t>(groups.size());
  Rcpp::List columns(n_group);
  Rcpp::List keep(n_group);
  Rcpp::IntegerVector count(n_group);
  Rcpp::IntegerVector n_keep(n_group);
  Rcpp::CharacterVector key(n_group);
  for (R_xlen_t i = 0; i < n_group; ++i) {
    const MaskGroup& group = groups[static_cast<std::size_t>(i)];
    columns[i] = Rcpp::wrap(group.columns);
    keep[i] = Rcpp::wrap(group.keep);
    count[i] = static_cast<int>(group.columns.size());
    n_keep[i] = static_cast<int>(group.keep.size());
    key[i] = printable_key(group.words, n_rows);
  }

  return Rcpp::List::create(
    Rcpp::_["group_id"] = group_id,
    Rcpp::_["columns"] = columns,
    Rcpp::_["keep"] = keep,
    Rcpp::_["count"] = count,
    Rcpp::_["n_keep"] = n_keep,
    Rcpp::_["key"] = key,
    Rcpp::_["n_group"] = static_cast<int>(n_group),
    Rcpp::_["n_word"] = static_cast<int>(n_word)
  );
}

