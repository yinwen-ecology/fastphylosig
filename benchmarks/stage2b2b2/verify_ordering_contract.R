#!/usr/bin/env Rscript

# Candidate 2 bounded ordering-contract smoke.  This script is deliberately
# not a benchmark and writes only to the caller-supplied output directory.

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1L && nzchar(args[[1L]])) {
  normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
if (length(args) < 2L || !nzchar(args[[2L]])) {
  stop("usage: verify_ordering_contract.R <repo> <ascii-output-dir>",
       call. = FALSE)
}
out_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
locale_set <- tryCatch(
  Sys.setlocale("LC_COLLATE", "C"),
  warning = function(w) conditionMessage(w),
  error = function(e) conditionMessage(e)
)
if (!identical(locale_set, "C")) {
  stop(sprintf("LC_COLLATE could not be set to C (got %s).", locale_set),
       call. = FALSE)
}
options(stringsAsFactors = FALSE)

installed_lib <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_LIB", "")
if (nzchar(installed_lib) &&
    dir.exists(file.path(installed_lib, "fastphylosig"))) {
  suppressPackageStartupMessages(
    library("fastphylosig", lib.loc = installed_lib, character.only = TRUE)
  )
  package_load <- paste0("installed library: ",
                         normalizePath(installed_lib, winslash = "/"))
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  package_load <- "pkgload::load_all(source)"
} else {
  stop("pkgload is required unless FASTPHYLOSIG_STAGE2B2B2_LIB is set.",
       call. = FALSE)
}
suppressPackageStartupMessages(library(ape))
ns <- asNamespace("fastphylosig")
get_internal <- function(name) get(name, envir = ns, inherits = FALSE)

# Four unique, non-empty labels are partitioned into two non-root clades.
# Under C collation both collapsed keys are exactly a\rb\rc.
tip_labels <- c("a", "b\rc", "a\rb", "c")
clade_a <- 5L
clade_b <- 6L
root <- 7L
tree <- structure(
  list(
    edge = matrix(
      c(root, clade_a,
        clade_a, 1L,
        clade_a, 2L,
        root, clade_b,
        clade_b, 3L,
        clade_b, 4L),
      ncol = 2L,
      byrow = TRUE
    ),
    edge.length = as.numeric(seq_len(6L)),
    tip.label = tip_labels,
    Nnode = 3L
  ),
  class = "phylo"
)

input_before <- serialize(tree, connection = NULL, version = 2L)
inspect <- get_internal(".inspect_tree_core")(
  tree, signal = c("K", "lambda", "D", "Delta")
)
children <- get_internal(".edge_children_index")(tree$edge)
keys <- get_internal(".descendant_keys_iterative")(
  c(clade_a, clade_b), children, as.character(tree$tip.label),
  length(tree$tip.label)
)
expected_key <- "a\rb\rc"
normalized <- get_internal(".safe_canonicalize_core")(tree)
metadata <- attr(normalized, "fastphylosig_canonicalization", exact = TRUE)
mapping <- if (is.list(metadata)) metadata$mapping else NULL
old_order <- if (is.list(mapping)) as.integer(mapping$old_internal_order) else
  integer()
new_order <- if (is.list(mapping)) as.integer(mapping$new_internal_order) else
  integer()

same_integer <- function(x, y) identical(as.integer(x), as.integer(y))
escape_text <- function(x) {
  if (!length(x) || is.na(x[[1L]])) return("NA")
  value <- as.character(x[[1L]])
  value <- gsub("\r", "\\r", value, fixed = TRUE)
  value <- gsub("\n", "\\n", value, fixed = TRUE)
  value
}
collapse_integer <- function(x) paste(as.integer(x), collapse = ",")

tree_valid <- isTRUE(inspect$tree_summary$topology_valid) &&
  isTRUE(inspect$tree_summary$structural_rooted) &&
  isTRUE(inspect$tree_summary$fully_dichotomous) &&
  isTRUE(inspect$tree_summary$branch_lengths_positive) &&
  !anyDuplicated(tree$tip.label)
keys_equal <- length(keys) == 2L && identical(keys[[1L]], keys[[2L]])
keys_expected <- length(keys) == 2L &&
  identical(as.character(keys[[1L]]), expected_key) &&
  identical(as.character(keys[[2L]]), expected_key)
metadata_safe <- is.list(metadata) && isTRUE(metadata$safe)
expected_old_order <- c(root, clade_a, clade_b)
expected_new_order <- c(5L, 6L, 7L)
mapping_order <- same_integer(old_order, expected_old_order) &&
  same_integer(new_order, expected_new_order)
numeric_tie_break <- keys_equal && mapping_order &&
  identical(old_order[2:3], c(clade_a, clade_b))
input_immutable <- identical(
  serialize(tree, connection = NULL, version = 2L), input_before
)

source_commit <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_SOURCE_COMMIT", "unknown")
if (!nzchar(source_commit)) source_commit <- "unknown"
package_version <- as.character(utils::packageVersion("fastphylosig"))

checks <- data.frame(
  check = c(
    "LC_COLLATE",
    "tree_is_rooted_binary_and_positive",
    "clade_A_labels",
    "clade_B_labels",
    "clade_A_key",
    "clade_B_key",
    "collapsed_keys_identical",
    "collapsed_key_matches_expected",
    "safe_canonicalization",
    "old_internal_order",
    "new_internal_order",
    "numeric_node_id_tie_break",
    "input_immutable"
  ),
  observed = c(
    as.character(Sys.getlocale("LC_COLLATE")),
    as.character(tree_valid),
    paste(vapply(tip_labels[1:2], escape_text, character(1L)), collapse = ","),
    paste(vapply(tip_labels[3:4], escape_text, character(1L)), collapse = ","),
    escape_text(if (length(keys) >= 1L) keys[[1L]] else NA_character_),
    escape_text(if (length(keys) >= 2L) keys[[2L]] else NA_character_),
    as.character(keys_equal),
    as.character(keys_expected),
    as.character(metadata_safe),
    collapse_integer(old_order),
    collapse_integer(new_order),
    as.character(numeric_tie_break),
    as.character(input_immutable)
  ),
  expected = c(
    "C",
    "TRUE",
    "a,b\\rc",
    "a\\rb,c",
    "a\\rb\\rc",
    "a\\rb\\rc",
    "TRUE",
    "TRUE",
    "TRUE",
    "7,5,6",
    "5,6,7",
    "TRUE",
    "TRUE"
  ),
  pass = c(
    identical(as.character(Sys.getlocale("LC_COLLATE")), "C"),
    tree_valid,
    identical(tip_labels[1:2], c("a", "b\rc")),
    identical(tip_labels[3:4], c("a\rb", "c")),
    length(keys) >= 1L && identical(as.character(keys[[1L]]), expected_key),
    length(keys) >= 2L && identical(as.character(keys[[2L]]), expected_key),
    keys_equal,
    keys_expected,
    metadata_safe,
    same_integer(old_order, expected_old_order),
    same_integer(new_order, expected_new_order),
    numeric_tie_break,
    input_immutable
  ),
  stringsAsFactors = FALSE
)

provenance <- rbind(
  data.frame(
    check = c("source_commit", "package_version", "R_version", "platform",
              "package_load", "protocol"),
    observed = c(
      source_commit,
      package_version,
      R.version.string,
      R.version$platform,
      package_load,
      "non-benchmark smoke; frozen descendant keys; safe canonicalization; serialized input check"
    ),
    expected = c("recorded", "0.2.0.9000", "recorded", "recorded", "recorded",
                 "recorded"),
    pass = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  ),
  checks
)
csv_path <- file.path(out_dir, "ordering_contract.csv")
txt_path <- file.path(out_dir, "ordering_contract.txt")
write.csv(provenance, csv_path, row.names = FALSE, na = "")

overall_pass <- all(provenance$pass)
txt <- c(
  "fastphylosig 0.2.0 Stage 2B2B Candidate 2 ordering-contract smoke",
  paste0("status=", if (overall_pass) "PASS" else "FAIL"),
  paste0("source_commit=", source_commit),
  paste0("package_version=", package_version),
  paste0("LC_COLLATE=", Sys.getlocale("LC_COLLATE")),
  "clade_A={a,b\\rc}",
  "clade_B={a\\rb,c}",
  paste0("clade_A_key=", escape_text(if (length(keys) >= 1L) keys[[1L]] else NA_character_)),
  paste0("clade_B_key=", escape_text(if (length(keys) >= 2L) keys[[2L]] else NA_character_)),
  paste0("collapsed_keys_identical=", keys_equal),
  paste0("old_internal_order=", collapse_integer(old_order)),
  paste0("new_internal_order=", collapse_integer(new_order)),
  paste0("numeric_node_id_tie_break=", numeric_tie_break),
  paste0("input_immutable=", input_immutable),
  paste0("csv=", normalizePath(csv_path, winslash = "/", mustWork = FALSE)),
  paste0("txt=", normalizePath(txt_path, winslash = "/", mustWork = FALSE))
)
writeLines(txt, txt_path, useBytes = TRUE)

if (!overall_pass) {
  stop("ordering-contract smoke failed; inspect ordering_contract.csv.",
       call. = FALSE)
}
message("ordering-contract smoke PASS")
message("CSV: ", normalizePath(csv_path, winslash = "/", mustWork = FALSE))
message("TXT: ", normalizePath(txt_path, winslash = "/", mustWork = FALSE))
