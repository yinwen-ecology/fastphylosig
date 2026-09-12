#!/usr/bin/env Rscript

# Stage 2E1A is an audit-only D profiling runner.  It never edits package
# source, tests, generated package artifacts, or public bindings.
#
# The runner is deliberately fail-closed: the controlled Brownian fixture and
# the complete D accounting path must agree with production before any formal
# timing is accepted.

argv <- commandArgs(trailingOnly = TRUE)

parse_args <- function(values) {
  modes <- intersect(values, c("--smoke", "--formal"))
  if (length(modes) != 1L) {
    stop("Use exactly one of --smoke or --formal.", call. = FALSE)
  }
  positional <- values[!values %in% modes]
  if (length(positional) > 2L) {
    stop("usage: run_stage2e1a_benchmark.R [--smoke|--formal] [repo] [out].",
         call. = FALSE)
  }
  repo <- if (length(positional)) positional[[1L]] else getwd()
  out <- if (length(positional) >= 2L) positional[[2L]] else
    file.path(repo, "benchmarks", "stage2e1a", "results", sub("^--", "", modes[[1L]]))
  list(
    mode = sub("^--", "", modes[[1L]]),
    repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

args <- parse_args(argv)
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_ints <- function(name, default, minimum = 1L) {
  value <- env_or(name, paste(default, collapse = ","))
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum, ".",
         call. = FALSE)
  }
  unique(out)
}

parse_reals <- function(name, default, minimum = 0, maximum = 1) {
  value <- env_or(name, paste(default, collapse = ","))
  out <- suppressWarnings(as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(!is.finite(out)) ||
      any(out < minimum) || any(out > maximum)) {
    stop(name, " must contain comma-separated finite values in [",
         minimum, ",", maximum, "].", call. = FALSE)
  }
  unique(out)
}

parse_shapes <- function(default) {
  value <- env_or("FASTPHYLOSIG_STAGE2E1A_SHAPES", paste(default, collapse = ","))
  out <- unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!nzchar(out)) ||
      any(!out %in% c("balanced", "random", "pectinate"))) {
    stop("FASTPHYLOSIG_STAGE2E1A_SHAPES must contain only balanced, random, or pectinate.",
         call. = FALSE)
  }
  out
}

parse_repeats <- function() {
  value <- suppressWarnings(as.integer(env_or(
    "FASTPHYLOSIG_STAGE2E1A_REPEATS",
    if (args$mode == "formal") "10" else "2"
  )))
  if (length(value) != 1L || is.na(value) || value < 1L ||
      (args$mode == "formal" && value < 10L)) {
    stop("Formal Stage 2E1A timing requires at least 10 paired repeats.",
         call. = FALSE)
  }
  value
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(paste(as.character(openssl::sha256(file(path))), collapse = ""))
  }
  NA_character_
}

sha256_text <- function(text) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(text, algo = "sha256"))
  }
  NA_character_
}

git_head <- function(repo) {
  value <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

clock <- function() unname(as.numeric(proc.time()[["elapsed"]]))

with_seed <- function(seed, fun) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(fun)()
}

capture <- function(fun, timed = FALSE) {
  warnings <- character()
  value <- NULL
  status <- "ok"
  error_message <- ""
  started <- if (timed) clock() else NA_real_
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- if (timed) max(0, clock() - started) else NA_real_
  list(
    value = value, elapsed_s = elapsed, status = status,
    warning_text = paste(unique(warnings), collapse = " | "),
    warning_count = length(warnings), error_message = error_message
  )
}

timed_call <- function(fun) {
  gc(FALSE)
  capture(fun, timed = TRUE)
}

failure_row <- function(label, message) {
  data.frame(check = label, status = "FAIL", detail = as.character(message),
             stringsAsFactors = FALSE)
}

load_package <- function(repo) {
  library_dir <- env_or("FASTPHYLOSIG_STAGE2E1A_LIBRARY")
  if (nzchar(library_dir)) {
    library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    suppressPackageStartupMessages(library(
      "fastphylosig", lib.loc = library_dir, character.only = TRUE
    ))
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required when FASTPHYLOSIG_STAGE2E1A_LIBRARY is unset.",
           call. = FALSE)
    }
    pkgload::load_all(repo, quiet = TRUE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required for Stage 2E1A fixtures.", call. = FALSE)
  }
  asNamespace("fastphylosig")
}

required_symbols <- c(
  "fast_d", "prepare_tree", ".prepare_analysis", ".binary_state",
  ".phylo_d_random_states", ".phylo_d_brownian_states", ".d_null_summary",
  "phylo_d_sums_cpp", "brownian_tree_threshold_cpp"
)

check_required_symbols <- function(ns) {
  missing <- required_symbols[!vapply(required_symbols, function(nm) {
    exists(nm, envir = ns, inherits = FALSE) &&
      (is.function(get(nm, envir = ns, inherits = FALSE)) || nm %in%
         c("phylo_d_sums_cpp", "brownian_tree_threshold_cpp"))
  }, logical(1))]
  if (length(missing)) stop("Missing production symbols: ", paste(missing, collapse = ", "),
                            call. = FALSE)
  invisible(TRUE)
}

# This translation unit is audit-only.  The traversal, RNG calls, full sort,
# and type-7 arithmetic intentionally mirror the current production helper;
# it exposes only phase-separated probes and is compiled in a temporary file.
audit_cpp_code <- paste(c(
  "// [[Rcpp::depends(RcppArmadillo)]]",
  "// [[Rcpp::plugins(cpp11)]]",
  "#include <RcppArmadillo.h>",
  "#include <algorithm>",
  "#include <cmath>",
  "#include <limits>",
  "#include <vector>",
  "using namespace Rcpp;",
  "namespace stage2e1a {",
  "inline double type7(const std::vector<double>& values, double probability) {",
  "  const int n = static_cast<int>(values.size());",
  "  if (n < 1) return std::numeric_limits<double>::quiet_NaN();",
  "  const double h = 1.0 + static_cast<double>(n - 1) * probability;",
  "  const int lo = static_cast<int>(std::floor(h));",
  "  const double gamma = h - static_cast<double>(lo);",
  "  if (lo <= 1) return values.front();",
  "  if (lo >= n) return values.back();",
  "  return (1.0 - gamma) * values[static_cast<std::size_t>(lo - 1)] +",
  "    gamma * values[static_cast<std::size_t>(lo)];",
  "}",
  "struct PreparedTree {",
  "  int n_tip; int max_node; int root;",
  "  std::vector< std::vector<int> > children;",
  "  std::vector< std::vector<double> > lengths;",
  "};",
  "PreparedTree prepare(const IntegerMatrix& edge, const arma::vec& edge_length, int n_tip) {",
  "  const int n_edge = edge.nrow();",
  "  if (edge.ncol() != 2 || edge_length.n_elem != static_cast<arma::uword>(n_edge))",
  "    stop(\"edge must be a two-column matrix matching edge_length.\");",
  "  if (n_tip < 1) stop(\"n_tip must be positive.\");",
  "  int max_node = n_tip;",
  "  for (int i = 0; i < n_edge; ++i) {",
  "    max_node = std::max(max_node, edge(i, 0));",
  "    max_node = std::max(max_node, edge(i, 1));",
  "    if (edge_length(i) < 0.0 || !std::isfinite(edge_length(i)))",
  "      stop(\"edge_length must be finite and non-negative.\");",
  "  }",
  "  std::vector< std::vector<int> > children(max_node + 1);",
  "  std::vector< std::vector<double> > lengths(max_node + 1);",
  "  std::vector<int> is_child(max_node + 1, 0);",
  "  for (int i = 0; i < n_edge; ++i) {",
  "    const int parent = edge(i, 0); const int child = edge(i, 1);",
  "    if (parent < 1 || child < 1 || parent > max_node || child > max_node)",
  "      stop(\"edge contains node numbers outside the expected range.\");",
  "    children[parent].push_back(child); lengths[parent].push_back(edge_length(i));",
  "    is_child[child] = 1;",
  "  }",
  "  int root = n_tip + 1;",
  "  if (root > max_node || children[root].empty() || is_child[root]) {",
  "    root = -1;",
  "    for (int node = n_tip + 1; node <= max_node; ++node)",
  "      if (!children[node].empty() && !is_child[node]) { root = node; break; }",
  "  }",
  "  if (root < 1) stop(\"Could not identify the root node.\");",
  "  PreparedTree out; out.n_tip = n_tip; out.max_node = max_node; out.root = root;",
  "  out.children.swap(children); out.lengths.swap(lengths); return out;",
  "}",
  "}",
  "",
  "// [[Rcpp::export]]",
  "arma::mat stage2e1a_brownian_continuous_cpp(const IntegerMatrix& edge,",
  "    const arma::vec& edge_length, const int n_tip, const int nsim) {",
  "  if (nsim < 1) stop(\"nsim must be positive.\");",
  "  stage2e1a::PreparedTree tree = stage2e1a::prepare(edge, edge_length, n_tip);",
  "  arma::mat out(n_tip, nsim); std::vector<double> node_value(tree.max_node + 1, 0.0);",
  "  std::vector<int> stack; stack.reserve(tree.max_node);",
  "  for (int sim = 0; sim < nsim; ++sim) {",
  "    std::fill(node_value.begin(), node_value.end(), 0.0); stack.clear(); stack.push_back(tree.root);",
  "    while (!stack.empty()) {",
  "      const int parent = stack.back(); stack.pop_back();",
  "      for (std::size_t k = 0; k < tree.children[parent].size(); ++k) {",
  "        const int child = tree.children[parent][k]; const double bl = tree.lengths[parent][k];",
  "        node_value[child] = node_value[parent] + std::sqrt(bl) * R::rnorm(0.0, 1.0);",
  "        if (!tree.children[child].empty()) stack.push_back(child);",
  "      }",
  "    }",
  "    for (int tip = 0; tip < n_tip; ++tip) out(tip, sim) = node_value[tip + 1];",
  "  }",
  "  return out;",
  "}",
  "",
  "// [[Rcpp::export]]",
  "IntegerMatrix stage2e1a_brownian_threshold_cpp(const arma::mat& samples, double prop_state1) {",
  "  const int n = static_cast<int>(samples.n_rows); const int p = static_cast<int>(samples.n_cols);",
  "  if (prop_state1 < 0.0 || prop_state1 > 1.0) stop(\"prop_state1 must be in [0, 1].\");",
  "  IntegerMatrix out(n, p); std::vector<double> values(n);",
  "  for (int j = 0; j < p; ++j) {",
  "    for (int i = 0; i < n; ++i) values[i] = samples(i, j);",
  "    std::sort(values.begin(), values.end());",
  "    const double threshold = stage2e1a::type7(values, prop_state1);",
  "    for (int i = 0; i < n; ++i) out(i, j) = samples(i, j) < threshold ? 1 : 0;",
  "  }",
  "  return out;",
  "}",
  "",
  "// [[Rcpp::export]]",
  "arma::mat stage2e1a_sort_only_cpp(const arma::mat& samples) {",
  "  const int n = static_cast<int>(samples.n_rows); const int p = static_cast<int>(samples.n_cols);",
  "  arma::mat out(n, p); std::vector<double> values(n);",
  "  for (int j = 0; j < p; ++j) {",
  "    for (int i = 0; i < n; ++i) values[i] = samples(i, j);",
  "    std::sort(values.begin(), values.end());",
  "    for (int i = 0; i < n; ++i) out(i, j) = values[i];",
  "  }",
  "  return out;",
  "}",
  "",
  "// [[Rcpp::export]]",
  "Rcpp::NumericVector stage2e1a_sort_threshold_cpp(const arma::mat& samples, double prop_state1) {",
  "  const int n = static_cast<int>(samples.n_rows); const int p = static_cast<int>(samples.n_cols);",
  "  if (prop_state1 < 0.0 || prop_state1 > 1.0) stop(\"prop_state1 must be in [0, 1].\");",
  "  NumericVector out(p); std::vector<double> values(n);",
  "  for (int j = 0; j < p; ++j) {",
  "    for (int i = 0; i < n; ++i) values[i] = samples(i, j);",
  "    std::sort(values.begin(), values.end()); out[j] = stage2e1a::type7(values, prop_state1);",
  "  }",
  "  return out;",
  "}",
  "",
  "// [[Rcpp::export]]",
  "IntegerMatrix stage2e1a_threshold_binary_cpp(const arma::mat& samples, const NumericVector& thresholds) {",
  "  const int n = static_cast<int>(samples.n_rows); const int p = static_cast<int>(samples.n_cols);",
  "  if (thresholds.size() != p) stop(\"thresholds must have one value per column.\");",
  "  IntegerMatrix out(n, p);",
  "  for (int j = 0; j < p; ++j) for (int i = 0; i < n; ++i)",
  "    out(i, j) = samples(i, j) < thresholds[j] ? 1 : 0;",
  "  return out;",
  "}"
), sep = "\n")

load_audit_cpp <- function() {
  if (!requireNamespace("Rcpp", quietly = TRUE) ||
      !requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Rcpp and RcppArmadillo are required for the audit harness.", call. = FALSE)
  }
  cpp_file <- tempfile("fastphylosig-stage2e1a-", fileext = ".cpp")
  writeLines(audit_cpp_code, cpp_file, useBytes = TRUE)
  hook_env <- new.env(parent = globalenv())
  Rcpp::sourceCpp(cpp_file, env = hook_env, rebuild = TRUE,
                  showOutput = FALSE, verbose = FALSE)
  exports <- c(
    "stage2e1a_brownian_continuous_cpp", "stage2e1a_brownian_threshold_cpp",
    "stage2e1a_sort_only_cpp", "stage2e1a_sort_threshold_cpp",
    "stage2e1a_threshold_binary_cpp"
  )
  if (!all(vapply(exports, exists, logical(1), envir = hook_env, inherits = FALSE))) {
    stop("The Stage 2E1A audit translation unit did not expose all probes.", call. = FALSE)
  }
  list(env = hook_env, cpp_file = cpp_file, sha256 = sha256_text(audit_cpp_code),
       continuous = get("stage2e1a_brownian_continuous_cpp", hook_env),
       threshold = get("stage2e1a_brownian_threshold_cpp", hook_env),
       sort_only = get("stage2e1a_sort_only_cpp", hook_env),
       sort_threshold = get("stage2e1a_sort_threshold_cpp", hook_env),
       threshold_binary = get("stage2e1a_threshold_binary_cpp", hook_env))
}

make_balanced_tree <- function(n) {
  root <- n + 1L
  internal <- n - 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  q_node <- q_lo <- q_hi <- integer(internal)
  q_node[[1L]] <- root
  q_lo[[1L]] <- 1L
  q_hi[[1L]] <- n
  head <- 1L
  tail <- 1L
  next_internal <- root + 1L
  edge_i <- 0L
  while (head <= tail) {
    node <- q_node[[head]]
    lo <- q_lo[[head]]
    hi <- q_hi[[head]]
    head <- head + 1L
    mid <- floor((lo + hi) / 2L)
    child_lo <- c(lo, mid + 1L)
    child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) child_lo[[j]] else {
        new_node <- next_internal
        next_internal <- next_internal + 1L
        tail <- tail + 1L
        q_node[[tail]] <- new_node
        q_lo[[tail]] <- child_lo[[j]]
        q_hi[[tail]] <- child_hi[[j]]
        new_node
      }
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(node, child)
    }
  }
  tree <- list(edge = edge[seq_len(edge_i), , drop = FALSE],
               tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(edge_i) / edge_i,
               Nnode = internal)
  class(tree) <- "phylo"
  tree
}

make_pectinate_tree <- function(n) {
  root <- n + 1L
  if (n == 2L) {
    edge <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  } else {
    edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
    current <- root
    next_internal <- root + 1L
    edge_i <- 0L
    for (tip in seq_len(n - 2L)) {
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(current, tip)
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(current, next_internal)
      current <- next_internal
      next_internal <- next_internal + 1L
    }
    edge_i <- edge_i + 1L
    edge[edge_i, ] <- c(current, n - 1L)
    edge_i <- edge_i + 1L
    edge[edge_i, ] <- c(current, n)
    edge <- edge[seq_len(edge_i), , drop = FALSE]
  }
  tree <- list(edge = edge, tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
               Nnode = n - 1L)
  class(tree) <- "phylo"
  tree
}

make_tree <- function(n, shape) {
  tree <- if (shape == "balanced") {
    make_balanced_tree(n)
  } else if (shape == "pectinate") {
    make_pectinate_tree(n)
  } else {
    with_seed(20260908L + as.integer(n), function() ape::rtree(n))
  }
  tree <- ape::reorder.phylo(tree, order = "postorder")
  tree$tip.label <- paste0("sp", seq_len(n))
  tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
  tree
}

make_trait <- function(tree, prevalence, seed) {
  n <- length(tree$tip.label)
  n_first <- max(1L, min(n - 1L, as.integer(round(n * prevalence))))
  values <- c(rep(0, n_first), rep(1, n - n_first))
  values <- with_seed(seed, function() values[sample.int(n)])
  stats::setNames(values, tree$tip.label)
}

make_permutations <- function(n, nsim, seed) {
  with_seed(seed, function() {
    out <- matrix(NA_integer_, nrow = nsim, ncol = n)
    out[1L, ] <- seq_len(n)
    if (nsim >= 2L) out[2L, ] <- rev(seq_len(n))
    if (nsim >= 3L) for (i in 3:nsim) out[i, ] <- sample.int(n)
    out
  })
}

make_direct_fixture <- function(tree, x, ns) {
  ctx <- ns$prepare_tree(tree)
  X <- matrix(as.numeric(x), ncol = 1L,
              dimnames = list(names(x), "trait1"))
  analysis <- ns$.prepare_analysis(ctx, X, signal = "D", data_kind = "binary",
                                   verbose = FALSE)
  if (!length(analysis$groups) || is.null(analysis$groups[[1L]]$group)) {
    stop("The prepared D fixture did not produce a retained group.", call. = FALSE)
  }
  group <- analysis$groups[[1L]]$group
  d <- ns$.binary_state(stats::setNames(X[, 1L], rownames(X)))
  phy <- group$d_tree
  edge <- matrix(as.integer(phy$edge), ncol = 2L)
  edge_length <- as.numeric(phy$edge.length)
  ds <- as.numeric(d$values[phy$tip.label])
  list(tree = tree, ctx = ctx, x = x, X = X, analysis = analysis,
       group = group, phy = phy, edge = edge, edge_length = edge_length,
       n_tip = length(phy$tip.label), ds = ds, prop_state1 = d$prop_state1,
       observed = matrix(ds, ncol = 1L), prepare_time_s = NA_real_)
}

prepare_fixture_timed <- function(tree, x, ns) {
  started <- clock()
  out <- make_direct_fixture(tree, x, ns)
  out$prepare_time_s <- max(0, clock() - started)
  out
}

stable_d <- function(value) {
  if (is.null(value)) return(NULL)
  list(
    DEstimate = value$DEstimate, Pval1 = value$Pval1, Pval0 = value$Pval0,
    P_random = value$P_random, P_Brownian = value$P_Brownian,
    MCSE_P_random = value$MCSE_P_random, MCSE_P_Brownian = value$MCSE_P_Brownian,
    Parameters = value$Parameters, nsim_requested = value$nsim_requested,
    nsim_successful_random = value$nsim_successful_random,
    nsim_successful_brownian = value$nsim_successful_brownian,
    nsim_failed_random = value$nsim_failed_random,
    nsim_failed_brownian = value$nsim_failed_brownian,
    status = value$status, note = value$note
  )
}

run_correctness <- function(ns, hooks, out_dir) {
  rows <- list()
  add <- function(label, pass, detail = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = label, status = if (isTRUE(pass)) "PASS" else "FAIL",
      detail = as.character(detail), stringsAsFactors = FALSE
    )
  }
  safe_identical <- function(a, b) isTRUE(identical(a, b))
  # ape::stree(type = "balanced") requires a power-of-two tip count.
  tree <- make_tree(16L, "balanced")
  x <- make_trait(tree, 0.5, 20260911L)
  fx <- tryCatch(make_direct_fixture(tree, x, ns), error = function(e) e)
  if (inherits(fx, "error")) {
    add("prepared_fixture", FALSE, conditionMessage(fx))
    result <- do.call(rbind, rows)
    status <- data.frame(status = "FAIL", stage = "Stage 2E1A correctness gate",
                         checks = nrow(result), failed = sum(result$status == "FAIL"),
                         production_code_changed = "NO", reason = "fixture construction failed",
                         stringsAsFactors = FALSE)
    write_csv(result, file.path(out_dir, "stage2e1a_correctness_checks.csv"))
    write_csv(status, file.path(out_dir, "stage2e1a_correctness_status.csv"))
    stop("Stage 2E1A correctness gate failed.", call. = FALSE)
  }

  nsim <- 37L
  permutations <- make_permutations(fx$n_tip, nsim, 20260912L)
  brown_seed <- 20260913L
  continuous <- with_seed(brown_seed, function() hooks$continuous(
    fx$edge, fx$edge_length, fx$n_tip, nsim
  ))
  production_brownian <- with_seed(brown_seed, function() ns$brownian_tree_threshold_cpp(
    edge = fx$edge, edge_length = fx$edge_length, n_tip = fx$n_tip,
    nsim = nsim, prop_state1 = fx$prop_state1
  ))
  audit_brownian <- hooks$threshold(continuous, fx$prop_state1)
  add("brownian_continuous_rng_replay", is.matrix(continuous) &&
        nrow(continuous) == fx$n_tip && ncol(continuous) == nsim)
  add("brownian_threshold_exact", safe_identical(production_brownian, audit_brownian))

  sorted <- hooks$sort_only(continuous)
  thresholds <- hooks$sort_threshold(continuous, fx$prop_state1)
  converted <- hooks$threshold_binary(continuous, thresholds)
  add("sort_only_shape", is.matrix(sorted) && all(dim(sorted) == dim(continuous)))
  add("sort_threshold_shape", length(thresholds) == nsim && all(is.finite(thresholds)))
  add("threshold_binary_exact", safe_identical(converted, audit_brownian))

  random_states <- ns$.phylo_d_random_states(
    fx$ds, nsim, rnd.bias = NULL, permutations = permutations,
    random_states = NULL
  )
  brownian_states <- ns$.phylo_d_brownian_states(
    fx$group$tree, nsim, fx$prop_state1, production_brownian
  )
  sums <- ns$phylo_d_sums_cpp(
    cbind(fx$observed, random_states, brownian_states), edge = fx$edge,
    edge_length = fx$edge_length, n_tip = fx$n_tip, n_threads = 1L
  )
  random_summary <- ns$.d_null_summary(
    sums[seq.int(2L, nsim + 1L)], sums[[1L]], direction = "less"
  )
  brownian_summary <- ns$.d_null_summary(
    sums[seq.int(nsim + 2L, 2L * nsim + 1L)], sums[[1L]], direction = "greater"
  )
  expected_d <- (sums[[1L]] - brownian_summary$mean) /
    (random_summary$mean - brownian_summary$mean)

  controlled <- capture(function() ns$fast_d(
    fx$ctx, x, test = TRUE, nsim = nsim, permutations = permutations,
    brownian_states = production_brownian, return_sim = TRUE, keep_null = TRUE,
    verbose = FALSE, progress = FALSE, ncores = 1L
  ))
  add("controlled_fast_d_status", controlled$status == "ok",
      paste(controlled$error_message, controlled$warning_text))
  actual <- controlled$value
  if (controlled$status == "ok") {
    add("observed_exact", safe_identical(actual$Parameters$Observed, sums[[1L]]))
    add("random_mean_exact", safe_identical(actual$Parameters$MeanRandom, random_summary$mean))
    add("brownian_mean_exact", safe_identical(actual$Parameters$MeanBrownian, brownian_summary$mean))
    add("D_exact", safe_identical(actual$DEstimate, expected_d))
    add("P_random_exact", safe_identical(actual$P_random, random_summary$p))
    add("P_brownian_exact", safe_identical(actual$P_Brownian, brownian_summary$p))
    add("MCSE_random_exact", safe_identical(actual$MCSE_P_random, random_summary$mcse))
    add("MCSE_brownian_exact", safe_identical(actual$MCSE_P_Brownian, brownian_summary$mcse))
    add("successful_accounting_exact",
        safe_identical(actual$nsim_successful_random, nsim) &&
          safe_identical(actual$nsim_successful_brownian, nsim))
    add("failed_accounting_exact",
        safe_identical(actual$nsim_failed_random, 0L) &&
          safe_identical(actual$nsim_failed_brownian, 0L))
    add("random_null_exact", safe_identical(as.numeric(actual$Permutations$random),
                                             as.numeric(sums[seq.int(2L, nsim + 1L)])))
    add("brownian_null_exact", safe_identical(as.numeric(actual$Permutations$brownian),
                                               as.numeric(sums[seq.int(nsim + 2L, 2L * nsim + 1L)])))
  }

  replay_one <- capture(function() with_seed(20260914L, function() ns$fast_d(
    fx$ctx, x, test = TRUE, nsim = 19L, return_sim = FALSE, keep_null = FALSE,
    verbose = FALSE, progress = FALSE, ncores = 1L
  )))
  replay_two <- capture(function() with_seed(20260914L, function() ns$fast_d(
    fx$ctx, x, test = TRUE, nsim = 19L, return_sim = FALSE, keep_null = FALSE,
    verbose = FALSE, progress = FALSE, ncores = 1L
  )))
  add("rng_replay_status", replay_one$status == "ok" && replay_two$status == "ok",
      paste(replay_one$error_message, replay_two$error_message))
  if (replay_one$status == "ok" && replay_two$status == "ok") {
    add("rng_replay_exact", safe_identical(stable_d(replay_one$value),
                                            stable_d(replay_two$value)))
  }

  strict <- ns$.d_null_summary(c(0, 2, 2), 1, direction = "greater")
  add("strict_tail_definition", safe_identical(strict$p, 2 / 3) &&
        safe_identical(strict$n_successful, 3L) &&
        safe_identical(strict$n_failed, 0L))
  checks <- do.call(rbind, rows)
  status_value <- if (all(checks$status == "PASS")) "PASS" else "FAIL"
  status <- data.frame(
    status = status_value, stage = "Stage 2E1A correctness gate",
    checks = nrow(checks), failed = sum(checks$status == "FAIL"),
    production_code_changed = "NO",
    reason = if (status_value == "PASS") "controlled production-equivalence gate passed" else
      "one or more controlled equivalence checks failed",
    stringsAsFactors = FALSE
  )
  write_csv(checks, file.path(out_dir, "stage2e1a_correctness_checks.csv"))
  write_csv(status, file.path(out_dir, "stage2e1a_correctness_status.csv"))
  if (status_value != "PASS") stop("Stage 2E1A correctness gate failed; timings are not authorized.",
                                    call. = FALSE)
  invisible(status)
}

source_type7_audit <- function(repo) {
  tryCatch({
    d_stream <- file.path(repo, "src", "d_stream.cpp")
    fast_signal <- file.path(repo, "src", "fast_signal.cpp")
    numeric_utils <- file.path(repo, "src", "numeric_utils.h")
    if (!file.exists(d_stream) || !file.exists(fast_signal) ||
        !file.exists(numeric_utils)) {
      return(list(status = "NOT_RUN", detail = "production source snapshot unavailable"))
    }
    exported_block <- function(path, symbol) {
      lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
      start <- grep(paste0(symbol, "("), lines, fixed = TRUE)
      start <- if (length(start)) start[[1L]] else NA_integer_
      next_export <- grep("// [[Rcpp::export]]", lines, fixed = TRUE)
      next_export <- next_export[next_export > start]
      finish <- if (length(next_export)) next_export[[1L]] - 1L else length(lines)
      if (is.finite(start)) paste(lines[start:finish], collapse = "\n") else ""
    }
    stream_block <- exported_block(d_stream, "phylo_d_stream_cpp")
    fallback_block <- exported_block(fast_signal, "brownian_tree_threshold_cpp")
    q_lines <- readLines(numeric_utils, warn = FALSE, encoding = "UTF-8")
    q_block <- paste(q_lines, collapse = "\n")
    count_fixed <- function(block, pattern) {
      sum(grepl(pattern, strsplit(block, "\n", fixed = TRUE)[[1L]], fixed = TRUE))
    }
    stream_sort <- count_fixed(stream_block, "std::sort(sorted_values")
    stream_quantile <- count_fixed(stream_block, "quantile_type7_sorted")
    fallback_sort <- count_fixed(fallback_block, "std::sort(sorted_values")
    fallback_quantile <- count_fixed(fallback_block, "quantile_type7_sorted")
    two_adjacent <- grepl("values\\[static_cast<std::size_t>\\(lo - 1\\)\\]", q_block) &&
      grepl("values\\[static_cast<std::size_t>\\(lo\\)\\]", q_block)
    stream_strict <- grepl("brownian_tip_values[static_cast<std::size_t>(tip)] < threshold",
                           stream_block, fixed = TRUE)
    fallback_strict <- grepl("tip_values[tip] < threshold", fallback_block,
                             fixed = TRUE)
    pass <- stream_sort == 1L && stream_quantile == 1L &&
      fallback_sort == 1L && fallback_quantile == 1L && two_adjacent &&
      stream_strict && fallback_strict
    list(
      status = if (pass) "PASS" else "FAIL",
      detail = if (pass) {
        "Streaming D and rooted-polytomy fallback each sort once and consume only the two adjacent type-7 order statistics before strict-threshold binary conversion."
      } else {
        "Static type-7/order-statistic contract could not be established from the source snapshot."
      }, stream_sort_once = stream_sort == 1L,
      stream_quantile_once = stream_quantile == 1L,
      fallback_sort_once = fallback_sort == 1L,
      fallback_quantile_once = fallback_quantile == 1L,
      strict_binary_threshold = stream_strict && fallback_strict,
      adjacent_order_statistics = two_adjacent,
      d_stream_sha256 = sha256_file(d_stream),
      fast_signal_sha256 = sha256_file(fast_signal),
      numeric_utils_sha256 = sha256_file(numeric_utils)
    )
  }, error = function(e) list(status = "NOT_RUN",
                               detail = paste0("source audit unavailable: ", conditionMessage(e))))
}

chunk_lengths <- function(nsim, chunk_size) {
  starts <- seq.int(1L, nsim, by = chunk_size)
  pmin(chunk_size, nsim - starts + 1L)
}

phase_capture_sum <- function(calls) {
  elapsed <- 0
  statuses <- character()
  warnings <- character()
  errors <- character()
  for (call in calls) {
    one <- timed_call(call)
    elapsed <- elapsed + if (is.finite(one$elapsed_s)) one$elapsed_s else 0
    statuses <- c(statuses, one$status)
    if (nzchar(one$warning_text)) warnings <- c(warnings, one$warning_text)
    if (nzchar(one$error_message)) errors <- c(errors, one$error_message)
  }
  list(elapsed_s = elapsed, status = if (all(statuses == "ok")) "ok" else "error",
       warning_text = paste(unique(warnings), collapse = " | "),
       error_message = paste(unique(errors), collapse = " | "))
}

phase_random <- function(fx, ns, nsim, chunk_size, seed) {
  with_seed(seed, function() {
    lens <- chunk_lengths(nsim, chunk_size)
    phase_capture_sum(lapply(lens, function(k) function() {
      invisible(ns$.phylo_d_random_states(fx$ds, k, NULL, NULL, NULL))
    }))
  })
}

phase_brownian_generation <- function(fx, hooks, nsim, chunk_size, seed) {
  with_seed(seed, function() {
    lens <- chunk_lengths(nsim, chunk_size)
    phase_capture_sum(lapply(lens, function(k) function() {
      invisible(hooks$continuous(fx$edge, fx$edge_length, fx$n_tip, k))
    }))
  })
}

phase_on_samples <- function(fx, hooks, nsim, chunk_size, seed, operation) {
  with_seed(seed, function() {
    lens <- chunk_lengths(nsim, chunk_size)
    elapsed <- 0
    statuses <- character()
    warnings <- character()
    errors <- character()
    for (k in lens) {
      generated <- capture(function() hooks$continuous(
        fx$edge, fx$edge_length, fx$n_tip, k
      ))
      if (generated$status != "ok") {
        statuses <- c(statuses, "error")
        errors <- c(errors, generated$error_message)
        next
      }
      samples <- generated$value
      one <- if (operation == "sort") {
        timed_call(function() invisible(hooks$sort_only(samples)))
      } else if (operation == "sort_threshold") {
        timed_call(function() invisible(hooks$sort_threshold(samples, fx$prop_state1)))
      } else {
        thresholds <- hooks$sort_threshold(samples, fx$prop_state1)
        timed_call(function() {
          brownian <- hooks$threshold_binary(samples, thresholds)
          invisible(ns_global$phylo_d_sums_cpp(
            brownian, edge = fx$edge, edge_length = fx$edge_length,
            n_tip = fx$n_tip, n_threads = 1L
          ))
        })
      }
      elapsed <- elapsed + if (is.finite(one$elapsed_s)) one$elapsed_s else 0
      statuses <- c(statuses, one$status)
      if (nzchar(one$warning_text)) warnings <- c(warnings, one$warning_text)
      if (nzchar(one$error_message)) errors <- c(errors, one$error_message)
    }
    list(elapsed_s = elapsed, status = if (length(statuses) && all(statuses == "ok")) "ok" else "error",
         warning_text = paste(unique(warnings), collapse = " | "),
         error_message = paste(unique(errors), collapse = " | "))
  })
}

phase_complete_brownian <- function(fx, ns, nsim, seed) {
  with_seed(seed, function() timed_call(function() invisible(ns$.phylo_d_brownian_states(
    fx$group$tree, nsim, fx$prop_state1, brownian_states = NULL
  ))))
}

phase_complete_d <- function(fx, ns, nsim, ncores, seed) {
  with_seed(seed, function() timed_call(function() invisible(ns$fast_d(
    fx$ctx, fx$x, test = TRUE, nsim = nsim, return_sim = FALSE,
    keep_null = FALSE, verbose = FALSE, progress = FALSE, ncores = ncores
  ))))
}

warmup_cell <- function(fx, ns, hooks, nsim, chunk_size) {
  warm_nsim <- min(as.integer(nsim), max(5L, min(19L, as.integer(chunk_size))))
  invisible(phase_random(fx, ns, warm_nsim, chunk_size, 31001L))
  invisible(phase_brownian_generation(fx, hooks, warm_nsim, chunk_size, 31002L))
  invisible(phase_on_samples(fx, hooks, warm_nsim, chunk_size, 31003L, "sort"))
  invisible(phase_on_samples(fx, hooks, warm_nsim, chunk_size, 31004L, "sort_threshold"))
  invisible(phase_on_samples(fx, hooks, warm_nsim, chunk_size, 31005L, "binary_d"))
  invisible(phase_complete_brownian(fx, ns, warm_nsim, 31006L))
  invisible(phase_complete_d(fx, ns, warm_nsim, 1L, 31007L))
  invisible(phase_complete_d(fx, ns, warm_nsim, 2L, 31008L))
  invisible(NULL)
}

run_cell <- function(tree, x, shape, n, prevalence, ns, hooks,
                     nsim, repeats, chunk_size, cell_index) {
  fx <- prepare_fixture_timed(tree, x, ns)
  warmup_cell(fx, ns, hooks, nsim, chunk_size)
  phase_names <- c("random", "brownian_generation", "brownian_sort",
                   "brownian_sort_threshold", "binary_D", "complete_brownian",
                   "complete_D_ncores1", "complete_D_ncores2")
  rows <- vector("list", repeats)
  phase_fun <- list(
    random = function(i) phase_random(fx, ns, nsim, chunk_size,
                                      410000L + cell_index * 1000L + i),
    brownian_generation = function(i) phase_brownian_generation(
      fx, hooks, nsim, chunk_size, 420000L + cell_index * 1000L + i
    ),
    brownian_sort = function(i) phase_on_samples(
      fx, hooks, nsim, chunk_size, 430000L + cell_index * 1000L + i, "sort"
    ),
    brownian_sort_threshold = function(i) phase_on_samples(
      fx, hooks, nsim, chunk_size, 440000L + cell_index * 1000L + i, "sort_threshold"
    ),
    binary_D = function(i) phase_on_samples(
      fx, hooks, nsim, chunk_size, 450000L + cell_index * 1000L + i, "binary_d"
    ),
    complete_brownian = function(i) phase_complete_brownian(
      fx, ns, nsim, 460000L + cell_index * 1000L + i
    ),
    complete_D_ncores1 = function(i) phase_complete_d(
      fx, ns, nsim, 1L, 470000L + cell_index * 1000L + i
    ),
    complete_D_ncores2 = function(i) phase_complete_d(
      fx, ns, nsim, 2L, 480000L + cell_index * 1000L + i
    )
  )
  for (i in seq_len(repeats)) {
    order <- if (i %% 2L) phase_names else rev(phase_names)
    values <- setNames(as.list(rep(NA_real_, length(phase_names))), phase_names)
    status <- character()
    warnings <- character()
    errors <- character()
    for (name in order) {
      one <- phase_fun[[name]](i)
      values[[name]] <- one$elapsed_s
      status <- c(status, one$status)
      if (nzchar(one$warning_text)) warnings <- c(warnings, one$warning_text)
      if (nzchar(one$error_message)) errors <- c(errors, one$error_message)
    }
    rows[[i]] <- data.frame(
      cell_id = paste(shape, n, sprintf("%.2f", prevalence), nsim, sep = "_"),
      shape = shape, n = n, prevalence = prevalence, nsim = nsim,
      `repeat` = i, order = paste(order, collapse = ">"),
      prepare_time_s = fx$prepare_time_s,
      random_null_s = values$random,
      brownian_generation_s = values$brownian_generation,
      brownian_sort_s = values$brownian_sort,
      brownian_sort_threshold_s = values$brownian_sort_threshold,
      binary_D_s = values$binary_D,
      complete_brownian_s = values$complete_brownian,
      complete_D_ncores1_s = values$complete_D_ncores1,
      complete_D_ncores2_s = values$complete_D_ncores2,
      status = if (all(status == "ok")) "PASS" else "FAIL",
      warning_text = paste(unique(warnings), collapse = " | "),
      error_message = paste(unique(errors), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

median_iqr <- function(x) {
  x <- as.numeric(x)
  c(median = if (length(x)) stats::median(x, na.rm = TRUE) else NA_real_,
    IQR = if (length(x)) stats::IQR(x, na.rm = TRUE) else NA_real_)
}

safe_div <- function(a, b) {
  if (!is.finite(a) || !is.finite(b) || b <= 0) NA_real_ else a / b
}

summarize_timings <- function(raw) {
  numeric_names <- c("prepare_time_s", "random_null_s", "brownian_generation_s",
                     "brownian_sort_s", "brownian_sort_threshold_s", "binary_D_s",
                     "complete_brownian_s", "complete_D_ncores1_s", "complete_D_ncores2_s")
  groups <- split(raw, raw[c("shape", "n", "prevalence", "nsim")], drop = TRUE)
  rows <- lapply(groups, function(z) {
    row <- z[1L, c("cell_id", "shape", "n", "prevalence", "nsim"), drop = FALSE]
    for (nm in numeric_names) {
      stats <- median_iqr(z[[nm]])
      row[[paste0(nm, "_median_s")]] <- stats[["median"]]
      row[[paste0(nm, "_iqr_s")]] <- stats[["IQR"]]
    }
    total <- row$complete_D_ncores1_s_median
    row$sort_fraction_of_complete_Brownian <- safe_div(
      row$brownian_sort_s_median, row$complete_brownian_s_median
    )
    row$sort_threshold_fraction_of_total_D <- safe_div(
      row$brownian_sort_threshold_s_median, total
    )
    row$random_fraction_of_total_D <- safe_div(row$random_null_s_median, total)
    row$Brownian_generation_fraction <- safe_div(row$brownian_generation_s_median, total)
    row$binary_D_fraction <- safe_div(row$binary_D_s_median, total)
    row$independent_phase_sum_s <- sum(
      row$random_null_s_median, row$brownian_generation_s_median,
      row$brownian_sort_s_median, row$brownian_sort_threshold_s_median,
      row$binary_D_s_median
    )
    row$sort_total_fraction_for_ceiling <- safe_div(row$brownian_sort_s_median, total)
    sort_share <- row$sort_total_fraction_for_ceiling
    row$theoretical_free_sort_speedup <- if (is.finite(sort_share) && sort_share < 1) {
      1 / (1 - sort_share)
    } else NA_real_
    for (remaining in c(0.2, 0.3, 0.5)) {
      row[[paste0("speedup_old_sort_", remaining * 100, "pct")]] <-
        if (is.finite(sort_share) && sort_share < 1) {
          1 / (1 - sort_share * (1 - remaining))
        } else NA_real_
    }
    row$ncores2_ratio <- safe_div(
      row$complete_D_ncores2_s_median, row$complete_D_ncores1_s_median
    )
    row$timing_status <- if (all(z$status == "PASS")) "PASS" else "FAIL"
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

make_decision <- function(summary, source_audit, raw, repeats, args,
                          formal_grid_complete) {
  heavy <- summary[summary$n >= 2000L & summary$nsim == 9999L &
                     summary$timing_status == "PASS", , drop = FALSE]
  valid_heavy <- nrow(heavy)
  sort_threshold_pass <- if (valid_heavy) {
    sum(heavy$sort_threshold_fraction_of_total_D >= 0.30, na.rm = TRUE)
  } else 0L
  sort_threshold_any <- if (valid_heavy) {
    any(heavy$sort_threshold_fraction_of_total_D >= 0.20, na.rm = TRUE)
  } else FALSE
  source_pass <- identical(source_audit$status, "PASS")
  # No exact-selection implementation is authorized in Stage 2E1A.  Even a
  # large share therefore cannot become GO in this profiling-only stage.
  decision <- if (!isTRUE(formal_grid_complete) || !nrow(summary) ||
                 any(raw$status != "PASS")) {
    "NEED_MORE_EVIDENCE"
  } else if (sort_threshold_pass >= 2L && source_pass) {
    "NEED_MORE_EVIDENCE"
  } else if (sort_threshold_any) {
    "NEED_MORE_EVIDENCE"
  } else {
    "NO_GO"
  }
  optimization <- if (decision == "NO_GO") "CLOSED" else "CONTINUE"
  next_stage <- if (decision == "NO_GO") "FINAL_0_2_0_VALIDATION" else
    "FINAL_DESIGN_REVIEW_NO_AUTO_PROTOTYPE"
  min_speed_50 <- if (nrow(heavy)) min(heavy$speedup_old_sort_50pct, na.rm = TRUE) else NA_real_
  max_free <- if (nrow(heavy)) max(heavy$theoretical_free_sort_speedup, na.rm = TRUE) else NA_real_
  data.frame(
    stage = "Stage 2E1A Final D Candidate Qualification",
    mode = args$mode, source_commit = git_head(args$repo), repeats = repeats,
    formal_grid_complete = identical(args$mode, "formal"),
    correctness_gate = "PASS", source_type7_order_statistic_audit = source_audit$status,
    heavy_valid_cells = valid_heavy, heavy_cells_sort_threshold_ge_30pct = sort_threshold_pass,
    heavy_cells_sort_threshold_ge_20pct = sum(heavy$sort_threshold_fraction_of_total_D >= 0.20,
                                               na.rm = TRUE),
    exact_selection_design = "NOT_IMPLEMENTED_IN_STAGE_2E1A",
    theoretical_free_sort_speedup_max = max_free,
    conservative_speedup_old_sort_50pct_min_heavy = min_speed_50,
    D_EXACT_ORDER_STATISTIC = decision,
    D_OPTIMIZATION_0_2_0 = optimization, NEXT_STAGE = next_stage,
    production_code_changed = "NO", stringsAsFactors = FALSE
  )
}

ns_global <- NULL

ns <- tryCatch(load_package(args$repo), error = function(e) {
  stop("Cannot load audited fastphylosig package: ", conditionMessage(e), call. = FALSE)
})
check_required_symbols(ns)
ns_global <- ns
hooks <- load_audit_cpp()

correctness_dir <- file.path(args$out, "correctness")
dir.create(correctness_dir, recursive = TRUE, showWarnings = FALSE)
run_correctness(ns, hooks, correctness_dir)

source_repo <- normalizePath(
  env_or("FASTPHYLOSIG_STAGE2E1A_SOURCE_REPO", args$repo),
  winslash = "/", mustWork = TRUE
)
source_audit <- source_type7_audit(source_repo)
write_csv(as.data.frame(source_audit, stringsAsFactors = FALSE),
          file.path(args$out, "stage2e1a_source_type7_audit.csv"))

if (args$mode == "formal") {
  n_grid <- parse_ints("FASTPHYLOSIG_STAGE2E1A_N", c(500L, 2000L, 5000L), minimum = 2L)
  if (tolower(env_or("FASTPHYLOSIG_STAGE2E1A_INCLUDE_10000", "false")) %in% c("true", "1", "yes")) {
    n_grid <- unique(c(n_grid, 10000L))
  }
  nsim_grid <- parse_ints("FASTPHYLOSIG_STAGE2E1A_NSIM", c(999L, 9999L), minimum = 1L)
  prevalence_grid <- parse_reals("FASTPHYLOSIG_STAGE2E1A_PREVALENCE", c(0.10, 0.50, 0.90))
  shapes <- parse_shapes(c("balanced", "random", "pectinate"))
} else {
  n_grid <- parse_ints("FASTPHYLOSIG_STAGE2E1A_N", c(50L, 200L), minimum = 2L)
  nsim_grid <- parse_ints("FASTPHYLOSIG_STAGE2E1A_NSIM", c(19L, 99L), minimum = 1L)
  prevalence_grid <- parse_reals("FASTPHYLOSIG_STAGE2E1A_PREVALENCE", c(0.10, 0.50, 0.90))
  shapes <- parse_shapes("balanced")
}
repeats <- parse_repeats()
chunk_size <- suppressWarnings(as.integer(env_or("FASTPHYLOSIG_STAGE2E1A_CHUNK_SIZE", "128")))
if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
  stop("FASTPHYLOSIG_STAGE2E1A_CHUNK_SIZE must be a positive integer.", call. = FALSE)
}

# A narrowed formal run is useful as a diagnostic, but it is not evidence for
# the requested qualification grid.  Keep that distinction explicit in both
# the decision and provenance rows.
formal_grid_complete <- identical(args$mode, "formal") &&
  all(c(500L, 2000L, 5000L) %in% n_grid) &&
  all(c(999L, 9999L) %in% nsim_grid) &&
  all(c(0.10, 0.50, 0.90) %in% prevalence_grid) &&
  all(c("balanced", "random", "pectinate") %in% shapes) &&
  repeats >= 10L

cells <- expand.grid(shape = shapes, n = n_grid, prevalence = prevalence_grid,
                     nsim = nsim_grid, KEEP.OUT.ATTRS = FALSE,
                     stringsAsFactors = FALSE)
raw_rows <- list()
checkpoint_path <- file.path(args$out, "stage2e1a_phase_timings.partial.csv")
checkpoint <- if (file.exists(checkpoint_path)) {
  tryCatch(utils::read.csv(checkpoint_path, stringsAsFactors = FALSE,
                           check.names = FALSE),
           error = function(e) NULL)
} else NULL
cell_index <- 0L
for (i in seq_len(nrow(cells))) {
  cell_index <- cell_index + 1L
  cell <- cells[i, , drop = FALSE]
  cell_id <- paste(cell$shape, cell$n, sprintf("%.2f", cell$prevalence),
                   cell$nsim, sep = "_")
  prior <- if (!is.null(checkpoint) && nrow(checkpoint)) {
    checkpoint[checkpoint$cell_id == cell_id, , drop = FALSE]
  } else NULL
  if (!is.null(prior) && nrow(prior) == repeats &&
      all(prior$status == "PASS")) {
    message(sprintf("[stage2e1a] %s/%s: resume completed %s",
                    cell_index, nrow(cells), cell_id))
    raw_rows[[length(raw_rows) + 1L]] <- prior
    next
  }
  message(sprintf("[stage2e1a] %s/%s: %s n=%s prevalence=%.2f nsim=%s",
                  cell_index, nrow(cells), cell$shape, cell$n,
                  cell$prevalence, cell$nsim))
  tree <- make_tree(as.integer(cell$n), cell$shape)
  x <- make_trait(tree, as.numeric(cell$prevalence),
                  220000L + cell_index)
  one <- tryCatch(run_cell(
    tree, x, cell$shape, as.integer(cell$n), as.numeric(cell$prevalence),
    ns, hooks, as.integer(cell$nsim), repeats, chunk_size, cell_index
  ), error = function(e) {
    data.frame(
      cell_id = paste(cell$shape, cell$n, sprintf("%.2f", cell$prevalence), cell$nsim, sep = "_"),
      shape = cell$shape, n = as.integer(cell$n), prevalence = as.numeric(cell$prevalence),
      nsim = as.integer(cell$nsim), `repeat` = NA_integer_, order = "",
      prepare_time_s = NA_real_, random_null_s = NA_real_, brownian_generation_s = NA_real_,
      brownian_sort_s = NA_real_, brownian_sort_threshold_s = NA_real_, binary_D_s = NA_real_,
      complete_brownian_s = NA_real_, complete_D_ncores1_s = NA_real_, complete_D_ncores2_s = NA_real_,
      status = "FAIL", warning_text = "", error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  raw_rows[[length(raw_rows) + 1L]] <- one
  write_csv(do.call(rbind, raw_rows), checkpoint_path)
  gc(FALSE)
}
raw <- do.call(rbind, raw_rows)
write_csv(raw, file.path(args$out, "stage2e1a_phase_timings.csv"))
summary <- summarize_timings(raw)
write_csv(summary, file.path(args$out, "stage2e1a_phase_summary.csv"))
decision <- make_decision(summary, source_audit, raw, repeats, args,
                          formal_grid_complete)
write_csv(decision, file.path(args$out, "stage2e1a_decision.csv"))

provenance <- data.frame(
  stage = "Stage 2E1A Final D Candidate Qualification", mode = args$mode,
  source_commit = git_head(args$repo), package_version = tryCatch(
    as.character(utils::packageVersion("fastphylosig")), error = function(e) NA_character_
  ), R_version = R.version.string, platform = R.version$platform,
  os = Sys.info()[["sysname"]], machine = Sys.info()[["machine"]],
  source_snapshot_repo = source_repo,
  locale_collate = Sys.getlocale("LC_COLLATE"), omp_num_threads = env_or("OMP_NUM_THREADS"),
  audit_cpp_sha256 = hooks$sha256, audit_cpp_file = hooks$cpp_file,
  correctness_status = "PASS", correctness_status_file = file.path(correctness_dir,
                                                                     "stage2e1a_correctness_status.csv"),
  fixed_fixture = TRUE, warmup = TRUE, serialized = TRUE, alternating_order = TRUE,
  direct_phases_use_bounded_chunks = TRUE, chunk_size = chunk_size,
  n_grid = paste(n_grid, collapse = ","), nsim_grid = paste(nsim_grid, collapse = ","),
  prevalence_grid = paste(prevalence_grid, collapse = ","), shapes = paste(shapes, collapse = ","),
  repeats = repeats, formal_grid_complete = formal_grid_complete,
  authoritative_ncores = 1L, observation_ncores = 2L,
  production_code_changed = "NO", stringsAsFactors = FALSE
)
write_csv(provenance, file.path(args$out, "stage2e1a_provenance.csv"))

overall <- if (nrow(raw) && all(raw$status == "PASS")) "PASS" else "FAIL"
status <- data.frame(
  status = overall, stage = "Stage 2E1A benchmark", mode = args$mode,
  correctness_gate = "PASS", production_code_changed = "NO",
  timing_rows = nrow(raw), timing_failures = sum(raw$status != "PASS"),
  decision = decision$D_EXACT_ORDER_STATISTIC[[1L]],
  optimization = decision$D_OPTIMIZATION_0_2_0[[1L]],
  next_stage = decision$NEXT_STAGE[[1L]],
  reason = if (overall == "PASS") "all requested timing cells completed without warning or error" else
    "one or more timing cells failed; do not interpret affected fractions",
  stringsAsFactors = FALSE
)
write_csv(status, file.path(args$out, "stage2e1a_benchmark_status.csv"))
message("[stage2e1a] complete; status=", overall,
        "; D_EXACT_ORDER_STATISTIC=", decision$D_EXACT_ORDER_STATISTIC[[1L]],
        "; production_code_changed=NO")

if (overall != "PASS") quit(save = "no", status = 1L)
