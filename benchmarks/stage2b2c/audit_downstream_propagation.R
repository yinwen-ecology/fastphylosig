#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2C downstream propagation audit.
#
# This is an audit-only script.  It compares representation-equivalent trees
# through the public diagnostics and estimators; it never edits package code,
# tests, or the supplied tree/data objects.  The output directory is supplied
# explicitly so formal evidence can be kept separate from source files.

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1L && nzchar(args[[1L]])) {
  normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
if (length(args) < 2L || !nzchar(args[[2L]])) {
  stop("usage: audit_downstream_propagation.R <repo> <output-dir>",
       call. = FALSE)
}
out_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
locale_before <- Sys.getlocale("LC_COLLATE")
locale_c <- tryCatch(
  Sys.setlocale("LC_COLLATE", "C"),
  warning = function(w) conditionMessage(w),
  error = function(e) conditionMessage(e)
)
if (!identical(locale_c, "C")) {
  stop(sprintf("LC_COLLATE could not be set to C (got %s).", locale_c),
       call. = FALSE)
}

if (requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
} else {
  stop("pkgload is required for this source audit.", call. = FALSE)
}
suppressPackageStartupMessages(library(ape))

source_commit <- if (length(args) >= 3L && nzchar(args[[3L]])) {
  args[[3L]]
} else {
  Sys.getenv("FASTPHYLOSIG_STAGE2B2C_SOURCE_COMMIT", "unknown")
}
if (!nzchar(source_commit) || identical(source_commit, "unknown")) {
  detected_commit <- tryCatch(
    system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(detected_commit) == 1L && nzchar(detected_commit)) {
    source_commit <- detected_commit
  } else {
    source_commit <- "unknown"
  }
}
package_version <- as.character(utils::packageVersion("fastphylosig"))

clone <- function(x) unserialize(serialize(x, connection = NULL, version = 2L))

# Delimiter-containing labels are legal under the current package contract.
# The two non-root clades deliberately have the same frozen collapsed key.
tip_labels <- c("a", "b\rc", "a\rb", "c")
base_tree <- structure(
  list(
    edge = matrix(
      c(5L, 6L, 6L, 1L, 6L, 2L, 5L, 7L, 7L, 3L, 7L, 4L),
      ncol = 2L, byrow = TRUE
    ),
    edge.length = c(1.10, 0.70, 0.90, 1.20, 0.80, 1.05),
    tip.label = tip_labels,
    Nnode = 3L
  ),
  class = "phylo"
)

renumber_internal <- function(tree) {
  out <- clone(tree)
  n_tip <- length(out$tip.label)
  # Keep the structural root (n_tip + 1) fixed and exchange the two clades.
  mapping <- c("6" = 7L, "7" = 6L)
  hit <- out$edge %in% c(6L, 7L)
  out$edge[hit] <- unname(mapping[as.character(out$edge[hit])])
  out
}

shuffle_edge_rows <- function(tree) {
  out <- clone(tree)
  rows <- c(3L, 1L, 6L, 4L, 2L, 5L)
  out$edge <- out$edge[rows, , drop = FALSE]
  out$edge.length <- out$edge.length[rows]
  out
}

trees <- list(
  base = base_tree,
  internal_renumbered = renumber_internal(base_tree),
  edge_rows_shuffled = shuffle_edge_rows(base_tree),
  renumbered_and_shuffled = shuffle_edge_rows(renumber_internal(base_tree))
)

# Fixed stochastic inputs are expressed in each tree's tip order.  The
# delimiter labels are in the same order for every representation.
permutations <- rbind(
  c(1L, 2L, 3L, 4L),
  c(2L, 1L, 3L, 4L),
  c(1L, 2L, 4L, 3L),
  c(4L, 3L, 2L, 1L)
)
brownian_states <- matrix(
  c(
    0.10, 0.20, 0.35, 0.40,
    0.25, 0.55, 0.15, 0.75,
    0.60, 0.30, 0.80, 0.45,
    0.90, 0.65, 0.50, 0.05
  ), nrow = 4L, ncol = 4L
)

trait_cont <- stats::setNames(c(0.55, 1.20, 2.10, 3.70), tip_labels)
trait_binary <- stats::setNames(c(0, 0, 1, 1), tip_labels)
trait_delta <- stats::setNames(
  factor(c("A", "A", "B", "B"), levels = c("A", "B")), tip_labels
)

capture <- function(fun) {
  warnings <- list()
  err <- NULL
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- list(
          message = conditionMessage(w), class = class(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- list(message = conditionMessage(e), class = class(e))
      NULL
    }
  )
  list(value = value, warnings = warnings, error = err)
}

warning_text <- function(x) {
  if (!length(x)) return("")
  paste(vapply(x, function(z) z$message, character(1L)), collapse = " || ")
}

warning_classes <- function(x) {
  if (!length(x)) return("")
  paste(vapply(x, function(z) paste(z$class, collapse = "/"), character(1L)),
        collapse = " || ")
}

error_text <- function(x) {
  if (is.null(x)) "" else x$message
}

error_classes <- function(x) {
  if (is.null(x)) "" else paste(x$class, collapse = "/")
}

escape_text <- function(x) {
  if (length(x) != 1L || is.na(x)) return("NA")
  value <- as.character(x)
  value <- gsub("\r", "\\r", value, fixed = TRUE)
  value <- gsub("\n", "\\n", value, fixed = TRUE)
  value
}

safe_equal <- function(x, y) isTRUE(identical(x, y))

all_equal_zero <- function(x, y) {
  if (length(x) != length(y)) return(FALSE)
  if (!length(x)) return(TRUE)
  if (!is.numeric(x) || !is.numeric(y)) return(safe_equal(x, y))
  all(is.finite(x) == is.finite(y)) &&
    all(!is.finite(x) | !is.finite(y) | abs(x - y) == 0)
}

relative_difference <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      !length(x)) return(NA_real_)
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  max(abs(x[ok] - y[ok]) / pmax(1, abs(y[ok])))
}

numeric_difference <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      !length(x)) return(NA_real_)
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  max(abs(x[ok] - y[ok]))
}

scalar_metric <- function(value, name) {
  if (is.null(value)) return(NA_character_)
  if (length(value) == 1L && (is.atomic(value) || is.null(dim(value)))) {
    if (is.logical(value)) return(as.character(value))
    if (is.numeric(value)) return(if (is.finite(value)) format(value, digits = 17) else as.character(value))
    return(escape_text(value))
  }
  if (is.numeric(value)) {
    return(paste(format(as.numeric(value), digits = 17), collapse = ","))
  }
  if (is.logical(value)) return(paste(as.logical(value), collapse = ","))
  if (is.character(value)) return(paste(vapply(value, escape_text, character(1L)), collapse = ","))
  paste0("<", name, ":", class(value)[[1L]], ">")
}

object_field <- function(object, name) {
  if (is.null(object) || !is.list(object) || !name %in% names(object)) NULL else object[[name]]
}

result_metric <- function(result, method, metric) {
  if (is.null(result)) return(NULL)
  if (metric == "matched_species" || metric == "retained_species") {
    metadata <- attr(result, "analysis_metadata", exact = TRUE)
    matching <- if (is.list(metadata)) metadata$matching else NULL
    if (is.list(matching)) {
      value <- matching[[metric]]
      if (!is.null(value)) return(value)
    }
    report <- attr(result, "match_report", exact = TRUE)
    if (is.list(report) && !is.null(report[[metric]])) return(report[[metric]])
    return(NULL)
  }
  if (method == "K") {
    if (metric == "estimate" && is.numeric(result) && is.null(dim(result))) return(result)
    if (metric == "estimate") return(object_field(result, "K"))
    if (metric == "P") return(object_field(result, "P"))
    if (metric == "MCSE") return(object_field(result, "MCSE_P"))
    if (metric == "nsim_successful") return(object_field(result, "nsim_successful"))
    if (metric == "status") return(attr(result, "status", exact = TRUE) %||% object_field(result, "status"))
  }
  if (method == "lambda") {
    if (metric == "estimate") return(object_field(result, "lambda"))
    if (metric == "logLik") return(object_field(result, "logL"))
    if (metric == "LR") return(object_field(result, "LR"))
    if (metric == "P") return(object_field(result, "P"))
    if (metric == "status") return(object_field(result, "status"))
  }
  if (method == "D") {
    if (metric == "estimate") return(object_field(result, "DEstimate"))
    if (metric == "P_random") return(object_field(result, "Pval1"))
    if (metric == "P_Brownian") return(object_field(result, "Pval0"))
    if (metric == "MCSE_random") return(object_field(result, "MCSE_P_random"))
    if (metric == "MCSE_Brownian") return(object_field(result, "MCSE_P_Brownian"))
    if (metric == "status") return(object_field(result, "status"))
    if (metric == "n_successful_random") return(object_field(result, "nsim_successful_random"))
    if (metric == "n_successful_Brownian") return(object_field(result, "nsim_successful_brownian"))
  }
  if (method == "Delta") {
    if (metric == "estimate") return(object_field(result, "delta"))
    if (metric == "P") return(object_field(result, "P"))
    if (metric == "MCSE") return(object_field(result, "MCSE_Delta"))
    if (metric == "ESS_alpha") return(object_field(result, "ESS_alpha"))
    if (metric == "ESS_beta") return(object_field(result, "ESS_beta"))
    if (metric == "Rhat_alpha") return(object_field(result, "split_Rhat_alpha"))
    if (metric == "Rhat_beta") return(object_field(result, "split_Rhat_beta"))
    if (metric == "status") return(object_field(result, "status"))
    if (metric == "n_successful") return(object_field(result, "successful_simulations"))
  }
  if (method == "ACE") {
    if (metric == "logLik") return(object_field(result, "loglik"))
    if (metric == "lik_anc") return(object_field(result, "lik.anc"))
    if (metric == "ace") return(object_field(result, "ace"))
    if (metric == "rates") return(object_field(result, "rates"))
    if (metric == "se") return(object_field(result, "se"))
    if (metric == "convergence") return(object_field(result, "convergence"))
    if (metric == "message") return(object_field(result, "message"))
  }
  NULL
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ACE reports one row per internal node.  Internal IDs are representation
# details, so align ancestral likelihood rows by a length-prefixed encoding of
# the sorted descendant tip labels.  Length prefixes avoid reintroducing the
# delimiter collision that motivated this audit.
descendant_sets_by_id <- function(tree) {
  n_tip <- length(tree$tip.label)
  edge <- as.matrix(tree$edge)
  children <- split(edge[, 2L], edge[, 1L])
  memo <- new.env(hash = TRUE, parent = emptyenv())
  descendants <- function(node) {
    key <- as.character(node)
    if (exists(key, memo, inherits = FALSE)) return(get(key, memo, inherits = FALSE))
    out <- if (node <= n_tip) tree$tip.label[[node]] else {
      kids <- children[[key]]
      sort(unlist(lapply(kids, descendants), use.names = FALSE))
    }
    assign(key, out, memo)
    out
  }
  ids <- sort(unique(c(edge[, 1L], edge[, 2L])))
  ids <- ids[ids > n_tip]
  stats::setNames(lapply(ids, descendants), as.character(ids))
}

clade_key <- function(labels) {
  labels <- as.character(labels)
  paste(sprintf("%d:%s", nchar(labels, type = "chars"), labels), collapse = ";")
}

aligned_ace_lik_anc <- function(result, tree) {
  lik <- if (is.list(result)) result$lik.anc else NULL
  if (is.null(lik) || is.null(dim(lik)) || !is.matrix(lik)) return(lik)
  by_id <- descendant_sets_by_id(tree)
  row_ids <- suppressWarnings(as.integer(rownames(lik)))
  if (length(row_ids) != nrow(lik) || anyNA(row_ids) ||
      !all(as.character(row_ids) %in% names(by_id))) {
    row_ids <- seq.int(nrow(lik)) + length(tree$tip.label)
  }
  keys <- vapply(row_ids, function(id) {
    if (!as.character(id) %in% names(by_id)) return(sprintf("<row:%d>", id))
    clade_key(by_id[[as.character(id)]])
  }, character(1L))
  if (anyDuplicated(keys)) {
    # Duplicate clade keys are not expected in a valid tree; retain a stable
    # row-index suffix so this audit never silently drops a row.
    keys <- paste0(keys, sprintf("#%08d", seq_along(keys)))
  }
  ord <- order(keys)
  out <- lik[ord, , drop = FALSE]
  rownames(out) <- keys[ord]
  out
}

ace_metric_value <- function(result, metric, tree) {
  if (identical(metric, "lik_anc")) return(aligned_ace_lik_anc(result, tree))
  result_metric(result, "ACE", metric)
}

comparison_rows <- list()
add_comparison <- function(...) comparison_rows[[length(comparison_rows) + 1L]] <<- data.frame(...)

# Biology-level equivalence is assessed independently of package
# canonicalization, using labeled patristic distances and the rooted edge
# graph after node IDs are abstracted by each node's descendant tip set.
descendant_sets <- function(tree) {
  n_tip <- length(tree$tip.label)
  edge <- as.matrix(tree$edge)
  parent <- split(edge[, 2L], edge[, 1L])
  memo <- new.env(hash = TRUE, parent = emptyenv())
  rec <- function(node) {
    key <- as.character(node)
    if (exists(key, memo, inherits = FALSE)) return(get(key, memo, inherits = FALSE))
    if (node <= n_tip) out <- tree$tip.label[[node]] else {
      kids <- parent[[key]]
      out <- sort(unlist(lapply(kids, rec), use.names = FALSE))
    }
    assign(key, out, memo)
    out
  }
  ids <- sort(unique(edge[, 1L]))
  internals <- ids[ids > n_tip]
  sets <- lapply(internals, rec)
  names(sets) <- vapply(sets, function(z) paste(z, collapse = "\r"), character(1L))
  sets
}

patristic <- function(tree) {
  as.matrix(ape::cophenetic.phylo(tree)[tip_labels, tip_labels, drop = FALSE])
}

tree_equivalence <- data.frame()
for (nm in names(trees)) {
  tr <- trees[[nm]]
  pd <- tryCatch(patristic(tr), error = function(e) NULL)
  tree_equivalence <- rbind(
    tree_equivalence,
    data.frame(
      fixture = nm,
      rooted = tryCatch(isTRUE(ape::is.rooted(tr)), error = function(e) NA),
      patristic_vs_base = if (nm == "base") TRUE else {
        p0 <- patristic(base_tree)
        !is.null(pd) && identical(dim(pd), dim(p0)) && all(pd == p0)
      },
      edge_rows = nrow(tr$edge),
      tip_label_count = length(tr$tip.label),
      Nnode = tr$Nnode,
      structural_descendant_sets = length(descendant_sets(tr)) == length(descendant_sets(base_tree)),
      stringsAsFactors = FALSE
    )
  )
}
write.csv(tree_equivalence, file.path(out_dir, "fixture_equivalence.csv"), row.names = FALSE, na = "")

prepared <- list()
for (nm in names(trees)) {
  prepared[[nm]] <- capture(function() prepare_tree(trees[[nm]]))
  c <- prepared[[nm]]
  if (!is.null(c$error)) {
    add_comparison(
      scope = "prepare_tree", fixture = nm, route = "raw", method = "all",
      metric = "error", value = error_text(c$error), warning = warning_text(c$warnings),
      warning_class = warning_classes(c$warnings), error_class = error_classes(c$error),
      stringsAsFactors = FALSE
    )
  }
}

prepared_fingerprint <- data.frame()
for (nm in names(trees)) {
  c <- prepared[[nm]]
  ctx <- c$value
  if (is.null(ctx) || !is.list(ctx)) next
  cache_keys <- if (is.environment(ctx$structural_cache)) ls(ctx$structural_cache, all.names = TRUE) else character()
  prepared_fingerprint <- rbind(
    prepared_fingerprint,
    data.frame(
      fixture = nm,
      route = "prepared_context",
      fingerprint = escape_text(ctx$fingerprint %||% NA_character_),
      cache_key_count = length(cache_keys),
      cache_keys = paste(cache_keys, collapse = ";"),
      inspection_present = is.list(ctx$inspection),
      canonical_mapping_present = is.list(ctx$canonical_mapping),
      error = "",
      stringsAsFactors = FALSE
    )
  )
}
write.csv(prepared_fingerprint, file.path(out_dir, "prepared_fingerprint_cache.csv"), row.names = FALSE, na = "")

# Public check_tree and prepare_tree propagation.  Check results are compared
# by complete object identity where possible; prepared contexts are compared by
# their contract-bearing fields because environments/compiled pointers are not
# identity-comparable across calls.
check_rows <- data.frame()
for (nm in names(trees)) {
  got <- capture(function() check_tree(trees[[nm]], signal = c("K", "lambda", "D", "Delta")))
  value <- got$value
  check_rows <- rbind(
    check_rows,
    data.frame(
      fixture = nm, status = if (is.null(got$error)) "ok" else "error",
      ready = if (is.list(value)) scalar_metric(value$ready, "ready") else "NA",
      ready_by_signal = if (is.list(value)) scalar_metric(value$ready_by_signal, "ready_by_signal") else "NA",
      tree_summary = if (is.list(value)) scalar_metric(value$tree_summary, "tree_summary") else "NA",
      issues = if (is.list(value)) scalar_metric(value$issues, "issues") else "NA",
      warning = warning_text(got$warnings), warning_class = warning_classes(got$warnings),
      error = error_text(got$error), error_class = error_classes(got$error),
      stringsAsFactors = FALSE
    )
  )
}
write.csv(check_rows, file.path(out_dir, "check_tree_results.csv"), row.names = FALSE, na = "")

check_base <- capture(function() check_tree(base_tree, signal = c("K", "lambda", "D", "Delta")))
for (nm in names(trees)[-1L]) {
  got <- capture(function() check_tree(trees[[nm]], signal = c("K", "lambda", "D", "Delta")))
  add_comparison(
    scope = "check_tree", fixture = nm, route = "raw", method = "all",
    metric = "complete_object_identical_to_base",
    value = safe_equal(got$value, check_base$value),
    warning = warning_text(got$warnings), warning_class = warning_classes(got$warnings),
    error_class = error_classes(got$error), error = error_text(got$error),
    stringsAsFactors = FALSE
  )
}

continuous_specs <- list(
  K_test_false = list(method = "K", test = FALSE),
  K_test_true_fixed = list(method = "K", test = TRUE, nsim = nrow(permutations), permutations = permutations),
  lambda_test_true = list(method = "lambda", test = TRUE, lambda_profile = FALSE)
)

run_continuous <- function(tree, ctx, spec, prepared_route) {
  args <- list(tree = if (prepared_route) ctx else tree, x = trait_cont,
               verbose = FALSE, progress = FALSE)
  args <- c(args, spec)
  args$method <- NULL
  if (identical(spec$method, "K")) {
    return(capture(function() do.call(fast_k, args)))
  }
  capture(function() do.call(fast_lambda, args))
}

continuous_rows <- data.frame()
continuous_results <- list()
for (nm in names(trees)) {
  ctx <- if (is.list(prepared[[nm]])) prepared[[nm]]$value else NULL
  for (spec_name in names(continuous_specs)) {
    spec <- continuous_specs[[spec_name]]
    for (route in c("raw", "prepared")) {
      set.seed(20260906)
      got <- run_continuous(trees[[nm]], ctx, spec, route == "prepared")
      key <- paste(nm, spec_name, route, sep = "|")
      continuous_results[[key]] <- got
      for (metric in if (spec$method == "K" && isTRUE(spec$test)) {
        c("estimate", "P", "MCSE", "nsim_successful", "matched_species", "retained_species", "status")
      } else if (spec$method == "K") {
        c("estimate", "matched_species", "retained_species", "status")
      } else {
        c("estimate", "logLik", "LR", "P", "matched_species", "retained_species", "status")
      }) {
        value <- result_metric(got$value, spec$method, metric)
        continuous_rows <- rbind(
          continuous_rows,
          data.frame(
            fixture = nm, route = route, method = spec$method,
            workload = spec_name, metric = metric,
            value = scalar_metric(value, metric),
            value_type = if (is.null(value)) "NULL" else paste(class(value), collapse = "/"),
            warning = warning_text(got$warnings), warning_class = warning_classes(got$warnings),
            error = error_text(got$error), error_class = error_classes(got$error),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}
write.csv(continuous_rows, file.path(out_dir, "continuous_results.csv"), row.names = FALSE, na = "")

for (spec_name in names(continuous_specs)) {
  spec <- continuous_specs[[spec_name]]
  base_key <- paste("base", spec_name, "raw", sep = "|")
  for (nm in names(trees)[-1L]) {
    for (route in c("raw", "prepared")) {
      key <- paste(nm, spec_name, route, sep = "|")
      lhs <- continuous_results[[key]]
      rhs <- continuous_results[[base_key]]
      for (metric in if (spec$method == "K" && isTRUE(spec$test)) c("estimate", "P", "MCSE", "nsim_successful", "matched_species", "retained_species", "status") else if (spec$method == "K") c("estimate", "matched_species", "retained_species", "status") else c("estimate", "logLik", "LR", "P", "matched_species", "retained_species", "status")) {
        x <- result_metric(lhs$value, spec$method, metric)
        y <- result_metric(rhs$value, spec$method, metric)
        equal <- if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y)
        add_comparison(
          scope = "estimator_vs_base", fixture = nm, route = route,
          method = spec$method, metric = paste(spec_name, metric, sep = ":"),
          value = equal, absolute_difference = numeric_difference(x, y),
          relative_difference = relative_difference(x, y),
          warning = warning_text(lhs$warnings), warning_class = warning_classes(lhs$warnings),
          error = error_text(lhs$error), error_class = error_classes(lhs$error),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

run_d <- function(tree, ctx, route, test_value) {
  args <- list(
    tree = if (route == "prepared") ctx else tree, x = trait_binary,
    test = test_value, nsim = nrow(permutations), permutations = permutations,
    brownian_states = brownian_states, return_sim = FALSE,
    verbose = FALSE, progress = FALSE, ncores = 1L
  )
  set.seed(20260906)
  capture(function() do.call(fast_d, args))
}

d_results <- list()
d_rows <- data.frame()
for (nm in names(trees)) {
  ctx <- if (is.list(prepared[[nm]])) prepared[[nm]]$value else NULL
  for (test_value in c(FALSE, TRUE)) {
    for (route in c("raw", "prepared")) {
      got <- run_d(trees[[nm]], ctx, route, test_value)
      key <- paste(nm, test_value, route, sep = "|")
      d_results[[key]] <- got
      metrics <- if (test_value) c("estimate", "P_random", "P_Brownian", "MCSE_random", "MCSE_Brownian", "matched_species", "retained_species", "status", "n_successful_random", "n_successful_Brownian") else c("estimate", "matched_species", "retained_species", "status")
      for (metric in metrics) {
        value <- result_metric(got$value, "D", metric)
        d_rows <- rbind(
          d_rows,
          data.frame(
            fixture = nm, route = route, test = test_value, method = "D",
            metric = metric, value = scalar_metric(value, metric),
            warning = warning_text(got$warnings), warning_class = warning_classes(got$warnings),
            error = error_text(got$error), error_class = error_classes(got$error),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}
write.csv(d_rows, file.path(out_dir, "d_results.csv"), row.names = FALSE, na = "")
for (test_value in c(FALSE, TRUE)) {
  for (nm in names(trees)[-1L]) {
    for (route in c("raw", "prepared")) {
      lhs <- d_results[[paste(nm, test_value, route, sep = "|")]]
      rhs <- d_results[[paste("base", test_value, "raw", sep = "|")]]
      metrics <- if (test_value) c("estimate", "P_random", "P_Brownian", "MCSE_random", "MCSE_Brownian", "matched_species", "retained_species", "status", "n_successful_random", "n_successful_Brownian") else c("estimate", "matched_species", "retained_species", "status")
      for (metric in metrics) {
        x <- result_metric(lhs$value, "D", metric)
        y <- result_metric(rhs$value, "D", metric)
        add_comparison(
          scope = "estimator_vs_base", fixture = nm, route = route,
          method = "D", metric = paste0("test=", test_value, ":", metric),
          value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
          absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
          warning = warning_text(lhs$warnings), warning_class = warning_classes(lhs$warnings),
          error = error_text(lhs$error), error_class = error_classes(lhs$error),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

run_delta <- function(tree, ctx, route, test_value) {
  args <- list(
    tree = if (route == "prepared") ctx else tree, x = trait_delta,
    test = test_value, nsim = nrow(permutations), permutations = permutations,
    mcmc_sim = 100L, thin = 2L, burn = 10L, lambda0 = 0.1,
    proposal_sd = 0.5, model = "ER", entropy = "LSE", return_sim = FALSE,
    verbose = FALSE, progress = FALSE, ncores = 1L
  )
  set.seed(20260906)
  capture(function() do.call(fast_delta, args))
}

delta_results <- list()
delta_rows <- data.frame()
for (nm in names(trees)) {
  ctx <- if (is.list(prepared[[nm]])) prepared[[nm]]$value else NULL
  for (test_value in c(FALSE, TRUE)) {
    for (route in c("raw", "prepared")) {
      got <- run_delta(trees[[nm]], ctx, route, test_value)
      key <- paste(nm, test_value, route, sep = "|")
      delta_results[[key]] <- got
      metrics <- if (test_value) c("estimate", "P", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status", "n_successful") else c("estimate", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status")
      for (metric in metrics) {
        value <- result_metric(got$value, "Delta", metric)
        delta_rows <- rbind(
          delta_rows,
          data.frame(
            fixture = nm, route = route, test = test_value, method = "Delta",
            metric = metric, value = scalar_metric(value, metric),
            warning = warning_text(got$warnings), warning_class = warning_classes(got$warnings),
            error = error_text(got$error), error_class = error_classes(got$error),
            rng_kind = paste(RNGkind(), collapse = "|"),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}
write.csv(delta_rows, file.path(out_dir, "delta_results.csv"), row.names = FALSE, na = "")
for (test_value in c(FALSE, TRUE)) {
  for (nm in names(trees)[-1L]) {
    for (route in c("raw", "prepared")) {
      lhs <- delta_results[[paste(nm, test_value, route, sep = "|")]]
      rhs <- delta_results[[paste("base", test_value, "raw", sep = "|")]]
      metrics <- if (test_value) c("estimate", "P", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status", "n_successful") else c("estimate", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status")
      for (metric in metrics) {
        x <- result_metric(lhs$value, "Delta", metric)
        y <- result_metric(rhs$value, "Delta", metric)
        add_comparison(
          scope = "estimator_vs_base", fixture = nm, route = route,
          method = "Delta", metric = paste0("test=", test_value, ":", metric),
          value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
          absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
          warning = warning_text(lhs$warnings), warning_class = warning_classes(lhs$warnings),
          error = error_text(lhs$error), error_class = error_classes(lhs$error),
          rng_same_seed = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

run_ace <- function(tree, ctx, route) {
  args <- list(
    x = trait_delta, phy = if (route == "prepared") NULL else tree,
    prepared = if (route == "prepared") ctx else NULL,
    type = "discrete", method = "ML", model = "ER", CI = TRUE,
    marginal = FALSE, progress = FALSE
  )
  capture(function() do.call(fast_ace, args))
}

ace_results <- list()
ace_rows <- data.frame()
for (nm in names(trees)) {
  ctx <- if (is.list(prepared[[nm]])) prepared[[nm]]$value else NULL
  for (route in c("raw", "prepared")) {
    got <- run_ace(trees[[nm]], ctx, route)
    key <- paste(nm, route, sep = "|")
    ace_results[[key]] <- got
    actual_tree <- if (route == "prepared" && is.list(ctx)) ctx$tree else trees[[nm]]
    for (metric in c("logLik", "lik_anc", "ace", "rates", "se", "convergence", "message")) {
      value <- ace_metric_value(got$value, metric, actual_tree)
      ace_rows <- rbind(
        ace_rows,
        data.frame(
          fixture = nm, route = route, method = "ACE", metric = metric,
          value = scalar_metric(value, metric),
          warning = warning_text(got$warnings), warning_class = warning_classes(got$warnings),
          error = error_text(got$error), error_class = error_classes(got$error),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}
write.csv(ace_rows, file.path(out_dir, "ace_results.csv"), row.names = FALSE, na = "")
for (nm in names(trees)[-1L]) {
  for (route in c("raw", "prepared")) {
    lhs <- ace_results[[paste(nm, route, sep = "|")]]
    rhs <- ace_results[["base|raw"]]
    lhs_tree <- if (route == "prepared" && is.list(prepared[[nm]]$value)) {
      prepared[[nm]]$value$tree
    } else trees[[nm]]
    for (metric in c("logLik", "lik_anc", "ace", "rates", "se", "convergence", "message")) {
      x <- ace_metric_value(lhs$value, metric, lhs_tree)
      y <- ace_metric_value(rhs$value, metric, base_tree)
      add_comparison(
        scope = "estimator_vs_base", fixture = nm, route = route,
        method = "ACE", metric = metric,
        value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
        absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
        warning = warning_text(lhs$warnings), warning_class = warning_classes(lhs$warnings),
        error = error_text(lhs$error), error_class = error_classes(lhs$error),
        stringsAsFactors = FALSE
      )
    }
  }
}

# Compare raw versus prepared within every representation.  Public prepared
# calls deliberately validate their context once; this table records whether
# their contract-bearing outputs match the raw route.
for (nm in names(trees)) {
  for (spec_name in names(continuous_specs)) {
    spec <- continuous_specs[[spec_name]]
    raw <- continuous_results[[paste(nm, spec_name, "raw", sep = "|")]]
    prep <- continuous_results[[paste(nm, spec_name, "prepared", sep = "|")]]
    metrics <- if (spec$method == "K" && isTRUE(spec$test)) c("estimate", "P", "MCSE", "nsim_successful", "matched_species", "retained_species", "status") else if (spec$method == "K") c("estimate", "matched_species", "retained_species", "status") else c("estimate", "logLik", "LR", "P", "matched_species", "retained_species", "status")
    for (metric in metrics) {
      x <- result_metric(raw$value, spec$method, metric)
      y <- result_metric(prep$value, spec$method, metric)
      add_comparison(
        scope = "raw_vs_prepared", fixture = nm, route = "within_fixture",
        method = spec$method, metric = paste(spec_name, metric, sep = ":"),
        value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
        absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
        warning = warning_text(prep$warnings), warning_class = warning_classes(prep$warnings),
        error = error_text(prep$error), error_class = error_classes(prep$error),
        stringsAsFactors = FALSE
      )
    }
  }
  for (test_value in c(FALSE, TRUE)) {
    raw <- d_results[[paste(nm, test_value, "raw", sep = "|")]]
    prep <- d_results[[paste(nm, test_value, "prepared", sep = "|")]]
    metrics <- if (test_value) c("estimate", "P_random", "P_Brownian", "MCSE_random", "MCSE_Brownian", "matched_species", "retained_species", "status", "n_successful_random", "n_successful_Brownian") else c("estimate", "matched_species", "retained_species", "status")
    for (metric in metrics) {
      x <- result_metric(raw$value, "D", metric)
      y <- result_metric(prep$value, "D", metric)
      add_comparison(
        scope = "raw_vs_prepared", fixture = nm, route = "within_fixture",
        method = "D", metric = paste0("test=", test_value, ":", metric),
        value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
        absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
        warning = warning_text(prep$warnings), warning_class = warning_classes(prep$warnings),
        error = error_text(prep$error), error_class = error_classes(prep$error),
        stringsAsFactors = FALSE
      )
    }
  }
  for (test_value in c(FALSE, TRUE)) {
    raw <- delta_results[[paste(nm, test_value, "raw", sep = "|")]]
    prep <- delta_results[[paste(nm, test_value, "prepared", sep = "|")]]
    metrics <- if (test_value) c("estimate", "P", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status", "n_successful") else c("estimate", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status")
    for (metric in metrics) {
      x <- result_metric(raw$value, "Delta", metric)
      y <- result_metric(prep$value, "Delta", metric)
      add_comparison(
        scope = "raw_vs_prepared", fixture = nm, route = "within_fixture",
        method = "Delta", metric = paste0("test=", test_value, ":", metric),
        value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
        absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
        warning = warning_text(prep$warnings), warning_class = warning_classes(prep$warnings),
        error = error_text(prep$error), error_class = error_classes(prep$error),
        rng_same_seed = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  raw <- ace_results[[paste(nm, "raw", sep = "|")]]
  prep <- ace_results[[paste(nm, "prepared", sep = "|")]]
  raw_tree <- trees[[nm]]
  prep_tree <- if (is.list(prepared[[nm]]$value)) prepared[[nm]]$value$tree else raw_tree
  for (metric in c("logLik", "lik_anc", "ace", "rates", "se", "convergence", "message")) {
    x <- ace_metric_value(raw$value, metric, raw_tree)
    y <- ace_metric_value(prep$value, metric, prep_tree)
    add_comparison(
      scope = "raw_vs_prepared", fixture = nm, route = "within_fixture",
      method = "ACE", metric = metric,
      value = if (is.numeric(x) && is.numeric(y)) all_equal_zero(x, y) else safe_equal(x, y),
      absolute_difference = numeric_difference(x, y), relative_difference = relative_difference(x, y),
      warning = warning_text(prep$warnings), warning_class = warning_classes(prep$warnings),
      error = error_text(prep$error), error_class = error_classes(prep$error),
      stringsAsFactors = FALSE
    )
  }
}

# Warning classes and failure behavior are separate contract dimensions from
# numerical values.  Record them explicitly for every representation and for
# raw/prepared parity; this prevents an equal estimate from masking a changed
# public failure or warning path.
add_condition_rows <- function(scope, fixture, route, method, workload,
                               lhs, rhs, rng_same_seed = NA) {
  add_comparison(
    scope = scope, fixture = fixture, route = route, method = method,
    metric = paste(workload, "warning_contract", sep = ":"),
    value = safe_equal(lhs$warnings, rhs$warnings),
    warning = warning_text(lhs$warnings), warning_class = warning_classes(lhs$warnings),
    error = error_text(lhs$error), error_class = error_classes(lhs$error),
    rng_same_seed = rng_same_seed, stringsAsFactors = FALSE
  )
  add_comparison(
    scope = scope, fixture = fixture, route = route, method = method,
    metric = paste(workload, "error_contract", sep = ":"),
    value = safe_equal(lhs$error, rhs$error),
    warning = warning_text(lhs$warnings), warning_class = warning_classes(lhs$warnings),
    error = error_text(lhs$error), error_class = error_classes(lhs$error),
    rng_same_seed = rng_same_seed, stringsAsFactors = FALSE
  )
}

for (spec_name in names(continuous_specs)) {
  spec <- continuous_specs[[spec_name]]
  for (nm in names(trees)[-1L]) {
    for (route in c("raw", "prepared")) {
      lhs <- continuous_results[[paste(nm, spec_name, route, sep = "|")]]
      rhs <- continuous_results[[paste("base", spec_name, "raw", sep = "|")]]
      add_condition_rows("estimator_vs_base", nm, route, spec$method,
                         spec_name, lhs, rhs, rng_same_seed = TRUE)
    }
  }
  for (nm in names(trees)) {
    add_condition_rows(
      "raw_vs_prepared", nm, "within_fixture", spec$method, spec_name,
      continuous_results[[paste(nm, spec_name, "raw", sep = "|")]],
      continuous_results[[paste(nm, spec_name, "prepared", sep = "|")]],
      rng_same_seed = TRUE
    )
  }
}
for (test_value in c(FALSE, TRUE)) {
  workload <- paste0("test=", test_value)
  for (nm in names(trees)[-1L]) {
    for (route in c("raw", "prepared")) {
      add_condition_rows(
        "estimator_vs_base", nm, route, "D", workload,
        d_results[[paste(nm, test_value, route, sep = "|")]],
        d_results[[paste("base", test_value, "raw", sep = "|")]],
        rng_same_seed = TRUE
      )
    }
  }
  for (nm in names(trees)) {
    add_condition_rows(
      "raw_vs_prepared", nm, "within_fixture", "D", workload,
      d_results[[paste(nm, test_value, "raw", sep = "|")]],
      d_results[[paste(nm, test_value, "prepared", sep = "|")]],
      rng_same_seed = TRUE
    )
  }
  for (nm in names(trees)[-1L]) {
    for (route in c("raw", "prepared")) {
      add_condition_rows(
        "estimator_vs_base", nm, route, "Delta", workload,
        delta_results[[paste(nm, test_value, route, sep = "|")]],
        delta_results[[paste("base", test_value, "raw", sep = "|")]],
        rng_same_seed = TRUE
      )
    }
  }
  for (nm in names(trees)) {
    add_condition_rows(
      "raw_vs_prepared", nm, "within_fixture", "Delta", workload,
      delta_results[[paste(nm, test_value, "raw", sep = "|")]],
      delta_results[[paste(nm, test_value, "prepared", sep = "|")]],
      rng_same_seed = TRUE
    )
  }
}
for (nm in names(trees)[-1L]) {
  for (route in c("raw", "prepared")) {
    add_condition_rows(
      "estimator_vs_base", nm, route, "ACE", "CI=TRUE", 
      ace_results[[paste(nm, route, sep = "|")]],
      ace_results[["base|raw"]], rng_same_seed = NA
    )
  }
}
for (nm in names(trees)) {
  add_condition_rows(
    "raw_vs_prepared", nm, "within_fixture", "ACE", "CI=TRUE",
    ace_results[[paste(nm, "raw", sep = "|")]],
    ace_results[[paste(nm, "prepared", sep = "|")]], rng_same_seed = NA
  )
}

# Cross-locale propagation audit.  This is intentionally a small contract
# audit, not a performance run.  The same legal delimiter-containing tree,
# traits, permutations, Brownian states, seed, and single-thread Delta path
# are evaluated under C and under the requested Windows English locale.
locale_a <- "C"
locale_b_requested <- "English_United States.1252"
locale_b_actual <- tryCatch(
  suppressWarnings(Sys.setlocale("LC_COLLATE", locale_b_requested)),
  error = function(e) NA_character_
)
locale_b_usable <- is.character(locale_b_actual) && length(locale_b_actual) == 1L &&
  !is.na(locale_b_actual) && nzchar(locale_b_actual) &&
  grepl("English_United States", locale_b_actual, fixed = TRUE)
locale_b_setup <- if (locale_b_usable) "RUN" else "NOT_RUN_ENVIRONMENT"

# Reuse the exact 12-tip locale-sensitive fixture from the canonicalization
# audit. Its six non-root clades combine case, prefixes, accents, punctuation,
# numeric-looking labels, and a separator-containing label. The fixed inputs
# are resized only after the representation-comparison audit above completes.
locale_labels <- c(
  "A", "a0", "a1", "aa", "e\u00e9", "e\u00e8",
  "a-1", "a_1", "2", "10", "a\rb", "c"
)
canonicalize_core <- getFromNamespace(".safe_canonicalize_core", "fastphylosig")
fingerprint_core <- getFromNamespace(".tree_fingerprint", "fastphylosig")
locale_tree <- NULL
locale_fixture_seed <- NA_integer_
locale_canonical_c <- NULL
locale_canonical_b <- NULL
if (locale_b_usable) {
  for (candidate_seed in 20260906L + 0:127L) {
    set.seed(candidate_seed)
    candidate <- ape::rtree(length(locale_labels))
    candidate$tip.label <- locale_labels
    suppressWarnings(Sys.setlocale("LC_COLLATE", "C"))
    canonical_c <- canonicalize_core(candidate)
    suppressWarnings(Sys.setlocale("LC_COLLATE", locale_b_actual))
    canonical_b <- canonicalize_core(candidate)
    if (!identical(fingerprint_core(canonical_c), fingerprint_core(canonical_b))) {
      locale_tree <- candidate
      locale_fixture_seed <- candidate_seed
      locale_canonical_c <- canonical_c
      locale_canonical_b <- canonical_b
      break
    }
  }
} else {
  locale_fixture_seed <- 20260906L
  set.seed(locale_fixture_seed)
  locale_tree <- ape::rtree(length(locale_labels))
  locale_tree$tip.label <- locale_labels
}
suppressWarnings(Sys.setlocale("LC_COLLATE", "C"))
if (locale_b_usable && is.null(locale_tree)) {
  stop("no deterministic locale-sensitive fully binary fixture was found",
       call. = FALSE)
}
write.csv(data.frame(
  fixture_seed = locale_fixture_seed,
  locale_a = "C",
  locale_b = if (locale_b_usable) locale_b_actual else locale_b_requested,
  rooted = isTRUE(ape::is.rooted(locale_tree)),
  fully_binary = isTRUE(ape::is.binary.phylo(locale_tree)),
  canonical_fingerprint_a = if (is.list(locale_canonical_c))
    escape_text(fingerprint_core(locale_canonical_c)) else "",
  canonical_fingerprint_b = if (is.list(locale_canonical_b))
    escape_text(fingerprint_core(locale_canonical_b)) else "",
  canonical_fingerprint_identical = if (is.list(locale_canonical_c) &&
                                        is.list(locale_canonical_b))
    identical(fingerprint_core(locale_canonical_c),
              fingerprint_core(locale_canonical_b)) else NA,
  stringsAsFactors = FALSE
), file.path(out_dir, "locale_fixture.csv"), row.names = FALSE, na = "")
trait_cont <- stats::setNames(seq(0.25, 3.00, length.out = 12L), locale_labels)
trait_binary <- stats::setNames(rep(c(0, 1), 6L), locale_labels)
trait_delta <- stats::setNames(
  factor(rep(c("A", "B", "C"), 4L), levels = c("A", "B", "C")),
  locale_labels
)
permutations <- rbind(
  seq_len(12L),
  c(2:12, 1L),
  rev(seq_len(12L)),
  c(7:12, 1:6)
)
continuous_specs$K_test_true_fixed$nsim <- nrow(permutations)
continuous_specs$K_test_true_fixed$permutations <- permutations
brownian_states <- matrix(
  seq(0.01, 0.99, length.out = 48L), nrow = 12L, ncol = 4L
)

stable_label_code <- function(value) {
  ints <- utf8ToInt(enc2utf8(as.character(value)))
  if (!length(ints)) "<empty>" else paste(sprintf("%08d", ints), collapse = ".")
}

stable_clade_key <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[order(vapply(labels, stable_label_code, character(1L)))]
  paste(sprintf("%d:%s", nchar(labels, type = "chars"), labels), collapse = ";")
}

aligned_ace_lik_anc_stable <- function(result, tree) {
  lik <- if (is.list(result)) result$lik.anc else NULL
  if (is.null(lik) || !is.matrix(lik)) return(lik)
  by_id <- descendant_sets_by_id(tree)
  row_ids <- suppressWarnings(as.integer(rownames(lik)))
  if (length(row_ids) != nrow(lik) || anyNA(row_ids) ||
      !all(as.character(row_ids) %in% names(by_id))) {
    row_ids <- seq.int(nrow(lik)) + length(tree$tip.label)
  }
  keys <- vapply(row_ids, function(id) {
    ids <- as.character(id)
    if (!ids %in% names(by_id)) return(sprintf("<row:%d>", id))
    stable_clade_key(by_id[[ids]])
  }, character(1L))
  if (anyDuplicated(keys)) keys <- paste0(keys, sprintf("#%08d", seq_along(keys)))
  ord <- order(keys)
  out <- lik[ord, , drop = FALSE]
  rownames(out) <- keys[ord]
  out
}

locale_suite <- function(label, requested, setup_status = "RUN") {
  actual <- tryCatch(
    suppressWarnings(Sys.setlocale("LC_COLLATE", requested)),
    error = function(e) NA_character_
  )
  usable <- identical(setup_status, "RUN") && is.character(actual) &&
    length(actual) == 1L && !is.na(actual) && nzchar(actual)
  suite <- list(label = label, requested = requested, actual = actual,
                status = if (usable) "RUN" else "NOT_RUN_ENVIRONMENT")
  if (!usable) return(suite)
  suite$check <- capture(function() check_tree(locale_tree, signal = c("K", "lambda", "D", "Delta")))
  suite$prepare <- capture(function() prepare_tree(locale_tree))
  suite$ctx <- suite$prepare$value
  suite$continuous <- list()
  for (spec_name in names(continuous_specs)) {
    set.seed(20260906)
    suite$continuous[[paste(spec_name, "raw", sep = "|")]] <-
      run_continuous(locale_tree, suite$ctx, continuous_specs[[spec_name]], FALSE)
    set.seed(20260906)
    suite$continuous[[paste(spec_name, "prepared", sep = "|")]] <-
      run_continuous(locale_tree, suite$ctx, continuous_specs[[spec_name]], TRUE)
  }
  suite$d <- list()
  for (test_value in c(FALSE, TRUE)) {
    suite$d[[paste(test_value, "raw", sep = "|")]] <-
      run_d(locale_tree, suite$ctx, "raw", test_value)
    suite$d[[paste(test_value, "prepared", sep = "|")]] <-
      run_d(locale_tree, suite$ctx, "prepared", test_value)
  }
  suite$delta <- list()
  for (test_value in c(FALSE, TRUE)) {
    suite$delta[[paste(test_value, "raw", sep = "|")]] <-
      run_delta(locale_tree, suite$ctx, "raw", test_value)
    suite$delta[[paste(test_value, "prepared", sep = "|")]] <-
      run_delta(locale_tree, suite$ctx, "prepared", test_value)
  }
  suite$ace <- list(
    raw = run_ace(locale_tree, suite$ctx, "raw"),
    prepared = run_ace(locale_tree, suite$ctx, "prepared")
  )
  suite
}

locale_c_suite <- locale_suite("C", "C", setup_status = "RUN")
locale_b_suite <- if (locale_b_usable) {
  locale_suite("English_United States.1252", locale_b_actual, setup_status = "RUN")
} else {
  list(label = locale_b_requested, requested = locale_b_requested,
       actual = locale_b_actual, status = "NOT_RUN_ENVIRONMENT",
       check = NULL, prepare = NULL, ctx = NULL,
       continuous = list(), d = list(), delta = list(), ace = list())
}

locale_rows <- list()
add_locale_row <- function(method, route, workload, metric, value_a, value_b,
                            result_a = NULL, result_b = NULL,
                            ctx_a = NULL, ctx_b = NULL,
                            tree_a = locale_tree, tree_b = locale_tree,
                            rng_same_seed = NA) {
  numeric_values <- is.numeric(value_a) && is.numeric(value_b)
  equal <- if (is.null(value_b) || (length(value_b) == 1L && is.na(value_b))) {
    NA
  } else if (numeric_values) {
    all_equal_zero(value_a, value_b)
  } else {
    safe_equal(value_a, value_b)
  }
  abs_diff <- if (numeric_values) numeric_difference(value_a, value_b) else NA_real_
  rel_diff <- if (numeric_values) relative_difference(value_a, value_b) else NA_real_
  status_value <- function(result) {
    if (is.null(result)) return("NOT_RUN_ENVIRONMENT")
    if (is.list(result) && !is.null(result$error)) return("error")
    if (is.list(result) && !is.null(result$value)) {
      z <- result_metric(result$value, method, "status")
      if (!is.null(z)) return(scalar_metric(z, "status"))
      if (method == "check_tree" && is.list(result$value)) return(scalar_metric(result$value$ready, "ready"))
    }
    "ok"
  }
  warning_value <- function(result) {
    if (is.null(result)) "" else warning_classes(result$warnings)
  }
  error_value <- function(result) {
    if (is.null(result)) "" else error_text(result$error)
  }
  retained_value <- function(result) {
    if (is.null(result) || is.null(result$value)) return(NA_character_)
    scalar_metric(result_metric(result$value, method, "retained_species"), "retained_species")
  }
  fingerprint_value <- function(ctx) {
    if (!is.list(ctx)) "" else escape_text(ctx$fingerprint %||% NA_character_)
  }
  locale_rows[[length(locale_rows) + 1L]] <<- data.frame(
    method = method, route = route, workload = workload, metric = metric,
    locale_fixture_seed = locale_fixture_seed,
    locale_a = locale_c_suite$label, locale_b = locale_b_suite$label,
    locale_a_status = locale_c_suite$status, locale_b_status = locale_b_suite$status,
    value_a = scalar_metric(value_a, metric),
    value_b = if (is.null(value_b)) NA_character_ else scalar_metric(value_b, metric),
    absolute_difference = abs_diff, relative_difference = rel_diff,
    equal = equal,
    status_a = status_value(result_a), status_b = status_value(result_b),
    warning_class_a = warning_value(result_a), warning_class_b = warning_value(result_b),
    failure_a = error_value(result_a), failure_b = error_value(result_b),
    retained_species_a = retained_value(result_a), retained_species_b = retained_value(result_b),
    fingerprint_a = fingerprint_value(ctx_a), fingerprint_b = fingerprint_value(ctx_b),
    rng_same_seed = rng_same_seed,
    stringsAsFactors = FALSE
  )
}

if (locale_c_suite$status == "RUN") {
  # check_tree and prepare_tree are public boundary observations.
  check_a <- locale_c_suite$check
  check_b <- locale_b_suite$check
  check_a_value <- if (is.list(check_a)) check_a$value else NULL
  check_b_value <- if (is.list(check_b)) check_b$value else NULL
  add_locale_row(
    "check_tree", "raw", "signal=K,lambda,D,Delta", "complete_object",
    if (is.list(check_a_value)) paste(scalar_metric(check_a_value$ready, "ready"), scalar_metric(check_a_value$ready_by_signal, "ready_by_signal"), collapse = "|") else NULL,
    if (is.list(check_b_value)) paste(scalar_metric(check_b_value$ready, "ready"), scalar_metric(check_b_value$ready_by_signal, "ready_by_signal"), collapse = "|") else NULL,
    check_a, check_b, rng_same_seed = NA
  )
  prep_a <- locale_c_suite$prepare
  prep_b <- locale_b_suite$prepare
  ctx_a <- locale_c_suite$ctx
  ctx_b <- locale_b_suite$ctx
  add_locale_row(
    "prepare_tree", "raw", "context", "fingerprint",
    if (is.list(ctx_a)) ctx_a$fingerprint else NULL,
    if (is.list(ctx_b)) ctx_b$fingerprint else NULL,
    prep_a, prep_b, ctx_a, ctx_b, rng_same_seed = NA
  )

  for (spec_name in names(continuous_specs)) {
    spec <- continuous_specs[[spec_name]]
    metrics <- if (spec$method == "K" && isTRUE(spec$test)) {
      c("estimate", "P", "MCSE", "matched_species", "retained_species", "status")
    } else if (spec$method == "K") {
      c("estimate", "matched_species", "retained_species", "status")
    } else {
      c("estimate", "logLik", "LR", "P", "matched_species", "retained_species", "status")
    }
    for (route in c("raw", "prepared")) {
      a <- locale_c_suite$continuous[[paste(spec_name, route, sep = "|")]]
      b <- locale_b_suite$continuous[[paste(spec_name, route, sep = "|")]]
      for (metric in metrics) {
        va <- result_metric(a$value, spec$method, metric)
        vb <- if (is.list(b)) result_metric(b$value, spec$method, metric) else NULL
        add_locale_row(spec$method, route, spec_name, metric, va, vb, a, b,
                       locale_c_suite$ctx, locale_b_suite$ctx,
                       rng_same_seed = TRUE)
      }
    }
  }

  for (test_value in c(FALSE, TRUE)) {
    metrics <- if (test_value) {
      c("estimate", "P_random", "P_Brownian", "MCSE_random", "MCSE_Brownian", "matched_species", "retained_species", "status")
    } else c("estimate", "matched_species", "retained_species", "status")
    for (route in c("raw", "prepared")) {
      a <- locale_c_suite$d[[paste(test_value, route, sep = "|")]]
      b <- locale_b_suite$d[[paste(test_value, route, sep = "|")]]
      for (metric in metrics) {
        va <- result_metric(a$value, "D", metric)
        vb <- if (is.list(b)) result_metric(b$value, "D", metric) else NULL
        add_locale_row("D", route, paste0("test=", test_value), metric, va, vb,
                       a, b, locale_c_suite$ctx, locale_b_suite$ctx,
                       rng_same_seed = TRUE)
      }
    }
  }

  for (test_value in c(FALSE, TRUE)) {
    metrics <- if (test_value) {
      c("estimate", "P", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status")
    } else c("estimate", "MCSE", "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "matched_species", "retained_species", "status")
    for (route in c("raw", "prepared")) {
      a <- locale_c_suite$delta[[paste(test_value, route, sep = "|")]]
      b <- locale_b_suite$delta[[paste(test_value, route, sep = "|")]]
      for (metric in metrics) {
        va <- result_metric(a$value, "Delta", metric)
        vb <- if (is.list(b)) result_metric(b$value, "Delta", metric) else NULL
        add_locale_row("Delta", route, paste0("test=", test_value), metric, va, vb,
                       a, b, locale_c_suite$ctx, locale_b_suite$ctx,
                       rng_same_seed = TRUE)
      }
    }
  }

  for (route in c("raw", "prepared")) {
    a <- locale_c_suite$ace[[route]]
    b <- locale_b_suite$ace[[route]]
    tree_a <- if (route == "prepared" && is.list(locale_c_suite$ctx)) locale_c_suite$ctx$tree else locale_tree
    tree_b <- if (route == "prepared" && is.list(locale_b_suite$ctx)) locale_b_suite$ctx$tree else locale_tree
    for (metric in c("logLik", "lik_anc", "rates", "se", "convergence", "message")) {
      va <- if (metric == "lik_anc") aligned_ace_lik_anc_stable(a$value, tree_a) else result_metric(a$value, "ACE", metric)
      vb <- if (is.list(b)) {
        if (metric == "lik_anc") aligned_ace_lik_anc_stable(b$value, tree_b) else result_metric(b$value, "ACE", metric)
      } else NULL
      add_locale_row("ACE", route, "CI=TRUE", metric, va, vb, a, b,
                     locale_c_suite$ctx, locale_b_suite$ctx,
                     tree_a = tree_a, tree_b = tree_b, rng_same_seed = NA)
    }
  }
}

locale_propagation <- if (length(locale_rows)) {
  keys <- unique(unlist(lapply(locale_rows, names)))
  do.call(rbind, lapply(locale_rows, function(x) {
    missing <- setdiff(keys, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[keys]
  }))
} else data.frame()
write.csv(locale_propagation, file.path(out_dir, "locale_propagation.csv"),
          row.names = FALSE, na = "")
try(Sys.setlocale("LC_COLLATE", locale_before), silent = TRUE)

comparison <- if (length(comparison_rows)) {
  keys <- unique(unlist(lapply(comparison_rows, names)))
  do.call(rbind, lapply(comparison_rows, function(x) {
    missing <- setdiff(keys, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[keys]
  }))
} else data.frame()
write.csv(comparison, file.path(out_dir, "propagation_comparisons.csv"), row.names = FALSE, na = "")

# Capture contract identity and downstream summary in a compact text file.
metric_pass <- if (nrow(comparison) && "value" %in% names(comparison)) {
  all(as.character(comparison$value) %in% c("TRUE", "True", "true"))
} else FALSE
error_count <- sum(nzchar(as.character(comparison$error %||% character())))
summary_lines <- c(
  "fastphylosig 0.2.0 Stage 2B2C downstream propagation audit",
  paste0("source_commit=", source_commit),
  paste0("package_version=", package_version),
  paste0("LC_COLLATE_before=", locale_before),
  paste0("LC_COLLATE_used=", Sys.getlocale("LC_COLLATE")),
  paste0("fixed_permutation_rows=", nrow(permutations)),
  paste0("fixed_brownian_states=", nrow(brownian_states), "x", ncol(brownian_states)),
  "delta_seed=20260906; delta ncores=1; ACE CI=TRUE",
  paste0("comparison_rows=", nrow(comparison)),
  paste0("all_comparison_rows_pass=", metric_pass),
  paste0("comparison_error_rows=", error_count),
  "D/Delta/ACE are recorded as executed results or explicit errors; no method was silently omitted.",
  paste0("output_dir=", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
)
writeLines(summary_lines, file.path(out_dir, "audit_summary.txt"), useBytes = TRUE)

provenance <- data.frame(
  item = c("source_commit", "package_version", "R_version", "platform", "LC_COLLATE_before", "LC_COLLATE_used", "RNGkind", "fixed_seed", "fixed_permutations", "fixed_brownian_states", "mcmc_settings", "production_code_changed", "tests_changed"),
  value = c(source_commit, package_version, R.version.string, R.version$platform,
            locale_before, Sys.getlocale("LC_COLLATE"), paste(RNGkind(), collapse = "|"),
            "20260906 for D/Delta; K fixed permutations", "4 rows x 4 indices", "4 x 4 finite matrix", "mcmc_sim=100, thin=2, burn=10", "NO", "NO"),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "provenance.csv"), row.names = FALSE, na = "")

message("Stage 2B2C downstream propagation audit completed.")
message("Output: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
