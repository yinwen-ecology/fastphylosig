#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2C.  Read-only canonicalization contract audit.
# This script writes only audit evidence and never changes package source or
# the existing test suite.

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args) >= 1L) args[[1L]] else getwd(),
                     winslash = "/", mustWork = TRUE)
out_dir <- if (length(args) >= 2L && nzchar(args[[2L]])) {
  normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(repo, "benchmarks", "stage2b2c", "results", "canonical")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
installed_lib <- if (length(args) >= 3L && nzchar(args[[3L]])) {
  normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
} else {
  Sys.getenv("FASTPHYLOSIG_STAGE2B2C_LIB", "")
}
if (nzchar(installed_lib) && dir.exists(file.path(installed_lib, "fastphylosig"))) {
  suppressPackageStartupMessages(
    library("fastphylosig", lib.loc = installed_lib, character.only = TRUE)
  )
  package_load <- paste0("installed library: ",
                         normalizePath(installed_lib, winslash = "/"))
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to load the source package.", call. = FALSE)
  }
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  package_load <- "pkgload::load_all(source)"
}
suppressPackageStartupMessages(library(ape))
ns <- asNamespace("fastphylosig")
get_internal <- function(name) get(name, envir = ns, inherits = FALSE)

source_commit <- if (length(args) >= 4L && nzchar(args[[4L]])) {
  args[[4L]]
} else {
  Sys.getenv("FASTPHYLOSIG_STAGE2B2C_SOURCE_COMMIT", "")
}
if (!nzchar(source_commit)) {
  source_commit <- tryCatch(
    system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE,
            stderr = FALSE)[[1L]],
    error = function(e) "unknown"
  )
}
if (!length(source_commit) || !nzchar(source_commit)) source_commit <- "unknown"

escape_text <- function(x) {
  if (!length(x) || is.na(x[[1L]])) return("NA")
  value <- as.character(x[[1L]])
  value <- gsub("\r", "\\r", value, fixed = TRUE)
  value <- gsub("\n", "\\n", value, fixed = TRUE)
  value
}
collapse_text <- function(x) {
  if (!length(x)) return("")
  paste(vapply(x, escape_text, character(1L)), collapse = ";")
}
collapse_integer <- function(x) {
  if (!length(x)) return("")
  paste(as.integer(x), collapse = ",")
}
collapse_numeric <- function(x) {
  if (!length(x)) return("")
  paste(sprintf("%.17g", as.numeric(x)), collapse = ",")
}
safe_serialize <- function(x) {
  tryCatch(serialize(x, NULL, version = 2L), error = function(e) NULL)
}

# Capture warning classes and errors without suppressing them in the evidence.
capture_call <- function(fun) {
  warnings <- character()
  warning_classes <- character()
  value <- NULL
  error_message <- ""
  error_classes <- character()
  status <- "ok"
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        warning_classes <<- c(warning_classes, class(w)[[1L]])
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      error_classes <<- class(e)
      NULL
    }
  )
  list(
    value = value,
    status = status,
    warnings = warnings,
    warning_classes = unique(warning_classes),
    error_message = error_message,
    error_classes = unique(error_classes)
  )
}

field_text <- function(x, name) {
  value <- if (is.list(x) && name %in% names(x)) x[[name]] else NULL
  if (is.null(value)) return("")
  if (is.numeric(value) || is.integer(value)) return(collapse_numeric(value))
  if (is.logical(value)) return(paste(as.character(value), collapse = ","))
  if (is.character(value)) return(collapse_text(value))
  "<complex>"
}

signature_text <- function(x) {
  if (!is.list(x)) return("NULL")
  paste(
    paste0("tip=", collapse_text(x$tip_label)),
    paste0("root=", collapse_integer(x$root)),
    paste0("root_desc=", escape_text(x$root_descendants)),
    paste0("root_degree=", collapse_integer(x$root_degree)),
    paste0("polytomies=", collapse_integer(x$polytomies)),
    paste0("edge=", collapse_integer(as.vector(x$edge))),
    paste0("length=", collapse_numeric(x$edge_length)),
    sep = "|"
  )
}

fingerprint <- get_internal(".tree_fingerprint")
canonical_signature <- get_internal(".canonical_tree_signature")
canonicalize <- get_internal(".safe_canonicalize_core")
inspect <- get_internal(".inspect_tree_core")
children_index <- get_internal(".edge_children_index")
descendant_keys <- get_internal(".descendant_keys_iterative")

tip_labels <- c("a", "b\rc", "a\rb", "c")

# Edge records are created from biological role names, so every representation
# uses the same branch length on the same biological edge.
make_collision_tree <- function(root = 5L, clade_a = 6L, clade_b = 7L,
                                order = c("root_a", "a_1", "a_2",
                                          "root_b", "b_3", "b_4")) {
  role <- list(
    root_a = c(root, clade_a, 1),
    a_1 = c(clade_a, 1L, 2),
    a_2 = c(clade_a, 2L, 3),
    root_b = c(root, clade_b, 4),
    b_3 = c(clade_b, 3L, 5),
    b_4 = c(clade_b, 4L, 6)
  )
  rows <- lapply(order, function(key) role[[key]])
  edge <- do.call(rbind, lapply(rows, function(z) as.integer(z[1:2])))
  edge.length <- vapply(rows, function(z) as.numeric(z[[3L]]), numeric(1L))
  structure(list(edge = edge, edge.length = edge.length,
                 tip.label = tip_labels, Nnode = 3L), class = "phylo")
}

trees <- list(
  base = make_collision_tree(),
  internal_renumbered = make_collision_tree(
    clade_a = 7L, clade_b = 6L,
    order = c("b_4", "root_a", "a_2", "root_b", "b_3", "a_1")
  ),
  edge_rows_shuffled = make_collision_tree(
    order = c("a_2", "root_b", "b_3", "root_a", "b_4", "a_1")
  )
)

tree_before <- lapply(trees, safe_serialize)

biological_equivalence <- function(x, y) {
  cx <- tryCatch(ape::cophenetic.phylo(x), error = function(e) NULL)
  cy <- tryCatch(ape::cophenetic.phylo(y), error = function(e) NULL)
  if (is.null(cx) || is.null(cy)) return(FALSE)
  isTRUE(all.equal(cx, cy, tolerance = 0, check.attributes = TRUE)) &&
    identical(as.character(x$tip.label), as.character(y$tip.label)) &&
    identical(sort(as.numeric(x$edge.length)),
              sort(as.numeric(y$edge.length)))
}

edge_records <- function(tree) {
  len <- if (is.null(tree$edge.length)) rep("<none>", nrow(tree$edge)) else
    sprintf("%.17g", as.numeric(tree$edge.length))
  sort(paste(tree$edge[, 1L], tree$edge[, 2L], len, sep = "\r"))
}

canonical_info <- lapply(trees, function(tree) {
  before <- canonical_signature(tree)
  out <- capture_call(function() canonicalize(tree))
  normalized <- out$value
  info <- if (is.list(normalized)) {
    attr(normalized, "fastphylosig_canonicalization", exact = TRUE)
  } else NULL
  after <- if (is.list(normalized)) canonical_signature(normalized) else NULL
  list(
    call = out,
    before = before,
    tree = normalized,
    info = info,
    after = after,
    before_fingerprint = fingerprint(tree),
    after_fingerprint = if (is.list(normalized)) fingerprint(normalized) else ""
  )
})

canonical_field_row <- function(pair_name, lhs_name, rhs_name) {
  lhs <- canonical_info[[lhs_name]]
  rhs <- canonical_info[[rhs_name]]
  left_tree <- lhs$tree
  right_tree <- rhs$tree
  left_info <- lhs$info
  right_info <- rhs$info
  flags <- c("changed", "tip_identity", "branch_identity", "edge_isomorphism",
             "cophenetic_identity", "root_identity", "polytomy_identity",
             "safe")
  flag_equal <- if (is.list(left_info) && is.list(right_info)) {
    identical(left_info[intersect(flags, names(left_info))],
              right_info[intersect(flags, names(right_info))])
  } else FALSE
  data.frame(
    pair = pair_name,
    lhs = lhs_name,
    rhs = rhs_name,
    biology_equivalent = biological_equivalence(trees[[lhs_name]],
                                                 trees[[rhs_name]]),
    canonical_edge_identical = is.list(left_tree) && is.list(right_tree) &&
      identical(left_tree$edge, right_tree$edge),
    canonical_edge_length_identical = is.list(left_tree) && is.list(right_tree) &&
      identical(as.numeric(left_tree$edge.length), as.numeric(right_tree$edge.length)),
    canonical_edge_length_association_identical = is.list(left_tree) &&
      is.list(right_tree) && identical(edge_records(left_tree), edge_records(right_tree)),
    canonical_tip_label_identical = is.list(left_tree) && is.list(right_tree) &&
      identical(as.character(left_tree$tip.label), as.character(right_tree$tip.label)),
    canonical_Nnode_identical = is.list(left_tree) && is.list(right_tree) &&
      identical(as.integer(left_tree$Nnode), as.integer(right_tree$Nnode)),
    mapping_identical = is.list(left_info) && is.list(right_info) &&
      identical(left_info$mapping, right_info$mapping),
    metadata_identical = is.list(left_info) && is.list(right_info) &&
      identical(left_info, right_info),
    metadata_safety_flags_identical = flag_equal,
    before_signature_identical = identical(lhs$before, rhs$before),
    after_signature_identical = identical(lhs$after, rhs$after),
    before_fingerprint_identical = identical(lhs$before_fingerprint,
                                             rhs$before_fingerprint),
    after_fingerprint_identical = identical(lhs$after_fingerprint,
                                            rhs$after_fingerprint),
    lhs_safe = is.list(left_info) && isTRUE(left_info$safe),
    rhs_safe = is.list(right_info) && isTRUE(right_info$safe),
    stringsAsFactors = FALSE
  )
}

pair_rows <- rbind(
  canonical_field_row("base_vs_internal_renumbered", "base", "internal_renumbered"),
  canonical_field_row("base_vs_edge_rows_shuffled", "base", "edge_rows_shuffled")
)
write.csv(pair_rows, file.path(out_dir, "representation_comparisons.csv"),
          row.names = FALSE, na = "")

canonical_detail <- do.call(rbind, lapply(names(canonical_info), function(nm) {
  z <- canonical_info[[nm]]
  info <- z$info
  data.frame(
    representation = nm,
    input_root = if (is.list(z$before)) z$before$root else NA_integer_,
    input_edge = if (is.list(z$before)) collapse_integer(as.vector(z$before$edge)) else "",
    input_lengths = if (is.list(z$before)) collapse_numeric(z$before$edge_length) else "",
    input_fingerprint = z$before_fingerprint,
    canonical_edge = if (is.list(z$tree)) collapse_integer(as.vector(z$tree$edge)) else "",
    canonical_lengths = if (is.list(z$tree)) collapse_numeric(z$tree$edge.length) else "",
    canonical_tip_label = if (is.list(z$tree)) collapse_text(z$tree$tip.label) else "",
    canonical_Nnode = if (is.list(z$tree)) as.integer(z$tree$Nnode) else NA_integer_,
    old_internal_order = if (is.list(info) && is.list(info$mapping))
      collapse_integer(info$mapping$old_internal_order) else "",
    new_internal_order = if (is.list(info) && is.list(info$mapping))
      collapse_integer(info$mapping$new_internal_order) else "",
    safe = is.list(info) && isTRUE(info$safe),
    changed = is.list(info) && isTRUE(info$changed),
    before_signature = signature_text(z$before),
    after_signature = signature_text(z$after),
    after_fingerprint = z$after_fingerprint,
    canonicalization_status = z$call$status,
    warning_classes = collapse_text(z$call$warning_classes),
    error_message = z$call$error_message,
    stringsAsFactors = FALSE
  )
}))
write.csv(canonical_detail, file.path(out_dir, "canonical_details.csv"),
          row.names = FALSE, na = "")

# Confirm the known collapsed-key collision and the legacy numeric tie-break
# without introducing a replacement ordering.
key_rows <- do.call(rbind, lapply(names(trees), function(nm) {
  tree <- trees[[nm]]
  roots <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
  root <- if (length(roots) == 1L) roots[[1L]] else NA_integer_
  internal_nodes <- setdiff(sort(unique(tree$edge[, 1L])), root)
  keys <- descendant_keys(
    internal_nodes,
    children_index(tree$edge), as.character(tree$tip.label), length(tree$tip.label)
  )
  data.frame(
    representation = nm,
    clade_keys = collapse_text(keys),
    clade_keys_identical = length(keys) == 2L && identical(keys[[1L]], keys[[2L]]),
    expected_key = "a\\rb\\rc",
    stringsAsFactors = FALSE
  )
}))
write.csv(key_rows, file.path(out_dir, "delimiter_collision_keys.csv"),
          row.names = FALSE, na = "")

# Locale audit.  Only LC_COLLATE is changed and it is restored on exit.
normal_locale <- Sys.getlocale("LC_COLLATE")
on.exit(try(Sys.setlocale("LC_COLLATE", normal_locale), silent = TRUE), add = TRUE)
locale_runs <- list()
run_locale <- function(label, requested) {
  set_result <- suppressWarnings(tryCatch(
    Sys.setlocale("LC_COLLATE", requested),
    error = function(e) NA_character_
  ))
  actual <- Sys.getlocale("LC_COLLATE")
  supported <- is.character(set_result) && length(set_result) == 1L &&
    !is.na(set_result) && identical(actual, set_result)
  if (!supported) {
    return(data.frame(
      locale = label, requested = requested, actual = actual,
      status = "NOT_RUN_ENVIRONMENT", labels = "", keys = "",
      old_internal_order = "", canonical_edge = "",
      canonical_fingerprint = "",
      stringsAsFactors = FALSE
    ))
  }
  locale_tree <- make_collision_tree()
  locale_tree$tip.label <- c(
    "A", "a0", "a1", "aa", "e\u00e9", "e\u00e8",
    "a-1", "a_1", "2", "10", "a\rb", "c"
  )
  # Expand a binary six-clade tree with a root.  The same labels are retained
  # in every locale; only R's documented ordering primitive is varied.
  pair_ids <- 13:18
  root <- 19L
  edge <- do.call(rbind, lapply(seq_along(pair_ids), function(i) {
    rbind(c(root, pair_ids[[i]]), c(pair_ids[[i]], 2L * i - 1L),
          c(pair_ids[[i]], 2L * i))
  }))
  locale_tree$edge <- edge
  locale_tree$edge.length <- as.numeric(seq_len(nrow(edge)))
  locale_tree$Nnode <- 7L
  children <- children_index(edge)
  internals <- pair_ids
  keys <- descendant_keys(internals, children, locale_tree$tip.label, 12L)
  normalized_call <- capture_call(function() canonicalize(locale_tree))
  normalized <- normalized_call$value
  normalized_info <- if (is.list(normalized)) {
    attr(normalized, "fastphylosig_canonicalization", exact = TRUE)
  } else NULL
  old_internal_order <- if (is.list(normalized_info) &&
                            is.list(normalized_info$mapping)) {
    normalized_info$mapping$old_internal_order
  } else integer()
  data.frame(
    locale = label, requested = requested, actual = actual,
    status = if (identical(label, "normal") || identical(actual, requested)) "PASS" else "NOT_RUN_ENVIRONMENT",
    labels = collapse_text(locale_tree$tip.label),
    keys = collapse_text(keys),
    old_internal_order = collapse_integer(old_internal_order),
    canonical_edge = if (is.list(normalized)) collapse_integer(as.vector(normalized$edge)) else "",
    canonical_fingerprint = if (is.list(normalized)) fingerprint(normalized) else "",
    stringsAsFactors = FALSE
  )
}

locale_runs[["C"]] <- run_locale("C", "C")
if (identical(normal_locale, "C")) {
  locale_runs[["normal"]] <- locale_runs[["C"]]
  locale_runs[["normal"]]$locale <- "normal"
  locale_runs[["normal"]]$status <- "SAME_AS_C"
} else {
  locale_runs[["normal"]] <- run_locale("normal", normal_locale)
}
optional_locales <- c(
  "English_United States.1252",
  "German_Germany.1252",
  "Chinese (Simplified)_China.936",
  "en_US.UTF-8"
)
for (candidate_locale in optional_locales) {
  label <- paste0("optional_", length(locale_runs) + 1L)
  locale_runs[[label]] <- run_locale(label, candidate_locale)
}
locale_rows <- do.call(rbind, locale_runs)
rownames(locale_rows) <- NULL
write.csv(locale_rows, file.path(out_dir, "locale_comparisons.csv"),
          row.names = FALSE, na = "")

locale_ok <- all(locale_rows$status %in%
                   c("PASS", "SAME_AS_C", "NOT_RUN_ENVIRONMENT"))
supported_actual <- unique(locale_rows$actual[locale_rows$status == "PASS"])
locale_has_two <- length(supported_actual) >= 2L
locale_dependent <- if (!locale_has_two || !locale_ok) {
  "PARTIALLY_TESTED"
} else {
  c_row <- locale_rows[locale_rows$locale == "C", , drop = FALSE][1L, ]
  alternatives <- locale_rows[
    locale_rows$status == "PASS" & locale_rows$actual != c_row$actual,
    , drop = FALSE
  ]
  differs <- vapply(seq_len(nrow(alternatives)), function(i) {
    row <- alternatives[i, , drop = FALSE]
    !identical(c_row$keys, row$keys) ||
      !identical(c_row$canonical_edge, row$canonical_edge) ||
      !identical(c_row$canonical_fingerprint, row$canonical_fingerprint)
  }, logical(1L))
  as.character(any(differs))
}

# Fixed inputs make K and D null accounting comparable across representations.
trait_cont <- c(0.2, 1.1, 2.4, 3.7)
names(trait_cont) <- tip_labels
trait_binary <- c(0, 1, 0, 1)
names(trait_binary) <- tip_labels
trait_cat <- c("A", "B", "A", "B")
names(trait_cat) <- tip_labels
permutations <- rbind(c(1L, 2L, 3L, 4L), c(2L, 1L, 4L, 3L),
                      c(4L, 3L, 2L, 1L), c(3L, 4L, 1L, 2L))
random_states <- cbind(c(0, 1, 0, 1), c(1, 0, 1, 0))
brownian_states <- cbind(c(0.1, 1.0, 2.0, 3.0), c(1.5, 0.5, 2.5, 3.5))

result_numeric <- function(value, names_to_keep) {
  if (!is.list(value)) return(setNames(numeric(), character()))
  keep <- intersect(names_to_keep, names(value))
  out <- unlist(lapply(keep, function(nm) {
    z <- value[[nm]]
    if (is.numeric(z) && length(z) == 1L) as.numeric(z) else NA_real_
  }), use.names = TRUE)
  out[!is.na(out)]
}
result_text <- function(value, name) {
  if (!is.list(value) || !name %in% names(value)) return("")
  z <- value[[name]]
  if (length(z) != 1L || !is.atomic(z)) return("")
  as.character(z[[1L]])
}
metadata_text <- function(value) {
  metadata <- attr(value, "analysis_metadata", exact = TRUE)
  if (!is.list(metadata)) return("")
  paste(
    paste0("matched=", field_text(metadata, "matched_species")),
    paste0("tree_removed=", field_text(metadata, "tree_tips_removed")),
    paste0("data_removed=", field_text(metadata, "data_rows_removed")),
    paste0("canonical=", field_text(metadata, "tree_auto_normalized")),
    sep = ";"
  )
}

message("audit checkpoint: before prepared contexts")
prepared_contexts <- lapply(names(trees), function(nm) {
  message("audit checkpoint: prepare ", nm)
  capture_call(function() prepare_tree(trees[[nm]]))
})
names(prepared_contexts) <- names(trees)
message("audit checkpoint: after prepared contexts")
context_rows <- do.call(rbind, lapply(names(prepared_contexts), function(nm) {
  z <- prepared_contexts[[nm]]
  ctx <- z$value
  ci <- if (is.list(ctx)) capture_call(function() cache_info(ctx)) else NULL
  data.frame(
    representation = nm,
    prepare_status = z$status,
    prepare_warnings = collapse_text(z$warning_classes),
    fingerprint = if (is.list(ctx)) ctx$fingerprint else "",
    cache_entries = if (is.list(ci) && is.list(ci$value)) ci$value$n_structural_entries else NA_integer_,
    cache_info_status = if (is.null(ci)) "not_run" else ci$status,
    stringsAsFactors = FALSE
  )
}))
write.csv(context_rows, file.path(out_dir, "prepared_contexts.csv"),
          row.names = FALSE, na = "")

# Exercise a protected-field mutation whose delimiter-collapsed tip-label
# bytes are unchanged. This tests context integrity directly, separately from
# comparing independently prepared representation-equivalent trees.
base_ctx <- prepared_contexts[["base"]]$value
mutated_ctx <- base_ctx
mutated_ctx$tree$tip.label <- c("a\rb", "c", "a", "b\rc")
mutation_fp <- if (is.list(mutated_ctx)) fingerprint(mutated_ctx$tree) else ""
mutation_validation <- capture_call(function() {
  get_internal(".validate_prepared_context")(mutated_ctx)
})
mutation_prepared_fit <- capture_call(function() {
  fast_k(mutated_ctx, trait_cont, test = FALSE, verbose = FALSE,
         progress = FALSE, ncores = 1L)
})
mutation_raw_fit <- capture_call(function() {
  fast_k(mutated_ctx$tree, trait_cont, test = FALSE, verbose = FALSE,
         progress = FALSE, ncores = 1L)
})
k_value <- function(z) {
  if (!identical(z$status, "ok") || is.null(z$value)) return(NA_real_)
  value <- suppressWarnings(as.numeric(z$value))
  if (length(value)) value[[1L]] else NA_real_
}
mutation_prepared_k <- k_value(mutation_prepared_fit)
mutation_raw_k <- k_value(mutation_raw_fit)
mutation_abs_diff <- abs(mutation_prepared_k - mutation_raw_k)
mutation_rel_diff <- mutation_abs_diff /
  max(abs(mutation_prepared_k), abs(mutation_raw_k), .Machine$double.eps)
mutation_rows <- data.frame(
  check = c(
    "tip_labels_changed",
    "tip_labels_remain_unique_nonempty",
    "fingerprint_collision",
    "prepared_context_validation_rejected",
    "prepared_vs_raw_mutated_K_identical"
  ),
  observed = c(
    as.character(!identical(base_ctx$tree$tip.label,
                            mutated_ctx$tree$tip.label)),
    as.character(!anyDuplicated(mutated_ctx$tree$tip.label) &&
                   all(nzchar(mutated_ctx$tree$tip.label))),
    as.character(identical(base_ctx$fingerprint, mutation_fp)),
    as.character(!identical(mutation_validation$status, "ok")),
    as.character(identical(mutation_prepared_k, mutation_raw_k))
  ),
  base_fingerprint = base_ctx$fingerprint,
  mutated_fingerprint = mutation_fp,
  validation_status = mutation_validation$status,
  validation_error = mutation_validation$error_message,
  prepared_K = mutation_prepared_k,
  raw_mutated_K = mutation_raw_k,
  absolute_difference = mutation_abs_diff,
  relative_difference = mutation_rel_diff,
  stringsAsFactors = FALSE
)
write.csv(mutation_rows,
          file.path(out_dir, "fingerprint_mutation_collision.csv"),
          row.names = FALSE, na = "")

public_specs <- list(
  fast_k = list(
    names = c("K", "P", "MCSE_P", "exceedance_count", "nsim_successful"),
    call = function(tree) fast_k(tree, trait_cont, test = TRUE, nsim = 4L,
                                  permutations = permutations, return_sim = TRUE,
                                  verbose = FALSE, progress = FALSE, ncores = 1L)
  ),
  fast_lambda = list(
    names = c("lambda", "logL", "LR", "P", "gls_mean", "sig2"),
    call = function(tree) fast_lambda(tree, trait_cont, test = TRUE,
                                       lambda_profile = FALSE, verbose = FALSE,
                                       progress = FALSE, ncores = 1L)
  ),
  fast_d = list(
    names = c("DEstimate", "Pval1", "Pval0", "P_random", "P_Brownian",
              "MCSE_P_random", "MCSE_P_Brownian", "observed", "mean_random",
              "mean_brownian", "nsim_successful_random",
              "nsim_successful_brownian"),
    call = function(tree) fast_d(tree, trait_binary, test = TRUE, nsim = 2L,
                                  random_states = random_states,
                                  brownian_states = brownian_states,
                                  return_sim = TRUE, verbose = FALSE,
                                  progress = FALSE, ncores = 1L)
  ),
  fast_delta = list(
    names = c("Delta", "P_fast", "P_MCSE", "MCSE_Delta", "alpha_mean",
              "beta_mean", "n_saved", "n_saved_successful"),
    call = function(tree) {
      set.seed(99173)
      fast_delta(tree, trait_cat, test = FALSE, mcmc_sim = 12L,
                 thin = 2L, burn = 2L, model = "ER", verbose = FALSE,
                 progress = FALSE, ncores = 1L)
    }
  ),
  fast_ace = list(
    names = c("ace", "loglik", "rates"),
    call = function(tree) fast_ace(trait_cat, phy = tree, CI = FALSE,
                                    marginal = FALSE, progress = FALSE)
  )
)

public_rows <- list()
public_values <- list()
for (method in names(public_specs)) {
  spec <- public_specs[[method]]
  method_values <- list()
  for (nm in names(trees)) {
    z <- if (method == "fast_ace") {
      normalized <- canonical_info[[nm]]$tree
      capture_call(function() spec$call(normalized))
    } else {
      capture_call(function() spec$call(trees[[nm]]))
    }
    method_values[[nm]] <- z
    public_rows[[length(public_rows) + 1L]] <- data.frame(
      method = method,
      representation = nm,
      input_route = if (method == "fast_ace") "canonicalized_tree" else "raw_tree",
      status = z$status,
      warning_classes = collapse_text(z$warning_classes),
      error_classes = collapse_text(z$error_classes),
      error_message = z$error_message,
      estimate_fields = paste(spec$names, collapse = ","),
      estimates = collapse_numeric(result_numeric(z$value, spec$names)),
      status_field = result_text(z$value, "status"),
      retained_species = result_text(z$value, "matched_species"),
      metadata = metadata_text(z$value),
      stringsAsFactors = FALSE
    )
  }
  public_values[[method]] <- method_values
}
public_rows <- do.call(rbind, public_rows)
write.csv(public_rows, file.path(out_dir, "downstream_propagation.csv"),
          row.names = FALSE, na = "")

propagation_rows <- list()
for (method in names(public_values)) {
  spec <- public_specs[[method]]
  ref <- public_values[[method]][["base"]]
  for (nm in setdiff(names(trees), "base")) {
    cur <- public_values[[method]][[nm]]
    lhs <- result_numeric(ref$value, spec$names)
    rhs <- result_numeric(cur$value, spec$names)
    all_names <- union(names(lhs), names(rhs))
    lhs <- lhs[all_names]
    rhs <- rhs[all_names]
    diffs <- abs(lhs - rhs)
    diffs[!is.finite(diffs)] <- NA_real_
    propagation_rows[[length(propagation_rows) + 1L]] <- data.frame(
      method = method,
      comparison = paste0("base_vs_", nm),
      base_status = ref$status,
      representation_status = cur$status,
      status_identical = identical(ref$status, cur$status),
      warning_classes_identical = identical(ref$warning_classes, cur$warning_classes),
      error_semantics_identical = identical(ref$error_classes, cur$error_classes) &&
        identical(ref$error_message, cur$error_message),
      numeric_fields = collapse_text(all_names),
      max_absolute_difference = if (length(diffs) && any(is.finite(diffs))) max(diffs, na.rm = TRUE) else NA_real_,
      numeric_exact_identical = identical(lhs, rhs),
      base_metadata = metadata_text(ref$value),
      representation_metadata = metadata_text(cur$value),
      stringsAsFactors = FALSE
    )
  }
}
propagation_rows <- do.call(rbind, propagation_rows)
write.csv(propagation_rows, file.path(out_dir, "downstream_comparisons.csv"),
          row.names = FALSE, na = "")

# Public check and prepared-context identity propagation are kept separate from
# estimator values because representation-specific fingerprints are part of the
# current safety contract rather than a biological estimate.
check_rows <- do.call(rbind, lapply(names(trees), function(nm) {
  z <- capture_call(function() check_tree(trees[[nm]],
                                          signal = c("K", "lambda", "D", "Delta")))
  v <- z$value
  data.frame(
    representation = nm,
    status = z$status,
    ready = if (is.list(v)) as.character(v$ready) else "",
    ready_by_signal = if (is.list(v)) paste(names(v$ready_by_signal),
                                            v$ready_by_signal, sep = "=", collapse = ";") else "",
    issue_codes = if (is.list(v) && is.data.frame(v$issues) && nrow(v$issues))
      paste(v$issues$code, collapse = ";") else "",
    warning_classes = collapse_text(z$warning_classes),
    error_message = z$error_message,
    stringsAsFactors = FALSE
  )
}))
write.csv(check_rows, file.path(out_dir, "check_tree_propagation.csv"),
          row.names = FALSE, na = "")

ctx_identity_rows <- do.call(rbind, lapply(setdiff(names(trees), "base"), function(nm) {
  base_ctx <- prepared_contexts[["base"]]$value
  other_ctx <- prepared_contexts[[nm]]$value
  data.frame(
    comparison = paste0("base_vs_", nm),
    prepared_status_identical = identical(prepared_contexts[["base"]]$status,
                                          prepared_contexts[[nm]]$status),
    fingerprint_identical = is.list(base_ctx) && is.list(other_ctx) &&
      identical(base_ctx$fingerprint, other_ctx$fingerprint),
    cache_key_names_identical = is.list(base_ctx) && is.list(other_ctx) &&
      identical(ls(base_ctx$structural_cache, all.names = TRUE),
                ls(other_ctx$structural_cache, all.names = TRUE)),
    retained_species_contract = "all four unique labels retained",
    stringsAsFactors = FALSE
  )
}))
write.csv(ctx_identity_rows, file.path(out_dir, "fingerprint_cache_identity.csv"),
          row.names = FALSE, na = "")

# Input immutability is checked after every canonicalization call.
immutability_rows <- do.call(rbind, lapply(names(trees), function(nm) {
  data.frame(
    representation = nm,
    input_immutable = identical(tree_before[[nm]], safe_serialize(trees[[nm]])),
    stringsAsFactors = FALSE
  )
}))
write.csv(immutability_rows, file.path(out_dir, "input_immutability.csv"),
          row.names = FALSE, na = "")

canonical_invariance <- pair_rows$canonical_edge_identical &
  pair_rows$canonical_edge_length_identical &
  pair_rows$canonical_edge_length_association_identical &
  pair_rows$canonical_tip_label_identical & pair_rows$canonical_Nnode_identical &
  pair_rows$after_signature_identical & pair_rows$after_fingerprint_identical
edge_order_row <- pair_rows$pair == "base_vs_edge_rows_shuffled"
renumber_row <- pair_rows$pair == "base_vs_internal_renumbered"
canonical_bug <- any(!canonical_invariance[renumber_row | edge_order_row])
fingerprint_mutation_collision <- identical(base_ctx$fingerprint, mutation_fp)
fingerprint_mutation_rejected <- !identical(mutation_validation$status, "ok")
stale_cache_propagated <- is.finite(mutation_abs_diff) && mutation_abs_diff > 0
canonical_bug <- canonical_bug ||
  (fingerprint_mutation_collision && !fingerprint_mutation_rejected)
numeric_propagation <- propagation_rows[,
  c("method", "comparison", "status_identical", "warning_classes_identical",
    "error_semantics_identical", "max_absolute_difference",
    "numeric_exact_identical"), drop = FALSE]

summary <- data.frame(
  item = c(
    "source_commit", "package_version", "delimiter_collision_fixture",
    "CANONICAL_RENUMBERING_INVARIANCE", "CANONICAL_EDGE_ORDER_INVARIANCE",
    "CANONICALIZATION_LOCALE_DEPENDENT", "fingerprint_contract",
    "FINGERPRINT_MUTATION_COLLISION", "FINGERPRINT_MUTATION_REJECTED",
    "STALE_CACHE_NUMERICAL_PROPAGATION",
    "downstream_fixed_input_audit", "input_immutability",
    "CANONICALIZATION_CONTRACT_BUG", "PUBLIC_RELEASE_BLOCKER_0_2_0",
    "NEXT_STEP"
  ),
  value = c(
    source_commit,
    as.character(utils::packageVersion("fastphylosig")),
    as.character(all(key_rows$clade_keys_identical)),
    if (all(canonical_invariance[renumber_row])) "PASS" else "FAIL",
    if (all(canonical_invariance[edge_order_row])) "PASS" else "FAIL",
    locale_dependent,
    "raw tip/edge/length/Nnode fingerprint is representation-sensitive",
    as.character(fingerprint_mutation_collision),
    as.character(fingerprint_mutation_rejected),
    as.character(stale_cache_propagated),
    if (all(propagation_rows$status_identical)) "PASS" else "PARTIAL",
    if (all(immutability_rows$input_immutable)) "PASS" else "FAIL",
    if (canonical_bug) "YES" else "NO",
    if (canonical_bug) "YES" else "NO",
    if (canonical_bug) "Canonicalization Contract v2" else "Stage 2C"
  ),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, "audit_summary.csv"), row.names = FALSE, na = "")

provenance <- data.frame(
  item = c("source_commit", "R_version", "platform", "locale_normal",
           "package_load", "script", "scope"),
  value = c(source_commit, R.version.string, R.version$platform, normal_locale,
            package_load, "audit_canonical_contract.R",
            "read-only canonicalization/fingerprint contract audit"),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "provenance.csv"), row.names = FALSE)

writeLines(c(
  "fastphylosig 0.2.0 Stage 2B2C Canonicalization Contract Safety Audit",
  paste0("source_commit=", source_commit),
  paste0("CANONICALIZATION_CONTRACT_BUG=", if (canonical_bug) "YES" else "NO"),
  paste0("PUBLIC_RELEASE_BLOCKER_0_2_0=", if (canonical_bug) "YES" else "NO"),
  paste0("NEXT_STEP=", if (canonical_bug) "Canonicalization Contract v2" else "Stage 2C"),
  paste0("CANONICAL_RENUMBERING_INVARIANCE=",
         if (all(canonical_invariance[renumber_row])) "PASS" else "FAIL"),
  paste0("CANONICAL_EDGE_ORDER_INVARIANCE=",
         if (all(canonical_invariance[edge_order_row])) "PASS" else "FAIL"),
  paste0("CANONICALIZATION_LOCALE_DEPENDENT=", locale_dependent),
  paste0("FINGERPRINT_MUTATION_COLLISION=", fingerprint_mutation_collision),
  paste0("FINGERPRINT_MUTATION_REJECTED=", fingerprint_mutation_rejected),
  paste0("STALE_CACHE_NUMERICAL_PROPAGATION=", stale_cache_propagated),
  paste0("evidence_dir=", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
), file.path(out_dir, "audit_status.txt"), useBytes = TRUE)

message("Stage 2B2C canonicalization audit completed")
message("Summary: ", normalizePath(file.path(out_dir, "audit_summary.csv"),
                                  winslash = "/", mustWork = FALSE))
