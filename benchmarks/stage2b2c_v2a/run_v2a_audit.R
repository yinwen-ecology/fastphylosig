#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Canonicalization Contract V2-A audit.
#
# This file is deliberately audit-only.  It sources a private prototype from
# this directory, compares exact canonical fields, and exercises the frozen
# production package as an estimator oracle.  It never edits production code,
# tests, or the supplied tree/trait objects.

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(
  if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else getwd(),
  winslash = "/", mustWork = TRUE
)
script_dir <- normalizePath(
  file.path(repo, "benchmarks", "stage2b2c_v2a"),
  winslash = "/", mustWork = TRUE
)
if (length(args) < 2L || !nzchar(args[[2L]])) {
  stop(
    "usage: run_v2a_audit.R <repo> <output-dir> [installed-lib] [source-commit]",
    call. = FALSE
  )
}
out_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

`%||%` <- function(x, y) if (is.null(x)) y else x

installed_lib <- if (length(args) >= 3L && nzchar(args[[3L]])) {
  normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
} else {
  Sys.getenv("FASTPHYLOSIG_V2A_LIB", "")
}
if (nzchar(installed_lib) && dir.exists(file.path(installed_lib, "fastphylosig"))) {
  if (!requireNamespace("fastphylosig", quietly = TRUE,
                       lib.loc = installed_lib)) {
    suppressPackageStartupMessages(
      library("fastphylosig", lib.loc = installed_lib, character.only = TRUE)
    )
  } else {
    suppressPackageStartupMessages(
      library("fastphylosig", lib.loc = installed_lib, character.only = TRUE)
    )
  }
  package_load <- paste0("installed library: ", installed_lib)
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to load the source package.", call. = FALSE)
  }
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  package_load <- "pkgload::load_all(source)"
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("ape is required for the V2-A audit fixtures.", call. = FALSE)
}
suppressPackageStartupMessages(library(ape))

prototype_file <- file.path(script_dir, "v2_prototype.R")
if (!file.exists(prototype_file)) {
  stop(
    paste0(
      "V2 prototype is missing: ", prototype_file,
      "\nThe audit is fail-closed; no evidence was generated."
    ),
    call. = FALSE
  )
}
v2_env <- new.env(parent = globalenv())
sys.source(prototype_file, envir = v2_env)
required_v2 <- c(
  "v2_canonicalize", "v2_fingerprint", "v2_protected_snapshot",
  "v2_prepare_context", "v2_validate_context"
)
missing_v2 <- required_v2[!vapply(
  required_v2, exists, logical(1L), envir = v2_env, inherits = FALSE
)]
if (length(missing_v2)) {
  stop(
    paste0("v2_prototype.R is missing required interface(s): ",
           paste(missing_v2, collapse = ", ")),
    call. = FALSE
  )
}

source_commit <- if (length(args) >= 4L && nzchar(args[[4L]])) {
  args[[4L]]
} else {
  detected <- tryCatch(
    system2("git", c("-C", repo, "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(detected) == 1L && nzchar(detected)) detected else "unknown"
}

clone <- function(x) unserialize(serialize(x, connection = NULL, version = 2L))

capture <- function(fun) {
  warnings <- list()
  captured_error <- NULL
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(condition) {
        warnings[[length(warnings) + 1L]] <<- list(
          message = conditionMessage(condition), class = class(condition)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      captured_error <<- list(
        message = conditionMessage(condition), class = class(condition)
      )
      NULL
    }
  )
  list(
    value = value,
    status = if (is.null(captured_error)) "ok" else "error",
    warnings = warnings,
    warning_text = if (length(warnings)) {
      paste(vapply(warnings, `[[`, character(1L), "message"), collapse = " || ")
    } else "",
    warning_classes = if (length(warnings)) {
      paste(vapply(warnings, function(x) paste(x$class, collapse = "/"),
                  character(1L)), collapse = " || ")
    } else "",
    error = captured_error,
    error_text = if (is.null(captured_error)) "" else captured_error$message,
    error_classes = if (is.null(captured_error)) "" else
      paste(captured_error$class, collapse = "/")
  )
}

collapse_text <- function(x) {
  if (is.null(x) || !length(x)) return("")
  value <- as.character(x)
  value <- gsub("\\r", "\\\\r", value, fixed = TRUE)
  value <- gsub("\\n", "\\\\n", value, fixed = TRUE)
  paste(value, collapse = ";")
}

scalar_text <- function(x) {
  if (is.null(x)) return("")
  if (is.logical(x)) return(paste(as.logical(x), collapse = ","))
  if (is.integer(x) || is.numeric(x)) {
    return(paste(sprintf("%.17g", as.numeric(x)), collapse = ","))
  }
  if (is.character(x)) return(collapse_text(x))
  paste0("<", paste(class(x), collapse = "/"), ">")
}

safe_identical <- function(x, y) isTRUE(identical(x, y))

first_numeric <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (!length(value) || is.na(value)) NA_real_ else value[[1L]]
}

numeric_difference <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      !length(x)) return(NA_real_)
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  max(abs(as.numeric(x[ok]) - as.numeric(y[ok])))
}

relative_difference <- function(x, y) {
  d <- numeric_difference(x, y)
  if (!is.finite(d) || !is.numeric(y) || !length(y)) return(NA_real_)
  max(abs(as.numeric(x) - as.numeric(y)) /
        pmax(1, abs(as.numeric(y))), na.rm = TRUE)
}

parity_equal <- function(x, y, tolerance = 1e-8) {
  if (safe_identical(x, y)) return(TRUE)
  if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
    ok <- is.finite(x) & is.finite(y)
    same_nonfinite <- identical(is.finite(x), is.finite(y)) &&
      all(is.na(x) == is.na(y))
    if (!same_nonfinite) return(FALSE)
    if (!any(ok)) return(TRUE)
    return(all(abs(as.numeric(x[ok]) - as.numeric(y[ok])) <=
                 tolerance * pmax(1, abs(as.numeric(y[ok])))))
  }
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}

make_pectinate <- function(n) {
  stopifnot(length(n) == 1L, n >= 2L)
  n <- as.integer(n)
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      parent <- n + i
      edge[2L * i - 1L, ] <- c(parent, i)
      edge[2L * i, ] <- c(parent, parent + 1L)
    }
  }
  root <- 2L * n - 1L
  edge[2L * n - 3L, ] <- c(root, n - 1L)
  edge[2L * n - 2L, ] <- c(root, n)
  structure(
    list(edge = edge, edge.length = rep(1, nrow(edge)),
         tip.label = paste0("tip", seq_len(n)), Nnode = n - 1L),
    class = "phylo"
  )
}

set_lengths <- function(tree, offset = 0) {
  tree$edge.length <- as.numeric(seq_len(nrow(tree$edge)) + offset) /
    (nrow(tree$edge) + 3)
  tree
}

shuffle_edges <- function(tree, seed = 1L) {
  out <- clone(tree)
  set.seed(seed)
  rows <- sample(seq_len(nrow(out$edge)))
  out$edge <- out$edge[rows, , drop = FALSE]
  if (!is.null(out$edge.length)) out$edge.length <- out$edge.length[rows]
  out
}

renumber_internals <- function(tree, permutation) {
  out <- clone(tree)
  n_tip <- length(out$tip.label)
  old <- seq.int(n_tip + 1L, n_tip + as.integer(out$Nnode))
  stopifnot(length(permutation) == length(old), setequal(permutation, old))
  map <- stats::setNames(as.integer(permutation), as.character(old))
  internal <- out$edge > n_tip
  out$edge[internal] <- unname(map[as.character(out$edge[internal])])
  if (!is.null(out$node.label)) {
    labels <- character(length(old))
    labels[as.integer(map) - n_tip] <- out$node.label
    out$node.label <- labels
  }
  out
}

make_collision_tree <- function() {
  structure(
    list(
      edge = matrix(c(
        5L, 6L, 6L, 1L, 6L, 2L,
        5L, 7L, 7L, 3L, 7L, 4L
      ), ncol = 2L, byrow = TRUE),
      edge.length = c(1.10, 0.70, 0.90, 1.20, 0.80, 1.05),
      tip.label = c("a", "b\rc", "a\rb", "c"), Nnode = 3L
    ),
    class = "phylo"
  )
}

make_label_tree <- function() {
  tree <- set_lengths(ape::stree(16L, type = "balanced"), offset = 11L)
  tree$tip.label <- enc2utf8(c(
    "Alpha", "alpha", "alpha0", "alpha-1", "alpha_1", "2", "10",
    "e\u00e9", "e\u00e8", "\u4e2d", "\u03b1", "a.prefix", "a prefix",
    "[alpha]", "01", "z"
  ))
  tree
}

make_random_tree <- function(n = 17L) {
  set.seed(20260906 + as.integer(n))
  tree <- ape::rtree(n)
  tree$tip.label <- paste0("tip", sprintf("%03d", seq_len(n)))
  tree$node.label <- paste0("node", seq_len(tree$Nnode))
  set_lengths(tree, offset = 7L)
}

random_base <- make_random_tree(17L)
random_internal <- seq.int(length(random_base$tip.label) + 1L,
                           length(random_base$tip.label) + random_base$Nnode)
random_permutation <- c(
  random_internal[[length(random_internal)]],
  random_internal[-length(random_internal)]
)
random_renumbered <- renumber_internals(random_base, random_permutation)
random_edge_shuffled <- shuffle_edges(random_base, seed = 17L)
random_renumbered_shuffled <- shuffle_edges(random_renumbered, seed = 18L)

collision_base <- make_collision_tree()
collision_internal <- seq.int(5L, 7L)
collision_renumbered <- renumber_internals(
  collision_base, c(7L, 5L, 6L)
)
collision_edge_shuffled <- shuffle_edges(collision_base, seed = 19L)
collision_both <- shuffle_edges(collision_renumbered, seed = 20L)

valid_fixtures <- list(
  balanced = set_lengths(ape::stree(16L, type = "balanced"), offset = 1L),
  random = random_base,
  pectinate = set_lengths(make_pectinate(40L), offset = 2L),
  polytomy = set_lengths(ape::stree(12L, type = "star"), offset = 3L),
  two_tip = set_lengths(ape::read.tree(text = "(a:1,b:2);"), offset = 4L),
  delimiter = collision_base,
  delimiter_edge_shuffled = collision_edge_shuffled,
  delimiter_internal_renumbered = collision_renumbered,
  delimiter_renumbered_shuffled = collision_both,
  labels_all_classes = make_label_tree(),
  edge_shuffled = random_edge_shuffled,
  internal_renumbered = random_renumbered,
  alternate_root_numbering = random_renumbered,
  renumbered_and_shuffled = random_renumbered_shuffled
)

equivalence_groups <- list(
  random_representations = c(
    "random", "edge_shuffled", "internal_renumbered",
    "alternate_root_numbering", "renumbered_and_shuffled"
  ),
  delimiter_representations = c(
    "delimiter", "delimiter_edge_shuffled", "delimiter_internal_renumbered",
    "delimiter_renumbered_shuffled"
  )
)

failure_fixtures <- list(
  not_phylo = list(edge = 1),
  missing_labels = {
    z <- clone(random_base); z$tip.label <- NULL; z
  },
  duplicate_labels = {
    z <- clone(random_base); z$tip.label[[2L]] <- z$tip.label[[1L]]; z
  },
  missing_label_value = {
    z <- clone(random_base); z$tip.label[[2L]] <- NA_character_; z
  },
  empty_label = {
    z <- clone(random_base); z$tip.label[[2L]] <- ""; z
  },
  malformed_edge = {
    z <- clone(random_base); z$edge <- matrix(seq_len(9L), nrow = 3L); z
  },
  invalid_nnode = {
    z <- clone(random_base); z$Nnode <- NA_integer_; z
  },
  negative_branch = {
    z <- clone(random_base); z$edge.length[[1L]] <- -1; z
  },
  nonfinite_branch = {
    z <- clone(random_base); z$edge.length[[1L]] <- Inf; z
  },
  cycle = {
    z <- clone(random_base); z$edge[1L, 2L] <- z$edge[1L, 1L]; z
  },
  disconnected = {
    structure(list(
      edge = matrix(c(5L, 1L, 5L, 2L, 6L, 3L, 6L, 4L), ncol = 2L,
                    byrow = TRUE), edge.length = rep(1, 4L),
      tip.label = paste0("tip", 1:4), Nnode = 2L), class = "phylo")
  }
)

get_v2 <- function(name) get(name, envir = v2_env, inherits = FALSE)
call_v2 <- function(name, argument) get_v2(name)(argument)

is_tree_like <- function(x) {
  is.list(x) && !is.null(x$edge) && !is.null(x$tip.label) &&
    !is.null(x$Nnode) && is.matrix(x$edge) && ncol(x$edge) == 2L
}

as_tree_candidate <- function(x) {
  if (inherits(x, "phylo") && is_tree_like(x)) return(x)
  if (!is.list(x)) return(NULL)
  direct <- c("tree", "canonical_tree", "phylo", "canonical")
  for (nm in direct) {
    if (!is.null(x[[nm]])) {
      candidate <- x[[nm]]
      if (inherits(candidate, "phylo") && is_tree_like(candidate)) {
        return(candidate)
      }
      if (is_tree_like(candidate)) {
        class(candidate) <- unique(c("phylo", class(candidate)))
        return(candidate)
      }
    }
  }
  if (is_tree_like(x)) {
    class(x) <- unique(c("phylo", class(x)))
    return(x)
  }
  NULL
}

metadata_candidates <- function(x) {
  if (!is.list(x)) return(list())
  out <- list(x)
  for (nm in c("canonical_fields", "fields", "representation",
               "canonical_metadata", "metadata", "provenance")) {
    if (is.list(x[[nm]])) out[[length(out) + 1L]] <- x[[nm]]
  }
  for (nm in c("fastphylosig_v2_contract", "fastphylosig_v2_metrics")) {
    value <- attr(x, nm, exact = TRUE)
    if (is.list(value)) out[[length(out) + 1L]] <- value
  }
  out
}

find_field <- function(x, aliases) {
  for (candidate in metadata_candidates(x)) {
    for (nm in aliases) {
      if (nm %in% names(candidate)) return(candidate[[nm]])
    }
  }
  NULL
}

extract_canonical <- function(capture_result) {
  raw <- capture_result$value
  tree <- as_tree_candidate(raw)
  field_aliases <- list(
    edge = c("edge", "canonical_edge"),
    edge.length = c("edge.length", "edge_length", "canonical_edge_length"),
    tip.label = c("tip.label", "tip_label", "canonical_tip_label"),
    Nnode = c("Nnode", "n_node", "nnode", "canonical_Nnode"),
    internal_mapping = c(
      "internal_mapping", "canonical_internal_mapping",
      "canonical_ids", "canonical_internal_ids"
    ),
    subtree_descriptors = c(
      "subtree_descriptors", "descriptor_table", "canonical_descriptors"
    ),
    edge_order = c("edge_order", "canonical_edge_order"),
    root = c("root", "canonical_root", "root_id")
  )
  fields <- lapply(field_aliases, function(aliases) find_field(raw, aliases))
  if (!is.null(tree)) {
    if (is.null(fields$edge)) fields$edge <- tree$edge
    if (is.null(fields$edge.length)) fields$edge.length <- tree$edge.length
    if (is.null(fields$tip.label)) fields$tip.label <- tree$tip.label
    if (is.null(fields$Nnode)) fields$Nnode <- tree$Nnode
    if (is.null(fields$root)) {
      fields$root <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
    }
    if (is.null(fields$edge_order)) fields$edge_order <- seq_len(nrow(tree$edge))
  }
  fields$edge <- if (is.null(fields$edge)) NULL else as.matrix(fields$edge)
  fields$edge.length <- if (is.null(fields$edge.length)) NULL else
    as.numeric(fields$edge.length)
  fields$tip.label <- if (is.null(fields$tip.label)) NULL else
    as.character(fields$tip.label)
  fields$Nnode <- if (is.null(fields$Nnode)) NULL else as.integer(fields$Nnode)
  fields$edge_order <- if (is.null(fields$edge_order)) NULL else
    as.integer(fields$edge_order)
  fields$root <- if (is.null(fields$root)) NULL else as.integer(fields$root)
  complexity <- list(
    child_rank_tuple_entries = find_field(raw, c(
      "child_rank_tuple_entries", "child_rank_entries",
      "n_child_rank_tuple_entries", "n_child_rank_entries"
    )),
    descriptor_entries = find_field(raw, c(
      "descriptor_entries", "n_descriptor_entries", "stored_descriptor_entries"
    )),
    peak_internal_object_proxy = find_field(raw, c(
      "peak_internal_object_proxy", "peak_descriptor_entries",
      "peak_internal_proxy", "peak_object_proxy"
    ))
  )
  list(
    raw = raw, tree = tree, fields = fields, complexity = complexity,
    missing_required = names(fields)[vapply(fields, is.null, logical(1L))]
  )
}

canonical_fields_identical <- function(lhs, rhs) {
  if (is.null(lhs) || is.null(rhs)) return(FALSE)
  if (length(lhs$missing_required) || length(rhs$missing_required)) return(FALSE)
  names_to_compare <- union(names(lhs$fields), names(rhs$fields))
  all(vapply(names_to_compare, function(nm) {
    nm %in% names(lhs$fields) && nm %in% names(rhs$fields) &&
      safe_identical(lhs$fields[[nm]], rhs$fields[[nm]])
  }, logical(1L)))
}

canonical_field_flags <- function(lhs, rhs) {
  names_to_compare <- union(names(lhs$fields), names(rhs$fields))
  flags <- vapply(names_to_compare, function(nm) {
    nm %in% names(lhs$fields) && nm %in% names(rhs$fields) &&
      safe_identical(lhs$fields[[nm]], rhs$fields[[nm]])
  }, logical(1L))
  list(
    all = length(lhs$missing_required) == 0L &&
      length(rhs$missing_required) == 0L && all(flags),
    fields = flags,
    missing_lhs = lhs$missing_required,
    missing_rhs = rhs$missing_required
  )
}

fingerprint_value <- function(capture_result) {
  if (is.null(capture_result) || !identical(capture_result$status, "ok")) {
    return(NULL)
  }
  value <- capture_result$value
  if (is.list(value) && !is.null(value$fingerprint)) value <- value$fingerprint
  value
}

fingerprint_text <- function(capture_result) {
  if (!identical(capture_result$status, "ok")) return("")
  value <- fingerprint_value(capture_result)
  if (is.raw(value)) return(paste(sprintf("%02x", as.integer(value)), collapse = ""))
  if (length(value) == 1L && (is.character(value) || is.numeric(value))) {
    return(scalar_text(value))
  }
  encoded <- serialize(value, NULL, version = 2L)
  paste(sprintf("%02x", as.integer(encoded)), collapse = "")
}

fingerprint_identical <- function(lhs, rhs) {
  safe_identical(fingerprint_value(lhs), fingerprint_value(rhs))
}

run_canonical <- function(tree) {
  local_tree <- clone(tree)
  before <- serialize(local_tree, NULL, version = 2L)
  z <- capture(function() call_v2("v2_canonicalize", local_tree))
  input_immutable <- identical(
    before, serialize(local_tree, NULL, version = 2L)
  )
  c <- if (identical(z$status, "ok")) extract_canonical(z) else NULL
  fp <- capture(function() call_v2("v2_fingerprint", clone(tree)))
  stored_tip_label_identical <- !is.null(c) && !is.null(c$tree) &&
    identical(c$tree$tip.label, tree$tip.label)
  list(
    call = z, canonical = c, fingerprint = fp,
    input_immutable = input_immutable,
    stored_tip_label_identical = stored_tip_label_identical
  )
}

canonical_runs <- lapply(valid_fixtures, run_canonical)

canonical_detail <- do.call(rbind, lapply(names(canonical_runs), function(nm) {
  z <- canonical_runs[[nm]]
  c <- z$canonical
  data.frame(
    fixture = nm,
    canonical_status = z$call$status,
    canonical_warnings = z$call$warning_classes,
    canonical_error = z$call$error_text,
    fingerprint_status = z$fingerprint$status,
    fingerprint = fingerprint_text(z$fingerprint),
    canonical_tree_present = !is.null(c) && !is.null(c$tree),
    input_immutable = isTRUE(z$input_immutable),
    stored_tip_label_identical = isTRUE(z$stored_tip_label_identical),
    missing_required_fields = if (is.null(c)) "" else
      paste(c$missing_required, collapse = ","),
    child_rank_tuple_entries = if (is.null(c)) NA_real_ else
      first_numeric(c$complexity$child_rank_tuple_entries),
    descriptor_entries = if (is.null(c)) NA_real_ else
      first_numeric(c$complexity$descriptor_entries),
    peak_internal_object_proxy = if (is.null(c)) NA_real_ else
      first_numeric(c$complexity$peak_internal_object_proxy),
    stringsAsFactors = FALSE
  )
}))
write.csv(canonical_detail, file.path(out_dir, "v2_canonical_details.csv"),
          row.names = FALSE, na = "")

invariance_rows <- list()
for (group_name in names(equivalence_groups)) {
  members <- equivalence_groups[[group_name]]
  reference <- canonical_runs[[members[[1L]]]]
  for (member in members[-1L]) {
    candidate <- canonical_runs[[member]]
    flags <- if (!is.null(reference$canonical) && !is.null(candidate$canonical)) {
      canonical_field_flags(reference$canonical, candidate$canonical)
    } else list(all = FALSE, fields = logical(), missing_lhs = character(),
                missing_rhs = character())
    fp_lhs <- fingerprint_text(reference$fingerprint)
    fp_rhs <- fingerprint_text(candidate$fingerprint)
    invariance_rows[[length(invariance_rows) + 1L]] <- data.frame(
      group = group_name, reference = members[[1L]], candidate = member,
      biology_fixture = TRUE,
      canonical_fields_identical = flags$all,
      edge_identical = isTRUE(flags$fields[["edge"]]),
      edge_length_identical = isTRUE(flags$fields[["edge.length"]]),
      tip_label_identical = isTRUE(flags$fields[["tip.label"]]),
      Nnode_identical = isTRUE(flags$fields[["Nnode"]]),
      internal_mapping_identical = isTRUE(flags$fields[["internal_mapping"]]),
      edge_order_identical = isTRUE(flags$fields[["edge_order"]]),
      root_identical = isTRUE(flags$fields[["root"]]),
      fingerprint_identical = fingerprint_identical(
        reference$fingerprint, candidate$fingerprint
      ),
      reference_missing = paste(flags$missing_lhs, collapse = ","),
      candidate_missing = paste(flags$missing_rhs, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
}
invariance_rows <- do.call(rbind, invariance_rows)
write.csv(invariance_rows,
          file.path(out_dir, "v2_representation_invariance.csv"),
          row.names = FALSE, na = "")

# Locale audit: only LC_COLLATE is changed.  Unsupported locale requests are
# retained as explicit NOT_RUN_ENVIRONMENT rows and are never treated as pass.
normal_locale <- Sys.getlocale("LC_COLLATE")
on.exit(try(Sys.setlocale("LC_COLLATE", normal_locale), silent = TRUE), add = TRUE)
locale_requests <- c(
  C = "C",
  normal = normal_locale,
  English_United_States = "English_United States.1252",
  German_Germany = "German_Germany.1252",
  Chinese_Simplified = "Chinese (Simplified)_China.936",
  en_US_UTF8 = "en_US.UTF-8"
)
locale_rows <- list()
locale_reference <- list()
for (locale_name in names(locale_requests)) {
  requested <- locale_requests[[locale_name]]
  set_result <- suppressWarnings(tryCatch(
    Sys.setlocale("LC_COLLATE", requested), error = function(e) NA_character_
  ))
  actual <- Sys.getlocale("LC_COLLATE")
  supported <- is.character(set_result) && length(set_result) == 1L &&
    !is.na(set_result) && identical(actual, set_result)
  if (!supported) {
    locale_rows[[length(locale_rows) + 1L]] <- data.frame(
      locale = locale_name, requested = requested, actual = actual,
      status = "NOT_RUN_ENVIRONMENT", fixture = "labels_all_classes",
      canonical_fields_identical_to_C = NA, fingerprint_identical_to_C = NA,
      missing_required_fields = "", error = "", stringsAsFactors = FALSE
    )
    next
  }
  runs <- run_canonical(valid_fixtures[["labels_all_classes"]])
  locale_reference[[locale_name]] <- runs
  c_ref <- locale_reference[["C"]]
  fields_same <- if (is.null(c_ref$canonical) || is.null(runs$canonical)) FALSE else
    canonical_fields_identical(c_ref$canonical, runs$canonical)
  fp_same <- fingerprint_identical(c_ref$fingerprint, runs$fingerprint)
  locale_rows[[length(locale_rows) + 1L]] <- data.frame(
    locale = locale_name, requested = requested, actual = actual,
    status = "PASS", fixture = "labels_all_classes",
    canonical_fields_identical_to_C = fields_same,
    fingerprint_identical_to_C = fp_same,
    missing_required_fields = if (is.null(runs$canonical)) "" else
      paste(runs$canonical$missing_required, collapse = ","),
    error = runs$call$error_text, stringsAsFactors = FALSE
  )
}
# The C row is the reference even if the current locale was processed first.
if (!is.null(locale_reference[["C"]])) {
  c_run <- locale_reference[["C"]]
  for (i in seq_along(locale_rows)) {
    row <- locale_rows[[i]]
    if (identical(row$status, "PASS")) {
      run <- locale_reference[[row$locale]]
      row$canonical_fields_identical_to_C <- if (is.null(run$canonical) ||
                                                   is.null(c_run$canonical)) FALSE else
        canonical_fields_identical(c_run$canonical, run$canonical)
      row$fingerprint_identical_to_C <- fingerprint_identical(
        c_run$fingerprint, run$fingerprint
      )
      locale_rows[[i]] <- row
    }
  }
}
locale_rows <- do.call(rbind, locale_rows)
write.csv(locale_rows, file.path(out_dir, "v2_locale_invariance.csv"),
          row.names = FALSE, na = "")

# Complexity evidence is counter-based; no large timing benchmark is run.
complexity_rows <- list()
for (shape in c("balanced", "random", "pectinate")) {
  for (n in c(32L, 64L, 128L, 256L, 512L, 1024L)) {
    tree <- switch(
      shape,
      balanced = set_lengths(ape::stree(n, type = "balanced"), offset = n),
      random = make_random_tree(n),
      pectinate = set_lengths(make_pectinate(n), offset = n * 2L)
    )
    got <- run_canonical(tree)
    c <- got$canonical
    values <- if (is.null(c)) list(child_rank_tuple_entries = NA_real_,
                                   descriptor_entries = NA_real_,
                                   peak_internal_object_proxy = NA_real_) else
      lapply(c$complexity, function(x) {
        first_numeric(x)
      })
    payload <- unlist(values, use.names = FALSE)
    finite_payload <- all(is.finite(payload))
    complexity_rows[[length(complexity_rows) + 1L]] <- data.frame(
      shape = shape, n_tip = n,
      canonical_status = got$call$status,
      child_rank_tuple_entries = values$child_rank_tuple_entries,
      descriptor_entries = values$descriptor_entries,
      peak_internal_object_proxy = values$peak_internal_object_proxy,
      payload_max = if (finite_payload) max(payload) else NA_real_,
      payload_over_n2 = if (finite_payload) max(payload) / (n * n) else NA_real_,
      linear_bound_128n = if (finite_payload) max(payload) <= 128 * n else NA,
      stringsAsFactors = FALSE
    )
  }
}
complexity_rows <- do.call(rbind, complexity_rows)
write.csv(complexity_rows, file.path(out_dir, "v2_complexity.csv"),
          row.names = FALSE, na = "")

# Protected snapshots and contexts.  The helper functions below mutate only
# serialized local copies, never the source fixture or a production context.
find_path <- function(x, wanted, prefix = character(), depth = 0L) {
  if (depth > 4L || !is.list(x)) return(NULL)
  hit <- intersect(wanted, names(x))
  if (length(hit)) return(c(prefix, hit[[1L]]))
  for (nm in names(x)) {
    value <- x[[nm]]
    if (is.list(value) && !is.environment(value)) {
      result <- find_path(value, wanted, c(prefix, nm), depth + 1L)
      if (!is.null(result)) return(result)
    }
  }
  NULL
}

find_phylo_path <- function(x, prefix = character(), depth = 0L) {
  if (depth > 4L) return(NULL)
  if (inherits(x, "phylo") && is_tree_like(x)) return(prefix)
  if (!is.list(x)) return(NULL)
  for (nm in names(x)) {
    value <- x[[nm]]
    if (is.environment(value)) next
    result <- find_phylo_path(value, c(prefix, nm), depth + 1L)
    if (!is.null(result)) return(result)
  }
  NULL
}

set_path <- function(x, path, value) {
  if (!length(path)) return(value)
  nm <- path[[1L]]
  x[[nm]] <- set_path(x[[nm]], path[-1L], value)
  x
}

get_path <- function(x, path) {
  if (!length(path)) return(x)
  if (!is.list(x) || !path[[1L]] %in% names(x)) return(NULL)
  get_path(x[[path[[1L]]]], path[-1L])
}

mutate_context_field <- function(context, field) {
  path <- find_phylo_path(context)
  if (is.null(path)) return(NULL)
  tree <- get_path(context, path)
  if (identical(field, "tip.label")) {
    if (identical(tree$tip.label, c("a", "b\rc", "a\rb", "c"))) {
      tree$tip.label <- c("a\rb", "c", "a", "b\rc")
    } else {
      tree$tip.label[[1L]] <- paste0(tree$tip.label[[1L]], "_mutated")
    }
  } else if (identical(field, "edge")) {
    tree$edge[[1L]] <- if (tree$edge[[1L]] == 1L) 2L else tree$edge[[1L]] + 1L
  } else if (identical(field, "edge.length")) {
    tree$edge.length[[1L]] <- tree$edge.length[[1L]] + 0.125
  } else if (identical(field, "Nnode")) {
    tree$Nnode <- as.integer(tree$Nnode) + 1L
  }
  set_path(context, path, tree)
}

context_schema_path <- function(context) {
  find_path(context, c(
    "context_schema_version", "schema_version", "fastphylosig_schema_version",
    "schema"
  ))
}

mutate_old_schema <- function(context) {
  path <- context_schema_path(context)
  if (is.null(path)) return(NULL)
  marker <- get_path(context, path)
  old <- if (is.numeric(marker)) 1L else "1"
  set_path(context, path, old)
}

base_context_call <- capture(function() call_v2("v2_prepare_context", clone(collision_base)))
base_context <- base_context_call$value
valid_context_call <- if (identical(base_context_call$status, "ok")) {
  capture(function() call_v2("v2_validate_context", clone(base_context)))
} else NULL

mutation_rows <- list()
for (field in c("tip.label", "edge", "edge.length", "Nnode")) {
  mutated <- if (is.null(base_context)) NULL else
    mutate_context_field(clone(base_context), field)
  snapshot_before <- capture(function() call_v2(
    "v2_protected_snapshot", clone(collision_base)
  ))
  mutated_tree <- if (is.null(mutated)) NULL else {
    path <- find_phylo_path(mutated)
    if (is.null(path)) NULL else get_path(mutated, path)
  }
  snapshot_after <- if (is.null(mutated_tree)) NULL else capture(function() {
    call_v2("v2_protected_snapshot", clone(mutated_tree))
  })
  validation <- if (is.null(mutated)) NULL else capture(function() {
    call_v2("v2_validate_context", mutated)
  })
  mutated_fingerprint <- if (is.null(mutated_tree)) NULL else capture(function() {
    call_v2("v2_fingerprint", clone(mutated_tree))
  })
  rejected <- !is.null(validation) && (
    identical(validation$status, "error") ||
      (is.logical(validation$value) && length(validation$value) == 1L &&
         isFALSE(validation$value))
  )
  snapshot_changed <- !is.null(snapshot_after) &&
    !safe_identical(snapshot_before$value, snapshot_after$value)
  mutation_rows[[length(mutation_rows) + 1L]] <- data.frame(
    field = field,
    context_prepared = identical(base_context_call$status, "ok"),
    mutation_applied = !is.null(mutated),
    snapshot_changed = snapshot_changed,
    fingerprint_changed = !is.null(mutated_fingerprint) &&
      identical(mutated_fingerprint$status, "ok") &&
      !safe_identical(mutated_fingerprint$value, base_context$fingerprint),
    validation_rejected = rejected,
    validation_status = if (is.null(validation)) "NOT_RUN_INTERFACE" else
      validation$status,
    validation_error = if (is.null(validation)) "" else validation$error_text,
    stringsAsFactors = FALSE
  )
}

old_schema_context <- if (is.null(base_context)) NULL else
  mutate_old_schema(clone(base_context))
old_schema_validation <- if (is.null(old_schema_context)) NULL else capture(function() {
  call_v2("v2_validate_context", old_schema_context)
})
mutation_rows[[length(mutation_rows) + 1L]] <- data.frame(
  field = "old_schema",
  context_prepared = identical(base_context_call$status, "ok"),
  mutation_applied = !is.null(old_schema_context),
  snapshot_changed = NA,
  fingerprint_changed = NA,
  validation_rejected = !is.null(old_schema_validation) && (
    identical(old_schema_validation$status, "error") ||
      (is.logical(old_schema_validation$value) &&
         length(old_schema_validation$value) == 1L &&
         isFALSE(old_schema_validation$value))
  ),
  validation_status = if (is.null(old_schema_validation)) "NOT_RUN_INTERFACE" else
    old_schema_validation$status,
  validation_error = if (is.null(old_schema_validation)) "" else
    old_schema_validation$error_text,
  stringsAsFactors = FALSE
)
mutation_rows <- do.call(rbind, mutation_rows)
write.csv(mutation_rows, file.path(out_dir, "v2_mutation_safety.csv"),
          row.names = FALSE, na = "")

# Failure semantics: invalid biology must be rejected or remain visibly
# invalid; the prototype may choose an error class/message that differs from
# production, but it may not silently return a valid canonical tree.
failure_rows <- list()
for (nm in names(failure_fixtures)) {
  tree <- failure_fixtures[[nm]]
  can <- capture(function() call_v2("v2_canonicalize", clone(tree)))
  prep <- capture(function() call_v2("v2_prepare_context", clone(tree)))
  can_tree <- if (identical(can$status, "ok")) as_tree_candidate(can$value) else NULL
  prep_tree <- if (identical(prep$status, "ok")) find_phylo_path(prep$value) else NULL
  failure_rows[[length(failure_rows) + 1L]] <- data.frame(
    fixture = nm,
    canonical_status = can$status,
    canonical_error_class = can$error_classes,
    canonical_error = can$error_text,
    canonical_valid_tree_returned = !is.null(can_tree),
    prepare_status = prep$status,
    prepare_error_class = prep$error_classes,
    prepare_error = prep$error_text,
    prepare_tree_returned = !is.null(prep_tree),
    rejected_without_repair = identical(can$status, "error") &&
      identical(prep$status, "error"),
    stringsAsFactors = FALSE
  )
}
failure_rows <- do.call(rbind, failure_rows)
write.csv(failure_rows, file.path(out_dir, "v2_failure_semantics.csv"),
          row.names = FALSE, na = "")

# Frozen estimator oracle.  The V2 prototype is never used by the production
# entry points here; only its canonical tree output is passed to production.
prod_ns <- asNamespace("fastphylosig")
prod_fun <- function(name) get(name, envir = prod_ns, inherits = FALSE)
trait_for <- function(tree) {
  labels <- as.character(tree$tip.label)
  list(
    continuous = stats::setNames(seq_along(labels) / 7, labels),
    binary = stats::setNames(as.numeric(seq_along(labels) %% 2L), labels),
    categorical = stats::setNames(
      factor(ifelse(seq_along(labels) %% 2L, "A", "B"), levels = c("A", "B")),
      labels
    )
  )
}

fixed_permutations <- function(n) {
  if (n < 2L) return(matrix(seq_len(n), nrow = 1L))
  rbind(seq_len(n), rev(seq_len(n)),
        c(2L:n, 1L), c(n, seq_len(n - 1L)))
}

result_field <- function(value, method, field) {
  if (is.null(value)) return(NULL)
  if (method == "K" && field == "estimate" && is.numeric(value) &&
      is.null(dim(value))) return(value)
  if (!is.list(value)) return(NULL)
  aliases <- switch(
    method,
    K = c(estimate = "K", P = "P", MCSE = "MCSE_P",
          status = "status", successful = "nsim_successful"),
    lambda = c(estimate = "lambda", logLik = "logL", LR = "LR",
               status = "status"),
    D = c(estimate = "DEstimate", P_random = "Pval1", P_Brownian = "Pval0",
          status = "status", successful_random = "nsim_successful_random",
          successful_brownian = "nsim_successful_brownian"),
    Delta = c(estimate = "delta", P = "P", MCSE = "MCSE_Delta",
              status = "status", successful = "successful_simulations"),
    ACE = c(logLik = "loglik", lik_anc = "lik.anc", ace = "ace",
            rates = "rates", se = "se", convergence = "convergence")
  )
  nm <- unname(aliases[[field]])
  if (is.null(nm) || !nm %in% names(value)) NULL else value[[nm]]
}

estimator_specs <- list(
  K = function(tree, trait, n) {
    prod_fun("fast_k")(tree, trait$continuous, test = TRUE, nsim = 4L,
                         permutations = fixed_permutations(n), return_sim = FALSE,
                         verbose = FALSE, progress = FALSE, ncores = 1L)
  },
  lambda = function(tree, trait, n) {
    prod_fun("fast_lambda")(tree, trait$continuous, test = FALSE,
                              lambda_profile = FALSE, verbose = FALSE,
                              progress = FALSE, ncores = 1L)
  },
  D = function(tree, trait, n) {
    random_states <- cbind(
      as.numeric(seq_len(n) %% 2L), as.numeric((seq_len(n) + 1L) %% 2L),
      as.numeric((seq_len(n) + 2L) %% 2L),
      rev(as.numeric(seq_len(n) %% 2L))
    )
    brownian_states <- cbind(
      seq_len(n) / 11, rev(seq_len(n) / 13),
      (seq_len(n) + 2) / 17, rep(0.25, n)
    )
    prod_fun("fast_d")(
      tree, trait$binary, test = TRUE, nsim = 4L,
      random_states = random_states, brownian_states = brownian_states,
      return_sim = FALSE, verbose = FALSE, progress = FALSE, ncores = 1L
    )
  },
  Delta = function(tree, trait, n) {
    set.seed(20260906)
    prod_fun("fast_delta")(
      tree, trait$categorical, test = FALSE, mcmc_sim = 40L, thin = 2L,
      burn = 10L, lambda0 = 0.1, proposal_sd = 0.5, model = "ER",
      entropy = "LSE", verbose = FALSE, progress = FALSE, ncores = 1L
    )
  },
  ACE = function(tree, trait, n) {
    prod_fun("fast_ace")(
      trait$categorical, phy = tree, type = "discrete", method = "ML",
      model = "ER", CI = FALSE, marginal = FALSE, progress = FALSE
    )
  }
)

estimator_metrics <- list(
  K = c("estimate", "P", "MCSE", "status", "successful"),
  lambda = c("estimate", "logLik", "LR", "status"),
  D = c("estimate", "P_random", "P_Brownian", "status",
        "successful_random", "successful_brownian"),
  Delta = c("estimate", "P", "MCSE", "status", "successful"),
  ACE = c("logLik", "lik_anc", "ace", "rates", "se", "convergence")
)

estimator_rows <- list()
parity_groups <- list(
  random = c("random", "edge_shuffled", "internal_renumbered",
             "alternate_root_numbering", "renumbered_and_shuffled"),
  delimiter = c("delimiter", "delimiter_edge_shuffled",
                "delimiter_internal_renumbered", "delimiter_renumbered_shuffled")
)
for (method in names(estimator_specs)) {
  for (fixture in c("random", "delimiter")) {
    reference_tree <- valid_fixtures[[fixture]]
    ref_trait <- trait_for(reference_tree)
    ref_call <- capture(function() estimator_specs[[method]](
      reference_tree, ref_trait, length(reference_tree$tip.label)
    ))
    for (route in c("frozen_raw", "v2_canonical")) {
      ref_value <- ref_call$value
      for (target in parity_groups[[fixture]]) {
        target_tree <- valid_fixtures[[target]]
        target_can <- canonical_runs[[target]]$canonical
        target_value <- if (route == "frozen_raw") {
          capture(function() estimator_specs[[method]](
            target_tree, ref_trait, length(target_tree$tip.label)
          ))
        } else if (is.null(target_can)) {
          list(value = NULL, status = "error", error_text = "V2 canonical output missing")
        } else {
          capture(function() estimator_specs[[method]](
            target_can$tree, ref_trait,
            length(target_can$tree$tip.label)
          ))
        }
        for (metric in estimator_metrics[[method]]) {
          lhs <- result_field(ref_value, method, metric)
          rhs <- result_field(target_value$value, method, metric)
          estimator_rows[[length(estimator_rows) + 1L]] <- data.frame(
            method = method, reference = fixture, target = target, route = route,
            metric = metric,
            reference_status = ref_call$status,
            target_status = target_value$status,
            exact_identical = safe_identical(lhs, rhs),
            parity_within_existing_tolerance = parity_equal(lhs, rhs),
            absolute_difference = numeric_difference(lhs, rhs),
            relative_difference = relative_difference(lhs, rhs),
            target_warning_classes = target_value$warning_classes %||% "",
            target_error = target_value$error_text %||% "",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
estimator_rows <- do.call(rbind, estimator_rows)
write.csv(estimator_rows, file.path(out_dir, "v2_estimator_parity.csv"),
          row.names = FALSE, na = "")

# `%||%` is kept local to the audit script so missing metadata never aborts a
# row after the actual estimator call has completed.

gate_from_logical <- function(values, missing = FALSE) {
  if (isTRUE(missing) || !length(values) || anyNA(values)) return("FAIL")
  if (all(values)) "PASS" else "FAIL"
}

renumber_rows <- invariance_rows[
  invariance_rows$group == "random_representations", , drop = FALSE
]
all_representation_rows <- invariance_rows[
  invariance_rows$group %in% c("random_representations",
                               "delimiter_representations"), , drop = FALSE
]
renumber_gate <- gate_from_logical(
  renumber_rows$canonical_fields_identical & renumber_rows$fingerprint_identical
)
edge_gate <- gate_from_logical(
  all_representation_rows$canonical_fields_identical &
    all_representation_rows$fingerprint_identical
)
delimiter_gate <- gate_from_logical(
  invariance_rows$canonical_fields_identical[
    invariance_rows$group == "delimiter_representations"
  ] & invariance_rows$fingerprint_identical[
    invariance_rows$group == "delimiter_representations"
  ]
)
delimiter_mutation <- mutation_rows[
  mutation_rows$field == "tip.label", , drop = FALSE
]
if (!nrow(delimiter_mutation) ||
    !isTRUE(delimiter_mutation$fingerprint_changed[[1L]])) {
  delimiter_gate <- "FAIL"
}
locale_supported <- locale_rows$status == "PASS"
locale_gate <- if (!any(locale_supported)) "FAIL" else gate_from_logical(
  locale_rows$canonical_fields_identical_to_C[locale_supported] &
    locale_rows$fingerprint_identical_to_C[locale_supported]
)
branch_gate <- gate_from_logical(
  invariance_rows$edge_length_identical &
    invariance_rows$canonical_fields_identical
)
mutation_gate <- gate_from_logical(mutation_rows$validation_rejected,
                                   missing = any(!mutation_rows$mutation_applied))
schema_row <- mutation_rows[mutation_rows$field == "old_schema", , drop = FALSE]
schema_gate <- if (!nrow(schema_row)) "FAIL" else
  gate_from_logical(schema_row$validation_rejected,
                    missing = !isTRUE(schema_row$mutation_applied[[1L]]))
stale_gate <- if (identical(mutation_gate, "PASS") &&
                   identical(schema_gate, "PASS")) {
  "IMPOSSIBLE_BY_VALIDATION"
} else "POSSIBLE_OR_NOT_ESTABLISHED"
failure_gate <- gate_from_logical(failure_rows$rejected_without_repair)
complexity_gate <- if (!nrow(complexity_rows) ||
                       anyNA(complexity_rows$linear_bound_128n) ||
                       any(!as.logical(complexity_rows$linear_bound_128n)) ||
                       any(!is.finite(complexity_rows$payload_over_n2))) {
  "FAIL"
} else "PASS"
immutability_gate <- gate_from_logical(vapply(
  canonical_runs, function(x) isTRUE(x$input_immutable), logical(1L)
))
stored_label_gate <- gate_from_logical(vapply(
  canonical_runs, function(x) isTRUE(x$stored_tip_label_identical), logical(1L)
))
v2_estimator_rows <- estimator_rows[
  estimator_rows$route == "v2_canonical", , drop = FALSE
]
estimator_gate <- if (!nrow(v2_estimator_rows)) "FAIL" else
  gate_from_logical(v2_estimator_rows$parity_within_existing_tolerance)

summary <- data.frame(
  gate = c(
    "V2_RENUMBERING_INVARIANCE", "V2_EDGE_ORDER_INVARIANCE",
    "V2_LOCALE_INVARIANCE", "V2_DELIMITER_SAFETY",
    "V2_BRANCH_ASSOCIATION", "V2_MUTATION_REJECTION",
    "STALE_CACHE_PROPAGATION", "V2_ESTIMATOR_PARITY",
    "V2_FAILURE_SEMANTICS", "V2_COMPLEXITY",
    "V2_INPUT_IMMUTABILITY", "V2_STORED_LABEL_PRESERVATION"
  ),
  status = c(
    renumber_gate, edge_gate, locale_gate, delimiter_gate, branch_gate,
    mutation_gate, stale_gate, estimator_gate, failure_gate, complexity_gate,
    immutability_gate, stored_label_gate
  ),
  stringsAsFactors = FALSE
)
prototype_gate <- if (all(summary$status %in% c("PASS", "IMPOSSIBLE_BY_VALIDATION"))) {
  "PASS"
} else "FAIL"
summary <- rbind(
  summary,
  data.frame(gate = "V2_PROTOTYPE", status = prototype_gate,
             stringsAsFactors = FALSE)
)
write.csv(summary, file.path(out_dir, "v2_gate_summary.csv"),
          row.names = FALSE, na = "")

provenance <- data.frame(
  field = c(
    "source_commit", "package_load", "package_version", "R_version",
    "platform", "LC_COLLATE_at_start", "prototype_file", "audit_scope"
  ),
  value = c(
    source_commit, package_load, as.character(utils::packageVersion("fastphylosig")),
    R.version.string, R.version$platform, normal_locale,
    normalizePath(prototype_file, winslash = "/"),
    "V2-A exact contract/prototype audit; no production integration or large benchmark"
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "v2_provenance.csv"),
          row.names = FALSE, na = "")

message("V2-A audit completed: ", prototype_gate)
message("Evidence directory: ", out_dir)
