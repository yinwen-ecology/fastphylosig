#!/usr/bin/env Rscript

# Stage 2B1 raw helper call-count audit.  This script changes only function
# bindings in an isolated child R process, restores them before exit, and
# never edits package source.  It is intended to be run after the paired
# timing harness, but has no dependency on timing output.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste("usage: raw_helper_call_audit.R BEFORE_LIB AFTER_LIB OUTPUT_DIR"),
    call. = FALSE
  )
}

before_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
after_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
options(stringsAsFactors = FALSE)

if (!requireNamespace("ape", quietly = TRUE) ||
    !requireNamespace("callr", quietly = TRUE)) {
  stop("the call-count audit requires ape and callr.", call. = FALSE)
}

make_fixture <- function(n, seed = 20260904L) {
  set.seed(as.integer(seed + n))
  ape::reorder.phylo(ape::rtree(n), order = "postorder")
}

helpers <- c(
  ".inspect_tree_core", ".safe_canonicalize_core",
  ".canonical_tree_signature", ".tree_fingerprint",
  ".prepare_tree_core", ".prepare_tree_subset", ".prepared_tree_subset",
  ".validated_context", ".validate_prepared_context",
  ".analysis_data_table", ".match_tree_data_core", ".analysis_na_groups",
  ".as_trait_matrix"
)

audit_child <- function(lib, route, tree, trait) {
  callr::r(
    function(lib, route, tree, trait, helpers) {
      library("fastphylosig", lib.loc = lib, character.only = TRUE)
      ns <- asNamespace("fastphylosig")
      available <- vapply(
        helpers, exists, logical(1), envir = ns, inherits = FALSE
      )
      active <- helpers[available]
      context <- if (identical(route, "prepared")) {
        # Preparation is deliberately outside the counted prepared public
        # call, matching the prepared timing boundary.
        fastphylosig::prepare_tree(tree)
      } else NULL
      originals <- lapply(active, get, envir = ns, inherits = FALSE)
      names(originals) <- active
      locked <- vapply(active, bindingIsLocked, logical(1), env = ns)
      counts <- setNames(as.list(rep(0L, length(active))), active)
      wrap <- function(name, fun) {
        force(name)
        force(fun)
        function(...) {
          counts[[name]] <<- counts[[name]] + 1L
          fun(...)
        }
      }
      for (name in active) {
        if (locked[[name]]) unlockBinding(name, ns)
        assign(name, wrap(name, originals[[name]]), envir = ns)
      }
      on.exit({
        for (name in active) {
          if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
          assign(name, originals[[name]], envir = ns)
          if (locked[[name]]) lockBinding(name, ns)
        }
      }, add = TRUE)

      result <- if (identical(route, "raw")) {
        fastphylosig::fast_k(
          tree, trait, test = FALSE, verbose = FALSE, progress = FALSE
        )
      } else {
        fastphylosig::fast_k(
          context, trait, test = FALSE, verbose = FALSE, progress = FALSE
        )
      }
      invisible(result)
      data.frame(
        helper = helpers,
        calls = vapply(seq_along(helpers), function(i) {
          if (!available[[i]]) NA_integer_ else
            as.integer(counts[[helpers[[i]]]])
        }, integer(1)),
        available = available,
        stringsAsFactors = FALSE
      )
    },
    args = list(lib = lib, route = route, tree = tree, trait = trait,
                helpers = helpers),
    libpath = c(lib, .libPaths()), spinner = FALSE
  )
}

audit_ns <- c(500L, 5000L)
rows <- list()
row_id <- 0L
message("[stage2b1] helper call-count audit: starting")
for (version in c("before", "after")) {
  lib <- if (identical(version, "before")) before_lib else after_lib
  for (n in audit_ns) {
    tree <- make_fixture(n)
    set.seed(as.integer(97531L + n))
    trait <- stats::setNames(stats::rnorm(n), tree$tip.label)
    for (route in c("raw", "prepared")) {
      counted <- audit_child(lib, route, tree, trait)
      for (i in seq_len(nrow(counted))) {
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          version = version, n = n, route = route,
          helper = counted$helper[[i]], calls = counted$calls[[i]],
          available = counted$available[[i]],
          stringsAsFactors = FALSE
        )
      }
    }
    rm(tree, trait)
    gc(FALSE)
  }
}
counts <- do.call(rbind, rows)
utils::write.csv(counts,
                 file.path(out_dir, "stage2b1_raw_helper_call_counts.csv"),
                 row.names = FALSE)

# Stage 2A's authoritative raw counts are included as a comparison contract;
# the actual before values are still read from the output above by reviewers.
contract <- data.frame(
  helper = helpers,
  stage2a_raw_calls = c(3L, 2L, 4L, 2L, NA, 1L, 1L, NA, NA,
                        1L, 1L, 1L, 1L),
  stage2b1_target_calls = c(1L, 1L, 2L, 1L, 1L, 1L, 1L, NA, NA,
                            1L, 1L, 1L, 1L),
  remaining_call_reason = c(
    "one post-canonicalization readiness inspection",
    "one representation-safe canonicalization",
    "before/after canonical signatures inside that canonicalization",
    "one fingerprint for the newly created context/integrity contract",
    "one private compilation boundary",
    "one full subset cache seed",
    "one retained-mask lookup in common analysis",
    "private validation token only; not an external boundary",
    "prepared boundary remains independently checked; early-return token is allowed",
    "one shared trait-table conversion",
    "one shared species matching pass",
    "one shared NA-mask grouping pass",
    "one shared numeric matrix coercion"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(contract,
                 file.path(out_dir, "stage2b1_helper_call_contract.csv"),
                 row.names = FALSE)

provenance <- data.frame(
  key = c("before_commit", "after_commit", "R_version", "R_platform",
          "before_library", "after_library", "audit_scope"),
  value = c(
    Sys.getenv("FASTPHYLOSIG_BEFORE_COMMIT", unset = "<unspecified>"),
    Sys.getenv("FASTPHYLOSIG_AFTER_COMMIT", unset = "<unspecified>"),
    R.version.string, R.version$platform, before_lib, after_lib,
    "one raw and one prepared fast_k call at n=500 and n=5000; package bindings restored in child"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(provenance,
                 file.path(out_dir, "stage2b1_call_audit_provenance.csv"),
                 row.names = FALSE)
manifest <- data.frame(
  file = list.files(out_dir, pattern = "^stage2b1_", full.names = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(out_dir, "stage2b1_call_audit_manifest.csv"),
                 row.names = FALSE)
print(counts)
