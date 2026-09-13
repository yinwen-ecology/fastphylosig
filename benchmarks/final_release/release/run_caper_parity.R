# Optional final-release observed-D parity check against caper.
#
# Usage:
#   Rscript --vanilla run_caper_parity.R <fresh-library> <output>
#
# The check is intentionally small and deterministic.  It compares only the
# observed contrast sum reported by caper::phylo.d() and fast_d().  The one
# fixed null column passed to fast_d() is required by its call contract, but
# null summaries and permutation results are deliberately out of scope here.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_caper_parity.R <fresh-library> <output>", call. = FALSE)
}

library_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
if (!nzchar(output)) stop("The output path must not be empty.", call. = FALSE)
output_dir <- dirname(output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_dir)) {
  stop("Cannot create the output directory.", call. = FALSE)
}

write_report <- function(lines) {
  writeLines(lines, output, useBytes = TRUE)
}

package_path <- find.package(
  "fastphylosig", lib.loc = library_path, quiet = TRUE
)
if (length(package_path) != 1L) {
  write_report(c(
    "fastphylosig 0.2.0 final release caper observed-D parity",
    paste0("LIBRARY_PATH=", library_path),
    "CAPER_REFERENCE_PARITY=ERROR",
    "CAPER_REFERENCE_DETAIL=fastphylosig is not installed in the supplied fresh library",
    "COMPARISON=observed D only; comparison was not attempted",
    "TOLERANCE_ABS=1e-10",
    "TOLERANCE_REL=1e-10"
  ))
  quit(save = "no", status = 1L)
}
package_path <- normalizePath(package_path, winslash = "/", mustWork = TRUE)
.libPaths(unique(c(library_path, .libPaths())))
suppressPackageStartupMessages(
  library("fastphylosig", character.only = TRUE, lib.loc = library_path)
)

loaded_path <- tryCatch(
  normalizePath(
    getNamespaceInfo(asNamespace("fastphylosig"), "path"),
    winslash = "/", mustWork = TRUE
  ),
  error = function(e) ""
)
if (!identical(loaded_path, package_path)) {
  write_report(c(
    "fastphylosig 0.2.0 final release caper observed-D parity",
    paste0("LIBRARY_PATH=", library_path),
    paste0("PACKAGE_PATH=", package_path),
    paste0("LOADED_PATH=", loaded_path),
    "CAPER_REFERENCE_PARITY=ERROR",
    "CAPER_REFERENCE_DETAIL=fastphylosig did not load from the supplied fresh library",
    "COMPARISON=observed D only; comparison was not attempted",
    "TOLERANCE_ABS=1e-10",
    "TOLERANCE_REL=1e-10"
  ))
  quit(save = "no", status = 1L)
}

package_version <- as.character(utils::packageVersion("fastphylosig"))
tolerance_abs <- 1e-10
tolerance_rel <- 1e-10

report_header <- c(
  "fastphylosig 0.2.0 final release caper observed-D parity",
  paste0("LIBRARY_PATH=", library_path),
  paste0("PACKAGE_PATH=", loaded_path),
  paste0("FASTPHYLOSIG_VERSION=", package_version),
  paste0("R_VERSION=", R.version.string),
  paste0("TOLERANCE_ABS=", format(tolerance_abs, scientific = TRUE)),
  paste0("TOLERANCE_REL=", format(tolerance_rel, scientific = TRUE)),
  "TOLERANCE_REASON=1e-10 follows the repository observed-contrast checks and allows harmless floating-point traversal differences",
  "COMPARISON=observed D only; P values, null means, permutations, and timing are not compared",
  "FIXTURE=inline fixed 4-tip ultrametric tree with named binary states",
  "CAPER_PERMUT=1",
  "CAPER_PERMUT_REASON=the tested caper interface rejects permut=0; one required draw is not used in the observed-D comparison",
  "CAPER_SEED=20260913",
  "FAST_NSIM=1",
  "FAST_NULL_INPUT=fixed one-column matrices supplied only to satisfy fast_d call contract"
)

if (!requireNamespace("caper", quietly = TRUE)) {
  write_report(c(
    report_header,
    "CAPER_VERSION=NOT_INSTALLED",
    "CAPER_REFERENCE_PARITY=NOT_RUN_DEPENDENCY",
    "CAPER_REFERENCE_DETAIL=caper is not installed; observed-D comparison was not attempted"
  ))
  quit(save = "no", status = 0L)
}

capture_call <- function(fun) {
  warnings <- character()
  error_message <- NULL
  value <- withCallingHandlers(
    tryCatch(fun(), error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, error = error_message, warnings = warnings)
}

extract_observed <- function(value) {
  if (!is.list(value) || is.null(value$Parameters)) return(NA_real_)
  parameters <- value$Parameters
  observed <- parameters$Observed
  if (is.null(observed)) observed <- parameters$observed
  if (length(observed) != 1L) return(NA_real_)
  observed <- tryCatch(as.numeric(observed[[1L]]), error = function(e) NA_real_)
  if (length(observed) != 1L) NA_real_ else observed
}

finite_scalar <- function(value) {
  is.numeric(value) && length(value) == 1L && is.finite(value)
}

within_tolerance <- function(x, y) {
  finite_scalar(x) && finite_scalar(y) &&
    abs(x - y) <= tolerance_abs + tolerance_rel * max(1, abs(x), abs(y))
}

tree <- ape::read.tree(text =
  "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);"
)
trait <- stats::setNames(c(0L, 1L, 0L, 1L), tree$tip.label)

# These fixed columns avoid RNG and keep the audit focused on the observed
# contrast.  They are not used in the parity decision.
random_states <- matrix(
  c(1, 0, 1, 0), nrow = length(tree$tip.label), ncol = 1L,
  dimnames = list(tree$tip.label, "draw1")
)
brownian_states <- matrix(
  c(0, 0, 1, 1), nrow = length(tree$tip.label), ncol = 1L,
  dimnames = list(tree$tip.label, "draw1")
)

reference_call <- capture_call(function() {
  set.seed(20260913L)
  data <- data.frame(
    sp = tree$tip.label,
    trait = as.numeric(trait),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  comparative <- caper::comparative.data(
    tree,
    data,
    names.col = sp,
    warn.dropped = FALSE
  )
  caper::phylo.d(
    comparative,
    binvar = trait,
    permut = 1L
  )
})

fast_call <- capture_call(function() {
  fastphylosig::fast_d(
    tree,
    trait,
    test = FALSE,
    nsim = 1L,
    random_states = random_states,
    brownian_states = brownian_states,
    return_sim = FALSE,
    verbose = FALSE,
    progress = FALSE,
    ncores = 1L
  )
})

reference_observed <- extract_observed(reference_call$value)
fast_observed <- extract_observed(fast_call$value)
detail <- character()
status <- "PASS"

if (!is.null(reference_call$error) || !is.null(fast_call$error)) {
  status <- "FAIL"
  detail <- c(
    "one or both observed-D calls failed",
    paste0("caper_error=", if (is.null(reference_call$error)) "none" else reference_call$error),
    paste0("fast_error=", if (is.null(fast_call$error)) "none" else fast_call$error)
  )
} else if (length(fast_call$warnings)) {
  status <- "FAIL"
  detail <- c(
    "fast_d emitted an unexpected warning",
    paste0("fast_warning=", paste(fast_call$warnings, collapse = " | "))
  )
} else if (!finite_scalar(reference_observed) || !finite_scalar(fast_observed)) {
  status <- "FAIL"
  detail <- c(
    "one or both calls returned a missing or non-finite Parameters$Observed",
    paste0("caper_observed=", format(reference_observed)),
    paste0("fast_observed=", format(fast_observed))
  )
} else {
  observed_difference <- fast_observed - reference_observed
  if (!within_tolerance(fast_observed, reference_observed)) {
    status <- "FAIL"
    detail <- c(
      "observed D contrast sums differ beyond the declared tolerance",
      paste0("caper_observed=", format(reference_observed, digits = 17)),
      paste0("fast_observed=", format(fast_observed, digits = 17)),
      paste0("difference_fast_minus_caper=", format(observed_difference, digits = 17))
    )
  } else {
    detail <- c(
      "observed D contrast sums agree within the declared tolerance",
      paste0("caper_observed=", format(reference_observed, digits = 17)),
      paste0("fast_observed=", format(fast_observed, digits = 17)),
      paste0("difference_fast_minus_caper=", format(observed_difference, digits = 17))
    )
  }
}

report <- c(
  report_header,
  paste0("CAPER_VERSION=", as.character(utils::packageVersion("caper"))),
  paste0("CAPER_WARNINGS=", paste(reference_call$warnings, collapse = " | ")),
  paste0("FAST_WARNINGS=", paste(fast_call$warnings, collapse = " | ")),
  paste0("CAPER_REFERENCE_PARITY=", status),
  paste0("CAPER_REFERENCE_DETAIL=", paste(detail, collapse = "; "))
)
write_report(report)

if (identical(status, "FAIL")) quit(save = "no", status = 1L)
