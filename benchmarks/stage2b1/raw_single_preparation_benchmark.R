#!/usr/bin/env Rscript

# Stage 2B1 raw single-preparation boundary benchmark.
#
# This is an audit harness only.  It compares two installed package libraries
# and never edits package source.  A timed raw call includes the package's
# preparation work.  A timed prepared call constructs its context before the
# timer and therefore measures the prepared public boundary.  prepare_tree()
# is measured independently.  callr keeps the two package namespaces isolated;
# process startup and serialization are outside the child timer.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "usage: raw_single_preparation_benchmark.R BEFORE_LIB AFTER_LIB",
      "OUTPUT_DIR [quick]"
    ),
    call. = FALSE
  )
}

before_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
after_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
quick_mode <- length(args) >= 4L && identical(args[[4L]], "quick")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
options(stringsAsFactors = FALSE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("the benchmark requires the ape package.", call. = FALSE)
}
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("the benchmark requires the callr package.", call. = FALSE)
}

parse_n_values <- function() {
  supplied <- Sys.getenv("FASTPHYLOSIG_STAGE2B1_NS", unset = "")
  if (!nzchar(supplied)) {
    return(if (quick_mode) c(500L, 1000L, 2000L) else
      c(500L, 1000L, 2000L, 5000L, 10000L, 20000L))
  }
  value <- suppressWarnings(as.integer(strsplit(supplied, ",", fixed = TRUE)[[1L]]))
  if (!length(value) || anyNA(value) || any(value < 2L)) {
    stop("FASTPHYLOSIG_STAGE2B1_NS must be comma-separated positive integers.",
         call. = FALSE)
  }
  value
}

n_values <- parse_n_values()
repeat_override <- suppressWarnings(as.integer(
  Sys.getenv("FASTPHYLOSIG_STAGE2B1_REPEATS", unset = "")
))
if (length(repeat_override) && !is.na(repeat_override) && repeat_override > 0L) {
  repeats_for <- function(n) repeat_override
} else if (quick_mode) {
  repeats_for <- function(n) 2L
} else {
  # This is the Stage 2A protocol: 20 paired observations below n=10,000 and
  # 10 paired observations for the two largest fixtures.
  repeats_for <- function(n) if (n >= 10000L) 10L else 20L
}

inner_iterations <- function(n) {
  # One timed child observation reports the average of this many calls.  This
  # reduces the Windows timer quantum for short prepared calls while keeping
  # the high-n raw route at one full call per observation.
  if (quick_mode) return(1L)
  if (n <= 500L) return(20L)
  if (n <= 1000L) return(10L)
  if (n <= 2000L) return(5L)
  if (n <= 5000L) return(2L)
  1L
}

make_fixture <- function(n, seed = 20260904L) {
  set.seed(as.integer(seed + n))
  ape::reorder.phylo(ape::rtree(n), order = "postorder")
}

make_trait <- function(tree, seed = 97531L) {
  set.seed(as.integer(seed + ape::Ntip(tree)))
  stats::setNames(stats::rnorm(ape::Ntip(tree)), tree$tip.label)
}

package_version_at <- function(lib) {
  desc <- utils::packageDescription("fastphylosig", lib.loc = lib)
  if (is.null(desc$Version)) "<unknown>" else as.character(desc$Version)
}

call_timed <- function(lib, route, tree, trait, iterations) {
  callr::r(
    function(lib, route, tree, trait, iterations) {
      library("fastphylosig", lib.loc = lib, character.only = TRUE)
      iterations <- as.integer(iterations)
      if (identical(route, "prepared_fast_k")) {
        context <- fastphylosig::prepare_tree(tree)
        run <- function() fastphylosig::fast_k(
          context, trait, test = FALSE, verbose = FALSE, progress = FALSE
        )
      } else if (identical(route, "raw_fast_k")) {
        run <- function() fastphylosig::fast_k(
          tree, trait, test = FALSE, verbose = FALSE, progress = FALSE
        )
      } else if (identical(route, "prepare_tree")) {
        run <- function() fastphylosig::prepare_tree(tree)
      } else {
        stop("unknown benchmark route", call. = FALSE)
      }
      invisible(run())
      gc(FALSE)
      start <- proc.time()[["elapsed"]]
      value <- NULL
      for (i in seq_len(iterations)) value <- run()
      elapsed <- proc.time()[["elapsed"]] - start
      statistic <- if (identical(route, "raw_fast_k") ||
                       identical(route, "prepared_fast_k")) {
        if (is.data.frame(value) && "K_fast" %in% names(value)) {
          as.numeric(value$K_fast[[1L]])
        } else if (is.list(value) && !is.null(value$K)) {
          as.numeric(value$K[[1L]])
        } else {
          as.numeric(value[[1L]])
        }
      } else {
        NA_real_
      }
      list(elapsed = max(0, as.numeric(elapsed)) / iterations,
           statistic = statistic)
    },
    args = list(lib = lib, route = route, tree = tree, trait = trait,
                iterations = as.integer(iterations)),
    libpath = c(lib, .libPaths()), spinner = FALSE
  )
}

commit_for <- function(which) {
  key <- if (identical(which, "before")) {
    "FASTPHYLOSIG_BEFORE_COMMIT"
  } else {
    "FASTPHYLOSIG_AFTER_COMMIT"
  }
  Sys.getenv(key, unset = "<unspecified>")
}

provenance <- data.frame(
  key = c(
    "before_commit", "after_commit", "before_version", "after_version",
    "R_version", "R_platform", "R_arch", "R_home", "compiler",
    "OpenMP", "before_library", "after_library", "fixture",
    "formal_repeats", "warmup", "alternation", "timer_scope",
    "source_change_scope"
  ),
  value = c(
    commit_for("before"), commit_for("after"),
    package_version_at(before_lib), package_version_at(after_lib),
    R.version.string, R.version$platform, R.version$arch, R.home(),
    paste0(R.version$compiler, "; CXX=", Sys.getenv("CXX", unset = "<unset>")),
    "not inferable from R capabilities(); no build claim made",
    before_lib, after_lib,
    "ape::rtree(seed=20260904+n), postorder; trait seed=97531+n_tip",
    "20 paired repeats for n<10000; 10 paired repeats for n>=10000",
    "one route-specific warmup before formal pairs",
    "odd pairs before>after; even pairs after>before",
    "child process startup/serialization excluded from child timer; gc before timer",
    "audit harness only; no production source change"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(provenance,
                 file.path(out_dir, "stage2b1_provenance.csv"),
                 row.names = FALSE)

routes <- c("raw_fast_k", "prepare_tree", "prepared_fast_k")
timing_rows <- list()
row_id <- 0L

message("[stage2b1] paired timing run: starting")
for (n in n_values) {
  message("[stage2b1] n=", n)
  tree_a <- make_fixture(n)
  tree_b <- unserialize(serialize(tree_a, NULL))
  trait_a <- make_trait(tree_a)
  trait_b <- unserialize(serialize(trait_a, NULL))
  iterations <- inner_iterations(n)
  repeats <- repeats_for(n)

  for (route in routes) {
    # Equal route-specific warmup for both installed versions.  No warmup
    # value is included in the formal output.
    invisible(call_timed(before_lib, route, tree_a, trait_a, iterations))
    invisible(call_timed(after_lib, route, tree_b, trait_b, iterations))
    for (pair in seq_len(repeats)) {
      versions <- if (pair %% 2L) c("before", "after") else
        c("after", "before")
      for (position in seq_along(versions)) {
        version <- versions[[position]]
        lib <- if (identical(version, "before")) before_lib else after_lib
        tree <- if (identical(version, "before")) tree_a else tree_b
        trait <- if (identical(version, "before")) trait_a else trait_b
        measured <- call_timed(lib, route, tree, trait, iterations)
        row_id <- row_id + 1L
        timing_rows[[row_id]] <- data.frame(
          n = as.integer(n), route = route, pair = pair, position = position,
          version = version, iterations = iterations,
          elapsed_ms = 1000 * as.numeric(measured$elapsed),
          statistic = as.numeric(measured$statistic),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  rm(tree_a, tree_b, trait_a, trait_b)
  gc(FALSE)
}

timings <- do.call(rbind, timing_rows)
utils::write.csv(timings,
                 file.path(out_dir, "stage2b1_raw_preparation_timings.csv"),
                 row.names = FALSE)
message("[stage2b1] paired timing run: finished")

safe_max_abs <- function(x, y) {
  delta <- abs(as.numeric(x) - as.numeric(y))
  if (!length(delta) || all(is.na(delta))) NA_real_ else max(delta, na.rm = TRUE)
}

summarise_version <- function(z) {
  values <- function(route) z$elapsed_ms[z$route == route]
  raw <- values("raw_fast_k")
  prep <- values("prepare_tree")
  prepared <- values("prepared_fast_k")
  raw_med <- stats::median(raw)
  prep_med <- stats::median(prep)
  prepared_med <- stats::median(prepared)
  denominator <- prep_med + prepared_med
  residual <- raw_med - denominator
  data.frame(
    version = as.character(z$version[[1L]]), n = as.integer(z$n[[1L]]),
    repeats = length(raw),
    raw_median_ms = raw_med, raw_IQR_ms = stats::IQR(raw),
    prepare_median_ms = prep_med, prepare_IQR_ms = stats::IQR(prep),
    prepared_median_ms = prepared_med,
    prepared_IQR_ms = stats::IQR(prepared),
    duplicate_residual_ms = residual,
    duplicate_residual_pct = 100 * residual / raw_med,
    raw_overhead_ratio = raw_med / denominator,
    max_K_statistic_abs_diff = safe_max_abs(
      z$statistic[z$route == "raw_fast_k"],
      z$statistic[z$route == "prepared_fast_k"]
    ),
    stringsAsFactors = FALSE
  )
}

summary_rows <- lapply(split(
  timings, interaction(timings$version, timings$n, drop = TRUE,
                       lex.order = TRUE)
), summarise_version)
version_summary <- do.call(rbind, summary_rows)
version_summary <- version_summary[order(version_summary$version,
                                         version_summary$n), , drop = FALSE]
utils::write.csv(version_summary,
                 file.path(out_dir, "stage2b1_raw_preparation_summary.csv"),
                 row.names = FALSE)

get_row <- function(z, version, n) {
  out <- z[z$version == version & z$n == n, , drop = FALSE]
  if (nrow(out) != 1L) stop("missing summary row", call. = FALSE)
  out[1L, , drop = FALSE]
}

paired_rows <- lapply(n_values, function(n) {
  before <- get_row(version_summary, "before", n)
  after <- get_row(version_summary, "after", n)
  raw_speedup <- before$raw_median_ms / after$raw_median_ms
  prepare_speedup <- before$prepare_median_ms / after$prepare_median_ms
  prepared_speedup <- before$prepared_median_ms / after$prepared_median_ms
  speedup_gate <- if (n >= 5000L) raw_speedup >= 1.7 else NA
  n20000_gate <- if (n == 20000L) raw_speedup >= 2.0 else NA
  residual_gate <- after$duplicate_residual_pct <= 15
  overhead_gate <- after$raw_overhead_ratio <= 1.20
  required <- c(speedup_gate, n20000_gate, residual_gate, overhead_gate)
  accepted <- if (all(is.na(required))) NA else all(required[!is.na(required)])
  data.frame(
    n = as.integer(n),
    before_raw_median_ms = before$raw_median_ms,
    after_raw_median_ms = after$raw_median_ms,
    raw_speedup = raw_speedup,
    before_prepare_median_ms = before$prepare_median_ms,
    after_prepare_median_ms = after$prepare_median_ms,
    prepare_speedup = prepare_speedup,
    before_prepared_median_ms = before$prepared_median_ms,
    after_prepared_median_ms = after$prepared_median_ms,
    prepared_speedup = prepared_speedup,
    after_duplicate_residual_ms = after$duplicate_residual_ms,
    after_duplicate_residual_pct = after$duplicate_residual_pct,
    after_raw_overhead_ratio = after$raw_overhead_ratio,
    speedup_ge_1_7x = speedup_gate,
    n20000_speedup_ge_2x = n20000_gate,
    duplicate_residual_le_15pct = residual_gate,
    raw_overhead_le_1_20 = overhead_gate,
    stage2b1_gate = accepted,
    stringsAsFactors = FALSE
  )
})
paired_summary <- do.call(rbind, paired_rows)
utils::write.csv(paired_summary,
                 file.path(out_dir, "stage2b1_before_after_summary.csv"),
                 row.names = FALSE)

manifest <- data.frame(
  file = list.files(out_dir, pattern = "^stage2b1_", full.names = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(out_dir, "stage2b1_manifest.csv"),
                 row.names = FALSE)

print(version_summary)
print(paired_summary)
