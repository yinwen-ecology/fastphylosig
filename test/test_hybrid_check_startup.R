# Development regression: R CMD check startup state must not leak into
# sourceCpp's build subprocesses. Run from a checkout with the package installed:
# Rscript --vanilla test/test_hybrid_check_startup.R <repo> <library> [baseline.R]
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".",
                      winslash = "/", mustWork = TRUE)
if (length(args) >= 2L) .libPaths(c(args[[2L]], .libPaths()))
script <- file.path(root, "tests", "test_hybrid_equivalence.R")
baseline <- normalizePath(if (length(args) >= 3L) args[[3L]] else script,
                          winslash = "/", mustWork = TRUE)
out_dir <- Sys.getenv("FASTPHYLOSIG_STARTUP_AUDIT", tempfile("startup-audit-"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
run_dir <- tempfile("startup-parent-", tmpdir = file.path(root, "docs"))
dir.create(run_dir)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
stopifnot(startsWith(run_dir, paste0(root, "/docs/")))
old_wd <- setwd(run_dir)
old_tests <- Sys.getenv("R_TESTS", unset = NA_character_)
writeLines("stop('The parent check startup hook leaked into a child')", "startup.Rs")
startup <- normalizePath("startup.Rs", winslash = "/", mustWork = TRUE)

# Record each existing native result, without adding another kernel call or draw.
run_equivalence <- function(path) {
  code <- readLines(path, warn = FALSE)
  call <- "out <- test_env$stage4_test_k_four_modes(compiled, X, permutations)"
  stopifnot(sum(grepl(call, code, fixed = TRUE)) == 1L)
  code <- sub(call, paste0(call, "\n    all_modes[[tree_kind]] <- out"),
              code, fixed = TRUE)
  scope <- new.env(parent = globalenv())
  scope$all_modes <- list()
  # Keep tree discovery in the checkout even though the startup hook lives in
  # another directory, as it does for the check/fastphylosig.Rcheck/tests layout.
  setwd(run_dir)
  set.seed(81371)
  eval(parse(text = code), envir = scope)
  stopifnot(nrow(scope$results) == 27L, length(scope$all_modes) == 3L)
  list(comparisons = scope$results, modes = scope$all_modes,
       seed = .Random.seed, compile = scope$compile_stage4_harness)
}

tryCatch({
  Sys.unsetenv("R_TESTS")
  reference <- run_equivalence(baseline)
  rows <- list()
  for (mode in c("unset", "empty", "relative", "absolute")) {
    if (mode == "unset") Sys.unsetenv("R_TESTS") else
      Sys.setenv(R_TESTS = switch(mode, empty = "", relative = "startup.Rs",
                                 absolute = startup))
    inherited <- Sys.getenv("R_TESTS", unset = NA_character_)
    result <- run_equivalence(script)
    stopifnot(identical(result$comparisons, reference$comparisons),
              identical(result$modes, reference$modes),
              identical(result$seed, reference$seed),
              identical(Sys.getenv("R_TESTS", unset = NA_character_), inherited))

    # A genuine compiler failure must also restore the parent's environment.
    bad <- tempfile(fileext = ".cpp")
    writeLines("#error intentional_startup_regression_failure", bad)
    seed_before <- .Random.seed
    failure <- tryCatch({
      result$compile(bad, new.env(parent = globalenv()))
      NULL
    }, error = conditionMessage)
    unlink(bad)
    stopifnot(!is.null(failure),
              identical(Sys.getenv("R_TESTS", unset = NA_character_), inherited),
              identical(.Random.seed, seed_before))
    rows[[mode]] <- data.frame(
      startup = mode, topology_trait_mode_comparisons = nrow(result$comparisons),
      all_native_outputs_identical = identical(result$modes, reference$modes),
      final_rng_identical = identical(result$seed, reference$seed),
      environment_restored_on_success = TRUE,
      environment_and_rng_restored_on_error = TRUE
    )
  }
  results <- do.call(rbind, rows)
  utils::write.csv(results, file.path(out_dir, "startup-regression.csv"),
                   row.names = FALSE)
  saveRDS(list(reference = reference[c("comparisons", "modes", "seed")],
               checks = results), file.path(out_dir, "startup-regression.rds"))
  print(results, row.names = FALSE)
  cat("STARTUP SUBPROCESS / RESULT / RNG CONTRACT PASS\n")
}, finally = {
  setwd(old_wd)
  if (is.na(old_tests)) Sys.unsetenv("R_TESTS") else Sys.setenv(R_TESTS = old_tests)
  unlink(run_dir, recursive = TRUE)
})
