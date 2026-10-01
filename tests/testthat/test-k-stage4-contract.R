# The full 0.2.0-vs-candidate comparison runs in isolated Rscript processes so
# both package versions can be loaded without namespace collision.  Set these
# library paths to include it in a package test run; without them the contract
# remains available as `Rscript test/test_rng_contract.R <v020-lib> <v030-lib>`.
test_that("the installed 0.2.0 RNG contract holds when both libraries are supplied", {
  baseline_lib <- Sys.getenv("FASTPHYLOSIG_V020_LIB", unset = "")
  candidate_lib <- Sys.getenv("FASTPHYLOSIG_CANDIDATE_LIB", unset = "")
  if (!nzchar(baseline_lib) || !nzchar(candidate_lib)) {
    skip("set FASTPHYLOSIG_V020_LIB and FASTPHYLOSIG_CANDIDATE_LIB for the cross-install contract")
  }

  script <- testthat::test_path("..", "..", "test", "test_rng_contract.R")
  if (!file.exists(script)) {
    skip("test/test_rng_contract.R is not included in this package source tree")
  }

  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (.Platform$OS.type != "windows") {
    rscript <- file.path(R.home("bin"), "Rscript")
  }
  output <- suppressWarnings(system2(
    rscript,
    args = c("--vanilla", shQuote(normalizePath(script, winslash = "/")),
             shQuote(normalizePath(baseline_lib, winslash = "/")),
             shQuote(normalizePath(candidate_lib, winslash = "/"))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L,
                   info = paste(tail(output, 20L), collapse = "\n"))
  expect_true(any(grepl("STAGE4 RNG CONTRACT PASS", output, fixed = TRUE)))
})
