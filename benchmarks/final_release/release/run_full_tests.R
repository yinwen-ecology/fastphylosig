# Final release full-test evidence runner.

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else "testthat.Rout"
status_file <- if (length(args) >= 2L) args[[2L]] else "testthat_status.csv"
commit <- if (length(args) >= 3L) args[[3L]] else NA_character_
package_root <- if (length(args) >= 4L) args[[4L]] else "."
package_root <- normalizePath(package_root, winslash = "/", mustWork = TRUE)
output <- normalizePath(output, winslash = "/", mustWork = FALSE)
status_file <- normalizePath(status_file, winslash = "/", mustWork = FALSE)
setwd(package_root)

connection <- file(output, open = "wt", encoding = "UTF-8")
sink(connection)
sink(connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(connection)
}, add = TRUE)

started <- Sys.time()
testthat::set_max_fails(Inf)
results <- devtools::test(".", stop_on_failure = FALSE)
summary <- as.data.frame(results)
totals <- c(
  FAIL = sum(summary$failed),
  WARN = sum(summary$warning),
  SKIP = sum(summary$skipped),
  PASS = sum(summary$passed)
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

cat("\nFINAL_TEST_TOTALS\n")
cat(paste(names(totals), totals, sep = "="), sep = "\n")
cat("\nELAPSED_SECONDS=", format(elapsed, digits = 12), "\n", sep = "")
cat("COMMIT=", commit, "\n", sep = "")

status <- data.frame(
  commit = commit,
  package_version = as.character(desc::desc_get_version(".")),
  pass = unname(totals[["PASS"]]),
  fail = unname(totals[["FAIL"]]),
  warning = unname(totals[["WARN"]]),
  skip = unname(totals[["SKIP"]]),
  elapsed_seconds = elapsed,
  status = if (all(totals[c("FAIL", "WARN", "SKIP")] == 0L)) "PASS" else "FAIL",
  stringsAsFactors = FALSE
)
utils::write.csv(status, status_file, row.names = FALSE, na = "")

if (status$status != "PASS") {
  quit(status = 1L)
}
