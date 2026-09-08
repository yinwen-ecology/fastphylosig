# Canonicalization Contract V2-B full testthat evidence runner.

output <- Sys.getenv("FASTPHYLOSIG_TEST_LOG", unset = "v2b_full_test.log")
connection <- file(output, open = "wt", encoding = "UTF-8")
sink(connection)
sink(connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(connection)
}, add = TRUE)

testthat::set_max_fails(Inf)
results <- devtools::test(".", stop_on_failure = FALSE)
summary <- as.data.frame(results)
totals <- c(
  FAIL = sum(summary$failed),
  WARN = sum(summary$warning),
  SKIP = sum(summary$skipped),
  PASS = sum(summary$passed)
)

cat("\nV2B_TEST_TOTALS\n")
cat(paste(names(totals), totals, sep = "="), sep = "\n")
cat("\n")

if (any(totals[c("FAIL", "WARN", "SKIP")] != 0L)) {
  quit(status = 1L)
}
