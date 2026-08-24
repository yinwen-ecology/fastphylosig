test_that("warning/status stability audit has a valid machine-readable schema", {
  path <- testthat::test_path("fixtures", "warning_status_stability.csv")
  skip_if_not(file.exists(path), "warning/status audit artifact is unavailable")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  expect_true(nrow(x) >= 1L)
  expect_true(all(c("method", "audit_status", "result_class", "fields",
                    "status_field", "requested_field", "successful_field",
                    "failed_field", "warning_count", "warning_classes",
                    "package_version") %in% names(x)))
  expect_true(all(x$audit_status %in% c("PASS", "NOT_RUN")))
  expect_true(all(is.na(x$warning_count) | x$warning_count >= 0))
})
