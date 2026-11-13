# Compare different floating-point implementations, not repeated calls to the
# same implementation. Shapes, field names, missingness, infinities and the
# invalid-objective sentinel remain exact; only finite arithmetic has a bound.
.expect_ace_numerical_equivalent <- function(actual, expected, tolerance = 1e-12) {
  expect_identical(typeof(actual), typeof(expected))
  expect_identical(attributes(actual), attributes(expected))
  expect_identical(length(actual), length(expected))
  if (!identical(typeof(actual), typeof(expected)) ||
      !identical(attributes(actual), attributes(expected)) ||
      length(actual) != length(expected)) return(invisible(NULL))
  if (is.list(actual)) {
    for (i in seq_along(actual))
      .expect_ace_numerical_equivalent(actual[[i]], expected[[i]], tolerance)
  } else if (is.double(actual) && !is.complex(actual)) {
    finite <- is.finite(expected)
    expect_identical(is.finite(actual), finite)
    expect_identical(is.na(actual), is.na(expected))
    expect_identical(is.nan(actual), is.nan(expected))
    expect_identical(actual[!finite], expected[!finite])
    if (length(actual) == 1L && identical(as.numeric(expected), 1e50)) {
      expect_identical(actual, expected)
    } else if (identical(is.finite(actual), finite) && any(finite)) {
      scaled <- abs(actual[finite] - expected[finite]) /
        pmax(1, abs(actual[finite]), abs(expected[finite]))
      expect_true(all(scaled <= tolerance),
        info = paste("maximum componentwise scaled error:", max(scaled),
                     "; bound:", tolerance))
    }
  } else {
    expect_identical(actual, expected)
  }
  invisible(NULL)
}
