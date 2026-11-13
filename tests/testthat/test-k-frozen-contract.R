test_that("public K complete objects and RNG retain the frozen contract", {
  frozen <- readRDS(test_path("fixtures", "k-public-contract.rds"))
  actual <- .compat_k_capture()
  expect_identical(names(actual), names(frozen$cases))
  for (name in names(actual)) {
    a <- actual[[name]]
    e <- frozen$cases[[name]]
    expect_identical(a$before, e$before, info=name)
    expect_identical(a$after, e$after, info=name)
    expect_identical(names(a$result), names(e$result), info=name)
    expect_identical(class(a$result), class(e$result), info=name)
    expect_equal(a$result, e$result, tolerance=1e-12, info=name)
    if (grepl("TRUE$", name)) expect_identical(a$before, a$after)
  }
})
