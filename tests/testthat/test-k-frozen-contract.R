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

test_that("frozen K capture restores the caller's RNG kind and seed presence", {
  run <- function() {
    original_kind <- RNGkind()
    original_present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (original_present) original_seed <- get(".Random.seed", .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(original_kind))
      if (original_present) assign(".Random.seed", original_seed, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    })
    set.seed(924L)
    before <- .Random.seed
    invisible(.compat_k_capture())
    expect_identical(RNGkind(), original_kind)
    expect_identical(.Random.seed, before)
    rm(".Random.seed", envir = .GlobalEnv)
    invisible(.compat_k_capture())
    expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    expect_identical(RNGkind(), original_kind)
  }
  run()
})
