.review_delta_shape_pair <- function(tree, x, ...) {
  set.seed(9127)
  scalar <- suppressWarnings(fast_delta(tree, x, model = "ER",
    mcmc_sim = 80, burn = 10, thin = 2, progress = FALSE, ...))
  scalar_seed <- .Random.seed
  set.seed(9127)
  batch <- suppressWarnings(fast_delta(tree,
    matrix(x, dimnames = list(names(x), "x")), model = "ER",
    mcmc_sim = 80, burn = 10, thin = 2, progress = FALSE, ...))
  expect_identical(.Random.seed, scalar_seed)
  fields <- c(names(.delta_calculation_flags()), "status", "message")
  expect_true(all(fields %in% names(scalar)))
  for (nm in fields) expect_identical(scalar[[nm]], batch[[nm]][[1L]], info = nm)
  list(scalar = scalar, batch = batch)
}

test_that("scalar Delta retains healthy calculation audit fields", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  for (test in c(FALSE, TRUE)) {
    got <- .review_delta_shape_pair(tree, x, test = test, nsim = 2,
      permutations = rbind(1:4, 4:1), return_sim = TRUE)
    expect_identical(got$scalar$calculation_failed, FALSE)
    expect_identical(got$scalar$failure_stage, NA_character_)
  }
})

test_that("scalar Delta retains failures raised by ancestral estimation", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  testthat::local_mocked_bindings(
    .fast_ace_lik_anc_workspace = function(...) stop("injected ancestral failure"))
  got <- .review_delta_shape_pair(tree, x)
  expect_identical(got$scalar$calculation_failed, TRUE)
  expect_identical(got$scalar$failure_stage, "ACE")
  expect_identical(got$scalar$status, "failed")
  expect_match(got$scalar$message, "injected ancestral failure")
})

test_that("scalar Delta retains sampler failure and preceding warnings", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  testthat::local_mocked_bindings(
    .fast_ace_lik_anc_workspace = function(...) {
      warning("injected ancillary warning")
      matrix(.5, nrow = 3, ncol = 2)
    },
    delta_mcmc_cpp = function(...) stop("injected sampler failure"))
  got <- .review_delta_shape_pair(tree, x)
  expect_identical(got$scalar$calculation_failed, TRUE)
  expect_identical(got$scalar$failure_stage, "MCMC")
  expect_identical(got$scalar$other_warning, TRUE)
  expect_identical(got$scalar$status, "failed")
  expect_match(got$scalar$message, "injected sampler failure")
})

test_that("scalar Delta retains the fallback calculation failure stage", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  testthat::local_mocked_bindings(.delta_one = function(...) stop("injected calculation failure"))
  got <- .review_delta_shape_pair(tree, x)
  expect_identical(got$scalar$calculation_failed, TRUE)
  expect_identical(got$scalar$failure_stage, "calculation")
  expect_identical(got$scalar$status, "failed")
})
test_that("unstructured Delta failures do not delete later trait audit slots", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- cbind(failed = c(a="A", b="A", c="B", d="B"),
             healthy = c(a="B", b="A", c="B", d="A"))
  original <- .delta_one
  testthat::local_mocked_bindings(.delta_one = function(tree, y, ...) {
    if (as.character(y[[1L]]) == "A") stop("injected first-trait failure")
    original(tree, y, ...)
  })
  set.seed(9128)
  got <- suppressWarnings(fast_delta(tree, x, model="ER", mcmc_sim=80,
    burn=10, thin=2, progress=FALSE))
  expect_identical(nrow(got), 2L)
  expect_identical(got$calculation_failed, c(TRUE, FALSE))
  expect_identical(got$failure_stage, c("calculation", NA_character_))
  expect_identical(got$status, c("failed", "ok"))
  expect_true(is.finite(got$Delta_fast[[2L]]))
})