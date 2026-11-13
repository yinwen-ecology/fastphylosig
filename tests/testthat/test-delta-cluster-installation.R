test_that("Delta PSOCK pins lib.loc installation without changing caller libraries", {
  skip_if_not_installed("ape")
  expected <- normalizePath(getNamespaceInfo(asNamespace("fastphylosig"), "path"),
                            winslash = "/", mustWork = TRUE)
  prior <- .libPaths()
  # Reproduce library(..., lib.loc=...) usage: its installation is loaded but
  # need not be in .libPaths(). Workers must pin that actual namespace.
  on.exit(.libPaths(prior), add = TRUE)
  .libPaths(prior[normalizePath(prior, winslash = "/") != dirname(expected)])
  caller_paths <- .libPaths()
  checked <- FALSE
  original <- fastphylosig:::.delta_initialize_cluster
  testthat::local_mocked_bindings(
    .delta_initialize_cluster = function(cl) {
      info <- original(cl)
      for (worker in info) {
        expect_identical(worker$path, expected)
        expect_identical(worker$version,
          as.character(getNamespaceVersion("fastphylosig")))
        expect_identical(worker$library_paths[[1]], dirname(expected))
        expect_true(all(caller_paths %in% worker$library_paths))
      }
      checked <<- TRUE
      info
    }, .package = "fastphylosig"
  )
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  set.seed(20261112)
  fit <- suppressWarnings(fastphylosig::fast_delta(
    tree, x, test = TRUE, nsim = 2, ncores = 2, model = "ER",
    mcmc_sim = 100, burn = 10, thin = 10, return_sim = TRUE,
    progress = FALSE, verbose = FALSE))
  expect_true(checked)
  expect_equal(fit$successful_simulations, 2)
  expect_equal(fit$n_failed_sim, 0)
  expect_identical(.libPaths(), caller_paths)
})

test_that("Delta PSOCK rejects an already loaded different installation", {
  expected <- normalizePath(getNamespaceInfo(asNamespace("fastphylosig"), "path"),
                            winslash = "/", mustWork = TRUE)
  alternate <- tempfile("delta-alternate-library-")
  dir.create(alternate)
  on.exit(unlink(alternate, recursive = TRUE), add = TRUE)
  expect_true(file.copy(expected, alternate, recursive = TRUE))
  cl <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cl), add = TRUE, after = FALSE)
  preload <- function(lib) {
    ns <- loadNamespace("fastphylosig", lib.loc = lib)
    normalizePath(getNamespaceInfo(ns, "path"), winslash = "/")
  }
  environment(preload) <- baseenv()
  loaded <- parallel::clusterCall(cl, preload, alternate)[[1]]
  expect_false(identical(loaded, expected))
  expect_error(fastphylosig:::.delta_initialize_cluster(cl),
               "different fastphylosig installation", fixed = TRUE)
})
test_that("Delta PSOCK rejects a version mismatch before computation", {
  cl <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  original_call <- parallel::clusterCall
  testthat::local_mocked_bindings(
    clusterCall = function(cl, fun, ...) {
      args <- list(...)
      args[[2L]] <- "0.0.0-unexpected"
      do.call(original_call, c(list(cl = cl, fun = fun), args))
    }, .package = "parallel"
  )
  expect_error(fastphylosig:::.delta_initialize_cluster(cl),
               "different fastphylosig installation", fixed = TRUE)
})
test_that("Delta PSOCK removes namespace-bearing source references", {
  original <- fastphylosig:::.delta_initialize_cluster
  code <- deparse(original, control = "all")
  source_file <- srcfilecopy("delta-keep-source-regression.R", code)
  parent.env(source_file) <- asNamespace("fastphylosig")
  with_source <- eval(parse(text = code, srcfile = source_file,
                            keep.source = TRUE),
                      envir = asNamespace("fastphylosig"))
  expect_false(is.null(attr(with_source, "srcref")))
  cl <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  info <- with_source(cl)[[1L]]
  expected <- normalizePath(getNamespaceInfo(asNamespace("fastphylosig"), "path"),
                            winslash = "/", mustWork = TRUE)
  expect_identical(info$path, expected)
  expect_identical(info$version,
                   as.character(getNamespaceVersion("fastphylosig")))
})
