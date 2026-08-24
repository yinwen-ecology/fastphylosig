# Development-only PSOCK failure harness. It deliberately uses local closures
# instead of monkey-patching production functions, whose public API has no
# fault-injection seam.

.worker_audit_closed <- function(cluster) {
  inherits(try(parallel::clusterCall(cluster, function() TRUE), silent = TRUE),
           "try-error")
}

.worker_audit_run <- function(kind, requested = 4L) {
  skip_if(parallel::detectCores(logical = TRUE) < 2L,
          "PSOCK worker audit requires two logical cores")
  cluster <- parallel::makeCluster(2L)
  on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
  values <- parallel::parLapply(cluster, seq_len(requested), function(i, kind) {
    tryCatch({
      if (identical(kind, "computation_error")) stop("injected computation error")
      if (identical(kind, "one_permutation_failure") && i == 2L) {
        stop("injected permutation error")
      }
      list(ok = TRUE, value = as.numeric(i))
    }, error = function(e) list(ok = FALSE, value = NA_real_))
  }, kind = kind)
  parallel::stopCluster(cluster)
  closed <- .worker_audit_closed(cluster)
  cluster <- NULL
  successful <- sum(vapply(values, function(x) isTRUE(x$ok) && is.finite(x$value), logical(1)))
  list(successful = successful, failed = requested - successful, closed = closed)
}

test_that("injected worker computation and one-permutation failures retain counts and cleanup", {
  for (kind in c("computation_error", "one_permutation_failure")) {
    got <- .worker_audit_run(kind)
    expect_equal(got$successful + got$failed, 4L)
    expect_true(got$closed)
  }
  expect_equal(.worker_audit_run("computation_error")$successful, 0L)
  expect_equal(.worker_audit_run("one_permutation_failure")$successful, 3L)
})

test_that("an injected initialization error closes its PSOCK cluster", {
  skip_if(parallel::detectCores(logical = TRUE) < 2L,
          "PSOCK initialization audit requires two logical cores")
  cluster <- parallel::makeCluster(2L)
  on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
  got <- try(parallel::clusterCall(cluster, function() {
    stop("injected initialization error")
  }), silent = TRUE)
  expect_s3_class(got, "try-error")
  parallel::stopCluster(cluster)
  expect_true(.worker_audit_closed(cluster))
})

test_that("a normal parallel Delta call remains available after the injected audit", {
  skip_if_not_installed("ape")
  skip_if(parallel::detectCores(logical = TRUE) < 2L,
          "PSOCK recovery check requires two logical cores")
  invisible(.worker_audit_run("one_permutation_failure"))
  set.seed(20260813L)
  tree <- ape::rtree(12L)
  x <- stats::setNames(rep(c("A", "B", "C"), length.out = 12L), tree$tip.label)
  permutations <- rbind(seq_len(12L), 12:1L, c(2:12, 1L), c(12L, 1:11))
  fit <- suppressWarnings(fastphylosig::fast_delta(
    tree, x, test = TRUE, nsim = 4L, permutations = permutations,
    mcmc_sim = 20L, burn = 5L, thin = 5L, model = "ER", ncores = 2L,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))
  expect_equal(fit$requested_simulations, 4L)
  expect_equal(fit$successful_simulations + fit$n_failed_sim, 4L)
  expect_equal(fit$successful_simulations, 4L)
  expect_equal(fit$n_failed_sim, 0L)
  expect_identical(fit$status, "ok")
})
