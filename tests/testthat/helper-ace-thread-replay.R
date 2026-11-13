# PORT can dispatch MKL dnrm2 differently with different work-array alignment.
# Replay actual optimizer inputs to test our kernel, not PORT's bitwise path.
.expect_ace_thread_replay <- function(w, state, runs) {
  for (threads in c(1L, 2L, 4L)) {
    ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, state,
      w$rate$rate_index, threads)
    tryCatch({
      for (run in runs) {
        seed <- .Random.seed
        for (call in run$calls) {
          expect_identical(fast_ace_fit_deviance_cpp(ptr, call$par), call$value)
        }
        expect_identical(fast_ace_fit_deviance_cpp(ptr, run$fit$par),
          run$fit$objective)
        expect_identical(.Random.seed, seed)
      }
    }, finally = fast_ace_fit_release_cpp(ptr))
  }
  for (run in runs[-1L]) {
    expect_identical(run$seed, runs[[1L]]$seed)
    expect_identical(run$next_normal, runs[[1L]]$next_normal)
    expect_identical(run$next_uniform, runs[[1L]]$next_uniform)
  }
}