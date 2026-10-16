# Independent R oracle: evaluates the full joint log posterior, without the
# production conditional simplification or accepted-state gamma cache.
.delta_posterior_reference <- function(x, rate, sd, sim, thin, burn) {
  target <- function(a, b) {
    length(x) * (lgamma(a + b) - lgamma(a) - lgamma(b)) +
      (a - 1) * sum(log(x)) + (b - 1) * sum(log1p(-x)) - rate * (a + b)
  }
  update <- function(a, b, which) {
    old <- if (which == 1L) a else b
    log_new <- rnorm(1, log(old), sd)
    new <- exp(log_new)
    if (!is.finite(new) || new <= 0) return(old)
    proposed <- if (which == 1L) target(new, b) else target(a, new)
    ratio <- proposed - target(a, b) + log_new - log(old)
    if (is.finite(ratio) && (ratio >= 0 || log(runif(1)) < ratio)) new else old
  }
  alpha <- beta <- matrix(NA_real_, 1L + (sim - burn) %/% thin, 2)
  for (chain in 1:2) {
    a <- rexp(1); b <- rexp(1); saved <- 0L
    for (iter in seq_len(sim)) {
      a <- update(a, b, 1L); b <- update(a, b, 2L)
      if (iter >= burn && (iter - burn) %% thin == 0) {
        saved <- saved + 1L; alpha[saved, chain] <- a; beta[saved, chain] <- b
      }
    }
  }
  list(alpha = alpha, beta = beta)
}
