# Saved-trace postprocessing only; never samples or changes the estimator.
.delta_rank_fields <- function() {
  c(as.vector(outer(c("rank_Rhat", "ESS_bulk", "ESS_tail"),
                    c("alpha", "beta", "Delta"), paste, sep = "_")),
    "relative_MCSE_Delta")
}

.delta_split_chains <- function(x) {
  half <- nrow(x) %/% 2L
  do.call(cbind, lapply(seq_len(ncol(x)), function(j)
    cbind(x[seq_len(half), j], x[seq.int(nrow(x)-half+1L, nrow(x)), j])))
}

.delta_rank_normalize <- function(x) {
  matrix(stats::qnorm((rank(x, ties.method = "average")-.375)/(length(x)+.25)),
         nrow(x), ncol(x))
}

.delta_basic_rhat <- function(x) {
  if (nrow(x) < 2L || ncol(x) < 2L || any(!is.finite(x)) ||
      diff(range(x)) <= .Machine$double.eps) return(NA_real_)
  W <- mean(apply(x, 2L, stats::var))
  if (W <= 0) return(NA_real_)
  sqrt((nrow(x)*stats::var(colMeans(x))/W + nrow(x)-1)/nrow(x))
}

.delta_multichain_ess <- function(x) {
  n <- nrow(x)
  if (n < 3L || any(!is.finite(x)) ||
      diff(range(x)) <= .Machine$double.eps) return(NA_real_)
  # R fallback uses direct acf rather than the native FFT implementation.
  acov <- rowMeans(vapply(seq_len(ncol(x)), function(j) {
    z <- x[,j]-mean(x[,j])
    if (sum(z*z) == 0) return(rep(0, n))
    as.numeric(stats::acf(z, type = "covariance", lag.max = n-1L,
                          plot = FALSE)$acf)
  }, numeric(n)))
  W <- acov[1L]*n/(n-1L)
  V <- acov[1L] + if (ncol(x)>1L) stats::var(colMeans(x)) else 0
  if (!is.finite(V) || V <= 0) return(NA_real_)
  rho <- numeric(n); rho[1L] <- even <- 1
  rho[2L] <- odd <- 1-(W-acov[2L])/V
  lag <- 0L
  while (lag+5L < n && is.finite(even+odd) && even+odd>0) {
    lag <- lag+2L
    even <- 1-(W-acov[lag+1L])/V
    odd <- 1-(W-acov[lag+2L])/V
    if (even+odd>=0) rho[lag+c(1L,2L)] <- c(even,odd)
  }
  if (even>0) rho[lag+1L] <- even
  if (lag>=4L) for (t in seq.int(2L,lag-2L,by=2L)) {
    prev <- sum(rho[c(t-1L,t)])
    if (sum(rho[t+c(1L,2L)])>prev) rho[t+c(1L,2L)] <- prev/2
  }
  sum_rho <- if (lag) sum(rho[seq_len(lag)]) else rho[1L]
  length(x)/max(1/log10(length(x)), -1+2*sum_rho+rho[lag+1L])
}

.delta_rank_diagnostics <- function(x) {
  out <- c(rank_Rhat = NA_real_, ESS_bulk = NA_real_, ESS_tail = NA_real_)
  if (nrow(x)<4L || any(!is.finite(x)) ||
      diff(range(x))<=.Machine$double.eps) return(out)
  z <- .delta_rank_normalize(.delta_split_chains(x))
  folded <- .delta_rank_normalize(.delta_split_chains(abs(x-stats::median(x))))
  out["rank_Rhat"] <- max(.delta_basic_rhat(z), .delta_basic_rhat(folded))
  out["ESS_bulk"] <- .delta_multichain_ess(z)
  out["ESS_tail"] <- min(vapply(c(.05,.95), function(p)
    .delta_multichain_ess(.delta_split_chains(
      1*(x<=stats::quantile(x,p,names=FALSE)))), numeric(1)))
  out
}

.delta_modern_diagnostics <- function(alpha, beta, mcse) {
  values <- list(alpha = alpha, beta = beta, Delta = beta/alpha)
  out <- list()
  for (nm in names(values)) {
    d <- .delta_rank_diagnostics(values[[nm]])
    for (field in names(d)) out[[paste(field,nm,sep="_")]] <- unname(d[field])
  }
  out$relative_MCSE_Delta <- mcse/mean(values$Delta)
  out
}
