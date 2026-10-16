# Independent likelihood oracle for the ARD recovery regression fixtures.
# Evaluate a GIVEN rate vector: an optimizer may stop at different valid points
# across R / BLAS / LAPACK configurations. Never compare a newly optimized
# value with an oracle evaluated at a different, historical parameter vector.
.ace_recovery_reference_deviance <- function(a, rates) {
  w <- a$workspace
  k <- w$nl
  q <- matrix(0, k, k)
  index <- w$rate$rate_index
  q[index > 0L] <- rates[index[index > 0L]]
  diag(q) <- -rowSums(q)

  likelihood <- matrix(1, w$nb_tip + w$nb_node, k)
  for (tip in seq_len(w$nb_tip)) {
    if (!is.na(a$tip_state[tip])) {
      likelihood[tip, ] <- as.numeric(seq_len(k) == a$tip_state[tip])
    }
  }
  log_scale <- numeric(nrow(likelihood))
  # The captured workspace edges are in postorder, children before parents.
  for (parent in unique(w$edge[, 1L])) {
    rows <- which(w$edge[, 1L] == parent)
    value <- rep(1, k)
    for (row in rows) {
      child <- w$edge[row, 2L]
      transition <- expm::expm(q * w$edge_length[row])
      value <- value * as.vector(transition %*% likelihood[child, ])
      log_scale[parent] <- log_scale[parent] + log_scale[child]
    }
    normalizer <- sum(value)
    stopifnot(is.finite(normalizer), normalizer > 0)
    likelihood[parent, ] <- value / normalizer
    log_scale[parent] <- log_scale[parent] + log(normalizer)
  }
  root <- setdiff(w$edge[, 1L], w$edge[, 2L])
  stopifnot(length(root) == 1L)
  # ACE uses the root-sum convention, not a normalized uniform-root prior.
  -2 * (log_scale[root] + log(sum(likelihood[root, ])))
}
