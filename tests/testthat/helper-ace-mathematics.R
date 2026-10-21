# Mathematical references: analytic binary CTMC / independent matrix exponent,
# exhaustive latent-state sums, and upward likelihood curvature. No ACE downpass.
.ace_math_transition <- function(Q, t) {
  if (nrow(Q) == 2L) {
    a <- Q[1L, 2L]; b <- Q[2L, 1L]; s <- a + b
    if (s == 0) return(diag(2L))
    z <- exp(-s * t)
    return(matrix(c(b + a*z, a*(1-z), b*(1-z), a + b*z),
                  2L, byrow = TRUE) / s)
  }
  as.matrix(expm::expm(Q * t))
}

.ace_math_q <- function(index, rates) {
  Q <- matrix(0, nrow(index), ncol(index))
  off <- !is.na(index) & index > 0L
  Q[off] <- rates[index[off]]
  diag(Q) <- -rowSums(Q)
  Q
}

.ace_math_enumerate <- function(tree, states, Q) {
  n <- length(tree$tip.label); k <- nrow(Q)
  nodes <- n + seq_len(tree$Nnode)
  P <- lapply(tree$edge.length, function(t) .ace_math_transition(Q, t))
  grid <- as.matrix(expand.grid(rep(list(seq_len(k)), tree$Nnode)))
  joint <- matrix(0, length(nodes), k, dimnames = list(nodes, seq_len(k)))
  likelihood <- 0
  for (r in seq_len(nrow(grid))) {
    assignment <- c(states, grid[r, ])
    mass <- 1 / k
    for (e in seq_len(nrow(tree$edge))) {
      child <- tree$edge[e, 2L]
      # Unknown observed tips integrate to a transition-row sum of one.
      if (!is.na(assignment[child]))
        mass <- mass * P[[e]][assignment[tree$edge[e, 1L]], assignment[child]]
    }
    likelihood <- likelihood + mass
    for (j in seq_along(nodes))
      joint[j, assignment[nodes[j]]] <- joint[j, assignment[nodes[j]]] + mass
  }
  list(probability = joint / likelihood, loglik = log(likelihood))
}

.ace_math_loglik <- function(tree, states, Q) {
  phy <- ape::reorder.phylo(tree, "postorder")
  n <- length(phy$tip.label); k <- nrow(Q)
  lik <- matrix(1, n + phy$Nnode, k)
  for (i in seq_len(n)) if (!is.na(states[i]))
    lik[i, ] <- as.numeric(seq_len(k) == states[i])
  log_scale <- 0
  for (parent in unique(phy$edge[, 1L])) {
    value <- rep(1, k)
    for (e in which(phy$edge[, 1L] == parent))
      value <- value * drop(.ace_math_transition(Q, phy$edge.length[e]) %*%
                              lik[phy$edge[e, 2L], ])
    s <- sum(value)
    lik[parent, ] <- value / s
    log_scale <- log_scale + log(s)
  }
  log_scale - log(k)
}
