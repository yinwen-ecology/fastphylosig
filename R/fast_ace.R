# Fast discrete ancestral-state reconstruction -------------------------------

fast_ace <- function(x, phy, type = "discrete", method = "ML", CI = TRUE,
                     model = c("ER", "ARD"), kappa = 1, ip = 0.1,
                     marginal = FALSE) {
  if (!inherits(phy, "phylo")) {
    stop("object \"phy\" is not of class \"phylo\".", call. = FALSE)
  }
  if (is.null(phy$edge.length)) {
    stop("tree has no branch lengths.", call. = FALSE)
  }
  if (!is.logical(CI) || length(CI) != 1L || is.na(CI) ||
      !is.logical(marginal) || length(marginal) != 1L || is.na(marginal)) {
    stop("CI and marginal must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      !is.numeric(ip) || length(ip) != 1L || !is.finite(ip) || ip < 0) {
    stop("kappa must be finite and ip must be finite and non-negative.",
         call. = FALSE)
  }
  type <- match.arg(type, c("discrete"))
  if (method != "ML") {
    stop("fast_ace currently supports only method = \"ML\".", call. = FALSE)
  }
  model <- match.arg(model)

  nb_tip <- length(phy$tip.label)
  nb_node <- phy$Nnode
  if (nb_node != nb_tip - 1L) {
    stop("\"phy\" must be rooted and fully dichotomous.", call. = FALSE)
  }
  if (length(x) != nb_tip) {
    stop("length of trait data and tree tips do not match.", call. = FALSE)
  }
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x))) {
      stop("names of x must be unique.", call. = FALSE)
    }
    if (!setequal(names(x), phy$tip.label)) {
      stop("names of x must match phy$tip.label exactly.", call. = FALSE)
    }
    x <- x[phy$tip.label]
  }
  if (kappa != 1) {
    phy$edge.length <- phy$edge.length^kappa
  }
  if (any(phy$edge.length < 0)) {
    stop("some branches have negative length.", call. = FALSE)
  }

  x <- droplevels(factor(x))
  lvls <- levels(x)
  nl <- nlevels(x)
  if (nl < 2L) {
    stop("x must contain at least two observed states.", call. = FALSE)
  }
  tip_state <- as.integer(x)
  rate <- .fast_ace_rate_index(model, nl)

  phy <- ape::reorder.phylo(phy, "postorder")
  edge <- matrix(as.integer(phy$edge), ncol = 2L)
  edge_length <- as.numeric(phy$edge.length)
  par0 <- rep(ip, length.out = rate$np)

  dev <- function(p) {
    fast_ace_discrete_deviance_cpp(
      edge = edge, edge_length = edge_length, tip_state = tip_state,
      rate_index = rate$rate_index, par = p
    )
  }

  fit <- stats::nlminb(
    start = par0, objective = dev,
    lower = rep(0, rate$np), upper = rep(1e50, rate$np)
  )

  hessian_fit <- try(
    suppressWarnings(stats::nlm(
      f = dev, p = fit$par, iterlim = 1, stepmax = 0, hessian = TRUE
    )),
    silent = TRUE
  )
  rate_se <- if (inherits(hessian_fit, "try-error") ||
      any(diag(hessian_fit$hessian) == 0)) {
    rep(NaN, rate$np)
  } else {
    solved <- try(solve(hessian_fit$hessian), silent = TRUE)
    if (inherits(solved, "try-error")) {
      rep(NaN, rate$np)
    } else {
      sqrt(diag(solved))
    }
  }

  obj <- list(
    loglik = -0.5 * fit$objective,
    rates = fit$par,
    se = rate_se,
    index.matrix = rate$index_matrix,
    convergence = fit$convergence,
    message = fit$message
  )
  if (fit$convergence != 0L) {
    warning("transition-rate optimization did not converge: ", fit$message,
            call. = FALSE)
  }

  if (isTRUE(CI)) {
    lik <- fast_ace_discrete_liks_cpp(
      edge = edge, edge_length = edge_length, tip_state = tip_state,
      rate_index = rate$rate_index, par = fit$par, marginal = marginal
    )
    lik_anc <- lik$lik.anc
    rownames(lik_anc) <- nb_tip + seq_len(nb_node)
    colnames(lik_anc) <- lvls
    obj$lik.anc <- lik_anc
  }

  obj$call <- match.call()
  obj$engine <- "fastphylosig::fast_ace"
  class(obj) <- "ace"
  obj
}

.fast_ace_rate_index <- function(model, nl) {
  rate_index <- matrix(0L, nl, nl)
  if (model == "ER") {
    np <- 1L
    rate_index[col(rate_index) != row(rate_index)] <- 1L
  } else if (model == "ARD") {
    np <- nl * (nl - 1L)
    rate_index[col(rate_index) != row(rate_index)] <- seq_len(np)
  } else {
    stop("fast_ace currently supports only model = \"ER\" or \"ARD\".",
         call. = FALSE)
  }

  index_matrix <- rate_index
  diag(index_matrix) <- NA_integer_
  list(rate_index = rate_index, index_matrix = index_matrix, np = np)
}
