# Continuous-trait signal API -------------------------------------------------

fast_signal <- function(tree, x = NULL, method = c("K", "lambda"),
                        test = FALSE, nsim = 1000, se = NULL,
                        start = NULL, control = list(),
                        permutations = NULL, return_sim = test,
                        verbose = TRUE, ncores = 1, X = NULL,
                        lambda_profile = NULL,
                        lambda_profile_points = 101) {
  method <- match.arg(method)
  ncores <- .normalize_ncores(ncores)
  if (!is.logical(test) || length(test) != 1L || is.na(test)) {
    stop("test must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_sim) || length(return_sim) != 1L ||
      is.na(return_sim)) {
    stop("return_sim must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(x)) {
    x <- X
  }
  if (is.null(x)) {
    stop("x must be a named vector, matrix, or data.frame.", call. = FALSE)
  }
  if (!is.null(se)) {
    stop("fast_signal currently supports se = NULL only.", call. = FALSE)
  }

  vector_input <- is.null(dim(x))
  if (is.null(lambda_profile)) {
    lambda_profile <- method == "lambda" && vector_input
  } else {
    lambda_profile <- isTRUE(lambda_profile)
  }
  lambda_profile_points <- as.integer(lambda_profile_points[[1L]])
  if (is.na(lambda_profile_points) || lambda_profile_points < 5L) {
    stop("lambda_profile_points must be an integer >= 5.", call. = FALSE)
  }

  out <- .fast_signal_batch(
    tree = tree, X = x, method = method, test = test, nsim = nsim,
    permutations = permutations, return_sim = return_sim, verbose = verbose,
    ncores = ncores, lambda_profile = lambda_profile,
    lambda_profile_points = lambda_profile_points
  )

  if (!vector_input) {
    return(out)
  }

  if (method == "K") {
    object <- if (!test) {
      out$K_fast[[1]]
    } else {
      z <- list(K = out$K_fast[[1]], P = out$P_fast[[1]])
      if (return_sim && "sim.K_fast" %in% names(out)) {
        z$sim.K <- out$sim.K_fast[[1]]
      }
      z
    }
  } else {
    object <- list(lambda = out$lambda_fast[[1]], logL = out$logL_fast[[1]])
    if (test) {
      object$logL0 <- out$logL0_fast[[1]]
      object$P <- out$P_fast[[1]]
    }
    if (lambda_profile && "lambda_profile_fast" %in% names(out)) {
      object$lambda_profile <- out$lambda_profile_fast[[1]]
      object$lambda_CI <- c(
        lower = out$lambda_CI_lower_fast[[1]],
        upper = out$lambda_CI_upper_fast[[1]]
      )
      object$lambda_CI_level <- 0.95
      object$lambda_CI_cutoff <- out$lambda_CI_cutoff_fast[[1]]
    }
  }

  attr(object, "class") <- "phylosig"
  attr(object, "method") <- method
  attr(object, "test") <- test
  attr(object, "se") <- FALSE
  object
}


# One internal batch engine ----------------------------------------------------

.fast_signal_batch <- function(tree, X, method, test, nsim, permutations,
                               return_sim, verbose, ncores,
                               lambda_profile, lambda_profile_points) {
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  if (method == "K" && test) {
    if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
        nsim < 1 || nsim != floor(nsim)) {
      stop("nsim must be a positive integer when test = TRUE.",
           call. = FALSE)
    }
    nsim <- as.integer(nsim)
  }

  X_raw <- .as_trait_matrix(X, tree)
  tree_tips <- tree$tip.label
  matched <- tree_tips[tree_tips %in% rownames(X_raw)]
  removed_tree <- setdiff(tree_tips, rownames(X_raw))
  removed_data <- setdiff(rownames(X_raw), tree_tips)
  if (length(matched) < 2) {
    stop("Fewer than two species are shared by tree and x.", call. = FALSE)
  }

  base_tree <- if (length(removed_tree)) {
    ape::drop.tip(tree, removed_tree)
  } else {
    tree
  }
  X0 <- X_raw[base_tree$tip.label, , drop = FALSE]
  n0 <- nrow(X0)
  traits <- colnames(X0)
  p <- ncol(X0)

  present <- !is.na(X0)
  n_species <- as.integer(colSums(present))
  n_removed_na <- n0 - n_species
  complete <- n_species == n0
  keys <- character(p)
  keys[complete] <- ".all"
  if (any(!complete)) {
    keys[!complete] <- apply(
      present[, !complete, drop = FALSE], 2,
      function(z) paste(which(z), collapse = ",")
    )
  }
  groups <- split(seq_len(p), keys)

  if (isTRUE(verbose)) {
    message(sprintf(
      paste(
        "Matched %d species; removed %d tree tips and %d data rows before",
        "trait-wise NA pruning."
      ),
      length(matched), length(removed_tree), length(removed_data)
    ))
  }

  K_fast <- P_fast <- lambda_fast <- logL_fast <- logL0_fast <- rep(NA_real_, p)
  lambda_CI_lower_fast <- lambda_CI_upper_fast <- lambda_CI_cutoff_fast <-
    rep(NA_real_, p)
  sim_list <- vector("list", p)
  lambda_profile_list <- vector("list", p)
  notes <- rep(NA_character_, p)
  for (key in names(groups)) {
    idx <- groups[[key]]
    keep <- if (key == ".all") seq_len(n0) else which(present[, idx[[1L]]])
    if (length(keep) < 2) {
      notes[idx] <- "fewer than 2 non-NA matched species"
      next
    }

    group_tree <- if (length(keep) == n0) {
      base_tree
    } else {
      ape::drop.tip(base_tree, base_tree$tip.label[-keep])
    }
    C <- ape::vcv.phylo(group_tree)
    Xg <- X0[rownames(C), idx, drop = FALSE]

    if (method == "K") {
      cholC <- chol(C)
      traceC <- sum(diag(C))

      if (!test) {
        K_fast[idx] <- fast_k_chol_batch_cpp(Xg, cholC, traceC)
      } else {
        n <- nrow(C)
        perms <- .permutation_matrix(
          permutations, n = n, nsim = nsim, include_observed = TRUE
        )
        kres <- if (return_sim) {
          fast_k_chol_permutation_cpp(
            Xg, cholC, traceC, perms, n_threads = ncores
          )
        } else {
          fast_k_chol_permutation_p_cpp(
            Xg, cholC, traceC, perms, n_threads = ncores
          )
        }
        K_fast[idx] <- as.numeric(kres$K)
        P_fast[idx] <- as.numeric(kres$P)
        if (return_sim) {
          for (j in seq_along(idx)) {
            sim_list[[idx[[j]]]] <- kres$sim_K[, j]
          }
        }
      }
    } else {
      maxlam <- .max_lambda(group_tree)
      for (j in seq_along(idx)) {
        col <- idx[[j]]
        y <- Xg[, j]
        lik <- function(lambda) lambda_loglik_cpp(lambda, C, y)
        opt <- stats::optimize(lik, c(0, maxlam), maximum = TRUE)
        lambda_fast[col] <- opt$maximum
        logL_fast[col] <- opt$objective
        if (lambda_profile) {
          profile <- .lambda_profile_data(
            lik = lik, max_lambda = maxlam,
            lambda_hat = opt$maximum, logL_hat = opt$objective,
            n_points = lambda_profile_points
          )
          ci <- .lambda_profile_ci(
            profile = profile, lambda_hat = opt$maximum,
            logL_hat = opt$objective
          )
          lambda_profile_list[[col]] <- profile
          lambda_CI_lower_fast[col] <- ci[["lower"]]
          lambda_CI_upper_fast[col] <- ci[["upper"]]
          lambda_CI_cutoff_fast[col] <- ci[["cutoff"]]
        }
        if (test) {
          logL0_fast[col] <- lik(0)
          P_fast[col] <- stats::pchisq(
            2 * (opt$objective - logL0_fast[col]),
            df = 1, lower.tail = FALSE
          )
        }
      }
    }
  }

  common <- data.frame(
    trait = traits,
    n_species = n_species,
    n_removed_na = n_removed_na,
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    note = notes,
    stringsAsFactors = FALSE
  )
  out <- if (method == "K") {
    data.frame(trait = traits, K_fast = K_fast, common[-1],
               stringsAsFactors = FALSE)
  } else {
    data.frame(trait = traits, lambda_fast = lambda_fast,
               logL_fast = logL_fast, common[-1], stringsAsFactors = FALSE)
  }

  if (method == "K" && test) {
    out$P_fast <- P_fast
    if (return_sim) out$sim.K_fast <- I(sim_list)
  }
  if (method == "lambda" && test) {
    out$logL0_fast <- logL0_fast
    out$P_fast <- P_fast
  }
  if (method == "lambda" && lambda_profile) {
    out$lambda_CI_lower_fast <- lambda_CI_lower_fast
    out$lambda_CI_upper_fast <- lambda_CI_upper_fast
    out$lambda_CI_cutoff_fast <- lambda_CI_cutoff_fast
    out$lambda_profile_fast <- I(lambda_profile_list)
  }

  attr(out, "match_report") <- list(
    original_tree_tips = ape::Ntip(tree),
    input_rows = nrow(X_raw),
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    tree_tips_removed = removed_tree,
    data_rows_removed = removed_data
  )
  out
}

.max_lambda <- function(tree) {
  if (!ape::is.ultrametric(tree)) {
    return(1)
  }
  depth <- ape::node.depth.edgelength(tree)
  parent_height <- depth[tree$edge[, 1]]
  child_height <- depth[tree$edge[, 2]]
  max(child_height) / max(parent_height)
}

.lambda_profile_data <- function(lik, max_lambda, lambda_hat, logL_hat,
                                 n_points) {
  grid <- seq(0, max_lambda, length.out = n_points)
  extra <- c(0, lambda_hat)
  if (max_lambda >= 1) {
    extra <- c(extra, 1)
  }
  lambda <- sort(unique(c(grid, extra)))
  logL <- vapply(lambda, lik, numeric(1))
  hit <- which.min(abs(lambda - lambda_hat))
  if (length(hit) == 1L && is.finite(logL_hat)) {
    logL[hit] <- logL_hat
  }
  data.frame(lambda = lambda, logL = logL, stringsAsFactors = FALSE)
}

.lambda_profile_ci <- function(profile, lambda_hat, logL_hat, level = 0.95) {
  profile_max <- max(c(profile$logL, logL_hat), na.rm = TRUE)
  if (!is.finite(profile_max)) {
    profile_max <- logL_hat
  }
  cutoff <- profile_max - 0.5 * stats::qchisq(level, df = 1)
  lambda <- profile$lambda
  delta <- profile$logL - cutoff
  finite <- is.finite(lambda) & is.finite(delta)
  lambda <- lambda[finite]
  delta <- delta[finite]
  if (length(lambda) < 2L || !is.finite(lambda_hat)) {
    return(c(lower = NA_real_, upper = NA_real_, cutoff = cutoff))
  }
  left <- which(lambda <= lambda_hat)
  right <- which(lambda >= lambda_hat)
  lower <- .lambda_profile_bound(lambda[left], delta[left],
                                 side = "left")
  upper <- .lambda_profile_bound(lambda[right], delta[right],
                                 side = "right")
  c(lower = lower, upper = upper, cutoff = cutoff)
}

.lambda_profile_bound <- function(lambda, delta, side = c("left", "right")) {
  side <- match.arg(side)
  if (length(lambda) < 1L) return(NA_real_)
  if (length(lambda) == 1L) return(lambda[[1L]])

  if (side == "left") {
    if (!any(delta < 0) || !any(delta >= 0)) {
      return(min(lambda))
    }
    idx <- which(delta[-length(delta)] < 0 & delta[-1L] >= 0)
    if (!length(idx)) {
      return(min(lambda))
    }
    i <- idx[[length(idx)]]
  } else {
    if (!any(delta < 0) || !any(delta >= 0)) {
      return(max(lambda))
    }
    idx <- which(delta[-length(delta)] >= 0 & delta[-1L] < 0)
    if (!length(idx)) {
      return(max(lambda))
    }
    i <- idx[[1L]]
  }

  d1 <- delta[[i]]
  d2 <- delta[[i + 1L]]
  if (d1 == d2) {
    return(mean(lambda[c(i, i + 1L)]))
  }
  lambda[[i]] + (0 - d1) * (lambda[[i + 1L]] - lambda[[i]]) / (d2 - d1)
}


# Small utilities --------------------------------------------------------------

.as_trait_matrix <- function(X, tree) {
  if (is.data.frame(X)) {
    if (!all(vapply(X, is.numeric, logical(1)))) {
      stop("All columns of x must be numeric.", call. = FALSE)
    }
    rn <- rownames(X)
    X <- data.matrix(X)
    rownames(X) <- rn
  } else if (is.null(dim(X))) {
    nm <- names(X)
    X <- matrix(as.numeric(X), ncol = 1L)
    rownames(X) <- nm
    colnames(X) <- "x"
  } else {
    if (!is.numeric(X)) {
      stop("x must be numeric.", call. = FALSE)
    }
    X <- as.matrix(X)
    storage.mode(X) <- "double"
  }

  default_rownames <- identical(rownames(X), as.character(seq_len(nrow(X))))
  if ((is.null(rownames(X)) || default_rownames) &&
      nrow(X) == ape::Ntip(tree)) {
    message("x has no names; assuming x is in tree$tip.label order")
    rownames(X) <- tree$tip.label
    default_rownames <- FALSE
  }
  if (is.null(rownames(X)) || default_rownames) {
    stop("x must have species names as names/rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(X))) {
    stop("x species names must be unique.", call. = FALSE)
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("trait_", seq_len(ncol(X)))
  }
  X
}

.permutation_matrix <- function(permutations, n, nsim,
                                include_observed = FALSE) {
  if (is.null(permutations)) {
    out <- matrix(NA_integer_, nrow = nsim, ncol = n)
    first_random <- 1L
    if (isTRUE(include_observed)) {
      out[1L, ] <- seq_len(n)
      first_random <- 2L
    }
    if (first_random <= nsim) {
      for (i in seq.int(first_random, nsim)) {
        out[i, ] <- sample.int(n)
      }
    }
    return(out)
  }
  permutations <- as.matrix(permutations)
  storage.mode(permutations) <- "integer"
  if (nrow(permutations) != nsim) {
    stop("permutations must have exactly nsim rows.", call. = FALSE)
  }
  if (ncol(permutations) != n) {
    stop("permutations must have one column per retained species.",
         call. = FALSE)
  }
  if (anyNA(permutations) || any(permutations < 1L) ||
      any(permutations > n)) {
    stop("permutations must contain 1-based indices in 1:n.",
         call. = FALSE)
  }
  duplicated_row <- vapply(
    seq_len(nrow(permutations)),
    function(i) anyDuplicated(permutations[i, ]) > 0L,
    logical(1)
  )
  if (any(duplicated_row)) {
    stop("each row of permutations must contain every index in 1:n once.",
         call. = FALSE)
  }
  permutations
}

.normalize_ncores <- function(ncores) {
  if (!is.numeric(ncores) || length(ncores) != 1L || is.na(ncores) ||
      ncores < 1) {
    stop("ncores must be a positive integer.", call. = FALSE)
  }
  as.integer(floor(ncores))
}
