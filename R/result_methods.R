# Stable result classes and display methods -----------------------------------

# The numerical engines intentionally keep their historical return values
# (including names used by phytools/caper callers).  This file is a thin
# compatibility layer: it adds a small set of common aliases and a shared S3
# class without changing any estimate or simulation value.

.result_first <- function(x, fields, default = NA_real_) {
  if (is.null(x) || !length(fields)) return(default)
  for (field in fields) {
    value <- NULL
    if (is.list(x)) value <- x[[field]]
    if (is.null(value)) value <- attr(x, field, exact = TRUE)
    if (!is.null(value) && length(value)) return(value)
  }
  default
}

.result_scalar <- function(value, default = NA_real_) {
  if (is.null(value) || !length(value)) return(default)
  value[[1L]]
}

.result_set <- function(x, name, value) {
  if (is.null(x[[name]])) x[[name]] <- value
  x
}

.result_table_value <- function(x, fields, default = NA_real_) {
  n <- nrow(x)
  for (field in fields) {
    if (field %in% names(x)) return(x[[field]])
  }
  rep(default, n)
}

.result_table_set <- function(x, name, value) {
  if (!(name %in% names(x))) x[[name]] <- value
  x
}

.result_method <- function(x, method = NULL) {
  if (!is.null(method) && length(method)) {
    method <- as.character(method[[1L]])
  } else {
    method <- attr(x, "method", exact = TRUE)
    if (!is.null(method) && length(method)) method <- as.character(method[[1L]])
  }
  if (!is.null(method) && length(method)) {
    key <- tolower(method[[1L]])
    if (key %in% c("k", "blomberg", "blomberg's k")) return("K")
    if (key %in% c("lambda", "pagel", "pagel's lambda")) return("lambda")
    if (key %in% c("d", "phylo.d", "phylo_d")) return("D")
    if (key %in% c("delta", "phylo_delta", "phylo.delta")) return("Delta")
  }

  cls <- class(x)
  if (any(c("fastphylosig_d", "phylo.d") %in% cls)) return("D")
  if (any(c("fastphylosig_delta", "phylo_delta") %in% cls)) return("Delta")
  if ("fastphylosig_signal" %in% cls) {
    # K is the only signal result that can be an atomic numeric.  Lists carry
    # an explicit method attribute in all public constructors; when a caller
    # strips that attribute, inspect the legacy fields as a fallback.
    if (is.list(x) && any(c("lambda", "logL", "logLik") %in% names(x))) {
      return("lambda")
    }
    return("K")
  }
  if ("phylosig" %in% cls) {
    if (is.list(x) && any(c("lambda", "logL", "logLik") %in% names(x))) {
      return("lambda")
    }
    return("K")
  }
  if (any(c("phylo_delta") %in% cls)) return("Delta")
  if (any(c("phylo.d") %in% cls)) return("D")
  if (is.list(x)) {
    nms <- names(x)
    if (any(c("DEstimate", "Pval1", "Pval0", "P_random", "p_random") %in% nms)) {
      return("D")
    }
    if (any(c("Delta_fast", "MCSE_Delta", "ESS_alpha", "estimate_mcse") %in% nms)) {
      return("Delta")
    }
    if (any(c("lambda", "logL", "logLik") %in% nms)) return("lambda")
    if (any(c("K", "K_fast") %in% nms)) return("K")
  }
  if (is.data.frame(x)) {
    nms <- names(x)
    if (any(c("D_fast", "Pval1_fast", "P_random", "p_random") %in% nms)) {
      return("D")
    }
    if (any(c("Delta_fast", "MCSE_Delta", "estimate_mcse") %in% nms)) {
      return("Delta")
    }
    if (any(c("lambda_fast", "logL_fast", "logLik") %in% nms)) {
      return("lambda")
    }
    if (any(c("K_fast", "K", "p_mcse") %in% nms)) return("K")
  }
  NULL
}

.result_special_class <- function(method) {
  switch(method,
    K = "fastphylosig_signal",
    lambda = "fastphylosig_signal",
    D = "fastphylosig_d",
    Delta = "fastphylosig_delta",
    "fastphylosig_result"
  )
}

.result_set_classes <- function(x, method, vector_input) {
  old <- attr(x, "class", exact = TRUE)
  if (is.null(old)) old <- character()
  old <- as.character(old)
  if (!isTRUE(vector_input)) {
    if (!length(old)) old <- "data.frame"
    if (!any(old == "data.frame")) old <- c(old, "data.frame")
    class(x) <- unique(c("fastphylosig_table", "fastphylosig_result", old))
  } else {
    class(x) <- unique(c(.result_special_class(method),
                         "fastphylosig_result", old))
  }
  x
}

.decorate_signal_vector <- function(x, method) {
  if (!is.list(x) || is.data.frame(x)) return(x)

  if (identical(method, "K")) {
    estimate <- .result_first(x, c("estimate", "K", "K_fast"))
    p_value <- .result_first(x, c("p_value", "P", "P_fast"))
    p_mcse <- .result_first(x, c("p_mcse", "MCSE_P", "P_MCSE"))
    n_requested <- .result_first(
      x, c("n_sim_requested", "nsim_requested", "nsim")
    )
    n_successful <- .result_first(
      x, c("n_sim_successful", "nsim_successful")
    )
    x <- .result_set(x, "estimate", .result_scalar(estimate))
    x <- .result_set(x, "p_value", .result_scalar(p_value))
    x <- .result_set(x, "p_mcse", .result_scalar(p_mcse))
    x <- .result_set(x, "n_sim_requested", .result_scalar(n_requested))
    x <- .result_set(x, "n_sim_successful", .result_scalar(n_successful))
    x <- .result_set(x, "trait", "x")
  } else if (identical(method, "lambda")) {
    estimate <- .result_first(x, c("estimate", "lambda", "lambda_fast"))
    log_lik <- .result_first(x, c("logLik", "logL", "logL_fast"))
    lr <- .result_first(x, c("LR", "LR_fast"))
    p_value <- .result_first(x, c("p_value", "P", "P_fast"))
    x <- .result_set(x, "estimate", .result_scalar(estimate))
    x <- .result_set(x, "logLik", .result_scalar(log_lik))
    x <- .result_set(x, "LR", .result_scalar(lr))
    x <- .result_set(x, "p_value", .result_scalar(p_value))
    x <- .result_set(x, "trait", "x")
  }
  x
}

.decorate_d_vector <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  estimate <- .result_first(x, c("estimate", "DEstimate", "D", "D_fast"))
  p_random <- .result_first(
    x, c("P_random", "p_random", "Pval1", "Pval1_fast")
  )
  p_brownian <- .result_first(
    x, c("P_Brownian", "p_brownian", "Pval0", "Pval0_fast")
  )
  mcse_random <- .result_first(
    x, c("MCSE_P_random", "p_mcse_random", "mcse_p_random")
  )
  mcse_brownian <- .result_first(
    x, c("MCSE_P_Brownian", "p_mcse_brownian", "mcse_p_brownian")
  )
  n_requested <- .result_first(
    x, c("n_sim_requested", "nsim_requested", "nPermut")
  )
  n_random <- .result_first(
    x, c("n_sim_successful_random", "nsim_successful_random")
  )
  n_brownian <- .result_first(
    x, c("n_sim_successful_brownian", "nsim_successful_brownian")
  )

  x <- .result_set(x, "estimate", .result_scalar(estimate))
  # Keep the two null-specific P values separate.  In particular, do not add
  # a p_value field for D: no scientifically meaningful combined test exists.
  x <- .result_set(x, "P_random", .result_scalar(p_random))
  x <- .result_set(x, "P_Brownian", .result_scalar(p_brownian))
  x <- .result_set(x, "p_random", .result_scalar(p_random))
  x <- .result_set(x, "p_brownian", .result_scalar(p_brownian))
  x <- .result_set(x, "MCSE_P_random", .result_scalar(mcse_random))
  x <- .result_set(x, "MCSE_P_Brownian", .result_scalar(mcse_brownian))
  x <- .result_set(x, "p_mcse_random", .result_scalar(mcse_random))
  x <- .result_set(x, "p_mcse_brownian", .result_scalar(mcse_brownian))
  x <- .result_set(x, "MCSE_random", .result_scalar(mcse_random))
  x <- .result_set(x, "MCSE_Brownian", .result_scalar(mcse_brownian))
  x <- .result_set(x, "mcse_random", .result_scalar(mcse_random))
  x <- .result_set(x, "mcse_brownian", .result_scalar(mcse_brownian))
  x <- .result_set(x, "n_sim_requested", .result_scalar(n_requested))
  x <- .result_set(x, "n_sim_successful_random", .result_scalar(n_random))
  x <- .result_set(x, "n_sim_successful_brownian", .result_scalar(n_brownian))
  x <- .result_set(x, "trait", .result_scalar(.result_first(
    x, c("trait", "binvar"), default = "x"
  ), "x"))
  x
}

.decorate_delta_vector <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  estimate <- .result_first(x, c("estimate", "delta", "Delta", "Delta_fast"))
  p_value <- .result_first(x, c("p_value", "P", "P_fast"))
  p_mcse <- .result_first(x, c("p_mcse", "P_MCSE", "P_mcse", "MCSE_P"))
  estimate_mcse <- .result_first(
    x, c("estimate_mcse", "MCSE_Delta", "delta_mcse", "mc_se")
  )
  ess_alpha <- .result_first(x, c("ESS_alpha", "ess_alpha"))
  ess_beta <- .result_first(x, c("ESS_beta", "ess_beta"))
  rhat_alpha <- .result_first(
    x, c("Rhat_alpha", "split_Rhat_alpha", "rhat_alpha")
  )
  rhat_beta <- .result_first(
    x, c("Rhat_beta", "split_Rhat_beta", "rhat_beta")
  )
  n_sim_requested <- .result_first(
    x, c("n_sim_requested", "requested_simulations",
         "requested_permutations", "nPermut")
  )
  n_sim_successful <- .result_first(
    x, c("n_sim_successful", "successful_simulations",
         "successful_permutations")
  )
  n_iter_requested <- .result_first(x, c("n_iter_requested", "requested_iterations"))
  n_iter_successful <- .result_first(x, c("n_iter_successful", "successful_iterations"))

  x <- .result_set(x, "estimate", .result_scalar(estimate))
  x <- .result_set(x, "p_value", .result_scalar(p_value))
  x <- .result_set(x, "p_mcse", .result_scalar(p_mcse))
  x <- .result_set(x, "estimate_mcse", .result_scalar(estimate_mcse))
  x <- .result_set(x, "ess_alpha", .result_scalar(ess_alpha))
  x <- .result_set(x, "ess_beta", .result_scalar(ess_beta))
  x <- .result_set(x, "rhat_alpha", .result_scalar(rhat_alpha))
  x <- .result_set(x, "rhat_beta", .result_scalar(rhat_beta))
  x <- .result_set(x, "ESS_alpha", .result_scalar(ess_alpha))
  x <- .result_set(x, "ESS_beta", .result_scalar(ess_beta))
  x <- .result_set(x, "Rhat_alpha", .result_scalar(rhat_alpha))
  x <- .result_set(x, "Rhat_beta", .result_scalar(rhat_beta))
  # ESS/R-hat are chain-specific diagnostics.  The aggregate aliases are
  # conservative summaries (minimum ESS and maximum R-hat), not new model
  # estimates, and are useful for compact reporting.
  if (is.null(x[["ESS"]])) {
    x[["ESS"]] <- suppressWarnings(min(
      c(.result_scalar(ess_alpha), .result_scalar(ess_beta)), na.rm = TRUE
    ))
    if (!is.finite(x[["ESS"]])) x[["ESS"]] <- NA_real_
  }
  x <- .result_set(x, "ess", x[["ESS"]])
  if (is.null(x[["Rhat"]])) {
    x[["Rhat"]] <- suppressWarnings(max(
      c(.result_scalar(rhat_alpha), .result_scalar(rhat_beta)), na.rm = TRUE
    ))
    if (!is.finite(x[["Rhat"]])) x[["Rhat"]] <- NA_real_
  }
  x <- .result_set(x, "rhat", x[["Rhat"]])
  x <- .result_set(x, "n_sim_requested", .result_scalar(n_sim_requested))
  x <- .result_set(x, "n_sim_successful", .result_scalar(n_sim_successful))
  x <- .result_set(x, "n_iter_requested", .result_scalar(n_iter_requested))
  x <- .result_set(x, "n_iter_successful", .result_scalar(n_iter_successful))
  x <- .result_set(x, "trait", .result_scalar(
    .result_first(x, c("trait"), default = "x"), "x"
  ))
  x
}

.decorate_signal_table <- function(x, method) {
  estimate <- .result_table_value(x, if (identical(method, "K")) {
    c("estimate", "K_fast", "K")
  } else {
    c("estimate", "lambda_fast", "lambda")
  })
  x <- .result_table_set(x, "estimate", estimate)
  if (identical(method, "K")) {
    p_value <- .result_table_value(x, c("p_value", "P_fast", "P"))
    p_mcse <- .result_table_value(x, c("p_mcse", "MCSE_P", "P_MCSE"))
    n_requested <- .result_table_value(
      x, c("n_sim_requested", "nsim_requested", "nsim")
    )
    n_successful <- .result_table_value(
      x, c("n_sim_successful", "nsim_successful")
    )
    x <- .result_table_set(x, "p_value", p_value)
    x <- .result_table_set(x, "p_mcse", p_mcse)
    x <- .result_table_set(x, "n_sim_requested", n_requested)
    x <- .result_table_set(x, "n_sim_successful", n_successful)
  } else {
    x <- .result_table_set(x, "logLik", .result_table_value(
      x, c("logLik", "logL_fast", "logL")
    ))
    x <- .result_table_set(x, "LR", .result_table_value(x, c("LR", "LR_fast")))
    x <- .result_table_set(x, "p_value", .result_table_value(
      x, c("p_value", "P_fast", "P")
    ))
  }
  x
}

.decorate_d_table <- function(x) {
  estimate <- .result_table_value(x, c("estimate", "D_fast", "D"))
  p_random <- .result_table_value(
    x, c("P_random", "p_random", "Pval1", "Pval1_fast")
  )
  p_brownian <- .result_table_value(
    x, c("P_Brownian", "p_brownian", "Pval0", "Pval0_fast")
  )
  mcse_random <- .result_table_value(
    x, c("MCSE_P_random", "p_mcse_random", "mcse_p_random")
  )
  mcse_brownian <- .result_table_value(
    x, c("MCSE_P_Brownian", "p_mcse_brownian", "mcse_p_brownian")
  )
  x <- .result_table_set(x, "estimate", estimate)
  x <- .result_table_set(x, "P_random", p_random)
  x <- .result_table_set(x, "P_Brownian", p_brownian)
  x <- .result_table_set(x, "p_random", p_random)
  x <- .result_table_set(x, "p_brownian", p_brownian)
  x <- .result_table_set(x, "MCSE_P_random", mcse_random)
  x <- .result_table_set(x, "MCSE_P_Brownian", mcse_brownian)
  x <- .result_table_set(x, "p_mcse_random", mcse_random)
  x <- .result_table_set(x, "p_mcse_brownian", mcse_brownian)
  x <- .result_table_set(x, "MCSE_random", mcse_random)
  x <- .result_table_set(x, "MCSE_Brownian", mcse_brownian)
  x <- .result_table_set(x, "mcse_random", mcse_random)
  x <- .result_table_set(x, "mcse_brownian", mcse_brownian)
  x <- .result_table_set(x, "n_sim_requested", .result_table_value(
    x, c("n_sim_requested", "nsim_requested")
  ))
  x <- .result_table_set(x, "n_sim_successful_random", .result_table_value(
    x, c("n_sim_successful_random", "nsim_successful_random")
  ))
  x <- .result_table_set(x, "n_sim_successful_brownian", .result_table_value(
    x, c("n_sim_successful_brownian", "nsim_successful_brownian")
  ))
  x
}

.decorate_delta_table <- function(x) {
  estimate <- .result_table_value(x, c("estimate", "Delta_fast", "delta"))
  p_value <- .result_table_value(x, c("p_value", "P_fast", "P"))
  p_mcse <- .result_table_value(x, c("p_mcse", "P_MCSE", "P_mcse", "MCSE_P"))
  estimate_mcse <- .result_table_value(
    x, c("estimate_mcse", "MCSE_Delta", "delta_mcse", "mc_se")
  )
  ess_alpha <- .result_table_value(x, c("ESS_alpha", "ess_alpha"))
  ess_beta <- .result_table_value(x, c("ESS_beta", "ess_beta"))
  rhat_alpha <- .result_table_value(
    x, c("Rhat_alpha", "split_Rhat_alpha", "rhat_alpha")
  )
  rhat_beta <- .result_table_value(
    x, c("Rhat_beta", "split_Rhat_beta", "rhat_beta")
  )
  x <- .result_table_set(x, "estimate", estimate)
  x <- .result_table_set(x, "p_value", p_value)
  x <- .result_table_set(x, "p_mcse", p_mcse)
  x <- .result_table_set(x, "estimate_mcse", estimate_mcse)
  x <- .result_table_set(x, "ess_alpha", ess_alpha)
  x <- .result_table_set(x, "ess_beta", ess_beta)
  x <- .result_table_set(x, "rhat_alpha", rhat_alpha)
  x <- .result_table_set(x, "rhat_beta", rhat_beta)
  x <- .result_table_set(x, "ESS_alpha", ess_alpha)
  x <- .result_table_set(x, "ESS_beta", ess_beta)
  x <- .result_table_set(x, "Rhat_alpha", rhat_alpha)
  x <- .result_table_set(x, "Rhat_beta", rhat_beta)
  if (!"ESS" %in% names(x)) {
    x$ESS <- apply(cbind(ess_alpha, ess_beta), 1L, function(z) {
      z <- z[is.finite(z)]
      if (!length(z)) NA_real_ else min(z)
    })
  }
  x <- .result_table_set(x, "ess", x$ESS)
  if (!"Rhat" %in% names(x)) {
    x$Rhat <- apply(cbind(rhat_alpha, rhat_beta), 1L, function(z) {
      z <- z[is.finite(z)]
      if (!length(z)) NA_real_ else max(z)
    })
  }
  x <- .result_table_set(x, "rhat", x$Rhat)
  x <- .result_table_set(x, "n_sim_requested", .result_table_value(
    x, c("n_sim_requested", "requested_simulations", "requested_permutations")
  ))
  x <- .result_table_set(x, "n_sim_successful", .result_table_value(
    x, c("n_sim_successful", "successful_simulations", "successful_permutations")
  ))
  x <- .result_table_set(x, "n_iter_requested", .result_table_value(
    x, c("n_iter_requested", "requested_iterations")
  ))
  x <- .result_table_set(x, "n_iter_successful", .result_table_value(
    x, c("n_iter_successful", "successful_iterations")
  ))
  x
}

#' Decorate a fastphylosig result with stable aliases and S3 classes.
#'
#' This is intentionally internal.  Public engines call it after constructing
#' their legacy object so numerical code remains independent of presentation.
.decorate_fastphylosig_result <- function(x, method = NULL,
                                           vector_input = NULL) {
  if (is.null(vector_input)) vector_input <- !is.data.frame(x)
  vector_input <- isTRUE(vector_input)
  method <- .result_method(x, method)
  if (is.null(method)) {
    stop("method must identify K, lambda, D, or Delta.", call. = FALSE)
  }

  if (vector_input) {
    if (identical(method, "K") && is.numeric(x) && !is.list(x)) {
      # A no-test K result is historically an atomic numeric with class
      # `phylosig`.  Keep it atomic; print/summary/data.frame derive estimate
      # directly from its value rather than forcing a list representation.
      attr(x, "method") <- if (is.null(attr(x, "method", exact = TRUE))) {
        method
      } else attr(x, "method", exact = TRUE)
    } else if (is.list(x) && !is.data.frame(x)) {
      x <- switch(method,
        K = .decorate_signal_vector(x, method),
        lambda = .decorate_signal_vector(x, method),
        D = .decorate_d_vector(x),
        Delta = .decorate_delta_vector(x)
      )
      if (is.null(x[["method"]])) x[["method"]] <- method
    }
  } else {
    if (!is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
    x <- switch(method,
      K = .decorate_signal_table(x, method),
      lambda = .decorate_signal_table(x, method),
      D = .decorate_d_table(x),
      Delta = .decorate_delta_table(x)
    )
    if (!"method" %in% names(x)) x$method <- method
    attr(x, "method") <- method
  }
  attr(x, "result_method") <- method
  attr(x, "vector_input") <- vector_input
  .result_set_classes(x, method, vector_input)
}

.result_simulation_columns <- function(nms) {
  nms[grepl(
    "^(sim([.]|_)|random_fast$|brownian_fast$|Permutations$|sim\\.)",
    nms, ignore.case = TRUE
  )]
}

.result_display_row <- function(x, method) {
  if (is.data.frame(x)) {
    out <- x
    drop <- unique(c(
      .result_simulation_columns(names(out)),
      names(out)[vapply(out, is.list, logical(1))]
    ))
    if (length(drop)) out[drop] <- NULL
    class(out) <- "data.frame"
    return(out)
  }

  # Pull scalar fields through the same aliases used by the decorator.  The
  # output is deliberately one row and never contains simulation vectors.
  get <- function(fields, default = NA_real_) {
    .result_scalar(.result_first(x, fields, default), default)
  }
  common <- list(
    trait = .result_scalar(
      .result_first(x, c("trait", "binvar"), default = "x"), "x"
    ),
    method = method,
    estimate = if (identical(method, "K") && is.numeric(x) && !is.list(x)) {
      .result_scalar(as.numeric(x))
    } else get(c("estimate", "K", "lambda", "DEstimate", "delta")),
    status = get(c("status"), NA_character_),
    message = get(c("message", "note", "diagnostics_note"), NA_character_),
    note = get(c("note", "diagnostics_note"), NA_character_)
  )
  if (identical(method, "K")) {
    common$p_value <- get(c("p_value", "P"))
    common$p_mcse <- get(c("p_mcse", "MCSE_P", "P_MCSE"))
    common$n_sim_requested <- get(c("n_sim_requested", "nsim_requested", "nsim"))
    common$n_sim_successful <- get(c("n_sim_successful", "nsim_successful"))
  } else if (identical(method, "lambda")) {
    common$logLik <- get(c("logLik", "logL"))
    common$LR <- get(c("LR", "LR_fast"))
    common$p_value <- get(c("p_value", "P"))
  } else if (identical(method, "D")) {
    common$P_random <- get(c("P_random", "Pval1"))
    common$P_Brownian <- get(c("P_Brownian", "Pval0"))
    common$MCSE_P_random <- get(c("MCSE_P_random", "p_mcse_random"))
    common$MCSE_P_Brownian <- get(c("MCSE_P_Brownian", "p_mcse_brownian"))
    common$n_sim_requested <- get(c("n_sim_requested", "nsim_requested", "nPermut"))
    common$n_sim_successful_random <- get(
      c("n_sim_successful_random", "nsim_successful_random")
    )
    common$n_sim_successful_brownian <- get(
      c("n_sim_successful_brownian", "nsim_successful_brownian")
    )
  } else if (identical(method, "Delta")) {
    common$p_value <- get(c("p_value", "P"))
    common$p_mcse <- get(c("p_mcse", "P_MCSE", "P_mcse", "MCSE_P"))
    common$estimate_mcse <- get(c("estimate_mcse", "MCSE_Delta"))
    common$ESS_alpha <- get(c("ESS_alpha", "ess_alpha"))
    common$ESS_beta <- get(c("ESS_beta", "ess_beta"))
    common$Rhat_alpha <- get(c("Rhat_alpha", "split_Rhat_alpha", "rhat_alpha"))
    common$Rhat_beta <- get(c("Rhat_beta", "split_Rhat_beta", "rhat_beta"))
    common$ESS <- get(c("ESS", "ess"))
    common$Rhat <- get(c("Rhat", "rhat"))
    common$n_sim_requested <- get(
      c("n_sim_requested", "requested_simulations", "requested_permutations")
    )
    common$n_sim_successful <- get(
      c("n_sim_successful", "successful_simulations", "successful_permutations")
    )
    common$n_saved_requested <- get(c("n_saved_requested"))
    common$n_saved_successful <- get(c("n_saved_successful"))
    common$n_iter_requested <- get(c("n_iter_requested", "requested_iterations"))
    common$n_iter_successful <- get(c("n_iter_successful", "successful_iterations"))
  }
  as.data.frame(common, stringsAsFactors = FALSE, check.names = FALSE)
}

.result_compact_frame <- function(x, method) {
  out <- .result_display_row(x, method)
  if (!is.data.frame(x)) return(out)
  wanted <- switch(method,
    K = c("trait", "method", "estimate", "p_value", "p_mcse",
          "n_sim_requested", "n_sim_successful", "n_species",
          "n_removed_na", "status", "message", "note"),
    lambda = c("trait", "method", "estimate", "logLik", "LR", "p_value",
               "n_species", "n_removed_na", "status", "message", "note"),
    D = c("trait", "method", "estimate", "P_random", "P_Brownian",
          "MCSE_P_random", "MCSE_P_Brownian", "n_sim_requested",
          "n_species", "n_removed_na", "status", "message", "note"),
    Delta = c("trait", "method", "estimate", "p_value", "p_mcse",
              "estimate_mcse", "ess_alpha", "ess_beta", "rhat_alpha",
              "rhat_beta", "n_sim_requested", "n_sim_successful",
              "n_species", "n_removed_na", "status", "message", "note"),
    names(out)
  )
  out[intersect(wanted, names(out))]
}

#' @export
print.fastphylosig_result <- function(x, ..., digits = getOption("digits")) {
  method <- .result_method(x)
  if (is.null(method)) method <- attr(x, "result_method", exact = TRUE)
  out <- .result_compact_frame(x, method)
  cat(sprintf("fastphylosig %s result\n", method %||% "result"))
  print.data.frame(out, ..., digits = digits, row.names = FALSE)
  invisible(x)
}

# `%||%` is kept local to avoid importing a dependency for a two-token helper.
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

#' @export
summary.fastphylosig_result <- function(object, ...) {
  method <- .result_method(object)
  out <- .result_compact_frame(object, method)
  class(out) <- c("summary.fastphylosig_result", "data.frame")
  attr(out, "method") <- method
  out
}

#' @export
print.summary.fastphylosig_result <- function(x, ..., digits = getOption("digits")) {
  out <- x
  class(out) <- "data.frame"
  print.data.frame(out, ..., digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.fastphylosig_result <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  method <- .result_method(x)
  out <- .result_display_row(x, method)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

# A table inherits these methods in normal S3 dispatch, but explicit wrappers
# make the registrations robust to classes supplied by downstream packages.
#' @export
print.fastphylosig_table <- print.fastphylosig_result

#' @export
summary.fastphylosig_table <- summary.fastphylosig_result

#' @export
as.data.frame.fastphylosig_table <- as.data.frame.fastphylosig_result
