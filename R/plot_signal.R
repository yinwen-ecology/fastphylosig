# Lightweight plotting helpers ------------------------------------------------

plot_signal <- function(result, p_col = NULL, alpha = 0.05, main = NULL,
                        xlab = NULL, ylab = NULL,
                        null = c("auto", "random", "brownian"),
                        label = c("none", "significant", "top", "all"),
                        top_n = 10, ridge_scale = 0.85, col = NULL,
                        pch = 19, cex = 0.7, ...) {
  null <- match.arg(null)
  label <- match.arg(label)
  base_dat <- .signal_plot_data(result, p_col = p_col)

  if (identical(base_dat$method[[1L]], "lambda")) {
    return(.plot_lambda_signal(
      base_dat, alpha = alpha, main = main, xlab = xlab, ylab = ylab,
      col = col, pch = pch, cex = cex, ...
    ))
  }

  if (identical(base_dat$method[[1L]], "D") &&
      nrow(base_dat) == 1L && null == "auto") {
    d_overlay <- .signal_d_overlay_data(result)
    return(.plot_d_overlay(
      d_overlay, alpha = alpha, main = main, xlab = xlab, ylab = ylab,
      col = col, ...
    ))
  }

  dat <- .signal_distribution_data(result, p_col = p_col, null = null)
  .plot_signal_ridge(
    dat, alpha = alpha, main = main, xlab = xlab, ylab = ylab,
    label = label, top_n = top_n, ridge_scale = ridge_scale,
    col = col, ...
  )
}

.plot_signal_ridge <- function(dat, alpha = 0.05, main = NULL,
                               xlab = NULL, ylab = NULL,
                               label = c("none", "significant", "top", "all"),
                               top_n = 10, ridge_scale = 0.85,
                               col = NULL, ...) {
  label <- match.arg(label)
  sim <- lapply(dat$sim, function(z) z[is.finite(z)])
  n_sim <- lengths(sim)
  ok <- is.finite(dat$estimate) & n_sim > 1L
  if (!any(ok)) {
    stop("plot_signal() needs simulated signal values. Run the fit with ",
         "test = TRUE and return_sim = TRUE.", call. = FALSE)
  }

  dat <- dat[ok, , drop = FALSE]
  sim <- sim[ok]
  n_sim <- n_sim[ok]
  n_ridge <- nrow(dat)
  sig <- !is.na(dat$p_value) & dat$p_value < alpha
  right_tail <- .is_right_tail_signal(dat$method[[1L]])
  obs_col <- "#D55E00"
  if (is.null(col)) {
    col <- if (right_tail) {
      rep("#4C78A8", n_ridge)
    } else {
      ifelse(sig, "#D55E00", "#4C78A8")
    }
  }
  if (length(col) == 1L) col <- rep(col, n_ridge)
  if (length(col) < n_ridge) col <- rep(col, length.out = n_ridge)
  if (is.null(main)) {
    main <- if (right_tail) {
      paste(dat$method[[1]], "permutation null distribution")
    } else {
      paste(dat$method[[1]], "signal ridge plot")
    }
  }
  if (is.null(xlab)) xlab <- paste(dat$method[[1]], "value")
  if (is.null(ylab)) {
    ylab <- if (n_ridge == 1L) "density" else "density, stacked by trait"
  }

  xlim <- range(c(dat$estimate, unlist(sim, use.names = FALSE)), finite = TRUE)
  if (diff(xlim) == 0) {
    xlim <- xlim + c(-0.5, 0.5)
  } else {
    pad <- 0.05 * diff(xlim)
    xlim <- xlim + c(-pad, pad)
  }

  if (n_ridge == 1L) {
    den <- stats::density(sim[[1L]], from = xlim[1], to = xlim[2], n = 256)
    graphics::plot(
      den$x, den$y, type = "n", xlim = xlim,
      ylim = c(0, max(den$y) * 1.15), xlab = xlab, ylab = ylab,
      main = main, ...
    )
    graphics::polygon(
      c(den$x, rev(den$x)), c(rep(0, length(den$x)), rev(den$y)),
      col = grDevices::adjustcolor(col[[1L]], alpha.f = 0.35),
      border = col[[1L]]
    )
    if (right_tail) {
      .draw_density_tail(
        den$x, den$y, threshold = dat$estimate[[1L]], baseline = 0,
        col = grDevices::adjustcolor(obs_col, alpha.f = 0.45)
      )
    }
    graphics::abline(v = dat$estimate[[1L]], col = obs_col, lwd = 2)
    legend <- c(
        paste0("distribution = ", dat$null_model[[1L]]),
      if (right_tail) "P tail: simulated >= observed" else NULL,
        paste0("observed ", dat$method[[1L]], " = ",
               signif(dat$estimate[[1L]], 4)),
      paste0(if (right_tail) "right-tail P = " else "P = ",
             .format_p_value(dat$p_value[[1L]])),
        paste0("nsim = ", n_sim[[1L]])
    )
    graphics::legend(
      "topright", bty = "n", legend = legend,
      fill = c(
        grDevices::adjustcolor(col[[1L]], alpha.f = 0.35),
        if (right_tail) grDevices::adjustcolor(obs_col, alpha.f = 0.45) else NULL,
        NA, NA, NA
      ),
      border = c(col[[1L]], if (right_tail) NA else NULL, NA, NA, NA),
      lwd = c(NA, if (right_tail) NA else NULL, 2, NA, NA),
      col = c(NA, if (right_tail) NA else NULL, obs_col, "black", "black")
    )
    dat$n_sim <- n_sim
    dat$tail_direction <- if (right_tail) "right" else NA_character_
    return(invisible(dat))
  }

  graphics::plot.new()
  graphics::plot.window(
    xlim = xlim, ylim = c(0.7, n_ridge + ridge_scale + 0.35), ...
  )
  graphics::axis(1)
  graphics::axis(2, at = seq_len(n_ridge), labels = dat$trait, las = 1)
  graphics::box()
  graphics::title(main = main, xlab = xlab, ylab = ylab)

  for (i in seq_len(n_ridge)) {
    z <- sim[[i]]
    graphics::segments(xlim[1], i, xlim[2], i, col = "gray80")
    if (length(unique(z)) > 1L && length(z) > 1L) {
      den <- stats::density(
        z, from = xlim[1], to = xlim[2], n = 256
      )
      height <- den$y / max(den$y) * ridge_scale
      graphics::polygon(
        c(den$x, rev(den$x)),
        c(rep(i, length(den$x)), rev(i + height)),
        col = grDevices::adjustcolor(col[[i]], alpha.f = 0.35),
        border = col[[i]]
      )
      if (right_tail) {
        .draw_density_tail(
          den$x, i + height, threshold = dat$estimate[[i]], baseline = i,
          col = grDevices::adjustcolor(obs_col, alpha.f = 0.45)
        )
      }
    }
    graphics::segments(
      dat$estimate[[i]], i, dat$estimate[[i]], i + ridge_scale,
      col = obs_col, lwd = 1.5
    )
    graphics::text(
      xlim[2], i + ridge_scale * 0.55,
      labels = paste0(
        dat$null_model[[i]], "; observed ", dat$method[[i]], "=",
        signif(dat$estimate[[i]], 4),
        if (right_tail) ", right-tail P=" else ", P=",
        .format_p_value(dat$p_value[[i]]),
        ", nsim=", n_sim[[i]]
      ),
      pos = 2, cex = 0.65, col = "gray20"
    )
  }

  label_idx <- .signal_label_indices(dat, rep(TRUE, n_ridge), sig, label, top_n)
  if (length(label_idx)) {
    graphics::text(
      dat$estimate[label_idx], label_idx + ridge_scale + 0.05,
      labels = dat$trait[label_idx], pos = 3, cex = 0.7,
      col = "gray20"
    )
  }

  dat$n_sim <- n_sim
  dat$tail_direction <- if (right_tail) "right" else NA_character_
  invisible(dat)
}

.plot_d_overlay <- function(dat, alpha = 0.05, main = NULL, xlab = NULL,
                            ylab = NULL, col = NULL, ...) {
  random <- dat$sim_random[[1L]]
  brownian <- dat$sim_brownian[[1L]]
  random <- random[is.finite(random)]
  brownian <- brownian[is.finite(brownian)]
  if (length(random) < 2L || length(brownian) < 2L) {
    stop("D plots require random and Brownian simulated D values.",
         call. = FALSE)
  }
  if (is.null(main)) main <- "Fritz & Purvis D calibration plot"
  if (is.null(xlab)) xlab <- "Fritz & Purvis D"
  if (is.null(ylab)) ylab <- "density"
  if (is.null(col)) col <- c(brownian = "#54A24B", random = "#4C78A8")
  if (!is.null(names(col)) && all(c("brownian", "random") %in% names(col))) {
    col <- col[c("brownian", "random")]
  }
  if (length(col) == 1L) col <- rep(col, 2L)
  obs_col <- "#D55E00"

  xlim <- range(c(0, 1, dat$estimate, random, brownian), finite = TRUE)
  pad <- if (diff(xlim) == 0) 0.5 else 0.05 * diff(xlim)
  xlim <- xlim + c(-pad, pad)
  den_random <- stats::density(random, from = xlim[1], to = xlim[2], n = 256)
  den_brownian <- stats::density(brownian, from = xlim[1], to = xlim[2], n = 256)
  ymax <- max(den_random$y, den_brownian$y) * 1.2
  extreme_random <- sum(random < dat$estimate[[1L]], na.rm = TRUE)
  extreme_brownian <- sum(brownian > dat$estimate[[1L]], na.rm = TRUE)
  n_random <- length(random)
  n_brownian <- length(brownian)
  p_random_label <- .format_sim_p_value(extreme_random, n_random)
  p_brownian_label <- .format_sim_p_value(extreme_brownian, n_brownian)

  graphics::plot(
    den_brownian$x, den_brownian$y, type = "n",
    xlim = xlim, ylim = c(0, ymax), xlab = xlab, ylab = ylab,
    main = main, ...
  )
  graphics::polygon(
    c(den_brownian$x, rev(den_brownian$x)),
    c(rep(0, length(den_brownian$x)), rev(den_brownian$y)),
    col = grDevices::adjustcolor(col[[1L]], alpha.f = 0.42),
    border = col[[1L]]
  )
  graphics::polygon(
    c(den_random$x, rev(den_random$x)),
    c(rep(0, length(den_random$x)), rev(den_random$y)),
    col = grDevices::adjustcolor(col[[2L]], alpha.f = 0.42),
    border = col[[2L]]
  )
  graphics::abline(v = 0, col = "gray35", lty = 3, lwd = 1.2)
  graphics::abline(v = 1, col = "gray35", lty = 3, lwd = 1.2)
  graphics::text(
    0, ymax * 0.98, labels = "D = 0\nBrownian expectation",
    pos = 4, cex = 0.7, col = "gray25"
  )
  graphics::text(
    1, ymax * 0.98, labels = "D = 1\nRandom expectation",
    pos = 4, cex = 0.7, col = "gray25"
  )
  graphics::abline(v = dat$estimate[[1L]], col = obs_col, lwd = 2)
  usr <- graphics::par("usr")
  graphics::text(
    usr[2], ymax * 0.04,
    labels = paste0(
      "Observed D = ", signif(dat$estimate[[1L]], 4),
      "\nP_random = ", p_random_label,
      "\nP_Brownian = ", p_brownian_label,
      "\nnsim = ", min(n_random, n_brownian)
    ),
    adj = c(1, 0), cex = 0.82, col = "gray15"
  )
  graphics::legend(
    "topleft", bty = "n",
    legend = c(
      "Brownian threshold null",
      "Random association null",
      "Observed D",
      "D = 0 / D = 1 references"
    ),
    fill = c(
      grDevices::adjustcolor(col[[1L]], alpha.f = 0.42),
      grDevices::adjustcolor(col[[2L]], alpha.f = 0.42),
      NA, NA
    ),
    border = c(col[[1L]], col[[2L]], NA, NA),
    lty = c(NA, NA, 1, 3),
    lwd = c(NA, NA, 2, 1.2),
    col = c(NA, NA, obs_col, "gray35")
  )
  dat$plot_type <- "D_calibration"
  dat$extreme_random <- extreme_random
  dat$extreme_brownian <- extreme_brownian
  dat$P_random_display <- p_random_label
  dat$P_Brownian_display <- p_brownian_label
  invisible(dat)
}

.plot_lambda_signal <- function(dat, alpha = 0.05, main = NULL,
                                xlab = NULL, ylab = NULL, col = NULL,
                                pch = 19, cex = 0.9, ...) {
  ok <- is.finite(dat$estimate)
  if (!any(ok)) {
    stop("No finite lambda estimates are available to plot.", call. = FALSE)
  }
  dat <- dat[ok, , drop = FALSE]
  n <- nrow(dat)
  sig <- !is.na(dat$p_value) & dat$p_value < alpha
  if (is.null(col)) col <- ifelse(sig, "#D55E00", "#4C78A8")
  if (length(col) == 1L) col <- rep(col, n)
  main_was_null <- is.null(main)
  if (is.null(xlab)) xlab <- "lambda value"
  if (is.null(ylab)) ylab <- if (n == 1L) "" else "trait"

  if (n == 1L && "profile" %in% names(dat) &&
      length(dat$profile[[1L]]) > 0L && is.data.frame(dat$profile[[1L]])) {
    return(.plot_lambda_profile(
      dat, alpha = alpha, main = if (main_was_null) NULL else main,
      xlab = xlab, ylab = ylab,
      col = col, pch = pch, cex = cex, ...
    ))
  }

  if (main_was_null) main <- "lambda likelihood-ratio summary"
  xlim <- range(c(0, 1, dat$estimate), finite = TRUE)
  pad <- if (diff(xlim) == 0) 0.1 else 0.06 * diff(xlim)
  xlim <- xlim + c(-pad, pad)

  if (n == 1L) {
    lr <- if ("LR" %in% names(dat)) dat$LR[[1L]] else NA_real_
    graphics::plot(
      NA, NA, xlim = xlim, ylim = c(0, 1), axes = FALSE,
      xlab = xlab, ylab = ylab, main = main, ...
    )
    graphics::axis(1)
    graphics::box()
    graphics::abline(v = c(0, 1), lty = 3, col = "gray70")
    graphics::segments(dat$estimate[[1L]], 0.25,
                       dat$estimate[[1L]], 0.75,
                       col = col[[1L]], lwd = 2)
    graphics::points(dat$estimate[[1L]], 0.5, pch = pch,
                     cex = cex, col = col[[1L]])
    graphics::legend(
      "topright", bty = "n",
      legend = c(
        "deterministic LR summary",
        paste0("observed lambda = ", signif(dat$estimate[[1L]], 4)),
        if (is.finite(lr)) paste0("LR = ", signif(lr, 4)) else NULL,
        paste0("P = ", .format_p_value(dat$p_value[[1L]]))
      ),
      pch = c(NA, pch, if (is.finite(lr)) NA else NULL, NA),
      col = c("black", col[[1L]],
              if (is.finite(lr)) "black" else NULL, "black")
    )
  } else {
    y <- seq_len(n)
    graphics::plot(
      NA, NA, xlim = xlim, ylim = c(0.5, n + 0.5), axes = FALSE,
      xlab = xlab, ylab = ylab, main = main, ...
    )
    graphics::axis(1)
    graphics::axis(2, at = y, labels = dat$trait, las = 1)
    graphics::box()
    graphics::abline(v = c(0, 1), lty = 3, col = "gray70")
    graphics::points(dat$estimate, y, pch = pch, cex = cex, col = col)
    graphics::text(
      xlim[2], y,
      labels = paste0(
        "lambda=", signif(dat$estimate, 4),
        if ("LR" %in% names(dat)) {
          ifelse(is.finite(dat$LR), paste0(", LR=", signif(dat$LR, 4)), "")
        } else "",
        ", P=", vapply(dat$p_value, .format_p_value, character(1))
      ),
      pos = 2, cex = 0.65, col = "gray20"
    )
  }
  dat$plot_type <- "lambda_lr_summary"
  invisible(dat)
}

.plot_lambda_profile <- function(dat, alpha = 0.05, main = NULL,
                                 xlab = NULL, ylab = NULL, col = NULL,
                                 pch = 19, cex = 0.9, ...) {
  profile <- dat$profile[[1L]]
  profile <- profile[is.finite(profile$lambda) & is.finite(profile$logL), ,
                     drop = FALSE]
  if (nrow(profile) < 2L) {
    stop("lambda profile plot requires at least two finite profile points.",
         call. = FALSE)
  }
  lambda_hat <- dat$estimate[[1L]]
  logL_hat <- if ("logL" %in% names(dat) && is.finite(dat$logL[[1L]])) {
    dat$logL[[1L]]
  } else {
    max(profile$logL)
  }
  ci_lower <- if ("ci_lower" %in% names(dat)) dat$ci_lower[[1L]] else NA_real_
  ci_upper <- if ("ci_upper" %in% names(dat)) dat$ci_upper[[1L]] else NA_real_
  cutoff <- if ("ci_cutoff" %in% names(dat) &&
      is.finite(dat$ci_cutoff[[1L]])) {
    dat$ci_cutoff[[1L]]
  } else {
    max(c(profile$logL, logL_hat), na.rm = TRUE) -
      0.5 * stats::qchisq(0.95, df = 1)
  }
  lr <- if ("LR" %in% names(dat)) dat$LR[[1L]] else NA_real_

  if (is.null(main)) main <- "lambda profile likelihood"
  if (is.null(xlab)) xlab <- "lambda value"
  if (is.null(ylab) || identical(ylab, "")) ylab <- "log-likelihood"
  if (is.null(col)) col <- "#4C78A8"
  obs_col <- "#D55E00"

  xlim <- range(c(0, 1, profile$lambda, lambda_hat, ci_lower, ci_upper),
                finite = TRUE)
  xpad <- if (diff(xlim) == 0) 0.1 else 0.04 * diff(xlim)
  xlim <- xlim + c(-xpad, xpad)
  ylim <- range(c(profile$logL, cutoff, logL_hat), finite = TRUE)
  ypad <- if (diff(ylim) == 0) 0.1 else 0.08 * diff(ylim)
  ylim <- ylim + c(-ypad, ypad)

  graphics::plot(
    profile$lambda, profile$logL, type = "l", lwd = 2,
    col = col[[1L]], xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, main = main, ...
  )
  graphics::abline(v = c(0, 1), lty = 3, col = "gray60")
  graphics::abline(h = cutoff, lty = 2, col = "gray45")
  if (is.finite(ci_lower) && is.finite(ci_upper)) {
    graphics::segments(ci_lower, cutoff, ci_upper, cutoff,
                       col = obs_col, lwd = 3)
    graphics::abline(v = c(ci_lower, ci_upper), lty = 2,
                     col = grDevices::adjustcolor(obs_col, alpha.f = 0.75))
  }
  graphics::abline(v = lambda_hat, col = obs_col, lwd = 2)
  graphics::points(lambda_hat, logL_hat, pch = pch, cex = cex,
                   col = obs_col)
  graphics::legend(
    "bottomleft", bty = "n",
    legend = c(
      "profile log-likelihood",
      paste0("lambda_hat = ", signif(lambda_hat, 4)),
      "lambda = 0 / 1",
      if (is.finite(ci_lower) && is.finite(ci_upper)) {
        paste0("approx. 95% CI = [", signif(ci_lower, 4), ", ",
               signif(ci_upper, 4), "]")
      } else {
        "approx. 95% CI unavailable"
      },
      paste0("P = ", .format_p_value(dat$p_value[[1L]])),
      if (is.finite(lr)) paste0("LR = ", signif(lr, 4)) else NULL
    ),
    lty = c(1, 1, 3, 1, NA, if (is.finite(lr)) NA else NULL),
    lwd = c(2, 2, 1, 3, NA, if (is.finite(lr)) NA else NULL),
    pch = c(NA, pch, NA, NA, NA, if (is.finite(lr)) NA else NULL),
    col = c(col[[1L]], obs_col, "gray60", obs_col, "black",
            if (is.finite(lr)) "black" else NULL)
  )

  dat$plot_type <- "lambda_profile"
  invisible(dat)
}

.format_p_value <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("0")
  format.pval(x, digits = 3, eps = .Machine$double.xmin)
}

.format_sim_p_value <- function(extreme, nsim) {
  if (!is.finite(extreme) || !is.finite(nsim) || nsim < 1) {
    return("NA")
  }
  extreme <- as.integer(extreme)
  nsim <- as.integer(nsim)
  if (extreme <= 0L) {
    return(paste0("< ", .format_p_value(1 / nsim)))
  }
  .format_p_value(extreme / nsim)
}

.is_right_tail_signal <- function(method) {
  method %in% c("K", "Delta")
}

.draw_density_tail <- function(x, y, threshold, baseline, col) {
  poly <- .density_tail_polygon(x, y, threshold, baseline)
  if (!is.null(poly)) {
    graphics::polygon(poly$x, poly$y, col = col, border = NA)
  }
  invisible(poly)
}

.density_tail_polygon <- function(x, y, threshold, baseline) {
  if (!is.finite(threshold) || length(x) < 2L || length(y) < 2L) {
    return(NULL)
  }
  keep <- x >= threshold
  if (!any(keep)) {
    return(NULL)
  }
  y_at_threshold <- stats::approx(x, y, xout = threshold, rule = 2)$y
  xt <- c(threshold, x[keep])
  yt <- c(y_at_threshold, y[keep])
  list(
    x = c(xt, rev(xt)),
    y = c(rep(baseline, length(xt)), rev(yt))
  )
}

.signal_label_indices <- function(dat, ok, sig, label, top_n) {
  if (label == "none") return(integer(0))
  ok_idx <- which(ok)
  if (label == "all") return(ok_idx)
  if (label == "significant") return(which(ok & sig))
  top_n <- min(as.integer(top_n), length(ok_idx))
  if (top_n < 1L) return(integer(0))
  ok_idx[order(dat$p_value[ok_idx], decreasing = FALSE)[seq_len(top_n)]]
}

.signal_distribution_data <- function(result, p_col = NULL,
                                      null = c("auto", "random",
                                               "brownian")) {
  null <- match.arg(null)
  dat <- .signal_plot_data(result, p_col = p_col)
  sim <- vector("list", nrow(dat))
  null_model <- rep(NA_character_, nrow(dat))

  if (is.data.frame(result)) {
    if ("K_fast" %in% names(result)) {
      if (!"sim.K_fast" %in% names(result)) {
        stop("K ridge plots require return_sim = TRUE.", call. = FALSE)
      }
      sim <- result$sim.K_fast
      null_model[] <- "permuted K"
    } else if ("D_fast" %in% names(result)) {
      d_null <- .resolve_d_null(null, p_col)
      sim_col <- if (d_null == "brownian") "brownian_fast" else "random_fast"
      if (!sim_col %in% names(result)) {
        stop("D plots require return_sim = TRUE.", call. = FALSE)
      }
      for (i in seq_len(nrow(result))) {
        denom <- result$mean_random[[i]] - result$mean_brownian[[i]]
        sim[[i]] <- (result[[sim_col]][[i]] - result$mean_brownian[[i]]) /
          denom
      }
      null_model[] <- paste0(d_null, " D")
    } else if ("Delta_fast" %in% names(result)) {
      if (!"sim.Delta_fast" %in% names(result)) {
        stop("Delta ridge plots require return_sim = TRUE.", call. = FALSE)
      }
      sim <- result$sim.Delta_fast
      null_model[] <- "permuted Delta"
    } else {
      stop("Could not identify simulated signal values for this result.",
           call. = FALSE)
    }
  } else if (inherits(result, "phylo.d")) {
    if (is.null(result$Permutations)) {
      stop("D plots require return_sim = TRUE.", call. = FALSE)
    }
    d_null <- .resolve_d_null(null, p_col)
    sums <- result$Permutations[[d_null]]
    denom <- result$Parameters$MeanRandom - result$Parameters$MeanBrownian
    sim[[1L]] <- (sums - result$Parameters$MeanBrownian) / denom
    null_model[[1L]] <- paste0(d_null, " D")
  } else if (inherits(result, "phylo_delta")) {
    if (is.null(result$sim.delta)) {
      stop("Delta ridge plots require return_sim = TRUE.", call. = FALSE)
    }
    sim[[1L]] <- result$sim.delta
    null_model[[1L]] <- "permuted Delta"
  } else if (inherits(result, "phylosig")) {
    if (!is.list(result) || is.null(result$sim.K)) {
      stop("K ridge plots require test = TRUE and return_sim = TRUE.",
           call. = FALSE)
    }
    sim[[1L]] <- result$sim.K
    null_model[[1L]] <- "permuted K"
  } else {
    stop("result must be produced by fast_signal(), fast_d(), or fast_delta().",
         call. = FALSE)
  }

  dat$sim <- I(lapply(sim, as.numeric))
  dat$null_model <- null_model
  dat
}

.signal_d_overlay_data <- function(result) {
  dat <- .signal_plot_data(result, p_col = NULL)
  if (!identical(dat$method[[1L]], "D") || nrow(dat) != 1L) {
    stop("D overlay plots require a single D result.", call. = FALSE)
  }

  if (is.data.frame(result)) {
    required <- c(
      "random_fast", "brownian_fast", "mean_random", "mean_brownian",
      "Pval1_fast", "Pval0_fast"
    )
    missing <- setdiff(required, names(result))
    if (length(missing)) {
      stop("D overlay plots require return_sim = TRUE.", call. = FALSE)
    }
    mean_random <- result$mean_random[[1L]]
    mean_brownian <- result$mean_brownian[[1L]]
    random_sums <- result$random_fast[[1L]]
    brownian_sums <- result$brownian_fast[[1L]]
    p_random <- result$Pval1_fast[[1L]]
    p_brownian <- result$Pval0_fast[[1L]]
  } else if (inherits(result, "phylo.d")) {
    if (is.null(result$Permutations)) {
      stop("D overlay plots require return_sim = TRUE.", call. = FALSE)
    }
    mean_random <- result$Parameters$MeanRandom
    mean_brownian <- result$Parameters$MeanBrownian
    random_sums <- result$Permutations$random
    brownian_sums <- result$Permutations$brownian
    p_random <- result$Pval1
    p_brownian <- result$Pval0
  } else {
    stop("D overlay plots require a fast_d() result.", call. = FALSE)
  }

  denom <- mean_random - mean_brownian
  if (!is.finite(denom) || denom == 0) {
    stop("D null distributions cannot be scaled because the random and ",
         "Brownian means are identical or non-finite.", call. = FALSE)
  }
  random_d <- (as.numeric(random_sums) - mean_brownian) / denom
  brownian_d <- (as.numeric(brownian_sums) - mean_brownian) / denom

  out <- data.frame(
    trait = dat$trait[[1L]],
    method = "D",
    estimate = dat$estimate[[1L]],
    p_value = dat$p_value[[1L]],
    Pval1 = p_random,
    Pval0 = p_brownian,
    n_sim = length(random_d),
    stringsAsFactors = FALSE
  )
  out$sim_random <- I(list(random_d))
  out$sim_brownian <- I(list(brownian_d))
  out
}

.resolve_d_null <- function(null, p_col) {
  if (null != "auto") return(null)
  if (!is.null(p_col) && grepl("Pval0", p_col, fixed = TRUE)) {
    return("brownian")
  }
  "random"
}

.signal_plot_data <- function(result, p_col = NULL) {
  if (is.data.frame(result)) {
    trait <- if ("trait" %in% names(result)) {
      as.character(result$trait)
    } else {
      as.character(seq_len(nrow(result)))
    }
    if ("K_fast" %in% names(result)) {
      method <- "K"
      estimate <- result$K_fast
      default_p <- "P_fast"
    } else if ("lambda_fast" %in% names(result)) {
      method <- "lambda"
      estimate <- result$lambda_fast
      default_p <- "P_fast"
    } else if ("D_fast" %in% names(result)) {
      method <- "D"
      estimate <- result$D_fast
      default_p <- "Pval1_fast"
    } else if ("Delta_fast" %in% names(result)) {
      method <- "Delta"
      estimate <- result$Delta_fast
      default_p <- "P_fast"
    } else {
      stop("Could not identify a supported signal estimate column.",
           call. = FALSE)
    }
    p_name <- if (is.null(p_col)) default_p else p_col
    p_value <- if (p_name %in% names(result)) result[[p_name]] else NA_real_
    out <- data.frame(
      trait = trait,
      method = method,
      estimate = as.numeric(estimate),
      p_value = as.numeric(p_value),
      stringsAsFactors = FALSE
    )
    if (identical(method, "lambda")) {
      out$logL <- if ("logL_fast" %in% names(result)) {
        as.numeric(result$logL_fast)
      } else {
        NA_real_
      }
      out$logL0 <- if ("logL0_fast" %in% names(result)) {
        as.numeric(result$logL0_fast)
      } else {
        NA_real_
      }
      out$LR <- ifelse(
        is.finite(out$logL) & is.finite(out$logL0),
        2 * (out$logL - out$logL0),
        NA_real_
      )
      if ("lambda_profile_fast" %in% names(result)) {
        out$profile <- I(result$lambda_profile_fast)
        out$ci_lower <- if ("lambda_CI_lower_fast" %in% names(result)) {
          as.numeric(result$lambda_CI_lower_fast)
        } else {
          NA_real_
        }
        out$ci_upper <- if ("lambda_CI_upper_fast" %in% names(result)) {
          as.numeric(result$lambda_CI_upper_fast)
        } else {
          NA_real_
        }
        out$ci_cutoff <- if ("lambda_CI_cutoff_fast" %in% names(result)) {
          as.numeric(result$lambda_CI_cutoff_fast)
        } else {
          NA_real_
        }
      }
    }
    return(out)
  }

  if (inherits(result, "phylo.d")) {
    return(data.frame(
      trait = result$binvar %||% "x",
      method = "D",
      estimate = as.numeric(result$DEstimate),
      p_value = as.numeric(result$Pval1),
      stringsAsFactors = FALSE
    ))
  }

  if (inherits(result, "phylo_delta")) {
    return(data.frame(
      trait = "x",
      method = "Delta",
      estimate = as.numeric(result$delta),
      p_value = if (!is.null(result$P)) as.numeric(result$P) else NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  if (inherits(result, "phylosig")) {
    method <- attr(result, "method")
    if (is.null(method)) method <- "signal"
    if (is.list(result)) {
      estimate <- if (!is.null(result$K)) result$K else result$lambda
      p_value <- if (!is.null(result$P)) result$P else NA_real_
    } else {
      estimate <- as.numeric(result)
      p_value <- NA_real_
    }
    out <- data.frame(
      trait = "x",
      method = method,
      estimate = as.numeric(estimate),
      p_value = as.numeric(p_value),
      stringsAsFactors = FALSE
    )
    if (identical(method, "lambda") && is.list(result)) {
      out$logL <- if (!is.null(result$logL)) as.numeric(result$logL) else NA_real_
      out$logL0 <- if (!is.null(result$logL0)) as.numeric(result$logL0) else NA_real_
      out$LR <- ifelse(
        is.finite(out$logL) & is.finite(out$logL0),
        2 * (out$logL - out$logL0),
        NA_real_
      )
      if (!is.null(result$lambda_profile)) {
        out$profile <- I(list(result$lambda_profile))
        out$ci_lower <- if (!is.null(result$lambda_CI)) {
          as.numeric(result$lambda_CI[["lower"]])
        } else {
          NA_real_
        }
        out$ci_upper <- if (!is.null(result$lambda_CI)) {
          as.numeric(result$lambda_CI[["upper"]])
        } else {
          NA_real_
        }
        out$ci_cutoff <- if (!is.null(result$lambda_CI_cutoff)) {
          as.numeric(result$lambda_CI_cutoff)
        } else {
          NA_real_
        }
      }
    }
    return(out)
  }

  stop("result must be produced by fast_signal(), fast_d(), or fast_delta().",
       call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}
