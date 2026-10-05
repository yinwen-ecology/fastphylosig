d_plot_fixture <- function(table = FALSE, canonical = TRUE, denominator = -2) {
  random <- c(0, 1, 2, 2, 3)
  brownian <- c(1, 2, 2, 3, 4)
  observed <- 2
  mb <- 3
  mr <- mb + denominator
  if (table) {
    fit <- data.frame(trait='x', D_fast=(observed-mb)/denominator,
      observed=observed, mean_random=mr, mean_brownian=mb)
    fit$random_fast <- I(list(random))
    fit$brownian_fast <- I(list(brownian))
    if (canonical) { fit$P_random <- .4; fit$P_Brownian <- .2 } else {
      fit$Pval1_fast <- .4; fit$Pval0_fast <- .2
    }
  } else {
    fit <- structure(list(DEstimate=(observed-mb)/denominator, binvar='x',
      Parameters=list(Observed=observed,MeanRandom=mr,MeanBrownian=mb),
      Permutations=list(random=random,brownian=brownian)), class='phylo.d')
    if (canonical) { fit$P_random <- .4; fit$P_Brownian <- .2 } else {
      fit$Pval1 <- .4; fit$Pval0 <- .2
    }
  }
  fit
}

test_that('D selected null and probability agree for scalar/table and all aliases', {
  for (table in c(FALSE,TRUE)) for (canonical in c(FALSE,TRUE)) {
    fit <- d_plot_fixture(table,canonical)
    for (null in c('random','brownian')) {
      dat <- fastphylosig:::.signal_distribution_data(fit,null=null)
      expect_equal(dat$p_value, if (null=='random') .4 else .2)
      expect_identical(dat$null_model, paste(null,'D'))
      sums <- if (null=='random') c(0,1,2,2,3) else c(1,2,2,3,4)
      expect_equal(dat$sim[[1]], (sums-3)/(-2))
      # An explicit null takes precedence over an opposite p_col alias.
      dat <- fastphylosig:::.signal_distribution_data(fit,null=null,
        p_col=if(null=='random') 'P_Brownian' else 'P_random')
      expect_equal(dat$p_value, if(null=='random') .4 else .2)
    }
    for (alias in c('P_Brownian','Pval0','Pval0_fast')) {
      dat <- fastphylosig:::.signal_distribution_data(fit,p_col=alias)
      expect_equal(dat$p_value,.2)
      expect_identical(dat$null_model,'brownian D')
      expect_equal(fastphylosig:::.signal_plot_data(fit,p_col=alias)$p_value,.2)
    }
    for (alias in c('P_random','Pval1','Pval1_fast')) {
      expect_equal(fastphylosig:::.signal_distribution_data(fit,p_col=alias)$p_value,.4)
    }
  }
})

test_that('D overlays retain contrast-scale strict counts and fitted probabilities', {
  path <- tempfile(fileext='.pdf')
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add=TRUE)
  for (table in c(FALSE,TRUE)) for (denominator in c(-2,2)) {
    fit <- d_plot_fixture(table,canonical=FALSE,denominator=denominator)
    dat <- fastphylosig:::.signal_d_overlay_data(fit)
    # Ties at observed=2 are excluded on both original contrast tails.
    expect_equal(dat$extreme_random,2)
    expect_equal(dat$extreme_brownian,2)
    plotted <- fastphylosig::plot_signal(fit)
    expect_equal(plotted$extreme_random,2)
    expect_equal(plotted$extreme_brownian,2)
    # Fitted P can differ from counts in a partially retained simulation set.
    expect_equal(plotted$P_random_display, fastphylosig:::.format_p_value(.4))
    expect_equal(plotted$P_Brownian_display, fastphylosig:::.format_p_value(.2))
  }
  fit <- d_plot_fixture()
  # Scaling at a huge offset can erase a strict original-scale difference.
  fit$Parameters$MeanBrownian <- -1e16
  fit$Parameters$MeanRandom <- -1e16 + 2
  fit$Parameters$Observed <- 1
  fit$DEstimate <- (1 + 1e16)/2
  fit$Permutations$random <- c(0,1,2)
  fit$Permutations$brownian <- c(0,1,2)
  dat <- fastphylosig:::.signal_d_overlay_data(fit)
  expect_equal(dat$extreme_random,1)
  expect_equal(dat$extreme_brownian,1)
  fit$P_random <- 0
  fit$P_Brownian <- 0
  plotted <- fastphylosig::plot_signal(fit)
  expect_match(plotted$P_random_display,'^< ')
})
