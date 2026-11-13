.with_review_plot_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

test_that("K scalar and ridge plots accept explicit graphics limits", {
  set.seed(91)
  tree <- ape::rtree(12)
  x <- stats::setNames(stats::rnorm(12), tree$tip.label)
  fits <- list(
    fast_k(tree, x, test = TRUE, nsim = 20, return_sim = TRUE,
           progress = FALSE),
    fast_k(tree, cbind(a = x, b = rev(x)), test = TRUE, nsim = 20,
           return_sim = TRUE, progress = FALSE)
  )
  .with_review_plot_device({
    for (fit in fits) {
      baseline <- plot_signal(fit)
      expect_silent(actual <- plot_signal(
        fit, xlim = c(0, 5), ylim = c(0, 4), lwd = 3,
        xaxs = "i", yaxs = "i"
      ))
      expect_equal(graphics::par("usr"), c(0, 5, 0, 4))
      expect_equal(actual, baseline)
    }
  })
})

test_that("D overlay and single-null plots accept graphics overrides", {
  set.seed(92)
  tree <- ape::rtree(24)
  x <- stats::setNames(rep(0:1, 12), tree$tip.label)
  fit <- fast_d(tree, x, nsim = 50, return_sim = TRUE, progress = FALSE)
  .with_review_plot_device({
    for (null in c("auto", "random", "brownian")) {
      baseline <- plot_signal(fit, null = null)
      expect_silent(actual <- plot_signal(
        fit, null = null, xlim = c(-4, 4), ylim = c(0, 4), lwd = 3,
        xaxs = "i", yaxs = "i"
      ))
      expect_equal(graphics::par("usr"), c(-4, 4, 0, 4))
      expect_equal(actual, baseline)
    }
  })
})

test_that("lambda profile and summaries accept graphics overrides", {
  set.seed(93)
  tree <- ape::rtree(12)
  x <- stats::setNames(stats::rnorm(12), tree$tip.label)
  fits <- list(
    fast_lambda(tree, x, test = TRUE, progress = FALSE),
    fast_lambda(tree, x, test = TRUE, lambda_profile = FALSE,
                progress = FALSE),
    fast_lambda(tree, cbind(a = x, b = rev(x)), test = TRUE,
                progress = FALSE)
  )
  .with_review_plot_device({
    for (fit in fits) {
      baseline <- plot_signal(fit)
      expect_silent(actual <- plot_signal(
        fit, xlim = c(-1, 2), ylim = c(-50, 5), lwd = 3,
        xaxs = "i", yaxs = "i"
      ))
      expect_equal(graphics::par("usr"), c(-1, 2, -50, 5))
      expect_equal(actual, baseline)
    }
  })
})


test_that("graphics defaults merge without hiding duplicate caller arguments", {
  capture <- function(x, lwd = 1) list(x = x, lwd = lwd)
  plot_call <- fastphylosig:::.signal_plot_call
  expect_identical(plot_call(capture, list(x = 1, lwd = 2)),
                   list(x = 1, lwd = 2))
  expect_identical(plot_call(capture, list(x = 1, lwd = 2), lwd = 3),
                   list(x = 1, lwd = 3))
  expect_error(plot_call(capture, list(x = 1, lwd = 2), lwd = 3, lwd = 4),
               "matched by multiple actual arguments")
})
