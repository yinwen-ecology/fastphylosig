# ACE release-blocker diagnosis.  This is an audit-only script.
# It does not modify package sources, tests, or statistical settings.

args <- commandArgs(trailingOnly = TRUE)
library_path <- Sys.getenv("FASTPHYLOSIG_LIBRARY", "")
if (nzchar(library_path)) {
  .libPaths(unique(c(library_path, .libPaths())))
}
suppressPackageStartupMessages(library(fastphylosig))
suppressPackageStartupMessages(library(ape))

make_tree <- function(n = 40L, heterogeneous = TRUE) {
  set.seed(20260913L + n + 2L)
  tr <- ape::rtree(n)
  tr$tip.label <- paste0("sp", seq_len(ape::Ntip(tr)))
  tr$edge.length <- if (heterogeneous) {
    10^seq(-2, 1, length.out = nrow(tr$edge))
  } else {
    pmax(as.numeric(tr$edge.length), 1e-4)
  }
  ape::reorder.phylo(tr, "cladewise")
}

capture <- function(fun) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(list(message = conditionMessage(e)),
                                  class = "ace_diag_error")
  )
  list(value = value, warnings = warnings)
}

tree <- make_tree()
x <- stats::setNames(
  factor(rep(c("A", "B", "C"), length.out = length(tree$tip.label)),
         levels = c("A", "B", "C")),
  tree$tip.label
)

cat("ACE release blocker diagnosis\n")
cat("R:", R.version.string, "\n")
cat("package:", as.character(packageVersion("fastphylosig")), "\n")
cat("tips:", length(x), "states:", paste(levels(x), collapse = ","),
    "branch_range:", paste(range(tree$edge.length), collapse = ","), "\n")

for (CI in c(FALSE, TRUE)) {
  result <- capture(function() fast_ace(
    x, tree, model = "ER", CI = CI, marginal = CI, progress = FALSE
  ))
  cat("\npublic fast_ace CI=", CI, "\n", sep = "")
  cat("warnings:", paste(result$warnings, collapse = " | "), "\n")
  if (inherits(result$value, "ace_diag_error")) {
    cat("error:", result$value$message, "\n")
  } else {
    cat("convergence:", result$value$convergence, "message:",
        result$value$message, "\n")
    cat("loglik:", format(result$value$loglik, digits = 17),
        "rates:", paste(format(result$value$rates, digits = 17), collapse = ","),
        "\n")
    cat("finite_rates:", all(is.finite(result$value$rates)),
        "finite_lik_anc:", if (is.null(result$value$lik.anc)) NA else
          all(is.finite(result$value$lik.anc)), "\n")
  }
}

phy <- fastphylosig:::.safe_canonicalize_core(tree)
workspace <- fastphylosig:::.fast_ace_workspace(
  fastphylosig:::.fast_ace_tree_workspace(phy), model = "ER", nl = 3L
)
tip_state <- as.integer(x)
dev <- function(p) fastphylosig:::fast_ace_discrete_deviance_cpp(
  edge = workspace$edge, edge_length = workspace$edge_length,
  tip_state = tip_state, rate_index = workspace$rate$rate_index, par = p
)

fit_one <- function(control = list()) {
  do.call(stats::nlminb, c(
    list(start = 0.1, objective = dev, lower = 0, upper = 1e50),
    if (length(control)) list(control = control) else list()
  ))
}

controls <- list(
  default = list(),
  more_iterations = list(eval.max = 1000L, iter.max = 1000L),
  relaxed_x = list(x.tol = 1e-7, xf.tol = 1e-7),
  strict_x = list(x.tol = 1e-10, xf.tol = 1e-10),
  strict_relative = list(rel.tol = 1e-12)
)

cat("\nmanual nlminb objective diagnostics\n")
cat("name,convergence,message,iterations,evaluations,par,objective,gradient,",
    "d_at_minus,d_at_plus\n", sep = "")
for (nm in names(controls)) {
  fit <- capture(function() fit_one(controls[[nm]]))
  if (inherits(fit$value, "ace_diag_error")) {
    cat(nm, ",ERROR,", fit$value$message, "\n", sep = "")
    next
  }
  p <- fit$value$par[[1L]]
  eps <- max(1e-8, abs(p) * 1e-6)
  dminus <- dev(max(0, p - eps))
  dplus <- dev(p + eps)
  gradient <- if (length(fit$value$gradient)) fit$value$gradient[[1L]] else NA_real_
  cat(nm, ",", fit$value$convergence, ",", fit$value$message, ",",
      fit$value$iterations, ",", paste(fit$value$evaluations, collapse = "/"), ",",
      format(p, digits = 17), ",", format(fit$value$objective, digits = 17), ",",
      format(gradient, digits = 17), ",", format(dminus, digits = 17), ",",
      format(dplus, digits = 17), "\n", sep = "")
}

ape_fit <- capture(function() ape::ace(
  x, tree, type = "discrete", method = "ML", model = "ER"
))
cat("\nape::ace\n")
cat("warnings:", paste(ape_fit$warnings, collapse = " | "), "\n")
if (!inherits(ape_fit$value, "ace_diag_error")) {
  cat("loglik:", format(ape_fit$value$loglik, digits = 17),
      "rates:", paste(format(ape_fit$value$rates, digits = 17), collapse = ","),
      "\n")
  cat("se:", paste(format(ape_fit$value$se, digits = 17), collapse = ","), "\n")
}

cat("\nseed=31 parity fixture\n")
set.seed(31)
tree31 <- ape::rtree(40L)
x31 <- factor(ape::rTraitDisc(tree31, model = "ER", k = 3,
                               rate = 1.5, states = letters[1:3]),
               levels = letters[1:3])
x31 <- x31[tree31$tip.label]
fit31 <- capture(function() fast_ace(x31, tree31, model = "ER", progress = FALSE))
ape31 <- capture(function() ape::ace(x31, tree31, type = "discrete",
                                     method = "ML", model = "ER"))
cat("fast warnings:", paste(fit31$warnings, collapse = " | "),
    "\nape warnings:", paste(ape31$warnings, collapse = " | "), "\n")
if (!inherits(fit31$value, "ace_diag_error") && !inherits(ape31$value, "ace_diag_error")) {
  cat("fast convergence:", fit31$value$convergence,
      "loglik_diff:", format(fit31$value$loglik - ape31$value$loglik,
                              digits = 17),
      "rate_diff:", format(fit31$value$rates - ape31$value$rates,
                            digits = 17), "\n")
}

run_variant <- function(label, tr, y) {
  fast <- capture(function() fast_ace(
    y, tr, model = "ER", CI = TRUE, marginal = TRUE, progress = FALSE
  ))
  ref <- capture(function() ape::ace(
    y, tr, type = "discrete", method = "ML", model = "ER"
  ))
  cat(label, ",fast_warning=", paste(fast$warnings, collapse = " | "),
      ",ape_warning=", paste(ref$warnings, collapse = " | "), sep = "")
  if (!inherits(fast$value, "ace_diag_error") &&
      !inherits(ref$value, "ace_diag_error")) {
    cat(",convergence=", fast$value$convergence,
        ",rate=", format(fast$value$rates, digits = 17),
        ",rate_diff=", format(fast$value$rates - ref$value$rates,
                                digits = 17),
        ",loglik_diff=", format(fast$value$loglik - ref$value$loglik,
                                 digits = 17),
        ",se_finite=", all(is.finite(fast$value$se)), sep = "")
  }
  cat("\n")
}

cat("\nfixture conditioning audit\n")
tree_homogeneous <- tree
tree_homogeneous$edge.length <- 1
tree_moderate <- tree
tree_moderate$edge.length <- 10^seq(-1, 1, length.out = nrow(tree$edge))
tree_original_order <- tree
tree_original_order$edge.length <- pmax(as.numeric(tree$edge.length), 1e-4)
y_cycle <- x
set.seed(3101)
y_random <- factor(ape::rTraitDisc(tree, model = "ER", k = 3,
                                    rate = 1.5, states = letters[1:3]),
                    levels = letters[1:3])
y_random <- y_random[tree$tip.label]
run_variant("cycle_heterogeneous", tree, y_cycle)
run_variant("cycle_moderate_heterogeneity", tree_moderate, y_cycle)
run_variant("cycle_homogeneous", tree_homogeneous, y_cycle)
run_variant("cycle_original_ape_lengths", tree_original_order, y_cycle)
run_variant("random_heterogeneous", tree, y_random)

cat("\nstarting-value audit for heterogeneous cycle fixture\n")
cat("start,convergence,message,par,objective,finite_rate\n")
for (start in c(1e-6, 1e-3, 0.01, 0.1, 1, 10, 100, 1000)) {
  fit <- capture(function() stats::nlminb(
    start = start, objective = dev, lower = 0, upper = 1e50
  ))
  if (inherits(fit$value, "ace_diag_error")) {
    cat(format(start, digits = 17), ",ERROR,", fit$value$message, "\n", sep = "")
  } else {
    cat(format(start, digits = 17), ",", fit$value$convergence, ",",
        fit$value$message, ",", format(fit$value$par, digits = 17), ",",
        format(fit$value$objective, digits = 17), ",",
        is.finite(fit$value$par), "\n", sep = "")
  }
}

cat("\nobjective local derivative/curvature at heterogeneous fit\n")
fit <- fit_one()
p <- fit$par[[1L]]
for (scale in c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)) {
  h <- max(scale, abs(p) * scale)
  left <- dev(p - h)
  center <- dev(p)
  right <- dev(p + h)
  derivative <- (right - left) / (2 * h)
  curvature <- (right - 2 * center + left) / (h * h)
  cat("relative_step=", format(scale, digits = 4),
      ",h=", format(h, digits = 8),
      ",derivative=", format(derivative, digits = 17),
      ",curvature=", format(curvature, digits = 17),
      ",delta_objective=", format(max(abs(left - center),
                                      abs(right - center)), digits = 17), "\n",
      sep = "")
}

ape_like_dev <- local({
  tr <- ape::reorder.phylo(tree, "postorder")
  n_tip <- length(tr$tip.label)
  n_node <- tr$Nnode
  e1 <- tr$edge[, 1L]
  e2 <- tr$edge[, 2L]
  el <- tr$edge.length
  function(p) {
    Q <- matrix(p, 3L, 3L)
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    decomp <- eigen(Q)
    P <- lapply(el, function(t) Re(
      decomp$vectors %*% diag(exp(decomp$values * t), 3L, 3L) %*%
        solve(decomp$vectors)
    ))
    liks <- matrix(0, n_tip + n_node, 3L)
    liks[cbind(seq_len(n_tip), tip_state)] <- 1
    comp <- numeric(n_tip + n_node)
    for (e in seq.int(1L, length.out = n_node, by = 2L)) {
      left <- P[[e]] %*% liks[e2[e], ]
      right <- P[[e + 1L]] %*% liks[e2[e + 1L], ]
      v <- as.numeric(left) * as.numeric(right)
      comp[e1[e]] <- sum(v)
      liks[e1[e], ] <- v / comp[e1[e]]
    }
    -2 * sum(log(comp[-seq_len(n_tip)]))
  }
})

cat("\nC++ versus independent R eigen objective\n")
cat("rate,cpp_dev,r_eigen_dev,absolute_difference\n")
for (rate in c(0.1, 1, 10, 100, 200, 218.6295090836829,
               218.62976768174769, 230, 256, 1000)) {
  cpp_value <- dev(rate)
  r_value <- ape_like_dev(rate)
  cat(format(rate, digits = 17), ",", format(cpp_value, digits = 17), ",",
      format(r_value, digits = 17), ",",
      format(cpp_value - r_value, digits = 17), "\n", sep = "")
}

cat("\nCI=false direct SE bypass check\n")
fit_no_ci <- capture(function() fast_ace(x, tree, model = "ER", CI = FALSE,
                                          progress = FALSE))
cat("warnings:", paste(fit_no_ci$warnings, collapse = " | "),
    "se_present:", !inherits(fit_no_ci$value, "ace_diag_error") &&
      !is.null(fit_no_ci$value$se), "se_nan:",
    if (inherits(fit_no_ci$value, "ace_diag_error")) NA else
      all(is.nan(fit_no_ci$value$se)), "\n")
