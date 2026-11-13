# Recheck the fit with dimensionless optimized candidate starts.
# ip retains its documented units and remains the original optimizer start.
# Candidate optimization is deterministic and never consumes RNG.
.fast_ace_recheck_rate_scale <- function(dev, first, start, workspace) {
  branch_scale <- stats::median(workspace$edge_length)
  if (!is.finite(branch_scale) || branch_scale <= 0 ||
      !is.finite(1 / branch_scale)) return(first)

  valid <- function(z) !inherits(z, "error") && is.finite(z$objective) &&
    z$objective < 1e49 && all(is.finite(z$par)) && all(z$par >= 0)
  baseline_valid <- valid(first)
  baseline_objective <- if (baseline_valid) first$objective else Inf
  tolerance <- if (baseline_valid) 1e-7 * (1 + abs(first$objective)) else 0
  intensities <- c(.01, .1, 1)
  rates <- intensities / branch_scale
  usable <- is.finite(rates) & rates > 0 & rates <= 1e50
  intensities <- intensities[usable]
  rates <- rates[usable]
  if (!length(rates)) return(first)
  # Optimize each dimensionless start before judging it. An equal-rate
  # point can have poor likelihood while leading to a good asymmetric fit.
  upper <- min(.Machine$double.xmax, 1e50 * branch_scale)
  candidates <- lapply(intensities, function(intensity) {
    z <- tryCatch(stats::nlminb(
      start = rep(intensity, length(start)),
      objective = function(p) dev(p / branch_scale),
      lower = rep(0, length(start)), upper = rep(upper, length(start))),
      error = identity)
    if (!inherits(z, "error")) z$par <- z$par / branch_scale
    z
  })
  labels <- paste0("branch_scale_", intensities)
  # Replay unfinished candidates from the same start with the reference-search
  # budget. Restarting from the endpoint would reset the quasi-Newton history;
  # replay preserves the original path beyond nlminb's default iteration limit.
  unfinished <- which(vapply(candidates, function(z)
    valid(z) && z$convergence != 0L, logical(1)))
  for (i in unfinished) {
    z <- tryCatch(stats::nlminb(
      start = rep(intensities[i], length(start)),
      objective = function(p) dev(p / branch_scale),
      lower = rep(0, length(start)), upper = rep(upper, length(start)),
      control = list(iter.max = 1500L, eval.max = 4000L)), error = identity)
    if (!inherits(z, "error")) z$par <- z$par / branch_scale
    candidates[[length(candidates) + 1L]] <- z
    intensities <- c(intensities, intensities[i])
    labels <- c(labels, paste0(labels[i], "_extended"))
  }
  if (length(start) == 1L) {
    # ER is scalar. Partition at the existing scale starts and fitted rates:
    # refining each interval avoids letting one broad Brent search skip a
    # low-rate optimum in favour of the saturated high-rate plateau.
    fitted <- vapply(candidates, function(z)
      if (valid(z)) z$par * branch_scale else NA_real_, numeric(1))
    knots <- sort(unique(c(0, intensities,
      if (baseline_valid) first$par * branch_scale, fitted)))
    knots <- knots[is.finite(knots) & knots >= 0 & knots <= upper]
    for (i in seq_len(length(knots) - 1L)) {
      lo <- knots[i]; hi <- knots[i + 1L]
      if (hi <= lo) next
      calls <- 0L
      objective <- function(p) {
        calls <<- calls + 1L
        dev(p / branch_scale)
      }
      z <- tryCatch({
        scalar <- stats::optimize(objective, c(lo, hi), tol = 1e-10)
        values <- c(objective(lo), scalar$objective, objective(hi))
        pars <- c(lo, scalar$minimum, hi)
        pick <- which.min(values)
        list(par = pars[pick] / branch_scale, objective = values[pick],
          convergence = 0L, iterations = NA_integer_,
          evaluations = stats::setNames(calls, "function"), message = "bounded scalar refinement")
      }, error = identity)
      candidates[[length(candidates) + 1L]] <- z
      intensities <- c(intensities, NA_real_)
      labels <- c(labels, paste0("scalar_interval_", i))
    }
  }
  probes <- vapply(candidates, function(z)
    if (valid(z)) z$objective else Inf, numeric(1))
  best <- which.min(probes)
  if (!length(best) || !is.finite(probes[best]) ||
      probes[best] >= baseline_objective - tolerance) return(first)
  candidate <- candidates[[best]]
  attempt_row <- function(fit, label) data.frame(
    start = label, valid = valid(fit), objective = fit$objective,
    convergence = fit$convergence,
    iterations = if (is.null(fit$iterations)) NA_integer_ else fit$iterations,
    evaluations = if (is.null(fit$evaluations)) NA_integer_ else
      unname(fit$evaluations[["function"]]),
    message = if (is.null(fit$message)) "" else fit$message,
    stringsAsFactors = FALSE)
  attempts <- first$diagnostics$attempts
  if (!is.data.frame(attempts)) attempts <- attempt_row(first, "baseline")
  baseline_rows <- nrow(attempts)
  for (i in seq_along(candidates)) {
    z <- candidates[[i]]
    if (inherits(z, "error")) z <- list(objective = Inf, par = NA_real_,
      convergence = NA_integer_, message = conditionMessage(z))
    attempts <- rbind(attempts, attempt_row(z, labels[i]))
  }
  candidate$diagnostics <- list(
    method = "branch-scale optimizer recovery",
    selected_attempt = baseline_rows + best, attempts = attempts,
    baseline_deviance = first$objective,
    best_deviance = candidate$objective,
    deviance_improvement = first$objective - candidate$objective,
    rate_scale_recovery = list(branch_scale = branch_scale,
      candidate_intensities = intensities, candidate_deviances = probes,
      selected_start = intensities[best]),
    optimizer_retry = list(original = first,
      recovered = candidate$convergence == 0L))
  candidate
}
