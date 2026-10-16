# Opt-in finite-sample nuisance bounds. Expensive Delta fits are performed once
# per proposal draw; continuous-Q search uses native likelihood intervals only.
.heterogeneity_nuisance_control <- function(z) {
  allowed <- c("gamma", "tolerance", "max_boxes", "rate_bounds",
    "proposal_intensities", "max_bank_bytes", "boundary_power", "proposal_family", "correlation", "proposal_counts")
  if (!is.list(z) || (length(z) && (is.null(names(z)) || anyNA(names(z)) ||
      any(!nzchar(names(z))) || anyDuplicated(names(z)))) ||
      length(setdiff(names(z), allowed)))
    stop("nuisance_control must be a uniquely named list of supported controls.", call. = FALSE)
  out <- utils::modifyList(list(gamma = .001, tolerance = .01,
    max_boxes = 128L, rate_bounds = NULL,
    proposal_intensities = c(.01, .1, 1, 10, Inf), max_bank_bytes = 64 * 1024^2,
    boundary_power = 24L, proposal_family = "auto", correlation = "auto"), z)
  out$proposal_family <- match.arg(out$proposal_family,c("auto","ER","ARD"))
  if("proposal_counts"%in%names(out)&&(!is.logical(out$proposal_counts)||length(out$proposal_counts)!=1L||is.na(out$proposal_counts)))
    stop("proposal_counts must be TRUE or FALSE.",call.=FALSE)
  out$correlation <- match.arg(out$correlation,c("auto","none","taylor","forest","transient"))
  out$boundary_power <- .heterogeneity_integer(out$boundary_power,"boundary_power",8L)
  if(out$boundary_power>30L)stop("boundary_power must not exceed 30.",call.=FALSE)
  for (nm in c("gamma", "tolerance"))
    if (!is.numeric(out[[nm]]) || length(out[[nm]]) != 1L ||
        !is.finite(out[[nm]]) || out[[nm]] < 0 || out[[nm]] >= 1)
      stop("nuisance_control$", nm, " must be in [0,1).", call. = FALSE)
  out$max_boxes <- .heterogeneity_integer(out$max_boxes, "max_boxes", 1L)
  if (!is.numeric(out$max_bank_bytes) || length(out$max_bank_bytes) != 1L ||
      !is.finite(out$max_bank_bytes) || out$max_bank_bytes <= 0)
    stop("max_bank_bytes must be positive and finite.", call. = FALSE)
  q <- out$proposal_intensities
  if (!is.numeric(q) || !length(q) || length(q) > 16L || anyNA(q) ||
      any(q <= 0) || !any(is.infinite(q)) || anyDuplicated(q))
    stop("proposal_intensities must be distinct positive numbers and include Inf (at most 16).", call. = FALSE)
  out
}

.heterogeneity_delta_state_space <- function(X, control) {
  spaces <- control$state_space
  if (is.list(spaces)) {
    if (is.null(names(spaces)) || anyNA(names(spaces)) || anyDuplicated(names(spaces)) ||
        !setequal(names(spaces), colnames(X)))
      stop("A list state_space must name every trait exactly once.", call. = FALSE)
    spaces <- spaces[colnames(X)]
  } else spaces <- rep(list(spaces), ncol(X))
  original <- attr(X, "state_labels")
  for (j in seq_len(ncol(X))) {
    s <- spaces[[j]]
    if ((!is.character(s) && !is.numeric(s)) || length(s) < 2L ||
        length(s) > 8L || anyNA(s) || any(!nzchar(as.character(s))) ||
        anyDuplicated(as.character(s)) || (is.numeric(s) && any(!is.finite(s))))
      stop("state_space must explicitly declare 2--8 distinct finite state labels.", call. = FALSE)
    s <- as.character(s)
    if (control$model == "ARD" && length(s) != 3L)
      stop("Nuisance-controlled ARD currently requires exactly three declared states.", call. = FALSE)
    if (length(setdiff(original[[j]], s)))
      stop("Observed states are missing from state_space for trait '", colnames(X)[j], "'.", call. = FALSE)
    used <- !is.na(X[,j])
    X[used,j] <- match(original[[j]][X[used,j]], s)
    spaces[[j]] <- s
  }
  attr(X, "state_labels") <- spaces
  X
}

.heterogeneity_nuisance_domain <- function(control, states) {
  index <- .fast_ace_rate_index(control$model, length(states))$rate_index
  d <- max(index)
  bounds <- control$nuisance_control$rate_bounds
  if (is.null(bounds)) bounds <- cbind(lower = rep(0, d), upper = rep(Inf, d))
  if (is.numeric(bounds) && is.null(dim(bounds)) && length(bounds) == 2L && d == 1L)
    bounds <- matrix(bounds, 1L, 2L)
  if (!is.numeric(bounds) || !is.matrix(bounds) || !identical(dim(bounds), c(d, 2L)) ||
      anyNA(bounds) || any(!is.finite(bounds[,1])) || any(bounds[,1] < 0) ||
      any(bounds[,2] < bounds[,1]))
    stop("rate_bounds must have one lower/upper row per native rate, with finite nonnegative lower bounds.", call. = FALSE)
  labels <- character(d)
  if (control$model == "ER") labels[] <- "q" else for (i in seq_along(states))
    for (j in seq_along(states)) if (i != j) labels[index[i,j]] <- paste0(states[i], " -> ", states[j])
  dimnames(bounds) <- list(labels, c("lower", "upper"))
  list(index = index, bounds = bounds, states = states, root = "uniform",
    model = control$model)
}

.heterogeneity_nuisance_search <- function(pointer, domain, control, scale) {
  b <- domain$bounds
  bank <- attr(pointer,"bank_diagnostics",exact=TRUE)
  compact <- function(x) { z <- 1/(1+scale/x); z[x==0] <- 0; z }
  correlation_evaluations <- correlation_fallbacks <- 0L
  forest_evaluations <- forest_fallbacks <- forest_subcells <- 0L
  transient_evaluations <- transient_fallbacks <- transient_subcells <- witness_evaluations <- 0L
  transient_messages <- transient_reuses <- mass_evaluations <- mass_fallbacks <- 0
  evaluate <- function(lo, hi) {
    fit <- heterogeneity_nuisance_box_cpp(pointer,as.numeric(lo),as.numeric(hi),control$gamma,
      boundary_power=control$boundary_power,correlation=if(domain$model=="ARD")
        match(control$correlation,c("none","auto","taylor","forest","transient"))-1L else 0L)
    if(isTRUE(fit$correlation_used))correlation_evaluations <<- correlation_evaluations+1L
    if(isTRUE(fit$correlation_fallback))correlation_fallbacks <<- correlation_fallbacks+1L
    if(isTRUE(fit$forest_used)) {
      forest_evaluations <<- forest_evaluations+1L
      forest_subcells <<- forest_subcells+fit$forest_cells
    }
    if(isTRUE(fit$forest_fallback))forest_fallbacks <<- forest_fallbacks+1L
    if(isTRUE(fit$transient_used)) {
      transient_evaluations <<- transient_evaluations+1L
      transient_subcells <<- transient_subcells+fit$transient_cells
      transient_messages <<- transient_messages+fit$transient_messages
      transient_reuses <<- transient_reuses+fit$transient_reuses
    }
    if(isTRUE(fit$transient_fallback))transient_fallbacks <<- transient_fallbacks+1L
    if(isTRUE(fit$mass_used))mass_evaluations <<- mass_evaluations+1L
    if(isTRUE(fit$mass_fallback))mass_fallbacks <<- mass_fallbacks+1L
    fit
  }
  midpoint <- function(lo, hi) {
    mid <- lo/2 + hi/2
    infinite <- is.infinite(hi)
    mid[infinite] <- lo[infinite] + pmax(scale, lo[infinite])
    mid[!is.finite(mid)] <- lo[!is.finite(mid)]
    if(domain$model=="ARD") {
      z <- (compact(lo)+compact(hi))/2
      candidate <- scale*z/(1-z)
      # Rounded compact coordinates can coincide near a very large endpoint.
      use <- is.finite(candidate) & candidate>=lo & candidate<=hi
      mid[use] <- candidate[use]
    }
    mid
  }
  lower <- 0; best_rates <- rep(NA_real_, nrow(b)); evaluations <- 0L
  best_observed_weight <- best_ESS <- NA_real_
  point <- function(lo, hi, p=midpoint(lo, hi)) {
    fit <- evaluate(p, p); evaluations <<- evaluations + 1L
    if (fit$included && (fit$lower > lower || anyNA(best_rates))) {
      lower <<- fit$lower; best_rates <<- p
      best_observed_weight <<- fit$observed_weight_upper; best_ESS <<- fit$ESS_lower
    }
  }
  initial <- evaluate(b[,1], b[,2]); evaluations <- evaluations + 1L
  cells <- list(list(lo = b[,1], hi = b[,2], fit = initial))
  if (!initial$excluded) point(b[,1], b[,2])
  if(!initial$excluded&&domain$model=="ARD"&&control$correlation%in%c("auto","transient")) {
    # Fixed witnesses improve only the certified lower bound of the supremum.
    # They are checked against C_gamma; no proposal or random draw is changed.
    seen <- list(midpoint(b[,1],b[,2]))
    for(intensity in c(1,10,100)) {
      p <- pmax(b[,1],pmin(b[,2],rep(intensity*scale,nrow(b))))
      if(any(!is.finite(p))||any(vapply(seen,function(x)identical(as.numeric(x),as.numeric(p)),logical(1))))next
      point(b[,1],b[,2],p);seen[[length(seen)+1L]]<-p
      witness_evaluations <- witness_evaluations+1L
    }
    # A closed pair plus one isolated state can assign zero likelihood to
    # every non-tail phenotype. Native support checks only select candidates;
    # the existing point engine must certify both the rank and membership.
    for(isolated in if(is.matrix(domain$index))bank$support_isolated_states else integer()) {
      pair <- setdiff(1:3,isolated)
      edges <- domain$index[cbind(pair,rev(pair))]
      if(any(b[-edges,1]>0))next
      for(intensity in c(.1,1,10)) {
        p <- numeric(6);p[edges] <- intensity*scale
        p <- pmax(b[,1],pmin(b[,2],p))
        point(b[,1],b[,2],p);witness_evaluations <- witness_evaluations+1L
        if(heterogeneity_nuisance_add_gamma_cpp(lower,1,control$gamma)[1]>=1)break
      }
      if(heterogeneity_nuisance_add_gamma_cpp(lower,1,control$gamma)[1]>=1)break
    }
  }
  unbounded_ard <- domain$model == "ARD" && any(is.infinite(b[,2]))
  stop_reason <- "box_budget"
  while (length(cells) < control$max_boxes) {
    if(heterogeneity_nuisance_add_gamma_cpp(lower,1,control$gamma)[1]>=1) {
      stop_reason <- "certified_P_one";break
    }
    upper <- vapply(cells, function(x) x$fit$upper, numeric(1))
    if (max(upper) - lower <= control$tolerance) {stop_reason <- "gap_target";break}
    # Among equal bounds, refine the largest compact-coordinate cell first.
    candidates <- which(upper==max(upper))
    sizes <- vapply(cells[candidates],function(x)
      prod(compact(x$hi)-compact(x$lo)),numeric(1))
    i <- candidates[which.max(sizes)]; cell <- cells[[i]]
    widths <- if(domain$model=="ARD") compact(cell$hi)-compact(cell$lo) else
      log1p(cell$hi/scale) - log1p(cell$lo/scale)
    dimension <- which.max(widths)
    if (cell$hi[dimension] <= cell$lo[dimension]) {stop_reason <- "unsplittable";break}
    mid <- midpoint(cell$lo, cell$hi)[dimension]
    if(domain$model=="ARD") {
      z <- (compact(cell$lo[dimension])+compact(cell$hi[dimension]))/2
      mid <- scale*z/(1-z)
    }
    if (domain$model == "ER") {
      # Correlation-scale midpoint resolves both zero and infinite endpoints.
      a <- exp(-cell$lo/scale); c <- exp(-cell$hi/scale)
      if (a > 0) mid <- -scale * log((a + c)/2)
    }
    if (!is.finite(mid) || mid <= cell$lo[dimension] || mid >= cell$hi[dimension]) {
      stop_reason <- "unsplittable";break
    }
    left_hi <- cell$hi; left_hi[dimension] <- mid
    right_lo <- cell$lo; right_lo[dimension] <- mid
    first <- evaluate(cell$lo, left_hi); second <- evaluate(right_lo, cell$hi)
    evaluations <- evaluations + 2L
    cells[[i]] <- list(lo = cell$lo, hi = left_hi, fit = first)
    cells[[length(cells)+1L]] <- list(lo = right_lo, hi = cell$hi, fit = second)
    if (!first$excluded) point(cell$lo, left_hi)
    if (!second$excluded) point(right_lo, cell$hi)
  }
  upper <- max(vapply(cells, function(x) x$fit$upper, numeric(1)))
  # The native helper encloses the final addition too; no R rounding can make
  # the reported bound smaller than gamma + the covered-domain supremum.
  bounds <- heterogeneity_nuisance_add_gamma_cpp(lower, upper, control$gamma)
  if(upper-lower<=control$tolerance)stop_reason <- "gap_target"
  if(bounds[1]>=1)stop_reason <- "certified_P_one"
  one_reason <- if(bounds[2]<1) "not_one" else if(!is.null(bank)&&bank$non_tail_count==0L)
    "all_upper_tail" else if(bounds[1]>=1) "certified_saturation" else if(upper-lower>control$tolerance)
    "unresolved_envelope" else if(upper<1)
    "gamma_ceiling" else "unresolved_envelope"
  list(lower = bounds[1], upper = bounds[2], gap = upper - lower,
    converged = upper - lower <= control$tolerance,
    cells = length(cells), evaluations = evaluations, best_rates = best_rates,
    best_point_observed_weight_upper = best_observed_weight, best_point_ESS_lower = best_ESS,
    n_unique_patterns = initial$n_unique_patterns,
    correlation_evaluations=correlation_evaluations,correlation_fallbacks=correlation_fallbacks,
    forest_evaluations=forest_evaluations,forest_fallbacks=forest_fallbacks,forest_subcells=forest_subcells,
    transient_evaluations=transient_evaluations,transient_fallbacks=transient_fallbacks,
    transient_subcells=transient_subcells,witness_evaluations=witness_evaluations,
    transient_messages=transient_messages,transient_reuses=transient_reuses,
    mass_evaluations=mass_evaluations,mass_fallbacks=mass_fallbacks,
    raw_lower=lower,raw_upper=upper,final_width=bounds[2]-bounds[1],
    gamma_clipped=upper<1&&bounds[2]>=1,
    stop_reason=stop_reason,P_one_reason=one_reason,bank_diagnostics=bank,
    unbounded_ard = unbounded_ard, boundary_covered=TRUE,
    boundary_Euler_error=if(unbounded_ard)32/2^control$boundary_power else 0,
    remaining = lapply(cells, function(x)
      list(lower = x$lo, upper = x$hi, P_upper = x$fit$upper, excluded = x$fit$excluded)))
}

.heterogeneity_delta_nuisance <- function(tree, core, global_ctx, contexts,
    observed_x, observed, states, control, nsim, keep_null) {
  settings <- control$nuisance_control
  domain <- .heterogeneity_nuisance_domain(control, states)
  n <- length(observed_x)
  if (4 * n * (nsim + 1) > settings$max_bank_bytes)
    stop("Nuisance reference bank exceeds max_bank_bytes; specify a smaller preplanned budget or a larger memory limit.", call. = FALSE)
  scale <- 1/(length(states) * stats::median(tree$edge.length))
  proposal_rates <- settings$proposal_intensities * scale
  er_index <- matrix(1L, length(states), length(states)); diag(er_index) <- 0L
  proposals <- list(); pointer <- NULL
  on.exit({
    if (!is.null(pointer)) heterogeneity_nuisance_release_cpp(pointer)
    for (p in proposals) heterogeneity_mk_release_cpp(p)
  }, add = TRUE)
  for (q in proposal_rates) proposals[[length(proposals)+1L]] <-
    heterogeneity_mk_prepare_cpp(core, er_index, q)
  extra <- .heterogeneity_ard_proposals(domain,settings,scale)
  for(r in extra)proposals[[length(proposals)+1L]] <- heterogeneity_mk_prepare_cpp(core,domain$index,r)
  count_proposal <- isTRUE(settings$proposal_counts)
  between <- identical(control$statistic,"between_clades")
  n_components <- length(proposals)
  patterns <- matrix(NA_integer_, n, nsim+1L); patterns[,1L] <- observed_x
  draws <- rep(NA_real_, nsim); tail <- rep(TRUE, nsim+1L)
  fields <- c("nsim_MCMC_warning", "nsim_optimizer_warning", "nsim_numerical_recovery",
    "nsim_entropy_recovery", "nsim_other_warning", "nsim_calculation_failure",
    "nsim_unavailable", "nsim_single_state")
  counts <- stats::setNames(integer(length(fields)), fields)
  retries <- c(nsim_optimizer_retry = 0L, nsim_optimizer_recovered = 0L)
  diagnostics <- matrix(0L, nsim, length(fields)+length(retries),
    dimnames = list(NULL,c(fields,names(retries))))
  failures <- vector("list", nsim)
  started <- proc.time()[[3L]]
  for (b in seq_len(nsim)) {
    component <- sample.int(if(count_proposal)2L*n_components else n_components, 1L)
    x <- if(component>n_components)heterogeneity_count_simulate_cpp(length(states),!is.na(observed_x)) else
      heterogeneity_mk_simulate_cpp(proposals[[component]], 1L)[,1L]
    x[is.na(observed_x)] <- NA_integer_; patterns[,b+1L] <- x
    full <- .heterogeneity_delta_evaluate(global_ctx, x)
    local <- lapply(contexts, .heterogeneity_delta_evaluate, x = x)
    flags <- heterogeneity_delta_flags_cpp(full, local)
    retry <- attr(flags, "optimizer_retry", exact = TRUE)
    if (!is.null(retry)) {
      retries <- retries + retry; diagnostics[b,names(retry)] <- retry
      attr(flags, "optimizer_retry") <- NULL
    }
    diagnostics[b,fields] <- flags
    if (flags[["nsim_unavailable"]] > 0L) failures[[b]] <- list(full = full, local = local)
    counts <- counts + flags
    values <- vapply(local, `[[`, numeric(1), "estimate")
    if (is.finite(full$estimate) && all(is.finite(values))) {
      draws[b] <- if(between)heterogeneity_delta_between_cpp(values) else
        sqrt(mean((log(values) - log(full$estimate))^2))
      tail[b+1L] <- draws[b] >= observed
    }
  }
  bank_seconds <- proc.time()[[3L]] - started
  summary <- as.list(counts)
  if (retries[1L] > 0L) summary <- c(summary, as.list(retries))
  summary$nsim_fit_warning <- counts[["nsim_optimizer_warning"]]
  summary$nsim_successful <- sum(is.finite(draws))
  summary$P_heterogeneity_kind <- "nuisance_upper_bound"
  summary$nuisance_reference_attempted <- nsim
  summary$reference_unavailable <- "single_state_upper_tail"
  summary$exceedance_count <- NA_integer_
  summary$finite_exceedance_count <- sum(draws >= observed, na.rm = TRUE)
  summary$MCSE_P_heterogeneity <- NA_real_
  model <- c(domain, list(proposal_rates = proposal_rates, proposal_ARD_rates=extra,
    proposal = if(length(extra))"fixed equal-weight ER/ARD mixture" else "fixed equal-weight ER mixture",
    correlation=settings$correlation,boundary_power=settings$boundary_power,
    gamma = settings$gamma, rate_bounds_supplied = !is.null(settings$rate_bounds),
    requires_predeclared_domain = TRUE, reference_draws_attempted = nsim,
    reference_diagnostics = diagnostics, reference_unavailable_details = failures))
  if(count_proposal) {
    model$proposal <- paste0("fixed half Mk mixture, half count-universal sequential proposal")
    model$proposal_counts <- list(alpha=1,weight=.5,density="actual sequential CDF kernel",observed_counts_used=FALSE)
  }
  if(between)model$statistic <- "SD(log(Delta_local))"
  if (counts[["nsim_calculation_failure"]] > 0L ||
      counts[["nsim_unavailable"]] != counts[["nsim_single_state"]]) {
    summary$status <- "simulation_failure"
    summary$message <- "Actual calculation failures retained; nuisance probability unavailable, no draws replaced"
    return(list(summary = summary, model = model, null = if (keep_null) draws else NULL))
  }
  started <- proc.time()[[3L]]
  components <- lapply(proposals, function(p)heterogeneity_mk_loglik_bounds_cpp(p, core, patterns))
  if(count_proposal)components <- c(components,rep(list(heterogeneity_count_loglik_cpp(patterns,length(states))),n_components))
  g <- heterogeneity_nuisance_mixture_cpp(components)
  density_seconds <- proc.time()[[3L]] - started
  pointer <- heterogeneity_nuisance_prepare_cpp(core, patterns, domain$index, g,
    tail, control$model == "ER")
  started <- proc.time()[[3L]]
  search <- .heterogeneity_nuisance_search(pointer, domain, settings, scale)
  model$timing_seconds <- c(reference_Delta = bank_seconds,
    proposal_density = density_seconds, nuisance_search = proc.time()[[3L]] - started)
  summary$P_heterogeneity <- search$upper
  summary$P_nuisance_lower <- search$lower; summary$P_nuisance_upper <- search$upper
  summary$nuisance_gap <- search$gap; summary$nuisance_converged <- search$converged
  summary$nuisance_cells <- search$cells; summary$nuisance_evaluations <- search$evaluations
  summary$nuisance_P_one_reason <- search$P_one_reason
  summary$nuisance_stop_reason <- search$stop_reason
  summary$nuisance_final_width <- search$final_width
  summary$status <- if (search$converged) "nuisance_bound" else "nuisance_budget_limited"
  summary$message <- if (!search$converged)
    "Search stopped before the gap target: complete-domain upper bound retained, candidate maximum is not used as P" else
    "Conservative nuisance-domain probability bound; finite-sample validity assumes the declared model, state space and preplanned proposal/domain"
  model$search <- search
  model$proposal_log_density_bounds <- g
  if (keep_null) { model$reference_patterns <- patterns; model$upper_tail <- tail }
  list(summary = summary, model = model, null = if (keep_null) draws else NULL)
}

.heterogeneity_ard_proposals <- function(domain,settings,scale) {
  enabled <- settings$proposal_family=="ARD" ||
    (settings$proposal_family=="auto" && domain$model=="ARD")
  if(!enabled)return(list())
  if(domain$model!="ARD" || length(domain$states)!=3L)
    stop("proposal_family='ARD' requires the three-state ARD model.",call.=FALSE)
  out <- list();index <- domain$index;off <- index>0
  pack <- function(Q) { r<-numeric(6);r[index[off]]<-Q[off];r*scale }
  for(rare in 1:3) {
    pi <- rep(.475,3);pi[rare]<-.05
    Q <- matrix(rep(3*pi,each=3),3,3);diag(Q)<-0
    out[[length(out)+1L]]<-pack(Q)
  }
  for(pair in list(c(1L,2L),c(1L,3L),c(2L,3L))) {
    Q <- matrix(.1,3,3);diag(Q)<-0
    Q[pair[1],pair[2]]<-Q[pair[2],pair[1]]<-10
    out[[length(out)+1L]]<-pack(Q)
  }
  out
}
