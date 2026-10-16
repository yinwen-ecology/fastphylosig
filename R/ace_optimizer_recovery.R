# Deterministic recovery for an unfinished Delta ARD fit. All objectives share
# the caller's existing C++ pointer; no new tree, random start or RNG scope.
.fast_ace_recover_ard <- function(dev, first, recheck = FALSE) {
  valid <- function(z) !inherits(z,"error") && is.finite(z$objective) &&
    z$objective < 1e49 && all(is.finite(z$par)) && all(z$par >= 0)
  if(first$convergence == 0L || !valid(first)) return(first)
  original <- first
  if(isTRUE(recheck)) {
    first$objective <- dev(first$par)
    if(!valid(first))stop("ARD recovery has no numerically valid initial deviance.",call.=FALSE)
  }
  attempts <- list(list(label="baseline",fit=first,valid=TRUE))
  run <- function(start, label, solver) {
    z <- tryCatch(solver(start),error=identity)
    attempts[[length(attempts)+1L]] <<- list(label=label,fit=z,valid=valid(z))
  }
  run(first$par,"scaled_continuation",function(p) stats::nlminb(
    p,dev,lower=rep(0,length(p)),upper=rep(1e50,length(p)),
    scale=1/pmax(abs(p),.1),control=list(iter.max=150L,eval.max=300L)))
  choice <- .fast_ace_select_candidate(attempts)
  if(attempts[[choice$selected]]$fit$convergence != 0L) {
    run(attempts[[choice$selected]]$fit$par,"bounded_lbfgsb",function(p) {
      z <- stats::optim(p,dev,method="L-BFGS-B",lower=rep(0,length(p)),
        upper=rep(1e50,length(p)),control=list(maxit=150L,ndeps=rep(1e-6,length(p))))
      list(par=z$par,objective=z$value,convergence=z$convergence,message=z$message,
        iterations=NA_integer_,evaluations=z$counts)
    })
    choice <- .fast_ace_select_candidate(attempts)
  }
  fit <- attempts[[choice$selected]]$fit
  table <- do.call(rbind,lapply(attempts,function(z) {
    bad <- inherits(z$fit,"error")
    data.frame(start=z$label,valid=z$valid,
      objective=if(bad)NA_real_ else z$fit$objective,
      convergence=if(bad)NA_integer_ else z$fit$convergence,
      iterations=if(bad)NA_integer_ else z$fit$iterations,
      evaluations=if(bad)NA_integer_ else unname(z$fit$evaluations[["function"]]),
      message=if(bad)conditionMessage(z$fit) else z$fit$message,
      stringsAsFactors=FALSE)
  }))
  fit$diagnostics <- list(method="bounded scaled ARD recovery",
    selected_attempt=choice$selected,attempts=table,baseline_deviance=first$objective,
    best_deviance=fit$objective,deviance_improvement=first$objective-fit$objective,
    optimizer_retry=list(original=original,recovered=fit$convergence==0L,
      scaled_iterations=150L,scaled_evaluations=300L,lbfgsb_iterations=150L))
  if(!identical(original$objective,first$objective))
    fit$diagnostics$initial_numerical_correction <- original$objective-first$objective
  if(choice$selected != choice$minimum)fit$diagnostics$selection_recovery <- list(
    minimum_attempt=choice$minimum,minimum_deviance=table$objective[choice$minimum],
    tolerance=choice$tolerance,reason="numerically equivalent converged candidate")
  fit
}

# Shared by the optional common multistart and bounded local recovery. Prefer
# actual convergence only when its objective is numerically equivalent.
.fast_ace_select_candidate <- function(attempts) {
  valid <- which(vapply(attempts,`[[`,logical(1),"valid"))
  objectives <- vapply(attempts[valid],function(z)z$fit$objective,numeric(1))
  selected <- minimum <- valid[which.min(objectives)]
  tolerance <- NULL
  if(attempts[[selected]]$fit$convergence != 0L) {
    tolerance <- 1e-7*(1+abs(attempts[[selected]]$fit$objective))
    converged <- valid[vapply(attempts[valid],function(z)z$fit$convergence==0L,logical(1))]
    if(length(converged)) {
      best <- converged[which.min(vapply(attempts[converged],function(z)z$fit$objective,numeric(1)))]
      if(attempts[[best]]$fit$objective <= attempts[[selected]]$fit$objective+tolerance)
        selected <- best
    }
  }
  list(selected=selected,minimum=minimum,tolerance=tolerance)
}
