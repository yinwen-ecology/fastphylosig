# Runtime/progress helpers ----------------------------------------------------

# Public engines use one small runtime context so progress output and elapsed
# time are identical across scalar and batch calls.  The context is an
# environment deliberately: stages are recorded by reference without passing
# mutable state through numerical helpers.

.runtime_now <- function() {
  unname(as.numeric(proc.time()[["elapsed"]]))
}

.runtime_begin <- function(method, progress = interactive(), verbose = TRUE) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  rt <- new.env(parent = emptyenv())
  rt$method <- as.character(method[[1L]])
  # `progress` is the explicit master switch for method-level status output.
  # The legacy `verbose` argument only controls the historical matching
  # notices (public wrappers suppress those whenever progress is FALSE).
  rt$enabled <- isTRUE(progress)
  rt$started <- .runtime_now()
  rt$finished <- FALSE
  rt$success <- FALSE
  rt$stages <- list()
  rt$timing <- NULL
  rt
}

.runtime_elapsed <- function(rt) {
  elapsed <- .runtime_now() - as.numeric(rt$started)
  if (!is.finite(elapsed) || elapsed < 0) 0 else unname(elapsed)
}

.runtime_format_elapsed <- function(seconds) {
  seconds <- as.numeric(seconds[[1L]])
  if (!is.finite(seconds) || seconds < 0) return("unknown time")
  if (seconds < 60) return(sprintf("%.2f s", seconds))
  if (seconds < 3600) {
    minutes <- floor(seconds / 60)
    return(sprintf("%d min %.1f s", minutes, seconds - 60 * minutes))
  }
  hours <- floor(seconds / 3600)
  remainder <- seconds - 3600 * hours
  minutes <- floor(remainder / 60)
  sprintf("%d h %d min %.1f s", hours, minutes,
          remainder - 60 * minutes)
}

.runtime_stage <- function(rt, text) {
  if (is.null(rt) || !is.environment(rt) || isTRUE(rt$finished)) {
    return(invisible(NULL))
  }
  text <- as.character(text[[1L]])
  elapsed <- .runtime_elapsed(rt)
  rt$stages[[length(rt$stages) + 1L]] <- list(
    name = text,
    elapsed = elapsed
  )
  if (isTRUE(rt$enabled)) {
    message(sprintf("[fastphylosig: %s] %s", rt$method, text))
  }
  invisible(NULL)
}

.runtime_timing <- function(rt, status = "success") {
  stages <- rt$stages
  stage_names <- if (length(stages)) {
    vapply(stages, `[[`, character(1), "name")
  } else {
    character()
  }
  stage_elapsed <- if (length(stages)) {
    vapply(stages, `[[`, numeric(1), "elapsed")
  } else {
    numeric()
  }
  names(stage_elapsed) <- stage_names
  list(
    total_elapsed = .runtime_elapsed(rt),
    units = "seconds",
    status = as.character(status[[1L]]),
    stages = stage_elapsed
  )
}

.runtime_close <- function(rt, success = TRUE) {
  if (is.null(rt) || !is.environment(rt) || isTRUE(rt$finished)) {
    return(if (is.null(rt)) NULL else rt$timing)
  }
  rt$success <- isTRUE(success)
  rt$finished <- TRUE
  timing <- .runtime_timing(rt, status = if (rt$success) "success" else "failure")
  rt$timing <- timing
  if (isTRUE(rt$enabled)) {
    if (rt$success) {
      message(sprintf(
        "[fastphylosig: %s] Done | %s",
        rt$method, .runtime_format_elapsed(timing$total_elapsed)
      ))
    } else {
      message(sprintf(
        paste0("[fastphylosig: %s] Failed | ",
               "elapsed_before_failure: %s"),
        rt$method, .runtime_format_elapsed(timing$total_elapsed)
      ))
    }
  }
  timing
}

.runtime_on_exit <- function(rt) {
  if (!is.null(rt) && is.environment(rt) && !isTRUE(rt$finished)) {
    # Let the active error propagate unchanged.  This handler only reports
    # timing; it never catches, reclassifies, or replaces the error.
    .runtime_close(rt, success = FALSE)
  }
  invisible(NULL)
}

.runtime_attach <- function(x, timing) {
  if (is.null(timing)) return(x)
  # Keep the timing as an attribute on every return type.  Scalar list
  # results additionally expose `$timing` for convenient programmatic use;
  # atomic K (no test) remains an atomic numeric and therefore uses only the
  # attribute.  Data-frame shape and legacy columns are intentionally kept
  # unchanged.
  attr(x, "timing") <- timing
  if (is.list(x) && !is.data.frame(x) && is.null(x$timing)) {
    x$timing <- timing
    attr(x, "timing") <- timing
  }
  x
}
