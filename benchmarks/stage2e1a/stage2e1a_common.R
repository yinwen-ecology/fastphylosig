# Shared helpers for the private Stage 2E1A D candidate qualification audit.
# This file is audit-only: it never edits package sources or namespace state.

stage2e1a_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

stage2e1a_restore_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

stage2e1a_parse_ints <- function(value, name, minimum = 1L) {
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum, ".",
         call. = FALSE)
  }
  unique(out)
}

stage2e1a_grid <- function(name, formal, smoke, minimum = 1L,
                           is_formal = FALSE) {
  value <- stage2e1a_env(name)
  if (nzchar(value)) stage2e1a_parse_ints(value, name, minimum) else
    if (isTRUE(is_formal)) formal else smoke
}

stage2e1a_find_library <- function() {
  lib <- stage2e1a_env("FASTPHYLOSIG_STAGE2C_LIB")
  if (nzchar(lib) && dir.exists(file.path(lib, "fastphylosig"))) return(lib)
  guesses <- Sys.glob("C:/Users/wenyi/AppData/Local/Temp/fastphylosig-stage2c-*/library")
  guesses <- guesses[dir.exists(file.path(guesses, "fastphylosig"))]
  if (length(guesses)) return(guesses[[length(guesses)]])
  installed <- find.package("fastphylosig", quiet = TRUE)
  if (nzchar(installed)) return(dirname(installed))
  stop("Set FASTPHYLOSIG_STAGE2C_LIB to the audited fastphylosig library.",
       call. = FALSE)
}

stage2e1a_compile <- function(repo, cpp_file) {
  lib <- stage2e1a_find_library()
  suppressPackageStartupMessages({
    library(fastphylosig, lib.loc = lib)
    library(ape)
    library(Rcpp)
  })
  private <- new.env(parent = globalenv())
  Rcpp::sourceCpp(normalizePath(cpp_file, winslash = "/", mustWork = TRUE),
                  env = private, rebuild = TRUE, showOutput = FALSE,
                  verbose = FALSE)
  list(env = private, lib = lib)
}

stage2e1a_make_balanced <- function(n) {
  root <- n + 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  q_node <- q_lo <- q_hi <- integer(n - 1L)
  q_node[[1L]] <- root
  q_lo[[1L]] <- 1L
  q_hi[[1L]] <- n
  head <- 1L
  tail <- 1L
  next_internal <- root + 1L
  edge_i <- 0L
  while (head <= tail) {
    node <- q_node[[head]]
    lo <- q_lo[[head]]
    hi <- q_hi[[head]]
    head <- head + 1L
    mid <- floor((lo + hi) / 2L)
    child_lo <- c(lo, mid + 1L)
    child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) child_lo[[j]] else {
        new_node <- next_internal
        next_internal <- next_internal + 1L
        tail <- tail + 1L
        q_node[[tail]] <- new_node
        q_lo[[tail]] <- child_lo[[j]]
        q_hi[[tail]] <- child_hi[[j]]
        new_node
      }
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(node, child)
    }
  }
  tree <- list(edge = edge[seq_len(edge_i), , drop = FALSE],
               tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(edge_i) / edge_i,
               Nnode = n - 1L)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

stage2e1a_make_pectinate <- function(n) {
  root <- n + 1L
  if (n == 2L) {
    edge <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  } else {
    edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
    current <- root
    next_internal <- root + 1L
    edge_i <- 0L
    for (tip in seq_len(n - 2L)) {
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(current, tip)
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(current, next_internal)
      current <- next_internal
      next_internal <- next_internal + 1L
    }
    edge_i <- edge_i + 1L
    edge[edge_i, ] <- c(current, n - 1L)
    edge_i <- edge_i + 1L
    edge[edge_i, ] <- c(current, n)
    edge <- edge[seq_len(edge_i), , drop = FALSE]
  }
  tree <- list(edge = edge, tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
               Nnode = n - 1L)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

stage2e1a_make_random <- function(n, seed = 206020L) {
  stage2e1a_restore_seed(seed + as.integer(n), {
    tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
    tree$tip.label <- paste0("sp", seq_len(n))
    tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
    tree
  })
}

stage2e1a_make_tree <- function(n, shape) {
  switch(shape,
         balanced = stage2e1a_make_balanced(n),
         random = stage2e1a_make_random(n),
         pectinate = stage2e1a_make_pectinate(n),
         stop("unknown shape: ", shape, call. = FALSE))
}

stage2e1a_make_observed <- function(tree) {
  n <- length(tree$tip.label)
  values <- as.numeric(seq_len(n) %% 2L == 0L)
  names(values) <- tree$tip.label
  values
}

stage2e1a_make_continuous <- function(tree, nsim, seed) {
  stage2e1a_restore_seed(seed, {
    matrix(stats::rnorm(length(tree$tip.label) * nsim),
           nrow = length(tree$tip.label), ncol = nsim)
  })
}

stage2e1a_capture <- function(code) {
  warning_text <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(code, warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(message = conditionMessage(e)),
                                  class = "stage2e1a_error")
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "stage2e1a_error")) {
    return(list(status = "error", value = NULL, elapsed_s = elapsed,
                warning_text = paste(unique(warning_text), collapse = " | "),
                error_message = value$message))
  }
  list(status = "ok", value = value, elapsed_s = elapsed,
       warning_text = paste(unique(warning_text), collapse = " | "),
       error_message = "")
}

stage2e1a_exact <- function(x, y) {
  isTRUE(identical(x, y))
}

# Type-7 is exact in its order-statistic selection.  The final interpolation
# can differ by one rounding step between R and the C++ scalar expression;
# keep that distinction explicit instead of using all.equal's broad default.
stage2e1a_numeric_close <- function(x, y, ulps = 32) {
  if (length(x) != length(y)) return(FALSE)
  if (!length(x)) return(TRUE)
  x <- as.numeric(x)
  y <- as.numeric(y)
  same_na <- is.na(x) & is.na(y)
  finite <- is.finite(x) & is.finite(y)
  scale <- pmax(1, abs(x), abs(y))
  all(same_na | (finite & abs(x - y) <=
                   ulps * .Machine$double.eps * scale))
}

# Mirror the currently frozen production helper exactly.  In particular,
# production returns the first order statistic when floor(h) <= 1; this is
# intentionally not replaced by stats::quantile(type = 7) in an audit.
stage2e1a_production_type7 <- function(values, probability) {
  values <- sort(as.numeric(values))
  n <- length(values)
  if (n < 1L) return(NaN)
  h <- 1 + (n - 1) * probability
  lo <- floor(h)
  gamma <- h - lo
  if (lo <= 1) return(values[[1L]])
  if (lo >= n) return(values[[n]])
  (1 - gamma) * values[[lo]] + gamma * values[[lo + 1L]]
}

stage2e1a_compare <- function(check, lhs, rhs, details = "") {
  data.frame(check = check, pass = stage2e1a_exact(lhs, rhs),
             details = details, stringsAsFactors = FALSE)
}

stage2e1a_write_provenance <- function(repo, out_dir, cpp_file, args) {
  git <- function(...) {
    value <- tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE),
                      error = function(e) character())
    if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
  }
  lines <- c(
    paste0("repo=", repo),
    paste0("git_head=", git("-C", repo, "rev-parse", "HEAD")),
    paste0("git_status=", git("-C", repo, "status", "--short")),
    paste0("harness=", normalizePath(cpp_file, winslash = "/", mustWork = TRUE)),
    paste0("stage=Stage 2E1A Final D Candidate Qualification"),
    paste0("production_translation_unit=src/d_stream.cpp"),
    paste0("production_contract=fast_d binary streaming path; R sample path retained for polytomy"),
    paste0("formal=", if (isTRUE(args$formal)) "YES" else "NO_SMOKE"),
    paste0("generated_at=", format(Sys.time(), tz = "UTC"))
  )
  writeLines(lines, file.path(out_dir, "stage2e1a_provenance.txt"))
  invisible(lines)
}
