# fastphylosig 0.2.0 final release scientific validation
#
# This script is deliberately limited to release-contract and scientific
# parity checks.  It does not modify package sources, tests, or documentation.

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

repo <- normalizePath(
  arg_value("--repo", getwd()), winslash = "/", mustWork = TRUE
)
source_commit <- arg_value(
  "--source-commit", Sys.getenv("FASTPHYLOSIG_SOURCE_COMMIT", "not supplied")
)
out_dir <- normalizePath(
  arg_value("--out", file.path(repo, "benchmarks", "final_release",
                                 "scientific")),
  winslash = "/", mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The child mode is used for the saveRDS/readRDS gate.  It is kept in this
# same file so a fresh R process can rerun the exact checked code path.
child_rds <- arg_value("--child-rds", NULL)
child_out <- arg_value("--child-out", NULL)

if (nzchar(Sys.getenv("FASTPHYLOSIG_LIBRARY"))) {
  .libPaths(unique(c(Sys.getenv("FASTPHYLOSIG_LIBRARY"), .libPaths())))
}

suppressPackageStartupMessages(library(fastphylosig))
suppressPackageStartupMessages(library(ape))

if (!is.null(child_rds)) {
  log_path <- if (is.null(child_out)) {
    file.path(dirname(child_rds), "fresh_session.log")
  } else {
    paste0(child_out, ".log")
  }
  log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
  sink(log_connection, type = "output", split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
    try(close(log_connection), silent = TRUE)
  }, add = TRUE)

  ctx <- readRDS(child_rds)
  tips <- ctx$tree$tip.label
  continuous <- stats::setNames(seq_along(tips) / 3, tips)
  binary <- stats::setNames(rep(c(0L, 0L, 1L, 1L), length.out =
                                length(tips)), tips)
  k <- fast_k(ctx, continuous, test = FALSE, verbose = FALSE,
              progress = FALSE)
  lambda <- fast_lambda(ctx, continuous, lambda_profile = FALSE,
                        verbose = FALSE, progress = FALSE)
  d <- fast_d(ctx, binary, test = FALSE, return_sim = FALSE,
              verbose = FALSE, progress = FALSE)
  mutation_message <- tryCatch({
    ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 0.25
    fast_k(ctx, continuous, test = FALSE, verbose = FALSE, progress = FALSE)
    "NO_ERROR"
  }, error = function(e) conditionMessage(e))
  result <- list(
    context_schema_version = ctx$context_schema_version,
    canonical_contract_version = ctx$canonical_contract_version,
    protected_snapshot_version = ctx$protected_snapshot_version,
    k = as.numeric(k),
    lambda = unname(lambda$lambda),
    lambda_logL = unname(lambda$logL),
    d = unname(d$DEstimate),
    mutation_message = mutation_message
  )
  if (!is.null(child_out)) saveRDS(result, child_out, version = 2L)
  cat("FRESH_SESSION_OK\n")
  quit(save = "no", status = 0L)
}

log_path <- file.path(out_dir, "scientific_validation.log")
log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

cat("fastphylosig 0.2.0 scientific validation\n")
cat("Started:", format(Sys.time(), tz = "UTC"), "UTC\n")

records <- data.frame(
  id = character(), status = character(), warning_expected = logical(),
  warnings = character(), detail = character(), stringsAsFactors = FALSE
)

capture_call <- function(fun) {
  warnings <- character()
  err <- NULL
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  list(value = value, warnings = warnings, error = err)
}

same_value <- function(a, b, tolerance = 1e-10) {
  if (identical(a, b)) return(TRUE)
  if (length(a) != length(b)) return(FALSE)
  if (is.raw(a) || is.raw(b)) return(identical(a, b))
  if (is.numeric(a) && is.numeric(b)) {
    na <- is.na(a)
    nb <- is.na(b)
    if (!identical(na, nb)) return(FALSE)
    keep <- !na
    if (!any(keep)) return(TRUE)
    scale <- pmax(1, abs(a[keep]), abs(b[keep]))
    return(all(abs(a[keep] - b[keep]) <= tolerance * scale))
  }
  isTRUE(all.equal(a, b, check.attributes = TRUE, tolerance = tolerance))
}

field_parity <- function(a, b, fields, tolerance = 1e-10) {
  if (!all(fields %in% names(a)) || !all(fields %in% names(b))) {
    return(FALSE)
  }
  all(vapply(fields, function(field) {
    same_value(a[[field]], b[[field]], tolerance)
  }, logical(1)))
}

record_case <- function(id, fun, check, warning_expected = FALSE) {
  captured <- capture_call(fun)
  ok <- is.null(captured$error) && isTRUE(check(captured$value))
  status <- if (ok) "PASS" else "FAIL"
  detail <- if (!is.null(captured$error)) {
    paste0("error: ", captured$error)
  } else if (!ok) {
    "check returned FALSE"
  } else {
    ""
  }
  warning_text <- paste(captured$warnings, collapse = " | ")
  records <<- rbind(
    records,
    data.frame(
      id = id, status = status, warning_expected = warning_expected,
      warnings = warning_text, detail = detail,
      stringsAsFactors = FALSE
    )
  )
  cat(sprintf("[%s] %s", status, id))
  if (nzchar(detail)) cat(" -- ", detail, sep = "")
  if (nzchar(warning_text)) cat(" -- warning: ", warning_text, sep = "")
  cat("\n")
  invisible(captured$value)
}

record_skip <- function(id, reason) {
  records <<- rbind(
    records,
    data.frame(id = id, status = "SKIP", warning_expected = FALSE,
               warnings = "", detail = reason, stringsAsFactors = FALSE)
  )
  cat(sprintf("[SKIP] %s -- %s\n", id, reason))
}

tree_snapshot <- function(x) serialize(x, connection = NULL, version = 3L)

make_tree <- function(n = 8L, shape = c("balanced", "random", "pectinate"),
                      heterogeneous = FALSE) {
  shape <- match.arg(shape)
  set.seed(20260913L + n + match(shape, c("balanced", "random", "pectinate")))
  if (shape == "balanced") {
    # ape::stree() is most stable for powers of two in a release fixture.
    n_use <- 2^ceiling(log2(n))
    tr <- ape::stree(n_use, type = "balanced")
    if (n_use != n) {
      tr <- ape::drop.tip(tr, tr$tip.label[seq.int(n + 1L, n_use)])
    }
  } else if (shape == "pectinate") {
    tr <- ape::stree(n, type = "left")
  } else {
    tr <- ape::rtree(n)
  }
  tr$tip.label <- paste0("sp", seq_len(ape::Ntip(tr)))
  tr$edge.length <- if (heterogeneous) {
    10^seq(-2, 1, length.out = nrow(tr$edge))
  } else {
    pmax(as.numeric(tr$edge.length), 1e-4)
  }
  ape::reorder.phylo(tr, "cladewise")
}

make_four_tree <- function(labels = c("a", "b", "c", "d")) {
  tr <- list(
    edge = matrix(c(7L, 5L, 7L, 6L, 5L, 1L, 5L, 2L,
                   6L, 3L, 6L, 4L), ncol = 2L, byrow = TRUE),
    edge.length = c(10.1, 20.2, 1.1, 2.2, 3.3, 4.4),
    Nnode = 3L, tip.label = labels
  )
  class(tr) <- "phylo"
  tr
}

renumber_internals <- function(tree) {
  n <- length(tree$tip.label)
  internal <- sort(unique(as.integer(tree$edge[tree$edge > n])))
  replacement <- n + 10L + (seq_along(internal) - 1L) * 7L
  map <- stats::setNames(replacement, as.character(internal))
  out <- tree
  mapped <- out$edge
  hit <- out$edge > n
  mapped[hit] <- unname(map[as.character(out$edge[hit])])
  out$edge <- matrix(as.integer(mapped), ncol = 2L,
                     dimnames = dimnames(tree$edge))
  out
}

shuffle_edges <- function(tree) {
  order_rows <- rev(seq_len(nrow(tree$edge)))
  out <- tree
  out$edge <- tree$edge[order_rows, , drop = FALSE]
  out$edge.length <- tree$edge.length[order_rows]
  out
}

make_v2_tree <- function() {
  tr <- list(
    edge = matrix(c(
      11L, 9L, 11L, 10L, 9L, 7L, 9L, 8L,
      7L, 1L, 7L, 2L, 8L, 3L, 8L, 4L,
      10L, 5L, 10L, 6L
    ), ncol = 2L, byrow = TRUE),
    edge.length = c(0.8, 0.9, 0.7, 0.6, 0.5,
                    0.4, 0.3, 0.2, 1.1, 1.2),
    Nnode = 5L,
    tip.label = c("a", "b", "c", "d", "e", "f"),
    root.edge = 0.125
  )
  class(tr) <- "phylo"
  tr
}

make_continuous <- function(tree) {
  n <- length(tree$tip.label)
  stats::setNames(sin(seq_len(n) / 2) + seq_len(n) / n,
                  tree$tip.label)
}

make_binary <- function(tree) {
  stats::setNames(rep(c(0L, 0L, 1L, 1L, 1L, 0L),
                      length.out = length(tree$tip.label)),
                  tree$tip.label)
}

make_categorical <- function(tree) {
  stats::setNames(factor(rep(c("A", "B", "C"),
                             length.out = length(tree$tip.label)),
                         levels = c("A", "B", "C")),
                  tree$tip.label)
}

make_permutations <- function(n, nsim = 3L) {
  base <- seq_len(n)
  out <- matrix(base, nrow = nsim, ncol = n, byrow = TRUE)
  if (nsim >= 2L) out[2L, ] <- rev(base)
  if (nsim >= 3L) out[3L, ] <- c(base[2:n], base[[1L]])
  if (nsim >= 4L) out[4L, ] <- c(base[3:n], base[1:2])
  out
}

call_k <- function(tree, x, test = FALSE, nsim = 4L,
                   permutations = NULL, return_sim = FALSE) {
  fast_k(tree, x, test = test, nsim = nsim, permutations = permutations,
         return_sim = return_sim, keep_null = return_sim,
         verbose = FALSE, progress = FALSE, ncores = 1L)
}

call_lambda <- function(tree, x, test = FALSE) {
  fast_lambda(tree, x, test = test, lambda_profile = FALSE,
              verbose = FALSE, progress = FALSE, ncores = 1L)
}

call_d <- function(tree, x, test = FALSE, nsim = 3L,
                   random_states = NULL, brownian_states = NULL,
                   return_sim = FALSE) {
  fast_d(tree, x, test = test, nsim = nsim,
         random_states = random_states, brownian_states = brownian_states,
         return_sim = return_sim, keep_null = return_sim,
         verbose = FALSE, progress = FALSE, ncores = 1L)
}

call_delta <- function(tree, x, test = FALSE, nsim = 3L,
                       permutations = NULL, return_sim = FALSE) {
  fast_delta(tree, x, test = test, nsim = nsim,
             permutations = permutations, return_sim = return_sim,
             mcmc_sim = 24L, thin = 4L, burn = 8L, model = "ER",
             verbose = FALSE, progress = FALSE, ncores = 1L)
}

normal <- make_tree(8L, "balanced", heterogeneous = TRUE)
normal_x <- make_continuous(normal)
normal_binary <- make_binary(normal)
normal_cat <- make_categorical(normal)
normal_ctx <- prepare_tree(normal)

four <- make_four_tree()
four_x <- stats::setNames(c(0.3, 1.7, -0.8, 0.4), four$tip.label)
four_binary <- stats::setNames(c(0L, 0L, 1L, 1L), four$tip.label)
four_ctx <- prepare_tree(four)

cat("\n-- K contract --\n")
k_perms <- make_permutations(length(normal$tip.label), nsim = 4L)
k_raw <- record_case(
  "K estimate/P/MCSE and raw-prepared parity",
  function() {
    raw <- call_k(normal, normal_x, test = TRUE, nsim = nrow(k_perms),
                  permutations = k_perms, return_sim = TRUE)
    prepared <- call_k(normal_ctx, normal_x, test = TRUE,
                       nsim = nrow(k_perms), permutations = k_perms,
                       return_sim = TRUE)
    list(raw = raw, prepared = prepared)
  },
  function(z) {
    fields <- c("K", "P", "MCSE_P", "exceedance_count", "sim.K")
    if (!field_parity(z$raw, z$prepared, fields, 1e-12)) return(FALSE)
    if (!all(is.finite(c(z$raw$K, z$raw$P, z$raw$MCSE_P)))) return(FALSE)
    sim <- as.numeric(z$raw$sim.K)
    obs <- as.numeric(z$raw$K)
    scale <- pmax(1, abs(sim), abs(obs))
    inclusive <- sim >= obs | (is.finite(sim) & is.finite(obs) &
      obs - sim <= 8 * .Machine$double.eps * scale)
    expected <- sum(inclusive)
    identical(as.integer(z$raw$exceedance_count), as.integer(expected)) &&
      same_value(z$raw$P, expected / length(sim), 1e-14) &&
      same_value(z$raw$MCSE_P,
                 sqrt(z$raw$P * (1 - z$raw$P) / length(sim)), 1e-14)
  }
)

record_case(
  "K inclusive upper tail and thread-invariant controlled null",
  function() {
    one <- call_k(normal, normal_x, test = TRUE, nsim = nrow(k_perms),
                  permutations = k_perms, return_sim = TRUE)
    two <- fast_k(normal, normal_x, test = TRUE, nsim = nrow(k_perms),
                  permutations = k_perms, return_sim = TRUE,
                  keep_null = TRUE, verbose = FALSE, progress = FALSE,
                  ncores = 2L)
    list(one = one, two = two)
  },
  function(z) field_parity(z$one, z$two,
                           c("K", "P", "MCSE_P", "exceedance_count", "sim.K"),
                           1e-12)
)

record_case(
  "K large-offset and trait-wise NA accounting",
  function() {
    X <- cbind(
      large_offset = normal_x + 1e12,
      with_na = normal_x
    )
    X[2L, 2L] <- NA_real_
    rownames(X) <- normal$tip.label
    raw <- fast_k(normal, X, test = FALSE, verbose = FALSE, progress = FALSE)
    prepared <- fast_k(normal_ctx, X, test = FALSE,
                       verbose = FALSE, progress = FALSE)
    list(raw = raw, prepared = prepared)
  },
  function(z) {
    field_parity(z$raw, z$prepared,
                 c("K_fast", "n_species", "n_removed_na", "status"), 1e-10) &&
      all(z$raw$n_species == c(8L, 7L)) &&
      all(is.finite(z$raw$K_fast))
  }
)

cat("\n-- lambda contract --\n")
record_case(
  "lambda estimate/logLik/LR/P and raw-prepared parity",
  function() {
    raw <- call_lambda(normal, normal_x, test = TRUE)
    prepared <- call_lambda(normal_ctx, normal_x, test = TRUE)
    list(raw = raw, prepared = prepared)
  },
  function(z) {
    fields <- c("lambda", "logL", "gls_mean", "sig2", "logL0", "LR", "P")
    field_parity(z$raw, z$prepared, fields, 1e-10) &&
      all(is.finite(unlist(z$raw[fields])))
  }
)

record_case(
  "lambda NA matching and immutable continuous input",
  function() {
    tr_before <- tree_snapshot(normal)
    x_before <- tree_snapshot(normal_x)
    x <- normal_x
    x[3L] <- NA_real_
    raw <- call_lambda(normal, x, test = FALSE)
    prepared <- call_lambda(normal_ctx, x, test = FALSE)
    list(raw = raw, prepared = prepared, tree_same = identical(
      tree_snapshot(normal), tr_before), x_same = identical(tree_snapshot(x),
                                                             tree_snapshot(x_before)))
  },
  function(z) {
    field_parity(z$raw, z$prepared,
                 c("lambda", "logL", "gls_mean", "sig2", "status"), 1e-9) &&
      isTRUE(z$tree_same) && !isTRUE(z$x_same) &&
      is.finite(z$raw$lambda)
  }
)

cat("\n-- D contract --\n")
d_random <- cbind(
  four_binary,
  c(1L, 0L, 1L, 0L),
  c(0L, 1L, 1L, 0L)
)
d_brownian <- cbind(
  four_binary,
  c(1L, 1L, 0L, 0L),
  c(0L, 1L, 0L, 1L)
)

d_controlled <- record_case(
  "D estimate/random-Brownian P/MCSE/strict tails/polytomy",
  function() {
    raw <- call_d(four, four_binary, test = TRUE, nsim = 3L,
                  random_states = d_random, brownian_states = d_brownian,
                  return_sim = TRUE)
    prepared <- call_d(four_ctx, four_binary, test = TRUE, nsim = 3L,
                       random_states = d_random, brownian_states = d_brownian,
                       return_sim = TRUE)
    obs_fit <- call_d(four, four_binary, test = FALSE)
    d_for_states <- function(states) vapply(seq_len(ncol(states)), function(j) {
      call_d(four, stats::setNames(states[, j], four$tip.label),
             test = FALSE)$Parameters$Observed
    }, numeric(1))
    list(raw = raw, prepared = prepared,
         observed = raw$Parameters$Observed,
         random_d = d_for_states(d_random), brownian_d = d_for_states(d_brownian))
  },
  function(z) {
    fields <- c("DEstimate", "Pval1", "Pval0", "P_random", "P_Brownian",
                "MCSE_P_random", "MCSE_P_Brownian", "Parameters",
                "nsim_requested", "nsim_successful_random",
                "nsim_successful_brownian", "status")
    if (!field_parity(z$raw, z$prepared, fields, 1e-12)) return(FALSE)
    if (!all(is.finite(c(z$raw$DEstimate, z$raw$Pval1, z$raw$Pval0,
                         z$raw$MCSE_P_random, z$raw$MCSE_P_Brownian)))) {
      return(FALSE)
    }
    p_random <- mean(z$random_d < z$observed)
    p_brownian <- mean(z$brownian_d > z$observed)
    same_value(z$raw$Pval1, p_random, 1e-14) &&
      same_value(z$raw$Pval0, p_brownian, 1e-14) &&
      same_value(z$raw$MCSE_P_random,
                 sqrt(p_random * (1 - p_random) / 3), 1e-14) &&
      same_value(z$raw$MCSE_P_Brownian,
                 sqrt(p_brownian * (1 - p_brownian) / 3), 1e-14) &&
      z$raw$nsim_successful_random == 3L &&
      z$raw$nsim_successful_brownian == 3L
  }
)

record_case(
  "D input immutability and raw-prepared controlled-state parity",
  function() {
    tree_before <- tree_snapshot(four)
    x_before <- tree_snapshot(four_binary)
    raw <- call_d(four, four_binary, test = TRUE, nsim = 3L,
                  random_states = d_random, brownian_states = d_brownian,
                  return_sim = TRUE)
    prepared <- call_d(four_ctx, four_binary, test = TRUE, nsim = 3L,
                       random_states = d_random, brownian_states = d_brownian,
                       return_sim = TRUE)
    list(raw = raw, prepared = prepared,
         tree_same = identical(tree_snapshot(four), tree_before),
         x_same = identical(tree_snapshot(four_binary), x_before))
  },
  function(z) {
    field_parity(z$raw, z$prepared,
                 c("DEstimate", "Pval1", "Pval0", "Parameters", "Permutations"),
                 1e-12) && isTRUE(z$tree_same) && isTRUE(z$x_same)
  }
)

record_case(
  "D rooted polytomy check_tree readiness",
  function() {
    # Use the fixed 20-tip polytomy fixture from the strict-tail audit.  The
    # four-tip toy state can have an undefined contrast denominator even when
    # the tree is structurally ready, so it is not a valid estimator gate.
    branch <- function(text) paste0(text, ":1")
    binary_subtree <- function(labels) {
      if (length(labels) == 1L) return(branch(labels))
      branch(paste0("(", branch(labels[[1L]]), ",",
                    binary_subtree(labels[-1L]), ")"))
    }
    labels <- paste0("t", seq_len(20L))
    polytomy_text <- branch(paste0(
      "(", paste(c(branch(labels[[1L]]), branch(labels[[2L]]),
                    binary_subtree(labels[3L:16L])), collapse = ","), ")"
    ))
    polytomy <- ape::read.tree(text = paste0(
      "(", polytomy_text, ",", binary_subtree(labels[17L:20L]), ");"
    ))
    binary <- stats::setNames(integer(20L), labels)
    binary[c(15L, 8L, 13L, 16L)] <- 1L
    fit <- call_d(polytomy, binary, test = FALSE)
    list(check = check_tree(polytomy, signal = "D"), fit = fit)
  },
  function(z) isTRUE(z$check$ready_by_signal[["D"]]) &&
    !is.null(z$fit) && is.finite(z$fit$DEstimate)
)

cat("\n-- Delta and ACE contract --\n")
delta_tree <- make_tree(20L, "random", heterogeneous = TRUE)
delta_x <- make_categorical(delta_tree)
delta_ctx <- prepare_tree(delta_tree)
delta_perms <- make_permutations(length(delta_tree$tip.label), nsim = 3L)

record_case(
  "Delta estimate/MCMC diagnostics/P/MCSE/status",
  function() {
    set.seed(2026091301L)
    fit <- call_delta(delta_tree, delta_x, test = TRUE, nsim = 3L,
                      permutations = delta_perms, return_sim = TRUE)
    fit_prepared <- {
      set.seed(2026091301L)
      call_delta(delta_ctx, delta_x, test = TRUE, nsim = 3L,
                 permutations = delta_perms, return_sim = TRUE)
    }
    list(fit = fit, prepared = fit_prepared)
  },
  function(z) {
    fields <- c("delta", "alpha_mean", "beta_mean", "n_saved",
                "n_saved_successful", "ESS_alpha", "ESS_beta",
                "split_Rhat_alpha", "split_Rhat_beta", "MCSE_Delta",
                "diagnostics_available", "requested_simulations",
                "successful_simulations", "P", "P_MCSE", "sim.delta",
                "status")
    field_parity(z$fit, z$prepared, fields, 1e-10) &&
      is.finite(z$fit$delta) && is.finite(z$fit$MCSE_Delta) &&
      is.finite(z$fit$P) && is.finite(z$fit$P_MCSE) &&
      z$fit$successful_simulations == 3L &&
      z$fit$status %in% c("ok", "partial") &&
      same_value(z$fit$P_MCSE,
                 sqrt(z$fit$P * (1 - z$fit$P) / z$fit$successful_simulations),
                 1e-12)
  },
  warning_expected = TRUE
)

# Keep the release parity fixture aligned with tests/testthat/test-fast-signal.R:
# this is a stochastic ER draw at seed 31 with ordinary branch lengths.  The
# deliberately ill-conditioned cyclic/high-heterogeneity case is audited
# below as an explicit expected-warning case instead of being used as a
# zero-warning release gate.
set.seed(31L)
ace_tree <- ape::rtree(40L)
ace_x <- ape::rTraitDisc(
  ace_tree, model = "ER", k = 3, rate = 1.5,
  states = letters[1:3]
)
ace_x <- factor(ace_x[ace_tree$tip.label])
ace_ctx <- prepare_tree(ace_tree)
record_case(
  "ACE logLik/rates/SE/ancestral likelihoods and raw-prepared parity",
  function() {
    raw <- fast_ace(ace_x, phy = ace_tree, model = "ER", CI = TRUE,
                    marginal = TRUE, progress = FALSE)
    prepared <- fast_ace(ace_x, prepared = ace_ctx, model = "ER", CI = TRUE,
                         marginal = TRUE, progress = FALSE)
    list(raw = raw, prepared = prepared)
  },
  function(z) {
    fields <- c("loglik", "rates", "se", "index.matrix", "convergence",
                "message", "lik.anc")
    field_parity(z$raw, z$prepared, fields, 1e-9) &&
      is.finite(z$raw$loglik) && all(is.finite(z$raw$rates)) &&
      length(z$raw$se) == length(z$raw$rates) &&
      all(is.finite(z$raw$se) | is.na(z$raw$se)) &&
      is.matrix(z$raw$lik.anc) && all(is.finite(z$raw$lik.anc))
  }
)

record_case(
  "ACE CI=false likelihood-only path has no spurious warning",
  function() fast_ace(ace_x, phy = ace_tree, model = "ER", CI = FALSE,
                      marginal = FALSE, progress = FALSE),
  function(z) is.finite(z$loglik) && all(is.finite(z$rates))
)

ace_pathological_tree <- make_tree(40L, "random", heterogeneous = TRUE)
ace_pathological_x <- make_categorical(ace_pathological_tree)
record_case(
  "ACE ill-conditioned fixture reports explicit nonconvergence",
  function() fast_ace(
    ace_pathological_x, phy = ace_pathological_tree, model = "ER",
    CI = FALSE, marginal = FALSE, progress = FALSE
  ),
  function(z) is.finite(z$loglik) && all(is.finite(z$rates)) &&
    isTRUE(z$convergence != 0L) &&
    identical(z$message, "singular convergence (7)"),
  warning_expected = TRUE
)

ace_warning_provenance <- local({
  fast_ci <- capture_call(function() fast_ace(
    ace_x, phy = ace_tree, model = "ER", CI = TRUE,
    marginal = TRUE, progress = FALSE
  ))
  fast_no_ci <- capture_call(function() fast_ace(
    ace_x, phy = ace_tree, model = "ER", CI = FALSE,
    marginal = FALSE, progress = FALSE
  ))
  ape_ci <- capture_call(function() ape::ace(
    ace_x, ace_tree, type = "discrete", method = "ML", model = "ER"
  ))
  pathological_no_ci <- capture_call(function() fast_ace(
    ace_pathological_x, phy = ace_pathological_tree, model = "ER",
    CI = FALSE, marginal = FALSE, progress = FALSE
  ))
  cat("[ACE-WARNING-PROVENANCE] fast_ace_CI_TRUE=",
      paste(fast_ci$warnings, collapse = " | "),
      "; fast_ace_CI_FALSE=", paste(fast_no_ci$warnings, collapse = " | "),
      "; ape_ace=", paste(ape_ci$warnings, collapse = " | "),
      "; pathological_fast_ace_CI_FALSE=",
      paste(pathological_no_ci$warnings, collapse = " | "), "\n", sep = "")
  error_text <- function(x) {
    if (is.null(x$error)) "" else x$error
  }
  convergence_value <- function(x) {
    if (is.null(x$value) || is.null(x$value$convergence)) NA_integer_ else
      as.integer(x$value$convergence)
  }
  message_value <- function(x) {
    if (is.null(x$value) || is.null(x$value$message)) "" else
      as.character(x$value$message)
  }
  fit_value <- function(x, field) {
    if (is.null(x$value) || is.null(x$value[[field]])) NA_real_ else
      as.numeric(x$value[[field]][[1L]])
  }
  data.frame(
    path = c("stable_fast_ace_CI_TRUE", "stable_fast_ace_CI_FALSE",
             "stable_ape_ace", "pathological_fast_ace_CI_FALSE"),
    warnings = c(paste(fast_ci$warnings, collapse = " | "),
                 paste(fast_no_ci$warnings, collapse = " | "),
                 paste(ape_ci$warnings, collapse = " | "),
                 paste(pathological_no_ci$warnings, collapse = " | ")),
    error = c(error_text(fast_ci), error_text(fast_no_ci),
              error_text(ape_ci), error_text(pathological_no_ci)),
    convergence = c(convergence_value(fast_ci), convergence_value(fast_no_ci),
                    NA_integer_, convergence_value(pathological_no_ci)),
    message = c(message_value(fast_ci), message_value(fast_no_ci), "",
                message_value(pathological_no_ci)),
    loglik = c(fit_value(fast_ci, "loglik"), fit_value(fast_no_ci, "loglik"),
               fit_value(ape_ci, "loglik"), fit_value(pathological_no_ci, "loglik")),
    rate = c(fit_value(fast_ci, "rates"), fit_value(fast_no_ci, "rates"),
             fit_value(ape_ci, "rates"), fit_value(pathological_no_ci, "rates")),
    stringsAsFactors = FALSE
  )
})
write.csv(ace_warning_provenance,
          file.path(out_dir, "ace_warning_provenance.csv"), row.names = FALSE,
          na = "")

cat("\n-- V2 representation and retained-field contract --\n")
v2 <- make_v2_tree()
v2_variants <- list(
  original = v2,
  renumbered = renumber_internals(v2),
  shuffled = shuffle_edges(v2),
  renumbered_shuffled = shuffle_edges(renumber_internals(v2))
)
v2_x <- make_continuous(v2)
v2_binary <- make_binary(v2)
v2_cat <- make_categorical(v2)
v2_random <- cbind(v2_binary, rev(v2_binary),
                   c(0L, 1L, 1L, 0L, 1L, 0L))
v2_brownian <- cbind(v2_binary, c(1L, 1L, 0L, 0L, 1L, 0L),
                     c(0L, 1L, 0L, 1L, 0L, 1L))
v2_perms <- make_permutations(length(v2$tip.label), nsim = 3L)

record_case(
  "V2 renumbering/edge-shuffle fingerprint and branch association",
  function() {
    canonical <- lapply(v2_variants, function(tr) {
      ct <- fastphylosig:::.safe_canonicalize_core(tr)
      list(edge = unname(ct$edge), edge.length = unname(ct$edge.length),
           tip.label = ct$tip.label, Nnode = ct$Nnode,
           fingerprint = fastphylosig:::.tree_fingerprint(ct, .canonical = TRUE))
    })
    canonical
  },
  function(z) all(vapply(z[-1L], function(x) {
    same_value(x$edge, z[[1L]]$edge, 0) &&
      same_value(x$edge.length, z[[1L]]$edge.length, 0) &&
      identical(x$tip.label, z[[1L]]$tip.label) &&
      identical(x$Nnode, z[[1L]]$Nnode) &&
      identical(x$fingerprint, z[[1L]]$fingerprint)
  }, logical(1)))
)

record_case(
  "V2 K/lambda/D/Delta/ACE representation propagation",
  function() {
    result <- lapply(v2_variants, function(tr) {
      ctx <- prepare_tree(tr)
      k <- call_k(tr, v2_x, test = TRUE, nsim = 3L,
                  permutations = v2_perms, return_sim = TRUE)
      kp <- call_k(ctx, v2_x, test = TRUE, nsim = 3L,
                   permutations = v2_perms, return_sim = TRUE)
      lambda <- call_lambda(tr, v2_x, test = TRUE)
      lambdap <- call_lambda(ctx, v2_x, test = TRUE)
      d <- call_d(tr, v2_binary, test = TRUE, nsim = 3L,
                  random_states = v2_random, brownian_states = v2_brownian,
                  return_sim = TRUE)
      dp <- call_d(ctx, v2_binary, test = TRUE, nsim = 3L,
                   random_states = v2_random, brownian_states = v2_brownian,
                   return_sim = TRUE)
      set.seed(2026091302L)
      delta <- call_delta(tr, v2_cat, test = TRUE, nsim = 3L,
                          permutations = v2_perms, return_sim = TRUE)
      set.seed(2026091302L)
      deltap <- call_delta(ctx, v2_cat, test = TRUE, nsim = 3L,
                           permutations = v2_perms, return_sim = TRUE)
      ace <- fast_ace(v2_cat, phy = tr, model = "ER", CI = FALSE,
                      marginal = TRUE, progress = FALSE)
      acep <- fast_ace(v2_cat, prepared = ctx, model = "ER", CI = FALSE,
                       marginal = TRUE, progress = FALSE)
      list(k = k, kp = kp, lambda = lambda, lambdap = lambdap,
           d = d, dp = dp, delta = delta, deltap = deltap,
           ace = ace, acep = acep)
    })
    result
  },
  function(z) {
    checks <- vapply(z, function(x) {
      kcheck <- field_parity(x$k, x$kp,
                             c("K", "P", "MCSE_P", "exceedance_count", "sim.K"),
                             1e-11)
      lcheck <- field_parity(x$lambda, x$lambdap,
                             c("lambda", "logL", "logL0", "LR", "P"), 1e-8)
      dcheck <- field_parity(x$d, x$dp,
                             c("DEstimate", "Pval1", "Pval0", "Parameters",
                               "Permutations"), 1e-11)
      decheck <- field_parity(x$delta, x$deltap,
                              c("delta", "alpha_mean", "beta_mean", "n_saved",
                                "ESS_alpha", "ESS_beta", "MCSE_Delta", "P", "P_MCSE",
                                "sim.delta", "status"), 1e-8)
      acecheck <- field_parity(x$ace, x$acep,
                               c("loglik", "rates"), 1e-8) &&
        ((is.null(x$ace$lik.anc) && is.null(x$acep$lik.anc)) ||
         same_value(x$ace$lik.anc, x$acep$lik.anc, 1e-8))
      x$checks <- c(k = kcheck, lambda = lcheck, D = dcheck,
                    Delta = decheck, ACE = acecheck)
      all(x$checks)
    }, logical(1))
    all(checks)
  },
  warning_expected = TRUE
)

unicode_tree <- make_four_tree(c("alpha", "b\u000dchar", "\u03b1", "\u4e2d"))
unicode_x <- stats::setNames(c(1.2, -0.7, 2.3, 0.4), unicode_tree$tip.label)
record_case(
  "Unicode/delimiter labels and heterogeneous branch association",
  function() {
    tr2 <- shuffle_edges(renumber_internals(unicode_tree))
    a <- fastphylosig:::.safe_canonicalize_core(unicode_tree)
    b <- fastphylosig:::.safe_canonicalize_core(tr2)
    ka <- call_k(unicode_tree, unicode_x)
    kb <- call_k(tr2, unicode_x)
    la <- call_lambda(unicode_tree, unicode_x)
    lb <- call_lambda(tr2, unicode_x)
    list(a = a, b = b, ka = ka, kb = kb, la = la, lb = lb)
  },
  function(z) {
    same_value(z$a$edge, z$b$edge, 0) &&
      same_value(z$a$edge.length, z$b$edge.length, 0) &&
      identical(fastphylosig:::.tree_fingerprint(z$a, .canonical = TRUE),
                fastphylosig:::.tree_fingerprint(z$b, .canonical = TRUE)) &&
      same_value(z$ka, z$kb, 1e-9) &&
      field_parity(z$la, z$lb, c("lambda", "logL", "gls_mean", "sig2"), 1e-8)
  }
)

cat("\n-- NA, root, and failure contracts --\n")
record_case(
  "raw/prepared NA parity across K/lambda/D/Delta",
  function() {
    X <- cbind(a = normal_x, b = normal_x)
    X[2L, 1L] <- NA_real_
    rownames(X) <- normal$tip.label
    k <- fast_k(normal, X, test = FALSE, verbose = FALSE, progress = FALSE)
    kp <- fast_k(normal_ctx, X, test = FALSE,
                 verbose = FALSE, progress = FALSE)
    l <- call_lambda(normal, X[, 1L])
    lp <- call_lambda(normal_ctx, X[, 1L])
    db <- normal_binary
    db[4L] <- NA_integer_
    db_keep <- db[!is.na(db)]
    db_random <- cbind(db_keep, rev(db_keep),
                       c(0L, 1L, 1L, 0L, 1L, 0L, 1L))
    db_brownian <- cbind(db_keep, c(1L, 1L, 0L, 0L, 1L, 0L, 0L),
                         c(0L, 1L, 0L, 1L, 0L, 1L, 1L))
    d <- call_d(normal, db, test = TRUE, nsim = 3L,
                random_states = db_random, brownian_states = db_brownian,
                return_sim = TRUE)
    dp <- call_d(normal_ctx, db, test = TRUE, nsim = 3L,
                 random_states = db_random, brownian_states = db_brownian,
                 return_sim = TRUE)
    dc <- normal_cat
    delta_na <- delta_x
    delta_na[5L] <- NA
    set.seed(2026091303L)
    delta <- call_delta(delta_tree, delta_na, test = FALSE)
    set.seed(2026091303L)
    deltap <- call_delta(delta_ctx, delta_na, test = FALSE)
    list(k = k, kp = kp, l = l, lp = lp, d = d, dp = dp,
         delta = delta, deltap = deltap)
  },
  function(z) {
    parts <- c(
      K = field_parity(z$k, z$kp,
                       c("K_fast", "n_species", "n_removed_na", "status"), 1e-9),
      lambda = field_parity(z$l, z$lp,
                            c("lambda", "logL", "gls_mean", "sig2", "status"),
                            1e-8),
      D = field_parity(z$d, z$dp,
                       c("DEstimate", "Pval1", "Pval0", "status"),
                       1e-9),
      Delta = field_parity(z$delta, z$deltap,
                           c("delta", "n_species", "status"), 1e-8),
      Delta_n = identical(z$delta$n_species, 19L)
    )
    all(parts)
  },
  warning_expected = TRUE
)

record_case(
  "large representative raw/prepared K and lambda parity",
  function() {
    tr <- make_tree(128L, "random", heterogeneous = TRUE)
    x <- make_continuous(tr)
    ctx <- prepare_tree(tr)
    list(k = call_k(tr, x), kp = call_k(ctx, x),
         l = call_lambda(tr, x), lp = call_lambda(ctx, x))
  },
  function(z) same_value(z$k, z$kp, 1e-8) &&
    field_parity(z$l, z$lp, c("lambda", "logL", "gls_mean", "sig2"), 1e-7)
)

record_case(
  "alternate biological root remains explicit and root-sensitive",
  function() {
    tr <- make_tree(8L, "pectinate", heterogeneous = TRUE)
    rooted_a <- ape::root(tr, outgroup = tr$tip.label[[1L]],
                          resolve.root = FALSE)
    rooted_b <- ape::root(tr, outgroup = tr$tip.label[[2L]],
                          resolve.root = FALSE)
    x <- make_continuous(rooted_a)
    a <- call_k(rooted_a, x)
    b <- call_k(rooted_b, x[rooted_b$tip.label])
    list(a = a, b = b, tree_a = rooted_a, tree_b = rooted_b)
  },
  function(z) is.finite(z$a) && is.finite(z$b) &&
    !identical(z$tree_a$edge, z$tree_b$edge) &&
    isTRUE(abs(as.numeric(z$a) - as.numeric(z$b)) > 1e-8)
)

record_case(
  "protected structural field mutations are rejected before estimation",
  function() {
    fields <- c("edge", "edge.length", "tip.label", "Nnode")
    out <- vapply(fields, function(field) {
      ctx <- prepare_tree(normal)
      if (field == "edge") {
        ctx$tree$edge[1L, 2L] <- ctx$tree$edge[1L, 2L] + 1L
      } else if (field == "edge.length") {
        ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 0.1
      } else if (field == "tip.label") {
        ctx$tree$tip.label[[1L]] <- paste0(ctx$tree$tip.label[[1L]], "_changed")
      } else {
        ctx$tree$Nnode <- ctx$tree$Nnode + 1L
      }
      msg <- tryCatch({
        call_k(ctx, normal_x)
        "NO_ERROR"
      }, error = function(e) conditionMessage(e))
      grepl("prepare_tree|fingerprint|context", msg, ignore.case = TRUE)
    }, logical(1))
    out
  },
  function(z) all(z)
)

record_case(
  "V1 and unknown prepared schemas have clear rejection",
  function() {
    v1 <- normal_ctx
    v1$context_schema_version <- 1L
    unknown <- normal_ctx
    unknown$protected_snapshot_version <- 999L
    missing <- normal_ctx
    missing$context_schema_version <- NULL
    one <- tryCatch({ call_k(v1, normal_x); "NO_ERROR" },
                    error = function(e) conditionMessage(e))
    two <- tryCatch({ call_k(unknown, normal_x); "NO_ERROR" },
                    error = function(e) conditionMessage(e))
    three <- tryCatch({ call_k(missing, normal_x); "NO_ERROR" },
                      error = function(e) conditionMessage(e))
    c(v1 = one, unknown = two, missing = three)
  },
  function(z) all(grepl("prepare_tree", z, fixed = TRUE))
)

cat("\n-- persistence and replay --\n")
persist_rds <- file.path(out_dir, "scientific_context_v2.rds")
persist_result <- file.path(out_dir, "fresh_session_result.rds")
saveRDS(normal_ctx, persist_rds, version = 2L)

record_case(
  "V2 saveRDS/readRDS fresh-session K/lambda/D and mutation gate",
  function() {
    rscript <- file.path(R.home("bin"), "Rscript")
    if (.Platform$OS.type == "windows" && !file.exists(rscript)) {
      rscript <- paste0(rscript, ".exe")
    }
    child_script <- normalizePath(
      file.path(repo, "benchmarks", "final_release", "scientific",
                "run_scientific_validation.R"),
      winslash = "/", mustWork = TRUE
    )
    old_lib <- Sys.getenv("FASTPHYLOSIG_LIBRARY", unset = NA_character_)
    pkg_path <- getNamespaceInfo(asNamespace("fastphylosig"), "path")
    Sys.setenv(
      FASTPHYLOSIG_LIBRARY = dirname(pkg_path),
      R_LIBS_USER = paste(unique(c(dirname(pkg_path), .libPaths())),
                          collapse = .Platform$path.sep)
    )
    on.exit({
      if (is.na(old_lib)) Sys.unsetenv("FASTPHYLOSIG_LIBRARY") else
        Sys.setenv(FASTPHYLOSIG_LIBRARY = old_lib)
    }, add = TRUE)
    output <- system2(
      rscript,
      args = c("--vanilla", shQuote(child_script), "--child-rds",
               shQuote(persist_rds), "--child-out", shQuote(persist_result)),
      stdout = TRUE, stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    result <- if (file.exists(persist_result)) readRDS(persist_result) else NULL
    list(status = as.integer(status), output = output, result = result)
  },
  function(z) {
    !is.null(z$result) && z$status == 0L &&
      identical(z$result$context_schema_version, 2L) &&
      identical(z$result$canonical_contract_version, 2L) &&
      identical(z$result$protected_snapshot_version, 2L) &&
      is.finite(z$result$k) && is.finite(z$result$lambda) &&
      is.finite(z$result$lambda_logL) && is.finite(z$result$d) &&
      grepl("prepare_tree|fingerprint|context", z$result$mutation_message,
            ignore.case = TRUE)
  }
)

record_case(
  "same-seed replay K/D/Delta under identical configuration",
  function() {
    set.seed(2026091310L)
    k1 <- call_k(normal, normal_x, test = TRUE, nsim = 16L,
                 return_sim = TRUE)
    set.seed(2026091310L)
    k2 <- call_k(normal, normal_x, test = TRUE, nsim = 16L,
                 return_sim = TRUE)
    set.seed(2026091311L)
    d1 <- call_d(normal, normal_binary, test = TRUE, nsim = 8L,
                 return_sim = TRUE)
    set.seed(2026091311L)
    d2 <- call_d(normal, normal_binary, test = TRUE, nsim = 8L,
                 return_sim = TRUE)
    set.seed(2026091312L)
    de1 <- call_delta(delta_tree, delta_x, test = TRUE, nsim = 3L,
                      permutations = delta_perms, return_sim = TRUE)
    set.seed(2026091312L)
    de2 <- call_delta(delta_tree, delta_x, test = TRUE, nsim = 3L,
                      permutations = delta_perms, return_sim = TRUE)
    list(k1 = k1, k2 = k2, d1 = d1, d2 = d2, de1 = de1, de2 = de2)
  },
  function(z) field_parity(z$k1, z$k2,
                           c("K", "P", "MCSE_P", "exceedance_count", "sim.K"),
                           1e-12) &&
    field_parity(z$d1, z$d2,
                 c("DEstimate", "Pval1", "Pval0", "Parameters",
                   "Permutations"), 1e-12) &&
    field_parity(z$de1, z$de2,
                 c("delta", "alpha_mean", "beta_mean", "n_saved",
                   "ESS_alpha", "ESS_beta", "MCSE_Delta", "P", "P_MCSE",
                   "sim.delta", "status"), 1e-10)
  , warning_expected = TRUE
)

write.csv(records, file.path(out_dir, "status.csv"), row.names = FALSE,
          na = "")

env_lines <- c(
  paste0("commit=", source_commit),
  paste0("package_version=", as.character(packageVersion("fastphylosig"))),
  paste0("package_library=", Sys.getenv("FASTPHYLOSIG_LIBRARY", "not supplied")),
  paste0("loaded_package_path=", getNamespaceInfo(asNamespace("fastphylosig"), "path")),
  paste0("R=", R.version$version.string),
  paste0("platform=", R.version$platform),
  paste0("os_type=", .Platform$OS.type),
  paste0("ape=", as.character(packageVersion("ape"))),
  paste0("locale_collate=", Sys.getlocale("LC_COLLATE")),
  paste0("locale_ctype=", Sys.getlocale("LC_CTYPE")),
  paste0("openmp_capability=", if (isTRUE(capabilities("long.double")))
    "reported-by-R" else "not-reported"),
  paste0("working_directory=", normalizePath(getwd(), winslash = "/")),
  paste0("production_change_scope=none; cross_platform=not_run"),
  capture.output(sessionInfo())
)
writeLines(env_lines, file.path(out_dir, "environment.txt"), useBytes = TRUE)

summary <- data.frame(
  metric = c(
    "total_cases", "pass_cases", "fail_cases", "skip_cases",
    "unexpected_warning_cases", "scientific_validation"
  ),
  value = c(
    nrow(records), sum(records$status == "PASS"), sum(records$status == "FAIL"),
    sum(records$status == "SKIP"),
    sum(nzchar(records$warnings) & !records$warning_expected),
    if (all(records$status != "FAIL") &&
        sum(nzchar(records$warnings) & !records$warning_expected) == 0L) {
      "PASS"
    } else {
      "FAIL"
    }
  ), stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE,
          na = "")

cat("\nSummary:\n")
print(summary, row.names = FALSE)
cat("Completed:", format(Sys.time(), tz = "UTC"), "UTC\n")

unexpected_warning_cases <- sum(nzchar(records$warnings) &
                                !records$warning_expected)
if (any(records$status == "FAIL") || unexpected_warning_cases > 0L) {
  quit(save = "no", status = 1L)
}
