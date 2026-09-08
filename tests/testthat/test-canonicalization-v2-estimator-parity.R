# Canonicalization V2 estimator parity ---------------------------------------

# This file exercises only the public raw/prepared boundary.  The equivalent
# trees keep public tip IDs and branch-child associations while varying the
# internal IDs and edge-row order.

.v2_ep_make_tree <- function() {
  tree <- list(
    edge = matrix(
      c(
        11L, 9L, 11L, 10L,
        9L, 7L, 9L, 8L,
        7L, 1L, 7L, 2L,
        8L, 3L, 8L, 4L,
        10L, 5L, 10L, 6L
      ),
      ncol = 2L,
      byrow = TRUE
    ),
    edge.length = c(0.8, 0.9, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 1.1, 1.2),
    Nnode = 5L,
    tip.label = c("a", "b", "c", "d", "e", "f"),
    root.edge = 0.125,
    v2_metadata = list(source = "estimator-parity")
  )
  class(tree) <- "phylo"
  tree
}

.v2_ep_renumber_internals <- function(tree) {
  id_map <- stats::setNames(
    c(seq_len(6L), 17L, 13L, 21L, 19L, 25L),
    as.character(seq_len(11L))
  )
  out <- tree
  mapped <- unname(id_map[as.character(tree$edge)])
  if (anyNA(mapped)) stop("test fixture contains an unmapped node",
                          call. = FALSE)
  out$edge <- matrix(
    as.integer(mapped), nrow = nrow(tree$edge), ncol = ncol(tree$edge),
    dimnames = dimnames(tree$edge)
  )
  out
}

.v2_ep_shuffle_edges <- function(tree) {
  rows <- c(7L, 2L, 10L, 4L, 1L, 9L, 5L, 8L, 3L, 6L)
  out <- tree
  out$edge <- tree$edge[rows, , drop = FALSE]
  out$edge.length <- tree$edge.length[rows]
  out
}

.v2_ep_variants <- function() {
  tree <- .v2_ep_make_tree()
  list(
    original = tree,
    internal_renumbered = .v2_ep_renumber_internals(tree),
    edge_shuffled = .v2_ep_shuffle_edges(tree),
    renumbered_and_shuffled = .v2_ep_shuffle_edges(
      .v2_ep_renumber_internals(tree)
    )
  )
}

.v2_ep_continuous <- function(tree) {
  stats::setNames(c(0.3, 1.7, -0.8, 0.4, 2.2, -1.1), tree$tip.label)
}

.v2_ep_binary <- function(tree) {
  stats::setNames(c(0, 1, 0, 1, 1, 0), tree$tip.label)
}

.v2_ep_categorical <- function(tree) {
  stats::setNames(
    factor(c("A", "B", "C", "A", "B", "C"),
           levels = c("A", "B", "C")),
    tree$tip.label
  )
}

.v2_ep_permutations <- function() {
  rbind(
    c(1L, 2L, 3L, 4L, 5L, 6L),
    c(6L, 5L, 4L, 3L, 2L, 1L),
    c(2L, 1L, 4L, 3L, 6L, 5L),
    c(3L, 4L, 5L, 6L, 1L, 2L)
  )
}

.v2_ep_capture <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(condition) {
      warnings[[length(warnings) + 1L]] <<- list(
        class = class(condition), message = conditionMessage(condition)
      )
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

.v2_ep_with_seed <- function(seed, expr) {
  set.seed(seed)
  .v2_ep_capture(expr)
}

.v2_ep_warning_classes <- function(captured) {
  lapply(captured$warnings, `[[`, "class")
}

.v2_ep_check_warnings <- function(raw, prepared, info) {
  expect_identical(
    .v2_ep_warning_classes(raw), .v2_ep_warning_classes(prepared),
    info = paste(info, "warning classes")
  )
}

.v2_ep_check_fields <- function(raw, prepared, fields, tolerance = 0,
                                info = "parity") {
  expect_true(all(fields %in% names(raw)),
              info = paste(info, "raw fields"))
  expect_true(all(fields %in% names(prepared)),
              info = paste(info, "prepared fields"))
  for (field in fields) {
    expect_equal(
      raw[[field]], prepared[[field]], tolerance = tolerance,
      info = paste(info, field)
    )
  }
  expect_identical(raw$status, prepared$status,
                   info = paste(info, "status"))
  if ("message" %in% names(raw) && "message" %in% names(prepared)) {
    expect_identical(raw$message, prepared$message,
                     info = paste(info, "message"))
  }
}

.v2_ep_retained_signature <- function(result) {
  report <- attr(result, "match_report", exact = TRUE)
  if (is.list(report)) {
    return(report[c("matched_species", "removed_tree_tips",
                    "removed_data_rows")])
  }
  metadata <- attr(result, "analysis_metadata", exact = TRUE)
  if (is.list(metadata)) {
    return(metadata[c("matched_species", "tree_tips_removed",
                      "data_rows_removed")])
  }
  NULL
}

.v2_ep_check_retained <- function(raw, prepared, info) {
  raw_retained <- .v2_ep_retained_signature(raw)
  prepared_retained <- .v2_ep_retained_signature(prepared)
  expect_false(is.null(raw_retained), info = paste(info, "raw retained data"))
  expect_false(is.null(prepared_retained),
               info = paste(info, "prepared retained data"))
  expect_identical(raw_retained, prepared_retained,
                   info = paste(info, "retained species"))
}

.v2_ep_source_snapshot <- function(value) {
  serialize(value, connection = NULL, version = 3L)
}

.v2_ep_ace_by_subtree <- function(fit, tree) {
  children <- split(as.integer(tree$edge[, 2L]),
                    as.character(tree$edge[, 1L]))
  descendant_tips <- function(node) {
    if (node <= length(tree$tip.label)) return(tree$tip.label[[node]])
    unlist(lapply(children[[as.character(node)]], descendant_tips),
           use.names = FALSE)
  }
  source_ids <- as.integer(rownames(fit$lik.anc))
  keys <- vapply(source_ids, function(node) {
    paste(sort(descendant_tips(node), method = "radix"), collapse = "\n")
  }, character(1L))
  out <- fit$lik.anc
  rownames(out) <- keys
  out[order(keys, method = "radix"), , drop = FALSE]
}

test_that("V2 raw and prepared K/lambda estimators are representation invariant", {
  testthat::skip_if_not_installed("ape")
  variants <- .v2_ep_variants()
  trait <- .v2_ep_continuous(variants[[1L]])
  permutations <- .v2_ep_permutations()
  reference <- list(k = NULL, lambda = NULL)

  for (variant_name in names(variants)) {
    tree <- variants[[variant_name]]
    tree_before <- .v2_ep_source_snapshot(tree)
    trait_before <- .v2_ep_source_snapshot(trait)
    context <- prepare_tree(tree)
    expect_identical(.v2_ep_source_snapshot(tree), tree_before,
                     info = paste(variant_name, "prepare_tree immutability"))

    raw_k <- .v2_ep_capture(fast_k(
      tree, trait, test = TRUE, nsim = nrow(permutations),
      permutations = permutations, return_sim = TRUE, keep_null = TRUE,
      verbose = FALSE, progress = FALSE
    ))
    prepared_k <- .v2_ep_capture(fast_k(
      context, trait, test = TRUE, nsim = nrow(permutations),
      permutations = permutations, return_sim = TRUE, keep_null = TRUE,
      verbose = FALSE, progress = FALSE
    ))
    .v2_ep_check_fields(
      raw_k$value, prepared_k$value,
      c("K", "P", "MCSE_P", "nsim_requested", "nsim_successful",
        "exceedance_count", "sim.K"),
      tolerance = 1e-12, info = paste(variant_name, "K")
    )
    .v2_ep_check_warnings(raw_k, prepared_k, paste(variant_name, "K"))
    .v2_ep_check_retained(raw_k$value, prepared_k$value,
                          paste(variant_name, "K"))
    expect_identical(.v2_ep_source_snapshot(tree), tree_before,
                     info = paste(variant_name, "K tree immutability"))
    expect_identical(.v2_ep_source_snapshot(trait), trait_before,
                     info = paste(variant_name, "K trait immutability"))
    if (is.null(reference$k)) reference$k <- raw_k$value
    .v2_ep_check_fields(
      reference$k, raw_k$value,
      c("K", "P", "MCSE_P", "nsim_requested", "nsim_successful",
        "exceedance_count", "sim.K"),
      tolerance = 1e-12, info = paste(variant_name, "K representation")
    )

    raw_lambda <- .v2_ep_capture(fast_lambda(
      tree, trait, test = TRUE, lambda_profile = FALSE,
      verbose = FALSE, progress = FALSE
    ))
    prepared_lambda <- .v2_ep_capture(fast_lambda(
      context, trait, test = TRUE, lambda_profile = FALSE,
      verbose = FALSE, progress = FALSE
    ))
    .v2_ep_check_fields(
      raw_lambda$value, prepared_lambda$value,
      c("lambda", "logL", "gls_mean", "sig2", "logL0", "LR", "P"),
      tolerance = 1e-12, info = paste(variant_name, "lambda")
    )
    .v2_ep_check_warnings(raw_lambda, prepared_lambda,
                          paste(variant_name, "lambda"))
    .v2_ep_check_retained(raw_lambda$value, prepared_lambda$value,
                          paste(variant_name, "lambda"))
    expect_identical(.v2_ep_source_snapshot(tree), tree_before,
                     info = paste(variant_name, "lambda tree immutability"))
    expect_identical(.v2_ep_source_snapshot(trait), trait_before,
                     info = paste(variant_name, "lambda trait immutability"))
    if (is.null(reference$lambda)) reference$lambda <- raw_lambda$value
    .v2_ep_check_fields(
      reference$lambda, raw_lambda$value,
      c("lambda", "logL", "gls_mean", "sig2", "logL0", "LR", "P"),
      tolerance = 1e-12, info = paste(variant_name, "lambda representation")
    )
  }
})

test_that("V2 raw and prepared D estimators preserve controlled nulls", {
  testthat::skip_if_not_installed("ape")
  variants <- .v2_ep_variants()
  trait <- .v2_ep_binary(variants[[1L]])
  random_states <- matrix(
    c(0, 1, 0, 1, 0, 1,
      1, 0, 1, 0, 1, 0,
      0, 0, 1, 1, 0, 1),
    nrow = 6L, ncol = 3L
  )
  brownian_states <- matrix(
    c(0, 0, 1, 1, 0, 0,
      1, 1, 0, 0, 1, 1,
      0, 1, 1, 0, 0, 1),
    nrow = 6L, ncol = 3L
  )
  reference <- NULL

  for (variant_name in names(variants)) {
    tree <- variants[[variant_name]]
    tree_before <- .v2_ep_source_snapshot(tree)
    trait_before <- .v2_ep_source_snapshot(trait)
    context <- prepare_tree(tree)
    raw <- .v2_ep_capture(fast_d(
      tree, trait, test = TRUE, nsim = 3L,
      random_states = random_states, brownian_states = brownian_states,
      return_sim = TRUE, verbose = FALSE, progress = FALSE
    ))
    prepared <- .v2_ep_capture(fast_d(
      context, trait, test = TRUE, nsim = 3L,
      random_states = random_states, brownian_states = brownian_states,
      return_sim = TRUE, verbose = FALSE, progress = FALSE
    ))
    .v2_ep_check_fields(
      raw$value, prepared$value,
      c("DEstimate", "Pval1", "Pval0", "P_random", "P_Brownian",
        "MCSE_P_random", "MCSE_P_Brownian", "Parameters", "StatesTable",
        "nPermut", "nsim_requested", "nsim_successful_random",
        "nsim_successful_brownian", "Permutations"),
      tolerance = 1e-12, info = paste(variant_name, "D")
    )
    .v2_ep_check_warnings(raw, prepared, paste(variant_name, "D"))
    .v2_ep_check_retained(raw$value, prepared$value,
                          paste(variant_name, "D"))
    expect_identical(.v2_ep_source_snapshot(tree), tree_before,
                     info = paste(variant_name, "D tree immutability"))
    expect_identical(.v2_ep_source_snapshot(trait), trait_before,
                     info = paste(variant_name, "D trait immutability"))
    if (is.null(reference)) reference <- raw$value
    .v2_ep_check_fields(
      reference, raw$value,
      c("DEstimate", "Pval1", "Pval0", "P_random", "P_Brownian",
        "MCSE_P_random", "MCSE_P_Brownian", "Parameters", "StatesTable",
        "nPermut", "nsim_requested", "nsim_successful_random",
        "nsim_successful_brownian", "Permutations"),
      tolerance = 1e-12, info = paste(variant_name, "D representation")
    )
  }
})

test_that("V2 raw and prepared Delta estimators preserve fixed stochastic paths", {
  testthat::skip_if_not_installed("ape")
  variants <- .v2_ep_variants()
  trait <- .v2_ep_categorical(variants[[1L]])
  permutations <- .v2_ep_permutations()
  delta_args <- list(
    x = trait, test = TRUE, nsim = nrow(permutations),
    mcmc_sim = 24L, thin = 2L, burn = 4L, lambda0 = 0.1,
    proposal_sd = 0.5, model = "ARD", permutations = permutations,
    return_sim = TRUE, verbose = FALSE, progress = FALSE, ncores = 1L
  )
  fields <- c(
    "delta", "alpha_mean", "beta_mean", "n_saved", "n_saved_requested",
    "n_saved_successful", "alpha_sd", "beta_sd", "ESS_alpha", "ESS_beta",
    "split_Rhat_alpha", "split_Rhat_beta", "Rhat_alpha", "Rhat_beta",
    "alpha_beta_cov", "MCSE_Delta", "diagnostics_available",
    "diagnostics_note", "requested_simulations", "successful_simulations",
    "requested_iterations", "successful_iterations",
    "requested_permutations", "successful_permutations", "n_states",
    "n_species", "status", "note", "message", "parameters", "P",
    "P_MCSE", "P_mcse", "MCSE_P", "n_failed_sim", "sim.delta"
  )
  reference <- NULL

  for (variant_name in names(variants)) {
    tree <- variants[[variant_name]]
    tree_before <- .v2_ep_source_snapshot(tree)
    trait_before <- .v2_ep_source_snapshot(trait)
    context <- prepare_tree(tree)
    raw <- .v2_ep_with_seed(20260908L, do.call(
      fast_delta, c(list(tree = tree), delta_args)
    ))
    prepared <- .v2_ep_with_seed(20260908L, do.call(
      fast_delta, c(list(tree = context), delta_args)
    ))
    .v2_ep_check_fields(raw$value, prepared$value, fields,
                        tolerance = 1e-12, info = paste(variant_name, "Delta"))
    .v2_ep_check_warnings(raw, prepared, paste(variant_name, "Delta"))
    .v2_ep_check_retained(raw$value, prepared$value,
                          paste(variant_name, "Delta"))
    expect_identical(.v2_ep_source_snapshot(tree), tree_before,
                     info = paste(variant_name, "Delta tree immutability"))
    expect_identical(.v2_ep_source_snapshot(trait), trait_before,
                     info = paste(variant_name, "Delta trait immutability"))
    if (is.null(reference)) reference <- raw$value
    .v2_ep_check_fields(reference, raw$value, fields,
                        tolerance = 1e-12,
                        info = paste(variant_name, "Delta representation"))
  }
})

test_that("V2 raw and prepared ACE estimators preserve likelihood outputs", {
  testthat::skip_if_not_installed("ape")
  variants <- .v2_ep_variants()
  trait <- .v2_ep_categorical(variants[[1L]])
  fields <- c("loglik", "rates", "se", "index.matrix", "convergence",
              "message", "lik.anc")
  reference <- NULL

  for (variant_name in names(variants)) {
    tree <- variants[[variant_name]]
    tree_before <- .v2_ep_source_snapshot(tree)
    trait_before <- .v2_ep_source_snapshot(trait)
    context <- prepare_tree(tree)
    raw <- .v2_ep_capture(fast_ace(
      trait, phy = tree, model = "ER", CI = TRUE, marginal = FALSE,
      progress = FALSE
    ))
    prepared <- .v2_ep_capture(fast_ace(
      trait, prepared = context, model = "ER", CI = TRUE, marginal = FALSE,
      progress = FALSE
    ))
    .v2_ep_check_fields(raw$value, prepared$value, fields,
                        tolerance = 1e-12, info = paste(variant_name, "ACE"))
    .v2_ep_check_warnings(raw, prepared, paste(variant_name, "ACE"))
    expect_identical(.v2_ep_source_snapshot(tree), tree_before,
                     info = paste(variant_name, "ACE tree immutability"))
    expect_identical(.v2_ep_source_snapshot(trait), trait_before,
                     info = paste(variant_name, "ACE trait immutability"))
    if (is.null(reference)) {
      reference <- list(
        value = raw$value,
        lik.anc = .v2_ep_ace_by_subtree(raw$value, tree)
      )
    }
    .v2_ep_check_fields(reference$value, raw$value,
                        setdiff(fields, "lik.anc"),
                        tolerance = 1e-12,
                        info = paste(variant_name, "ACE representation"))
    expect_equal(
      reference$lik.anc,
      .v2_ep_ace_by_subtree(raw$value, tree),
      tolerance = 1e-12,
      info = paste(variant_name, "ACE representation lik.anc by subtree")
    )
  }
})
