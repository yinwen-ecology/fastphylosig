# Trait-level failure semantics audit ---------------------------------------

.audit_failure_tree <- function(n = 8L) {
  set.seed(20260811L + as.integer(n))
  ape::rtree(n)
}

.audit_named <- function(x, tree) {
  stats::setNames(x, tree$tip.label)
}

.audit_try <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = force(expr), error = NULL),
    error = function(e) list(ok = FALSE, value = NULL,
                              error = conditionMessage(e))
  )
}

.audit_has_diagnostic <- function(x) {
  is.data.frame(x) && any(c("status", "message", "note",
                            "diagnostics_note") %in% names(x))
}

.audit_status_levels <- c(
  "ok", "partial", "invalid_trait", "insufficient_data", "undefined",
  "failed"
)

test_that("constant and near-constant continuous traits remain batch-auditable", {
  tree <- .audit_failure_tree()
  n <- ape::Ntip(tree)
  X <- cbind(
    constant = rep(4, n),
    near_constant = 4 + seq(-1, 1, length.out = n) * 1e-9,
    finite = seq_len(n),
    nan_missing = c(seq_len(n - 1L), NaN),
    too_few = c(1, rep(NA_real_, n - 1L))
  )
  dimnames(X) <- list(tree$tip.label, colnames(X))

  k <- suppressMessages(.audit_try(fast_k(
    tree, X, test = FALSE, verbose = FALSE, progress = FALSE
  )))
  expect_true(k$ok)
  expect_true(.audit_has_diagnostic(k$value))
  expect_true(all(k$value$status %in% .audit_status_levels))
  expect_true(any(grepl("constant|degenerate", k$value$note,
                        ignore.case = TRUE), na.rm = TRUE))
  expect_true(is.na(k$value$K_fast[[1L]]) ||
              !is.finite(k$value$K_fast[[1L]]))
  expect_true(k$value$n_removed_na[[4L]] >= 1L)
  expect_equal(k$value$status[[5L]], "insufficient_data")
  expect_true(grepl("fewer than 2", k$value$note[[5L]], fixed = TRUE))
  expect_true(k$value$status[[1L]] %in% c("undefined", "invalid_trait"))
  expect_true(is.character(k$value$message[[1L]]) &&
              nzchar(k$value$message[[1L]]))

  lambda <- suppressMessages(.audit_try(fast_lambda(
    tree, X, verbose = FALSE, progress = FALSE,
    lambda_profile = FALSE
  )))
  expect_true(lambda$ok)
  expect_true(.audit_has_diagnostic(lambda$value))
  expect_true(all(lambda$value$status %in% .audit_status_levels))
  # A near-constant fit may be recovered by the dense fallback; the audit is
  # interested in preserving either a finite estimate or an explicit note.
  near_ok <- is.finite(lambda$value$lambda_fast[[2L]]) ||
    (length(lambda$value$note[[2L]]) &&
       !is.na(lambda$value$note[[2L]]) && nzchar(lambda$value$note[[2L]]))
  expect_true(near_ok)
})

test_that("continuous non-finite inputs are per-trait invalid_trait rows", {
  tree <- .audit_failure_tree()
  n <- ape::Ntip(tree)
  X <- cbind(
    good = seq_len(n),
    bad_inf = c(seq_len(n - 1L), Inf),
    bad_nan = c(seq_len(n - 1L), NaN),
    bad_neg_inf = c(seq_len(n - 1L), -Inf)
  )
  dimnames(X) <- list(tree$tip.label, colnames(X))

  k <- suppressMessages(fast_k(
    tree, X, verbose = FALSE, progress = FALSE
  ))
  lambda <- suppressMessages(fast_lambda(
    tree, X, verbose = FALSE, progress = FALSE,
    lambda_profile = FALSE
  ))
  expect_true(all(c("status", "message", "note") %in% names(k)))
  expect_true(all(c("status", "message", "note") %in% names(lambda)))
  expect_true(all(k$status %in% .audit_status_levels))
  expect_true(all(lambda$status %in% .audit_status_levels))
  expect_equal(k$status[2:4], rep("invalid_trait", 3L))
  expect_equal(lambda$status[2:4], rep("invalid_trait", 3L))
  expect_true(all(!is.finite(k$K_fast[2:4])))
  expect_true(all(!is.finite(lambda$lambda_fast[2:4])))
  expect_true(all(grepl("non-finite", k$message[2:4], fixed = TRUE)))
  expect_true(all(grepl("non-finite", lambda$message[2:4], fixed = TRUE)))

  # The no-test K scalar remains atomic; status and message are attributes.
  k_scalar <- suppressMessages(fast_k(
    tree, .audit_named(c(seq_len(n - 1L), Inf), tree),
    verbose = FALSE, progress = FALSE
  ))
  expect_true(is.numeric(k_scalar) && !is.list(k_scalar))
  expect_identical(attr(k_scalar, "status"), "invalid_trait")
  expect_true(nzchar(attr(k_scalar, "message")))

  k_scalar_test <- suppressMessages(fast_k(
    tree, .audit_named(c(seq_len(n - 1L), Inf), tree), test = TRUE,
    nsim = 2, verbose = FALSE, progress = FALSE
  ))
  expect_true(is.list(k_scalar_test))
  expect_identical(k_scalar_test$status, "invalid_trait")
  expect_true(nzchar(k_scalar_test$message))

  lambda_scalar <- suppressMessages(fast_lambda(
    tree, .audit_named(c(seq_len(n - 1L), Inf), tree),
    verbose = FALSE, progress = FALSE, lambda_profile = FALSE
  ))
  expect_true(is.list(lambda_scalar))
  expect_identical(lambda_scalar$status, "invalid_trait")
  expect_true(nzchar(lambda_scalar$message))
  expect_true(is.character(lambda_scalar$note) && nzchar(lambda_scalar$note))
})

test_that("binary single-state, non-finite, and too-few traits retain notes", {
  tree <- .audit_failure_tree()
  n <- ape::Ntip(tree)
  X <- cbind(
    single_state = rep(0, n),
    valid = rep(0:1, length.out = n),
    too_few = c(0, rep(NA_real_, n - 1L)),
    infinite = c(rep(0:1, length.out = n - 1L), Inf)
  )
  dimnames(X) <- list(tree$tip.label, colnames(X))

  fit <- suppressMessages(fast_d(
    tree, X, test = FALSE, nsim = 4, return_sim = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  expect_true(is.data.frame(fit))
  expect_true(.audit_has_diagnostic(fit))
  expect_true(all(c("status", "message", "note") %in% names(fit)))
  expect_true(all(fit$status %in% .audit_status_levels))
  expect_true(grepl("single state", fit$note[[1L]], fixed = TRUE))
  expect_true(grepl("fewer than 2", fit$note[[3L]], fixed = TRUE))
  expect_true(grepl("non-finite|infinite", fit$note[[4L]],
                    ignore.case = TRUE))
  expect_equal(fit$status[[1L]], "invalid_trait")
  expect_equal(fit$status[[3L]], "insufficient_data")
  expect_equal(fit$status[[4L]], "invalid_trait")
  expect_true(is.na(fit$D_fast[[1L]]))

  scalar <- suppressMessages(fast_d(
    tree, X[, "single_state"], test = FALSE, nsim = 4,
    return_sim = FALSE, verbose = FALSE, progress = FALSE
  ))
  expect_identical(scalar$status, "invalid_trait")
  expect_true(nzchar(scalar$message))
  expect_true(is.character(scalar$note) && nzchar(scalar$note))
})

test_that("categorical single-state and too-few traits are explicit in batches", {
  tree <- .audit_failure_tree()
  n <- ape::Ntip(tree)
  X <- data.frame(
    single_state = rep("a", n),
    valid = rep(c("a", "b", "c"), length.out = n),
    too_few = c("a", rep(NA_character_, n - 1L)),
    numeric_inf = c(seq_len(n - 1L), Inf),
    numeric_nan = c(seq_len(n - 1L), NaN),
    stringsAsFactors = FALSE,
    row.names = tree$tip.label
  )

  fit <- suppressWarnings(suppressMessages(fast_delta(
    tree, X, test = FALSE, mcmc_sim = 20, thin = 5, burn = 5,
    verbose = FALSE, progress = FALSE
  )))
  expect_true(is.data.frame(fit))
  expect_true(.audit_has_diagnostic(fit))
  expect_true(all(c("status", "message", "note") %in% names(fit)))
  expect_true(all(fit$status %in% .audit_status_levels))
  expect_true(grepl("single state", fit$note[[1L]], fixed = TRUE))
  expect_true(grepl("fewer than 2", fit$note[[3L]], fixed = TRUE))
  expect_equal(fit$status[[1L]], "invalid_trait")
  expect_equal(fit$status[[3L]], "insufficient_data")
  expect_equal(fit$status[[4L]], "invalid_trait")
  expect_equal(fit$status[[5L]], "invalid_trait")
  expect_true(grepl("non-finite", fit$note[[4L]], fixed = TRUE))
  expect_true(grepl("non-finite", fit$note[[5L]], fixed = TRUE))
  expect_true(is.na(fit$Delta_fast[[1L]]))
  expect_true(is.na(fit$Delta_fast[[3L]]))

  scalar <- suppressWarnings(suppressMessages(fast_delta(
    tree, X[, "single_state"], test = FALSE,
    mcmc_sim = 20, thin = 5, burn = 5,
    verbose = FALSE, progress = FALSE
  )))
  expect_identical(scalar$status, "invalid_trait")
  expect_true(nzchar(scalar$message))
  expect_true(is.character(scalar$note) && nzchar(scalar$note))
})
