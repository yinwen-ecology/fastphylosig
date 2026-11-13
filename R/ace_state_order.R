# Use observed tip partitions, not arbitrary state labels, for numerical order.
# This removes label-dependent optimizer coordinates; it is not a global-ML guarantee.
.fast_ace_state_order <- function(tip_state, workspace) {
  identity <- seq_len(workspace$nl)
  if (!identical(workspace$model, "ARD"))
    return(identity)
  # Use stable tip names so reordered tree representations use the same groups.
  tip_order <- order(workspace$tip_label, method = "radix")
  seen <- unique(tip_state[tip_order])
  seen <- seen[!is.na(seen)]
  c(seen, setdiff(identity, seen))
}

.fast_ace_remap_parameters <- function(x, inverse) {
  if (!is.list(x)) return(x)
  for (name in names(x)) {
    value <- x[[name]]
    if (name %in% c("par", "start") && is.numeric(value) &&
        length(value) == length(inverse)) {
      x[[name]] <- value[inverse]
    } else if (identical(name, "starts") && is.list(value)) {
      x[[name]] <- lapply(value, function(z)
        if (is.numeric(z) && length(z) == length(inverse)) z[inverse] else z)
    } else if (is.list(value) && !is.data.frame(value)) {
      x[[name]] <- .fast_ace_remap_parameters(value, inverse)
    }
  }
  x
}

.fast_ace_fit_workspace <- function(tip_state, lvls, workspace, ip,
                                    CI = TRUE, marginal = FALSE,
                                    estimate_se = TRUE,
                                    reuse_buffers = TRUE,
                                    fit_strategy = "default", recover_ard = FALSE) {
  # Existing minimal workspaces need no state metadata unless ARD reordering applies.
  state_order <- if (identical(workspace$model, "ARD"))
    .fast_ace_state_order(tip_state, workspace) else seq_along(lvls)
  if (identical(as.integer(state_order), seq_along(lvls))) {
    return(.fast_ace_fit_workspace_impl(tip_state, lvls, workspace, ip,
      CI, marginal, estimate_se, reuse_buffers, fit_strategy, recover_ard))
  }
  # ARD uses column-major off-diagonal parameter indices. Keep the external
  # matrix contract unchanged while permuting internal states and starts.
  rate_index <- workspace$rate$rate_index
  canonical_to_original <- rate_index[state_order, state_order][rate_index > 0L]
  inverse <- order(canonical_to_original)
  canonical_ip <- rep(ip, length.out = workspace$rate$np)[canonical_to_original]
  result <- withCallingHandlers(
    .fast_ace_fit_workspace_impl(match(tip_state, state_order), lvls[state_order],
      workspace, canonical_ip, CI, marginal, estimate_se, reuse_buffers,
      fit_strategy, recover_ard),
    warning = function(w) {
      if (!is.null(w$diagnostics)) {
        w$diagnostics <- .fast_ace_remap_parameters(w$diagnostics, inverse)
        warning(w)
        invokeRestart("muffleWarning")
      }
    })
  result$rates <- result$rates[inverse]
  result$se <- result$se[inverse]
  if (!is.null(result$lik.anc))
    result$lik.anc <- result$lik.anc[, order(state_order), drop = FALSE]
  if (!is.null(result$diagnostics))
    result$diagnostics <- .fast_ace_remap_parameters(result$diagnostics, inverse)
  result
}
