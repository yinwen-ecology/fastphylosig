# Shared tree requirements and issue metadata --------------------------------

# The registry is deliberately descriptive rather than executable.  Production
# kernels continue to own their numerical contracts; this table only provides a
# single vocabulary for the read-only preflight and for callers that need to
# explain why a method is (or is not) ready.
.tree_requirement_registry <- local({
  registry <- list(
    K = list(
      method = "K",
      min_tips = 2L,
      finite_branch_lengths = TRUE,
      zero_terminal = "error",
      zero_internal = "error",
      structural_root = TRUE,
      conventional_root = FALSE,
      polytomy = "supported",
      unary = "supported",
      binary = FALSE,
      non_ultrametric = "allowed"
    ),
    lambda = list(
      method = "lambda",
      min_tips = 2L,
      finite_branch_lengths = TRUE,
      zero_terminal = "user_action_required",
      zero_internal = "warning",
      structural_root = TRUE,
      conventional_root = FALSE,
      polytomy = "supported",
      unary = "supported",
      binary = FALSE,
      non_ultrametric = "allowed"
    ),
    D = list(
      method = "D",
      min_tips = 2L,
      finite_branch_lengths = TRUE,
      zero_terminal = "error",
      zero_internal = "error",
      structural_root = TRUE,
      conventional_root = FALSE,
      polytomy = "compatibility",
      unary = "error",
      binary = FALSE,
      non_ultrametric = "allowed"
    ),
    Delta = list(
      method = "Delta",
      min_tips = 2L,
      finite_branch_lengths = TRUE,
      zero_terminal = "error",
      zero_internal = "error",
      structural_root = TRUE,
      conventional_root = TRUE,
      polytomy = "error",
      unary = "error",
      binary = TRUE,
      non_ultrametric = "allowed"
    ),
    ACE = list(
      method = "ACE",
      min_tips = 2L,
      finite_branch_lengths = TRUE,
      zero_terminal = "error",
      zero_internal = "error",
      structural_root = TRUE,
      conventional_root = TRUE,
      polytomy = "error",
      unary = "error",
      binary = TRUE,
      non_ultrametric = "allowed"
    )
  )

  function() {
    # Return a fresh copy so a diagnostic consumer cannot mutate the package
    # level vocabulary for subsequent checks.
    lapply(registry, function(requirement) as.list(requirement))
  }
})

.tree_requirement <- function(signal) {
  if (missing(signal) || !is.character(signal) || length(signal) != 1L) {
    stop("signal must identify one registered method.", call. = FALSE)
  }
  key <- tolower(trimws(signal[[1L]]))
  key <- switch(key, k = "K", lambda = "lambda", d = "D",
                delta = "Delta", ace = "ACE", signal[[1L]])
  registry <- .tree_requirement_registry()
  if (!key %in% names(registry)) {
    stop("signal is not present in the tree requirement registry.", call. = FALSE)
  }
  registry[[key]]
}

# Keep this formatter in one place.  The fields and wording intentionally
# mirror the historical check_tree() output; moving the code here is a
# vocabulary refactor, not a change to the readiness or numerical rules.
.tree_issue_metadata <- function(code, severity, signal_name, message) {
  code <- as.character(code[[1L]])
  severity <- as.character(severity[[1L]])
  signal_name <- as.character(signal_name[[1L]])
  representation <- code %in% c(
    "d_requires_canonical_root", "delta_requires_canonical_root"
  )
  check <- if (representation) {
    "representation"
  } else if (grepl("tip_label", code, fixed = TRUE)) {
    "tip labels"
  } else if (code == "too_few_tips") {
    "tree size"
  } else if (grepl("root_distance|tip_height|ultrametric", code)) {
    "branch geometry"
  } else if (grepl("branch|zero", code)) {
    "branch lengths"
  } else if (grepl("root", code)) {
    "rooting"
  } else if (grepl("edge|nnode|topology|polytom|binary|single_child", code)) {
    "topology"
  } else {
    "tree input"
  }
  problem <- as.character(message[[1L]])
  action <- if (representation) {
    "canonicalize the internal representation with resolve_tree() and re-run check_tree()"
  } else if (grepl("tip_label", code, fixed = TRUE)) {
    "provide non-empty, non-missing, unique species names"
  } else if (code == "too_few_tips") {
    "provide a tree and matched data retaining at least two tips"
  } else if (grepl("requires_rooted|conventional_unrooted|conventional_root", code)) {
    "provide a biologically justified rooted phylogeny; specify an outgroup explicitly when appropriate"
  } else if (grepl("binary_tree|rooted_binary", code)) {
    "provide a biologically justified rooted binary phylogeny"
  } else if (grepl("single_child", code)) {
    "review unary nodes and collapse them only when that transformation is scientifically justified"
  } else if (grepl("branch|zero|tip_height|root_distance", code)) {
    "provide branch lengths that satisfy the selected method; do not replace zero or invalid values automatically"
  } else if (grepl("edge|nnode|topology", code)) {
    "provide a connected, acyclic phylo edge structure with consistent node metadata"
  } else if (severity %in% c("ERROR", "USER_ACTION_REQUIRED")) {
    "repair the reported tree issue, then re-run check_tree()"
  } else {
    "inspect this diagnostic with check_tree()"
  }
  list(check = check, problem = problem, action = action,
       auto_fixable = isTRUE(representation))
}
