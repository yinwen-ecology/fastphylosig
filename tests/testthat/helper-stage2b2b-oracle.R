# Stage 2B2B Candidate 1: baseline inspection oracle and fixture helpers.

.stage2b2b_oracle_env <- new.env(parent = asNamespace("fastphylosig"))
.stage2b2b_oracle_file <- testthat::test_path(
  "fixtures", "stage2b2b-old-inspect-tree.R"
)
sys.source(.stage2b2b_oracle_file, envir = .stage2b2b_oracle_env)
.stage2b2b_old_inspect_tree_core <- get(
  ".stage2b2b_old_inspect_tree_core", envir = .stage2b2b_oracle_env
)
.stage2b2b_current_inspect_tree_core <- get(
  ".inspect_tree_core", envir = asNamespace("fastphylosig")
)

.stage2b2b_capture <- function(fun, tree, signal) {
  warning_records <- list()
  captured_error <- NULL
  value <- tryCatch(
    withCallingHandlers(
      fun(tree, signal = signal),
      warning = function(condition) {
        warning_records[[length(warning_records) + 1L]] <<- list(
          message = conditionMessage(condition),
          class = class(condition)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      captured_error <<- list(
        message = conditionMessage(condition),
        class = class(condition)
      )
      NULL
    }
  )
  list(value = value, warnings = warning_records, error = captured_error)
}

.stage2b2b_compare_capture <- function(old, current, label) {
  testthat::expect_identical(
    current$error, old$error,
    info = paste(label, "error")
  )
  testthat::expect_identical(
    current$warnings, old$warnings,
    info = paste(label, "warnings")
  )
  if (is.null(old$error) && is.null(current$error)) {
    testthat::expect_identical(
      current$value, old$value,
      info = paste(label, "complete inspection object")
    )
    for (field in c("ready", "ready_by_signal", "tree_summary", "issues")) {
      testthat::expect_identical(
        current$value[[field]], old$value[[field]],
        info = paste(label, field)
      )
    }
  }
  invisible(NULL)
}

.stage2b2b_clone <- function(tree) {
  unserialize(serialize(tree, connection = NULL, version = 2L))
}

.stage2b2b_make_pectinate <- function(n) {
  stopifnot(length(n) == 1L, n >= 2L)
  n <- as.integer(n)
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      parent <- n + i
      edge[2L * i - 1L, ] <- c(parent, i)
      edge[2L * i, ] <- c(parent, parent + 1L)
    }
  }
  root <- 2L * n - 1L
  edge[2L * n - 3L, ] <- c(root, n - 1L)
  edge[2L * n - 2L, ] <- c(root, n)
  structure(
    list(
      edge = edge,
      edge.length = rep(1, nrow(edge)),
      tip.label = paste0("tip", seq_len(n)),
      Nnode = n - 1L
    ),
    class = "phylo"
  )
}

.stage2b2b_shuffle_edges <- function(tree) {
  out <- .stage2b2b_clone(tree)
  rows <- rev(seq_len(nrow(out$edge)))
  out$edge <- out$edge[rows, , drop = FALSE]
  if (!is.null(out$edge.length)) out$edge.length <- out$edge.length[rows]
  out
}

.stage2b2b_renumber_internals <- function(tree, permutation) {
  out <- .stage2b2b_clone(tree)
  n_tip <- length(out$tip.label)
  old <- seq.int(n_tip + 1L, n_tip + out$Nnode)
  stopifnot(length(permutation) == length(old),
            setequal(permutation, old))
  mapping <- stats::setNames(as.integer(permutation), as.character(old))
  internal <- out$edge > n_tip
  out$edge[internal] <- unname(mapping[as.character(out$edge[internal])])
  if (!is.null(out$node.label)) {
    labels <- character(length(old))
    labels[as.integer(mapping) - n_tip] <- out$node.label
    out$node.label <- labels
  }
  out
}

.stage2b2b_unary_tree <- function() {
  structure(
    list(
      edge = matrix(c(3L, 4L, 4L, 1L, 3L, 2L), ncol = 2L, byrow = TRUE),
      edge.length = c(1, 1, 1),
      tip.label = c("a", "b"),
      Nnode = 2L
    ),
    class = "phylo"
  )
}

.stage2b2b_disconnected_tree <- function() {
  structure(
    list(
      edge = matrix(c(5L, 1L, 5L, 2L, 6L, 3L, 6L, 4L),
                    ncol = 2L, byrow = TRUE),
      edge.length = rep(1, 4L),
      tip.label = paste0("tip", 1:4),
      Nnode = 2L
    ),
    class = "phylo"
  )
}

.stage2b2b_fixture_set <- function() {
  random <- ape::rtree(14L)
  random$node.label <- paste0("node", seq_len(random$Nnode))
  balanced <- ape::compute.brlen(ape::stree(16L, type = "balanced"))
  pectinate <- .stage2b2b_make_pectinate(40L)
  polytomy <- ape::compute.brlen(ape::stree(12L, type = "star"))
  shuffled <- .stage2b2b_shuffle_edges(random)
  internal <- seq.int(length(random$tip.label) + 1L,
                      length(random$tip.label) + random$Nnode)
  renumbered <- .stage2b2b_renumber_internals(random, rev(internal))
  root <- setdiff(unique(random$edge[, 1L]), unique(random$edge[, 2L]))
  other <- setdiff(internal, root)[[1L]]
  root_swap <- internal
  root_swap[match(root, internal)] <- other
  root_swap[match(other, internal)] <- root
  different_root_numbering <- .stage2b2b_renumber_internals(random, root_swap)
  two_tip <- ape::read.tree(text = "(a:1,b:2);")
  large <- .stage2b2b_make_pectinate(1000L)

  cycle <- .stage2b2b_clone(random)
  cycle$edge[1L, 2L] <- cycle$edge[1L, 1L]
  disconnected <- .stage2b2b_disconnected_tree()
  unary <- .stage2b2b_unary_tree()
  malformed_edge <- .stage2b2b_clone(random)
  malformed_edge$edge <- matrix(seq_len(9L), nrow = 3L, ncol = 3L)
  invalid_nnode <- .stage2b2b_clone(random)
  invalid_nnode$Nnode <- NA_real_
  negative_branch <- .stage2b2b_clone(random)
  negative_branch$edge.length[[1L]] <- -1
  zero_terminal <- .stage2b2b_clone(random)
  terminal <- which(zero_terminal$edge[, 2L] <= length(zero_terminal$tip.label))[[1L]]
  zero_terminal$edge.length[[terminal]] <- 0
  invalid_root <- ape::unroot(.stage2b2b_clone(random))

  list(
    balanced = balanced,
    random = random,
    pectinate = pectinate,
    polytomy = polytomy,
    shuffled_edge_order = shuffled,
    safe_internal_renumbering = renumbered,
    different_root_numbering = different_root_numbering,
    two_tip = two_tip,
    large_pectinate = large,
    cycle = cycle,
    disconnected = disconnected,
    unary = unary,
    malformed_edge = malformed_edge,
    invalid_nnode = invalid_nnode,
    negative_branch = negative_branch,
    zero_terminal_branch = zero_terminal,
    invalid_root_representation = invalid_root
  )
}

