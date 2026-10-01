# Direct production-kernel equivalence test for Stage 4.
#
# This is a standalone script (run by R CMD check from tests/).  The direct
# four-mode bridge is kept under benchmarks/stage4 and includes the production
# k_permutation.cpp translation unit.  That benchmark directory is intentionally
# excluded from source tarballs, so the package-check copy reports a visible
# SKIP when the bridge is absent; run this script from the repository checkout
# to exercise all four evaluator modes.

find_repo_root <- function(max_levels = 6L) {
  at <- "."
  for (i in seq_len(max_levels)) {
    if (file.exists(file.path(at, "src", "k_permutation.cpp")) &&
        file.exists(file.path(at, "benchmarks", "stage4",
                              "test_k_four_modes.cpp"))) {
      return(at)
    }
    at <- file.path(at, "..")
  }
  NA_character_
}

repo_root <- find_repo_root()
if (is.na(repo_root)) {
  cat(paste(
    "STAGE4 HYBRID EQUIVALENCE SKIPPED:",
    "benchmarks/stage4/test_k_four_modes.cpp is absent from this source tree.",
    "Run from the development checkout for direct four-mode coverage.\n"
  ))
} else {
  if (!requireNamespace("ape", quietly = TRUE) ||
      !requireNamespace("Rcpp", quietly = TRUE)) {
    stop("The hybrid equivalence test requires ape and Rcpp.", call. = FALSE)
  }
  if (!requireNamespace("fastphylosig", quietly = TRUE)) {
    stop("The hybrid equivalence test requires fastphylosig to be installed.",
         call. = FALSE)
  }

  # Compile from an ASCII-only temporary mirror.  Besides keeping the
  # production translation unit and its relative include layout intact, this
  # avoids toolchain locale failures when the checkout path contains CJK text.
  compile_root <- tempfile("fastphylosig-stage4-")
  tryCatch({
  dir.create(compile_root)
  source_dir <- file.path(repo_root, "src")
  compile_src <- file.path(compile_root, "src")
  dir.create(compile_src)
  source_files <- list.files(
    source_dir, pattern = "\\.(cpp|h|hpp)$", full.names = TRUE,
    recursive = TRUE
  )
  if (!length(source_files) ||
      !all(file.copy(source_files, compile_src, overwrite = TRUE))) {
    stop("Could not stage the production C++ sources for the test harness.",
         call. = FALSE)
  }
  compile_harness_dir <- file.path(compile_root, "benchmarks", "stage4")
  dir.create(compile_harness_dir, recursive = TRUE)
  harness <- file.path(compile_harness_dir, "test_k_four_modes.cpp")
  if (!file.copy(file.path(repo_root, "benchmarks", "stage4",
                           "test_k_four_modes.cpp"),
                 harness, overwrite = TRUE)) {
    stop("Could not stage the four-mode C++ harness.", call. = FALSE)
  }

  test_env <- new.env(parent = globalenv())
  Rcpp::sourceCpp(harness, env = test_env, rebuild = TRUE, showOutput = TRUE)

  make_tree <- function(kind, seed) {
    set.seed(seed)
    if (kind == "balanced") {
      tree <- ape::stree(n = 32L, type = "balanced")
    } else if (kind == "pectinate") {
      tree <- ape::stree(n = 32L, type = "left")
    } else {
      tree <- ape::rtree(32L)
    }
    if (is.null(tree$edge.length)) {
      tree$edge.length <- rep(1, nrow(tree$edge))
    }
    tree$edge.length <- pmax(as.numeric(tree$edge.length), 1e-10)
    tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
    ape::reorder.phylo(tree, "cladewise")
  }

  make_traits <- function(tree, seed) {
    set.seed(seed)
    n <- ape::Ntip(tree)
    X <- cbind(
      continuous = stats::rnorm(n),
      binary = as.numeric(sample.int(2L, n, replace = TRUE) - 1L),
      categorical = as.numeric(sample.int(4L, n, replace = TRUE))
    )
    rownames(X) <- tree$tip.label
    X
  }

  make_permutations <- function(n, nsim, seed) {
    set.seed(seed)
    out <- matrix(NA_integer_, nrow = nsim, ncol = n)
    out[1L, ] <- seq_len(n)
    if (nsim > 1L) {
      for (i in 2:nsim) out[i, ] <- sample.int(n)
    }
    out
  }

  bit_difference <- function(a, b) {
    a <- as.numeric(a)
    b <- as.numeric(b)
    if (length(a) != length(b)) return(Inf)
    ra <- writeBin(a, raw(), size = 8L, endian = .Platform$endian)
    rb <- writeBin(b, raw(), size = 8L, endian = .Platform$endian)
    xor <- bitwXor(as.integer(ra), as.integer(rb))
    sum(vapply(xor, function(value) {
      sum(as.integer(intToBits(value)))
    }, numeric(1)))
  }

  max_abs_difference <- function(a, b) {
    a <- as.numeric(a)
    b <- as.numeric(b)
    if (length(a) != length(b)) return(Inf)
    if (!length(a)) return(0)
    if (!identical(is.finite(a), is.finite(b))) return(Inf)
    finite <- is.finite(a) & is.finite(b)
    if (!any(finite)) return(0)
    max(abs(a[finite] - b[finite]))
  }

  compare_mode <- function(reference, candidate, tree, trait, mode) {
    k_ref <- c(reference$K_obs, as.numeric(reference$K_null))
    k_got <- c(candidate$K_obs, as.numeric(candidate$K_null))
    p_ref <- as.numeric(reference$P)
    p_got <- as.numeric(candidate$P)
    all_ref <- c(k_ref, p_ref)
    all_got <- c(k_got, p_got)
    data.frame(
      tree = tree,
      trait = trait,
      mode = mode,
      bit_difference = bit_difference(k_ref, k_got),
      max_numerical_difference = max_abs_difference(all_ref, all_got),
      P_difference = sum(p_ref != p_got),
      exceedance_difference = sum(reference$exceedance !=
                                    candidate$exceedance),
      tail_decision_difference = sum(reference$tail_decision !=
                                       candidate$tail_decision),
      stringsAsFactors = FALSE
    )
  }

  comparison_rows <- list()
  row_id <- 0L
  for (tree_kind in c("balanced", "pectinate", "random")) {
    tree <- make_tree(tree_kind, seed = 20260928L + match(
      tree_kind, c("balanced", "pectinate", "random")
    ))
    X <- make_traits(tree, seed = 20261001L + ape::Ntip(tree))
    permutations <- make_permutations(
      ape::Ntip(tree), nsim = 17L,
      seed = 20261023L + match(tree_kind, c("balanced", "pectinate", "random"))
    )
    compiled <- fastphylosig:::compile_tree_cpp(
      tree$edge, as.numeric(tree$edge.length), ape::Ntip(tree)
    )
    out <- test_env$stage4_test_k_four_modes(compiled, X, permutations)
    stopifnot(identical(out$permutations, permutations))

    reference <- out$legacy_block1
    for (mode in c("hybrid_block1", "legacy_block4", "hybrid_block4")) {
      candidate <- out[[mode]]
      for (trait_index in seq_len(ncol(X))) {
        trait_comparison <- compare_mode(
          reference = list(
            K_obs = reference$K_obs[trait_index],
            K_null = reference$K_null[, trait_index, drop = FALSE],
            P = reference$P[trait_index],
            exceedance = reference$exceedance[trait_index],
            tail_decision = reference$tail_decision[, trait_index, drop = FALSE]
          ),
          candidate = list(
            K_obs = candidate$K_obs[trait_index],
            K_null = candidate$K_null[, trait_index, drop = FALSE],
            P = candidate$P[trait_index],
            exceedance = candidate$exceedance[trait_index],
            tail_decision = candidate$tail_decision[, trait_index, drop = FALSE]
          ),
          tree = tree_kind,
          trait = colnames(X)[trait_index],
          mode = mode
        )
        row_id <- row_id + 1L
        comparison_rows[[row_id]] <- trait_comparison

        # Stage 3C acceptance: the hybrid may change the last few floating
        # point bits, but it must stay within 1e-10 relative; the exact upper
        # tail, integer exceedance count, and P value must remain unchanged.
        scale <- max(1, max(abs(c(
          reference$K_obs[trait_index],
          reference$K_null[, trait_index],
          reference$P[trait_index]
        ))))
        stopifnot(
          trait_comparison$max_numerical_difference <= 1e-10 * scale,
          trait_comparison$P_difference == 0L,
          trait_comparison$exceedance_difference == 0L,
          trait_comparison$tail_decision_difference == 0L
        )
      }
    }

    # Changing only the fusion block must remain bitwise invariant within each
    # evaluator implementation.
    legacy_big <- out$legacy_block4
    hybrid_big <- out$hybrid_block4
    stopifnot(
      bit_difference(out$legacy_block1$K_obs, legacy_big$K_obs) == 0,
      bit_difference(out$legacy_block1$K_null, legacy_big$K_null) == 0,
      identical(out$legacy_block1$P, legacy_big$P),
      identical(out$legacy_block1$exceedance, legacy_big$exceedance),
      identical(out$legacy_block1$tail_decision, legacy_big$tail_decision),
      bit_difference(out$hybrid_block1$K_obs, hybrid_big$K_obs) == 0,
      bit_difference(out$hybrid_block1$K_null, hybrid_big$K_null) == 0,
      identical(out$hybrid_block1$P, hybrid_big$P),
      identical(out$hybrid_block1$exceedance, hybrid_big$exceedance),
      identical(out$hybrid_block1$tail_decision, hybrid_big$tail_decision)
    )
  }

  results <- do.call(rbind, comparison_rows)
  cat("STAGE4 HYBRID EQUIVALENCE\n")
  print(results, row.names = FALSE)
  cat(sprintf(
    "summary: %d topology-trait-mode comparisons; max numerical difference=%.17g; total tail decision difference=%d\n",
    nrow(results), max(results$max_numerical_difference),
    sum(results$tail_decision_difference)
  ))
  cat("Four direct modes: hybrid=false/true x block=1/4 all executed.\n")
  cat("STAGE4 HYBRID EQUIVALENCE PASS\n")
  }, finally = unlink(compile_root, recursive = TRUE))
}
