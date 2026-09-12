#!/usr/bin/env Rscript

# Stage 2E1A correctness gate for the private D phase harness.
# No production source, package binding, or test expectation is modified.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2e1a_correctness.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
args <- commandArgs(trailingOnly = TRUE)
formal <- "--formal" %in% args
smoke <- "--smoke" %in% args || !formal
args <- args[!args %in% c("--formal", "--smoke")]
repo <- if (length(args)) normalizePath(args[[1L]], winslash = "/", mustWork = TRUE) else
  normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(repo, "benchmarks", "stage2e1a", "results", "correctness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "stage2e1a_common.R"), local = TRUE)

cpp_file <- stage2e1a_env(
  "FASTPHYLOSIG_STAGE2E1A_CPP",
  file.path(repo, "benchmarks", "stage2e1a", "prototype", "d_phase_harness.cpp")
)
compiled <- stage2e1a_compile(repo, cpp_file)
e <- compiled$env
if (!is.function(e$stage2e1a_brownian_generation) ||
    !is.function(e$stage2e1a_brownian_phase) ||
    !is.function(e$stage2e1a_random_null) ||
    !is.function(e$stage2e1a_complete_d)) {
  stop("Stage 2E1A harness exports are incomplete.", call. = FALSE)
}

n_grid <- stage2e1a_grid("FASTPHYLOSIG_STAGE2E1A_CORRECTNESS_N",
                         c(2L, 4L, 7L, 17L), c(2L, 4L), minimum = 2L,
                         is_formal = formal)
nsim_grid <- stage2e1a_grid("FASTPHYLOSIG_STAGE2E1A_CORRECTNESS_NSIM",
                            c(3L, 7L, 19L), c(3L), minimum = 1L,
                            is_formal = formal)
shapes <- strsplit(stage2e1a_env(
  "FASTPHYLOSIG_STAGE2E1A_CORRECTNESS_SHAPES",
  if (isTRUE(formal)) "balanced,random,pectinate" else "balanced,pectinate"
), ",", fixed = TRUE)[[1L]]
shapes <- unique(trimws(shapes[nzchar(trimws(shapes))]))
if (any(!shapes %in% c("balanced", "random", "pectinate"))) {
  stop("unknown Stage 2E1A correctness shape.", call. = FALSE)
}
base_seed <- as.integer(stage2e1a_env("FASTPHYLOSIG_STAGE2E1A_SEED", "206021"))
# Match the frozen public fast_d() default so the audit harness and the public
# route consume random and Brownian draws in the same chunk order.
chunk <- as.integer(stage2e1a_env("FASTPHYLOSIG_STAGE2E1A_CHUNK", "128"))
if (is.na(base_seed) || is.na(chunk) || chunk < 1L)
  stop("invalid Stage 2E1A seed/chunk.", call. = FALSE)

rows <- list()
add <- function(check, pass, shape, n, nsim, details = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    check = check, pass = isTRUE(pass), shape = shape,
    n = as.integer(n), nsim = as.integer(nsim), details = details,
    stringsAsFactors = FALSE
  )
}

for (shape in shapes) for (n in n_grid) {
  tree <- stage2e1a_make_tree(n, shape)
  x <- stats::setNames(stage2e1a_make_observed(tree), tree$tip.label)
  ctx <- fastphylosig::prepare_tree(tree)
  X <- matrix(as.numeric(x), ncol = 1L,
              dimnames = list(names(x), "trait1"))
  analysis <- getFromNamespace(".prepare_analysis", "fastphylosig")(
    ctx, X, signal = "D", data_kind = "binary", verbose = FALSE
  )
  group <- analysis$groups[[1L]]$group
  d_tree <- group$d_tree
  edge <- matrix(as.integer(d_tree$edge), ncol = 2L)
  edge_length <- as.numeric(d_tree$edge.length)
  binary <- getFromNamespace(".binary_state", "fastphylosig")(x)
  observed <- as.numeric(binary$values[d_tree$tip.label])
  prop <- binary$prop_state1
  for (nsim in nsim_grid) {
    seed <- base_seed + n + nsim + match(shape, shapes)
    message("[stage2e1a] correctness shape=", shape,
            " n=", n, " nsim=", nsim)

    generated <- stage2e1a_capture({
      stage2e1a_restore_seed(seed, e$stage2e1a_brownian_generation(
        edge, edge_length, n, nsim, prop, return_states = TRUE,
        n_threads = 1L
      ))
    })
    if (generated$status != "ok") {
      add("brownian_generation_call", FALSE, shape, n, nsim,
          generated$error_message)
      next
    }
    states <- generated$value$states
    add("brownian_generation_finite", all(is.finite(states)), shape, n, nsim)
    add("brownian_rng_draw_count",
        identical(as.integer(generated$value$rng_draws), as.integer(nrow(edge) * nsim)),
        shape, n, nsim)

    replay <- stage2e1a_restore_seed(seed, e$stage2e1a_brownian_generation(
      edge, edge_length, n, nsim, prop, return_states = TRUE, n_threads = 1L
    ))
    add("brownian_same_seed_replay", stage2e1a_exact(states, replay$states),
        shape, n, nsim)

    sorted <- e$stage2e1a_brownian_phase(
      states, prop, phase = "sort_threshold", thresholds = NULL,
      edge = edge, edge_length = edge_length, n_tip = n, n_threads = 1L
    )
    q_r <- apply(states, 2L, stage2e1a_production_type7,
                 probability = prop)
    threshold_delta <- max(abs(as.numeric(sorted$thresholds) - as.numeric(q_r)))
    add("type7_threshold", stage2e1a_numeric_close(
      sorted$thresholds, q_r, ulps = 32L
    ), shape, n, nsim, paste0(
      "max_abs_difference=", format(threshold_delta, scientific = TRUE),
      ";cpp=", paste(format(as.numeric(sorted$thresholds), digits = 17), collapse = ","),
      ";r=", paste(format(as.numeric(q_r), digits = 17), collapse = ",")
    ))

    binary <- e$stage2e1a_brownian_phase(
      states, prop, phase = "binary_d", thresholds = sorted$thresholds,
      edge = edge, edge_length = edge_length, n_tip = n, n_threads = 1L,
      return_payload = TRUE
    )
    expected_binary <- states < rep(as.numeric(sorted$thresholds), each = n)
    expected_binary <- matrix(as.integer(expected_binary), nrow = n, ncol = nsim)
    add("thresholded_binary_states", stage2e1a_exact(binary$binary, expected_binary),
        shape, n, nsim)
    production_threshold <- getFromNamespace(
      "brownian_threshold_cpp", "fastphylosig"
    )(states, prop)
    add("production_brownian_threshold_binary",
        stage2e1a_exact(production_threshold, expected_binary), shape, n, nsim)
    production_D <- getFromNamespace("phylo_d_sums_cpp", "fastphylosig")(
      expected_binary * 1, edge, edge_length, n, 1L
    )
    add("binary_D_matches_production",
        stage2e1a_exact(as.numeric(binary$D), as.numeric(production_D)),
        shape, n, nsim)

    complete <- stage2e1a_restore_seed(seed, e$stage2e1a_complete_d(
      edge, edge_length, observed, n, nsim, prop, n_threads = 1L,
      simulation_chunk = chunk, return_sim = TRUE
    ))
    # The complete pipeline uses a fresh Brownian stream after random draws;
    # compare independently generated Brownian states via a stream replay
    # with the exact same C++ call order below.
    add("complete_random_finite", all(is.finite(complete$random_D)), shape, n, nsim)
    add("complete_brownian_finite", all(is.finite(complete$brownian_D)), shape, n, nsim)
    add("complete_accounting_requested",
        identical(as.integer(complete$nsim_requested), as.integer(nsim)),
        shape, n, nsim)
    add("complete_accounting_success_random",
        identical(as.integer(complete$nsim_successful_random), as.integer(nsim)),
        shape, n, nsim)
    add("complete_accounting_success_brownian",
        identical(as.integer(complete$nsim_successful_brownian), as.integer(nsim)),
        shape, n, nsim)

    # A second same-seed complete call must replay both null streams and every
    # derived summary.  This is the authoritative RNG/replay gate.
    replay_complete <- stage2e1a_restore_seed(seed, e$stage2e1a_complete_d(
      edge, edge_length, observed, n, nsim, prop, n_threads = 1L,
      simulation_chunk = chunk, return_sim = TRUE
    ))
    add("complete_same_seed_random_replay",
        stage2e1a_exact(complete$random_D, replay_complete$random_D), shape, n, nsim)
    add("complete_same_seed_brownian_replay",
        stage2e1a_exact(complete$brownian_D, replay_complete$brownian_D), shape, n, nsim)
    for (field in c("mean_random", "mean_brownian", "p_random", "p_brownian",
                    "mcse_random", "mcse_brownian")) {
      add(paste0("complete_summary_", field),
          stage2e1a_exact(complete[[field]], replay_complete[[field]]),
          shape, n, nsim)
    }

    random_phase <- stage2e1a_restore_seed(seed, e$stage2e1a_random_null(
      edge, edge_length, observed, n, nsim, n_threads = 1L,
      simulation_chunk = chunk, return_sim = TRUE
    ))
    add("random_phase_replay_matches_complete",
        stage2e1a_exact(random_phase$D, complete$random_D), shape, n, nsim)

    # Public production path: same prepared binary tree, same seed, and the
    # same retained null vectors.  This checks P/MCSE/accounting and output
    # contract without asking the harness to approximate R packaging.
    public_online <- stage2e1a_restore_seed(seed, fastphylosig::fast_d(
      ctx, x, test = TRUE, nsim = nsim, return_sim = FALSE,
      keep_null = FALSE, ncores = 1L, progress = FALSE, verbose = FALSE
    ))
    complete_D <- (as.numeric(complete$observed) -
                   as.numeric(complete$mean_brownian)) /
      (as.numeric(complete$mean_random) -
       as.numeric(complete$mean_brownian))
    add("public_online_summary_exact", all(
      stage2e1a_exact(as.numeric(public_online$Parameters$Observed),
                      as.numeric(complete$observed)),
      stage2e1a_exact(as.numeric(public_online$Parameters$MeanRandom),
                      as.numeric(complete$mean_random)),
      stage2e1a_exact(as.numeric(public_online$Parameters$MeanBrownian),
                      as.numeric(complete$mean_brownian)),
      stage2e1a_exact(as.numeric(public_online$DEstimate), complete_D),
      stage2e1a_exact(as.numeric(public_online$P_random),
                      as.numeric(complete$p_random)),
      stage2e1a_exact(as.numeric(public_online$P_Brownian),
                      as.numeric(complete$p_brownian)),
      stage2e1a_exact(as.numeric(public_online$MCSE_P_random),
                      as.numeric(complete$mcse_random)),
      stage2e1a_exact(as.numeric(public_online$MCSE_P_Brownian),
                      as.numeric(complete$mcse_brownian))
    ), shape, n, nsim)

    public <- stage2e1a_restore_seed(seed, fastphylosig::fast_d(
      ctx, x, test = TRUE, nsim = nsim, return_sim = TRUE,
      keep_null = TRUE, ncores = 1L, progress = FALSE, verbose = FALSE
    ))
    d_null_summary <- getFromNamespace(".d_null_summary", "fastphylosig")
    vector_random <- d_null_summary(
      complete$random_D, complete$observed, direction = "less"
    )
    vector_brownian <- d_null_summary(
      complete$brownian_D, complete$observed, direction = "greater"
    )
    add("public_random_null_matches",
        stage2e1a_exact(as.numeric(public$Permutations$random),
                        as.numeric(complete$random_D)), shape, n, nsim)
    add("public_brownian_null_matches",
        stage2e1a_exact(as.numeric(public$Permutations$brownian),
                        as.numeric(complete$brownian_D)), shape, n, nsim)
    # `complete$observed` is the raw observed contrast used by the D
    # calibration.  The public `DEstimate` is the calibrated D statistic, so
    # compare the raw value with Parameters$Observed and check the calibrated
    # value against the same accounting equation separately.
    add("public_observed_D_matches",
        stage2e1a_exact(as.numeric(public$Parameters$Observed),
                        as.numeric(complete$observed)), shape, n, nsim)
    vector_D <- (as.numeric(complete$observed) - vector_brownian$mean) /
      (vector_random$mean - vector_brownian$mean)
    add("public_D_estimate_matches",
        stage2e1a_exact(as.numeric(public$DEstimate), vector_D),
        shape, n, nsim)
    add("public_mean_random_matches",
        stage2e1a_exact(as.numeric(public$Parameters$MeanRandom),
                        as.numeric(vector_random$mean)), shape, n, nsim)
    add("public_mean_brownian_matches",
        stage2e1a_exact(as.numeric(public$Parameters$MeanBrownian),
                        as.numeric(vector_brownian$mean)), shape, n, nsim)
    add("public_P_random_matches",
        stage2e1a_exact(as.numeric(public$P_random), as.numeric(complete$p_random)),
        shape, n, nsim)
    add("public_P_brownian_matches",
        stage2e1a_exact(as.numeric(public$P_Brownian), as.numeric(complete$p_brownian)),
        shape, n, nsim)
    add("public_MCSE_random_matches",
        stage2e1a_exact(as.numeric(public$MCSE_P_random), as.numeric(complete$mcse_random)),
        shape, n, nsim)
    add("public_MCSE_brownian_matches",
        stage2e1a_exact(as.numeric(public$MCSE_P_Brownian), as.numeric(complete$mcse_brownian)),
        shape, n, nsim)
    add("public_random_success_count",
        identical(as.integer(public$nsim_successful_random),
                  as.integer(complete$nsim_successful_random)), shape, n, nsim)
    add("public_brownian_success_count",
        identical(as.integer(public$nsim_successful_brownian),
                  as.integer(complete$nsim_successful_brownian)), shape, n, nsim)

    two <- stage2e1a_restore_seed(seed, e$stage2e1a_complete_d(
      edge, edge_length, observed, n, nsim, prop, n_threads = 2L,
      simulation_chunk = chunk, return_sim = TRUE
    ))
    add("two_thread_random_determinism", stage2e1a_exact(two$random_D,
        complete$random_D), shape, n, nsim)
    add("two_thread_brownian_determinism", stage2e1a_exact(two$brownian_D,
        complete$brownian_D), shape, n, nsim)
  }
}

# Frozen rooted-polytomy public compatibility smoke with supplied states.  The
# binary harness intentionally rejects polytomies; this check ensures the
# production public route remains the owner of that contract.
poly <- list(
  edge = matrix(c(5L, 1L, 5L, 2L, 5L, 3L, 5L, 4L), ncol = 2L, byrow = TRUE),
  edge.length = rep(1, 4), tip.label = paste0("sp", 1:4), Nnode = 1L
)
class(poly) <- "phylo"
poly_x <- stats::setNames(c(0, 1, 0, 1), poly$tip.label)
poly_random <- matrix(c(0, 1, 1, 0, 1, 0, 0, 1), nrow = 4L, ncol = 2L)
poly_brown <- matrix(c(0.1, -0.2, 0.3, 0.2, -0.4, 0.5, -0.1, 0.8),
                     nrow = 4L, ncol = 2L)
poly_public <- stage2e1a_capture(fastphylosig::fast_d(
  poly, poly_x, test = TRUE, nsim = 2L, random_states = poly_random,
  brownian_states = poly_brown, return_sim = TRUE, progress = FALSE,
  verbose = FALSE
))
add("rooted_polytomy_public_contract", poly_public$status == "ok",
    "polytomy", 4L, 2L, if (poly_public$status == "ok") "" else poly_public$error_message)

# Exercise the literal strict tails at one representable step around 1.0.
# Equality must remain outside both tails; no epsilon is introduced.
d_null_summary <- getFromNamespace(".d_null_summary", "fastphylosig")
ulp_up <- 1 + .Machine$double.eps
ulp_down <- 1 - .Machine$double.eps / 2
brownian_ulp <- d_null_summary(c(1, ulp_up), 1, direction = "greater")
random_ulp <- d_null_summary(c(1, ulp_down), 1, direction = "less")
add("one_ulp_brownian_strict_tail", identical(brownian_ulp$p, 0.5),
    "tail-boundary", 2L, 2L)
add("one_ulp_random_strict_tail", identical(random_ulp$p, 0.5),
    "tail-boundary", 2L, 2L)

result <- if (length(rows)) do.call(rbind, rows) else
  data.frame(check = character(), pass = logical(), shape = character(),
             n = integer(), nsim = integer(), details = character(),
             stringsAsFactors = FALSE)
utils::write.csv(result, file.path(out_dir, "stage2e1a_correctness.csv"),
                 row.names = FALSE)
status <- data.frame(
  formal = if (formal) "YES" else "NO_SMOKE",
  checks = nrow(result), failures = sum(!result$pass),
  status = if (nrow(result) && all(result$pass)) "PASS" else "FAIL",
  production_code_changed = "NO", stringsAsFactors = FALSE
)
utils::write.csv(status, file.path(out_dir, "stage2e1a_correctness_status.csv"),
                 row.names = FALSE)
stage2e1a_write_provenance(repo, out_dir, cpp_file, list(formal = formal))
print(status)
if (!nrow(result) || any(!result$pass)) {
  bad <- result[!result$pass, , drop = FALSE]
  print(utils::head(bad, 40L))
  stop("Stage 2E1A correctness gate failed.", call. = FALSE)
}
message("[stage2e1a] correctness PASS: ", nrow(result), " checks")
