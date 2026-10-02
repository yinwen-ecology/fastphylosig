# Matched single-call benchmark: fastphylosig 0.3.0 vs D/Delta references.
#
# Run from the repository root with Rscript --vanilla and an explicit
# --library PATH. The default profile is a small n=10/20 smoke test; request
# --profile full only after review/authorization. All generated trees, traits,
# fast tree contexts, and reference-side comparative.data objects are prepared
# before the timer. Each formal timer covers exactly one public call.

options(stringsAsFactors = FALSE, warn = 1)

parse_cli <- function(args) {
  profile <- "smoke"
  library_path <- NULL
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--profile")) {
      if (i == length(args)) stop("--profile needs smoke or full.", call. = FALSE)
      profile <- args[[i + 1L]]
      i <- i + 2L
    } else if (startsWith(arg, "--profile=")) {
      profile <- sub("^--profile=", "", arg)
      i <- i + 1L
    } else if (identical(arg, "--library")) {
      if (i == length(args)) stop("--library needs a path.", call. = FALSE)
      library_path <- args[[i + 1L]]
      i <- i + 2L
    } else if (startsWith(arg, "--library=")) {
      library_path <- sub("^--library=", "", arg)
      i <- i + 1L
    } else {
      stop("Unknown benchmark argument: ", arg, call. = FALSE)
    }
  }
  if (!profile %in% c("smoke", "full")) {
    stop("--profile must be smoke or full.", call. = FALSE)
  }
  if (is.null(library_path) || !nzchar(library_path)) {
    stop("Usage: Rscript --vanilla benchmark_reference_d_delta.R ",
         "--profile smoke|full --library PATH", call. = FALSE)
  }
  list(profile = profile, library_path = library_path)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
profile <- cli$profile
repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run this script from the fastphylosig repository root.", call. = FALSE)
}
out_dir <- file.path(repo_root, "docs", "release-0.3.0")
if (!dir.exists(out_dir)) stop("Missing benchmark output directory.", call. = FALSE)

install_lib <- normalizePath(cli$library_path, winslash = "/", mustWork = TRUE)
expected_package_dir <- normalizePath(
  file.path(install_lib, "fastphylosig"), winslash = "/", mustWork = TRUE
)
.libPaths(unique(c(install_lib, .libPaths())))

required_packages <- c("ape", "caper", "digest", "expm")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

tarball_raw <- Sys.getenv(
  "FASTPHYLOSIG_SOURCE_TARBALL",
  unset = file.path(
    repo_root, "benchmarks", "stage5", "documentation_fix_release_20261001",
    "fastphylosig_0.3.0.tar.gz"
  )
)
tarball <- normalizePath(tarball_raw, winslash = "/", mustWork = TRUE)
expected_tarball_sha256 <-
  "909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff"
observed_tarball_sha256 <- digest::digest(
  tarball, algo = "sha256", serialize = FALSE, file = TRUE
)
if (!identical(tolower(observed_tarball_sha256), expected_tarball_sha256)) {
  stop("Release tarball SHA-256 mismatch; no benchmark was run.", call. = FALSE)
}

if ("fastphylosig" %in% loadedNamespaces()) {
  stop("fastphylosig is already loaded; use a fresh Rscript --vanilla process.",
       call. = FALSE)
}
suppressPackageStartupMessages(library(
  "fastphylosig", character.only = TRUE, lib.loc = install_lib
))
package_path <- normalizePath(
  find.package("fastphylosig"), winslash = "/", mustWork = TRUE
)
if (!identical(tolower(package_path), tolower(expected_package_dir))) {
  stop("Loaded package path is not the explicitly selected installation.",
       call. = FALSE)
}
loaded_version <- as.character(utils::packageVersion("fastphylosig"))
if (!identical(loaded_version, "0.3.0")) {
  stop("Expected fastphylosig 0.3.0; loaded ", loaded_version, ".",
       call. = FALSE)
}
dll_list <- getLoadedDLLs()
dll_names <- names(dll_list)
dll_index <- which(grepl("^fastphylosig(\\.dll)?$", dll_names,
                         ignore.case = TRUE))
if (length(dll_index) != 1L) {
  stop("Could not identify exactly one loaded fastphylosig DLL.", call. = FALSE)
}
dll_path <- normalizePath(dll_list[[dll_index]][["path"]], winslash = "/",
                          mustWork = TRUE)
if (!startsWith(tolower(dll_path), paste0(tolower(package_path), "/"))) {
  stop("Loaded DLL is outside the selected package installation.", call. = FALSE)
}
if (!all(c("fast_d", "fast_delta", "prepare_tree") %in%
         getNamespaceExports("fastphylosig"))) {
  stop("Selected 0.3.0 library lacks a required public API.", call. = FALSE)
}

delta_source <- file.path(
  repo_root, "paper", "reference_sources", "delta_borges2019_code.R"
)
delta_source_sha256 <-
  "0dc1e987c13f0e972b0e7c3b4829b4d0fdc41ff0de3856bae4803f9ce0115175"
delta_source_url <- paste0(
  "https://raw.githubusercontent.com/mrborges23/",
  "delta_statistic/master/code.R"
)
if (!file.exists(delta_source)) {
  stop("Missing pinned Borges Delta reference source.", call. = FALSE)
}
observed_delta_source_sha256 <- digest::digest(
  delta_source, algo = "sha256", serialize = FALSE, file = TRUE
)
if (!identical(tolower(observed_delta_source_sha256), delta_source_sha256)) {
  stop("Borges Delta source SHA-256 mismatch; no benchmark was run.",
       call. = FALSE)
}
suppressPackageStartupMessages({
  delta_reference_env <- new.env(parent = globalenv())
  sys.source(delta_source, envir = delta_reference_env)
})
if (!exists("delta", envir = delta_reference_env, mode = "function",
            inherits = FALSE)) {
  stop("Pinned Borges source does not define delta().", call. = FALSE)
}
delta_reference <- get("delta", envir = delta_reference_env, inherits = FALSE)

package_path_label <- "isolated_release_library/fastphylosig"
dll_path_label <- "isolated_release_library/fastphylosig/libs/<platform>/fastphylosig"
package_path_sha256 <- digest::digest(package_path, algo = "sha256",
                                      serialize = FALSE)
dll_path_sha256 <- digest::digest(dll_path, algo = "sha256", serialize = FALSE)
dll_file_sha256 <- digest::digest(
  dll_path, algo = "sha256", serialize = FALSE, file = TRUE
)
repo_prefix <- paste0(tolower(repo_root), "/")
tarball_relative <- if (startsWith(tolower(tarball), repo_prefix)) {
  substring(tarball, nchar(repo_root) + 2L)
} else {
  paste0("<external>/", basename(tarball))
}

run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
output_prefix <- if (profile == "smoke") {
  paste0("benchmark_reference_d_delta_smoke_", run_id)
} else {
  "benchmark_reference_d_delta"
}
output_names <- paste0(output_prefix, c(
  "_runs.csv", "_pairs.csv", "_summary.csv", "_metadata.csv",
  "_warmups.csv", "_session_info.txt"
))
output_paths <- file.path(out_dir, output_names)
if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite existing D/Delta benchmark output(s): ",
       paste(output_names[file.exists(output_paths)], collapse = ", "),
       call. = FALSE)
}
run_file <- output_paths[[1L]]
warmup_file <- output_paths[[5L]]

full_n_grid <- c(50L, 100L, 500L, 2000L)
n_grid <- if (profile == "smoke") c(10L, 20L) else full_n_grid
formal_reps <- if (profile == "smoke") 1L else 10L
scenario_grid <- data.frame(
  scenario = c("none", "moderate", "strong"),
  signal_mixture_h2 = c(0, 0.5, 0.95),
  stringsAsFactors = FALSE
)
d_nsim <- if (profile == "smoke") 19L else 199L
delta_mcmc_sim <- if (profile == "smoke") 1000L else 10000L
delta_thin <- 10L
delta_burn <- 100L
delta_lambda0 <- 0.1
delta_proposal_sd <- 0.5
delta_entropy <- "LSE"
delta_model <- "ARD"
base_seed <- 20260901L
methods <- c("D", "Delta")
expected_pair_count <- length(n_grid) * length(methods) *
  nrow(scenario_grid) * formal_reps

set_reproducible_seed <- function(seed) {
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed)) {
    stop("seed must be one finite integer.", call. = FALSE)
  }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  invisible(as.integer(seed))
}

hash_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

make_tree <- function(n, seed) {
  set_reproducible_seed(seed)
  tree <- ape::rtree(n = n, rooted = TRUE)
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  branch_lengths <- as.numeric(tree$edge.length)
  branch_lengths[!is.finite(branch_lengths) | branch_lengths <= 0] <- 1
  tree$edge.length <- branch_lengths
  tree$tip.label <- paste0("sp_", seq_len(ape::Ntip(tree)))
  ape::reorder.phylo(tree, order = "cladewise")
}

standardize_component <- function(x, label) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x))) {
    stop(label, " component is not finite.", call. = FALSE)
  }
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx <= 0) {
    stop(label, " component has zero or undefined variance.", call. = FALSE)
  }
  (x - mean(x)) / sx
}

make_scenario_traits <- function(tree, bm_seed, noise_seed, h2) {
  set_reproducible_seed(bm_seed)
  bm <- standardize_component(
    ape::rTraitCont(tree, model = "BM", ancestor = FALSE), "Brownian"
  )
  names(bm) <- tree$tip.label
  bm <- bm[tree$tip.label]
  set_reproducible_seed(noise_seed)
  noise <- standardize_component(stats::rnorm(ape::Ntip(tree)), "iid noise")
  names(noise) <- tree$tip.label
  noise <- noise[tree$tip.label]
  continuous <- sqrt(h2) * bm + sqrt(1 - h2) * noise
  names(continuous) <- tree$tip.label
  binary <- as.integer(continuous > stats::median(continuous))
  names(binary) <- tree$tip.label
  ranks <- rank(continuous, ties.method = "first")
  categorical <- cut(
    ranks,
    breaks = c(0, length(ranks) / 3, 2 * length(ranks) / 3, length(ranks)),
    labels = c("a", "b", "c"), include.lowest = TRUE
  )
  names(categorical) <- tree$tip.label
  list(
    D = binary, Delta = categorical,
    bm_component = stats::setNames(bm, tree$tip.label),
    noise_component = stats::setNames(noise, tree$tip.label)
  )
}

assert_trait <- function(tree, trait, method, expected_n) {
  tree_tips <- tree$tip.label
  if (!inherits(tree, "phylo") || ape::Ntip(tree) != expected_n ||
      length(tree_tips) != expected_n) {
    stop(method, ": tree does not have requested tip count.", call. = FALSE)
  }
  if (length(trait) != expected_n || is.null(names(trait)) || anyNA(trait) ||
      !setequal(names(trait), tree_tips) ||
      !identical(names(trait)[match(tree_tips, names(trait))], tree_tips)) {
    stop(method, ": trait is not a complete named match to tree tips.",
         call. = FALSE)
  }
  trait <- trait[tree_tips]
  if (!identical(names(trait), tree_tips)) {
    stop(method, ": trait order does not match tree-tip order.", call. = FALSE)
  }
  if (method == "D" &&
      (any(!is.finite(as.numeric(trait))) || length(unique(trait)) != 2L)) {
    stop("D scenario must contain exactly two finite states.", call. = FALSE)
  }
  if (method == "Delta" && length(unique(trait)) != 3L) {
    stop("Delta scenario must contain exactly three states.", call. = FALSE)
  }
  trait
}

make_fast_input <- function(trait, method) {
  states <- if (method == "D") as.integer(trait) else as.integer(trait)
  matrix(states, ncol = 1L,
         dimnames = list(names(trait), "trait"))
}

fast_call <- function(method, context, trait_matrix) {
  if (method == "D") {
    return(fastphylosig::fast_d(
      tree = context$tree, x = trait_matrix, prepared = context,
      test = TRUE, nsim = d_nsim, return_sim = FALSE, keep_null = FALSE,
      verbose = FALSE, ncores = 1L, progress = FALSE
    ))
  }
  if (method == "Delta") {
    return(fastphylosig::fast_delta(
      tree = context$tree, x = trait_matrix, prepared = context,
      test = FALSE, mcmc_sim = delta_mcmc_sim, thin = delta_thin,
      burn = delta_burn, lambda0 = delta_lambda0,
      proposal_sd = delta_proposal_sd, entropy = delta_entropy,
      model = delta_model, return_sim = FALSE, verbose = FALSE,
      ncores = 1L, progress = FALSE
    ))
  }
  stop("Unknown method: ", method, call. = FALSE)
}

prepare_reference_input <- function(method, tree, trait) {
  if (method == "D") {
    dat <- data.frame(
      species = tree$tip.label,
      trait = as.integer(trait),
      stringsAsFactors = FALSE
    )
    comparative <- caper::comparative.data(
      tree, dat, names.col = "species", warn.dropped = FALSE
    )
    if (!is.data.frame(comparative$data) || is.null(comparative$phy)) {
      stop("caper::comparative.data returned an unexpected object.",
           call. = FALSE)
    }
    data_species <- if ("species" %in% names(comparative$data)) {
      as.character(comparative$data$species)
    } else {
      rownames(comparative$data)
    }
    actual_tree_n <- as.integer(ape::Ntip(comparative$phy))
    actual_data_n <- as.integer(nrow(comparative$data))
    matched <- if (is.null(data_species)) 0L else
      as.integer(length(intersect(tree$tip.label, data_species)))
    if (actual_tree_n != length(tree$tip.label) || actual_data_n !=
        length(tree$tip.label) || matched != length(tree$tip.label) ||
        !setequal(data_species, tree$tip.label)) {
      stop("caper reference preparation changed the species set.",
           call. = FALSE)
    }
    return(list(
      comparative = comparative, trait = dat$trait,
      n_species_used = actual_data_n, matched_species = matched,
      n_removed_na = as.integer(sum(is.na(dat$trait))),
      input_tree_tip_count = actual_tree_n,
      input_trait_length = actual_data_n,
      reference_data_rows = actual_data_n,
      reference_tree_tip_count = actual_tree_n
    ))
  }

  trait_names <- names(trait)
  matched <- if (is.null(trait_names)) 0L else
    as.integer(sum(tree$tip.label %in% trait_names))
  n_removed_na <- as.integer(sum(is.na(trait)))
  n_used <- as.integer(sum(tree$tip.label %in% trait_names &
                           !is.na(trait[match(tree$tip.label, trait_names)])))
  if (is.null(trait_names) || !identical(trait_names, tree$tip.label) ||
      length(trait) != ape::Ntip(tree) || n_removed_na != 0L ||
      n_used != ape::Ntip(tree)) {
    stop("Borges reference trait vector is not complete and tip-matched.",
         call. = FALSE)
  }
  list(
    tree = tree, trait = trait, n_species_used = n_used,
    matched_species = matched, n_removed_na = n_removed_na,
    input_tree_tip_count = as.integer(ape::Ntip(tree)),
    input_trait_length = as.integer(length(trait)),
    reference_data_rows = NA_integer_,
    reference_tree_tip_count = as.integer(ape::Ntip(tree))
  )
}

reference_phylo_d <- function(comparative, trait, permut) {
  # caper captures the unevaluated binvar symbol; the formal name `trait`
  # matches the pre-prepared comparative.data column.
  caper::phylo.d(comparative, binvar = trait, permut = permut)
}

reference_call <- function(method, input) {
  if (method == "D") {
    return(reference_phylo_d(input$comparative, input$trait, d_nsim))
  }
  if (method == "Delta") {
    return(delta_reference(
      input$trait, input$tree, delta_lambda0, delta_proposal_sd,
      delta_mcmc_sim, delta_thin, delta_burn
    ))
  }
  stop("Unknown method: ", method, call. = FALSE)
}

get_field <- function(value, candidates, default = NA) {
  for (field in candidates) {
    candidate <- NULL
    if (is.data.frame(value) && field %in% names(value) && nrow(value)) {
      candidate <- value[[field]][[1L]]
    } else if (is.list(value) && !is.null(value[[field]])) {
      candidate <- value[[field]][[1L]]
    }
    if (!is.null(candidate) && length(candidate)) return(candidate)
  }
  default
}

extract_estimate <- function(value, method) {
  if (is.numeric(value) && length(value)) {
    return(suppressWarnings(as.numeric(value[[1L]])))
  }
  candidates <- if (method == "D") {
    c("estimate", "DEstimate", "D_fast", "D")
  } else {
    c("estimate", "delta", "Delta_fast", "Delta")
  }
  candidate <- get_field(value, candidates, default = NA_real_)
  suppressWarnings(as.numeric(candidate)[[1L]])
}

result_fields <- function(value) {
  if (is.data.frame(value)) return(paste(names(value), collapse = ";"))
  if (is.list(value)) return(paste(names(value), collapse = ";"))
  if (is.numeric(value)) return("<numeric>")
  paste0("<", paste(class(value), collapse = "/"), ">")
}

capture_call <- function(fun) {
  warning_messages <- character()
  error_message <- NA_character_
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  list(value = value, warning = paste(unique(warning_messages), collapse = " || "),
       error = error_message)
}

timed_call_wall <- function(fun) {
  started <- Sys.time()
  captured <- capture_call(fun)
  stopped <- Sys.time()
  captured$elapsed_seconds <- as.numeric(
    difftime(stopped, started, units = "secs")
  )
  captured
}

audit_fast_result <- function(value) {
  required <- c("n_species", "matched_species", "n_removed_na", "status")
  if (!is.data.frame(value) || nrow(value) != 1L ||
      !all(required %in% names(value))) {
    return(list(
      n_species_used = NA_integer_, matched_species = NA_integer_,
      n_removed_na = NA_integer_, status = "missing_fields",
      validation_error = "fast result lacks actual one-row species/status fields"
    ))
  }
  n_species_used <- suppressWarnings(as.integer(value$n_species[[1L]]))
  matched_species <- suppressWarnings(as.integer(value$matched_species[[1L]]))
  n_removed_na <- suppressWarnings(as.integer(value$n_removed_na[[1L]]))
  status <- as.character(value$status[[1L]])
  valid <- all(is.finite(c(n_species_used, matched_species, n_removed_na))) &&
    n_species_used >= 0L && matched_species >= 0L && n_removed_na >= 0L
  list(
    n_species_used = n_species_used, matched_species = matched_species,
    n_removed_na = n_removed_na, status = status,
    validation_error = if (valid) NA_character_ else
      "fast result reports invalid species audit counts"
  )
}

audit_fast_delta_control <- function(value) {
  expected <- list(
    lambda0 = delta_lambda0,
    proposal_sd = delta_proposal_sd,
    mcmc_sim = as.integer(delta_mcmc_sim),
    thin = as.integer(delta_thin),
    burn = as.integer(delta_burn),
    entropy = delta_entropy,
    model = delta_model,
    ace_engine = "fast",
    ncores = 1L
  )
  control <- attr(value, "delta_control", exact = TRUE)
  errors <- character()
  if (!is.list(control)) {
    errors <- c(errors, "fast Delta result lacks delta_control attribute")
  } else {
    for (field in names(expected)) {
      observed <- control[[field]]
      target <- expected[[field]]
      equal <- if (is.numeric(target)) {
        length(observed) == 1L && is.numeric(observed) &&
          is.finite(observed) && as.numeric(observed) == as.numeric(target)
      } else {
        length(observed) == 1L && identical(as.character(observed), target)
      }
      if (!isTRUE(equal)) {
        errors <- c(errors, sprintf(
          "delta_control$%s did not match the requested value", field
        ))
      }
    }
  }
  expected_iterations <- as.numeric(2L * delta_mcmc_sim)
  requested_iterations <- suppressWarnings(as.numeric(get_field(
    value, "requested_iterations", default = NA_real_
  )))
  successful_iterations <- suppressWarnings(as.numeric(get_field(
    value, "successful_iterations", default = NA_real_
  )))
  requested_match <- length(requested_iterations) == 1L &&
    is.finite(requested_iterations) &&
    requested_iterations == expected_iterations
  successful_match <- length(successful_iterations) == 1L &&
    is.finite(successful_iterations) &&
    successful_iterations == expected_iterations
  if (!requested_match) {
    errors <- c(errors, "Delta requested_iterations did not equal 2*mcmc_sim")
  }
  if (!successful_match) {
    errors <- c(errors, "Delta successful_iterations did not equal 2*mcmc_sim")
  }
  reported <- lapply(names(expected), function(field) {
    value <- if (is.list(control)) control[[field]] else NULL
    if (is.null(value) || !length(value)) return(NA_character_)
    as.character(value[[1L]])
  })
  names(reported) <- paste0("reported_delta_control_", names(expected))
  c(reported, list(
    delta_expected_total_iterations = expected_iterations,
    delta_reported_requested_iterations = requested_iterations,
    delta_reported_successful_iterations = successful_iterations,
    delta_requested_iterations_match = requested_match,
    delta_successful_iterations_match = successful_match,
    validation_error = if (length(errors)) paste(errors, collapse = "; ") else
      NA_character_
  ))
}

append_csv_row <- function(row, path) {
  utils::write.table(
    row, file = path, sep = ",", row.names = FALSE,
    col.names = !file.exists(path), append = file.exists(path),
    quote = TRUE, qmethod = "double", na = "", fileEncoding = "UTF-8"
  )
}

parameter_value <- function(method, field) {
  if (field == "d_nsim" && method == "D") return(as.integer(d_nsim))
  if (method == "Delta") {
    return(switch(field,
      d_nsim = NA_integer_,
      delta_mcmc_sim = as.integer(delta_mcmc_sim),
      delta_thin = as.integer(delta_thin),
      delta_burn = as.integer(delta_burn),
      delta_lambda0 = delta_lambda0,
      delta_proposal_sd = delta_proposal_sd,
      delta_chains = 2L,
      delta_entropy = delta_entropy,
      delta_model = delta_model,
      NA
    ))
  }
  switch(field,
    d_nsim = as.integer(d_nsim),
    delta_mcmc_sim = NA_integer_,
    delta_thin = NA_integer_,
    delta_burn = NA_integer_,
    delta_lambda0 = NA_real_,
    delta_proposal_sd = NA_real_,
    delta_chains = NA_integer_,
    delta_entropy = NA_character_,
    delta_model = NA_character_,
    NA
  )
}

run_one <- function(implementation, method, context, tree, trait,
                    trait_matrix, reference_input, meta,
                    call_seed, first_implementation, call_order) {
  # RNG setup and garbage collection are intentionally outside the timer.
  set_reproducible_seed(call_seed)
  invisible(gc())
  timed <- if (implementation == "fastphylosig") {
    timed_call_wall(function() fast_call(method, context, trait_matrix))
  } else {
    timed_call_wall(function() reference_call(method, reference_input))
  }
  estimate <- extract_estimate(timed$value, method)
  if (implementation == "fastphylosig") {
    audit <- audit_fast_result(timed$value)
    parameter_audit <- if (method == "Delta") {
      audit_fast_delta_control(timed$value)
    } else {
      list(
        reported_delta_control_lambda0 = NA_character_,
        reported_delta_control_proposal_sd = NA_character_,
        reported_delta_control_mcmc_sim = NA_character_,
        reported_delta_control_thin = NA_character_,
        reported_delta_control_burn = NA_character_,
        reported_delta_control_entropy = NA_character_,
        reported_delta_control_model = NA_character_,
        reported_delta_control_ace_engine = NA_character_,
        reported_delta_control_ncores = NA_character_,
        delta_expected_total_iterations = NA_real_,
        delta_reported_requested_iterations = NA_real_,
        delta_reported_successful_iterations = NA_real_,
        delta_requested_iterations_match = NA,
        delta_successful_iterations_match = NA,
        validation_error = NA_character_
      )
    }
    reported_status <- audit$status
    status_source <- "fast_result"
    n_species_used <- audit$n_species_used
    matched_species <- audit$matched_species
    n_removed_na <- audit$n_removed_na
    input_tree_tip_count <- NA_integer_
    input_trait_length <- NA_integer_
    reference_data_rows <- NA_integer_
    reference_tree_tip_count <- NA_integer_
    validation_errors <- c(audit$validation_error,
                           parameter_audit$validation_error)
    validation_errors <- validation_errors[
      !is.na(validation_errors) & nzchar(validation_errors)
    ]
    validation_error <- if (length(validation_errors)) {
      paste(validation_errors, collapse = "; ")
    } else {
      NA_character_
    }
    status_ok <- identical(reported_status, "ok")
  } else {
    audit <- reference_input
    parameter_audit <- list(
      reported_delta_control_lambda0 = NA_character_,
      reported_delta_control_proposal_sd = NA_character_,
      reported_delta_control_mcmc_sim = NA_character_,
      reported_delta_control_thin = NA_character_,
      reported_delta_control_burn = NA_character_,
      reported_delta_control_entropy = NA_character_,
      reported_delta_control_model = NA_character_,
      reported_delta_control_ace_engine = NA_character_,
      reported_delta_control_ncores = NA_character_,
      delta_expected_total_iterations = NA_real_,
      delta_reported_requested_iterations = NA_real_,
      delta_reported_successful_iterations = NA_real_,
      delta_requested_iterations_match = NA,
      delta_successful_iterations_match = NA,
      validation_error = NA_character_
    )
    n_species_used <- audit$n_species_used
    matched_species <- audit$matched_species
    n_removed_na <- audit$n_removed_na
    input_tree_tip_count <- audit$input_tree_tip_count
    input_trait_length <- audit$input_trait_length
    reference_data_rows <- audit$reference_data_rows
    reference_tree_tip_count <- audit$reference_tree_tip_count
    reported_status <- "not_reported_by_reference"
    status_source <- "finite_reference_estimate_and_input_audit"
    validation_error <- NA_character_
    status_ok <- length(estimate) == 1L && is.finite(estimate)
  }
  call_ok <- is.na(timed$error) && length(estimate) == 1L &&
    is.finite(estimate) && status_ok && is.na(validation_error) &&
    is.finite(timed$elapsed_seconds) && timed$elapsed_seconds > 0
  prefix <- if (implementation == "fastphylosig") "fast_" else "reference_"
  run_row <- data.frame(
    run_id = run_id,
    profile = profile,
    comparison_id = meta$comparison_id,
    implementation = implementation,
    method = method,
    reference = if (method == "D") "caper::phylo.d" else
      "Borges2019 delta() [pinned source]",
    scenario = meta$scenario,
    signal_mixture_h2 = meta$signal_mixture_h2,
    n_tip = meta$n_tip,
    replicate = meta$replicate,
    n_species_used = n_species_used,
    matched_species = matched_species,
    n_removed_na = n_removed_na,
    input_tree_tip_count = input_tree_tip_count,
    input_trait_length = input_trait_length,
    reference_data_rows = reference_data_rows,
    reference_tree_tip_count = reference_tree_tip_count,
    tree_seed = meta$tree_seed,
    component_seed = meta$component_seed,
    bm_seed = meta$bm_seed,
    noise_seed = meta$noise_seed,
    call_seed = as.integer(call_seed),
    first_implementation = first_implementation,
    call_order = as.integer(call_order),
    test = method == "D",
    d_nsim = parameter_value(method, "d_nsim"),
    delta_mcmc_sim = parameter_value(method, "delta_mcmc_sim"),
    delta_thin = parameter_value(method, "delta_thin"),
    delta_burn = parameter_value(method, "delta_burn"),
    delta_lambda0 = parameter_value(method, "delta_lambda0"),
    delta_proposal_sd = parameter_value(method, "delta_proposal_sd"),
    delta_chains = parameter_value(method, "delta_chains"),
    delta_entropy = parameter_value(method, "delta_entropy"),
    delta_model = parameter_value(method, "delta_model"),
    reported_delta_control_lambda0 =
      parameter_audit$reported_delta_control_lambda0,
    reported_delta_control_proposal_sd =
      parameter_audit$reported_delta_control_proposal_sd,
    reported_delta_control_mcmc_sim =
      parameter_audit$reported_delta_control_mcmc_sim,
    reported_delta_control_thin =
      parameter_audit$reported_delta_control_thin,
    reported_delta_control_burn =
      parameter_audit$reported_delta_control_burn,
    reported_delta_control_entropy =
      parameter_audit$reported_delta_control_entropy,
    reported_delta_control_model =
      parameter_audit$reported_delta_control_model,
    reported_delta_control_ace_engine =
      parameter_audit$reported_delta_control_ace_engine,
    reported_delta_control_ncores =
      parameter_audit$reported_delta_control_ncores,
    delta_expected_total_iterations =
      parameter_audit$delta_expected_total_iterations,
    delta_reported_requested_iterations =
      parameter_audit$delta_reported_requested_iterations,
    delta_reported_successful_iterations =
      parameter_audit$delta_reported_successful_iterations,
    delta_requested_iterations_match =
      parameter_audit$delta_requested_iterations_match,
    delta_successful_iterations_match =
      parameter_audit$delta_successful_iterations_match,
    tree_hash = meta$tree_hash,
    bm_hash = meta$bm_hash,
    noise_hash = meta$noise_hash,
    trait_hash = meta$trait_hash,
    reference_input_hash = meta$reference_input_hash,
    package_version = loaded_version,
    package_path = package_path_label,
    package_path_sha256 = package_path_sha256,
    dll_path = dll_path_label,
    dll_path_sha256 = dll_path_sha256,
    dll_file_sha256 = dll_file_sha256,
    source_tarball = tarball_relative,
    source_tarball_sha256 = observed_tarball_sha256,
    reference_version = if (method == "D") {
      as.character(utils::packageVersion("caper"))
    } else {
      paste0("Borges2019-source-sha256:", delta_source_sha256)
    },
    elapsed_seconds = as.numeric(timed$elapsed_seconds),
    estimate = as.numeric(estimate),
    timing_clock = "Sys.time",
    timing_iterations = 1L,
    prepare_excluded = TRUE,
    reported_status = reported_status,
    status_source = status_source,
    success = call_ok,
    validation_error = validation_error,
    result_fields = result_fields(timed$value),
    warning = timed$warning,
    error = timed$error,
    # Optional method-specific diagnostics are returned by fastphylosig; keep
    # them as NA for reference implementations that do not expose them.
    n_states = if (implementation == "fastphylosig")
      suppressWarnings(as.integer(get_field(timed$value, "n_states"))) else
      as.integer(length(unique(trait))),
    nsim_requested = get_field(timed$value,
      c("nsim_requested", "n_sim_requested", "requested_permutations")),
    nsim_successful_random = get_field(timed$value,
      c("nsim_successful_random", "n_sim_successful_random")),
    nsim_successful_brownian = get_field(timed$value,
      c("nsim_successful_brownian", "n_sim_successful_brownian")),
    requested_iterations = get_field(timed$value,
      c("requested_iterations", "n_iter_requested")),
    successful_iterations = get_field(timed$value,
      c("successful_iterations", "n_iter_successful")),
    diagnostics_available = as.character(get_field(
      timed$value, c("diagnostics_available"), default = NA_character_
    )),
    diagnostics_note = as.character(get_field(
      timed$value, c("diagnostics_note"), default = NA_character_
    )),
    ESS_alpha = suppressWarnings(as.numeric(get_field(
      timed$value, c("ESS_alpha", "ess_alpha"), default = NA_real_
    ))),
    ESS_beta = suppressWarnings(as.numeric(get_field(
      timed$value, c("ESS_beta", "ess_beta"), default = NA_real_
    ))),
    Rhat_alpha = suppressWarnings(as.numeric(get_field(
      timed$value, c("Rhat_alpha", "split_Rhat_alpha"), default = NA_real_
    ))),
    Rhat_beta = suppressWarnings(as.numeric(get_field(
      timed$value, c("Rhat_beta", "split_Rhat_beta"), default = NA_real_
    ))),
    MCSE_Delta = suppressWarnings(as.numeric(get_field(
      timed$value, c("MCSE_Delta", "estimate_mcse"), default = NA_real_
    ))),
    stringsAsFactors = FALSE
  )
  append_csv_row(run_row, run_file)
  if (!call_ok) {
    reason <- if (!is.na(timed$error)) timed$error else
      if (!is.finite(timed$elapsed_seconds) || timed$elapsed_seconds <= 0)
        "non-positive/non-finite elapsed time" else
        if (!identical(status_ok, TRUE)) paste0("reported status is ",
          reported_status) else
        if (!is.na(validation_error)) validation_error else
          "estimate was not finite"
    stop("Call failed; raw row retained in ", basename(run_file), ": ",
         meta$comparison_id, " / ", implementation, ": ", reason,
         call. = FALSE)
  }
  invisible(run_row)
}

record_warmup <- function(implementation, method, n_tip, meta, fun) {
  captured <- capture_call(fun)
  estimate <- extract_estimate(captured$value, method)
  status <- if (implementation == "fastphylosig") {
    audit_fast_result(captured$value)$status
  } else {
    "not_reported_by_reference"
  }
  success <- is.na(captured$error) && is.finite(estimate) &&
    (implementation == "reference" || identical(status, "ok"))
  row <- data.frame(
    run_id = run_id, profile = profile, method = method,
    implementation = implementation, n_tip = as.integer(n_tip),
    warmup_id = meta$warmup_id, call_seed = as.integer(meta$call_seed),
    reported_status = status, success = success,
    estimate = as.numeric(estimate), warning = captured$warning,
    error = captured$error, result_fields = result_fields(captured$value),
    stringsAsFactors = FALSE
  )
  append_csv_row(row, warmup_file)
  if (!success) {
    stop("Untimed warm-up failed; warm-up row retained: ", method, "/n",
         n_tip, "/", implementation, call. = FALSE)
  }
  invisible(row)
}

warm_up <- function(n_tip) {
  tree_seed <- base_seed + n_tip
  tree_raw <- make_tree(n_tip, tree_seed)
  context <- fastphylosig::prepare_tree(tree_raw)
  tree <- context$tree
  component_seed <- base_seed + n_tip + 1L
  bm_seed <- component_seed + 1L
  noise_seed <- component_seed + 2L
  traits <- make_scenario_traits(tree, bm_seed, noise_seed, h2 = 0.5)
  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    trait <- assert_trait(tree, traits[[method]], method, n_tip)
    trait_matrix <- make_fast_input(trait, method)
    reference_input <- prepare_reference_input(method, tree, trait)
    order <- if ((n_tip + method_index) %% 2L) {
      c("reference", "fastphylosig")
    } else {
      c("fastphylosig", "reference")
    }
    for (j in seq_along(order)) {
      implementation <- order[[j]]
      call_seed <- base_seed + n_tip * 10L + method_index * 2L + j
      fun <- if (implementation == "fastphylosig") {
        function() fast_call(method, context, trait_matrix)
      } else {
        function() reference_call(method, reference_input)
      }
      record_warmup(
        implementation, method, n_tip,
        list(warmup_id = paste0(method, "_n", n_tip), call_seed = call_seed),
        function() {
          set_reproducible_seed(call_seed)
          fun()
        }
      )
    }
  }
}

for (n_tip in n_grid) warm_up(n_tip)

for (n_tip in n_grid) {
  for (replicate in seq_len(formal_reps)) {
    tree_seed <- base_seed + n_tip * 1000L + replicate * 10L
    tree_raw <- make_tree(n_tip, tree_seed)
    context <- fastphylosig::prepare_tree(tree_raw)
    tree <- context$tree
    if (ape::Ntip(tree) != n_tip || length(tree$tip.label) != n_tip) {
      stop("Prepared tree has unexpected tip count.", call. = FALSE)
    }
    tree_hash <- hash_object(tree)
    component_seed <- base_seed + n_tip * 100000L + replicate * 1000L
    bm_seed <- component_seed + 1L
    noise_seed <- component_seed + 2L
    for (scenario_index in seq_len(nrow(scenario_grid))) {
      scenario <- scenario_grid$scenario[[scenario_index]]
      h2 <- scenario_grid$signal_mixture_h2[[scenario_index]]
      traits <- make_scenario_traits(tree, bm_seed, noise_seed, h2)
      bm_hash <- hash_object(traits$bm_component)
      noise_hash <- hash_object(traits$noise_component)
      for (method_index in seq_along(methods)) {
        method <- methods[[method_index]]
        trait <- assert_trait(tree, traits[[method]], method, n_tip)
        trait_matrix <- make_fast_input(trait, method)
        canonical_states <- if (method == "D") as.integer(trait) else
          as.integer(trait)
        trait_hash <- hash_object(list(
          species = names(trait), states = canonical_states,
          state_labels = if (method == "D") c("0", "1") else c("a", "b", "c")
        ))
        reference_input <- prepare_reference_input(method, tree, trait)
        method_params <- if (method == "D") {
          list(test = TRUE, permut = d_nsim)
        } else {
          list(
            test = FALSE, mcmc_sim = delta_mcmc_sim,
            thin = delta_thin, burn = delta_burn, lambda0 = delta_lambda0,
            proposal_sd = delta_proposal_sd, entropy = delta_entropy,
            model = delta_model, chains = 2L
          )
        }
        reference_input_hash <- hash_object(list(
          method = method, tree = tree, trait_hash = trait_hash,
          parameters = method_params
        ))
        comparison_id <- sprintf(
          "%s_%s_n%04d_r%02d", method, scenario, n_tip, replicate
        )
        call_order_vector <- if (
          (replicate + method_index + scenario_index + n_tip) %% 2L
        ) c("reference", "fastphylosig") else c("fastphylosig", "reference")
        first_implementation <- call_order_vector[[1L]]
        meta <- list(
          comparison_id = comparison_id,
          scenario = scenario,
          signal_mixture_h2 = h2,
          n_tip = n_tip,
          replicate = replicate,
          tree_seed = tree_seed,
          component_seed = component_seed,
          bm_seed = bm_seed,
          noise_seed = noise_seed,
          tree_hash = tree_hash,
          bm_hash = bm_hash,
          noise_hash = noise_hash,
          trait_hash = trait_hash,
          reference_input_hash = reference_input_hash
        )
        pair_call_seed <- base_seed + n_tip * 1000L + replicate * 100L +
          scenario_index * 10000L + method_index * 100L
        for (call_order in seq_along(call_order_vector)) {
          implementation <- call_order_vector[[call_order]]
          call_seed <- pair_call_seed +
            if (implementation == "fastphylosig") 1L else 2L
          run_one(
            implementation, method, context, tree, trait, trait_matrix,
            reference_input, meta, call_seed, first_implementation, call_order
          )
        }
      }
    }
  }
}

runs <- utils::read.csv(run_file, stringsAsFactors = FALSE,
                        check.names = FALSE)
pair_ids <- split(seq_len(nrow(runs)), runs$comparison_id)
pair_rows <- lapply(pair_ids, function(index) {
  z <- runs[index, , drop = FALSE]
  fast <- z[z$implementation == "fastphylosig", , drop = FALSE]
  reference <- z[z$implementation == "reference", , drop = FALSE]
  if (nrow(fast) != 1L || nrow(reference) != 1L) {
    stop("Each comparison needs one fast and one reference row.", call. = FALSE)
  }
  out <- fast[1L, c(
    "run_id", "profile", "comparison_id", "first_implementation",
    "method", "reference", "scenario",
    "signal_mixture_h2", "n_tip", "replicate", "tree_seed", "component_seed",
    "bm_seed", "noise_seed", "tree_hash", "bm_hash", "noise_hash",
    "trait_hash", "reference_input_hash", "d_nsim", "delta_mcmc_sim",
    "delta_thin", "delta_burn", "delta_lambda0", "delta_proposal_sd",
    "delta_chains", "delta_expected_total_iterations",
    "delta_reported_requested_iterations",
    "delta_reported_successful_iterations",
    "delta_requested_iterations_match", "delta_successful_iterations_match",
    "delta_entropy", "delta_model", "package_version", "source_tarball_sha256",
    "dll_file_sha256", "reference_version", "test", "timing_clock",
    "timing_iterations"
  ), drop = FALSE]
  out$n_fast_species <- fast$n_species_used
  out$n_reference_species <- reference$n_species_used
  out$fast_reported_delta_control_mcmc_sim <-
    fast$reported_delta_control_mcmc_sim
  out$fast_reported_delta_control_thin <-
    fast$reported_delta_control_thin
  out$fast_reported_delta_control_burn <-
    fast$reported_delta_control_burn
  out$fast_reported_delta_control_lambda0 <-
    fast$reported_delta_control_lambda0
  out$fast_reported_delta_control_proposal_sd <-
    fast$reported_delta_control_proposal_sd
  out$fast_reported_delta_control_entropy <-
    fast$reported_delta_control_entropy
  out$fast_reported_delta_control_model <-
    fast$reported_delta_control_model
  out$fast_reported_delta_control_ace_engine <-
    fast$reported_delta_control_ace_engine
  out$fast_reported_delta_control_ncores <-
    fast$reported_delta_control_ncores
  out$fast_delta_requested_iterations_match <-
    fast$delta_requested_iterations_match
  out$fast_delta_successful_iterations_match <-
    fast$delta_successful_iterations_match
  out$fast_matched_species <- fast$matched_species
  out$reference_matched_species <- reference$matched_species
  out$fast_n_removed_na <- fast$n_removed_na
  out$reference_n_removed_na <- reference$n_removed_na
  out$fast_input_tree_tip_count <- fast$input_tree_tip_count
  out$reference_input_tree_tip_count <- reference$input_tree_tip_count
  out$fast_input_trait_length <- fast$input_trait_length
  out$reference_input_trait_length <- reference$input_trait_length
  out$reference_data_rows <- reference$reference_data_rows
  out$reference_tree_tip_count <- reference$reference_tree_tip_count
  out$fast_elapsed_seconds <- fast$elapsed_seconds
  out$reference_elapsed_seconds <- reference$elapsed_seconds
  out$fast_estimate <- fast$estimate
  out$reference_estimate <- reference$estimate
  out$fast_status <- fast$reported_status
  out$reference_status <- reference$reported_status
  out$fast_status_source <- fast$status_source
  out$reference_status_source <- reference$status_source
  out$fast_success <- fast$success
  out$reference_success <- reference$success
  out$complete_pair <- isTRUE(fast$success) && isTRUE(reference$success)
  out$same_tree_hash <- identical(fast$tree_hash, reference$tree_hash)
  out$same_trait_hash <- identical(fast$trait_hash, reference$trait_hash)
  out$same_reference_input_hash <- identical(
    fast$reference_input_hash, reference$reference_input_hash
  )
  out$species_match <-
    is.finite(out$n_fast_species) && is.finite(out$n_reference_species) &&
    out$n_fast_species == out$n_reference_species &&
    out$n_fast_species == fast$matched_species &&
    out$n_reference_species == reference$matched_species &&
    out$n_fast_species == fast$n_tip &&
    out$n_reference_species == reference$n_tip &&
    fast$n_removed_na == 0L && reference$n_removed_na == 0L &&
    (is.na(fast$input_tree_tip_count) ||
       fast$input_tree_tip_count == fast$n_species_used) &&
    (is.na(reference$input_tree_tip_count) ||
       reference$input_tree_tip_count == reference$n_species_used) &&
    (is.na(reference$input_trait_length) ||
       reference$input_trait_length == reference$n_species_used)
  out$fast_warning <- fast$warning
  out$reference_warning <- reference$warning
  out$fast_error <- fast$error
  out$reference_error <- reference$error
  out$reference_over_fast_ratio <-
    reference$elapsed_seconds / fast$elapsed_seconds
  out
})
pairs <- do.call(rbind, pair_rows)
pairs <- pairs[order(pairs$method, pairs$scenario, pairs$n_tip,
                     pairs$replicate), , drop = FALSE]
rownames(pairs) <- NULL
# Persist the derived pair audit before assertions so any failed audit remains
# inspectable alongside the append-only raw call rows.
utils::write.csv(pairs, output_paths[[2L]], row.names = FALSE, na = "",
                 fileEncoding = "UTF-8")

expected_cells <- expand.grid(
  method = methods, scenario = scenario_grid$scenario, n_tip = n_grid,
  stringsAsFactors = FALSE
)
observed_counts <- as.data.frame(table(
  factor(pairs$method, levels = methods),
  factor(pairs$scenario, levels = scenario_grid$scenario),
  factor(pairs$n_tip, levels = n_grid)
), stringsAsFactors = FALSE)
names(observed_counts) <- c(
  "method", "scenario", "n_tip", "n_pairs_observed"
)
observed_counts$n_tip <- as.integer(as.character(observed_counts$n_tip))
expected_cells$n_pairs_expected <- formal_reps
coverage <- merge(expected_cells, observed_counts,
                  by = c("method", "scenario", "n_tip"), all = TRUE,
                  suffixes = c("_expected", "_observed"))
audit_checks <- c(
  run_count = nrow(runs) == expected_pair_count * 2L,
  pair_count = nrow(pairs) == expected_pair_count,
  run_success = all(runs$success %in% TRUE),
  pair_complete = all(pairs$complete_pair %in% TRUE),
  species_match = all(pairs$species_match %in% TRUE),
  same_tree_hash = all(pairs$same_tree_hash %in% TRUE),
  same_trait_hash = all(pairs$same_trait_hash %in% TRUE),
  same_reference_input_hash = all(pairs$same_reference_input_hash %in% TRUE),
  coverage = all(coverage$n_pairs_observed == coverage$n_pairs_expected)
)
if (anyNA(audit_checks) || any(!audit_checks)) {
  failed_checks <- names(audit_checks)[is.na(audit_checks) | !audit_checks]
  stop(
    "D/Delta matched-input/species/coverage audit failed (",
    paste(failed_checks, collapse = ", "), "); raw rows retained.",
    call. = FALSE
  )
}
if (profile == "full" && !identical(sort(unique(pairs$n_tip)), full_n_grid)) {
  stop("Full profile must contain n=50, 100, 500 and 2000.", call. = FALSE)
}
if (profile == "full") {
  order_counts <- aggregate(
    comparison_id ~ method + scenario + n_tip + first_implementation,
    pairs, length
  )
  if (any(order_counts$comparison_id != formal_reps / 2L)) {
    stop("Full profile call order is not balanced within each cell.",
         call. = FALSE)
  }
}

summary_rows <- lapply(split(
  pairs, interaction(pairs$method, pairs$scenario, pairs$n_tip, drop = TRUE)
), function(z) {
  data.frame(
    package_version = loaded_version,
    reference_version = z$reference_version[[1L]],
    method = z$method[[1L]], scenario = z$scenario[[1L]],
    signal_mixture_h2 = z$signal_mixture_h2[[1L]],
    n_tip = z$n_tip[[1L]], n_pairs = nrow(z),
    fast_median_seconds = stats::median(z$fast_elapsed_seconds),
    fast_q1_seconds = unname(stats::quantile(z$fast_elapsed_seconds, 0.25)),
    fast_q3_seconds = unname(stats::quantile(z$fast_elapsed_seconds, 0.75)),
    reference_median_seconds = stats::median(z$reference_elapsed_seconds),
    reference_q1_seconds = unname(stats::quantile(z$reference_elapsed_seconds, 0.25)),
    reference_q3_seconds = unname(stats::quantile(z$reference_elapsed_seconds, 0.75)),
    reference_over_fast_ratio_of_medians =
      stats::median(z$reference_elapsed_seconds) /
      stats::median(z$fast_elapsed_seconds),
    paired_ratio_median = stats::median(z$reference_over_fast_ratio),
    paired_ratio_q1 = unname(stats::quantile(z$reference_over_fast_ratio, 0.25)),
    paired_ratio_q3 = unname(stats::quantile(z$reference_over_fast_ratio, 0.75)),
    test = z$test[[1L]], timing_clock = "Sys.time", timing_iterations = 1L,
    preparation_included = FALSE,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(summary$method, summary$scenario, summary$n_tip), ]
rownames(summary) <- NULL

script_path <- file.path(out_dir, "benchmark_reference_d_delta.R")
script_sha256 <- digest::digest(
  script_path, algo = "sha256", serialize = FALSE, file = TRUE
)
metadata <- data.frame(
  field = c(
    "run_id", "profile", "package_version", "source_tarball",
    "source_tarball_sha256", "expected_source_tarball_sha256",
    "package_install_library_label", "loaded_package_path_sha256",
    "loaded_dll_path_sha256", "loaded_dll_file_sha256", "R_version",
    "platform", "operating_system", "benchmark_script_sha256",
    "methods", "D_reference", "D_nsim", "D_permutations",
    "Delta_reference_source_sha256", "Delta_reference_source_url",
    "Delta_reference", "Delta_mcmc_sim", "Delta_thin", "Delta_burn",
    "Delta_lambda0", "Delta_proposal_sd", "Delta_chains",
    "Delta_entropy", "Delta_model", "n_grid", "scenario_grid",
    "formal_replicates_per_cell", "expected_pairs", "expected_run_rows",
    "warmup", "call_order", "timing_boundary", "timing_clock",
    "timing_iterations", "species_contract", "statistics",
    "source_data_files"
  ),
  value = c(
    run_id, profile, loaded_version, tarball_relative,
    observed_tarball_sha256, expected_tarball_sha256,
    "FASTPHYLOSIG_INSTALL_LIB", package_path_sha256, dll_path_sha256,
    dll_file_sha256, R.version.string, R.version$platform,
    as.character(Sys.info()[["sysname"]]), script_sha256,
    "D, Delta", "caper::phylo.d(DEstimate)", as.character(d_nsim),
    as.character(d_nsim), delta_source_sha256, delta_source_url,
    paste0(
      "pinned Borges2019 delta() source; ace(model='ARD'); nentropy() LSE-equivalent; ",
      "two emcmc chains"
    ), as.character(delta_mcmc_sim),
    as.character(delta_thin), as.character(delta_burn),
    as.character(delta_lambda0), as.character(delta_proposal_sd), "2",
    delta_entropy, delta_model, paste(n_grid, collapse = ","),
    "none=h2 0; moderate=h2 0.5; strong=h2 0.95",
    as.character(formal_reps), as.character(expected_pair_count),
    as.character(expected_pair_count * 2L),
    "one untimed reference/fast warm-up per method and n_tip; warnings/errors retained",
    if (profile == "full") "alternating; five fast-first and five reference-first per cell" else
      "alternating by replicate/method/scenario/n_tip parity; smoke has one pair per cell",
    paste(
      "one public fast_d/fast_delta or reference call; tree/trait generation,",
      "prepare_tree(), caper::comparative.data(), matching, RNG seeding and gc()",
      "are outside the timed interval; Delta ace/MCMC remain inside each method call"
    ),
    "Sys.time wall-clock", "one call per implementation per matched pair",
    "same prepared random rooted positive-branch tree and complete named trait set",
    "descriptive medians/Q1/Q3; reference median / fast median; no inference",
    paste(output_names[1:5], collapse = "; ")
  ), stringsAsFactors = FALSE
)

utils::write.csv(summary, output_paths[[3L]], row.names = FALSE, na = "",
                 fileEncoding = "UTF-8")
utils::write.csv(metadata, output_paths[[4L]], row.names = FALSE, na = "",
                 fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), output_paths[[6L]], useBytes = TRUE)
message(
  "Completed ", nrow(pairs), " D/Delta matched pairs (", nrow(runs),
  " timed calls) for profile ", profile,
  "; all input/species/status audits passed."
)
