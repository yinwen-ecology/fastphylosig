# Matched single-call benchmark: fastphylosig 0.3.0 vs phytools.
#
# Scope: point-estimate calls only (test = FALSE). For K this is the K
# estimator, not the permutation test. Lambda uses the maximum-likelihood
# estimate (lambda_profile = FALSE). Tree/data matching and preparation are
# completed before timing. Run from the repository root with --vanilla and pass
# the isolated release install library as --library PATH. No package is loaded
# from the user's default library by this script.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("ape", "digest", "phytools")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run this script from the fastphylosig repository root.", call. = FALSE)
}
out_dir <- file.path(repo_root, "docs", "release-0.3.0")
if (!dir.exists(out_dir)) stop("Missing benchmark output directory.", call. = FALSE)

cli_args <- commandArgs(trailingOnly = TRUE)
library_arg <- which(cli_args == "--library")
library_equals_arg <- grep("^--library=", cli_args)
if (length(library_arg) + length(library_equals_arg) != 1L) {
  stop("Usage: Rscript --vanilla benchmark_reference.R --library PATH",
       call. = FALSE)
}
if (length(library_equals_arg)) {
  install_lib_raw <- sub("^--library=", "", cli_args[[library_equals_arg]])
  accepted_args <- cli_args[[library_equals_arg]]
} else {
  if (library_arg == length(cli_args)) {
    stop("--library must be followed by an installation-library path.",
         call. = FALSE)
  }
  install_lib_raw <- cli_args[[library_arg + 1L]]
  if (startsWith(install_lib_raw, "--")) {
    stop("--library must be followed by an installation-library path.",
         call. = FALSE)
  }
  accepted_args <- c("--library", install_lib_raw)
}
unknown_args <- setdiff(cli_args, accepted_args)
if (length(unknown_args)) {
  stop("Unknown benchmark arguments: ", paste(unknown_args, collapse = ", "),
       call. = FALSE)
}
install_lib <- normalizePath(install_lib_raw, winslash = "/", mustWork = TRUE)
expected_package_dir <- normalizePath(
  file.path(install_lib, "fastphylosig"), winslash = "/", mustWork = TRUE
)

tarball_raw <- Sys.getenv(
  "FASTPHYLOSIG_SOURCE_TARBALL",
  unset = file.path(
    repo_root, "benchmarks", "stage5", "documentation_fix_release_20261001",
    "fastphylosig_0.3.0.tar.gz"
  )
)
tarball <- normalizePath(tarball_raw, winslash = "/", mustWork = TRUE)
expected_tarball_sha256 <- "909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff"
provided_sha256 <- Sys.getenv(
  "FASTPHYLOSIG_EXPECTED_TARBALL_SHA256", expected_tarball_sha256
)
if (!identical(tolower(provided_sha256), expected_tarball_sha256)) {
  stop("The expected source tarball SHA-256 does not match the approved value.",
       call. = FALSE)
}
observed_tarball_sha256 <- digest::digest(
  tarball, algo = "sha256", serialize = FALSE, file = TRUE
)
if (!identical(tolower(observed_tarball_sha256), expected_tarball_sha256)) {
  stop("Release tarball SHA-256 mismatch; no benchmark was run.", call. = FALSE)
}

if ("fastphylosig" %in% loadedNamespaces()) {
  stop("fastphylosig is already loaded; rerun in a fresh Rscript --vanilla process.",
       call. = FALSE)
}
suppressPackageStartupMessages(library(
  "fastphylosig", character.only = TRUE, lib.loc = install_lib
))
package_path <- normalizePath(
  find.package("fastphylosig"), winslash = "/", mustWork = TRUE
)
if (!identical(tolower(package_path), tolower(expected_package_dir))) {
  stop("Loaded package path is not the explicitly selected release library.",
       call. = FALSE)
}
loaded_version <- as.character(utils::packageVersion("fastphylosig"))
if (!identical(loaded_version, "0.3.0")) {
  stop("Expected fastphylosig 0.3.0; loaded ", loaded_version, ".",
       call. = FALSE)
}
dll_list <- getLoadedDLLs()
dll_names <- names(dll_list)
dll_index <- which(grepl("^fastphylosig(\\.dll)?$", dll_names, ignore.case = TRUE))
if (length(dll_index) != 1L) {
  stop("Could not identify exactly one loaded fastphylosig DLL.", call. = FALSE)
}
dll_path <- normalizePath(dll_list[[dll_index]][["path"]], winslash = "/",
                          mustWork = TRUE)
if (!startsWith(tolower(dll_path), paste0(tolower(package_path), "/"))) {
  stop("Loaded DLL is outside the selected package installation.", call. = FALSE)
}
if (!all(c("fast_k", "fast_lambda", "prepare_tree") %in%
         getNamespaceExports("fastphylosig"))) {
  stop("The selected 0.3.0 library lacks a required public API.", call. = FALSE)
}
phytools_version <- as.character(utils::packageVersion("phytools"))

# Keep private machine paths out of source-data files that may be published.
# The path labels plus SHA-256 fingerprints attest to the resolved locations.
package_path_label <- "isolated_release_library/fastphylosig"
dll_path_label <- "isolated_release_library/fastphylosig/libs/<platform>/fastphylosig"
package_path_sha256 <- digest::digest(
  package_path, algo = "sha256", serialize = FALSE
)
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
local_path_log <- file.path(
  repo_root, "benchmarks", "stage5",
  paste0("benchmark_reference_local_paths_", run_id, ".txt")
)
writeLines(c(
  "Local-only benchmark execution provenance; do not publish or commit.",
  paste0("run_id=", run_id),
  paste0("source_tarball_path=", tarball),
  paste0("source_tarball_sha256=", observed_tarball_sha256),
  paste0("install_library_path=", install_lib),
  paste0("loaded_package_path=", package_path),
  paste0("loaded_package_version=", loaded_version),
  paste0("loaded_dll_path=", dll_path)
), local_path_log, useBytes = TRUE)

output_names <- c(
  "benchmark_reference_runs.csv",
  "benchmark_reference_pairs.csv",
  "benchmark_reference_summary.csv",
  "benchmark_reference_metadata.csv",
  "benchmark_reference_session_info.txt"
)
existing_outputs <- file.path(out_dir, output_names)
if (any(file.exists(existing_outputs))) {
  stop(
    "Refusing to overwrite existing benchmark outputs: ",
    paste(basename(existing_outputs[file.exists(existing_outputs)]),
          collapse = ", "),
    call. = FALSE
  )
}
run_file <- file.path(out_dir, output_names[[1L]])

n_grid <- c(50L, 100L, 500L, 2000L)
methods <- c("K", "lambda")
scenario_grid <- data.frame(
  scenario = c("none", "moderate", "strong"),
  signal_mixture_h2 = c(0, 0.5, 0.95),
  stringsAsFactors = FALSE
)
formal_reps <- 10L
base_seed <- 20260914L
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

make_components <- function(tree, trait_seed) {
  set_reproducible_seed(trait_seed)
  bm <- ape::rTraitCont(tree, model = "BM", ancestor = FALSE)
  bm <- standardize_component(bm, "Brownian")
  names(bm) <- tree$tip.label
  bm <- bm[tree$tip.label]
  noise <- standardize_component(stats::rnorm(ape::Ntip(tree)), "iid noise")
  names(noise) <- tree$tip.label
  noise <- noise[tree$tip.label]
  list(bm = bm, noise = noise)
}

make_scenario_trait <- function(tree, components, h2) {
  trait <- sqrt(h2) * components$bm + sqrt(1 - h2) * components$noise
  names(trait) <- tree$tip.label
  trait[tree$tip.label]
}

assert_matched_trait <- function(tree, trait, expected_n, label) {
  if (!inherits(tree, "phylo") || ape::Ntip(tree) != expected_n ||
      length(tree$tip.label) != expected_n) {
    stop(label, ": unexpected tree tip count.", call. = FALSE)
  }
  if (length(trait) != expected_n || is.null(names(trait)) || anyNA(trait) ||
      any(!is.finite(trait)) || !setequal(names(trait), tree$tip.label)) {
    stop(label, ": trait does not match the complete tree species set.",
         call. = FALSE)
  }
  trait <- as.numeric(trait[tree$tip.label])
  names(trait) <- tree$tip.label
  if (!identical(names(trait), tree$tip.label)) {
    stop(label, ": trait order differs from tree tip order.", call. = FALSE)
  }
  trait
}

hash_object <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

fast_call <- function(method, context, trait) {
  # Use a pre-built one-column matrix so the public result exposes actual
  # species-audit counts without charging input conversion to call time.
  if (method == "K") {
    return(fastphylosig::fast_k(
      tree = context$tree, x = trait, prepared = context, test = FALSE,
      verbose = FALSE, ncores = 1L, progress = FALSE
    ))
  }
  if (method == "lambda") {
    return(fastphylosig::fast_lambda(
      tree = context$tree, x = trait, prepared = context, test = FALSE,
      verbose = FALSE, ncores = 1L, lambda_profile = FALSE, progress = FALSE
    ))
  }
  stop("Unsupported method: ", method, call. = FALSE)
}

make_trait_matrix <- function(trait) {
  matrix(
    as.numeric(trait), ncol = 1L,
    dimnames = list(names(trait), "trait")
  )
}

reference_call <- function(method, tree, trait) {
  phytools::phylosig(
    tree, trait, method = method, test = FALSE, se = NULL
  )
}

extract_estimate <- function(value, method, implementation) {
  if (is.null(value)) return(NA_real_)
  if (implementation == "reference" && is.numeric(value)) {
    return(suppressWarnings(as.numeric(value[[1L]])))
  }
  keys <- if (method == "K") c("K_fast", "K", "estimate") else
    c("lambda_fast", "lambda", "estimate")
  for (key in keys) {
    if (is.data.frame(value) && key %in% names(value)) {
      return(suppressWarnings(as.numeric(value[[key]][[1L]])))
    }
    if (is.list(value) && !is.null(value[[key]])) {
      return(suppressWarnings(as.numeric(value[[key]][[1L]])))
    }
  }
  if (is.numeric(value) && length(value)) {
    return(suppressWarnings(as.numeric(value[[1L]])))
  }
  NA_real_
}

extract_status <- function(value, implementation, estimate) {
  if (implementation == "reference") {
    return(if (length(estimate) == 1L && is.finite(estimate)) "ok" else "invalid_estimate")
  }
  if (is.list(value) && !is.null(value[["status"]])) {
    return(as.character(value[["status"]][[1L]]))
  }
  if (is.data.frame(value) && "status" %in% names(value) &&
      nrow(value) == 1L) {
    return(as.character(value$status[[1L]]))
  }
  if (length(estimate) == 1L && is.finite(estimate)) "not_reported" else
    "invalid_estimate"
}

audit_call_species <- function(value, implementation, tree, trait) {
  if (implementation == "fastphylosig") {
    required <- c("n_species", "matched_species", "n_removed_na")
    if (!is.data.frame(value) || nrow(value) != 1L ||
        !all(required %in% names(value))) {
      return(list(
        n_species_used = NA_integer_, matched_species = NA_integer_,
        n_removed_na = NA_integer_, tree_tip_count = NA_integer_,
        trait_length = NA_integer_, validation_error =
          "fastphylosig result lacks actual species-audit fields"
      ))
    }
    n_species_used <- suppressWarnings(as.integer(value$n_species[[1L]]))
    matched_species <- suppressWarnings(as.integer(value$matched_species[[1L]]))
    n_removed_na <- suppressWarnings(as.integer(value$n_removed_na[[1L]]))
    valid <- all(is.finite(c(n_species_used, matched_species, n_removed_na))) &&
      n_species_used >= 0L && matched_species >= 0L && n_removed_na >= 0L
    return(list(
      n_species_used = n_species_used,
      matched_species = matched_species,
      n_removed_na = n_removed_na,
      tree_tip_count = NA_integer_, trait_length = NA_integer_,
      validation_error = if (valid) NA_character_ else
        "fastphylosig returned invalid species-audit counts"
    ))
  }

  tree_tip_count <- as.integer(ape::Ntip(tree))
  trait_length <- as.integer(length(trait))
  trait_names <- names(trait)
  matched_species <- if (is.null(trait_names)) NA_integer_ else
    as.integer(sum(tree$tip.label %in% trait_names))
  n_removed_na <- as.integer(sum(is.na(trait)))
  valid <- inherits(tree, "phylo") &&
    length(tree$tip.label) == tree_tip_count &&
    trait_length == tree_tip_count && !is.null(trait_names) &&
    !anyNA(trait) && all(is.finite(trait)) &&
    identical(trait_names, tree$tip.label) &&
    matched_species == tree_tip_count && n_removed_na == 0L
  list(
    n_species_used = tree_tip_count,
    matched_species = matched_species,
    n_removed_na = n_removed_na,
    tree_tip_count = tree_tip_count,
    trait_length = trait_length,
    validation_error = if (valid) NA_character_ else
      "reference tree tips and complete named trait vector do not match"
  )
}

timed_call <- function(fun) {
  warning_messages <- character()
  error_message <- NA_character_
  value <- NULL
  started <- Sys.time()
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
  stopped <- Sys.time()
  list(
    value = value,
    elapsed_seconds = as.numeric(difftime(stopped, started, units = "secs")),
    warning = paste(unique(warning_messages), collapse = " || "),
    error = error_message
  )
}

append_run <- function(row) {
  utils::write.table(
    row, file = run_file, sep = ",", row.names = FALSE,
    col.names = !file.exists(run_file), append = file.exists(run_file),
    quote = TRUE, qmethod = "double", na = "", fileEncoding = "UTF-8"
  )
}

run_one <- function(implementation, method, context, trait, fast_trait_matrix,
                    tree, meta,
                    call_seed, first_implementation, call_order) {
  # State setup and garbage collection are outside the timed interval.
  set_reproducible_seed(call_seed)
  invisible(gc())
  timed <- if (implementation == "fastphylosig") {
    timed_call(function() fast_call(method, context, fast_trait_matrix))
  } else {
    timed_call(function() reference_call(method, tree, trait))
  }
  estimate <- extract_estimate(timed$value, method, implementation)
  status <- extract_status(timed$value, implementation, estimate)
  species_audit <- audit_call_species(timed$value, implementation, tree, trait)
  success <- is.na(timed$error) && length(estimate) == 1L &&
    is.finite(estimate) && identical(status, "ok") &&
    is.na(species_audit$validation_error)
  row <- data.frame(
    run_id = run_id,
    comparison_id = meta$comparison_id,
    implementation = implementation,
    method = method,
    reference = "phytools::phylosig",
    scenario = meta$scenario,
    signal_mixture_h2 = meta$signal_mixture_h2,
    n_tip = meta$n_tip,
    replicate = meta$replicate,
    n_species_used = species_audit$n_species_used,
    matched_species = species_audit$matched_species,
    n_removed_na = species_audit$n_removed_na,
    input_tree_tip_count = species_audit$tree_tip_count,
    input_trait_length = species_audit$trait_length,
    tree_seed = meta$tree_seed,
    trait_seed = meta$trait_seed,
    call_seed = as.integer(call_seed),
    first_implementation = first_implementation,
    call_order = as.integer(call_order),
    test = FALSE,
    lambda_profile = if (method == "lambda") FALSE else NA,
    prepare_excluded = TRUE,
    timing_clock = "Sys.time",
    timing_iterations = 1L,
    tree_hash = meta$tree_hash,
    bm_hash = meta$bm_hash,
    noise_hash = meta$noise_hash,
    trait_hash = meta$trait_hash,
    reference_input_hash = meta$reference_input_hash,
    prepared_fingerprint = meta$prepared_fingerprint,
    package_version = loaded_version,
    package_path = package_path_label,
    package_path_sha256 = package_path_sha256,
    dll_path = dll_path_label,
    dll_path_sha256 = dll_path_sha256,
    dll_file_sha256 = dll_file_sha256,
    source_tarball = tarball_relative,
    source_tarball_sha256 = observed_tarball_sha256,
    reference_version = phytools_version,
    elapsed_seconds = timed$elapsed_seconds,
    estimate = estimate,
    reported_status = status,
    success = success,
    validation_error = species_audit$validation_error,
    warning = timed$warning,
    error = timed$error,
    stringsAsFactors = FALSE
  )
  # Persist every call immediately; if one fails, keep its row and stop.
  append_run(row)
  if (!is.finite(row$elapsed_seconds) || row$elapsed_seconds <= 0 || !success) {
    reason <- if (!is.na(timed$error)) timed$error else
      if (!is.finite(row$elapsed_seconds) || row$elapsed_seconds <= 0)
        "non-positive or non-finite elapsed time" else
        if (!identical(status, "ok")) paste0("reported status is ", status) else
        if (!is.na(species_audit$validation_error))
          species_audit$validation_error else "non-finite estimate"
    stop(
      "Benchmark call failed; its raw row is retained in benchmark_reference_runs.csv. ",
      meta$comparison_id, " / ", implementation, ": ", reason,
      call. = FALSE
    )
  }
  invisible(row)
}

run_file <- file.path(out_dir, "benchmark_reference_runs.csv")

warm_up <- function(n_tip) {
  warm_seed <- base_seed + n_tip + 1L
  warm_tree_raw <- make_tree(n_tip, warm_seed)
  warm_context <- fastphylosig::prepare_tree(warm_tree_raw)
  warm_tree <- warm_context$tree
  warm_components <- make_components(warm_tree, base_seed + n_tip + 1L)
  warm_trait <- assert_matched_trait(
    warm_tree, make_scenario_trait(warm_tree, warm_components, 0.5),
    n_tip, "warm-up trait"
  )
  warm_trait_matrix <- make_trait_matrix(warm_trait)
  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    set_reproducible_seed(base_seed + n_tip * 10L + method_index)
    fast_value <- tryCatch(fast_call(method, warm_context, warm_trait_matrix),
                           error = identity)
    if (inherits(fast_value, "error")) {
      stop("Untimed fastphylosig warm-up failed for ", method, "/n", n_tip,
           ": ", conditionMessage(fast_value), call. = FALSE)
    }
    invisible(gc())
    reference_value <- tryCatch(
      reference_call(method, warm_tree, warm_trait), error = identity
    )
    if (inherits(reference_value, "error")) {
      stop("Untimed phytools warm-up failed for ", method, "/n", n_tip,
           ": ", conditionMessage(reference_value), call. = FALSE)
    }
    invisible(gc())
  }
}

for (n_tip in n_grid) {
  warm_up(n_tip)
  for (replicate in seq_len(formal_reps)) {
    tree_seed <- base_seed + n_tip * 1000L + replicate * 10L
    tree_raw <- make_tree(n_tip, tree_seed)
    context <- fastphylosig::prepare_tree(tree_raw)
    tree <- context$tree
    tree_hash <- hash_object(tree)
    context_fingerprint <- if (!is.null(context$fingerprint)) {
      hash_object(context$fingerprint)
    } else {
      NA_character_
    }
    trait_seed <- base_seed + n_tip * 100000L + replicate * 1000L
    components <- make_components(tree, trait_seed)
    bm_hash <- hash_object(components$bm)
    noise_hash <- hash_object(components$noise)

    for (scenario_index in seq_len(nrow(scenario_grid))) {
      scenario <- scenario_grid$scenario[[scenario_index]]
      h2 <- scenario_grid$signal_mixture_h2[[scenario_index]]
      trait <- assert_matched_trait(
        tree, make_scenario_trait(tree, components, h2), n_tip,
        paste0("trait/", scenario, "/n", n_tip, "/r", replicate)
      )
      trait_hash <- hash_object(trait)
      fast_trait_matrix <- make_trait_matrix(trait)
      reference_input_hash <- hash_object(list(tree = tree, trait = trait))

      for (method_index in seq_along(methods)) {
        method <- methods[[method_index]]
        comparison_id <- sprintf("%s_%s_n%04d_r%02d", method, scenario,
                                 n_tip, replicate)
        call_order_vector <- if (
          (replicate + method_index + scenario_index + n_tip) %% 2L
        ) c("reference", "fastphylosig") else c("fastphylosig", "reference")
        meta <- list(
          comparison_id = comparison_id,
          scenario = scenario,
          signal_mixture_h2 = h2,
          n_tip = n_tip,
          replicate = replicate,
          tree_seed = tree_seed,
          trait_seed = trait_seed,
          tree_hash = tree_hash,
          bm_hash = bm_hash,
          noise_hash = noise_hash,
          trait_hash = trait_hash,
          reference_input_hash = reference_input_hash,
          prepared_fingerprint = context_fingerprint
        )
        pair_call_seed <- base_seed + n_tip * 1000L + replicate * 100L +
          scenario_index * 10L + method_index * 2L
        for (call_order in seq_along(call_order_vector)) {
          run_one(
            implementation = call_order_vector[[call_order]],
            method = method,
            context = context,
            trait = trait,
            fast_trait_matrix = fast_trait_matrix,
            tree = tree,
            meta = meta,
            call_seed = pair_call_seed,
            first_implementation = call_order_vector[[1L]],
            call_order = call_order
          )
        }
      }
    }
  }
}

runs <- utils::read.csv(run_file, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(runs) != expected_pair_count * 2L || any(!runs$success) ||
    any(runs$test != FALSE) || any(!runs$prepare_excluded) ||
    any(runs$timing_iterations != 1L)) {
  stop("Raw-run completeness or timing-contract validation failed.", call. = FALSE)
}

pair_ids <- split(seq_len(nrow(runs)), runs$comparison_id)
pair_rows <- lapply(pair_ids, function(index) {
  z <- runs[index, , drop = FALSE]
  fast <- z[z$implementation == "fastphylosig", , drop = FALSE]
  reference <- z[z$implementation == "reference", , drop = FALSE]
  if (nrow(fast) != 1L || nrow(reference) != 1L) {
    stop("Each comparison must have one fast and one reference row.",
         call. = FALSE)
  }
  pair <- fast[1L, c(
    "run_id", "comparison_id", "method", "reference", "scenario",
    "signal_mixture_h2", "n_tip", "replicate", "tree_seed", "trait_seed",
    "call_seed", "first_implementation", "tree_hash", "bm_hash",
    "noise_hash", "trait_hash", "reference_input_hash",
    "prepared_fingerprint", "n_species_used", "matched_species",
    "n_removed_na", "package_version", "package_path",
    "package_path_sha256", "dll_path", "dll_path_sha256", "dll_file_sha256",
    "source_tarball", "source_tarball_sha256", "reference_version", "test",
    "lambda_profile", "prepare_excluded", "timing_clock",
    "timing_iterations"
  ), drop = FALSE]
  pair$n_fast_species <- fast$n_species_used
  pair$n_reference_species <- reference$n_species_used
  pair$matched_n_species <- fast$matched_species
  pair$fast_n_removed_na <- fast$n_removed_na
  pair$reference_n_removed_na <- reference$n_removed_na
  pair$fast_elapsed_seconds <- fast$elapsed_seconds
  pair$reference_elapsed_seconds <- reference$elapsed_seconds
  pair$fast_estimate <- fast$estimate
  pair$reference_estimate <- reference$estimate
  pair$fast_status <- fast$reported_status
  pair$reference_status <- reference$reported_status
  pair$fast_success <- fast$success
  pair$reference_success <- reference$success
  pair$fast_warning <- fast$warning
  pair$reference_warning <- reference$warning
  pair$fast_error <- fast$error
  pair$reference_error <- reference$error
  pair$same_tree_hash <- identical(fast$tree_hash, reference$tree_hash)
  pair$same_trait_hash <- identical(fast$trait_hash, reference$trait_hash)
  pair$same_reference_input_hash <- identical(
    fast$reference_input_hash, reference$reference_input_hash
  )
  pair$species_match <-
    fast$n_species_used == reference$n_species_used &&
    fast$n_species_used == fast$matched_species &&
    fast$n_removed_na == 0L && reference$n_removed_na == 0L &&
    reference$input_tree_tip_count == reference$n_species_used &&
    reference$input_trait_length == reference$n_species_used
  pair$complete_pair <- isTRUE(fast$success) && isTRUE(reference$success)
  pair$paired_reference_over_fast <-
    reference$elapsed_seconds / fast$elapsed_seconds
  pair
})
pairs <- do.call(rbind, pair_rows)
pairs <- pairs[order(pairs$method, pairs$scenario, pairs$n_tip,
                     pairs$replicate), , drop = FALSE]
rownames(pairs) <- NULL

expected_cells <- expand.grid(
  method = methods, scenario = scenario_grid$scenario, n_tip = n_grid,
  stringsAsFactors = FALSE
)
cell_key <- function(x) paste(x$method, x$scenario, x$n_tip, sep = "|")
observed_counts <- as.data.frame(table(
  factor(pairs$method, levels = methods),
  factor(pairs$scenario, levels = scenario_grid$scenario),
  factor(pairs$n_tip, levels = n_grid)
), stringsAsFactors = FALSE)
names(observed_counts) <- c("method", "scenario", "n_tip", "n_pairs")
observed_counts$n_tip <- as.integer(as.character(observed_counts$n_tip))
expected_cells$n_pairs <- formal_reps
coverage <- merge(expected_cells, observed_counts,
                  by = c("method", "scenario", "n_tip"),
                  all = TRUE, suffixes = c("_expected", "_observed"))
if (nrow(pairs) != expected_pair_count || any(!pairs$complete_pair) ||
    any(!pairs$species_match) || any(!pairs$same_tree_hash) ||
    any(!pairs$same_trait_hash) || any(!pairs$same_reference_input_hash) ||
    any(pairs$n_fast_species != pairs$n_reference_species) ||
    any(pairs$n_fast_species != pairs$n_tip) ||
    any(pairs$fast_n_removed_na != 0L) ||
    any(pairs$reference_n_removed_na != 0L) ||
    any(coverage$n_pairs_observed != coverage$n_pairs_expected)) {
  stop("Matched-input pair/coverage audit failed; inspect retained raw CSV.",
       call. = FALSE)
}
order_counts <- aggregate(
  comparison_id ~ method + scenario + n_tip + first_implementation,
  pairs, length
)
if (any(order_counts$comparison_id != formal_reps / 2L)) {
  stop("Call order is not balanced within every cell.", call. = FALSE)
}

summary_rows <- lapply(split(
  pairs, interaction(pairs$method, pairs$scenario, pairs$n_tip, drop = TRUE)
), function(z) {
  data.frame(
    package_version = loaded_version,
    reference_version = phytools_version,
    method = z$method[[1L]],
    scenario = z$scenario[[1L]],
    signal_mixture_h2 = z$signal_mixture_h2[[1L]],
    n_tip = z$n_tip[[1L]],
    n_pairs = nrow(z),
    fast_median_seconds = stats::median(z$fast_elapsed_seconds),
    fast_q1_seconds = unname(stats::quantile(z$fast_elapsed_seconds, 0.25)),
    fast_q3_seconds = unname(stats::quantile(z$fast_elapsed_seconds, 0.75)),
    reference_median_seconds = stats::median(z$reference_elapsed_seconds),
    reference_q1_seconds = unname(stats::quantile(z$reference_elapsed_seconds, 0.25)),
    reference_q3_seconds = unname(stats::quantile(z$reference_elapsed_seconds, 0.75)),
    reference_over_fast_ratio_of_medians =
      stats::median(z$reference_elapsed_seconds) /
      stats::median(z$fast_elapsed_seconds),
    paired_ratio_median = stats::median(z$paired_reference_over_fast),
    paired_ratio_q1 = unname(stats::quantile(z$paired_reference_over_fast, 0.25)),
    paired_ratio_q3 = unname(stats::quantile(z$paired_reference_over_fast, 0.75)),
    test = FALSE,
    timing_clock = "Sys.time",
    timing_iterations = 1L,
    preparation_included = FALSE,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(summary$method, summary$scenario, summary$n_tip), ]
rownames(summary) <- NULL

benchmark_script <- file.path(out_dir, "benchmark_reference.R")
plot_script <- file.path(out_dir, "plot_benchmark_reference.R")
script_sha <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE)
}
metadata <- data.frame(
  field = c(
    "run_id", "package_version", "source_tarball", "source_tarball_sha256",
    "source_tarball_expected_sha256", "package_install_library_label",
    "loaded_package_path_label", "loaded_package_path_sha256",
    "loaded_dll_path_label", "loaded_dll_path_sha256",
    "loaded_dll_file_sha256", "phytools_version",
    "R_version", "platform", "operating_system", "benchmark_script_sha256",
    "plot_script_sha256", "methods", "reference_methods", "n_grid",
    "scenario_grid", "formal_replicates_per_cell", "expected_pairs",
    "expected_run_rows", "warm_up", "call_order", "test", "estimator_scope",
    "timing_boundary", "timing_clock", "timing_iterations", "species_contract",
    "statistics", "source_data_files"
  ),
  value = c(
    run_id, loaded_version, tarball_relative, observed_tarball_sha256,
    expected_tarball_sha256, "FASTPHYLOSIG_INSTALL_LIB",
    package_path_label, package_path_sha256, dll_path_label, dll_path_sha256,
    dll_file_sha256,
    phytools_version, R.version.string, R.version$platform,
    as.character(Sys.info()[["sysname"]]), script_sha(benchmark_script),
    script_sha(plot_script), "K, lambda", "phytools::phylosig(method=K/lambda)",
    paste(n_grid, collapse = ","),
    "none=h2 0; moderate=h2 0.5; strong=h2 0.95",
    as.character(formal_reps), as.character(expected_pair_count),
    as.character(expected_pair_count * 2L),
    "one untimed fast/reference warm-up per method and n_tip",
    "alternating; 5 fast-first and 5 reference-first pairs per cell",
    "FALSE for both implementations; K point estimator only, no permutation test",
    "K point estimate and lambda maximum-likelihood estimate",
    paste(
      "one public fast_* or phytools::phylosig call; tree/data generation,",
      "matching, prepare_tree(), reference input preparation, gc(), and RNG",
      "seeding are outside the timed interval"
    ),
    "Sys.time wall clock", "1 per implementation per matched pair",
    "same canonical tree, trait values/order, and complete species set",
    paste(
      "descriptive medians and Q1/Q3; ratio of medians is reference median /",
      "fast median; paired call ratios reported separately; no inferential test"
    ),
    paste(output_names[[1L]], "benchmark_reference_pairs.csv",
          "benchmark_reference_summary.csv", "benchmark_reference_metadata.csv",
          sep = "; ")
  ), stringsAsFactors = FALSE
)

utils::write.csv(pairs, file.path(out_dir, output_names[[2L]]), row.names = FALSE,
                 na = "", fileEncoding = "UTF-8")
utils::write.csv(summary, file.path(out_dir, output_names[[3L]]), row.names = FALSE,
                 na = "", fileEncoding = "UTF-8")
utils::write.csv(metadata, file.path(out_dir, output_names[[4L]]), row.names = FALSE,
                 na = "", fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, output_names[[5L]]), useBytes = TRUE)

message(
  "Completed ", nrow(pairs), " matched pairs (", nrow(runs),
  " timed calls); all 24 cells contain ", formal_reps,
  " pairs with identical species/tree/trait hashes."
)
