.v2_persistence_child_source <- function() {
  c(
    "mode <- Sys.getenv('FASTPHYLOSIG_V2_MODE')",
    "rds_path <- Sys.getenv('FASTPHYLOSIG_V2_RDS')",
    "pkg_path <- Sys.getenv('FASTPHYLOSIG_V2_PACKAGE')",
    "pkg_lib <- Sys.getenv('FASTPHYLOSIG_V2_LIBRARY')",
    "if (nzchar(pkg_lib)) .libPaths(unique(c(pkg_lib, .libPaths())))",
    "source_pkg <- file.exists(file.path(pkg_path, 'R', 'prepare_tree.R'))",
    "if (source_pkg) {",
    "  if (!requireNamespace('pkgload', quietly = TRUE)) stop('pkgload is required to load the source package in the child session')",
    "  pkgload::load_all(pkg_path, quiet = TRUE)",
    "} else {",
    "  suppressPackageStartupMessages(library('fastphylosig'))",
    "}",
    "if (mode == 'prepare') {",
    "  tr <- ape::rtree(4L)",
    "  ctx <- prepare_tree(tr)",
    "  stopifnot(identical(ctx$context_schema_version, 2L))",
    "  stopifnot(identical(ctx$canonical_contract_version, 2L))",
    "  stopifnot(identical(ctx$protected_snapshot_version, 2L))",
    "  stopifnot(is.raw(ctx$protected_snapshot), is.raw(ctx$fingerprint))",
    "  saveRDS(ctx, rds_path, version = 2L)",
    "  cat('PREPARE_OK\\n')",
    "} else if (mode == 'reuse') {",
    "  ctx <- readRDS(rds_path)",
    "  stopifnot(identical(ctx$context_schema_version, 2L))",
    "  stopifnot(identical(ctx$canonical_contract_version, 2L))",
    "  stopifnot(identical(ctx$protected_snapshot_version, 2L))",
    "  stopifnot(is.raw(ctx$protected_snapshot), is.raw(ctx$fingerprint))",
    "  x <- stats::setNames(seq_len(ctx$n_tip), ctx$tree$tip.label)",
    "  fit_k <- fast_k(ctx, x, test = FALSE, verbose = FALSE, progress = FALSE)",
    "  fit_lambda <- fast_lambda(ctx, x, test = FALSE, verbose = FALSE, progress = FALSE)",
    "  stopifnot(is.finite(as.numeric(fit_k)), is.finite(fit_lambda$lambda))",
    "  cat('REUSE_OK\\n')",
    "} else if (mode == 'mutate') {",
    "  ctx <- readRDS(rds_path)",
    "  x <- stats::setNames(seq_len(ctx$n_tip), ctx$tree$tip.label)",
    "  hits_before <- ctx$cache_meta$hits",
    "  misses_before <- ctx$cache_meta$misses",
    "  ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 1",
    "  edge_error <- tryCatch({ fast_k(ctx, x, test = FALSE, verbose = FALSE, progress = FALSE); '__NO_ERROR__' }, error = function(e) conditionMessage(e))",
    "  stopifnot(grepl('prepare_tree', edge_error, fixed = TRUE))",
    "  stopifnot(identical(ctx$cache_meta$hits, hits_before), identical(ctx$cache_meta$misses, misses_before))",
    "  ctx <- readRDS(rds_path)",
    "  ctx$tree$tip.label <- c(paste0('a', rawToChar(as.raw(13)), 'b'), 'c', 'a', paste0('b', rawToChar(as.raw(13)), 'c'))",
    "  x <- stats::setNames(seq_len(ctx$n_tip), ctx$tree$tip.label)",
    "  hits_before <- ctx$cache_meta$hits",
    "  misses_before <- ctx$cache_meta$misses",
    "  label_error <- tryCatch({ fast_lambda(ctx, x, test = FALSE, verbose = FALSE, progress = FALSE); '__NO_ERROR__' }, error = function(e) conditionMessage(e))",
    "  stopifnot(grepl('prepare_tree', label_error, fixed = TRUE))",
    "  stopifnot(identical(ctx$cache_meta$hits, hits_before), identical(ctx$cache_meta$misses, misses_before))",
    "  cat('MUTATION_OK\\n')",
    "} else {",
    "  stop('unknown persistence child mode')",
    "}"
  )
}

.run_v2_persistence_child <- function(mode, rds_path, pkg_path, pkg_lib,
                                      script_path) {
  rscript <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows" && !file.exists(rscript)) {
    rscript <- paste0(rscript, ".exe")
  }
  env <- c(
    FASTPHYLOSIG_V2_MODE = mode,
    FASTPHYLOSIG_V2_RDS = rds_path,
    FASTPHYLOSIG_V2_PACKAGE = pkg_path,
    FASTPHYLOSIG_V2_LIBRARY = pkg_lib,
    R_LIBS_USER = paste(unique(c(pkg_lib, .libPaths())),
                        collapse = .Platform$path.sep)
  )
  old_env <- Sys.getenv(names(env), unset = NA_character_)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    for (i in seq_along(env)) {
      if (is.na(old_env[[i]])) {
        Sys.unsetenv(names(env)[[i]])
      } else {
        do.call(Sys.setenv, setNames(list(old_env[[i]]), names(env)[[i]]))
      }
    }
  }, add = TRUE)
  output <- system2(
    rscript,
    args = c("--vanilla", script_path),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(
    status = as.integer(status),
    output = output,
    text = paste(output, collapse = "\n")
  )
}

test_that("V2 prepared contexts persist and validate across R sessions", {
  skip_if_not_installed("ape")

  pkg_path <- system.file(package = "fastphylosig")
  if (!nzchar(pkg_path) || !file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    pkg_path <- getNamespaceInfo(asNamespace("fastphylosig"), "path")
  }
  pkg_path <- normalizePath(pkg_path, winslash = "/", mustWork = TRUE)
  pkg_lib <- dirname(pkg_path)
  rds_path <- tempfile("fastphylosig-v2-context-", fileext = ".rds")
  script_path <- tempfile("fastphylosig-v2-child-", fileext = ".R")
  writeLines(.v2_persistence_child_source(), script_path, useBytes = TRUE)
  on.exit(unlink(c(rds_path, script_path), force = TRUE), add = TRUE)

  prepared <- .run_v2_persistence_child(
    "prepare", rds_path, pkg_path, pkg_lib, script_path
  )
  expect_identical(prepared$status, 0L, info = prepared$text)
  expect_match(prepared$text, "PREPARE_OK", fixed = TRUE)

  reused <- .run_v2_persistence_child(
    "reuse", rds_path, pkg_path, pkg_lib, script_path
  )
  expect_identical(reused$status, 0L, info = reused$text)
  expect_match(reused$text, "REUSE_OK", fixed = TRUE)

  mutated <- .run_v2_persistence_child(
    "mutate", rds_path, pkg_path, pkg_lib, script_path
  )
  expect_identical(mutated$status, 0L, info = mutated$text)
  expect_match(mutated$text, "MUTATION_OK", fixed = TRUE)
})
