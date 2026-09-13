#!/usr/bin/env Rscript

# Fresh-process helper for the final memory characterization.  The parent
# runner supplies an RDS path, a CSV output path, and (optionally) the source
# repository when no private installed library is selected.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: memory_validate.R <rds> <output.csv> [repo]", call. = FALSE)
}

rds_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out_path <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
repo <- if (length(args) >= 3L) normalizePath(args[[3L]], winslash = "/",
                                               mustWork = FALSE) else getwd()

library_dir <- Sys.getenv("FASTPHYLOSIG_FINAL_LIBRARY", unset = "")
if (nzchar(library_dir) && dir.exists(file.path(library_dir, "fastphylosig"))) {
  .libPaths(c(normalizePath(library_dir, winslash = "/", mustWork = TRUE),
              .libPaths()))
  suppressPackageStartupMessages(library(fastphylosig))
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
} else {
  stop("fastphylosig is unavailable in the fresh validation process.",
       call. = FALSE)
}

started <- proc.time()[["elapsed"]]
status <- "PASS"
message_text <- ""
warning_text <- character()
value <- tryCatch(
  withCallingHandlers({
    ctx <- readRDS(rds_path)
    validated <- prepare_tree(ctx)
    if (!inherits(validated, "fastphylosig_tree")) {
      stop("readRDS object did not validate as fastphylosig_tree.",
           call. = FALSE)
    }
    invisible(validated)
    TRUE
  }, warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }),
  error = function(e) {
    status <<- "FAIL"
    message_text <<- conditionMessage(e)
    FALSE
  }
)
elapsed <- proc.time()[["elapsed"]] - started

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(
    status = status,
    elapsed_s = as.numeric(elapsed),
    warning_text = paste(unique(warning_text), collapse = " | "),
    error_message = message_text,
    r_version = R.version.string,
    package_version = if (requireNamespace("fastphylosig", quietly = TRUE))
      as.character(utils::packageVersion("fastphylosig")) else NA_character_,
    stringsAsFactors = FALSE
  ),
  out_path,
  row.names = FALSE,
  na = "NA"
)

if (!isTRUE(value) || !identical(status, "PASS")) quit(status = 1L)
