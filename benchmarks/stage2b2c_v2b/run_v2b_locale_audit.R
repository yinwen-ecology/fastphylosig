#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: run_v2b_locale_audit.R LIBRARY OUTPUT_CSV", call. = FALSE)
}

lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(fastphylosig, lib.loc = lib))

tree <- list(
  edge = matrix(c(
    7, 1, 7, 8, 8, 2, 8, 9, 9, 3,
    9, 10, 10, 4, 10, 5, 7, 6
  ), byrow = TRUE, ncol = 2L),
  edge.length = c(0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99),
  tip.label = c("a", "A", "a.1", "a10", "\u00e4", "10"),
  Nnode = 4L
)
class(tree) <- "phylo"

ns <- asNamespace("fastphylosig")
canonicalize <- get(".safe_canonicalize_core", envir = ns)
fingerprint <- get(".tree_fingerprint", envir = ns)

encode_result <- function() {
  canonical <- canonicalize(tree)
  list(
    edge = canonical$edge,
    edge.length = canonical$edge.length,
    tip.label = canonical$tip.label,
    Nnode = canonical$Nnode,
    fingerprint = fingerprint(canonical, .canonical = TRUE)
  )
}

original <- Sys.getlocale("LC_COLLATE")
on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", original)), add = TRUE)
candidates <- c(
  current = original,
  C = "C",
  English_US_CP1252 = "English_United States.1252",
  German_CP1252 = "German_Germany.1252",
  Chinese_Simplified_CP936 = "Chinese (Simplified)_China.936",
  en_US_UTF8 = "en_US.UTF-8"
)

baseline <- NULL
rows <- lapply(names(candidates), function(label) {
  requested <- candidates[[label]]
  actual <- suppressWarnings(Sys.setlocale("LC_COLLATE", requested))
  available <- is.character(actual) && length(actual) == 1L && nzchar(actual)
  value <- if (available) tryCatch(encode_result(), error = identity) else NULL
  if (identical(label, "C") && is.list(value)) baseline <<- value
  data.frame(
    locale_case = label,
    requested = requested,
    actual = if (available) actual else "",
    status = if (!available) {
      "NOT_RUN_ENVIRONMENT"
    } else if (inherits(value, "error")) {
      "ERROR"
    } else {
      "RUN"
    },
    identical_to_C = NA,
    message = if (inherits(value, "error")) conditionMessage(value) else "",
    stringsAsFactors = FALSE
  )
})

# Re-run available locales after the C baseline is known, so output equality is
# based on exact production objects rather than a character summary.
for (i in seq_along(rows)) {
  if (!identical(rows[[i]]$status, "RUN") || is.null(baseline)) next
  actual <- suppressWarnings(Sys.setlocale("LC_COLLATE", rows[[i]]$requested))
  if (!is.character(actual) || !nzchar(actual)) next
  value <- tryCatch(encode_result(), error = identity)
  if (is.list(value)) rows[[i]]$identical_to_C <- identical(value, baseline)
}

result <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, output, row.names = FALSE, fileEncoding = "UTF-8")
print(result, row.names = FALSE)
