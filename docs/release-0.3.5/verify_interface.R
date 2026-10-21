# Version/documentation-only preparation: preserve every export and its formals.
args <- commandArgs(TRUE)
.libPaths(c(args[[1L]], .libPaths()))
library(fastphylosig)
ns <- asNamespace("fastphylosig")
exports <- sort(getNamespaceExports(ns))
interface <- lapply(setNames(exports, exports), function(name) {
  value <- get(name, ns, inherits = FALSE)
  if (is.function(value)) formals(value) else typeof(value)
})
saveRDS(interface, args[[2L]])
if (length(args) > 2L) stopifnot(identical(interface, readRDS(args[[3L]])))
cat(length(exports), "exports and their argument defaults checked\n")
