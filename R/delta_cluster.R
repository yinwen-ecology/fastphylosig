# Load exactly the installation used by the calling process. lib.loc does not
# change .libPaths(), so forwarding only library paths can silently select an
# older package. A base-environment closure avoids namespace deserialization
# loading that older package before we have pinned the intended library.
.delta_initialize_cluster <- function(cl) {
  ns <- asNamespace("fastphylosig")
  package_path <- normalizePath(getNamespaceInfo(ns, "path"),
                                winslash = "/", mustWork = TRUE)
  package_version <- as.character(getNamespaceVersion(ns))
  initialize <- function(package_path, package_version, library_paths) {
    expected <- normalizePath(package_path, winslash = "/", mustWork = TRUE)
    .libPaths(unique(c(dirname(expected), library_paths)))
    ns <- loadNamespace("fastphylosig", lib.loc = dirname(expected))
    actual <- normalizePath(getNamespaceInfo(ns, "path"),
                            winslash = "/", mustWork = TRUE)
    version <- as.character(getNamespaceVersion(ns))
    if (!identical(actual, expected) || !identical(version, package_version)) {
      stop("Delta worker loaded a different fastphylosig installation: ",
           actual, " (", version, "); expected ", expected, " (",
           package_version, "). Restart the worker with the same installation.",
           call. = FALSE)
    }
    list(path = actual, version = version, library_paths = .libPaths())
  }
  # Source references can retain the package namespace even when the closure
  # environment is baseenv(); strip them before serialization to the worker.
  initialize <- utils::removeSource(initialize)
  environment(initialize) <- baseenv()
  invisible(parallel::clusterCall(cl, initialize, package_path,
                                 package_version, .libPaths()))
}