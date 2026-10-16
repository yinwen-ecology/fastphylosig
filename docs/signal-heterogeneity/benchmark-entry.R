# Historical baseline libraries predate the final API name. This lookup is
# confined to benchmark controls; the current package exposes fast_hetero only.
hetero_benchmark <- getExportedValue("fastphylosig",
  if ("fast_hetero" %in% getNamespaceExports("fastphylosig")) "fast_hetero" else
    "fast_signal_heterogeneity")
