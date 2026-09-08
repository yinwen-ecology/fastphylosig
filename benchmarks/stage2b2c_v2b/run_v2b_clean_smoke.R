# Canonicalization Contract V2-B clean-install smoke.

smoke_library <- Sys.getenv("FASTPHYLOSIG_SMOKE_LIBRARY", unset = "")
if (!nzchar(smoke_library) || !dir.exists(smoke_library)) {
  stop("FASTPHYLOSIG_SMOKE_LIBRARY must name the clean installation library.",
       call. = FALSE)
}
.libPaths(c(smoke_library, .libPaths()))

smoke_log <- Sys.getenv("FASTPHYLOSIG_SMOKE_LOG", unset = "")
if (nzchar(smoke_log)) {
  connection <- file(smoke_log, open = "wt", encoding = "UTF-8")
  sink(connection)
  sink(connection, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(connection)
  }, add = TRUE)
}

library(fastphylosig)

set.seed(20260908)
tree <- ape::rtree(40)
tips <- tree$tip.label
trait <- setNames(seq_len(40) / 5 + stats::rnorm(40, sd = 0.1), tips)
binary <- setNames(rep(c(0L, 1L), 20), tips)
categorical <- ape::rTraitDisc(
  tree, model = "ER", k = 3, rate = 1.5, states = letters[1:3]
)
categorical <- setNames(as.character(categorical[tips]), tips)

checked <- check_tree(tree)
context <- prepare_tree(tree)
results <- list(
  fast_k = fast_k(context, trait, test = FALSE, progress = FALSE,
                  verbose = FALSE),
  fast_lambda = fast_lambda(context, trait, test = FALSE, progress = FALSE,
                            verbose = FALSE),
  fast_d = fast_d(context, binary, test = FALSE, progress = FALSE,
                  verbose = FALSE),
  fast_delta = fast_delta(context, categorical, test = FALSE, model = "ER",
                          progress = FALSE, verbose = FALSE),
  fast_ace = fast_ace(categorical, prepared = context, CI = FALSE,
                      model = "ER", progress = FALSE),
  fast_signal = fast_signal(context, data = trait, method = "K",
                            test = FALSE, progress = FALSE, verbose = FALSE)
)

stopifnot(
  isTRUE(checked$ready),
  inherits(context, "fastphylosig_tree"),
  all(vapply(results, function(result) !is.null(result), logical(1)))
)

cat("CLEAN_INSTALL_SMOKE=PASS\n")
cat("PACKAGE_VERSION=", as.character(packageVersion("fastphylosig")), "\n",
    sep = "")
cat("ENTRY_POINTS=check_tree,prepare_tree,fast_k,fast_lambda,fast_d,",
    "fast_delta,fast_ace,fast_signal\n", sep = "")
