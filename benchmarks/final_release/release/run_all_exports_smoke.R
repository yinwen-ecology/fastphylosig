# Smoke-test every public fastphylosig export from a fresh library.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_all_exports_smoke.R <library> <output>", call. = FALSE)
}

smoke_library <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)

main <- function() {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("cannot create output directory.", call. = FALSE)
  }

  log_con <- file(output, open = "wt", encoding = "UTF-8")
  sink(log_con, type = "output", split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
    try(close(log_con), silent = TRUE)
  }, add = TRUE)

  capture_call <- function(fun) {
    warnings <- character()
    messages <- character()
    err <- NULL
    value <- tryCatch(
      withCallingHandlers(
        fun(),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      ),
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
    list(value = value, warnings = warnings, messages = messages, error = err)
  }

  collapse_text <- function(value) {
    if (!length(value)) return("")
    paste(gsub("[\\r\\n]+", " ", as.character(value)), collapse = " | ")
  }

  expected_deprecation <- function(value) {
    grepl("match_phylo_data", value, fixed = TRUE) &&
      grepl("deprecated", value, fixed = TRUE)
  }

  public_exports <- c(
    "fast_signal", "fast_k", "fast_lambda", "fast_d", "fast_delta",
    "fast_ace", "prepare_tree", "cache_info", "match_tree_data",
    "match_phylo_data", "plot_signal", "check_tree", "resolve_tree"
  )

  package_path <- find.package(
    "fastphylosig", lib.loc = smoke_library, quiet = TRUE
  )
  if (length(package_path) != 1L) {
    stop("fastphylosig is not installed in the supplied library.",
         call. = FALSE)
  }
  package_path <- normalizePath(package_path, winslash = "/", mustWork = TRUE)
  if ("fastphylosig" %in% loadedNamespaces()) {
    stop("fastphylosig is already loaded; run this script in a fresh R session.",
         call. = FALSE)
  }

  .libPaths(unique(c(smoke_library, .libPaths())))
  loaded <- capture_call(function() {
    suppressPackageStartupMessages(
      library("fastphylosig", character.only = TRUE, lib.loc = smoke_library)
    )
  })
  if (!is.null(loaded$error)) {
    stop(paste0("failed to load fastphylosig: ", loaded$error), call. = FALSE)
  }
  if (length(loaded$warnings)) {
    stop(paste0("package load emitted warning(s): ",
                collapse_text(loaded$warnings)), call. = FALSE)
  }
  loaded_path <- normalizePath(
    getNamespaceInfo(asNamespace("fastphylosig"), "path"),
    winslash = "/", mustWork = TRUE
  )
  if (!identical(loaded_path, package_path)) {
    stop("fastphylosig was not loaded from the supplied library.", call. = FALSE)
  }
  if (!setequal(getNamespaceExports("fastphylosig"), public_exports)) {
    stop("installed NAMESPACE does not expose the expected 13 public exports.",
         call. = FALSE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required by the installed fastphylosig package.", call. = FALSE)
  }

  records <- data.frame(
    export = character(), status = character(), warnings = character(),
    error = character(), stringsAsFactors = FALSE
  )

  run_export <- function(name, fun, check, allow_deprecation = FALSE) {
    captured <- capture_call(fun)
    check_error <- NULL
    value_ok <- tryCatch(
      isTRUE(check(captured$value)),
      error = function(e) {
        check_error <<- conditionMessage(e)
        FALSE
      }
    )
    warning_ok <- if (isTRUE(allow_deprecation)) {
      length(captured$warnings) > 0L &&
        all(vapply(captured$warnings, expected_deprecation, logical(1L)))
    } else {
      length(captured$warnings) == 0L
    }
    ok <- is.null(captured$error) && is.null(check_error) && value_ok && warning_ok
    error_text <- captured$error
    if (!is.null(check_error)) {
      error_text <- paste0("result check failed: ", check_error)
    }
    records <<- rbind(
      records,
      data.frame(
        export = name,
        status = if (ok) "PASS" else "FAIL",
        warnings = collapse_text(captured$warnings),
        error = if (is.null(error_text)) "" else error_text,
        stringsAsFactors = FALSE
      )
    )
    cat(sprintf("[%s] %s", if (ok) "PASS" else "FAIL", name))
    if (length(captured$warnings)) {
      cat(" -- warning: ", collapse_text(captured$warnings), sep = "")
    }
    if (!is.null(error_text)) cat(" -- error: ", error_text, sep = "")
    cat("\n")
    invisible(captured$value)
  }

  make_tree20 <- function() {
    set.seed(20260913L)
    tree <- ape::rtree(20L)
    tree$tip.label <- paste0("sp", seq_len(20L))
    tree$edge.length <- seq(0.4, 2.2, length.out = nrow(tree$edge))
    ape::reorder.phylo(tree, "cladewise")
  }

  make_ace_tree <- function() {
    ape::read.tree(
      text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"
    )
  }

  fixed_permutations <- function(n, nsim) {
    out <- vapply(seq_len(nsim), function(i) {
      ((seq_len(n) + i - 2L) %% n) + 1L
    }, integer(n))
    t(out)
  }

  smoke_tree <- make_tree20()
  smoke_tips <- smoke_tree$tip.label
  smoke_cont <- stats::setNames(seq_along(smoke_tips) / 5, smoke_tips)
  smoke_binary <- stats::setNames(
    rep(c(0, 1), length.out = length(smoke_tips)), smoke_tips
  )
  smoke_categorical <- stats::setNames(
    factor(rep(c("A", "B"), length.out = length(smoke_tips)),
           levels = c("A", "B")),
    smoke_tips
  )

  run_export(
    "fast_signal",
    function() fast_signal(
      smoke_tree, data = smoke_cont, method = "K", test = FALSE,
      verbose = FALSE, progress = FALSE
    ),
    function(value) is.numeric(value) && length(value) == 1L &&
      all(is.finite(as.numeric(value)))
  )
  run_export(
    "fast_k",
    function() fast_k(
      smoke_tree, smoke_cont, test = FALSE, verbose = FALSE, progress = FALSE
    ),
    function(value) is.numeric(value) && length(value) == 1L &&
      all(is.finite(as.numeric(value)))
  )
  run_export(
    "fast_lambda",
    function() fast_lambda(
      smoke_tree, smoke_cont, lambda_profile = FALSE,
      verbose = FALSE, progress = FALSE
    ),
    function(value) is.list(value) && length(value$lambda) == 1L &&
      is.finite(value$lambda) && length(value$logL) == 1L &&
      is.finite(value$logL)
  )
  run_export(
    "fast_d",
    function() fast_d(
      smoke_tree, smoke_binary, test = FALSE, nsim = 4L,
      return_sim = FALSE, verbose = FALSE, progress = FALSE
    ),
    function(value) is.list(value) && length(value$DEstimate) == 1L &&
      is.finite(value$DEstimate) && identical(value$status, "ok")
  )
  run_export(
    "fast_delta",
    function() fast_delta(
      smoke_tree, smoke_categorical, test = FALSE,
      mcmc_sim = 1000L, thin = 10L, burn = 100L,
      model = "ER", verbose = FALSE, progress = FALSE
    ),
    function(value) is.list(value) && length(value$delta) == 1L &&
      is.finite(value$delta) && identical(value$status, "ok")
  )

  ace_tree <- make_ace_tree()
  ace_x <- stats::setNames(
    factor(c("A", "A", "B", "B", "A", "A", "B", "B"),
           levels = c("A", "B")),
    ace_tree$tip.label
  )
  run_export(
    "fast_ace",
    function() fast_ace(
      ace_x, phy = ace_tree, model = "ER", CI = FALSE,
      marginal = FALSE, progress = FALSE
    ),
    function(value) is.list(value) && length(value$loglik) == 1L &&
      is.finite(value$loglik) && is.numeric(value$rates) &&
      length(value$rates) > 0L && all(is.finite(value$rates))
  )

  smoke_ctx <- NULL
  run_export(
    "prepare_tree",
    function() {
      smoke_ctx <<- prepare_tree(smoke_tree)
      smoke_ctx
    },
    function(value) inherits(value, "fastphylosig_tree")
  )
  run_export(
    "cache_info",
    function() cache_info(smoke_ctx),
    function(value) is.list(value) && length(value$bytes_used) == 1L &&
      is.numeric(value$bytes_used) && is.data.frame(value$entries)
  )
  run_export(
    "match_tree_data",
    function() match_tree_data(smoke_tree, data = smoke_cont, verbose = FALSE),
    function(value) is.list(value) &&
      identical(value$matched_species, smoke_tips)
  )
  run_export(
    "match_phylo_data",
    function() match_phylo_data(smoke_tree, smoke_cont, verbose = FALSE),
    function(value) is.list(value) &&
      identical(value$matched_species, smoke_tips),
    allow_deprecation = TRUE
  )
  run_export(
    "plot_signal",
    function() {
      plot_path <- tempfile("fastphylosig-export-smoke-", fileext = ".png")
      device_open <- FALSE
      on.exit({
        if (device_open && !is.null(grDevices::dev.list())) {
          try(grDevices::dev.off(), silent = TRUE)
        }
        if (file.exists(plot_path)) unlink(plot_path, force = TRUE)
      }, add = TRUE)
      grDevices::png(plot_path, width = 800L, height = 600L)
      device_open <- TRUE
      permutations <- fixed_permutations(length(smoke_tips), nsim = 4L)
      fit <- fast_k(
        smoke_tree, smoke_cont, test = TRUE, nsim = 4L,
        permutations = permutations, return_sim = TRUE,
        verbose = FALSE, progress = FALSE
      )
      value <- plot_signal(fit)
      grDevices::dev.off()
      device_open <- FALSE
      value
    },
    function(value) is.data.frame(value) && nrow(value) >= 1L
  )
  run_export(
    "check_tree",
    function() check_tree(smoke_tree),
    function(value) is.list(value) && is.logical(value$ready) &&
      all(c("K", "lambda", "D", "Delta") %in% names(value$ready_by_signal))
  )
  run_export(
    "resolve_tree",
    function() resolve_tree(smoke_tree, signal = "K"),
    function(value) inherits(value, "phylo") &&
      isTRUE(attr(value, "fastphylosig_resolution")$ready)
  )

  missing_exports <- setdiff(public_exports, records$export)
  duplicate_exports <- records$export[duplicated(records$export)]
  if (length(missing_exports) || length(duplicate_exports)) {
    stop("public export coverage is incomplete or duplicated.", call. = FALSE)
  }
  if (any(records$status != "PASS")) {
    failed <- records$export[records$status != "PASS"]
    stop(paste0("export smoke failed: ", paste(failed, collapse = ", ")),
         call. = FALSE)
  }

  cat("EXPORT_SMOKE=PASS\n")
  cat("PACKAGE_VERSION=", as.character(packageVersion("fastphylosig")), "\n",
      sep = "")
  cat("LIBRARY=", loaded_path, "\n", sep = "")
  cat("EXPORTS=", paste(public_exports, collapse = ","), "\n", sep = "")
  invisible(records)
}

main()
