#!/usr/bin/env Rscript

# Reproduce the fastphylosig 0.3.0 plot_signal style gallery.
# All plotted panels are drawn by the package's native base-R plot_signal().

args <- commandArgs(trailingOnly = TRUE)
library_args <- sub("^--library=", "", grep("^--library=", args, value = TRUE))
out_arg <- grep("^--outdir=", args, value = TRUE)
if (length(out_arg) > 1L) stop("Use --outdir only once.", call. = FALSE)
out_dir_arg <- if (length(out_arg)) sub("^--outdir=", "", out_arg) else NULL
unknown_args <- args[!grepl("^--library=", args) & !grepl("^--outdir=", args)]
if (length(unknown_args)) {
  stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "),
       call. = FALSE)
}

if (length(library_args)) {
  missing_libraries <- library_args[!dir.exists(library_args)]
  if (length(missing_libraries)) {
    stop("A supplied R library directory does not exist.", call. = FALSE)
  }
  .libPaths(c(normalizePath(library_args, winslash = "/", mustWork = TRUE),
              .libPaths()))
}

if (!requireNamespace("fastphylosig", quietly = TRUE)) {
  stop("Install fastphylosig 0.3.0 before running this script.",
       call. = FALSE)
}
installed_version <- as.character(utils::packageVersion("fastphylosig"))
if (!identical(installed_version, "0.3.0")) {
  stop("This gallery requires fastphylosig 0.3.0.", call. = FALSE)
}
required_exports <- c("fast_signal", "fast_d", "fast_delta", "plot_signal")
missing_exports <- setdiff(
  required_exports, getNamespaceExports("fastphylosig")
)
if (length(missing_exports)) {
  stop("The installed package is missing a required public function.",
       call. = FALSE)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Install the ape package before running this script.",
       call. = FALSE)
}

out_dir <- if (is.null(out_dir_arg)) {
  file.path("docs", "release-0.3.0")
} else {
  out_dir_arg
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

n_tip <- 32L
seed_data <- 20261001L
seed_k <- 20261002L
seed_d <- 20261003L
seed_delta <- 20261004L
nsim_k <- 99L
nsim_d <- 99L
nsim_delta <- 3L
mcmc_sim_delta <- 60L
thin_delta <- 10L
burn_delta <- 20L

standardize <- function(x) {
  as.numeric(scale(as.numeric(x)))
}

set.seed(seed_data)
tree <- ape::rtree(n_tip)
bm_trait <- ape::rTraitCont(tree, model = "BM", sigma = 1)
tip_names <- tree$tip.label

set.seed(seed_data + 1L)
k_noise <- stats::rnorm(n_tip)
names(k_noise) <- tip_names
k_moderate <- 0.6 * standardize(bm_trait) +
  0.8 * stats::rnorm(n_tip)
names(k_moderate) <- tip_names
k_strong <- bm_trait
names(k_strong) <- tip_names

lambda_trait <- 0.7 * standardize(bm_trait) +
  0.5 * stats::rnorm(n_tip)
names(lambda_trait) <- tip_names

d_trait <- as.integer(bm_trait >= stats::median(bm_trait))
names(d_trait) <- tip_names

delta_trait <- as.character(cut(
  bm_trait,
  breaks = stats::quantile(bm_trait, probs = c(0, 1 / 3, 2 / 3, 1)),
  include.lowest = TRUE,
  labels = c("low", "middle", "high")
))
names(delta_trait) <- tip_names

k_multi <- cbind(
  Trait_A = k_noise,
  Trait_B = k_moderate,
  Trait_C = k_strong
)
rownames(k_multi) <- tip_names

demo_data <- data.frame(
  species = tip_names,
  K_single = k_strong[tip_names],
  K_trait_A = k_noise[tip_names],
  K_trait_B = k_moderate[tip_names],
  K_trait_C = k_strong[tip_names],
  lambda_trait = lambda_trait[tip_names],
  D_binary = d_trait[tip_names],
  Delta_categorical = delta_trait[tip_names],
  stringsAsFactors = FALSE
)
utils::write.csv(
  demo_data,
  file = file.path(out_dir, "plot_signal_demo_traits.csv"),
  row.names = FALSE,
  na = ""
)
ape::write.tree(
  tree,
  file = file.path(out_dir, "plot_signal_demo_tree.nwk")
)
saveRDS(
  list(
    package_version = installed_version,
    seeds = c(data = seed_data, K = seed_k, D = seed_d, Delta = seed_delta),
    n_tip = n_tip,
    settings = list(
      nsim_K = nsim_k,
      nsim_D = nsim_d,
      nsim_Delta = nsim_delta,
      mcmc_sim_Delta = mcmc_sim_delta,
      thin_Delta = thin_delta,
      burn_Delta = burn_delta
    ),
    tree = tree,
    traits = list(
      K_single = k_strong,
      K_multi = k_multi,
      lambda = lambda_trait,
      D = d_trait,
      Delta = delta_trait
    )
  ),
  file = file.path(out_dir, "plot_signal_demo_data.rds"),
  version = 3
)

open_gallery_device <- function(path, type, overview = FALSE) {
  width_px <- if (overview) 1600L else 1800L
  height_px <- if (overview) 1200L else 1350L
  width_in <- 12
  height_in <- 9
  if (identical(type, "png")) {
    grDevices::png(
      filename = path, width = width_px, height = height_px,
      res = width_px / width_in,
      type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows",
      pointsize = if (overview) 14 else 18
    )
  } else if (identical(type, "svg")) {
    grDevices::svg(
      filename = path, width = width_in, height = height_in,
      pointsize = 14,
      family = "sans"
    )
  } else {
    grDevices::pdf(
      file = path, width = width_in, height = height_in, pointsize = 14,
      family = "Helvetica", useDingbats = FALSE
    )
  }
}

render_gallery_plot <- function(name, draw, par_overrides = list()) {
  return_value <- NULL
  for (type in c("png", "svg", "pdf")) {
    path <- file.path(out_dir, paste0(name, ".", type))
    open_gallery_device(path, type)
    result <- tryCatch(
      {
        par_settings <- utils::modifyList(
          list(
            family = "sans",
            mar = c(5, 5.5, 4.5, 2),
            oma = c(0, 0, 0.5, 0),
            mgp = c(2.7, 0.75, 0),
            cex = 1.15,
            cex.axis = 1.05,
            cex.lab = 1.1,
            cex.main = 1.2,
            lwd = 1.2
          ),
          par_overrides
        )
        do.call(graphics::par, par_settings)
        draw()
      },
      finally = grDevices::dev.off()
    )
    if (identical(type, "png")) return_value <- result
  }
  invisible(return_value)
}

set.seed(seed_k)
k_single_fit <- fastphylosig::fast_signal(
  tree, data = k_strong, method = "K", test = TRUE, nsim = nsim_k,
  return_sim = TRUE, verbose = FALSE, progress = FALSE
)
k_multi_fit <- fastphylosig::fast_signal(
  tree, data = k_multi, method = "K", test = TRUE, nsim = nsim_k,
  return_sim = TRUE, verbose = FALSE, progress = FALSE
)

set.seed(seed_k)
draw_k_single <- function() fastphylosig::plot_signal(
    k_single_fit,
    main = "K permutation null: single trait",
    xlab = "Blomberg's K",
    ylab = "Density"
  )
draw_k_multi <- function() fastphylosig::plot_signal(
  k_multi_fit,
  main = "K permutation nulls: multiple traits",
  xlab = "Blomberg's K",
  ylab = "Density",
  label = "none"
)
plot_k_single <- render_gallery_plot(
  "plot_signal_K_single_trait", draw_k_single
)
plot_k_multi <- render_gallery_plot(
  "plot_signal_K_multi_trait", draw_k_multi
)

lambda_fit <- fastphylosig::fast_signal(
  tree, data = lambda_trait, method = "lambda", test = TRUE,
  lambda_profile_points = 101L, verbose = FALSE, progress = FALSE
)
draw_lambda <- function() fastphylosig::plot_signal(
    lambda_fit,
    main = "Pagel's lambda profile likelihood",
    xlab = "Pagel's lambda",
    ylab = "Profile log-likelihood"
  )
plot_lambda <- render_gallery_plot(
  "plot_signal_lambda_profile", draw_lambda
)

set.seed(seed_d)
d_fit <- fastphylosig::fast_d(
  tree, x = d_trait, test = TRUE, nsim = nsim_d,
  return_sim = TRUE, verbose = FALSE, progress = FALSE
)
draw_d_overlay <- function() fastphylosig::plot_signal(
    d_fit,
    main = "D calibration: Brownian and random nulls",
    xlab = "Fritz & Purvis D",
    ylab = "Density"
  )
plot_d_overlay <- render_gallery_plot(
  "plot_signal_D_random_brownian_overlay",
  draw_d_overlay,
  par_overrides = list(
    cex = 0.82, cex.axis = 0.95, cex.lab = 1, cex.main = 1.1
  )
)
plot_d_single_nulls <- render_gallery_plot(
  "plot_signal_D_single_nulls",
  function() {
    graphics::par(
      mfrow = c(1, 2),
      mar = c(4.5, 4.1, 4.1, 0.7),
      mgp = c(2.3, 0.65, 0),
      cex = 0.9,
      cex.axis = 0.84,
      cex.lab = 0.9,
      cex.main = 1.05
    )
    random_result <- fastphylosig::plot_signal(
      d_fit, null = "random",
      main = "Random association null",
      xlab = "Fritz & Purvis D",
      ylab = "Density"
    )
    brownian_result <- fastphylosig::plot_signal(
      d_fit, null = "brownian",
      main = "Brownian threshold null",
      xlab = "Fritz & Purvis D",
      ylab = "Density"
    )
    list(random = random_result, brownian = brownian_result)
  }
)

set.seed(seed_delta)
delta_warnings <- character()
delta_fit <- withCallingHandlers(
  fastphylosig::fast_delta(
    tree, x = delta_trait, test = TRUE, nsim = nsim_delta,
    mcmc_sim = mcmc_sim_delta, thin = thin_delta, burn = burn_delta,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ),
  warning = function(w) {
    delta_warnings <<- c(delta_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
plot_delta <- render_gallery_plot(
  "plot_signal_Delta_categorical_demo",
  function() fastphylosig::plot_signal(
    delta_fit,
    main = "Delta permutation demo (not for inference)",
    xlab = "Delta",
    ylab = "Density"
  )
)

base_plot_rect_pt <- function(fig, plt, width_pt, height_pt) {
  if (length(fig) != 4L || length(plt) != 4L ||
      any(!is.finite(c(fig, plt)))) {
    stop("Could not measure the base-R plot-area rectangle.", call. = FALSE)
  }
  figure_width <- fig[[2L]] - fig[[1L]]
  figure_height <- fig[[4L]] - fig[[3L]]
  c(
    left = (fig[[1L]] + plt[[1L]] * figure_width) * width_pt,
    bottom = (fig[[3L]] + plt[[3L]] * figure_height) * height_pt,
    right = (fig[[1L]] + plt[[2L]] * figure_width) * width_pt,
    top = (fig[[3L]] + plt[[4L]] * figure_height) * height_pt
  )
}

write_base_alignment_manifest <- function(rectangles, path,
                                          width_in, height_in) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required to write alignment QA.",
         call. = FALSE)
  }
  if (length(rectangles) != 4L || is.null(names(rectangles))) {
    stop("The overview alignment manifest requires four named panels.",
         call. = FALSE)
  }
  panels <- lapply(seq_along(rectangles), function(index) {
    row <- as.integer((index - 1L) %/% 2L)
    column <- as.integer((index - 1L) %% 2L)
    list(
      id = names(rectangles)[[index]],
      bbox_pt = unname(as.numeric(rectangles[[index]])),
      grid_id = "overview-2x2",
      row_start = row,
      row_stop = row + 1L,
      col_start = column,
      col_stop = column + 1L
    )
  })
  manifest <- list(
    schema_version = 1L,
    backend = "r-base",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels,
    row_groups = list(c("a", "b"), c("c", "d")),
    column_groups = list(c("a", "c"), c("b", "d"))
  )
  jsonlite::write_json(
    manifest, path = path, auto_unbox = TRUE, pretty = TRUE, digits = NA
  )
  invisible(manifest)
}

render_gallery_overview <- function() {
  overview_rectangles <- list()
  capture_overview_rectangle <- function(panel_id) {
    fig <- graphics::par("fig")
    plt <- graphics::par("plt")
    overview_rectangles[[panel_id]] <<- base_plot_rect_pt(
      fig, plt, width_pt = 12 * 72, height_pt = 9 * 72
    )
  }
  overview_par <- list(
    mfrow = c(2, 2),
    family = "sans",
    mar = c(4.2, 4.8, 3.8, 1.1),
    oma = c(0, 0, 0, 0),
    mgp = c(2.15, 0.55, 0),
    cex = 0.88,
    cex.axis = 0.86,
    cex.lab = 0.9,
    cex.main = 0.95,
    lwd = 1
  )
  for (type in c("png", "svg", "pdf")) {
    path <- file.path(
      out_dir, paste0("plot_signal_gallery_overview.", type)
    )
    open_gallery_device(path, type, overview = TRUE)
    tryCatch(
      {
        do.call(graphics::par, overview_par)
        draw_k_multi()
        if (identical(type, "pdf")) capture_overview_rectangle("a")
        draw_lambda()
        if (identical(type, "pdf")) capture_overview_rectangle("b")
        graphics::par(
          cex = 0.66, cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.25,
          xpd = FALSE
        )
        draw_d_overlay()
        if (identical(type, "pdf")) capture_overview_rectangle("c")
        graphics::par(
          cex = overview_par$cex,
          cex.axis = overview_par$cex.axis,
          cex.lab = overview_par$cex.lab,
          cex.main = overview_par$cex.main
        )
        fastphylosig::plot_signal(
          delta_fit,
          main = "Delta demo (not for inference)",
          xlab = "Delta",
          ylab = "Density"
        )
        if (identical(type, "pdf")) capture_overview_rectangle("d")
      },
      finally = grDevices::dev.off()
    )
    if (identical(type, "pdf")) {
      write_base_alignment_manifest(
        overview_rectangles,
        path = file.path(
          out_dir, "plot_signal_gallery_overview_alignment-layout.json"
        ),
        width_in = 12,
        height_in = 9
      )
    }
  }
  invisible(NULL)
}
render_gallery_overview()

plot_data_rows <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    return(do.call(rbind, lapply(x, plot_data_rows)))
  }
  data.frame(
    method = as.character(x$method[[1L]]),
    trait = if ("trait" %in% names(x)) as.character(x$trait) else NA_character_,
    estimate = if ("estimate" %in% names(x)) as.numeric(x$estimate) else NA_real_,
    p_value = if ("p_value" %in% names(x)) as.numeric(x$p_value) else NA_real_,
    nsim = if ("n_sim" %in% names(x)) as.integer(x$n_sim) else NA_integer_,
    plot_type = if ("plot_type" %in% names(x)) {
      as.character(x$plot_type[[1L]])
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}

summary_table <- rbind(
  plot_data_rows(plot_k_single),
  plot_data_rows(plot_k_multi),
  plot_data_rows(plot_lambda),
  plot_data_rows(plot_d_overlay),
  plot_data_rows(plot_d_single_nulls),
  plot_data_rows(plot_delta)
)
summary_table$n_species <- n_tip
summary_table$package_version <- installed_version
summary_table$seed_data <- seed_data
summary_table$seed_K <- seed_k
summary_table$seed_D <- seed_d
summary_table$seed_Delta <- seed_delta
summary_table$delta_demo_nsim <- nsim_delta
summary_table$delta_demo_mcmc_sim <- mcmc_sim_delta
summary_table$delta_demo_thin <- thin_delta
summary_table$delta_demo_burn <- burn_delta
summary_table$delta_warnings <- ifelse(
  summary_table$method == "Delta",
  paste(unique(delta_warnings), collapse = " | "),
  ""
)
utils::write.csv(
  summary_table,
  file = file.path(out_dir, "plot_signal_gallery_summary.csv"),
  row.names = FALSE,
  na = ""
)

cat("fastphylosig version:", installed_version, "\n")
cat("tips:", n_tip, "\n")
cat("Delta warnings captured:", length(delta_warnings), "\n")
cat("Generated PNG, SVG, PDF, and reproducible demo data.\n")
