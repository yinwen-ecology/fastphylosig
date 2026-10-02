# Render the v0.3.0 matched K/lambda benchmark for the GitHub package page.
# Run in R from the repository root after benchmark_reference.R succeeds.

required_packages <- c(
  "digest", "ggplot2", "grid", "jsonlite", "ragg", "scales", "svglite"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required plotting packages: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}
figure_dir <- file.path(repo_root, "docs", "release-0.3.0")
summary_path <- file.path(figure_dir, "benchmark_reference_summary.csv")
pairs_path <- file.path(figure_dir, "benchmark_reference_pairs.csv")
metadata_path <- file.path(figure_dir, "benchmark_reference_metadata.csv")
for (path in c(summary_path, pairs_path, metadata_path)) {
  if (!file.exists(path)) stop("Missing benchmark source file: ", basename(path),
                               call. = FALSE)
}

summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
pairs <- utils::read.csv(pairs_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
metadata_value <- function(field) {
  value <- metadata$value[metadata$field == field]
  if (length(value) != 1L) NA_character_ else as.character(value[[1L]])
}
expected_tarball_sha256 <-
  "909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff"
if (!identical(metadata_value("package_version"), "0.3.0") ||
    !identical(metadata_value("source_tarball_sha256"), expected_tarball_sha256)) {
  stop("Benchmark provenance does not match the approved v0.3.0 tarball.",
       call. = FALSE)
}
required_summary <- c(
  "package_version", "method", "scenario", "n_tip", "n_pairs",
  "fast_median_seconds", "reference_median_seconds",
  "reference_over_fast_ratio_of_medians"
)
if (!all(required_summary %in% names(summary))) {
  stop("Benchmark summary is missing required columns.", call. = FALSE)
}
if (!all(c("same_tree_hash", "same_trait_hash", "species_match",
           "complete_pair", "test", "prepare_excluded") %in% names(pairs))) {
  stop("Paired source data lack required audit fields.", call. = FALSE)
}
if (any(!pairs$same_tree_hash) || any(!pairs$same_trait_hash) ||
    any(!pairs$species_match) || any(!pairs$complete_pair) ||
    any(pairs$test != FALSE) || any(!pairs$prepare_excluded)) {
  stop("Pair audit failed; figure generation stopped.", call. = FALSE)
}
if (any(pairs$n_fast_species != pairs$n_reference_species) ||
    any(pairs$n_fast_species != pairs$matched_n_species)) {
  stop("Fast/reference species counts do not match.", call. = FALSE)
}
expected_grid <- expand.grid(
  method = c("K", "lambda"),
  scenario = c("none", "moderate", "strong"),
  n_tip = c(50L, 100L, 500L, 2000L),
  stringsAsFactors = FALSE
)
key <- function(x) paste(x$method, x$scenario, x$n_tip, sep = "|")
if (nrow(summary) != nrow(expected_grid) ||
    !setequal(key(summary), key(expected_grid)) ||
    any(summary$n_pairs != 10L) ||
    !all(summary$package_version == "0.3.0") ||
    any(!is.finite(summary$reference_over_fast_ratio_of_medians)) ||
    any(summary$reference_over_fast_ratio_of_medians <= 0)) {
  stop("Summary does not match the expected 24-cell v0.3.0 grid.",
       call. = FALSE)
}
if (nrow(pairs) != 240L || nrow(pairs) != 10L * nrow(expected_grid)) {
  stop("Expected exactly 240 complete matched pairs.", call. = FALSE)
}
speed_ratio <- summary$reference_over_fast_ratio_of_medians
if (any(!is.finite(speed_ratio)) || !all(speed_ratio > 0)) {
  stop("Log-scaled speed ratios must all be finite and strictly positive.",
       call. = FALSE)
}

summary$method <- factor(summary$method, levels = c("K", "lambda"))
summary$scenario <- factor(
  summary$scenario, levels = c("none", "moderate", "strong"),
  labels = c("No signal", "Moderate signal", "Strong signal")
)
summary$n_tip <- factor(summary$n_tip, levels = c(50, 100, 500, 2000))

scenario_colours <- c(
  "No signal" = "#0072B2",
  "Moderate signal" = "#D89000",
  "Strong signal" = "#009E73"
)
scenario_shapes <- c("No signal" = 16, "Moderate signal" = 17,
                     "Strong signal" = 15)
scenario_linetypes <- c("No signal" = "solid", "Moderate signal" = "dashed",
                        "Strong signal" = "dotdash")

figure <- ggplot2::ggplot(
  summary,
  ggplot2::aes(
    x = n_tip, y = reference_over_fast_ratio_of_medians,
    colour = scenario, shape = scenario, linetype = scenario, group = scenario
  )
) +
  ggplot2::geom_hline(
    yintercept = 1, linewidth = 0.45, linetype = "dashed", colour = "#555555"
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.2, stroke = 0.3) +
  ggplot2::facet_wrap(
    ~method, nrow = 1, labeller = ggplot2::as_labeller(c(K = "Blomberg's K",
                                                        lambda = "Pagel's lambda"))
  ) +
  ggplot2::scale_colour_manual(values = scenario_colours, name = "Signal scenario") +
  ggplot2::scale_shape_manual(values = scenario_shapes, name = "Signal scenario") +
  ggplot2::scale_linetype_manual(values = scenario_linetypes,
                                 name = "Signal scenario") +
  ggplot2::scale_x_discrete(drop = FALSE, name = "Number of species") +
  ggplot2::scale_y_log10(
    name = "Reference / fastphylosig median runtime",
    breaks = scales::breaks_log(n = 5),
    labels = scales::label_number(accuracy = 0.1, trim = TRUE)
  ) +
  ggplot2::labs(
    title = "Matched single-call runtime: fastphylosig 0.3.0 vs phytools",
    subtitle = paste(
      "Speedup = median reference / median fastphylosig call time.",
      "Above 1 favours fastphylosig; below 1 favours phytools.",
      sep = "\n"
    ),
    caption = paste(
      "Each point summarizes 10 matched pairs; line and marker encode signal scenario.",
      "One public call timed per method; tree/data matching and preparation excluded.",
      "test = FALSE; K is the point estimator, not the permutation test.",
      sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    axis.title = ggplot2::element_text(size = 8),
    axis.text = ggplot2::element_text(size = 7.2, colour = "#222222"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8.2, face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 7.5),
    legend.text = ggplot2::element_text(size = 7.2),
    legend.key.width = grid::unit(11, "mm"),
    plot.title = ggplot2::element_text(size = 10, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.5, lineheight = 1.05),
    plot.caption = ggplot2::element_text(size = 6.5, hjust = 0,
                                        lineheight = 1.05),
    panel.spacing = grid::unit(8, "mm"),
    plot.margin = ggplot2::margin(5, 7, 4, 5, unit = "mm")
  )

output_files <- file.path(
  "docs", "release-0.3.0",
  c("benchmark_reference_speedup.png", "benchmark_reference_speedup.svg",
    "benchmark_reference_speedup.pdf", "benchmark_reference_qa.md",
    "benchmark_reference_speedup.alignment-layout.json")
)
if (any(file.exists(file.path(figure_dir, basename(output_files))))) {
  stop("Refusing to overwrite existing figure/QA outputs.", call. = FALSE)
}
width_mm <- 183
height_mm <- 115

# Measure the final ggplot facet panel rectangles on an R PDF device at the
# exact export dimensions. The JSON is backend-neutral and audited separately.
write_facet_alignment_manifest <- function(plot, manifest_path,
                                           width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grob <- ggplot2::ggplotGrob(plot)
  panel_rows <- grob$layout[
    grepl("^panel-[0-9]+-[0-9]+$", grob$layout$name), , drop = FALSE
  ]
  if (nrow(panel_rows) != 2L) {
    stop("Expected exactly two rendered metric facet panels.", call. = FALSE)
  }
  panel_rows <- panel_rows[order(panel_rows$t, panel_rows$l), , drop = FALSE]
  probe_path <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(probe_path, width = width_in, height = height_in,
                       family = "Arial")
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(probe_path)
  }, add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grid::grid.force()
  resolve_gtable_units <- function(units, extent_pt, axis) {
    unit_types <- grid::unitType(units)
    flexible <- unit_types == "null"
    resolved <- vapply(seq_along(units), function(index) {
      if (flexible[[index]]) return(0)
      if (axis == "width") {
        grid::convertWidth(units[index], "pt", valueOnly = TRUE)
      } else {
        grid::convertHeight(units[index], "pt", valueOnly = TRUE)
      }
    }, numeric(1))
    if (any(!is.finite(resolved))) return(rep(NA_real_, length(units)))
    if (any(flexible)) {
      weights <- as.numeric(units[flexible])
      remaining <- extent_pt - sum(resolved)
      if (!is.finite(remaining) || remaining <= 0 ||
          any(!is.finite(weights)) || any(weights <= 0)) {
        return(rep(NA_real_, length(units)))
      }
      resolved[flexible] <- remaining * weights / sum(weights)
    }
    resolved
  }
  widths_pt <- resolve_gtable_units(grob$widths, width_in * 72, "width")
  heights_pt <- resolve_gtable_units(grob$heights, height_in * 72, "height")
  if (any(!is.finite(widths_pt)) || any(!is.finite(heights_pt)) ||
      abs(sum(widths_pt) - width_in * 72) > 0.05 ||
      abs(sum(heights_pt) - height_in * 72) > 0.05) {
    stop("Could not measure final facet panel dimensions in R.", call. = FALSE)
  }
  x_edges <- c(0, cumsum(widths_pt))
  top_edges <- c(0, cumsum(heights_pt))
  total_height_pt <- sum(heights_pt)
  panels <- lapply(seq_len(nrow(panel_rows)), function(index) {
    row <- panel_rows[index, ]
    list(
      id = letters[index],
      bbox_pt = unname(c(
        x_edges[row$l],
        total_height_pt - top_edges[row$b + 1L],
        x_edges[row$r + 1L],
        total_height_pt - top_edges[row$t]
      )),
      grid_id = "ggplot-facet-grid",
      row_start = as.integer(row$t - 1L),
      row_stop = as.integer(row$b),
      col_start = as.integer(row$l - 1L),
      col_stop = as.integer(row$r)
    )
  })
  manifest <- list(
    schema_version = 1L,
    backend = "r-ggplot2-gtable",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels,
    row_groups = list(list(id = "metric-panels", panels = c("a", "b")))
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE,
                       pretty = TRUE, digits = NA)
  grDevices::dev.off()
  device_open <- FALSE
  unlink(probe_path)
  invisible(manifest)
}

write_facet_alignment_manifest(
  figure, output_files[[5L]], width_mm = width_mm, height_mm = height_mm
)
# ragg's direct device rejects the relative path in this Windows workspace;
# ggsave resolves the absolute target correctly for the PNG export.
ggplot2::ggsave(
  file.path(figure_dir, basename(output_files[[1L]])), figure,
  device = ragg::agg_png, width = width_mm, height = height_mm,
  units = "mm", dpi = 300, bg = "white"
)
svglite::svglite(
  output_files[[2L]], width = width_mm / 25.4, height = height_mm / 25.4,
  bg = "white"
)
print(figure)
grDevices::dev.off()
grDevices::cairo_pdf(
  output_files[[3L]], width = width_mm / 25.4,
  height = height_mm / 25.4, family = "Arial"
)
print(figure)
grDevices::dev.off()

ratio_range <- range(summary$reference_over_fast_ratio_of_medians)
method_ranges <- lapply(split(summary, summary$method), function(z) {
  range(z$reference_over_fast_ratio_of_medians)
})
fastest_cells <- lapply(split(summary, summary$method), function(z) {
  z[which.max(z$reference_over_fast_ratio_of_medians), , drop = FALSE]
})
format_cell <- function(z) {
  paste0(
    z$method[[1L]], "/", z$scenario[[1L]], "/n", z$n_tip[[1L]],
    " (", signif(z$reference_over_fast_ratio_of_medians[[1L]], 5), "x)"
  )
}
metadata_plot_sha <- metadata_value("plot_script_sha256")
final_plot_sha <- digest::digest(
  file.path("docs", "release-0.3.0", "plot_benchmark_reference.R"),
  algo = "sha256",
  serialize = FALSE, file = TRUE
)
qa_text <- c(
  "# Benchmark figure QA",
  "",
  "## Evidence and interpretation",
  "",
  "- Claim: matched point-estimate call times for fastphylosig 0.3.0 and phytools vary by metric, signal scenario, and species count.",
  "- Panels: Blomberg's K and Pagel's lambda; each contains all three signal scenarios and species counts 50, 100, 500, and 2,000.",
  "- Statistic plotted: ratio of the median phytools call time to the median fastphylosig call time. Ratios above 1 favour fastphylosig; values below 1 remain visible and favour phytools.",
  "- Each cell contains 10 paired replicates. The plot is descriptive; no inferential test or confidence interval is shown.",
  "- Scope: test = FALSE; K point estimation only, not the permutation K test. Tree, trait, matching and preparation are outside the timer.",
  "",
  "## Automated checks",
  "",
  paste0("- Grid coverage: PASS (", nrow(summary), " cells; 10 pairs per cell)."),
  paste0("- Matched inputs and species counts: PASS (", nrow(pairs), " pairs; all pair audit fields passed)."),
  paste0("- Ratio range retained: PASS (", signif(ratio_range[[1L]], 4), " to ",
         signif(ratio_range[[2L]], 4), "; no values below 1 were filtered)."),
  paste0("- K reference/fast ratio range: ", signif(method_ranges$K[[1L]], 4),
         " to ", signif(method_ranges$K[[2L]], 4), "."),
  paste0("- Lambda reference/fast ratio range: ",
         signif(method_ranges$lambda[[1L]], 4), " to ",
         signif(method_ranges$lambda[[2L]], 4), "."),
  paste0("- Fastest cell by K: ", format_cell(fastest_cells$K), "; by lambda: ",
         format_cell(fastest_cells$lambda), "."),
  "- Backend and exports: R/ggplot2; PNG, SVG, and PDF generated at 183 x 115 mm.",
  "- Scenario separation: redundant colour, marker-shape and line-type encodings.",
  "- Alignment: final facet plot-area rectangles measured in R; layout manifest emitted for backend-neutral JSON audit.",
  "- Editable text: SVG/PDF use vector text; verify selection after opening exports.",
  "- Visual review: PENDING; inspect each panel and full figure at final size.",
  "",
  "## Provenance",
  "",
  "- Source data: `benchmark_reference_summary.csv`; pair-level records: `benchmark_reference_pairs.csv`.",
  "- Benchmark metadata: `benchmark_reference_metadata.csv`.",
  paste0("- Metadata plot-script SHA-256 captured at benchmark time: ",
         metadata_plot_sha, "."),
  paste0("- Final plotting script SHA-256 used for this export: ",
         final_plot_sha, "."),
  paste0("- Summary CSV SHA-256: ", digest::digest(summary_path, algo = "sha256", serialize = FALSE, file = TRUE), "."),
  paste0("- Pair CSV SHA-256: ", digest::digest(pairs_path, algo = "sha256", serialize = FALSE, file = TRUE), "."),
  "- R session details: `benchmark_reference_session_info.txt`."
)
writeLines(qa_text, output_files[[4L]], useBytes = TRUE)

message("Wrote speedup figure: ", basename(output_files[[1L]]),
        ", ", basename(output_files[[2L]]), ", and ", basename(output_files[[3L]]))
