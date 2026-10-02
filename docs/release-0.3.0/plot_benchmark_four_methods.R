# Render the four-method matched runtime figure for the fastphylosig 0.3.0 page.
# Run from the repository root. This script reads audited summaries; it does not
# run benchmarks or rewrite their source data.

required_packages <- c("digest", "ggplot2", "grid", "jsonlite", "patchwork",
                       "ragg", "scales", "svglite")
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
k_summary_path <- file.path(figure_dir, "benchmark_reference_summary.csv")
d_summary_path <- file.path(figure_dir, "benchmark_reference_d_delta_summary.csv")
k_metadata_path <- file.path(figure_dir, "benchmark_reference_metadata.csv")
d_metadata_path <- file.path(figure_dir, "benchmark_reference_d_delta_metadata.csv")
required_files <- c(k_summary_path, d_summary_path, k_metadata_path, d_metadata_path)
if (any(!file.exists(required_files))) stop("Benchmark source file missing.", call. = FALSE)

k_summary <- utils::read.csv(k_summary_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
d_summary <- utils::read.csv(d_summary_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
k_metadata <- utils::read.csv(k_metadata_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
d_metadata <- utils::read.csv(d_metadata_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
metadata_value <- function(metadata, field) {
  value <- metadata$value[metadata$field == field]
  if (length(value) != 1L) NA_character_ else as.character(value[[1L]])
}
expected_tarball_sha <-
  "909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff"
if (!identical(metadata_value(k_metadata, "package_version"), "0.3.0") ||
    !identical(metadata_value(d_metadata, "package_version"), "0.3.0") ||
    !identical(metadata_value(k_metadata, "source_tarball_sha256"), expected_tarball_sha) ||
    !identical(metadata_value(d_metadata, "source_tarball_sha256"), expected_tarball_sha)) {
  stop("Benchmark provenance does not match the approved v0.3.0 source tarball.",
       call. = FALSE)
}

required_summary <- c(
  "method", "scenario", "n_tip", "n_pairs", "fast_median_seconds",
  "reference_median_seconds", "reference_over_fast_ratio_of_medians",
  "test", "preparation_included"
)
if (!all(required_summary %in% names(k_summary)) ||
    !all(required_summary %in% names(d_summary))) {
  stop("A benchmark summary lacks required columns.", call. = FALSE)
}
source <- rbind(
  k_summary[, required_summary, drop = FALSE],
  d_summary[, required_summary, drop = FALSE]
)
expected_grid <- expand.grid(
  method = c("K", "lambda", "D", "Delta"),
  scenario = c("none", "moderate", "strong"),
  n_tip = c(50L, 100L, 500L, 2000L), stringsAsFactors = FALSE
)
cell_key <- function(x) paste(x$method, x$scenario, x$n_tip, sep = "|")
if (nrow(source) != 48L || !setequal(cell_key(source), cell_key(expected_grid)) ||
    anyDuplicated(cell_key(source)) || any(source$n_pairs != 10L) ||
    any(!is.finite(source$reference_over_fast_ratio_of_medians)) ||
    any(source$reference_over_fast_ratio_of_medians <= 0) ||
    any(!is.finite(source$fast_median_seconds)) ||
    any(!is.finite(source$reference_median_seconds)) ||
    any(source$fast_median_seconds <= 0) || any(source$reference_median_seconds <= 0) ||
    any(source$preparation_included != FALSE)) {
  stop("The summaries do not match the complete 48-cell prepared-input grid.",
       call. = FALSE)
}
if (any(k_summary$test != FALSE) ||
    any(d_summary$test[d_summary$method == "Delta"] != FALSE) ||
    any(d_summary$test[d_summary$method == "D"] != TRUE)) {
  stop("Estimator/permutation settings do not match the documented methods.",
       call. = FALSE)
}

source$method <- factor(source$method, levels = c("K", "lambda", "D", "Delta"))
source$scenario <- factor(
  source$scenario, levels = c("none", "moderate", "strong"),
  labels = c("No signal", "Moderate signal", "Strong signal")
)
scenario_colours <- c("No signal" = "#0072B2", "Moderate signal" = "#D89000",
                      "Strong signal" = "#009E73")
scenario_shapes <- c("No signal" = 16, "Moderate signal" = 17,
                     "Strong signal" = 15)
scenario_linetypes <- c("No signal" = "solid", "Moderate signal" = "dashed",
                        "Strong signal" = "dotdash")
metric_titles <- c(
  K = "Blomberg's K (test = FALSE)",
  lambda = "Pagel's lambda (test = FALSE)",
  D = "D (199 permutations)",
  Delta = "Delta* (two-chain MCMC)"
)
metric_plot <- function(metric, show_x_title = FALSE) {
  data <- source[source$method == metric, , drop = FALSE]
  data$n_tip_factor <- factor(data$n_tip, levels = c(50, 100, 500, 2000))
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = n_tip_factor, y = reference_over_fast_ratio_of_medians,
                 colour = scenario, shape = scenario, linetype = scenario,
                 group = scenario)
  ) +
    ggplot2::geom_hline(yintercept = 1, linewidth = 0.45,
                        linetype = "dashed", colour = "#555555") +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 2, stroke = 0.25) +
    ggplot2::scale_colour_manual(values = scenario_colours, name = "Signal scenario") +
    ggplot2::scale_shape_manual(values = scenario_shapes, name = "Signal scenario") +
    ggplot2::scale_linetype_manual(values = scenario_linetypes,
                                   name = "Signal scenario") +
    ggplot2::scale_x_discrete(
      limits = c("50", "100", "500", "2000"),
      labels = c("50", "100", "500", "2,000"),
      name = if (show_x_title) "Number of species" else NULL
    ) +
    ggplot2::scale_y_log10(
      name = "Reference / fast runtime",
      breaks = scales::breaks_log(n = 5),
      labels = scales::label_number(accuracy = 0.1, trim = TRUE)
    ) +
    ggplot2::labs(title = unname(metric_titles[[metric]])) +
    ggplot2::theme_classic(base_size = 8, base_family = "sans") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = 7.2),
      axis.text = ggplot2::element_text(size = 6.8, colour = "#222222"),
      plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6.8),
      legend.key.width = grid::unit(10, "mm"),
      plot.margin = ggplot2::margin(2, 4, 2, 4, unit = "mm")
    )
}
panel_plots <- list(metric_plot("K"), metric_plot("lambda"),
                    metric_plot("D", TRUE), metric_plot("Delta", TRUE))
caption_sentences <- c(
  "Speedup = median(reference) / median(fastphylosig); >1 favours fastphylosig, <1 favours reference.",
  "Points summarize 10 paired one-call timings per cell; tree, trait, and analysis-context preparation was excluded.",
  "K/lambda use test=FALSE. D compares fast_d(test=TRUE, nsim=199) with caper::phylo.d(permut=199).",
  "Delta uses two MCMC chains (10,000 iterations/chain, thin=10, burn=100) and is compared with the pinned delta() reference.",
  "Delta caveat: reference warnings in 112/120 calls (NaNs, Inf replacement, complex coercion, non-finite gradients); fast Delta had one optimizer-limit warning. The reference returns no convergence object.",
  "Runtime only; no equal-accuracy or convergence claim. X positions are categorical/equally spaced; each panel has its own logarithmic y-axis."
)
caption_text <- paste(unlist(lapply(caption_sentences, strwrap, width = 112)),
                      collapse = "\n")
figure <- patchwork::wrap_plots(panel_plots, ncol = 2, guides = "collect") +
  patchwork::plot_annotation(
    title = "Matched single-call runtime: fastphylosig 0.3.0 vs reference implementations",
    subtitle = "Reference / fastphylosig ratio of median wall-clock time; horizontal rule marks parity.",
    caption = caption_text, tag_levels = "a",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.2, lineheight = 1),
      plot.caption = ggplot2::element_text(size = 6.1, hjust = 0, lineheight = 1.03),
      plot.tag = ggplot2::element_text(size = 8, face = "bold")
    )
  ) & ggplot2::theme(legend.position = "bottom")

output_stem <- "benchmark_reference_four_methods_speedup"
output_paths <- file.path(figure_dir, paste0(output_stem, c(".png", ".svg", ".pdf")))
source_path <- file.path(figure_dir, "benchmark_reference_four_methods_source.csv")
qa_path <- file.path(figure_dir, "benchmark_reference_four_methods_qa.md")
alignment_path <- file.path("docs", "release-0.3.0", paste0(output_stem, ".alignment.json"))
alignment_report_path <- file.path("docs", "release-0.3.0",
                                   paste0(output_stem, ".alignment-audit.json"))
utils::write.csv(source, source_path, row.names = FALSE, na = "")
width_mm <- 183
height_mm <- 142
alignment_status <- "NOT AUDITABLE: alignment helper/auditor was not supplied."
audit_script <- Sys.getenv("NATURE_FIGURE_ALIGNMENT_AUDITOR", unset = "")
python <- Sys.getenv("NATURE_FIGURE_PYTHON_EXE", unset = "")
if (nzchar(audit_script) && file.exists(audit_script) &&
    nzchar(python) && file.exists(python)) {
  # Measure final patchwork panel rectangles on a PDF device matching export
  # dimensions. Resolve flexible gtable units by their null-unit weights.
  grob <- patchwork::patchworkGrob(figure)
    panel_rows <- grob$layout[
      grepl("^panel(?:-[0-9]+)?$", grob$layout$name, perl = TRUE), , drop = FALSE
    ]
    if (nrow(panel_rows) != 4L) stop("Expected four patchwork panel cells.", call. = FALSE)
    panel_rows <- panel_rows[order(panel_rows$t, panel_rows$l), , drop = FALSE]
    width_in <- width_mm / 25.4
    height_in <- height_mm / 25.4
    probe <- tempfile(fileext = ".pdf")
    grDevices::pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
    grid::grid.newpage()
    grid::grid.draw(grob)
    grid::grid.force()
    resolve_units <- function(units, extent_pt, axis) {
      types <- grid::unitType(units)
      flexible <- types == "null"
      resolved <- vapply(seq_along(units), function(i) {
        if (flexible[[i]]) return(0)
        if (axis == "width") grid::convertWidth(units[i], "pt", valueOnly = TRUE)
        else grid::convertHeight(units[i], "pt", valueOnly = TRUE)
      }, numeric(1))
      if (any(!is.finite(resolved))) return(rep(NA_real_, length(units)))
      if (any(flexible)) {
        weights <- as.numeric(units[flexible])
        remaining <- extent_pt - sum(resolved)
        if (!is.finite(remaining) || remaining <= 0 || any(weights <= 0)) {
          return(rep(NA_real_, length(units)))
        }
        resolved[flexible] <- remaining * weights / sum(weights)
      }
      resolved
    }
    widths_pt <- resolve_units(grob$widths, width_in * 72, "width")
    heights_pt <- resolve_units(grob$heights, height_in * 72, "height")
    grDevices::dev.off()
    unlink(probe)
    if (any(!is.finite(widths_pt)) || any(!is.finite(heights_pt)) ||
        abs(sum(widths_pt) - width_in * 72) > 0.05 ||
        abs(sum(heights_pt) - height_in * 72) > 0.05) {
      stop("Could not measure final patchwork panel rectangles.", call. = FALSE)
    }
    x_edges <- c(0, cumsum(widths_pt))
    top_edges <- c(0, cumsum(heights_pt))
    total_height <- sum(heights_pt)
    panels <- lapply(seq_len(nrow(panel_rows)), function(i) {
      row <- panel_rows[i, ]
      list(
        id = letters[i],
        bbox_pt = unname(c(x_edges[row$l], total_height - top_edges[row$b + 1L],
                           x_edges[row$r + 1L], total_height - top_edges[row$t])),
        grid_id = "patchwork-grid",
        row_start = as.integer(row$t - 1L), row_stop = as.integer(row$b),
        col_start = as.integer(row$l - 1L), col_stop = as.integer(row$r)
      )
    })
    jsonlite::write_json(list(
      schema_version = 1L, backend = "r-patchwork-gtable-resolved-null-units",
      figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
      panels = panels,
      row_groups = list(list(id = "top", panels = c("a", "b")),
                        list(id = "bottom", panels = c("c", "d"))),
      column_groups = list(list(id = "left", panels = c("a", "c")),
                           list(id = "right", panels = c("b", "d")))
    ), alignment_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    audit_args <- c(shQuote(audit_script), shQuote(alignment_path),
                    "--json-out", shQuote(alignment_report_path), "--strict")
    audit_output <- suppressWarnings(system2(python, audit_args, stdout = TRUE,
                                             stderr = TRUE))
    audit_status <- attr(audit_output, "status")
    if (is.null(audit_status)) audit_status <- 0L
    if (length(audit_output)) message(paste(audit_output, collapse = "\n"))
    if (audit_status != 0L) {
      stop("Panel alignment JSON audit failed; inspect the saved audit report.",
           call. = FALSE)
    }
  alignment_status <- paste(
    "PASS: four panel rectangles measured in R at export dimensions and",
    "audited by the nature-figure alignment auditor."
  )
} else {
  alignment_status <- paste(
    "NOT AUDITABLE: optional alignment QA requires",
    "NATURE_FIGURE_ALIGNMENT_AUDITOR and NATURE_FIGURE_PYTHON_EXE."
  )
}

ggplot2::ggsave(output_paths[[1L]], figure, device = ragg::agg_png,
                width = width_mm, height = height_mm, units = "mm", dpi = 300,
                bg = "white")
svg_relative_path <- file.path("docs", "release-0.3.0", basename(output_paths[[2L]]))
pdf_relative_path <- file.path("docs", "release-0.3.0", basename(output_paths[[3L]]))
svglite::svglite(svg_relative_path, width = width_mm / 25.4,
                 height = height_mm / 25.4, bg = "white")
print(figure)
grDevices::dev.off()
grDevices::cairo_pdf(pdf_relative_path, width = width_mm / 25.4,
                     height = height_mm / 25.4, family = "sans")
print(figure)
grDevices::dev.off()

ratio_ranges <- lapply(split(source, source$method), function(x) {
  range(x$reference_over_fast_ratio_of_medians)
})
max_cells <- lapply(split(source, source$method), function(x) {
  x[which.max(x$reference_over_fast_ratio_of_medians), , drop = FALSE]
})
format_cell <- function(x) paste0(x$method[[1L]], "/", x$scenario[[1L]], "/n",
                                  x$n_tip[[1L]], " = ",
                                  signif(x$reference_over_fast_ratio_of_medians[[1L]], 5), "x")
final_plot_sha <- digest::digest(
  file.path(figure_dir, "plot_benchmark_four_methods.R"), algo = "sha256",
  serialize = FALSE, file = TRUE
)
qa <- c(
  "# Four-method benchmark figure QA", "", "## Figure and evidence", "",
  "- Outputs: `benchmark_reference_four_methods_speedup.png`, `.svg`, and `.pdf` (R/ggplot2; 183 x 142 mm).",
  "- These are direct fastphylosig 0.3.0 versus corresponding-reference timings; no v0.2-to-v0.3 ratio was multiplied or inferred.",
  "- All metrics include n = 50, 100, 500, and 2,000 and all three signal scenarios, with 10 paired calls per cell.",
  "- One public call per method was timed with tree, trait, and analysis context prepared beforehand; preparation was excluded.",
  "- D is a 199-randomization test; Delta is a two-chain MCMC call with 10,000 iterations/chain, thin 10, burn 100; K/lambda use test=FALSE.",
  "- Plotted value is median(reference time) / median(fast time), not median paired ratio. Ratios below 1 remain included. Species counts are equally spaced categorical x positions; panels have independent logarithmic y axes.",
  "- Runtime comparison only: the figure does not establish equal accuracy or estimator validity.",
  "- Delta: reference warnings in 112/120 calls (NaNs, Inf replacement, imaginary-part coercion, non-finite-gradient warnings); fast Delta had one optimizer iteration-limit warning. The reference returns no convergence object; no convergence certification is made.",
  "", "## Observed ratio ranges", "",
  paste0("- K: ", signif(ratio_ranges$K[[1L]], 4), " to ", signif(ratio_ranges$K[[2L]], 4), "x; fastest cell: ", format_cell(max_cells$K), "."),
  paste0("- Lambda: ", signif(ratio_ranges$lambda[[1L]], 4), " to ", signif(ratio_ranges$lambda[[2L]], 4), "x; fastest cell: ", format_cell(max_cells$lambda), "."),
  paste0("- D: ", signif(ratio_ranges$D[[1L]], 4), " to ", signif(ratio_ranges$D[[2L]], 4), "x; fastest cell: ", format_cell(max_cells$D), "."),
  paste0("- Delta: ", signif(ratio_ranges$Delta[[1L]], 4), " to ", signif(ratio_ranges$Delta[[2L]], 4), "x; fastest cell: ", format_cell(max_cells$Delta), "."),
  "", "## Automated QA", "",
  "- Grid completeness, 10 pairs/cell, positive finite ratios, and preparation-excluded contract: PASS (48 cells).",
  paste0("- Panel alignment: ", alignment_status),
  "- PDF text audit, collision audit, and visual inspection are recorded separately; successful rendering alone does not certify them.",
  "", "## Provenance", "",
  paste0("- Validated source tarball SHA-256: ", expected_tarball_sha, "."),
  paste0("- K/lambda summary SHA-256: ", digest::digest(k_summary_path, algo = "sha256", serialize = FALSE, file = TRUE), "."),
  paste0("- D/Delta summary SHA-256: ", digest::digest(d_summary_path, algo = "sha256", serialize = FALSE, file = TRUE), "."),
  paste0("- D/Delta raw calls SHA-256: ", metadata_value(d_metadata, "raw_runs_sha256"), "."),
  paste0("- Delta reference source SHA-256: ", metadata_value(d_metadata, "Delta_reference_source_sha256"), "."),
  paste0("- D/Delta benchmark process exit code: ", metadata_value(d_metadata, "original_benchmark_process_exit_code"), "; raw calls completed, but post-run pair-order aggregation exited 1 and was reconstructed from immutable raw rows (details in `benchmark_reference_d_delta_qa.md`)."),
  paste0("- Plot script SHA-256: ", final_plot_sha, "."),
  "- R/package versions are listed in `benchmark_reference_d_delta_smoke_20261001T181645Z_session_info.txt`; it is the immediately preceding matched session, not a session dump from the full run.",
  "- The pinned Delta reference source is not redistributed; reproducing that reference requires a user-provided source file and SHA verification as documented in the D/Delta QA record."
)
writeLines(qa, file.path(figure_dir, "benchmark_reference_four_methods_qa.md"),
           useBytes = TRUE)
message("Wrote four-method figure: ", paste(basename(output_paths), collapse = ", "))
