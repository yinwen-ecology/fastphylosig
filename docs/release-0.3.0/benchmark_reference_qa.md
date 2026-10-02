# Benchmark figure QA

## Evidence and interpretation

- Claim represented: matched single-call point-estimator runtimes for `fastphylosig` 0.3.0 versus `phytools` 2.5.2, shown separately for Blomberg's K and Pagel's lambda.
- Scope: `test = FALSE`; K is its point estimator, not the permutation test. The timing starts after tree/trait creation, matching, and preparation. One public call per method was timed with `Sys.time`; the 10 paired replicates per cell are summarized by the ratio of medians.
- Statistic: median `phytools` runtime / median `fastphylosig` runtime. Values above 1 favor `fastphylosig`; values below 1 remain visible and favor `phytools`.
- Coverage: 4 species counts (50, 100, 500, 2,000) × 3 signal scenarios × 2 metrics = 24 cells; 10 matched pairs per cell (240 pairs; 480 timed calls).
- This is a descriptive runtime figure only. It does not establish numerical-result equivalence or provide an inferential test/confidence interval.

## Benchmark result summary

- K reference/fast ratio range: 0.1434–151.3968. The lowest cell is n=50 (no signal): the reference median is faster; the largest is n=2,000 (moderate signal), 151.4×.
- Lambda reference/fast ratio range: 6.4802–23,591.4444. The largest is n=2,000 (moderate signal), 23,591.4×.
- No cells were filtered for falling below 1. In particular, the K results at n=50 and n=100 show slower `fastphylosig` calls; the plot preserves this regression visibly.

## Data and plot QA

- Matched-input audit: PASS (240 pairs; all 24 cells have 10 pairs; fast/reference species counts equal the target count; input hashes match; zero removed trait values; all timed-call statuses are `ok`).
- Data grid and ratio guard: PASS (all 24 cells present; finite positive ratios; values below 1 retained).
- Static source validator: 17 PASS, 4 WARN, 0 FAIL. Warnings: generic static syntax check asks for backend parsing (R `parse()` passed); no TIFF export; raster is 300 dpi rather than the 600-dpi default; static width detector does not infer the explicit 183-mm setting. This is a GitHub preview package, not a journal-raster submission package.
- R rendering and output: PASS; PNG, SVG, and PDF rendered at 183 × 115 mm. Final PNG was reviewed at full resolution; title/subtitle, both panels, legend, and caption are not clipped, and no visible text/mark collision was found. Some scenario points are numerically close and overlap in the panels; the legend explicitly identifies their color, marker, and line-type encodings.
- Facet alignment: PASS (R/ggplot2 gtable measurement; 1 comparison, 0 failures/warnings; tolerance 1.5 pt). Actual panel boxes are recorded in `benchmark_reference_speedup.alignment-layout.json`; full audit is `benchmark_reference_speedup.alignment.json`.
- PDF text size: the skill's `audit_pdf_text.py --min-pt 5` reports 9 text operators at 1 pt. Cross-checking the final PDF with PyMuPDF 1.28.2 finds 25 rendered text spans, minimum 6.5 pt, and 0 spans below 5 pt. The discrepancy is consistent with the primary auditor reading the Cairo PDF text matrix/operator size without its scale; therefore the script result is recorded as a limitation, not relabeled as a PASS. Span-level evidence is in `benchmark_reference_speedup.pdf-text-audit.json`.
- Collision audit: the automated collision tool returns `FIX BEFORE DELIVERY` (15 FAIL, 1 WARN). Inspection of its JSON findings against individual PyMuPDF text-span boxes and the full-size PNG indicates coarse union text-trace boxes: e.g., separate x ticks, facet titles, subtitle/caption lines and legend labels are grouped into overlapping bounding boxes; the reported legend stroke does not cross the actual “No signal” span. The audit output is preserved as `benchmark_reference_speedup.collision-audit.json` with its overlay PDF. Because the automated audit does not pass, this figure is suitable as a local/GitHub preview for user review, not certified here as publication-ready; resolve the auditor limitations or obtain an independent final collision review before journal submission.
- Editable text: SVG/PDF retain vector text; verify selection after opening the exports in the target editor.

## Provenance

- Source data: `benchmark_reference_summary.csv`; pair-level records: `benchmark_reference_pairs.csv`.
- Benchmark metadata: `benchmark_reference_metadata.csv`; R session details: `benchmark_reference_session_info.txt`.
- Metadata plot-script SHA-256 captured at benchmark time: `7524130660ae1d3c141d8c0f153b79ea0175584553b7fab07d0d7446d482f297`.
- Final plotting script SHA-256 used for this export: `5d6cfeeb615ef41d7f19f70de933f66f6eab4a2ff6e0b07e97f7b1f4baba964c`.
- Summary CSV SHA-256: `9e49227b2544d7eefbb2922afb8b4ac95cbd14ad3af60178a44bff4680206387`.
- Pair CSV SHA-256: `74df84f9b35081e3fdc7c79232c2ad82c744f65da6f3d8d410f438ba9f131797`.
