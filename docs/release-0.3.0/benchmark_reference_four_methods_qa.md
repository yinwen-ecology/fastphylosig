# Four-method benchmark figure QA

## Figure and evidence

- Outputs: `benchmark_reference_four_methods_speedup.png`, `.svg`, and `.pdf` (R/ggplot2; 183 x 142 mm).
- These are direct fastphylosig 0.3.0 versus corresponding-reference timings; no v0.2-to-v0.3 ratio was multiplied or inferred.
- All metrics include n = 50, 100, 500, and 2,000 and all three signal scenarios, with 10 paired calls per cell.
- One public call per method was timed with tree, trait, and analysis context prepared beforehand; preparation was excluded.
- D is a 199-randomization test; Delta is a two-chain MCMC call with 10,000 iterations/chain, thin 10, burn 100; K/lambda use test=FALSE.
- Plotted value is median(reference time) / median(fast time), not median paired ratio. Ratios below 1 remain included. Species counts are equally spaced categorical x positions; panels have independent logarithmic y axes.
- Runtime comparison only: the figure does not establish equal accuracy or estimator validity.
- Delta: reference warnings in 112/120 calls (NaNs, Inf replacement, imaginary-part coercion, non-finite-gradient warnings); fast Delta had one optimizer iteration-limit warning. The reference returns no convergence object; no convergence certification is made.

## Observed ratio ranges

- K: 0.1434 to 151.4x; fastest cell: K/Moderate signal/n2000 = 151.4x.
- Lambda: 6.48 to 23590x; fastest cell: lambda/Moderate signal/n2000 = 23591x.
- D: 5.27 to 224.1x; fastest cell: D/Strong signal/n2000 = 224.09x.
- Delta: 23.89 to 37.69x; fastest cell: Delta/Strong signal/n500 = 37.693x.

## Automated QA

- Grid completeness, 10 pairs/cell, positive finite ratios, and preparation-excluded contract: PASS (48 cells).
- Panel alignment: PASS: four panel rectangles measured in R at export dimensions and audited by the nature-figure alignment auditor.
- R-source static validation: 16 PASS, 5 WARN, 0 FAIL. The warnings concern website-preview-only PNG resolution/no TIFF, no static target width, and a conservative log-guard heuristic; an R render completed successfully.
- The PDF `Tf`-operator font audit reported 1 pt for every text run. This is a known transformation-matrix misread for this R PDF; a PyMuPDF text-span cross-check found 61 spans, 6.1 pt minimum and 10 pt maximum, with no span below 5 pt.
- The rendered-collision auditor reported FAIL (35 findings) and WARN (1). Its text boxes aggregate transformed R text runs (for example, every caption line is assigned one large shared box), and it also flags intended legend key strokes. This machine verdict is retained as a non-pass; final-size PNG inspection found no actual text collision or clipping.
- Visual QA: final-size R-rendered PNG inspected; panel labels/titles, legend, and wrapped caption are readable and not clipped. Figure rendering is not a claim of statistical accuracy.

## Provenance

- Validated source tarball SHA-256: 909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff.
- K/lambda summary SHA-256: 9e49227b2544d7eefbb2922afb8b4ac95cbd14ad3af60178a44bff4680206387.
- D/Delta summary SHA-256: 4f88293c0be9bc2cd914f2c6be8d216939bc135d122271fc34dff983b2a8115c.
- D/Delta raw calls SHA-256: e46a67a58f3e8519b1f50ae8ed31e88005053e947949517ec647e6fe99400318.
- Delta reference source SHA-256: 0dc1e987c13f0e972b0e7c3b4829b4d0fdc41ff0de3856bae4803f9ce0115175.
- D/Delta benchmark process exit code: 1; raw calls completed, but post-run pair-order aggregation exited 1 and was reconstructed from immutable raw rows (details in `benchmark_reference_d_delta_qa.md`).
- Plotting script SHA-256 used to render the current exports: 6e3ca80a9d538bc067bda29ea7e65f2ccf30ebc62af4207b4c687b6bd6ec95b5. The public script was then changed only to remove host-specific automatic QA-tool discovery in favor of optional `NATURE_FIGURE_ALIGNMENT_AUDITOR` and `NATURE_FIGURE_PYTHON_EXE` variables; plot geometry/exports were not changed. Current source SHA-256: 6ff391076fd773ccd9315b01bf25809ea3e855b875ac2cb18980b52445053325.
- Matched benchmark environment: R 4.6.1 (ucrt), platform x86_64-w64-mingw32/x64, Windows 11 x64 build 26300; fastphylosig 0.3.0, phytools 2.5.2, caper 1.0.4, ape 5.8-1, expm 1.0-0, Rcpp 1.1.2, mvtnorm 1.4-2, Matrix 1.7-5, and nlme 3.1-169. These versions come from the immediately preceding matched smoke session; the full-run session dump was unavailable after the post-run aggregation exit 1.
- The pinned Delta reference source is not redistributed; reproducing that reference requires a user-provided source file and SHA verification as documented in the D/Delta QA record.
