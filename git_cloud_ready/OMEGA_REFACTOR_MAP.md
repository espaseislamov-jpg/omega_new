# Omega Refactor Map

This file is the working map for replacing the current monolithic processing code with a cleaner, library-based engine.

## Current Pipeline

1. Load input.
   `omega_core.io.load_batches()` reads single chromatograms and multi-sample `CHROMTAB.CSV` files, then returns ordered batches.

2. Normalize each chromatogram.
   `omega_core.io.finalize_chromatogram_dataframe()` keeps the useful time region and produces `x_corrected` / `y`.

3. Baseline correction.
   `omega_core.signal.add_baseline()` uses a custom Chebyshev/quantile baseline.
   `omega_core.pipeline` can select `omega_core.signal.add_asls_baseline()` for internally detected wide C22 clusters.
   If cluster quality is bad, `omega_core.signal.add_arpls_baseline()` tries `pybaselines.arpls`.

4. Smoothing and derivatives.
   `omega_core.signal.add_smoothing_and_derivatives()` selects a Savitzky-Golay window and produces `y_smooth`, `dy`, `d2y`.

5. Peak detection.
   `omega_core.signal.detect_peak_candidates()` combines `scipy.find_peaks`, `peak_widths`, targeted C18/C20/C22 augmentation, and optional `pyOpenMS` peak assistance.

6. Target matching.
   `omega_core.matching.match_targets_to_peaks()` maps detected peaks to fatty-acid targets by RT, order, and cluster overrides.

7. Cluster repair.
   `omega_core.clusters` now owns local non-fit repair layers:
   `refine_c18_c20_cluster_matches`, `refine_overlapped_c22_cluster_areas`, `refine_cluster_areas_by_local_valleys`,
   `tighten_overwide_c22_cluster_tails`, `refine_small_peak_integrations`.
   The remaining delegated layers are fit/deconvolution based:
   `recover_missing_c22_components_with_fit`, `recover_underintegrated_c20_components_with_fit`,
   `recover_overlapped_c18_components_with_fit`, `refine_overwide_c22_cluster_with_pvfit`.
   These are isolated behind `omega_core.legacy_fit`.

8. Omega calculation.
   `omega_core.metrics.compute_omega()` calculates strict and corrected omega with several C20/C22/C18 correction models.

9. Confidence.
   `omega_core.metrics.assess_confidence()` scores how likely a sample needs manual review.

10. GUI.
   `ChromatogramApp` owns file selection, batch navigation, plotting, tables, confidence popup, and preview windows.

## What Is Currently Fragile

- Peak detection, integration boundaries, and cluster repair are too tightly coupled.
- Some rescue layers are chemically useful, but hidden inside generic processing flow.
- Baseline choice is mixed with cluster quality decisions.
- The confidence score is useful, but it is not yet a hard gate for expensive alternative algorithms.
- The GUI is now mostly a shell over `omega_core`, but the old monolithic file still contains compatibility wrappers and fit recovery code.

## Refactor Target

The production shape should be:

1. `omega_io`
   Batch loading, CHROMTAB parsing, reference target loading.

2. `omega_signal`
   Baseline, smoothing, derivative calculation, peak candidates.

3. `omega_clusters`
   C18/C20/C22/C24 cluster models, boundaries, fits, quality metrics.

4. `omega_matching`
   Fatty-acid target assignment, order constraints, RT priors.

5. `omega_metrics`
   Omega calculation, debug tables, confidence scoring.

6. `omega_gui`
   Tkinter only. No processing algorithms hidden here.

7. `omega_regression`
   Mandatory validation against all reference batches before any engine change is accepted.

## Replacement Map

### Baseline

Current:
Custom Chebyshev baseline plus optional `pybaselines.arpls`.

Replace with:
`pybaselines.Baseline` as the main baseline provider. Candidate algorithms:
`asls`, `iasls`, `airpls`, `arpls`, `drpls`, `aspls`, `psalsa`, `derpsalsa`.

Why:
`pybaselines` exposes many published 1D baseline algorithms behind one API.

Implementation:
Create a baseline selector that evaluates candidates by internal signal criteria, not reference answers:
negative corrected area, baseline slope, cluster completeness, and valley clarity.

Sources:
- https://pybaselines.readthedocs.io/en/stable/api/Baseline.html
- https://pypi.org/project/pybaselines/

### Peak Candidate Picking

Current:
`scipy.find_peaks`, custom derivative geometry, targeted cluster recovery, plus optional `pyOpenMS`.

Replace with:
Keep `scipy.find_peaks` as the fast first pass, but make `pyOpenMS.PeakPickerChromatogram` the second official candidate source.

Why:
`PeakPickerChromatogram` explicitly searches left/right borders and handles overlapping chromatographic peaks.

Sources:
- https://pyopenms.readthedocs.io/en/latest/apidocs/_autosummary/pyopenms/pyopenms.PeakPickerChromatogram.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html

### Peak Area Integration

Current:
`np.trapezoid` over boundaries chosen by custom heuristics.

Replace with:
`pyOpenMS.PeakIntegrator` for final area, background, and shape metrics when boundaries are known.

Why:
It supports trapezoid, Simpson, intensity-sum integration and baseline approaches such as `base_to_base`.
It also computes shape metrics useful for confidence.

Sources:
- https://pyopenms.readthedocs.io/en/latest/apidocs/_autosummary/pyopenms/pyopenms.PeakIntegrator.html
- https://openms.de/documentation/PeakIntegrator_8h.html

### Cluster Deconvolution

Current:
Many hand-written rescue layers for C18/C20/C22.

Experimental implementation:
`Omega_cluster_engine.py` now uses raw-direct local cluster processing. It does
not subtract a global baseline or reuse a smoothed signal for integration.
Pseudo-Voigt fitting estimates a two-parameter local linear background together
with the peak components. The valley fallback uses a single detection-only
Savitzky-Golay profile and integrates the unsmoothed local signal.

Replace with:
A constrained local cluster model:
fixed peak order, bounded center drift, shared/bounded widths, local baseline, and robust weighted residuals.

Libraries:
`lmfit` for `PseudoVoigtModel`, `SkewedVoigtModel`, and `ExponentialGaussianModel`.

Why:
`lmfit` provides composable fitting models and parameter bounds, which matches known peak order and expected RT windows.

Sources:
- https://lmfit.github.io/lmfit-py/builtin_models.html
- https://lmfit.github.io/lmfit-py/model.html

### ChromatoPy

Current status:
Installed with `--no-deps` because the normal installation fails on `hdbscan` under Python 3.13 without Microsoft C++ Build Tools.
The top-level `chromatopy` import also pulls `PyQt5`, so direct package import is not lightweight.

What is usable:
`chromatopy/FID/FID_Integration_functions.py` can be loaded directly.
It provides valley finding, derivative boundary finding, Gaussian/multi-Gaussian fitting, skewed Gaussian, and area ensembles.

Where it helps:
As an experimental candidate generator for cluster boundaries and fit diagnostics, not as a drop-in replacement for our whole application.

### Anchor-Gated Gaussian Candidate

Current status:
Implemented experimentally in `Omega_cluster_engine.py --mode anchored`.

Approach:
Use four independent stable anchor peaks to estimate sample-level RT warp first:
C16:1N7, C16:0, C24:1N9, C24:0. The model classifies the shift pattern as
quiet, coherent, gradient, noisy, or missing. Only coherent/gradient samples can
try Gaussian-only cluster deconvolution on raw local signal. The candidate is
not allowed to replace production unless anchor shift is meaningful, the omega
movement is analytically significant but bounded, and internal shape/area/order
gates pass.

Why it matters:
It removes part of the over-heavy Pseudo-Voigt stack for cases where the main
problem is RT drift, not peak profile precision.

Latest regression:
ALL MAE improved from 0.159349 to 0.157527 without reducing the +/-0.5 count.

Source:
- https://github.com/GerardOtiniano/chromatoPy
- https://pypi.org/project/chromatopy/

### hplc-py

Candidate:
`hplc-py` has a `Chromatogram` class with peak detection, baseline correction, deconvolution, and fitted peak tables.

Why it matters:
It is better documented than ChromatoPy and is built around a clean chromatogram-processing object.

Risks:
It is HPLC-oriented, may be slow on our full corpus, and has GPL licensing.

Source:
- https://github.com/cremerlab/hplc-py

### Vendor Loading

Current:
We parse exported CSV/CHROMTAB.

Optional future:
`rainbow-api` reads raw chromatography and MS vendor binary formats.

Use only if we later want direct vendor-file ingestion.

Source:
- https://rainbow-api.readthedocs.io/

## Next Implementation Steps

1. Replace `omega_core.legacy_fit` with native fit/deconvolution recovery code in `omega_core.cluster_fit`.
2. Keep `omega_regression.py` mandatory before accepting engine changes.
3. Replace baseline selection with a `pybaselines` selector.
4. Replace final area calculation with `pyOpenMS.PeakIntegrator`.
5. Replace C18/C20/C22 rescue layers with constrained local cluster models.
6. Keep further GUI changes behind the stable `omega_core` API.
