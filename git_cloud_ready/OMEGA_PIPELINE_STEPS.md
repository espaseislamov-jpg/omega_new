# Omega Pipeline Steps

This is the current production pipeline after the first modular refactor.

## 1. IO

Module: `omega_core.io`

Responsibilities:
- load reference targets from `reference_targets_reverted_c22fixed.json`;
- load single chromatogram CSV files;
- split CHROMTAB-style batch files into ordered sample batches.

Current implementation:
Implemented directly in `omega_core.io`.

## 2. Signal Preparation

Module: `omega_core.signal`

Responsibilities:
- baseline correction;
- shape-gated alternative baseline fallback;
- optional arPLS fallback baseline;
- Savitzky-Golay smoothing and derivatives;
- first-pass peak candidate detection.

Current implementation:
Implemented directly in `omega_core.signal` for Chebyshev baseline, shape-gated `pybaselines.asls`, optional `pybaselines.arpls`, Savitzky-Golay smoothing, primary `scipy.find_peaks` picking, and targeted C18/C20/C22 peak augmentation.
Optional pyOpenMS assistance is implemented directly in this module.

Next replacement target:
Replace baseline selection with a `pybaselines` selector and keep the smoothing rules behind this module.

## 3. Matching

Module: `omega_core.matching`

Responsibilities:
- assign detected peak candidates to fatty-acid targets;
- preserve known elution order;
- keep RT priors and target identity separate from integration details.

Current implementation:
Implemented directly in `omega_core.matching`, including RT-shift estimation, order-based assignment, and C18/C20/C22 target overrides.

## 4. Cluster Refinement

Module: `omega_core.clusters`

Responsibilities:
- refine C18/C20/C22 clusters;
- recover missing components;
- split overlapped components;
- tighten small or overwide peaks.

Current implementation:
Implemented directly in `omega_core.clusters` for:
- C18/C20 local cluster rematching;
- C22 overlapped split by local minima;
- C18/C20/C22 valley reintegration;
- C22 tail tightening;
- small-peak sharp reintegration.

Still delegated through `omega_core.legacy_fit` compatibility wrappers:
- C22 missing-component fit recovery;
- C20 underintegrated EPA fit recovery;
- C18 overlapped fit recovery;
- overwide C22 PV fit refinement.

Next replacement target:
Replace `omega_core.legacy_fit` with a real `omega_core.cluster_fit` implementation.

## 5. Metrics

Module: `omega_core.metrics`

Responsibilities:
- calculate omega metrics;
- calculate cluster quality;
- build confidence score;
- attach baseline mode and confidence metadata.

Current implementation:
Omega calculation and cluster quality are implemented directly in `omega_core.metrics`.
Confidence scoring is implemented directly in `omega_core.metrics`.

## 6. Pipeline

Module: `omega_core.pipeline`

Responsibilities:
- orchestrate signal preparation, peak detection, matching, cluster refinement, metrics, and fallback baseline;
- select ASLS shape fallback only for internally detected wide C22 clusters;
- expose `process_batch()` and `process_file()` as stable engine entry points.

This module is now the engine seam. GUI code calls `omega_core.process_batch()` and `omega_core.load_batches()` directly.

## 7. Regression

Module: `omega_regression.py`

Responsibilities:
- run all available `test_bigbatch_*.xlsx` references against matching CSV batches;
- report full statistics and outliers;
- save `omega_regression_current.xlsx`.

Rule:
Every engine change must run this before being accepted.

## Experimental Raw-Direct Cluster Engine

Module: `Omega_cluster_engine.py`

The experimental cluster path no longer applies a global baseline or a globally
smoothed signal. `prepare_raw_signal()` preserves the input intensity unchanged.
Each cluster is fitted directly against raw intensity with a local linear
background included in the Pseudo-Voigt model. Valley mode uses one
Savitzky-Golay profile only to locate apices and split valleys; areas are
integrated from the unsmoothed local signal.

## Experimental Anchor-Gated Gaussian Engine

Module: `Omega_cluster_engine.py --mode anchored`

This path keeps the production result as the safe base, then estimates local RT
warp from four independent stable anchor peaks: C16:1N7, C16:0, C24:1N9,
C24:0. C18:2N6C is intentionally not used as an anchor, because it belongs to
one of the difficult target clusters.

The anchor model classifies each sample as quiet, coherent, gradient, noisy, or
missing. Only coherent/gradient samples with meaningful local RT movement can
try a raw-signal Gaussian-only deconvolution candidate for C18/C20/C22.
Candidate replacement is accepted only after internal checks for order, width,
cluster area ratio, and bounded omega movement. Tiny omega changes are ignored,
because they add refit risk without meaningful analytical gain.

Current test result:
- production ALL MAE: 0.159349;
- anchored Gaussian ALL MAE: 0.157527;
- threshold counts stayed unchanged at +/-0.5: 132/136.

Interpretation:
This is a promising simplification path, but still experimental. It should stay
outside the GUI until the acceptance rules are stress-tested further.

## Experimental No-Smoothing Checks

Module: `Omega_cluster_engine.py`

Modes added for diagnostics:
- `--mode core-nosmooth`: keeps the normal Chebyshev baseline but sets
  `y_smooth = y_corrected`, so peak detection and boundary geometry operate on
  the unsmoothed corrected signal.
- `--mode core-raw-boundaries`: keeps Savitzky-Golay for peak search and
  derivatives, but sets peak/cluster boundary metrics to raw-corrected signal.
- `--mode anchored-nosmooth`: keeps production as the base and disables
  Savitzky-Golay only inside anchored Gaussian candidate detection.

Current result:
- full `core-nosmooth` is much worse: MAE 0.298297;
- raw-only boundaries are also worse: MAE 0.201479;
- lighter boundary-smoothing weights from 0.25 to 0.55 were worse than current
  production weights;
- current production weights remain best in this sweep: MAE 0.159349.

Interpretation:
Smoothing should not be removed globally. It is stabilizing boundary search more
than it is distorting area. The safer next direction is not "no smoothing", but
localized quality gates around cases where smoothed boundaries visibly drift.

## Final Boundary Judge V0

Module: `omega_core/clusters.py`

After all cluster-specific refinements, production now builds a conservative
`baseexpand` candidate. The candidate expands an integration boundary only when
the raw corrected signal reaches a low local support threshold within a short
search distance. It is bounded by neighboring apices, maximum width, and maximum
allowed area growth.

This candidate is no longer accepted blindly. Rule-based judge v0 compares
`current` vs `baseexpand_candidate` and rejects the candidate when it changes
too many peaks, moves final omega too much, or increases the strict/final omega
spread too much.

Current conservative settings:
- max boundary extension: 0.010 min;
- max accepted area growth: 4%;
- minimum accepted area growth: 1%;
- raw support threshold: max(local noise * 0.20, apex height * 0.003).
- max accepted changed peaks: 2;
- max accepted omega shift: 0.050;
- max accepted strict/final spread increase: 0.045.

Audit output:
- each processed batch carries `judge_decisions_df`;
- regression output includes `judge_accepted`, `judge_rejected`, `judge_codes`,
  and `judge_reasons`;
- the GUI confidence popup shows accepted/rejected judge decisions for the
  current sample.

Current regression:
- previous production ALL MAE: 0.159349;
- with judge v0 ALL MAE: 0.157595;
- RMSE improved from 0.220224 to 0.218849;
- mean delta is effectively neutral: -0.000324.

Interpretation:
The visual diagnosis is valid: some peaks are cut above the base. However,
broadly expanding every peak is harmful and creates a positive omega bias. The
safe version must stay conservative and should eventually become one candidate
inside a per-cluster arbitration model rather than an unconditional large
expansion.
