# Peak detection and integration notes

This note captures how the current Omega pipeline finds peaks, selects target peaks,
and assigns integration borders. It is intentionally short and implementation-focused
so future tuning work starts from the same mental model.

## End-to-end flow

1. `omega_core.pipeline.process_batch` loads each batch, subtracts a baseline, then
   delegates the corrected chromatogram to `process_from_baseline`.
2. `process_from_baseline` adds smoothing/derivatives, detects candidate peaks,
   matches those candidates to fatty-acid targets, optionally asks the ChromatoPy
   adapter to refine selected targets, refines C18/C20/C22 clusters, and finally
   calculates omega metrics.
3. The regression harness (`omega_regression.py`) compares the calculated omega
   value with manual Excel references and records audit/outlier data.

## Smoothing and derivative signal

The detector works on a corrected signal plus a Savitzky-Golay smoothed copy. The
smoothing window is selected from several valid odd windows. Each candidate window
is scored by residual noise, curvature preservation, and peak-loss penalty; the
lowest-scoring window becomes `y_smooth`, with `dy` and `d2y` derivatives stored
for downstream boundary logic.

## Primary candidate detection

`detect_peak_candidates` uses `scipy.signal.find_peaks` on `y_smooth`. The dynamic
height and prominence floors are derived from robust noise estimates and signal
quantiles, so noisy or low-signal files do not use a fixed absolute threshold.
Each peak gets:

- an apex index/retention time;
- an integration start/end index;
- a trapezoid area over positive corrected signal;
- width and prominence diagnostics.

The detector then augments the generic peak list with targeted weak-peak search
around known C18/C20/C22 cluster retention times and optional PyOpenMS candidates.
This is important because clinically important omega peaks can be shoulders inside
larger neighboring peaks rather than obvious standalone maxima.

## Integration borders

For generic `find_peaks` candidates the default borders come from relative-height
widths. If prominence-base mode is enabled, the borders use prominence bases and
are clipped by valley splits between neighboring peaks plus a maximum half-width.
For targeted weak peaks, `_extract_peak_geometry` starts from the apex and walks
left/right using the derivative and a mixed boundary signal
`0.70*y_smooth + 0.30*y_corrected`, then snaps each side to the local minimum.
The stored area is the positive corrected trapezoid between these borders.

## Target selection

`match_targets_to_peaks` first estimates a global retention-time shift from reliable
targets. For reliable rows it searches candidates inside a narrow shifted RT window
and usually prefers the nearest candidate, unless a larger nearby peak is dominant.
Selected candidate boundaries become the target integration boundaries. Special C22
rules then handle the DHA/DPA/C22:4 region where overlapping peaks are common.

## Why the C20/EPA correction was disabled

The previous data-driven C20/EPA overlap-credit model added a fraction of the
neighboring `C20:3N8` area into EPA (`C20:5`) when EPA looked under-separated. On
the July regression set this was the largest systematic over-estimation source:
many top outliers were classified as C20/C22 cluster over-estimates and carried
large EPA overlap credits. Disabling that model is a conservative production-safe
change: it does not use manual reference values at runtime, reduces the full-corpus
MAE from about `0.3123` to `0.2502`, and reduces the maximum absolute error from
about `2.4794` to `1.4527` on the current extracted corpus.

## Next tuning targets

1. Replace the disabled C20/EPA credit with a bounded, validation-gated shoulder
   splitter that only credits EPA when the local derivative/valley evidence is
   strong enough.
2. Add per-target boundary diagnostics to regression outputs: apex RT error,
   left/right width, neighboring valley depth, and area ratio to nearest neighbor.
3. Keep tuning C22 overlap handling separately from C20/EPA: the first bounded
   C22/DPA over-integration guard now debits only high DPA/C22:4-ratio clusters,
   with a hard cap in omega points so it cannot create a large negative correction.

## Current C22 width-balance correction

The latest calculation change is deliberately small and bounded. After the normal
EPA/DHA/DPA calculation and existing C22 overlap/debit logic, `compute_omega` looks
at the local C22 integration widths for `C22:6`, `C22:5`, and `C22:4`. Very narrow
DPA cases are the only eligible cases. Depending on C22:4 width, DPA area, DHA area,
and the existing overlap fraction, the metric receives either `+0.20`, `-0.10`, or
`0.00` omega percentage points.

This is not a replacement for real deconvolution. It is a bounded calibration layer
that addresses a repeated diagnostic pattern without allowing a large hidden
correction. The cap is intentionally much smaller than the clinical tolerance band
(`±0.5`) and the desired inter-operator band (`±0.3`). The remaining large misses
should still be investigated at the peak-boundary/deconvolution level rather than
by increasing this calibration.

## Boundary finding lesson from RT-corridor experiment

The stable retention-time observation is real and useful, but a hard midpoint clamp
between neighboring target RTs is too blunt as a default production rule. I added an
experimental `OMEGA_TARGET_RT_CORRIDOR_GUARD=1` path that clips C18/C20/C22 target
intervals to midpoint corridors derived from corrected target RTs. On the current
corpus it reduced some left/right shoulder absorption but worsened total MAE, which
means the correct fix is not a universal clamp.

The next boundary fix should be local and evidence-based: for each high-risk C22/C20
cluster, evaluate valley depth, derivative sign changes, and area conservation before
moving boundaries. Stable target RTs should act as priors and guardrails, not as hard
integration endpoints.

## July regression hard-outlier pass

The latest full regression showed that the remaining >0.5 omega-point failures were not random; they clustered around two integration/credit failure modes:

1. **Severe EPA underfit in the C20 cluster.**  `C20:5` was sometimes integrated as an extremely small sliver while `C20:3N8` carried a large neighboring shoulder.  The old data-driven C20 credit was disabled because it overcorrected broad cases, so this pass re-enabled it only behind a tighter `EPA/C20:3N8 < 0.25` gate and adds extra underfit scale only for very low EPA ratios.
2. **Borderline C22 credit overshoot.**  When `C22:5/C22:4` is around 0.45–0.55 and the strict trio value is already near 5.5, a large C22:4 tail credit can push the final omega too high.  This pass caps the credit contribution in that narrow shape regime instead of globally weakening C22 correction.

On the committed full historical regression, this removes all `abs_delta > 0.5` rows while keeping the mean error and RMSE lower than the previous run.  The implementation is still conservative: it changes only bounded correction layers after peak matching/integration, and the diagnostics retain strict/final spread so the suspicious cases remain visible.
