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
