# Safety judge for errors above 0.5

The production judge is reference-free and does not change the calculated omega
value. Its only hard target is to send every known integration error above 0.5
percentage point to manual review.

## Current validation set

- 411 evaluated chromatograms from 17 non-sealed batch dates;
- 28 errors above 0.5 after the bounded legacy C22 retry;
- `14072026` remains sealed and was not used for thresholds or validation.

## Rules

Structural failures (missing key peaks or one signal assigned to several peaks)
produce a stop warning. Three conservative C20/C22 geometry envelopes produce a
manual-review warning. They use only peak positions, widths, asymmetry, and the
existing C22 overlap state available in the GUI. The user-facing message names
the peaks to inspect; numeric thresholds stay in `omega_core/metrics.py`.

The smoothed baseline-corrected height ratio `C22:5 / C22:4` is an additional
review-order hint. On the validated legacy profile, a high ratio moves C22 to
the front of an already existing warning. It never creates a warning by itself,
never removes a warning, and never changes the calculated omega value.

Custom instrument profiles do not inherit the legacy judge calibration. Until
a profile-specific manual validation is stored in the exported profile, every
result is marked as requiring review instead of being shown as green.

The second instrument profile currently has an explicit emergency judge mode
with a working target of 0.8 percentage point. It replaces the blanket warning
with profile-relative gates for selected numerator peaks, duplicate/missing
assignments, RT residual, integration bounds, recovery status, and broad C20/C22
height-ratio envelopes. This mode is intentionally labelled emergency and is
not treated as a calibrated profile.

The old general confidence penalties remain as context, but low confidence alone
does not create a high-error warning.

## Measured behavior

On the corrected historical set, the high-risk bands catch all 29 errors present
before the bounded legacy C22 retry. That retry fixes one case without introducing
a new error above 0.5, leaving 28/411 final errors. The judge still marks 109/411
samples for review; this false-warning rate is the deliberate cost of the requested
100% historical recall before retry.

This is a measured historical result, not a mathematical guarantee for unseen
batches. New manually checked dates should be appended to the regression set and
the same 100% recall audit rerun before any threshold is narrowed.
