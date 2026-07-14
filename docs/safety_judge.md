# Safety judge for errors above 0.5

The production judge is reference-free and does not change the calculated omega
value. Its only hard target is to send every known integration error above 0.5
percentage point to manual review.

## Current validation set

- 411 evaluated chromatograms from 17 non-sealed batch dates;
- 35 errors above 0.5 with the current integration engine;
- `14072026` remains sealed and was not used for thresholds or validation.

## Rules

Structural failures (missing key peaks or one signal assigned to several peaks)
produce a stop warning. Three conservative C20/C22 geometry envelopes produce a
manual-review warning. They use only peak positions, widths, asymmetry, and the
existing C22 overlap state available in the GUI. The user-facing message names
the peaks to inspect; numeric thresholds stay in `omega_core/metrics.py`.

The old general confidence penalties remain as context, but low confidence alone
does not create a high-error warning.

## Measured behavior

On the corrected historical set, the high-risk bands catch 35/35 known errors
above 0.5. They mark 106/411 samples for review, including 71/376 samples that are
inside the tolerance band (18.9% of normal samples). This false-warning rate is
the deliberate cost of the requested 100% historical recall.

This is a measured historical result, not a mathematical guarantee for unseen
batches. New manually checked dates should be appended to the regression set and
the same 100% recall audit rerun before any threshold is narrowed.
