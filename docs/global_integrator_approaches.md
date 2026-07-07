# Global integrator improvement approaches

The current evidence says integration boundaries, not just omega post-correction,
are the main remaining problem. The strongest user observation is that after a few
samples in a batch, target RTs stabilize to the thousandths. That should be used as
a batch-level prior, but not as a hard midpoint clamp.

## Approach 1 — Batch RT/boundary profile

Build a robust profile from the whole batch: median RT, median left/right width,
robust scale, and chi-square-like stability score for each dangerous target. The
profile is saved as JSON and can become the zero/initial prior for the next run.

This is implemented as a read-only experiment in `scripts/build_batch_profile.py` and
`omega_core.batch_profile`. It currently profiles the two latest batches and a set
of old problem batches.

## Approach 2 — Soft RT-prior integrator

Instead of clipping to hard midpoints, use the profile as a soft prior when choosing
candidate boundaries. A boundary candidate should score well only if it satisfies
several signals at once:

- it stays near the learned RT/width distribution;
- it lands on a local valley or clear derivative transition;
- it does not create a large omega jump;
- it conserves cluster area within a plausible range;
- it improves the safety judge / confidence signal.

This should replace universal hard clamping, which worsened MAE in the first
corridor experiment.

## Approach 3 — Iterative batch convergence

For each batch, run the current pipeline once, build/update the profile, then rerun
high-risk samples with the profile. Repeat until the RT/width chi-square-like score
stops improving or a small iteration cap is reached. The profile JSON can then seed
the next launch as a warm start, while still allowing per-batch drift.

The key safety rule: the profile may propose priors and candidate windows, but the
final boundary should still require local signal evidence. Stable RT tells us where
a peak should be; it does not by itself prove where the integration baseline ends.
