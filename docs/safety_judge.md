# Safety judge for large omega misses

This layer is a reference-free evaluator for the current engine. It does **not**
change the omega number. Its job is to catch samples that look like the known
large failures (`abs(delta) > 0.5`) before they are silently accepted.

## Training target

The current regression corpus has 286 evaluated samples and 9 samples outside the
clinical tolerance band `±0.5`. The judge was tuned against those failures using
only runtime features: C22/DHA/DPA integration geometry, C22 area ratios,
confidence, baseline mode, and existing cluster/refinement status.

## Current rule family

The high-risk class is intentionally narrow. It fires on repeated C22/C20 geometry
patterns seen in the scary samples:

1. Low DPA-to-C22:4 ratio with asymmetric DHA integration.
2. High DPA-to-C22:4 ratio combined with low DHA area and omega above the low-end
   guard threshold.
3. A low-confidence C20 shape exception observed in one remaining over-estimate.

Small supporting penalties are added for fallback baseline, low/medium confidence,
and existing cluster fit/overlap reasons. These supporting signals cannot create a
high-risk decision alone; they only raise confidence once a learned geometry pattern
has already fired.

## Current measured behavior

On the current corpus, `HIGH_RISK_GT_0_5` catches all 9 known `>0.5` misses while
marking 16/286 samples for focused review. That is intentionally much narrower than
the older general `review_flag`, which is conservative but too broad to be useful
as a triage queue.

The low-risk band currently has zero `>0.5` misses in this corpus, but this is not a
formal guarantee on new unseen data. Treat it as a practical triage layer until we
have enough independent future batches to validate it.
