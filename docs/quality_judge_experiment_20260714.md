# Quality-judge experiment — 2026-07-14

## Dataset

- 411 matched labeled chromatograms from 17 complete batch dates.
- 152 samples with absolute error greater than 0.3.
- 86 samples with absolute error greater than 0.5.
- The greater-than-0.5 errors occur in only three batches: 8 in `02072026`,
  64 in `03072026`, and 14 in `06072026`.
- `14072026` remained sealed and was not exported, trained on, or evaluated.

## Validation design

Each error-bearing batch was held out in turn together with several normal
batches. The network therefore had to learn from the other two faulty dates and
generalize to a completely unseen faulty date. This is stricter and more useful
than a random sample split.

## Results

The small MLP using the final integration geometry achieved an aggregate ROC-AUC
of 0.691 for errors greater than 0.5. At an operating point with 50% precision,
it detected 25 of 86 errors and reviewed 50 of 411 samples. It detected none of
the 8 errors in the held-out `02072026` batch.

Adding disagreement across five baseline variants did not improve transfer:

- aggregate ROC-AUC: 0.580;
- average precision: 0.245;
- held-out batch ROC-AUC: 0.632 (`02072026`), 0.388 (`03072026`), and 0.599
  (`06072026`);
- the high-recall threshold detected 83 of 86 errors but marked 314 of 411
  samples (76.4%) for review.

Simple logistic regression, Extra Trees, and histogram gradient boosting were
also checked with the same held-out batches and did not transfer better. This
indicates a data/feature limitation rather than a need for a larger network.

## Decision

Do not integrate the trained model into the GUI. The current labeled set has
enough individual errors but not enough independent error-bearing batches. A
future experiment should add several new fully labeled dates containing both
good and bad integrations, then rerun the existing cloud workflow without using
the sealed `14072026` answers for feature design or model selection.
