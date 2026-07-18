# Quality-judge experiment — 2026-07-14

## Data correction

The original cloud dataset matched `03072026` by the printed sample ID. That is
invalid for this acquisition: after an empty first injection, the 75 non-empty
signals and the 76 workbook rows are offset by one. The resulting alternating
manual/calculated pattern created 64 artificial `>0.5` labels and invalidates
the neural-network metrics below as model-selection evidence.

After matching `03072026` by acquisition position, the 411 evaluated
chromatograms contain 32 real errors greater than 0.5: 8 in `02072026`, 10 in
`03072026`, and 14 in `06072026`. `14072026` remained sealed and was not
exported, trained on, or evaluated.

## Validation design

Each error-bearing batch was held out in turn together with several normal
batches. The network therefore had to learn from the other two faulty dates and
generalize to a completely unseen faulty date. This is stricter and more useful
than a random sample split.

## Superseded neural results

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
also checked with the same misaligned labels. These figures are retained only
to explain why the neural model was not integrated; they must not be compared
with future correctly aligned runs.

## Decision

Do not integrate the trained model into the GUI. The rule-based production judge
was recalibrated on the corrected labels with a hard 100% recall constraint. A
future neural experiment must first regenerate the dataset with the positional
`03072026` mapping and must not use the sealed `14072026` answers for feature
design or model selection.
