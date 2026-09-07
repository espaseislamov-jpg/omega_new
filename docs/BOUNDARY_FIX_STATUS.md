# Boundary continuity patch — 2026-09-07

Base: main 7ab0c5c129b7acdd0493a304ee4d94281b455e75.

Implemented:
- Empty primary detections and all-area-rejected results now reach existing targeted/OpenMS recovery.
- Supplemental candidates cannot displace a primary candidate merely by having an earlier RT within 0.006 min; primary detections retain their geometry.
- Exact same-instrument profile copies and import/export preserve processing mode and judge metadata. New profiles remain uncalibrated. Changes to calculation inputs invalidate calibration and emergency mode.
- Unedited four-decimal RT display cells preserve the original full precision.
- Painted footprint no longer extends past the numerical integration interval.
- Cluster refinement records per-stage changes to RT, bounds, area and status in result['boundary_history']. The GUI button exports this automatic-stage history, not a complete manual-edit or baseline-candidate history.

Verification:
- 42 unittest cases, including a synthetic bright narrow Gaussian at 7.650 min: old detector returns no peak, patched detector returns area 501.3256549262.
- Real CSV 25032026, two samples: old/new Omega values identical at 4.708780275475674 and 5.308608797634255 in the available runtime.
- These two samples are a non-regression smoke test, NOT manual-reference accuracy validation. Optional-library environments and Windows GUI/build have not been validated.

Not completed:
- Full paired CSV/workbook regression across all dates and instruments, including missing bright target assignments and final boundary errors.
- Validated replacement of the legacy relative-height and cluster boundary rules.
- Complete provenance across rejected baseline alternatives, joint retries and manual edits.
- Background GUI execution, comprehensive profiling, Windows installer verification.
- Any changes to judge training. Corrected historical acquisition-label alignment is accepted; this patch does not retrain or retune the judge.

Do not treat this branch as a fully validated production release. Keep main unchanged pending broader regression.
