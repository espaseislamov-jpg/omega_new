# New archive data diagnostics

Generated after extracting `Desktop.part1.rar`/`Desktop.part2.rar` with `unar` into `extracted_desktop/` and running the current regression with the diagnostic layer.

## Input coverage

| date | reference_rows | instrument_batches | matched_by_sample_id | missing_reference_ids |
| --- | --- | --- | --- | --- |
| 02072026 | 75 | 75 | 75 | 0 |
| 03072026 | 76 | 75 | 75 | 1 |

## Current-engine summary

| date | n | MAE | RMSE | median_abs | max_abs | within_0_3 | within_0_5 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 02072026 | 75 | 0.170360 | 0.241696 | 0.127291 | 1.271281 | 67 | 73 |
| 03072026 | 75 | 0.204861 | 0.274619 | 0.168970 | 1.002719 | 63 | 71 |

## Mapping notes

- `02072026` is fully extracted and matched by sample ID: 75 manual rows, 75 instrument batches, 75 sample-id matches.
- `03072026` is intentionally matched by row position rather than sample ID because this workbook is shifted: the first manual row corresponds to the first instrument batch `O2`, and the final manual row has no matching instrument batch.

## Algorithm note

The full corpus now has `ALL` MAE `0.1842`, `244/286` within `±0.3`, `277/286` within `±0.5`, and max absolute error `1.4367`. This run includes the bounded C22 width-balance calibration.

## July diagnostic issue summary

| diagnostic_bucket | n | MAE | max_abs | within_0_3 | within_0_5 |
| --- | --- | --- | --- | --- | --- |
| under_c22_cluster | 1 | 1.271281 | 1.271281 | 0 | 0 |
| over_c22_cluster | 5 | 0.747962 | 1.002719 | 0 | 0 |
| watch_within_0_5 | 14 | 0.403076 | 0.491258 | 0 | 14 |
| ok_within_0_3 | 130 | 0.134518 | 0.299952 | 130 | 130 |

## Errors / skipped rows

| date | sample_no | reference | match_method | error |
| --- | --- | --- | --- | --- |
| 03072026 | 1110012956 | 3.800000 | missing_position_date_override | No matching instrument batch |

## Top 20 July diagnostic samples

| date | sample_name | reference | calculated | delta | confidence | diagnostic_bucket | diagnostic_reasons |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 02072026 | O4_105067678302.D | 4.930000 | 3.658719 | -1.271281 | 70.000000 | under_c22_cluster | baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,wide_peak_window,asymmetric_peak_window,underestimated_gt_0_5 |
| 03072026 | O27_925663916002.D | 2.500000 | 3.502719 | 1.002719 | 50.000000 | over_c22_cluster | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,overestimated_gt_0_5 |
| 03072026 | O8_104837397699.D | 3.000000 | 3.801493 | 0.801493 | 74.000000 | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,overestimated_gt_0_5 |
| 03072026 | O69_1110012953.D | 2.700000 | 3.478291 | 0.778291 | 60.000000 | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,wide_peak_window,asymmetric_peak_window,overestimated_gt_0_5 |
| 02072026 | O71_903862928199.D | 3.950000 | 4.577730 | 0.627730 | 62.000000 | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,wide_peak_window,asymmetric_peak_window,overestimated_gt_0_5 |
| 03072026 | O10_929990895102.D | 3.020000 | 3.549577 | 0.529577 | 42.000000 | over_c22_cluster | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,overestimated_gt_0_5 |
| 02072026 | O26_105068232403.D | 6.430000 | 6.921258 | 0.491258 | 67.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 03072026 | O38_925380252607.D | 7.330000 | 7.818360 | 0.488360 | 68.000000 | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 03072026 | O72_1110012951.D | 3.400000 | 3.885498 | 0.485498 | 50.000000 | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | O51_105067906599.D | 6.340000 | 5.861877 | -0.478123 | 76.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 03072026 | O40_105070893202.D | 9.830000 | 9.380991 | -0.449009 | 60.000000 | watch_within_0_5 | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | O66_1110012914.D | 2.700000 | 3.137572 | 0.437572 | 64.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | O71_1110012917.D | 4.200000 | 3.782969 | -0.417031 | 94.000000 | watch_within_0_5 | c18_complex_boundaries,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 02072026 | O49_104895240699.D | 3.530000 | 3.925246 | 0.395246 | 72.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,wide_peak_window,asymmetric_peak_window,within_clinical_band |
| 02072026 | O30_104714605399.D | 2.460000 | 2.845791 | 0.385791 | 58.000000 | watch_within_0_5 | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 02072026 | O19_105066911701.D | 2.810000 | 3.179492 | 0.369492 | 62.000000 | watch_within_0_5 | c22_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,within_clinical_band |
| 02072026 | O67_1110012961.D | 3.200000 | 3.517048 | 0.317048 | 62.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | O39_104962405001.D | 6.630000 | 6.943750 | 0.313750 | 60.000000 | watch_within_0_5 | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | O30_910130574901.D | 5.860000 | 5.551396 | -0.308604 | 68.000000 | watch_within_0_5 | baseline_fallback,c22_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 02072026 | O15_922730889202.D | 4.240000 | 4.546285 | 0.306285 | 59.000000 | watch_within_0_5 | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,within_clinical_band |
