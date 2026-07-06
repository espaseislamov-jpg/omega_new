# Omega regression report

Generated with: `python omega_regression.py --data-dir . --out regression_outputs/omega_regression_current.xlsx`

Total evaluated samples: 286
Overall MAE: 0.192675
Overall RMSE: 0.271020
Overall max abs delta: 1.436689

## Summary

| scope | n | MAE | RMSE | mean_delta | median_abs | std_delta | within_0_2 | within_0_3 | within_0_4 | within_0_5 | within_0_6 | max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ALL | 286 | 0.192675 | 0.271020 | 0.062582 | 0.140797 | 0.263696 | 179 | 236 | 261 | 275 | 278 | 1.436689 |
| OLD45 | 45 | 0.189594 | 0.230734 | 0.049514 | 0.171403 | 0.225359 | 25 | 35 | 41 | 45 | 45 | 0.487020 |
| 01032026 | 4 | 0.542328 | 0.767982 | -0.443823 | 0.363310 | 0.626752 | 2 | 2 | 2 | 2 | 3 | 1.428173 |
| 02072026 | 75 | 0.175944 | 0.248407 | 0.094551 | 0.126947 | 0.229709 | 50 | 65 | 71 | 73 | 73 | 1.271281 |
| 03072026 | 75 | 0.201826 | 0.276635 | 0.078080 | 0.168388 | 0.265387 | 48 | 62 | 66 | 70 | 71 | 1.002719 |
| 06032026 | 15 | 0.224745 | 0.409470 | -0.112991 | 0.116387 | 0.393572 | 10 | 12 | 13 | 14 | 14 | 1.436689 |
| 09032026 | 6 | 0.178398 | 0.200993 | 0.006149 | 0.176312 | 0.200899 | 3 | 5 | 6 | 6 | 6 | 0.308157 |
| 13022026 | 11 | 0.203994 | 0.252392 | 0.130396 | 0.216004 | 0.216098 | 5 | 9 | 9 | 11 | 11 | 0.477590 |
| 13032026 | 13 | 0.158485 | 0.201482 | 0.086426 | 0.120046 | 0.182005 | 10 | 11 | 12 | 13 | 13 | 0.420672 |
| 14012026 | 18 | 0.207480 | 0.229020 | 0.112183 | 0.220086 | 0.199663 | 8 | 14 | 17 | 18 | 18 | 0.404052 |
| 16012026 | 17 | 0.174680 | 0.210700 | 0.066477 | 0.139156 | 0.199938 | 11 | 13 | 17 | 17 | 17 | 0.395057 |
| 18022026 | 9 | 0.185581 | 0.202088 | 0.030425 | 0.203797 | 0.199785 | 4 | 9 | 9 | 9 | 9 | 0.293322 |
| 20022026 | 7 | 0.124124 | 0.160613 | 0.031568 | 0.092139 | 0.157480 | 6 | 6 | 7 | 7 | 7 | 0.362190 |
| 20032026 | 14 | 0.195484 | 0.256887 | -0.065337 | 0.184920 | 0.248439 | 7 | 10 | 12 | 14 | 14 | 0.487020 |
| 23012026 | 7 | 0.149733 | 0.178849 | 0.109435 | 0.141044 | 0.141460 | 6 | 6 | 7 | 7 | 7 | 0.363465 |
| 27022026 | 15 | 0.188604 | 0.244956 | 0.155403 | 0.103298 | 0.189349 | 9 | 12 | 13 | 14 | 15 | 0.550894 |

## Input audit

| date | xlsx_path | csv_path | reference_rows | instrument_batches | matched_by_sample_id | reference_id_rows | position_rows | missing_reference_ids |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 02072026 | extracted_desktop/test_bigbatch_020726.xlsx | extracted_desktop/02072026.CSV | 75 | 75 | 75 | 75 | 0 | 0 |
| 03072026 | extracted_desktop/test_bigbatch_03072026.xlsx | extracted_desktop/03072026.CSV | 76 | 75 | 75 | 76 | 0 | 1 |
| 01032026 | test_bigbatch_01032026.xlsx | 01032026.CSV | 4 | 4 | 0 | 0 | 4 | 0 |
| 06032026 | test_bigbatch_06032026.xlsx | 06032026.CSV | 15 | 15 | 0 | 0 | 15 | 0 |
| 09032026 | test_bigbatch_09032026.xlsx | 09032026.CSV | 6 | 6 | 0 | 0 | 6 | 0 |
| 13022026 | test_bigbatch_13022026.xlsx | 13022026.CSV | 11 | 11 | 0 | 0 | 11 | 0 |
| 13032026 | test_bigbatch_13032026.xlsx | 13032026.CSV | 13 | 13 | 0 | 0 | 13 | 0 |
| 14012026 | test_bigbatch_14012026.xlsx | 14012026.CSV | 18 | 18 | 0 | 0 | 18 | 0 |
| 16012026 | test_bigbatch_16012026.xlsx | 16012026.CSV | 17 | 17 | 0 | 0 | 17 | 0 |
| 18022026 | test_bigbatch_18022026.xlsx | 18022026.CSV | 9 | 9 | 0 | 0 | 9 | 0 |
| 20022026 | test_bigbatch_20022026.xlsx | 20022026.CSV | 7 | 7 | 0 | 0 | 7 | 0 |
| 20032026 | test_bigbatch_20032026.xlsx | 20032026.CSV | 14 | 14 | 0 | 0 | 14 | 0 |
| 23012026 | test_bigbatch_23012026.xlsx | 23012026.CSV | 7 | 7 | 0 | 0 | 7 | 0 |
| 27022026 | test_bigbatch_27022026.xlsx | 27022026.CSV | 15 | 15 | 0 | 0 | 15 | 0 |

## Errors

| date | xlsx_path | csv_path | excel_row | sample_no | reference | match_method | error |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 03072026 | extracted_desktop/test_bigbatch_03072026.xlsx | extracted_desktop/03072026.CSV | 77 | 1110012956 | 3.800000 | missing_position_date_override | No matching instrument batch |

## Review / outlier classification

| review_flag | outlier_class | n |
| --- | --- | --- |
| OK | overestimated_unclassified | 78 |
| OK | underestimated_unclassified | 44 |
| REJECT | overestimated_cluster | 2 |
| REJECT | underestimated_cluster | 2 |
| REJECT | overestimated_low_confidence | 1 |
| REJECT | underestimated_low_confidence | 1 |
| REVIEW | overestimated_cluster | 86 |
| REVIEW | underestimated_cluster | 20 |
| REVIEW | overestimated_correction_spread | 14 |
| REVIEW | underestimated_unclassified | 12 |
| REVIEW | overestimated_low_confidence | 10 |
| REVIEW | underestimated_low_confidence | 8 |
| REVIEW | overestimated_unclassified | 4 |
| REVIEW | underestimated_correction_spread | 4 |

## Diagnostic issue summary

| diagnostic_bucket | n | MAE | RMSE | max_abs | within_0_3 | within_0_5 |
| --- | --- | --- | --- | --- | --- | --- |
| under_c22_cluster | 3 | 1.378714 | 1.380810 | 1.436689 | 0 | 0 |
| over_c22_cluster | 6 | 0.715117 | 0.733889 | 1.002719 | 0 | 0 |
| under_c20_cluster | 2 | 0.589309 | 0.592325 | 0.649009 | 0 | 0 |
| watch_within_0_5 | 39 | 0.384883 | 0.388877 | 0.491258 | 0 | 39 |
| ok_within_0_3 | 236 | 0.129192 | 0.152053 | 0.299952 | 236 | 236 |

## Top diagnostic samples

| date | sample_no | sample_name | reference | calculated | delta | confidence | diagnostic_bucket | diagnostic_reasons |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 06032026 | 2 | O2_1110012176.D | 5.200000 | 3.763311 | -1.436689 | 70.000000 | under_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,underestimated_gt_0_5 |
| 01032026 | 4 | O4_5555839154.D | 5.300000 | 3.871827 | -1.428173 | 42.000000 | under_c22_cluster | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,underestimated_gt_0_5 |
| 02072026 | 105067678302 | O4_105067678302.D | 4.930000 | 3.658719 | -1.271281 | 70.000000 | under_c22_cluster | baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,wide_peak_window,asymmetric_peak_window,underestimated_gt_0_5 |
| 03072026 | 105069916298 | O27_925663916002.D | 2.500000 | 3.502719 | 1.002719 | 50.000000 | over_c22_cluster | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,overestimated_gt_0_5 |
| 03072026 | 900512824605 | O8_104837397699.D | 3.000000 | 3.801493 | 0.801493 | 74.000000 | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,overestimated_gt_0_5 |
| 03072026 | 1110012947 | O69_1110012953.D | 2.700000 | 3.478291 | 0.778291 | 60.000000 | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,wide_peak_window,asymmetric_peak_window,overestimated_gt_0_5 |
| 03072026 | 104962405001 | O40_105070893202.D | 9.830000 | 9.180991 | -0.649009 | 60.000000 | under_c20_cluster | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,underestimated_gt_0_5 |
| 02072026 | 903862928199 | O71_903862928199.D | 3.950000 | 4.577730 | 0.627730 | 62.000000 | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,wide_peak_window,asymmetric_peak_window,overestimated_gt_0_5 |
| 27022026 | 5 | O5_1100021851.D | 5.200000 | 5.750894 | 0.550894 | 45.000000 | over_c22_cluster | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,overestimated_gt_0_5 |
| 01032026 | 3 | O3_1110012153.D | 12.500000 | 11.970391 | -0.529609 | 54.000000 | under_c20_cluster | low_or_medium_confidence,baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,underestimated_gt_0_5 |
| 03072026 | 105068695803 | O10_929990895102.D | 3.020000 | 3.549577 | 0.529577 | 42.000000 | over_c22_cluster | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,overestimated_gt_0_5 |
| 02072026 | 105068232403 | O26_105068232403.D | 6.430000 | 6.921258 | 0.491258 | 67.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 20032026 | 8 | O8_1100027096.D | 7.300000 | 6.812980 | -0.487020 | 69.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 1110012917 | O72_1110012951.D | 3.400000 | 3.885498 | 0.485498 | 50.000000 | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 105070691099 | O51_105067906599.D | 6.340000 | 5.856915 | -0.483085 | 60.000000 | watch_within_0_5 | baseline_fallback,c22_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 13022026 | 1 | O1_5551112961.D | 7.000000 | 7.477590 | 0.477590 | 69.000000 | watch_within_0_5 | c20_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 20032026 | 7 | O7_1100027106.D | 4.100000 | 3.640605 | -0.459395 | 95.000000 | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 1110012934 | O66_1110012914.D | 2.700000 | 3.137572 | 0.437572 | 64.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 02072026 | 903982571001 | O24_903982571001.D | 8.020000 | 7.590799 | -0.429201 | 49.000000 | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 13032026 | 8 | O8_1100026594.D | 5.900000 | 6.320672 | 0.420672 | 69.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 13022026 | 9 | O9_1100020744.D | 5.200000 | 5.616856 | 0.416856 | 64.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 27022026 | 14 | O14_1100021859.D | 8.500000 | 8.909283 | 0.409283 | 54.000000 | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 105070893202 | O41_900917265501.D | 7.800000 | 7.391116 | -0.408884 | 60.000000 | watch_within_0_5 | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 06032026 | 15 | O15_1100026030.D | 7.500000 | 7.904240 | 0.404240 | 69.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 14012026 | 12 | O12_1100020191.D | 4.800000 | 5.204052 | 0.404052 | 66.000000 | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |

## Outliers > 0.5

Count: 11

| date | sample_no | instrument_no | sample_id | sample_name | match_method | reference | calculated | delta | confidence | selected_variant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 06032026 | 2 | 2 | 1110012176 | O2_1110012176.D | instrument_no | 5.200000 | 3.763311 | -1.436689 | 70.000000 | current_pipeline |
| 01032026 | 4 | 4 | 5555839154 | O4_5555839154.D | instrument_no | 5.300000 | 3.871827 | -1.428173 | 42.000000 | current_pipeline |
| 02072026 | 105067678302 | 4 | 105067678302 | O4_105067678302.D | sample_id | 4.930000 | 3.658719 | -1.271281 | 70.000000 | current_pipeline |
| 03072026 | 105069916298 | 27 | 925663916002 | O27_925663916002.D | position_date_override | 2.500000 | 3.502719 | 1.002719 | 50.000000 | current_pipeline |
| 03072026 | 900512824605 | 8 | 104837397699 | O8_104837397699.D | position_date_override | 3.000000 | 3.801493 | 0.801493 | 74.000000 | current_pipeline |
| 03072026 | 1110012947 | 69 | 1110012953 | O69_1110012953.D | position_date_override | 2.700000 | 3.478291 | 0.778291 | 60.000000 | current_pipeline |
| 03072026 | 104962405001 | 40 | 105070893202 | O40_105070893202.D | position_date_override | 9.830000 | 9.180991 | -0.649009 | 60.000000 | current_pipeline |
| 02072026 | 903862928199 | 71 | 903862928199 | O71_903862928199.D | sample_id | 3.950000 | 4.577730 | 0.627730 | 62.000000 | current_pipeline |
| 27022026 | 5 | 5 | 1100021851 | O5_1100021851.D | instrument_no | 5.200000 | 5.750894 | 0.550894 | 45.000000 | current_pipeline |
| 01032026 | 3 | 3 | 1110012153 | O3_1110012153.D | instrument_no | 12.500000 | 11.970391 | -0.529609 | 54.000000 | current_pipeline |
| 03072026 | 105068695803 | 10 | 929990895102 | O10_929990895102.D | position_date_override | 3.020000 | 3.549577 | 0.529577 | 42.000000 | current_pipeline |