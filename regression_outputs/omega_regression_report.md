# Omega regression report

Generated with: `python omega_regression.py --data-dir . --out regression_outputs/omega_regression_current.xlsx`

Total evaluated samples: 286
Overall MAE: 0.184235
Overall RMSE: 0.262305
Overall max abs delta: 1.436689

## Summary

| scope | n | MAE | RMSE | mean_delta | median_abs | std_delta | within_0_2 | within_0_3 | within_0_4 | within_0_5 | within_0_6 | max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ALL | 286 | 0.184235 | 0.262305 | 0.058449 | 0.141806 | 0.255710 | 195 | 244 | 265 | 277 | 278 | 1.436689 |
| OLD45 | 45 | 0.164665 | 0.192390 | -0.007508 | 0.159049 | 0.192244 | 31 | 39 | 44 | 45 | 45 | 0.426502 |
| 01032026 | 4 | 0.523450 | 0.742707 | -0.355441 | 0.263310 | 0.652131 | 2 | 2 | 3 | 3 | 3 | 1.428173 |
| 02072026 | 75 | 0.170360 | 0.241696 | 0.091865 | 0.127291 | 0.223558 | 52 | 67 | 72 | 73 | 73 | 1.271281 |
| 03072026 | 75 | 0.204861 | 0.274619 | 0.083531 | 0.168970 | 0.261607 | 48 | 63 | 65 | 71 | 72 | 1.002719 |
| 06032026 | 15 | 0.259112 | 0.421432 | -0.035587 | 0.169599 | 0.419927 | 9 | 12 | 13 | 14 | 14 | 1.436689 |
| 09032026 | 6 | 0.199356 | 0.236387 | -0.101627 | 0.174575 | 0.213426 | 4 | 5 | 5 | 6 | 6 | 0.454813 |
| 13022026 | 11 | 0.143511 | 0.179485 | 0.118973 | 0.116004 | 0.134389 | 9 | 9 | 11 | 11 | 11 | 0.377590 |
| 13032026 | 13 | 0.142634 | 0.181305 | 0.005308 | 0.076126 | 0.181227 | 9 | 11 | 13 | 13 | 13 | 0.356720 |
| 14012026 | 18 | 0.166244 | 0.190334 | 0.043493 | 0.164229 | 0.185298 | 13 | 15 | 18 | 18 | 18 | 0.346602 |
| 16012026 | 17 | 0.169668 | 0.220203 | 0.062175 | 0.135054 | 0.211244 | 11 | 13 | 16 | 17 | 17 | 0.495057 |
| 18022026 | 9 | 0.150084 | 0.174717 | 0.098479 | 0.148469 | 0.144319 | 7 | 8 | 9 | 9 | 9 | 0.331283 |
| 20022026 | 7 | 0.112971 | 0.132301 | 0.045854 | 0.097620 | 0.124101 | 6 | 7 | 7 | 7 | 7 | 0.262190 |
| 20032026 | 14 | 0.183093 | 0.204628 | -0.084983 | 0.159997 | 0.186146 | 9 | 13 | 13 | 14 | 14 | 0.426502 |
| 23012026 | 7 | 0.070948 | 0.107332 | 0.050819 | 0.047603 | 0.094539 | 6 | 7 | 7 | 7 | 7 | 0.263465 |
| 27022026 | 15 | 0.190793 | 0.272323 | 0.169084 | 0.113817 | 0.213473 | 10 | 12 | 13 | 14 | 14 | 0.750894 |

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
| REJECT | overestimated_cluster | 55 |
| REJECT | underestimated_cluster | 16 |
| REJECT | overestimated_low_confidence | 10 |
| REJECT | underestimated_low_confidence | 7 |
| REJECT | overestimated_correction_spread | 3 |
| REJECT | underestimated_correction_spread | 3 |
| REJECT | underestimated_unclassified | 3 |
| REJECT | overestimated_unclassified | 2 |
| REVIEW | overestimated_unclassified | 76 |
| REVIEW | underestimated_unclassified | 58 |
| REVIEW | overestimated_cluster | 34 |
| REVIEW | overestimated_correction_spread | 9 |
| REVIEW | underestimated_cluster | 5 |
| REVIEW | underestimated_correction_spread | 5 |

## Safety judge

| metric | value |
| --- | --- |
| HIGH_RISK_GT_0_5 catch rate | 9/9 |
| HIGH_RISK_GT_0_5 review load | 16/286 |
| HIGH_RISK false positives | 7 |
| HIGH_RISK missed >0.5 | 0 |

### Safety judge by band

| safety_judge_band | n | MAE | max_abs | gt_0_5 | within_0_5 |
| --- | --- | --- | --- | --- | --- |
| HIGH_RISK_GT_0_5 | 16 | 0.625684 | 1.436689 | 9 | 7 |
| LOW_RISK | 270 | 0.158075 | 0.495057 | 0 | 270 |

## Diagnostic issue summary

| diagnostic_bucket | n | MAE | RMSE | max_abs | within_0_3 | within_0_5 |
| --- | --- | --- | --- | --- | --- | --- |
| under_c22_cluster | 3 | 1.378714 | 1.380810 | 1.436689 | 0 | 0 |
| over_c22_cluster | 6 | 0.748451 | 0.762880 | 1.002719 | 0 | 0 |
| watch_within_0_5 | 33 | 0.380665 | 0.386008 | 0.495057 | 0 | 33 |
| ok_within_0_3 | 244 | 0.129108 | 0.150804 | 0.299952 | 244 | 244 |

## Top diagnostic samples

| date | sample_no | sample_name | reference | calculated | delta | confidence | safety_judge_band | safety_judge_score | safety_judge_reasons | diagnostic_bucket | diagnostic_reasons |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 06032026 | 2 | O2_1110012176.D | 5.200000 | 3.763311 | -1.436689 | 70.000000 | HIGH_RISK_GT_0_5 | 94 | high_c22_ratio_low_dha_area,cluster_overlap_or_fit | under_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,underestimated_gt_0_5 |
| 01032026 | 4 | O4_5555839154.D | 5.300000 | 3.871827 | -1.428173 | 42.000000 | HIGH_RISK_GT_0_5 | 100 | high_c22_ratio_low_dha_area,baseline_fallback,low_or_medium_confidence,cluster_overlap_or_fit | under_c22_cluster | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,underestimated_gt_0_5 |
| 02072026 | 105067678302 | O4_105067678302.D | 4.930000 | 3.658719 | -1.271281 | 70.000000 | HIGH_RISK_GT_0_5 | 97 | c22_low_ratio_with_asymmetric_dha,baseline_fallback | under_c22_cluster | baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,wide_peak_window,asymmetric_peak_window,underestimated_gt_0_5 |
| 03072026 | 105069916298 | O27_925663916002.D | 2.500000 | 3.502719 | 1.002719 | 50.000000 | HIGH_RISK_GT_0_5 | 99 | high_c22_ratio_low_dha_area,low_or_medium_confidence,cluster_overlap_or_fit | over_c22_cluster | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,overestimated_gt_0_5 |
| 03072026 | 900512824605 | O8_104837397699.D | 3.000000 | 3.801493 | 0.801493 | 74.000000 | HIGH_RISK_GT_0_5 | 90 | high_c22_ratio_low_dha_area | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,overestimated_gt_0_5 |
| 03072026 | 1110012947 | O69_1110012953.D | 2.700000 | 3.478291 | 0.778291 | 60.000000 | HIGH_RISK_GT_0_5 | 94 | high_c22_ratio_low_dha_area,cluster_overlap_or_fit | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,wide_peak_window,asymmetric_peak_window,overestimated_gt_0_5 |
| 27022026 | 5 | O5_1100021851.D | 5.200000 | 5.950894 | 0.750894 | 45.000000 | HIGH_RISK_GT_0_5 | 91 | low_confidence_c20_shape_exception,low_or_medium_confidence,cluster_overlap_or_fit | over_c22_cluster | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,overestimated_gt_0_5 |
| 02072026 | 903862928199 | O71_903862928199.D | 3.950000 | 4.577730 | 0.627730 | 62.000000 | HIGH_RISK_GT_0_5 | 94 | high_c22_ratio_low_dha_area,cluster_overlap_or_fit | over_c22_cluster | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,wide_peak_window,asymmetric_peak_window,overestimated_gt_0_5 |
| 03072026 | 105068695803 | O10_929990895102.D | 3.020000 | 3.549577 | 0.529577 | 42.000000 | HIGH_RISK_GT_0_5 | 100 | high_c22_ratio_low_dha_area,baseline_fallback,low_or_medium_confidence,cluster_overlap_or_fit | over_c22_cluster | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,overestimated_gt_0_5 |
| 16012026 | 3 | O3_1100020357.D | 7.700000 | 7.204943 | -0.495057 | 65.500000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 02072026 | 105068232403 | O26_105068232403.D | 6.430000 | 6.921258 | 0.491258 | 67.000000 | LOW_RISK | 4 | cluster_overlap_or_fit | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 03072026 | 105070335002 | O38_925380252607.D | 7.330000 | 7.818360 | 0.488360 | 68.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 03072026 | 1110012917 | O72_1110012951.D | 3.400000 | 3.885498 | 0.485498 | 50.000000 | LOW_RISK | 14 | baseline_fallback,low_or_medium_confidence,cluster_overlap_or_fit | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 105070691099 | O51_105067906599.D | 6.340000 | 5.861877 | -0.478123 | 76.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 09032026 | 6 | O6_5573805910.D | 5.200000 | 4.745187 | -0.454813 | 90.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 104962405001 | O40_105070893202.D | 9.830000 | 9.380991 | -0.449009 | 60.000000 | LOW_RISK | 5 | baseline_fallback | watch_within_0_5 | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 1110012934 | O66_1110012914.D | 2.700000 | 3.137572 | 0.437572 | 64.000000 | LOW_RISK | 4 | cluster_overlap_or_fit | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 20032026 | 7 | O7_1100027106.D | 4.100000 | 3.673498 | -0.426502 | 85.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 06032026 | 9 | O9_1100026024.D | 7.500000 | 7.075867 | -0.424133 | 61.000000 | LOW_RISK | 5 | baseline_fallback | watch_within_0_5 | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 03072026 | 1110012955 | O71_1110012917.D | 4.200000 | 3.782969 | -0.417031 | 94.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c18_complex_boundaries,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 27022026 | 14 | O14_1100021859.D | 8.500000 | 8.909283 | 0.409283 | 54.000000 | LOW_RISK | 14 | baseline_fallback,low_or_medium_confidence,cluster_overlap_or_fit | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 27022026 | 8 | O8_1100021861.D | 5.100000 | 5.499491 | 0.399491 | 71.000000 | LOW_RISK | 4 | cluster_overlap_or_fit | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 02072026 | 104895240699 | O49_104895240699.D | 3.530000 | 3.925246 | 0.395246 | 72.000000 | LOW_RISK | 4 | cluster_overlap_or_fit | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,wide_peak_window,asymmetric_peak_window,within_clinical_band |
| 02072026 | 104714605399 | O30_104714605399.D | 2.460000 | 2.845791 | 0.385791 | 58.000000 | LOW_RISK | 9 | low_or_medium_confidence,cluster_overlap_or_fit | watch_within_0_5 | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 13022026 | 1 | O1_5551112961.D | 7.000000 | 7.377590 | 0.377590 | 69.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |

## Outliers > 0.5

Count: 9

| date | sample_no | instrument_no | sample_id | sample_name | match_method | reference | calculated | delta | confidence | selected_variant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 06032026 | 2 | 2 | 1110012176 | O2_1110012176.D | instrument_no | 5.200000 | 3.763311 | -1.436689 | 70.000000 | current_pipeline |
| 01032026 | 4 | 4 | 5555839154 | O4_5555839154.D | instrument_no | 5.300000 | 3.871827 | -1.428173 | 42.000000 | current_pipeline |
| 02072026 | 105067678302 | 4 | 105067678302 | O4_105067678302.D | sample_id | 4.930000 | 3.658719 | -1.271281 | 70.000000 | current_pipeline |
| 03072026 | 105069916298 | 27 | 925663916002 | O27_925663916002.D | position_date_override | 2.500000 | 3.502719 | 1.002719 | 50.000000 | current_pipeline |
| 03072026 | 900512824605 | 8 | 104837397699 | O8_104837397699.D | position_date_override | 3.000000 | 3.801493 | 0.801493 | 74.000000 | current_pipeline |
| 03072026 | 1110012947 | 69 | 1110012953 | O69_1110012953.D | position_date_override | 2.700000 | 3.478291 | 0.778291 | 60.000000 | current_pipeline |
| 27022026 | 5 | 5 | 1100021851 | O5_1100021851.D | instrument_no | 5.200000 | 5.950894 | 0.750894 | 45.000000 | current_pipeline |
| 02072026 | 903862928199 | 71 | 903862928199 | O71_903862928199.D | sample_id | 3.950000 | 4.577730 | 0.627730 | 62.000000 | current_pipeline |
| 03072026 | 105068695803 | 10 | 929990895102 | O10_929990895102.D | position_date_override | 3.020000 | 3.549577 | 0.529577 | 42.000000 | current_pipeline |