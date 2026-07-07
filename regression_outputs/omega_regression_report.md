# Omega regression report

Generated with: `python omega_regression.py --data-dir . --out regression_outputs/omega_regression_current.xlsx`

Total evaluated samples: 136
Overall MAE: 0.164266
Overall RMSE: 0.200762
Overall max abs delta: 0.497193

## Summary

| scope | n | MAE | RMSE | mean_delta | median_abs | std_delta | within_0_2 | within_0_3 | within_0_4 | within_0_5 | within_0_6 | max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ALL | 136 | 0.164266 | 0.200762 | 0.039956 | 0.145891 | 0.196745 | 96 | 115 | 130 | 136 | 136 | 0.497193 |
| OLD45 | 45 | 0.161213 | 0.184509 | -0.001961 | 0.160945 | 0.184499 | 32 | 41 | 45 | 45 | 45 | 0.356720 |
| 01032026 | 4 | 0.193276 | 0.211146 | -0.025268 | 0.168009 | 0.209628 | 3 | 3 | 4 | 4 | 4 | 0.329609 |
| 06032026 | 15 | 0.203635 | 0.240343 | 0.042566 | 0.199209 | 0.236543 | 8 | 11 | 13 | 15 | 15 | 0.435115 |
| 09032026 | 6 | 0.251329 | 0.276038 | -0.153601 | 0.204294 | 0.229356 | 3 | 4 | 5 | 6 | 6 | 0.454813 |
| 13022026 | 11 | 0.158042 | 0.188513 | 0.133504 | 0.134865 | 0.133093 | 9 | 9 | 11 | 11 | 11 | 0.377590 |
| 13032026 | 13 | 0.137209 | 0.172568 | -0.000014 | 0.076126 | 0.172568 | 9 | 12 | 13 | 13 | 13 | 0.356720 |
| 14012026 | 18 | 0.167639 | 0.191387 | 0.044889 | 0.164229 | 0.186049 | 13 | 15 | 18 | 18 | 18 | 0.346602 |
| 16012026 | 17 | 0.169668 | 0.220203 | 0.062175 | 0.135054 | 0.211244 | 11 | 13 | 16 | 17 | 17 | 0.495057 |
| 18022026 | 9 | 0.150087 | 0.174720 | 0.098482 | 0.148496 | 0.144320 | 7 | 8 | 9 | 9 | 9 | 0.331283 |
| 20022026 | 7 | 0.112971 | 0.132301 | 0.045854 | 0.097620 | 0.124101 | 6 | 7 | 7 | 7 | 7 | 0.262190 |
| 20032026 | 14 | 0.175241 | 0.186223 | -0.064005 | 0.161617 | 0.174878 | 10 | 14 | 14 | 14 | 14 | 0.295864 |
| 23012026 | 7 | 0.070907 | 0.107310 | 0.034015 | 0.047603 | 0.101776 | 6 | 7 | 7 | 7 | 7 | 0.263465 |
| 27022026 | 15 | 0.165945 | 0.224270 | 0.129033 | 0.113817 | 0.183432 | 11 | 12 | 13 | 15 | 15 | 0.497193 |

## Input audit

| date | xlsx_path | csv_path | reference_rows | instrument_batches | matched_by_sample_id | reference_id_rows | position_rows | missing_reference_ids |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 01032026 | data/regression/test_bigbatch_01032026.xlsx | data/regression/01032026.CSV | 4 | 4 | 0 | 0 | 4 | 0 |
| 06032026 | data/regression/test_bigbatch_06032026.xlsx | data/regression/06032026.CSV | 15 | 15 | 0 | 0 | 15 | 0 |
| 09032026 | data/regression/test_bigbatch_09032026.xlsx | data/regression/09032026.CSV | 6 | 6 | 0 | 0 | 6 | 0 |
| 13022026 | data/regression/test_bigbatch_13022026.xlsx | data/regression/13022026.CSV | 11 | 11 | 0 | 0 | 11 | 0 |
| 13032026 | data/regression/test_bigbatch_13032026.xlsx | data/regression/13032026.CSV | 13 | 13 | 0 | 0 | 13 | 0 |
| 14012026 | data/regression/test_bigbatch_14012026.xlsx | data/regression/14012026.CSV | 18 | 18 | 0 | 0 | 18 | 0 |
| 16012026 | data/regression/test_bigbatch_16012026.xlsx | data/regression/16012026.CSV | 17 | 17 | 0 | 0 | 17 | 0 |
| 18022026 | data/regression/test_bigbatch_18022026.xlsx | data/regression/18022026.CSV | 9 | 9 | 0 | 0 | 9 | 0 |
| 20022026 | data/regression/test_bigbatch_20022026.xlsx | data/regression/20022026.CSV | 7 | 7 | 0 | 0 | 7 | 0 |
| 20032026 | data/regression/test_bigbatch_20032026.xlsx | data/regression/20032026.CSV | 14 | 14 | 0 | 0 | 14 | 0 |
| 23012026 | data/regression/test_bigbatch_23012026.xlsx | data/regression/23012026.CSV | 7 | 7 | 0 | 0 | 7 | 0 |
| 27022026 | data/regression/test_bigbatch_27022026.xlsx | data/regression/27022026.CSV | 15 | 15 | 0 | 0 | 15 | 0 |

## Errors

_empty_

## Review / outlier classification

| review_flag | outlier_class | n |
| --- | --- | --- |
| REJECT | underestimated_low_confidence | 7 |
| REJECT | underestimated_cluster | 6 |
| REJECT | overestimated_cluster | 5 |
| REJECT | overestimated_low_confidence | 5 |
| REJECT | overestimated_correction_spread | 3 |
| REJECT | underestimated_unclassified | 3 |
| REJECT | underestimated_correction_spread | 2 |
| REJECT | overestimated_unclassified | 1 |
| REVIEW | overestimated_unclassified | 56 |
| REVIEW | underestimated_unclassified | 31 |
| REVIEW | overestimated_correction_spread | 7 |
| REVIEW | overestimated_cluster | 5 |
| REVIEW | underestimated_correction_spread | 5 |

## Safety judge

| metric | value |
| --- | --- |
| HIGH_RISK_GT_0_5 catch rate | 0/0 |
| HIGH_RISK_GT_0_5 review load | 3/136 |
| HIGH_RISK false positives | 3 |
| HIGH_RISK missed >0.5 | 0 |

### Safety judge by band

| safety_judge_band | n | MAE | max_abs | gt_0_5 | within_0_5 |
| --- | --- | --- | --- | --- | --- |
| LOW_RISK | 133 | 0.163210 | 0.497193 | 0 | 133 |
| HIGH_RISK_GT_0_5 | 3 | 0.211062 | 0.435115 | 0 | 3 |

## Diagnostic issue summary

| diagnostic_bucket | n | MAE | RMSE | max_abs | within_0_3 | within_0_5 |
| --- | --- | --- | --- | --- | --- | --- |
| watch_within_0_5 | 21 | 0.370251 | 0.375052 | 0.497193 | 0 | 21 |
| ok_within_0_3 | 115 | 0.126651 | 0.148253 | 0.296175 | 115 | 115 |

## Top diagnostic samples

| date | sample_no | sample_name | reference | calculated | delta | confidence | safety_judge_band | safety_judge_score | safety_judge_reasons | diagnostic_bucket | diagnostic_reasons |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 27022026 | 5 | O5_1100021851.D | 5.200000 | 5.697193 | 0.497193 | 49.000000 | LOW_RISK | 9 | low_or_medium_confidence,cluster_overlap_or_fit | watch_within_0_5 | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 16012026 | 3 | O3_1100020357.D | 7.700000 | 7.204943 | -0.495057 | 65.500000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 09032026 | 6 | O6_5573805910.D | 5.200000 | 4.745187 | -0.454813 | 90.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 06032026 | 2 | O2_1110012176.D | 5.200000 | 4.764885 | -0.435115 | 57.000000 | HIGH_RISK_GT_0_5 | 99 | high_c22_ratio_low_dha_area,low_or_medium_confidence,cluster_overlap_or_fit | watch_within_0_5 | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,within_clinical_band |
| 06032026 | 9 | O9_1100026024.D | 7.500000 | 7.075867 | -0.424133 | 61.000000 | LOW_RISK | 5 | baseline_fallback | watch_within_0_5 | baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 27022026 | 14 | O14_1100021859.D | 8.500000 | 8.909283 | 0.409283 | 54.000000 | LOW_RISK | 14 | baseline_fallback,low_or_medium_confidence,cluster_overlap_or_fit | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 27022026 | 8 | O8_1100021861.D | 5.100000 | 5.499491 | 0.399491 | 71.000000 | LOW_RISK | 4 | cluster_overlap_or_fit | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 13022026 | 1 | O1_5551112961.D | 7.000000 | 7.377590 | 0.377590 | 69.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 16012026 | 17 | O17_1100020383.D | 4.900000 | 5.264508 | 0.364508 | 75.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 13032026 | 2 | O2_1100026588.D | 7.300000 | 6.943280 | -0.356720 | 69.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 09032026 | 3 | O3_1100026395.D | 5.800000 | 5.448612 | -0.351388 | 86.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 14012026 | 2 | O2_5557347805.D | 6.500000 | 6.153398 | -0.346602 | 82.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 16012026 | 4 | O4_1100020363.D | 5.600000 | 5.937273 | 0.337273 | 68.150000 | LOW_RISK | 0 |  | watch_within_0_5 | missing_target,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,within_clinical_band |
| 18022026 | 5 | O5_1100020799.D | 8.000000 | 8.331283 | 0.331283 | 77.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 01032026 | 3 | O3_1110012153.D | 12.500000 | 12.170391 | -0.329609 | 54.000000 | LOW_RISK | 10 | baseline_fallback,low_or_medium_confidence | watch_within_0_5 | low_or_medium_confidence,baseline_fallback,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 06032026 | 4 | O4_1100026044.D | 2.700000 | 3.019879 | 0.319879 | 74.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 13022026 | 9 | O9_1100020744.D | 5.200000 | 5.516856 | 0.316856 | 72.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 16012026 | 2 | O2_1100020391.D | 5.200000 | 4.888017 | -0.311983 | 57.500000 | LOW_RISK | 5 | low_or_medium_confidence | watch_within_0_5 | low_or_medium_confidence,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 14012026 | 7 | O7_1100020307.D | 4.600000 | 4.908195 | 0.308195 | 79.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 06032026 | 15 | O15_1100026030.D | 7.500000 | 7.804240 | 0.304240 | 69.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_clinical_band |
| 14012026 | 12 | O12_1100020191.D | 4.800000 | 5.104052 | 0.304052 | 66.000000 | LOW_RISK | 0 |  | watch_within_0_5 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_clinical_band |
| 16012026 | 5 | O5_1100020377.D | 5.700000 | 5.996175 | 0.296175 | 69.500000 | LOW_RISK | 0 |  | ok_within_0_3 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_inter_operator_band |
| 20032026 | 7 | O7_1100027106.D | 4.100000 | 3.804136 | -0.295864 | 85.000000 | LOW_RISK | 0 |  | ok_within_0_3 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,high_dpa_to_c22_4_ratio,c22_debit_applied,large_rt_error,asymmetric_peak_window,within_inter_operator_band |
| 06032026 | 6 | O6_1100026035.D | 9.300000 | 9.588330 | 0.288330 | 57.000000 | LOW_RISK | 5 | low_or_medium_confidence | ok_within_0_3 | low_or_medium_confidence,c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,within_inter_operator_band |
| 20032026 | 8 | O8_1100027096.D | 7.300000 | 7.012980 | -0.287020 | 69.000000 | LOW_RISK | 0 |  | ok_within_0_3 | c22_complex_boundaries,c20_complex_boundaries,c18_complex_boundaries,c22_credit_applied,large_rt_error,asymmetric_peak_window,within_inter_operator_band |

## Outliers > 0.5

Count: 0

_empty_