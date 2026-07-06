# Omega regression report

Generated with: `python omega_regression.py --data-dir . --out regression_outputs/omega_regression_current.xlsx --debug-dir regression_debug --debug-threshold 0.5`

Total evaluated samples: 286
Overall MAE: 0.312341
Overall RMSE: 0.485320
Overall max abs delta: 2.479391

## Summary

| scope | n | MAE | RMSE | mean_delta | median_abs | std_delta | within_0_2 | within_0_3 | within_0_4 | within_0_5 | within_0_6 | max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ALL | 286 | 0.312341 | 0.485320 | 0.229466 | 0.197855 | 0.427645 | 148 | 195 | 222 | 238 | 248 | 2.479391 |
| OLD45 | 45 | 0.195044 | 0.241717 | 0.099509 | 0.171403 | 0.220284 | 26 | 34 | 40 | 43 | 45 | 0.528791 |
| 01032026 | 4 | 0.235608 | 0.300012 | -0.137103 | 0.199152 | 0.266852 | 2 | 3 | 3 | 3 | 4 | 0.529609 |
| 02072026 | 75 | 0.399631 | 0.591811 | 0.355696 | 0.217957 | 0.472992 | 37 | 46 | 51 | 56 | 57 | 1.679417 |
| 03072026 | 75 | 0.440805 | 0.663716 | 0.339217 | 0.254169 | 0.570483 | 29 | 43 | 49 | 52 | 56 | 2.479391 |
| 06032026 | 15 | 0.171176 | 0.206725 | 0.064329 | 0.126137 | 0.196461 | 10 | 11 | 14 | 15 | 15 | 0.404240 |
| 09032026 | 6 | 0.161212 | 0.174463 | 0.091683 | 0.136659 | 0.148431 | 4 | 6 | 6 | 6 | 6 | 0.255721 |
| 13022026 | 11 | 0.248205 | 0.293994 | 0.185844 | 0.254057 | 0.227804 | 4 | 7 | 8 | 11 | 11 | 0.477590 |
| 13032026 | 13 | 0.142829 | 0.175848 | 0.132202 | 0.134655 | 0.115954 | 11 | 12 | 12 | 13 | 13 | 0.420672 |
| 14012026 | 18 | 0.224947 | 0.252366 | 0.129650 | 0.220086 | 0.216516 | 8 | 13 | 17 | 17 | 18 | 0.528791 |
| 16012026 | 17 | 0.176841 | 0.213084 | 0.052169 | 0.135054 | 0.206599 | 11 | 13 | 17 | 17 | 17 | 0.395057 |
| 18022026 | 9 | 0.168313 | 0.183385 | 0.078340 | 0.132225 | 0.165810 | 5 | 9 | 9 | 9 | 9 | 0.268160 |
| 20022026 | 7 | 0.127424 | 0.162721 | 0.036445 | 0.097620 | 0.158587 | 6 | 6 | 7 | 7 | 7 | 0.362190 |
| 20032026 | 14 | 0.205082 | 0.277852 | 0.030400 | 0.157968 | 0.276184 | 7 | 9 | 11 | 13 | 14 | 0.500939 |
| 23012026 | 7 | 0.159233 | 0.186246 | 0.118935 | 0.141044 | 0.143325 | 6 | 6 | 7 | 7 | 7 | 0.363465 |
| 27022026 | 15 | 0.252300 | 0.345429 | 0.222701 | 0.138050 | 0.264055 | 8 | 11 | 11 | 12 | 14 | 0.859283 |

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
| OK | overestimated_unclassified | 84 |
| OK | underestimated_unclassified | 31 |
| REJECT | overestimated_cluster | 15 |
| REJECT | underestimated_cluster | 4 |
| REJECT | underestimated_low_confidence | 2 |
| REJECT | overestimated_low_confidence | 1 |
| REVIEW | overestimated_cluster | 82 |
| REVIEW | overestimated_correction_spread | 14 |
| REVIEW | overestimated_low_confidence | 12 |
| REVIEW | underestimated_unclassified | 10 |
| REVIEW | overestimated_unclassified | 9 |
| REVIEW | underestimated_cluster | 9 |
| REVIEW | underestimated_low_confidence | 9 |
| REVIEW | underestimated_correction_spread | 4 |

## Outliers > 0.5

Count: 48

| date | sample_no | instrument_no | sample_id | sample_name | match_method | reference | calculated | delta | confidence | selected_variant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 03072026 | 105069916298 | 27 | 925663916002 | O27_925663916002.D | position_date_override | 2.500000 | 4.979391 | 2.479391 | 45.000000 | current_pipeline |
| 03072026 | 1110012947 | 69 | 1110012953 | O69_1110012953.D | position_date_override | 2.700000 | 4.510912 | 1.810912 | 55.000000 | current_pipeline |
| 03072026 | 105068695803 | 10 | 929990895102 | O10_929990895102.D | position_date_override | 3.020000 | 4.712192 | 1.692192 | 37.000000 | current_pipeline |
| 02072026 | 105068236802 | 65 | 105068236802 | O65_105068236802.D | sample_id | 3.490000 | 5.169417 | 1.679417 | 31.000000 | current_pipeline |
| 03072026 | 1110012916 | 57 | 1110012957 | O57_1110012957.D | position_date_override | 3.000000 | 4.618952 | 1.618952 | 39.000000 | current_pipeline |
| 02072026 | 104714605399 | 30 | 104714605399 | O30_104714605399.D | sample_id | 2.460000 | 4.077767 | 1.617767 | 45.000000 | current_pipeline |
| 02072026 | 105069117602 | 16 | 105069117602 | O16_105069117602.D | sample_id | 3.220000 | 4.824124 | 1.604124 | 39.000000 | current_pipeline |
| 02072026 | 105066909103 | 68 | 105066909103 | O68_105066909103.D | sample_id | 3.420000 | 4.970112 | 1.550112 | 43.000000 | current_pipeline |
| 02072026 | 903862928199 | 71 | 903862928199 | O71_903862928199.D | sample_id | 3.950000 | 5.490953 | 1.540953 | 31.000000 | current_pipeline |
| 03072026 | 1110012914 | 67 | 1110012959 | O67_1110012959.D | position_date_override | 1.800000 | 3.237897 | 1.437897 | 39.000000 | current_pipeline |
| 03072026 | 1110012917 | 72 | 1110012951 | O72_1110012951.D | position_date_override | 3.400000 | 4.806425 | 1.406425 | 37.000000 | current_pipeline |
| 03072026 | 905564620603 | 16 | 905564619906 | O16_905564619906.D | position_date_override | 4.000000 | 5.269973 | 1.269973 | 47.000000 | current_pipeline |
| 03072026 | 900512824605 | 8 | 104837397699 | O8_104837397699.D | position_date_override | 3.000000 | 4.253579 | 1.253579 | 84.000000 | current_pipeline |
| 02072026 | 905564618106 | 75 | 905564618106 | O75_905564618106.D | sample_id | 4.260000 | 5.481613 | 1.221613 | 31.000000 | current_pipeline |
| 03072026 | 905564620803 | 45 | 104236538002 | O45_104236538002.D | position_date_override | 3.010000 | 4.214649 | 1.204649 | 39.000000 | current_pipeline |
| 02072026 | 105067123002 | 29 | 105067123002 | O29_105067123002.D | sample_id | 3.270000 | 4.379894 | 1.109894 | 49.000000 | current_pipeline |
| 02072026 | 905564617506 | 63 | 905564617506 | O63_905564617506.D | sample_id | 5.000000 | 6.090876 | 1.090876 | 39.000000 | current_pipeline |
| 02072026 | 929850052502 | 66 | 929850052502 | O66_929850052502.D | sample_id | 5.440000 | 6.455657 | 1.015657 | 82.000000 | current_pipeline |
| 03072026 | 105068114497 | 36 | 105071495901 | O36_105071495901.D | position_date_override | 5.130000 | 6.100919 | 0.970919 | 39.000000 | current_pipeline |
| 03072026 | 1110012934 | 66 | 1110012914 | O66_1110012914.D | position_date_override | 2.700000 | 3.667673 | 0.967673 | 77.000000 | current_pipeline |
| 03072026 | 1110012957 | 58 | 1110012945 | O58_1110012945.D | position_date_override | 3.200000 | 4.143902 | 0.943902 | 45.000000 | current_pipeline |
| 02072026 | 104895240699 | 49 | 104895240699 | O49_104895240699.D | sample_id | 3.530000 | 4.448388 | 0.918388 | 49.000000 | current_pipeline |
| 02072026 | 105068553703 | 33 | 105068553703 | O33_105068553703.D | sample_id | 3.690000 | 4.566451 | 0.876451 | 41.000000 | current_pipeline |
| 27022026 | 14 | 14 | 1100021859 | O14_1100021859.D | instrument_no | 8.500000 | 9.359283 | 0.859283 | 64.000000 | current_pipeline |
| 02072026 | 105065793597 | 1 | 105065793597 | O1_105065793597.D | sample_id | 3.300000 | 4.143971 | 0.843971 | 69.000000 | current_pipeline |
| 03072026 | 1110012949 | 56 | 1110012916 | O56_1110012916.D | position_date_override | 3.400000 | 4.236324 | 0.836324 | 45.000000 | current_pipeline |
| 02072026 | 922730889202 | 15 | 922730889202 | O15_922730889202.D | sample_id | 4.240000 | 5.075555 | 0.835555 | 77.000000 | current_pipeline |
| 02072026 | 105066911701 | 19 | 105066911701 | O19_105066911701.D | sample_id | 2.810000 | 3.573096 | 0.763096 | 72.000000 | current_pipeline |
| 02072026 | 105067678302 | 4 | 105067678302 | O4_105067678302.D | sample_id | 4.930000 | 4.167045 | -0.762955 | 47.000000 | current_pipeline |
| 02072026 | 1110012961 | 67 | 1110012961 | O67_1110012961.D | sample_id | 3.200000 | 3.918620 | 0.718620 | 72.000000 | current_pipeline |
| 02072026 | 905564617703 | 62 | 905564617703 | O62_905564617703.D | sample_id | 4.480000 | 5.174218 | 0.694218 | 77.000000 | current_pipeline |
| 03072026 | 105063823498 | 35 | 105068114497 | O35_105068114497.D | position_date_override | 2.920000 | 3.586143 | 0.666143 | 72.000000 | current_pipeline |
| 02072026 | 105015060401 | 31 | 105015060401 | O31_105015060401.D | sample_id | 2.990000 | 3.640050 | 0.650050 | 72.000000 | current_pipeline |
| 03072026 | 104962405001 | 40 | 105070893202 | O40_105070893202.D | position_date_override | 9.830000 | 9.180991 | -0.649009 | 60.000000 | current_pipeline |
| 03072026 | 905564620206 | 12 | 105071372802 | O12_105071372802.D | position_date_override | 2.720000 | 3.356295 | 0.636295 | 67.000000 | current_pipeline |
| 03072026 | 1110012954 | 55 | 1110012949 | O55_1110012949.D | position_date_override | 3.000000 | 3.621859 | 0.621859 | 72.000000 | current_pipeline |
| 03072026 | 105067908799 | 18 | 105067955299 | O18_105067955299.D | position_date_override | 3.120000 | 3.724187 | 0.604187 | 72.000000 | current_pipeline |
| 03072026 | 1110012958 | 60 | 1110012915 | O60_1110012915.D | position_date_override | 4.100000 | 4.702476 | 0.602476 | 41.000000 | current_pipeline |
| 03072026 | 905564619506 | 53 | 1110012933 | O53_1110012933.D | position_date_override | 2.960000 | 3.556348 | 0.596348 | 72.000000 | current_pipeline |
| 03072026 | 1110012960 | 62 | 1110012946 | O62_1110012946.D | position_date_override | 2.800000 | 3.381476 | 0.581476 | 82.000000 | current_pipeline |
| 27022026 | 5 | 5 | 1100021851 | O5_1100021851.D | instrument_no | 5.200000 | 5.750894 | 0.550894 | 45.000000 | current_pipeline |
| 27022026 | 9 | 9 | 1100021818 | O9_1100021818.D | instrument_no | 5.100000 | 5.647341 | 0.547341 | 70.000000 | current_pipeline |
| 01032026 | 3 | 3 | 1110012153 | O3_1110012153.D | instrument_no | 12.500000 | 11.970391 | -0.529609 | 54.000000 | current_pipeline |
| 14012026 | 6 | 6 | 1100020225 | O6_1100020225.D | instrument_no | 3.800000 | 4.328791 | 0.528791 | 90.000000 | current_pipeline |
| 02072026 | 905564618203 | 64 | 905564618203 | O64_905564618203.D | sample_id | 4.850000 | 5.363438 | 0.513438 | 39.000000 | current_pipeline |
| 03072026 | 1110012918 | 65 | 1110012934 | O65_1110012934.D | position_date_override | 3.300000 | 3.801917 | 0.501917 | 82.000000 | current_pipeline |
| 20032026 | 13 | 13 | 1100027101 | O13_1100027101.D | instrument_no | 3.300000 | 3.800939 | 0.500939 | 77.000000 | current_pipeline |
| 03072026 | 905564620406 | 44 | 905564620803 | O44_905564620803.D | position_date_override | 4.930000 | 5.430698 | 0.500698 | 74.000000 | current_pipeline |