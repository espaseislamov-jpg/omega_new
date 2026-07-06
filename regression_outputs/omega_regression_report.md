# Omega regression report

Generated with: `python omega_regression.py --data-dir . --out regression_outputs/omega_regression_current.xlsx`

Total evaluated samples: 286
Overall MAE: 0.250167
Overall RMSE: 0.350054
Overall max abs delta: 1.452719

## Summary

| scope | n | MAE | RMSE | mean_delta | median_abs | std_delta | within_0_2 | within_0_3 | within_0_4 | within_0_5 | within_0_6 | max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ALL | 286 | 0.250167 | 0.350054 | 0.143466 | 0.177972 | 0.319304 | 155 | 207 | 232 | 250 | 260 | 1.452719 |
| OLD45 | 45 | 0.187570 | 0.234535 | 0.091730 | 0.159049 | 0.215852 | 26 | 35 | 39 | 45 | 45 | 0.499520 |
| 01032026 | 4 | 0.429828 | 0.564874 | -0.331323 | 0.363310 | 0.457502 | 2 | 2 | 2 | 2 | 3 | 0.978173 |
| 02072026 | 75 | 0.279618 | 0.383487 | 0.212084 | 0.194682 | 0.319503 | 39 | 52 | 59 | 62 | 64 | 1.271281 |
| 03072026 | 75 | 0.318800 | 0.440543 | 0.203428 | 0.208604 | 0.390762 | 34 | 47 | 52 | 58 | 63 | 1.452719 |
| 06032026 | 15 | 0.199789 | 0.323813 | -0.076018 | 0.116387 | 0.314763 | 10 | 12 | 13 | 14 | 14 | 1.058912 |
| 09032026 | 6 | 0.149625 | 0.166230 | 0.080095 | 0.125339 | 0.145661 | 4 | 6 | 6 | 6 | 6 | 0.255721 |
| 13022026 | 11 | 0.228544 | 0.268469 | 0.154945 | 0.254057 | 0.219243 | 4 | 8 | 9 | 11 | 11 | 0.477590 |
| 13032026 | 13 | 0.131668 | 0.168901 | 0.121041 | 0.114541 | 0.117799 | 11 | 12 | 12 | 13 | 13 | 0.420672 |
| 14012026 | 18 | 0.227811 | 0.255581 | 0.132514 | 0.229453 | 0.218544 | 7 | 13 | 16 | 18 | 18 | 0.499520 |
| 16012026 | 17 | 0.174680 | 0.210700 | 0.066477 | 0.139156 | 0.199938 | 11 | 13 | 17 | 17 | 17 | 0.395057 |
| 18022026 | 9 | 0.167276 | 0.181980 | 0.077303 | 0.132225 | 0.164745 | 5 | 9 | 9 | 9 | 9 | 0.268160 |
| 20022026 | 7 | 0.124124 | 0.160613 | 0.031568 | 0.092139 | 0.157480 | 6 | 6 | 7 | 7 | 7 | 0.362190 |
| 20032026 | 14 | 0.187741 | 0.257550 | 0.012075 | 0.122987 | 0.257267 | 8 | 10 | 11 | 14 | 14 | 0.487020 |
| 23012026 | 7 | 0.149733 | 0.178849 | 0.109435 | 0.141044 | 0.141460 | 6 | 6 | 7 | 7 | 7 | 0.363465 |
| 27022026 | 15 | 0.247052 | 0.338623 | 0.213851 | 0.126101 | 0.262552 | 8 | 11 | 12 | 12 | 14 | 0.859283 |

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
| OK | overestimated_unclassified | 82 |
| OK | underestimated_unclassified | 39 |
| REJECT | underestimated_cluster | 2 |
| REJECT | overestimated_low_confidence | 1 |
| REJECT | underestimated_low_confidence | 1 |
| REVIEW | overestimated_cluster | 94 |
| REVIEW | overestimated_correction_spread | 14 |
| REVIEW | underestimated_cluster | 14 |
| REVIEW | overestimated_low_confidence | 10 |
| REVIEW | underestimated_unclassified | 10 |
| REVIEW | underestimated_low_confidence | 8 |
| REVIEW | overestimated_unclassified | 7 |
| REVIEW | underestimated_correction_spread | 4 |

## Outliers > 0.5

Count: 36

| date | sample_no | instrument_no | sample_id | sample_name | match_method | reference | calculated | delta | confidence | selected_variant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 03072026 | 105069916298 | 27 | 925663916002 | O27_925663916002.D | position_date_override | 2.500000 | 3.952719 | 1.452719 | 68.000000 | current_pipeline |
| 02072026 | 105067678302 | 4 | 105067678302 | O4_105067678302.D | sample_id | 4.930000 | 3.658719 | -1.271281 | 70.000000 | current_pipeline |
| 03072026 | 900512824605 | 8 | 104837397699 | O8_104837397699.D | position_date_override | 3.000000 | 4.251493 | 1.251493 | 84.000000 | current_pipeline |
| 03072026 | 1110012947 | 69 | 1110012953 | O69_1110012953.D | position_date_override | 2.700000 | 3.928291 | 1.228291 | 78.000000 | current_pipeline |
| 06032026 | 2 | 2 | 1110012176 | O2_1110012176.D | instrument_no | 5.200000 | 4.141088 | -1.058912 | 80.000000 | current_pipeline |
| 02072026 | 929850052502 | 66 | 929850052502 | O66_929850052502.D | sample_id | 5.440000 | 6.429513 | 0.989513 | 82.000000 | current_pipeline |
| 03072026 | 105068695803 | 10 | 929990895102 | O10_929990895102.D | position_date_override | 3.020000 | 3.999577 | 0.979577 | 60.000000 | current_pipeline |
| 01032026 | 4 | 4 | 5555839154 | O4_5555839154.D | instrument_no | 5.300000 | 4.321827 | -0.978173 | 60.000000 | current_pipeline |
| 02072026 | 903862928199 | 71 | 903862928199 | O71_903862928199.D | sample_id | 3.950000 | 4.845298 | 0.895298 | 54.000000 | current_pipeline |
| 03072026 | 1110012934 | 66 | 1110012914 | O66_1110012914.D | position_date_override | 2.700000 | 3.587572 | 0.887572 | 82.000000 | current_pipeline |
| 27022026 | 14 | 14 | 1100021859 | O14_1100021859.D | instrument_no | 8.500000 | 9.359283 | 0.859283 | 64.000000 | current_pipeline |
| 03072026 | 1110012917 | 72 | 1110012951 | O72_1110012951.D | position_date_override | 3.400000 | 4.242102 | 0.842102 | 60.000000 | current_pipeline |
| 02072026 | 104714605399 | 30 | 104714605399 | O30_104714605399.D | sample_id | 2.460000 | 3.235693 | 0.775693 | 68.000000 | current_pipeline |
| 02072026 | 105066911701 | 19 | 105066911701 | O19_105066911701.D | sample_id | 2.810000 | 3.568394 | 0.758394 | 72.000000 | current_pipeline |
| 02072026 | 922730889202 | 15 | 922730889202 | O15_922730889202.D | sample_id | 4.240000 | 4.996285 | 0.756285 | 77.000000 | current_pipeline |
| 02072026 | 105066909103 | 68 | 105066909103 | O68_105066909103.D | sample_id | 3.420000 | 4.157840 | 0.737840 | 60.000000 | current_pipeline |
| 02072026 | 1110012961 | 67 | 1110012961 | O67_1110012961.D | sample_id | 3.200000 | 3.918485 | 0.718485 | 72.000000 | current_pipeline |
| 02072026 | 105069117602 | 16 | 105069117602 | O16_105069117602.D | sample_id | 3.220000 | 3.913542 | 0.693542 | 62.000000 | current_pipeline |
| 03072026 | 105063823498 | 35 | 105068114497 | O35_105068114497.D | position_date_override | 2.920000 | 3.585613 | 0.665613 | 72.000000 | current_pipeline |
| 03072026 | 104962405001 | 40 | 105070893202 | O40_105070893202.D | position_date_override | 9.830000 | 9.180991 | -0.649009 | 60.000000 | current_pipeline |
| 02072026 | 105015060401 | 31 | 105015060401 | O31_105015060401.D | sample_id | 2.990000 | 3.636754 | 0.646754 | 72.000000 | current_pipeline |
| 02072026 | 105067123002 | 29 | 105067123002 | O29_105067123002.D | sample_id | 3.270000 | 3.911328 | 0.641328 | 72.000000 | current_pipeline |
| 03072026 | 905564620603 | 16 | 905564619906 | O16_905564619906.D | position_date_override | 4.000000 | 4.626670 | 0.626670 | 70.000000 | current_pipeline |
| 03072026 | 1110012954 | 55 | 1110012949 | O55_1110012949.D | position_date_override | 3.000000 | 3.620775 | 0.620775 | 72.000000 | current_pipeline |
| 03072026 | 105068114497 | 36 | 105071495901 | O36_105071495901.D | position_date_override | 5.130000 | 5.748970 | 0.618970 | 54.000000 | current_pipeline |
| 03072026 | 105067908799 | 18 | 105067955299 | O18_105067955299.D | position_date_override | 3.120000 | 3.720618 | 0.600618 | 72.000000 | current_pipeline |
| 03072026 | 905564619506 | 53 | 1110012933 | O53_1110012933.D | position_date_override | 2.960000 | 3.551859 | 0.591859 | 72.000000 | current_pipeline |
| 02072026 | 905564618106 | 75 | 905564618106 | O75_905564618106.D | sample_id | 4.260000 | 4.850326 | 0.590326 | 54.000000 | current_pipeline |
| 03072026 | 1110012960 | 62 | 1110012946 | O62_1110012946.D | position_date_override | 2.800000 | 3.375801 | 0.575801 | 82.000000 | current_pipeline |
| 03072026 | 1110012916 | 57 | 1110012957 | O57_1110012957.D | position_date_override | 3.000000 | 3.573486 | 0.573486 | 62.000000 | current_pipeline |
| 27022026 | 5 | 5 | 1100021851 | O5_1100021851.D | instrument_no | 5.200000 | 5.750894 | 0.550894 | 45.000000 | current_pipeline |
| 01032026 | 3 | 3 | 1110012153 | O3_1110012153.D | instrument_no | 12.500000 | 11.970391 | -0.529609 | 54.000000 | current_pipeline |
| 27022026 | 9 | 9 | 1100021818 | O9_1100021818.D | instrument_no | 5.100000 | 5.612960 | 0.512960 | 70.000000 | current_pipeline |
| 02072026 | 905564617506 | 63 | 905564617506 | O63_905564617506.D | sample_id | 5.000000 | 5.504731 | 0.504731 | 62.000000 | current_pipeline |
| 03072026 | 905564620206 | 12 | 105071372802 | O12_105071372802.D | position_date_override | 2.720000 | 3.224081 | 0.504081 | 77.000000 | current_pipeline |
| 03072026 | 1110012918 | 65 | 1110012934 | O65_1110012934.D | position_date_override | 3.300000 | 3.801776 | 0.501776 | 82.000000 | current_pipeline |