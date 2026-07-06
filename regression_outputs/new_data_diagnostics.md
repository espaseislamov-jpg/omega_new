# New archive data diagnostics

Generated after extracting `Desktop.part1.rar`/`Desktop.part2.rar` with `unar` into `extracted_desktop/` and running current regression.

## Input coverage

| date | reference_rows | instrument_batches | matched_by_sample_id | missing_reference_ids |
| --- | --- | --- | --- | --- |
| 02072026 | 75 | 75 | 75 | 0 |
| 03072026 | 76 | 75 | 75 | 1 |

## Current-engine summary

| date | n | MAE | RMSE | median_abs | max_abs | within_0_5 |
| --- | --- | --- | --- | --- | --- | --- |
| 02072026 | 75 | 0.399631 | 0.591811 | 0.217957 | 1.679417 | 56 |
| 03072026 | 75 | 0.440805 | 0.663716 | 0.254169 | 2.479391 | 52 |

## Mapping notes

- `02072026` is fully extracted and matched by sample ID: 75 manual rows, 75 instrument batches, 75 sample-id matches.
- `03072026` is intentionally matched by row position rather than sample ID because this workbook is shifted: the first manual row corresponds to the first instrument batch `O2`, and the final manual row has no matching instrument batch. This lowers the 03072026 MAE from 1.836145 to 0.440805.

## Errors / skipped rows

| date | sample_no | reference | match_method | error |
| --- | --- | --- | --- | --- |
| 03072026 | 1110012956 | 3.8 | missing_position_date_override | No matching instrument batch |

## Top 20 new-data outliers

| date | sample_name | reference | calculated | delta | confidence |
| --- | --- | --- | --- | --- | --- |
| 03072026 | O27_925663916002.D | 2.500000 | 4.979391 | 2.479391 | 45.00 |
| 03072026 | O69_1110012953.D | 2.700000 | 4.510912 | 1.810912 | 55.00 |
| 03072026 | O10_929990895102.D | 3.020000 | 4.712192 | 1.692192 | 37.00 |
| 02072026 | O65_105068236802.D | 3.490000 | 5.169417 | 1.679417 | 31.00 |
| 03072026 | O57_1110012957.D | 3.000000 | 4.618952 | 1.618952 | 39.00 |
| 02072026 | O30_104714605399.D | 2.460000 | 4.077767 | 1.617767 | 45.00 |
| 02072026 | O16_105069117602.D | 3.220000 | 4.824124 | 1.604124 | 39.00 |
| 02072026 | O68_105066909103.D | 3.420000 | 4.970112 | 1.550112 | 43.00 |
| 02072026 | O71_903862928199.D | 3.950000 | 5.490953 | 1.540953 | 31.00 |
| 03072026 | O67_1110012959.D | 1.800000 | 3.237897 | 1.437897 | 39.00 |
| 03072026 | O72_1110012951.D | 3.400000 | 4.806425 | 1.406425 | 37.00 |
| 03072026 | O16_905564619906.D | 4.000000 | 5.269973 | 1.269973 | 47.00 |
| 03072026 | O8_104837397699.D | 3.000000 | 4.253579 | 1.253579 | 84.00 |
| 02072026 | O75_905564618106.D | 4.260000 | 5.481613 | 1.221613 | 31.00 |
| 03072026 | O45_104236538002.D | 3.010000 | 4.214649 | 1.204649 | 39.00 |
| 02072026 | O29_105067123002.D | 3.270000 | 4.379894 | 1.109894 | 49.00 |
| 02072026 | O63_905564617506.D | 5.000000 | 6.090876 | 1.090876 | 39.00 |
| 02072026 | O66_929850052502.D | 5.440000 | 6.455657 | 1.015657 | 82.00 |
| 03072026 | O36_105071495901.D | 5.130000 | 6.100919 | 0.970919 | 39.00 |
| 03072026 | O66_1110012914.D | 2.700000 | 3.667673 | 0.967673 | 77.00 |