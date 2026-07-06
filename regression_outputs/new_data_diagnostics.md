# New archive data diagnostics

Generated after extracting `Desktop.part1.rar`/`Desktop.part2.rar` with `unar` into `extracted_desktop/` and running the current regression.

## Input coverage

| date | reference_rows | instrument_batches | matched_by_sample_id | missing_reference_ids |
| --- | --- | --- | --- | --- |
| 02072026 | 75 | 75 | 75 | 0 |
| 03072026 | 76 | 75 | 75 | 1 |

## Current-engine summary

| date | n | MAE | RMSE | median_abs | max_abs | within_0_5 |
| --- | --- | --- | --- | --- | --- | --- |
| 02072026 | 75 | 0.175944 | 0.248407 | 0.126947 | 1.271281 | 73 |
| 03072026 | 75 | 0.201826 | 0.276635 | 0.168388 | 1.002719 | 70 |

## Mapping notes

- `02072026` is fully extracted and matched by sample ID: 75 manual rows, 75 instrument batches, 75 sample-id matches.
- `03072026` is intentionally matched by row position rather than sample ID because this workbook is shifted: the first manual row corresponds to the first instrument batch `O2`, and the final manual row has no matching instrument batch.

## Algorithm note

The data-driven C20/EPA overlap-credit model remains disabled by default. This run also applies a bounded C22/DPA over-integration debit for high DPA/C22:4-ratio clusters. On the extracted full corpus, `ALL` MAE is `0.1927` and max absolute error is `1.4367`.

## Review / outlier classification for July data

| review_flag | outlier_class | n |
| --- | --- | --- |
| OK | overestimated_unclassified | 22 |
| OK | underestimated_unclassified | 26 |
| REJECT | overestimated_cluster | 2 |
| REJECT | overestimated_low_confidence | 1 |
| REJECT | underestimated_cluster | 2 |
| REVIEW | overestimated_cluster | 78 |
| REVIEW | overestimated_correction_spread | 1 |
| REVIEW | overestimated_low_confidence | 4 |
| REVIEW | underestimated_cluster | 12 |
| REVIEW | underestimated_low_confidence | 1 |
| REVIEW | underestimated_unclassified | 1 |

## Errors / skipped rows

| date | sample_no | reference | match_method | error |
| --- | --- | --- | --- | --- |
| 03072026 | 1110012956 | 3.800000 | missing_position_date_override | No matching instrument batch |

## Top 20 new-data outliers

| date | sample_name | reference | calculated | delta | confidence | review_flag | outlier_class |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 02072026 | O4_105067678302.D | 4.930000 | 3.658719 | -1.271281 | 70.000000 | OK | underestimated_unclassified |
| 03072026 | O27_925663916002.D | 2.500000 | 3.502719 | 1.002719 | 50.000000 | REVIEW | overestimated_cluster |
| 03072026 | O8_104837397699.D | 3.000000 | 3.801493 | 0.801493 | 74.000000 | OK | overestimated_unclassified |
| 03072026 | O69_1110012953.D | 2.700000 | 3.478291 | 0.778291 | 60.000000 | REVIEW | overestimated_cluster |
| 03072026 | O40_105070893202.D | 9.830000 | 9.180991 | -0.649009 | 60.000000 | OK | underestimated_unclassified |
| 02072026 | O71_903862928199.D | 3.950000 | 4.577730 | 0.627730 | 62.000000 | REVIEW | overestimated_cluster |
| 03072026 | O10_929990895102.D | 3.020000 | 3.549577 | 0.529577 | 42.000000 | REVIEW | overestimated_cluster |
