# Omega regression data

Each labeled batch uses a pair of files:

- `DDMMYYYY.CSV` contains the raw chromatograms exported by the instrument.
- `test_bigbatch_DDMMYYYY.xlsx` contains the corresponding manual omega values.

The first reference row in `test_bigbatch_03072026.xlsx` belongs to the blank
`O1` injection, which is not present in `03072026.CSV`. Exclude that row when
building the labeled sample table; the remaining 75 sample IDs match the 75
chromatograms exactly.

The `06072026.CSV` batch currently has no manual reference workbook.

The `14072026.CSV` batch is the sealed final evaluation set. Its manual results
must not be added to the repository or used during feature design, training, or
model selection. They may be opened only for the final one-time evaluation of a
frozen quality model.
