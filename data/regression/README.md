# Omega regression data

Each labeled batch uses a pair of files:

- `DDMMYYYY.CSV` contains the raw chromatograms exported by the instrument.
- `test_bigbatch_DDMMYYYY.xlsx` contains the corresponding manual omega values.

`03072026` is a special acquisition. The workbook has 76 reference rows, while
the CSV has 75 non-empty chromatograms (`O2` through `O76`) after an empty first
injection. The printed IDs are shifted relative to the actual signal sequence:
the first workbook result belongs to the first non-empty chromatogram, the
second result to the second chromatogram, and so on. Match this date by
acquisition position, not by the printed sample ID; the final workbook row has
no non-empty chromatogram. Direct ID matching produces an artificial alternating
error pattern (64 false `>0.5` labels).

The `14072026.CSV` batch is the sealed final evaluation set. Its manual results
must not be added to the repository or used during feature design, training, or
model selection. They may be opened only for the final one-time evaluation of a
frozen quality model.
