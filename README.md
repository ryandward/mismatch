# Mismatch Analysis

This repository contains a Python script `mismatch.py` that can be used to calculate the predicted `y_pred` values for a given set of parameters and original vs variant sequences.

## Usage

### To calculate the y_pred values and mismatches spread evenly across a range of values using existing parameters:

The [parameters](https://github.com/traeki/sgRNAs-for-mismatch-CRISPRi/blob/8b76b67f5dd9c69bf9f40917d8cc826890dbf15a/model_param.csv) provided here are calculated and posted by John Hawkins.

```
python mismatch.py mismatches --parameters parameters.csv --min 0 --max 1 --step 0.1 --spacers_file original_guides.tsv
```

This will output a table of original sequences, variant sequences, change descriptions, and predicted `y_pred` values.

| original             | variant              | change_description | y_pred |
|----------------------|----------------------|--------------------|--------|
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAA | C20A               | 0.0636 |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAG | C20G               | 0.2308 |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAT | C20T               | 0.2541 |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACAAGTAC | C15A               | 0.3170 |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTTC | A19T               | 0.3914 |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGGAC | T18G               | 0.5042 |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCTGTAC | A16T               | 0.6004 |
| AAAAACATGTCCACCAGTAC | AAAAAAATGTCCACCAGTAC | C6A                | 0.6957 |
| AAAAACATGTCCACCAGTAC | AAAACCATGTCCACCAGTAC | A5C                | 0.7896 |
| AAAAACATGTCCACCAGTAC | AAAAACACGTCCACCAGTAC | T8C                | 0.8879 |
| AAAAACATGTCCACCAGTAC | CAAAACATGTCCACCAGTAC | A1C                | 0.9918 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAA | C20A               | 0.1365 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAG | C20G               | 0.3037 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAT | C20T               | 0.3270 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCAGCTAC | C15A               | 0.3900 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTCC | A19C               | 0.4097 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGCCCGCTAC | G13C               | 0.5129 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTGC | A19G               | 0.5994 |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTAGCCGCTAC | G12A               | 0.6951 |


### To recalculate `y_pred` values for a set of parameters and existing mismatches:

```
python mismatch.py recalculate --parameters parameters.csv --existing_mismatches mismatches.tsv
```

This will output a table of original sequences, variant sequences, change descriptions, original `y_pred` values, and recalculated `y_pred_new` values.

| original             | variant              | change_description | y_pred | y_pred_new |
|----------------------|----------------------|--------------------|--------|------------|
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAA | C20A               | 0.0636 | 0.0636     |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAG | C20G               | 0.2308 | 0.2308     |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAT | C20T               | 0.2541 | 0.2541     |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACAAGTAC | C15A               | 0.317  | 0.317      |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTTC | A19T               | 0.3914 | 0.3914     |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGGAC | T18G               | 0.5042 | 0.5042     |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCTGTAC | A16T               | 0.6004 | 0.6004     |
| AAAAACATGTCCACCAGTAC | AAAAAAATGTCCACCAGTAC | C6A                | 0.6957 | 0.6957     |
| AAAAACATGTCCACCAGTAC | AAAACCATGTCCACCAGTAC | A5C                | 0.7896 | 0.7896     |
| AAAAACATGTCCACCAGTAC | AAAAACACGTCCACCAGTAC | T8C                | 0.8879 | 0.8879     |
| AAAAACATGTCCACCAGTAC | CAAAACATGTCCACCAGTAC | A1C                | 0.9918 | 0.9918     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAA | C20A               | 0.1365 | 0.1365     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAG | C20G               | 0.3037 | 0.3037     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAT | C20T               | 0.327  | 0.327      |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCAGCTAC | C15A               | 0.39   | 0.39       |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTCC | A19C               | 0.4097 | 0.4097     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGCCCGCTAC | G13C               | 0.5129 | 0.5129     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTGC | A19G               | 0.5994 | 0.5994     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTAGCCGCTAC | G12A               | 0.6951 | 0.6951     |



## Adjustments

Please adjust your input files `parameters.csv` and `original_guides.tsv` (for the `mismatches` command) or `mismatches.tsv` (for the `recalculate` command) accordingly to match your data.
