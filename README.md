Use this repo to calculate mismatches and y_preds aiming for an evenly spread target range, i.e. between 0 and 1, with 0.1 steps as below.

```bash
python make_mismatches.py --parameters parameters.csv --min 0 --max 1 --step 0.1 --spacers_file original_guides.tsv
```

Adjust your files accordingly
