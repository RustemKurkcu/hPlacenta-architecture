# output_share

This folder is tracked by git and intended for curated, shareable outputs only.

## Included
- `figures/`: selected `.png` figures copied from `output/figures/`
- `tables/`: selected small summary `.csv` / `.tsv` tables copied from `output/tables/`
- `logs/run_manifest.json`: run metadata (commit hash, date, publish parameters, copied files)

## Excluded
Large intermediates and heavy pipeline outputs should stay in `output/` (git-ignored).
