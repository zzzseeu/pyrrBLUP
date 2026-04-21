# Examples

This directory contains deterministic subsets built from the real genotype and phenotype sources used in this repository.

- `minimal_real_data/`
  - `32` samples
  - `64` markers
  - intended for README-scale examples

- `medium_real_data/`
  - `256` samples
  - `512` markers
  - intended for more realistic local runs

Both examples are derived from:

- location: `loc_BeiJ`
- genotype source: `../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv`
- phenotype source: `../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv`

Rebuild them with:

```bash
python scripts/build_examples.py
```
