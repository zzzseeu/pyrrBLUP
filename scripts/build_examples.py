#!/usr/bin/env python
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
import sys

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pyrrblup.io import resolve_project_data_path


GENOTYPE_PATH = resolve_project_data_path(
    "../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv"
)
PHENOTYPE_PATH = resolve_project_data_path(
    "../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv"
)
LOCATION = "loc_BeiJ"
EXAMPLE_SPECS = {
    "minimal_real_data": {"sample_count": 32, "marker_count": 64},
    "medium_real_data": {"sample_count": 256, "marker_count": 512},
}


def _selected_samples(genotype_df: pd.DataFrame, phenotype_df: pd.DataFrame, sample_count: int) -> list[str]:
    pheno_loc = phenotype_df.loc[phenotype_df[LOCATION] == 1, ["sample", "value", LOCATION]].copy()
    common = sorted(set(pheno_loc["sample"]).intersection(genotype_df["sample"]))
    if len(common) < sample_count:
        raise ValueError(f"Requested {sample_count} samples, but only {len(common)} overlap at {LOCATION}")
    return common[:sample_count]


def _selected_markers(genotype_subset: pd.DataFrame, marker_count: int) -> list[str]:
    marker_frame = genotype_subset.drop(columns=["sample"]).apply(pd.to_numeric)
    marker_matrix = marker_frame.to_numpy(dtype=float)
    n_lines = marker_matrix.shape[0]

    frac_missing = np.isnan(marker_matrix).mean(axis=0)
    freq = np.nanmean(marker_matrix + 1.0, axis=0) / 2.0
    maf = np.minimum(freq, 1.0 - freq)
    min_maf = 1.0 / (2.0 * n_lines)
    max_missing = 1.0 - 1.0 / (2.0 * n_lines)
    keep_mask = (maf >= min_maf) & (frac_missing <= max_missing)
    kept_markers = marker_frame.columns[keep_mask].tolist()

    if len(kept_markers) < marker_count:
        raise ValueError(
            f"Requested {marker_count} markers, but only {len(kept_markers)} survive rrBLUP-style filtering"
        )
    return kept_markers[:marker_count]


def _build_example(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    *,
    example_name: str,
    sample_count: int,
    marker_count: int,
) -> None:
    out_dir = ROOT / "examples" / example_name
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = _selected_samples(genotype_df, phenotype_df, sample_count)
    genotype_sample_subset = genotype_df.loc[genotype_df["sample"].isin(samples)].copy()
    genotype_sample_subset = genotype_sample_subset.sort_values("sample").reset_index(drop=True)
    marker_cols = _selected_markers(genotype_sample_subset, marker_count)
    genotype_subset = genotype_sample_subset.loc[:, ["sample", *marker_cols]].copy()

    phenotype_subset = phenotype_df.loc[
        (phenotype_df["sample"].isin(samples)) & (phenotype_df[LOCATION] == 1),
        ["sample", "value", LOCATION],
    ].copy()
    phenotype_subset = phenotype_subset.sort_values("sample").reset_index(drop=True)

    genotype_subset.to_csv(out_dir / "genotype_012.csv", index=False)
    phenotype_subset.to_csv(out_dir / "phenotype.csv", index=False)

    metadata = {
        "source_genotype": str(GENOTYPE_PATH),
        "source_phenotype": str(PHENOTYPE_PATH),
        "location": LOCATION,
        "sample_count": int(sample_count),
        "marker_count": int(marker_count),
        "sample_selection": "samples filtered by location, intersected with genotype, sorted by sample, first N kept",
        "marker_selection": "original genotype marker column order after rrBLUP-style filtering on the selected samples, first M kept",
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main() -> int:
    genotype_df = pd.read_csv(GENOTYPE_PATH)
    phenotype_df = pd.read_csv(PHENOTYPE_PATH)

    for example_name, spec in EXAMPLE_SPECS.items():
        _build_example(
            genotype_df,
            phenotype_df,
            example_name=example_name,
            sample_count=spec["sample_count"],
            marker_count=spec["marker_count"],
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
