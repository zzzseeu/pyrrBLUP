from __future__ import annotations

from pathlib import Path

import pandas as pd

from pyrrblup.models import AlignmentResult


REQUIRED_GENO_COLUMNS = {"sample"}
REQUIRED_PHENO_COLUMNS = {"sample", "value"}


def align_single_location(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    location_col: str,
) -> AlignmentResult:
    missing_geno = REQUIRED_GENO_COLUMNS.difference(genotype_df.columns)
    missing_pheno = REQUIRED_PHENO_COLUMNS.difference(phenotype_df.columns)
    if missing_geno:
        raise ValueError(f"Missing genotype columns: {sorted(missing_geno)}")
    if missing_pheno:
        raise ValueError(f"Missing phenotype columns: {sorted(missing_pheno)}")
    if location_col not in phenotype_df.columns:
        available = [col for col in phenotype_df.columns if col.startswith("loc_")]
        raise ValueError(f"Unknown location column: {location_col}. Available: {available}")

    pheno_loc = phenotype_df.loc[phenotype_df[location_col] == 1, ["sample", "value"]].copy()
    common = pd.Index(pheno_loc["sample"]).intersection(pd.Index(genotype_df["sample"]), sort=False)
    if common.empty:
        raise ValueError("No overlapping samples remain after location filtering")

    pheno_aligned = pheno_loc.set_index("sample").loc[common, "value"]
    geno_aligned = genotype_df.set_index("sample").loc[common]
    marker_values = geno_aligned.apply(pd.to_numeric)
    if not marker_values.isin([0, 1, 2]).all().all():
        raise ValueError("Genotype matrix contains values outside the supported 0/1/2 encoding")

    return AlignmentResult(samples=common, phenotype=pheno_aligned, markers=marker_values)


def project_root() -> Path:
    here = Path(__file__).resolve()
    if here.parents[3].name == ".worktrees":
        return here.parents[4]
    return here.parents[2]


def resolve_project_data_path(relative_path: str) -> Path:
    return (project_root() / relative_path).resolve()


def load_and_align(
    genotype_path: str,
    phenotype_path: str,
    location_col: str,
) -> AlignmentResult:
    genotype_df = pd.read_csv(genotype_path)
    phenotype_df = pd.read_csv(phenotype_path)
    return align_single_location(genotype_df, phenotype_df, location_col)
