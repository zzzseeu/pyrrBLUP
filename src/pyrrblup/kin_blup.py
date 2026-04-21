from __future__ import annotations

from pathlib import Path

import pandas as pd

from pyrrblup.amat import amat
from pyrrblup.io import align_single_location, load_and_align
from pyrrblup.mixed_solve import mixed_solve
from pyrrblup.models import KinBlupResult


def _fit_from_aligned(aligned) -> KinBlupResult:
    kinship = amat(aligned.markers).matrix
    fit = mixed_solve(aligned.phenotype.to_numpy(), K=kinship.to_numpy())
    predictions = pd.Series(fit.beta[0] + fit.u, index=aligned.samples, name="prediction")
    return KinBlupResult(predictions=predictions, mixed_solve=fit, kinship=kinship)


def kin_blup(
    genotype: str | Path | pd.DataFrame,
    phenotype: str | Path | pd.DataFrame,
    location_col: str,
) -> KinBlupResult:
    genotype_is_frame = isinstance(genotype, pd.DataFrame)
    phenotype_is_frame = isinstance(phenotype, pd.DataFrame)
    genotype_is_path = isinstance(genotype, (str, Path))
    phenotype_is_path = isinstance(phenotype, (str, Path))

    if genotype_is_frame and phenotype_is_frame:
        return kin_blup_from_frames(genotype, phenotype, location_col)
    if genotype_is_path and phenotype_is_path:
        return kin_blup_from_paths(genotype, phenotype, location_col)

    raise TypeError("genotype and phenotype must both be paths or both be pandas.DataFrame")


def kin_blup_from_frames(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    location_col: str,
) -> KinBlupResult:
    aligned = align_single_location(genotype_df, phenotype_df, location_col)
    return _fit_from_aligned(aligned)


def kin_blup_from_paths(genotype_path: str | Path, phenotype_path: str | Path, location_col: str) -> KinBlupResult:
    aligned = load_and_align(genotype_path, phenotype_path, location_col)
    return _fit_from_aligned(aligned)
