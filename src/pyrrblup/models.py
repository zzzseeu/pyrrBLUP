from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(slots=True)
class AlignmentResult:
    samples: pd.Index
    phenotype: pd.Series
    markers: pd.DataFrame


@dataclass(slots=True)
class AMatResult:
    matrix: pd.DataFrame


@dataclass(slots=True)
class MixedSolveResult:
    beta: np.ndarray
    u: np.ndarray
    Vu: float
    Ve: float
    beta_terms: list[str]
    u_terms: list[str]


@dataclass(slots=True)
class KinBlupResult:
    predictions: pd.Series
    mixed_solve: MixedSolveResult
    kinship: pd.DataFrame
