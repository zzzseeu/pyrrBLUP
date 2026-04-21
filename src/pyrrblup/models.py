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


@dataclass(slots=True)
class BatchRunStatus:
    location: str
    tool: str
    success: bool
    error_message: str
    started_at: str
    finished_at: str


@dataclass(slots=True)
class BatchComparisonResult:
    location: str
    r_success: bool
    python_success: bool
    sample_count: int
    amat_rmse: float | None
    prediction_rmse: float | None
    beta_rmse: float | None
    u_rmse: float | None
    vu_r: float | None
    vu_py: float | None
    ve_r: float | None
    ve_py: float | None
    error_stage: str
    error_message: str
