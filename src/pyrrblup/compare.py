from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from pyrrblup.models import BatchComparisonResult


def _load_and_align_beta(r_dir: Path, python_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    r_beta = pd.read_csv(r_dir / "mixed_solve_beta.csv")
    python_beta = pd.read_csv(python_dir / "mixed_solve_beta.csv")

    if "term" in r_beta.columns and "term" in python_beta.columns:
        python_beta = python_beta.set_index("term").loc[r_beta["term"]].reset_index()
    else:
        r_beta = r_beta.reset_index(names="row_id")
        python_beta = python_beta.reset_index(names="row_id")
        python_beta = python_beta.set_index("row_id").loc[r_beta["row_id"]].reset_index()

    return r_beta, python_beta


def compare_location_outputs(
    location: str,
    r_dir: Path,
    python_dir: Path,
    r_success: bool,
    python_success: bool,
) -> BatchComparisonResult:
    if not (r_success and python_success):
        return BatchComparisonResult(
            location=location,
            r_success=r_success,
            python_success=python_success,
            sample_count=0,
            amat_rmse=None,
            prediction_rmse=None,
            beta_rmse=None,
            u_rmse=None,
            vu_r=None,
            vu_py=None,
            ve_r=None,
            ve_py=None,
            error_stage="comparison",
            error_message="comparison skipped because one side failed",
        )

    r_amat = pd.read_csv(r_dir / "amat.csv").set_index("sample")
    python_amat = pd.read_csv(python_dir / "amat.csv")
    if python_amat.columns[0] != "sample":
        python_amat = python_amat.rename(columns={python_amat.columns[0]: "sample"})
    python_amat = python_amat.set_index("sample").loc[r_amat.index, r_amat.columns]
    amat_diff = r_amat.to_numpy(float) - python_amat.to_numpy(float)

    r_pred = pd.read_csv(r_dir / "kin_blup_predictions.csv")
    python_pred = pd.read_csv(python_dir / "kin_blup_predictions.csv")
    pred_merged = r_pred.merge(python_pred, on="sample")
    pred_diff = pred_merged["pred"] - pred_merged[python_pred.columns[-1]]

    r_beta, python_beta = _load_and_align_beta(r_dir, python_dir)
    beta_diff = r_beta["beta"].to_numpy(float) - python_beta["beta"].to_numpy(float)

    r_u = pd.read_csv(r_dir / "mixed_solve_u.csv")
    python_u = pd.read_csv(python_dir / "mixed_solve_u.csv")
    u_merged = r_u.merge(python_u, on="sample", suffixes=("_r", "_py"))
    u_diff = u_merged["u_r"] - u_merged["u_py"]

    r_var = pd.read_csv(r_dir / "mixed_solve_varcomp.csv").iloc[0]
    python_var = json.loads((python_dir / "mixed_solve_varcomp.json").read_text(encoding="utf-8"))

    return BatchComparisonResult(
        location=location,
        r_success=True,
        python_success=True,
        sample_count=int(len(pred_merged)),
        amat_rmse=float(np.sqrt(np.mean(amat_diff**2))),
        prediction_rmse=float(np.sqrt(np.mean(pred_diff**2))),
        beta_rmse=float(np.sqrt(np.mean(beta_diff**2))),
        u_rmse=float(np.sqrt(np.mean(u_diff**2))),
        vu_r=float(r_var["vu"]),
        vu_py=float(python_var["vu"]),
        ve_r=float(r_var["ve"]),
        ve_py=float(python_var["ve"]),
        error_stage="",
        error_message="",
    )
