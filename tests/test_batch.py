from pathlib import Path

import pandas as pd

from pyrrblup.batch import run_batch_locations, summarize_batch_results
from pyrrblup.models import BatchComparisonResult


def test_summarize_batch_results_writes_summary_for_mixed_success(tmp_path: Path) -> None:
    results = [
        BatchComparisonResult(
            location="loc_A",
            r_success=True,
            python_success=True,
            sample_count=10,
            amat_rmse=0.0,
            prediction_rmse=0.0,
            beta_rmse=0.0,
            u_rmse=0.0,
            vu_r=1.0,
            vu_py=1.0,
            ve_r=2.0,
            ve_py=2.0,
            error_stage="",
            error_message="",
        ),
        BatchComparisonResult(
            location="loc_B",
            r_success=False,
            python_success=True,
            sample_count=0,
            amat_rmse=None,
            prediction_rmse=None,
            beta_rmse=None,
            u_rmse=None,
            vu_r=None,
            vu_py=None,
            ve_r=None,
            ve_py=None,
            error_stage="r",
            error_message="R failed",
        ),
    ]

    summary_path = tmp_path / "runs" / "summary.csv"
    summarize_batch_results(results, summary_path)

    summary = pd.read_csv(summary_path)
    assert summary["location"].tolist() == ["loc_A", "loc_B"]
    assert summary.loc[1, "error_stage"] == "r"


def test_run_batch_locations_writes_summary_with_r_and_python_results(tmp_path: Path) -> None:
    def fake_r_runner(location: str, out_dir: Path) -> tuple[bool, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"sample": ["s1"], "s1": [1.0]}).to_csv(out_dir / "amat.csv", index=False)
        pd.DataFrame({"sample": ["s1"], "pred": [1.0]}).to_csv(out_dir / "kin_blup_predictions.csv", index=False)
        pd.DataFrame({"term": ["Intercept"], "beta": [0.75]}).to_csv(out_dir / "mixed_solve_beta.csv", index=False)
        pd.DataFrame({"sample": ["s1"], "u": [0.25]}).to_csv(out_dir / "mixed_solve_u.csv", index=False)
        pd.DataFrame({"vu": [1.0], "ve": [2.0]}).to_csv(out_dir / "mixed_solve_varcomp.csv", index=False)
        return True, ""

    def fake_python_runner(location: str, out_dir: Path) -> tuple[bool, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"sample": ["s1"], "s1": [1.0]}).to_csv(out_dir / "amat.csv", index=False)
        pd.DataFrame({"sample": ["s1"], "prediction": [1.0]}).to_csv(out_dir / "kin_blup_predictions.csv", index=False)
        pd.DataFrame({"term": ["Intercept"], "beta": [0.75]}).to_csv(out_dir / "mixed_solve_beta.csv", index=False)
        pd.DataFrame({"sample": ["s1"], "u": [0.25]}).to_csv(out_dir / "mixed_solve_u.csv", index=False)
        (out_dir / "mixed_solve_varcomp.json").write_text('{"vu": 1.0, "ve": 2.0}', encoding="utf-8")
        return True, ""

    results = run_batch_locations(
        locations=["loc_A"],
        runs_root=tmp_path / "runs",
        r_runner=fake_r_runner,
        python_runner=fake_python_runner,
    )

    assert results[0].amat_rmse == 0.0
    assert results[0].prediction_rmse == 0.0
    assert results[0].beta_rmse == 0.0
    assert results[0].u_rmse == 0.0


def test_batch_summary_contains_required_columns(tmp_path: Path) -> None:
    results = [
        BatchComparisonResult(
            location="loc_A",
            r_success=True,
            python_success=True,
            sample_count=10,
            amat_rmse=0.1,
            prediction_rmse=0.2,
            beta_rmse=0.25,
            u_rmse=0.3,
            vu_r=1.0,
            vu_py=1.1,
            ve_r=2.0,
            ve_py=2.1,
            error_stage="",
            error_message="",
        )
    ]
    summary_path = tmp_path / "runs" / "summary.csv"
    summarize_batch_results(results, summary_path)
    frame = pd.read_csv(summary_path)
    required = {
        "location",
        "r_success",
        "python_success",
        "sample_count",
        "amat_rmse",
        "prediction_rmse",
        "beta_rmse",
        "u_rmse",
        "vu_r",
        "vu_py",
        "ve_r",
        "ve_py",
        "error_stage",
        "error_message",
    }
    assert required.issubset(frame.columns)
