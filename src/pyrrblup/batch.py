from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import pandas as pd

from pyrrblup.compare import compare_location_outputs
from pyrrblup.models import BatchComparisonResult, BatchRunStatus
from pyrrblup.status import write_status_json

RRunner = Callable[[str, Path], tuple[bool, str]]
PythonRunner = Callable[[str, Path], tuple[bool, str]]


def summarize_batch_results(results: list[BatchComparisonResult], summary_path: Path) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([asdict(result) for result in results])
    frame.to_csv(summary_path, index=False)


def run_batch_locations(
    locations: list[str],
    runs_root: Path,
    r_runner: RRunner,
    python_runner: PythonRunner,
) -> list[BatchComparisonResult]:
    results: list[BatchComparisonResult] = []
    for location in locations:
        r_dir = runs_root / location / "r"
        python_dir = runs_root / location / "python"

        r_started = datetime.now(timezone.utc).isoformat()
        r_success, r_error = r_runner(location, r_dir)
        r_finished = datetime.now(timezone.utc).isoformat()
        write_status_json(
            BatchRunStatus(location, "r", r_success, r_error, r_started, r_finished),
            r_dir / "status.json",
        )

        python_started = datetime.now(timezone.utc).isoformat()
        python_success, python_error = python_runner(location, python_dir)
        python_finished = datetime.now(timezone.utc).isoformat()
        write_status_json(
            BatchRunStatus(location, "python", python_success, python_error, python_started, python_finished),
            python_dir / "status.json",
        )

        if r_success and python_success:
            results.append(compare_location_outputs(location, r_dir, python_dir, True, True))
        else:
            results.append(
                BatchComparisonResult(
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
                    error_stage="r" if not r_success else "python",
                    error_message=r_error if not r_success else python_error,
                )
            )
    return results
