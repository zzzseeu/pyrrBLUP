# pyrrBLUP Batch Multi-Location Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a batch orchestration layer that discovers phenotype locations, runs per-location R baselines and Python outputs into `runs/<location>/r/` and `runs/<location>/python/`, records per-tool status, and writes a cross-location summary.

**Architecture:** Keep the current single-location R and Python entrypoints as the stable core, then add a lightweight orchestration layer for location discovery, status tracking, per-location comparison, and summary generation. Batch mode should continue after per-location failures, preserve per-location artifacts, and return a non-zero exit code when any location fails.

**Tech Stack:** Python 3, numpy, pandas, scipy, pytest, R, rrBLUP, json, csv

---

### Task 1: Add Batch Models And Status Serialization

**Files:**
- Modify: `src/pyrrblup/models.py`
- Create: `src/pyrrblup/status.py`
- Create: `tests/test_status.py`

- [ ] **Step 1: Write the failing status serialization test**

```python
from pathlib import Path

from pyrrblup.models import BatchRunStatus
from pyrrblup.status import write_status_json


def test_write_status_json_creates_expected_fields(tmp_path: Path) -> None:
    status = BatchRunStatus(
        location="loc_A",
        tool="python",
        success=True,
        error_message="",
        started_at="2026-04-16T10:00:00",
        finished_at="2026-04-16T10:00:03",
    )

    out_path = tmp_path / "status.json"
    write_status_json(status, out_path)

    payload = out_path.read_text(encoding="utf-8")
    assert '"location": "loc_A"' in payload
    assert '"tool": "python"' in payload
    assert '"success": true' in payload
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_status.py -v`
Expected: FAIL with `ImportError` for `BatchRunStatus` or `write_status_json`

- [ ] **Step 3: Write minimal batch status models and serializer**

```python
# src/pyrrblup/models.py
from dataclasses import asdict, dataclass

# keep existing imports and classes above


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
    vu_r: float | None
    vu_py: float | None
    ve_r: float | None
    ve_py: float | None
    error_stage: str
    error_message: str
```

```python
# src/pyrrblup/status.py
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from pyrrblup.models import BatchRunStatus


def write_status_json(status: BatchRunStatus, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(asdict(status), indent=2), encoding="utf-8")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_status.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/models.py src/pyrrblup/status.py tests/test_status.py
git commit -m "feat: add batch status models"
```

### Task 2: Add Location Discovery And Filtering

**Files:**
- Create: `src/pyrrblup/locations.py`
- Create: `tests/test_locations.py`

- [ ] **Step 1: Write the failing location discovery tests**

```python
import pandas as pd
import pytest

from pyrrblup.locations import discover_locations, filter_locations


def test_discover_locations_returns_loc_columns_in_order() -> None:
    phenotype = pd.DataFrame(
        {
            "sample": ["s1"],
            "value": [1.0],
            "loc_B": [1],
            "misc": [0],
            "loc_A": [0],
        }
    )

    assert discover_locations(phenotype) == ["loc_B", "loc_A"]


def test_filter_locations_rejects_unknown_requested_location() -> None:
    with pytest.raises(ValueError, match="Unknown requested locations"):
        filter_locations(["loc_A", "loc_B"], ["loc_C"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_locations.py -v`
Expected: FAIL with `ImportError` for `pyrrblup.locations`

- [ ] **Step 3: Write minimal location helpers**

```python
# src/pyrrblup/locations.py
from __future__ import annotations

import pandas as pd



def discover_locations(phenotype_df: pd.DataFrame) -> list[str]:
    locations = [col for col in phenotype_df.columns if col.startswith("loc_")]
    if not locations:
        raise ValueError("No phenotype location columns with prefix 'loc_' were found")
    return locations



def filter_locations(discovered: list[str], requested: list[str] | None) -> list[str]:
    if not requested:
        return discovered
    missing = sorted(set(requested).difference(discovered))
    if missing:
        raise ValueError(f"Unknown requested locations: {missing}")
    filtered = [loc for loc in discovered if loc in requested]
    if not filtered:
        raise ValueError("No locations remain after applying the requested filter")
    return filtered
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_locations.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/locations.py tests/test_locations.py
git commit -m "feat: add location discovery helpers"
```

### Task 3: Add Per-Location Comparison Logic

**Files:**
- Create: `src/pyrrblup/compare.py`
- Create: `tests/test_compare.py`

- [ ] **Step 1: Write the failing comparison test**

```python
from pathlib import Path

import pandas as pd

from pyrrblup.compare import compare_location_outputs


def test_compare_location_outputs_returns_metrics(tmp_path: Path) -> None:
    r_dir = tmp_path / "loc_A" / "r"
    py_dir = tmp_path / "loc_A" / "python"
    r_dir.mkdir(parents=True)
    py_dir.mkdir(parents=True)

    pd.DataFrame({"sample": ["s1", "s2"], "s1": [1.0, 0.1], "s2": [0.1, 1.0]}).to_csv(r_dir / "amat.csv", index=False)
    pd.DataFrame({"sample": ["s1", "s2"], "pred": [1.0, 2.0]}).to_csv(r_dir / "kin_blup_predictions.csv", index=False)
    pd.DataFrame({"vu": [3.0], "ve": [4.0]}).to_csv(r_dir / "mixed_solve_varcomp.csv", index=False)

    pd.DataFrame({"sample": ["s1", "s2"], "s1": [1.0, 0.1], "s2": [0.1, 1.0]}).to_csv(py_dir / "amat.csv", index=False)
    pd.DataFrame({"sample": ["s1", "s2"], "prediction": [1.0, 2.0]}).to_csv(py_dir / "kin_blup_predictions.csv", index=False)
    (py_dir / "mixed_solve_varcomp.json").write_text('{"vu": 3.0, "ve": 4.0}', encoding="utf-8")

    result = compare_location_outputs("loc_A", r_dir, py_dir, r_success=True, python_success=True)

    assert result.location == "loc_A"
    assert result.amat_rmse == 0.0
    assert result.prediction_rmse == 0.0
    assert result.vu_r == 3.0
    assert result.vu_py == 3.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_compare.py -v`
Expected: FAIL with `ImportError` for `pyrrblup.compare`

- [ ] **Step 3: Write minimal comparison implementation**

```python
# src/pyrrblup/compare.py
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from pyrrblup.models import BatchComparisonResult



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
            vu_r=None,
            vu_py=None,
            ve_r=None,
            ve_py=None,
            error_stage="comparison",
            error_message="comparison skipped because one side failed",
        )

    r_amat = pd.read_csv(r_dir / "amat.csv").set_index("sample")
    p_amat = pd.read_csv(python_dir / "amat.csv")
    if p_amat.columns[0] != "sample":
        p_amat = p_amat.rename(columns={p_amat.columns[0]: "sample"})
    p_amat = p_amat.set_index("sample").loc[r_amat.index, r_amat.columns]
    amat_diff = r_amat.to_numpy(float) - p_amat.to_numpy(float)

    r_pred = pd.read_csv(r_dir / "kin_blup_predictions.csv")
    p_pred = pd.read_csv(python_dir / "kin_blup_predictions.csv")
    pred_merged = r_pred.merge(p_pred, on="sample")
    pred_diff = pred_merged["pred"] - pred_merged[p_pred.columns[-1]]

    r_var = pd.read_csv(r_dir / "mixed_solve_varcomp.csv").iloc[0]
    p_var = json.loads((python_dir / "mixed_solve_varcomp.json").read_text(encoding="utf-8"))

    return BatchComparisonResult(
        location=location,
        r_success=True,
        python_success=True,
        sample_count=int(len(pred_merged)),
        amat_rmse=float(np.sqrt(np.mean(amat_diff**2))),
        prediction_rmse=float(np.sqrt(np.mean(pred_diff**2))),
        vu_r=float(r_var["vu"]),
        vu_py=float(p_var["vu"]),
        ve_r=float(r_var["ve"]),
        ve_py=float(p_var["ve"]),
        error_stage="",
        error_message="",
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_compare.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/compare.py tests/test_compare.py
git commit -m "feat: add per-location comparison metrics"
```

### Task 4: Add Batch Orchestration Core

**Files:**
- Create: `src/pyrrblup/batch.py`
- Create: `tests/test_batch.py`

- [ ] **Step 1: Write the failing batch continuation test**

```python
from pathlib import Path

import pandas as pd

from pyrrblup.batch import summarize_batch_results
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: FAIL with `ImportError` for `pyrrblup.batch`

- [ ] **Step 3: Write minimal batch summary implementation**

```python
# src/pyrrblup/batch.py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pandas as pd

from pyrrblup.models import BatchComparisonResult



def summarize_batch_results(results: list[BatchComparisonResult], summary_path: Path) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([asdict(result) for result in results])
    frame.to_csv(summary_path, index=False)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/batch.py tests/test_batch.py
git commit -m "feat: add batch summary writer"
```

### Task 5: Add Python Batch Entry Point

**Files:**
- Modify: `scripts/run_python_pipeline.py`
- Create: `scripts/run_batch_pipeline.py`
- Modify: `src/pyrrblup/batch.py`
- Test: `tests/test_batch.py`

- [ ] **Step 1: Extend the failing batch test to exercise continue-on-failure behavior**

```python
def test_run_batch_locations_continues_after_python_failure(tmp_path: Path) -> None:
    calls: list[str] = []

    def fake_runner(location: str, out_dir: Path) -> tuple[bool, str]:
        calls.append(location)
        if location == "loc_B":
            return False, "python failed"
        return True, ""

    results = run_batch_locations(
        locations=["loc_A", "loc_B", "loc_C"],
        runs_root=tmp_path / "runs",
        python_runner=fake_runner,
    )

    assert calls == ["loc_A", "loc_B", "loc_C"]
    assert [result.python_success for result in results] == [True, False, True]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: FAIL because `run_batch_locations` is not implemented

- [ ] **Step 3: Write minimal batch runner and CLI**

```python
# append to src/pyrrblup/batch.py
from datetime import datetime, UTC
from pathlib import Path
from typing import Callable

from pyrrblup.compare import compare_location_outputs
from pyrrblup.models import BatchComparisonResult, BatchRunStatus
from pyrrblup.status import write_status_json


PythonRunner = Callable[[str, Path], tuple[bool, str]]


def run_batch_locations(
    locations: list[str],
    runs_root: Path,
    python_runner: PythonRunner,
) -> list[BatchComparisonResult]:
    results: list[BatchComparisonResult] = []
    for location in locations:
        python_dir = runs_root / location / "python"
        started_at = datetime.now(UTC).isoformat()
        success, error_message = python_runner(location, python_dir)
        finished_at = datetime.now(UTC).isoformat()
        write_status_json(
            BatchRunStatus(
                location=location,
                tool="python",
                success=success,
                error_message=error_message,
                started_at=started_at,
                finished_at=finished_at,
            ),
            python_dir / "status.json",
        )
        results.append(
            BatchComparisonResult(
                location=location,
                r_success=False,
                python_success=success,
                sample_count=0,
                amat_rmse=None,
                prediction_rmse=None,
                vu_r=None,
                vu_py=None,
                ve_r=None,
                ve_py=None,
                error_stage="python" if not success else "",
                error_message=error_message,
            )
        )
    return results
```

```python
# scripts/run_batch_pipeline.py
#!/usr/bin/env python
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pandas as pd

from pyrrblup.batch import run_batch_locations, summarize_batch_results
from pyrrblup.locations import discover_locations, filter_locations
from pyrrblup.models import BatchComparisonResult


PHENOTYPE_PATH = ROOT / ".." / ".." / "project" / "genomic_selection_project" / "data" / "processed" / "phenotype" / "phenotype_ETN_processed_20260227.csv"


def run_single_python(location: str, out_dir: Path) -> tuple[bool, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "run_python_pipeline.py"), location, str(out_dir)],
        capture_output=True,
        text=True,
    )
    if proc.returncode == 0:
        return True, ""
    return False, (proc.stderr or proc.stdout).strip()


def main() -> int:
    requested = sys.argv[1].split(",") if len(sys.argv) > 1 and sys.argv[1] else None
    runs_root = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "runs"
    phenotype_df = pd.read_csv(PHENOTYPE_PATH)
    discovered = discover_locations(phenotype_df)
    locations = filter_locations(discovered, requested)
    results = run_batch_locations(locations, runs_root, run_single_python)
    summarize_batch_results(results, runs_root / "summary.csv")
    return 1 if any(not result.python_success for result in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/batch.py scripts/run_batch_pipeline.py tests/test_batch.py
git commit -m "feat: add python batch runner"
```

### Task 6: Add R Batch Entry Point And Full Per-Location Comparison Flow

**Files:**
- Create: `scripts/build_r_baseline_batch.R`
- Modify: `src/pyrrblup/batch.py`
- Modify: `scripts/run_batch_pipeline.py`
- Modify: `tests/test_batch.py`

- [ ] **Step 1: Extend the failing batch test to require comparison and summary output**

```python
def test_run_batch_locations_writes_summary_with_r_and_python_results(tmp_path: Path) -> None:
    def fake_r_runner(location: str, out_dir: Path) -> tuple[bool, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"sample": ["s1"], "s1": [1.0]}).to_csv(out_dir / "amat.csv", index=False)
        pd.DataFrame({"sample": ["s1"], "pred": [1.0]}).to_csv(out_dir / "kin_blup_predictions.csv", index=False)
        pd.DataFrame({"vu": [1.0], "ve": [2.0]}).to_csv(out_dir / "mixed_solve_varcomp.csv", index=False)
        return True, ""

    def fake_python_runner(location: str, out_dir: Path) -> tuple[bool, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"sample": ["s1"], "s1": [1.0]}).to_csv(out_dir / "amat.csv", index=False)
        pd.DataFrame({"sample": ["s1"], "prediction": [1.0]}).to_csv(out_dir / "kin_blup_predictions.csv", index=False)
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: FAIL because `run_batch_locations` does not yet support both runners and comparison

- [ ] **Step 3: Implement combined batch execution and R batch script**

```python
# update src/pyrrblup/batch.py signatures and flow
from typing import Callable

RRunner = Callable[[str, Path], tuple[bool, str]]
PythonRunner = Callable[[str, Path], tuple[bool, str]]


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

        r_started = datetime.now(UTC).isoformat()
        r_success, r_error = r_runner(location, r_dir)
        r_finished = datetime.now(UTC).isoformat()
        write_status_json(BatchRunStatus(location, "r", r_success, r_error, r_started, r_finished), r_dir / "status.json")

        py_started = datetime.now(UTC).isoformat()
        py_success, py_error = python_runner(location, python_dir)
        py_finished = datetime.now(UTC).isoformat()
        write_status_json(BatchRunStatus(location, "python", py_success, py_error, py_started, py_finished), python_dir / "status.json")

        if r_success and py_success:
            result = compare_location_outputs(location, r_dir, python_dir, True, True)
        else:
            result = BatchComparisonResult(
                location=location,
                r_success=r_success,
                python_success=py_success,
                sample_count=0,
                amat_rmse=None,
                prediction_rmse=None,
                vu_r=None,
                vu_py=None,
                ve_r=None,
                ve_py=None,
                error_stage="r" if not r_success else "python",
                error_message=r_error if not r_success else py_error,
            )
        results.append(result)
    return results
```

```python
# update scripts/run_batch_pipeline.py

def run_single_r(location: str, out_dir: Path) -> tuple[bool, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        ["Rscript", str(ROOT / "scripts" / "build_r_baseline.R"), location, str(out_dir)],
        capture_output=True,
        text=True,
    )
    if proc.returncode == 0:
        return True, ""
    return False, (proc.stderr or proc.stdout).strip()

# in main():
results = run_batch_locations(locations, runs_root, run_single_r, run_single_python)
summarize_batch_results(results, runs_root / "summary.csv")
return 1 if any((not result.r_success) or (not result.python_success) for result in results) else 0
```

```r
# scripts/build_r_baseline_batch.R
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
location_cols <- strsplit(args[[1]], ",", fixed = TRUE)[[1]]
runs_root <- if (length(args) >= 2) args[[2]] else "runs"

for (location_col in location_cols) {
  out_dir <- file.path(runs_root, location_col, "r")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2("Rscript", c("scripts/build_r_baseline.R", location_col, out_dir))
  if (status != 0) {
    quit(status = status)
  }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: PASS

- [ ] **Step 5: Run a small real-data smoke batch for two locations**

Run: `source ~/.bashrc && conda activate evo2 && python scripts/run_batch_pipeline.py loc_BeiJ,loc_FuJ runs_smoke`
Expected: `runs_smoke/summary.csv` exists and contains both locations

- [ ] **Step 6: Commit**

```bash
git add src/pyrrblup/batch.py scripts/run_batch_pipeline.py scripts/build_r_baseline_batch.R tests/test_batch.py
git commit -m "feat: add multi-location batch orchestration"
```

### Task 7: Add Real-Data Batch Validation

**Files:**
- Modify: `tests/test_batch.py`
- Modify: `scripts/run_batch_pipeline.py`

- [ ] **Step 1: Write the failing real-data smoke validation test**

```python
from pathlib import Path


def test_batch_summary_contains_required_columns(tmp_path: Path) -> None:
    summary = Path("runs_smoke/summary.csv")
    assert summary.exists()
    frame = pd.read_csv(summary)
    required = {
        "location",
        "r_success",
        "python_success",
        "sample_count",
        "amat_rmse",
        "prediction_rmse",
        "vu_r",
        "vu_py",
        "ve_r",
        "ve_py",
        "error_stage",
        "error_message",
    }
    assert required.issubset(frame.columns)
```

- [ ] **Step 2: Run test to verify it fails if summary is incomplete**

Run: `source ~/.bashrc && conda activate evo2 && python -m pytest tests/test_batch.py -v`
Expected: FAIL until summary columns are complete

- [ ] **Step 3: Adjust summary generation if needed so all required columns are present**

```python
# ensure summarize_batch_results writes every field from BatchComparisonResult
# and preserves missing numeric values as blank or NaN in CSV output
```

- [ ] **Step 4: Run the full test suite**

Run: `source ~/.bashrc && conda activate evo2 && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_batch.py src/pyrrblup/batch.py scripts/run_batch_pipeline.py
git commit -m "test: validate batch summary output"
```
