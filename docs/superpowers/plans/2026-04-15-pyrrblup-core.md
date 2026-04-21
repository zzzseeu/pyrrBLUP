# pyrrBLUP Core Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Python rrBLUP core that can run `A.mat`, `mixed.solve`, and `kin.blup` on the provided real genotype and phenotype files, while comparing outputs against R `rrBLUP` baselines.

**Architecture:** Keep the first milestone small and explicit: one data-loading layer for sample alignment and single-location filtering, one R baseline script that emits stable artifacts, one Python algorithm layer for `A.mat`, `mixed.solve`, and `kin.blup`, and one verification layer for unit tests, smoke tests, and baseline regression tests. All implementation work should happen inside the dedicated worktree at `.worktrees/pyrrblup-core` on branch `feature/pyrrblup-core`.

**Tech Stack:** Python 3, numpy, pandas, scipy, pytest, R, rrBLUP

---

### Task 1: Scaffold The Python Package And Test Harness

**Files:**
- Create: `pyproject.toml`
- Create: `src/pyrrblup/__init__.py`
- Create: `src/pyrrblup/models.py`
- Create: `tests/conftest.py`
- Create: `tests/test_smoke_import.py`

- [ ] **Step 1: Write the failing import test**

```python
from pyrrblup import __version__


def test_package_imports() -> None:
    assert isinstance(__version__, str)
    assert __version__
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_smoke_import.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pyrrblup'`

- [ ] **Step 3: Write minimal project files**

```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "pyrrblup"
version = "0.1.0"
description = "Python rrBLUP alignment project"
requires-python = ">=3.10"
dependencies = [
  "numpy",
  "pandas",
  "scipy",
]

[tool.setuptools]
package-dir = {"" = "src"}

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

```python
# src/pyrrblup/__init__.py
__version__ = "0.1.0"
```

```python
# src/pyrrblup/models.py
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
    vu: float
    ve: float


@dataclass(slots=True)
class KinBlupResult:
    predictions: pd.Series
    mixed_solve: MixedSolveResult
    kinship: pd.DataFrame
```

```python
# tests/conftest.py
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_smoke_import.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/pyrrblup/__init__.py src/pyrrblup/models.py tests/conftest.py tests/test_smoke_import.py
git commit -m "chore: scaffold pyrrblup package"
```

### Task 2: Build The R Baseline Export Script

**Files:**
- Create: `scripts/build_r_baseline.R`
- Create: `baselines/.gitkeep`
- Create: `tests/test_baseline_artifacts_exist.py`

- [ ] **Step 1: Write the failing baseline artifact test**

```python
from pathlib import Path


def test_baseline_script_exists() -> None:
    assert Path("scripts/build_r_baseline.R").exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_baseline_artifacts_exist.py -v`
Expected: FAIL because `scripts/build_r_baseline.R` does not exist

- [ ] **Step 3: Write the baseline script and baseline directory marker**

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rrBLUP)
})

args <- commandArgs(trailingOnly = TRUE)
location_col <- if (length(args) >= 1) args[[1]] else "loc_BeiJ"
out_dir <- if (length(args) >= 2) args[[2]] else "baselines/current"

geno_path <- "../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv"
pheno_path <- "../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

geno <- read.csv(geno_path, check.names = FALSE)
pheno <- read.csv(pheno_path, check.names = FALSE)

stopifnot("sample" %in% names(geno))
stopifnot(all(c("sample", "value", location_col) %in% names(pheno)))

pheno_loc <- pheno[pheno[[location_col]] == 1, c("sample", "value")]
common <- intersect(pheno_loc$sample, geno$sample)
stopifnot(length(common) > 1)

pheno_aligned <- pheno_loc[match(common, pheno_loc$sample), , drop = FALSE]
geno_aligned <- geno[match(common, geno$sample), , drop = FALSE]
rownames(geno_aligned) <- geno_aligned$sample
M <- as.matrix(geno_aligned[, setdiff(names(geno_aligned), "sample"), drop = FALSE])
A <- A.mat(M)
ms <- mixed.solve(y = pheno_aligned$value, K = A)
kb <- kin.blup(data = pheno_aligned, geno = "sample", pheno = "value", K = A)

write.csv(data.frame(sample = common), file.path(out_dir, "samples.csv"), row.names = FALSE)
write.csv(data.frame(sample = pheno_aligned$sample, value = pheno_aligned$value), file.path(out_dir, "phenotype.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(A), A, check.names = FALSE), file.path(out_dir, "amat.csv"), row.names = FALSE)
write.csv(data.frame(beta = ms$beta), file.path(out_dir, "mixed_solve_beta.csv"), row.names = FALSE)
write.csv(data.frame(sample = names(ms$u), u = as.numeric(ms$u)), file.path(out_dir, "mixed_solve_u.csv"), row.names = FALSE)
write.csv(data.frame(vu = ms$Vu, ve = ms$Ve), file.path(out_dir, "mixed_solve_varcomp.csv"), row.names = FALSE)
write.csv(data.frame(sample = names(kb$g.pred), pred = as.numeric(kb$g.pred)), file.path(out_dir, "kin_blup_predictions.csv"), row.names = FALSE)
writeLines(location_col, file.path(out_dir, "location.txt"))
```

```text
# baselines/.gitkeep
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_baseline_artifacts_exist.py -v`
Expected: PASS

- [ ] **Step 5: Generate baseline artifacts in the R environment**

Run: `conda activate rrblup && Rscript scripts/build_r_baseline.R loc_BeiJ baselines/loc_BeiJ`
Expected: baseline CSV files appear under `baselines/loc_BeiJ/`

- [ ] **Step 6: Commit**

```bash
git add scripts/build_r_baseline.R baselines/.gitkeep tests/test_baseline_artifacts_exist.py
git commit -m "feat: add rrBLUP baseline export script"
```

### Task 3: Implement Data Loading, Location Filtering, And Sample Alignment

**Files:**
- Create: `src/pyrrblup/io.py`
- Create: `tests/test_io.py`

- [ ] **Step 1: Write the failing IO tests**

```python
import pandas as pd
import pytest

from pyrrblup.io import align_single_location


def test_align_single_location_filters_and_aligns() -> None:
    geno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3"],
            "m1": [0, 1, 2],
            "m2": [2, 1, 0],
        }
    )
    pheno = pd.DataFrame(
        {
            "sample": ["s2", "s1", "s1"],
            "value": [1.5, 2.5, 9.9],
            "loc_A": [1, 1, 0],
        }
    )

    result = align_single_location(geno, pheno, "loc_A")

    assert list(result.samples) == ["s2", "s1"]
    assert list(result.phenotype) == [1.5, 2.5]
    assert list(result.markers.index) == ["s2", "s1"]


def test_align_single_location_rejects_missing_location() -> None:
    geno = pd.DataFrame({"sample": ["s1"], "m1": [0]})
    pheno = pd.DataFrame({"sample": ["s1"], "value": [1.0], "loc_A": [1]})

    with pytest.raises(ValueError, match="Unknown location column"):
        align_single_location(geno, pheno, "loc_B")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_io.py -v`
Expected: FAIL with `ModuleNotFoundError` or missing function error

- [ ] **Step 3: Write the minimal implementation**

```python
from __future__ import annotations

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
        available = [c for c in phenotype_df.columns if c.startswith("loc_")]
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


def load_and_align(
    genotype_path: str,
    phenotype_path: str,
    location_col: str,
) -> AlignmentResult:
    genotype_df = pd.read_csv(genotype_path)
    phenotype_df = pd.read_csv(phenotype_path)
    return align_single_location(genotype_df, phenotype_df, location_col)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_io.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/io.py tests/test_io.py
git commit -m "feat: add genotype phenotype alignment"
```

### Task 4: Implement Python `A.mat`

**Files:**
- Create: `src/pyrrblup/amat.py`
- Create: `tests/test_amat.py`

- [ ] **Step 1: Write the failing `A.mat` tests**

```python
import numpy as np
import pandas as pd

from pyrrblup.amat import amat


def test_amat_returns_symmetric_dataframe() -> None:
    markers = pd.DataFrame(
        {"m1": [0, 1, 2], "m2": [2, 1, 0]},
        index=pd.Index(["s1", "s2", "s3"], name="sample"),
    )

    result = amat(markers)

    assert list(result.matrix.index) == ["s1", "s2", "s3"]
    assert np.allclose(result.matrix.values, result.matrix.values.T)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_amat.py -v`
Expected: FAIL because `pyrrblup.amat` does not exist

- [ ] **Step 3: Write the minimal implementation**

```python
from __future__ import annotations

import numpy as np
import pandas as pd

from pyrrblup.models import AMatResult


def amat(markers: pd.DataFrame) -> AMatResult:
    if markers.empty:
        raise ValueError("Marker matrix is empty")

    M = markers.astype(float).to_numpy(copy=True)
    allele_freq = M.mean(axis=0) / 2.0
    centered = M - (2.0 * allele_freq)
    denom = 2.0 * np.sum(allele_freq * (1.0 - allele_freq))
    if denom <= 0:
        raise ValueError("A.mat denominator is not positive")

    A = centered @ centered.T / denom
    matrix = pd.DataFrame(A, index=markers.index.copy(), columns=markers.index.copy())
    return AMatResult(matrix=matrix)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_amat.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/amat.py tests/test_amat.py
git commit -m "feat: implement additive relationship matrix"
```

### Task 5: Implement Python `mixed.solve`

**Files:**
- Create: `src/pyrrblup/mixed_solve.py`
- Create: `tests/test_mixed_solve.py`

- [ ] **Step 1: Write the failing mixed-solve tests**

```python
import numpy as np

from pyrrblup.mixed_solve import mixed_solve


def test_mixed_solve_returns_expected_shapes() -> None:
    y = np.array([1.0, 2.0, 3.0])
    k = np.eye(3)

    result = mixed_solve(y=y, k=k)

    assert result.beta.shape == (1,)
    assert result.u.shape == (3,)
    assert result.vu >= 0.0
    assert result.ve >= 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mixed_solve.py -v`
Expected: FAIL because `pyrrblup.mixed_solve` does not exist

- [ ] **Step 3: Write the minimal implementation**

```python
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

from pyrrblup.models import MixedSolveResult


def _neg_reml(log_lambda: float, y: np.ndarray, x: np.ndarray, k: np.ndarray) -> float:
    lam = float(np.exp(log_lambda))
    v = k + lam * np.eye(k.shape[0])
    sign_v, logdet_v = np.linalg.slogdet(v)
    if sign_v <= 0:
        return np.inf

    v_inv = np.linalg.inv(v)
    xt_vinv = x.T @ v_inv
    xt_vinv_x = xt_vinv @ x
    sign_x, logdet_x = np.linalg.slogdet(xt_vinv_x)
    if sign_x <= 0:
        return np.inf

    beta = np.linalg.solve(xt_vinv_x, xt_vinv @ y)
    resid = y - x @ beta
    quad = float(resid.T @ v_inv @ resid)
    n = y.shape[0]
    p = x.shape[1]
    return 0.5 * (logdet_v + logdet_x + (n - p) * np.log(quad))


def mixed_solve(y: np.ndarray, k: np.ndarray, x: np.ndarray | None = None) -> MixedSolveResult:
    y = np.asarray(y, dtype=float).reshape(-1)
    k = np.asarray(k, dtype=float)
    if x is None:
        x = np.ones((y.shape[0], 1), dtype=float)
    else:
        x = np.asarray(x, dtype=float)

    opt = minimize_scalar(_neg_reml, bounds=(-8.0, 8.0), method="bounded", args=(y, x, k))
    if not opt.success:
        raise RuntimeError(f"REML optimization failed: {opt.message}")

    lam = float(np.exp(opt.x))
    v = k + lam * np.eye(k.shape[0])
    v_inv = np.linalg.inv(v)
    beta = np.linalg.solve(x.T @ v_inv @ x, x.T @ v_inv @ y)
    resid = y - x @ beta
    vu = float((resid.T @ v_inv @ resid) / max(y.shape[0] - x.shape[1], 1))
    ve = float(vu * lam)
    u = vu * (k @ v_inv @ resid)
    return MixedSolveResult(beta=beta, u=u, vu=vu, ve=ve)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mixed_solve.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/mixed_solve.py tests/test_mixed_solve.py
git commit -m "feat: implement mixed solve core"
```

### Task 6: Implement Python `kin.blup` And CLI Pipeline

**Files:**
- Create: `src/pyrrblup/kin_blup.py`
- Create: `scripts/run_python_pipeline.py`
- Create: `tests/test_kin_blup.py`

- [ ] **Step 1: Write the failing `kin.blup` tests**

```python
import pandas as pd

from pyrrblup.kin_blup import kin_blup_from_frames


def test_kin_blup_returns_predictions_for_aligned_samples() -> None:
    geno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3"],
            "m1": [0, 1, 2],
            "m2": [2, 1, 0],
        }
    )
    pheno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3"],
            "value": [1.0, 2.0, 3.0],
            "loc_A": [1, 1, 1],
        }
    )

    result = kin_blup_from_frames(geno, pheno, "loc_A")

    assert list(result.predictions.index) == ["s1", "s2", "s3"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_kin_blup.py -v`
Expected: FAIL because `pyrrblup.kin_blup` does not exist

- [ ] **Step 3: Write the minimal implementation**

```python
from __future__ import annotations

from pathlib import Path

import pandas as pd

from pyrrblup.amat import amat
from pyrrblup.io import align_single_location, load_and_align
from pyrrblup.mixed_solve import mixed_solve
from pyrrblup.models import KinBlupResult


def kin_blup_from_frames(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    location_col: str,
) -> KinBlupResult:
    aligned = align_single_location(genotype_df, phenotype_df, location_col)
    kinship = amat(aligned.markers).matrix
    fit = mixed_solve(y=aligned.phenotype.to_numpy(), k=kinship.to_numpy())
    predictions = pd.Series(fit.beta[0] + fit.u, index=aligned.samples, name="prediction")
    return KinBlupResult(predictions=predictions, mixed_solve=fit, kinship=kinship)


def kin_blup_from_paths(genotype_path: str, phenotype_path: str, location_col: str) -> KinBlupResult:
    aligned = load_and_align(genotype_path, phenotype_path, location_col)
    kinship = amat(aligned.markers).matrix
    fit = mixed_solve(y=aligned.phenotype.to_numpy(), k=kinship.to_numpy())
    predictions = pd.Series(fit.beta[0] + fit.u, index=aligned.samples, name="prediction")
    return KinBlupResult(predictions=predictions, mixed_solve=fit, kinship=kinship)
```

```python
#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pyrrblup.kin_blup import kin_blup_from_paths


def main() -> int:
    location_col = sys.argv[1] if len(sys.argv) > 1 else "loc_BeiJ"
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "baselines" / "python_loc_BeiJ"
    out_dir.mkdir(parents=True, exist_ok=True)

    result = kin_blup_from_paths(
        "../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv",
        "../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv",
        location_col,
    )

    result.kinship.to_csv(out_dir / "amat.csv")
    result.predictions.to_csv(out_dir / "kin_blup_predictions.csv", header=True)
    with (out_dir / "mixed_solve_varcomp.json").open("w", encoding="utf-8") as handle:
        json.dump({"vu": result.mixed_solve.vu, "ve": result.mixed_solve.ve}, handle, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_kin_blup.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyrrblup/kin_blup.py scripts/run_python_pipeline.py tests/test_kin_blup.py
git commit -m "feat: add kin blup pipeline"
```

### Task 7: Add Real-Data Smoke Tests And R Baseline Regression Tests

**Files:**
- Create: `tests/test_against_r_baseline.py`

- [ ] **Step 1: Write the failing regression tests**

```python
from pathlib import Path

import numpy as np
import pandas as pd

from pyrrblup.kin_blup import kin_blup_from_paths


def test_r_baseline_files_exist() -> None:
    baseline_dir = Path("baselines/loc_BeiJ")
    assert (baseline_dir / "amat.csv").exists()
    assert (baseline_dir / "kin_blup_predictions.csv").exists()


def test_python_pipeline_matches_r_sample_order() -> None:
    baseline_samples = pd.read_csv("baselines/loc_BeiJ/samples.csv")
    result = kin_blup_from_paths(
        "../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv",
        "../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv",
        "loc_BeiJ",
    )
    assert list(result.predictions.index) == baseline_samples["sample"].tolist()


def test_python_pipeline_has_positive_prediction_correlation_with_r() -> None:
    baseline_pred = pd.read_csv("baselines/loc_BeiJ/kin_blup_predictions.csv")
    result = kin_blup_from_paths(
        "../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv",
        "../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv",
        "loc_BeiJ",
    )
    merged = baseline_pred.merge(
        result.predictions.rename("prediction").rename_axis("sample").reset_index(),
        on="sample",
        suffixes=("_r", "_py"),
    )
    corr = np.corrcoef(merged["pred_r"], merged["prediction"])[0, 1]
    assert corr > 0.5
```

- [ ] **Step 2: Run test to verify it fails before baseline generation or Python completion**

Run: `pytest tests/test_against_r_baseline.py -v`
Expected: FAIL due to missing baseline files or missing implementation details

- [ ] **Step 3: Adjust Python and baseline serialization until the tests pass**

```python
# Example normalization rule to apply if needed in tests or implementation:
# - always sort or reindex by baseline sample order before numeric comparison
# - compare shape first, then correlation, then elementwise values
```

- [ ] **Step 4: Run the targeted regression test suite**

Run: `pytest tests/test_against_r_baseline.py -v`
Expected: PASS for the initial loose regression criteria

- [ ] **Step 5: Run the full test suite**

Run: `pytest -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add tests/test_against_r_baseline.py scripts/run_python_pipeline.py src/pyrrblup/*.py
git commit -m "test: add rrBLUP baseline regression coverage"
```

### Task 8: Run End-To-End Baseline And Python Pipelines In Their Target Environments

**Files:**
- Modify: `scripts/build_r_baseline.R`
- Modify: `scripts/run_python_pipeline.py`
- Modify: `tests/test_against_r_baseline.py`

- [ ] **Step 1: Run the R baseline command in the `rrblup` environment**

Run: `conda activate rrblup && Rscript scripts/build_r_baseline.R loc_BeiJ baselines/loc_BeiJ`
Expected: baseline artifacts are refreshed successfully

- [ ] **Step 2: Run the Python pipeline command in the `evo2` environment**

Run: `conda activate evo2 && python scripts/run_python_pipeline.py loc_BeiJ baselines/python_loc_BeiJ`
Expected: Python output artifacts are written successfully

- [ ] **Step 3: Run the final regression test in the `evo2` environment**

Run: `conda activate evo2 && pytest tests/test_against_r_baseline.py -v`
Expected: PASS

- [ ] **Step 4: Commit any final alignment fixes**

```bash
git add scripts/build_r_baseline.R scripts/run_python_pipeline.py tests/test_against_r_baseline.py src/pyrrblup/*.py baselines/loc_BeiJ
git commit -m "feat: align python rrBLUP workflow with R baseline"
```
