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
from pyrrblup.io import resolve_project_data_path


PHENOTYPE_PATH = resolve_project_data_path("/share/org/YZWL/yzbsl_zhangchao/project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv")


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
    results = run_batch_locations(locations, runs_root, run_single_r, run_single_python)
    summarize_batch_results(results, runs_root / "summary.csv")
    return 1 if any((not result.r_success) or (not result.python_success) for result in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
