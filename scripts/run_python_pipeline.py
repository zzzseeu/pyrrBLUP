#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pyrrblup.io import resolve_project_data_path
from pyrrblup.kin_blup import kin_blup


def main() -> int:
    location_col = sys.argv[1] if len(sys.argv) > 1 else "loc_BeiJ"
    out_dir = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else ROOT / "baselines" / "python_loc_BeiJ"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    result = kin_blup(
        str(
            resolve_project_data_path(
                "/share/org/YZWL/yzbsl_zhangchao/project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv"
            )
        ),
        str(
            resolve_project_data_path(
                "/share/org/YZWL/yzbsl_zhangchao/project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv"
            )
        ),
        location_col,
    )

    result.kinship.to_csv(out_dir / "amat.csv")
    result.predictions.to_csv(out_dir / "kin_blup_predictions.csv", header=True)
    pd.DataFrame({"term": ["Intercept"], "beta": result.mixed_solve.beta}).to_csv(
        out_dir / "mixed_solve_beta.csv", index=False
    )
    pd.DataFrame({"sample": result.predictions.index, "u": result.mixed_solve.u}).to_csv(
        out_dir / "mixed_solve_u.csv", index=False
    )
    with (out_dir / "mixed_solve_varcomp.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {"vu": result.mixed_solve.Vu, "ve": result.mixed_solve.Ve}, handle, indent=2
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
