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
    pd.DataFrame({"term": ["Intercept"], "beta": [1.5]}).to_csv(r_dir / "mixed_solve_beta.csv", index=False)
    pd.DataFrame({"sample": ["s1", "s2"], "u": [0.5, -0.5]}).to_csv(r_dir / "mixed_solve_u.csv", index=False)
    pd.DataFrame({"vu": [3.0], "ve": [4.0]}).to_csv(r_dir / "mixed_solve_varcomp.csv", index=False)

    pd.DataFrame({"sample": ["s1", "s2"], "s1": [1.0, 0.1], "s2": [0.1, 1.0]}).to_csv(py_dir / "amat.csv", index=False)
    pd.DataFrame({"sample": ["s1", "s2"], "prediction": [1.0, 2.0]}).to_csv(py_dir / "kin_blup_predictions.csv", index=False)
    pd.DataFrame({"term": ["Intercept"], "beta": [1.5]}).to_csv(py_dir / "mixed_solve_beta.csv", index=False)
    pd.DataFrame({"sample": ["s1", "s2"], "u": [0.5, -0.5]}).to_csv(py_dir / "mixed_solve_u.csv", index=False)
    (py_dir / "mixed_solve_varcomp.json").write_text('{"vu": 3.0, "ve": 4.0}', encoding="utf-8")

    result = compare_location_outputs("loc_A", r_dir, py_dir, r_success=True, python_success=True)

    assert result.location == "loc_A"
    assert result.amat_rmse == 0.0
    assert result.prediction_rmse == 0.0
    assert result.beta_rmse == 0.0
    assert result.u_rmse == 0.0
    assert result.vu_r == 3.0
    assert result.vu_py == 3.0
