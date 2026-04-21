from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyrrblup.amat import amat
from pyrrblup.io import load_and_align, resolve_project_data_path
from pyrrblup.kin_blup import kin_blup_from_paths


GENOTYPE_PATH = str(resolve_project_data_path("../../project/genomic_selection_project/data/processed/model_dataset/agront_ds_202603/genotype_012.csv"))
PHENOTYPE_PATH = str(resolve_project_data_path("../../project/genomic_selection_project/data/processed/phenotype/phenotype_ETN_processed_20260227.csv"))

EXTERNAL_DATA_AVAILABLE = Path(GENOTYPE_PATH).exists() and Path(PHENOTYPE_PATH).exists()


pytestmark = pytest.mark.skipif(
    not EXTERNAL_DATA_AVAILABLE,
    reason="external genotype/phenotype data is not available in this environment",
)


def test_r_baseline_files_exist() -> None:
    baseline_dir = Path("baselines/loc_BeiJ")
    assert (baseline_dir / "amat.csv").exists()
    assert (baseline_dir / "kin_blup_predictions.csv").exists()


def test_python_amat_matches_r_baseline_structure_and_values() -> None:
    baseline_amat = pd.read_csv("baselines/loc_BeiJ/amat.csv")
    aligned = load_and_align(GENOTYPE_PATH, PHENOTYPE_PATH, "loc_BeiJ")
    result = amat(aligned.markers).matrix

    assert list(result.index) == baseline_amat["sample"].tolist()
    assert list(result.columns) == baseline_amat["sample"].tolist()

    baseline_values = baseline_amat.drop(columns=["sample"]).to_numpy()
    result_values = result.to_numpy()
    diff = baseline_values - result_values

    assert np.corrcoef(baseline_values.ravel(), result_values.ravel())[0, 1] > 0.99
    assert np.sqrt(np.mean(diff**2)) < 0.05


def test_python_pipeline_matches_r_sample_order() -> None:
    baseline_samples = pd.read_csv("baselines/loc_BeiJ/samples.csv")
    result = kin_blup_from_paths(GENOTYPE_PATH, PHENOTYPE_PATH, "loc_BeiJ")
    assert list(result.predictions.index) == baseline_samples["sample"].tolist()


def test_python_pipeline_matches_r_predictions_numerically() -> None:
    baseline_pred = pd.read_csv("baselines/loc_BeiJ/kin_blup_predictions.csv")
    result = kin_blup_from_paths(GENOTYPE_PATH, PHENOTYPE_PATH, "loc_BeiJ")
    merged = baseline_pred.merge(
        result.predictions.rename("prediction").rename_axis("sample").reset_index(),
        on="sample",
    )
    diff = merged["pred"] - merged["prediction"]
    assert np.corrcoef(merged["pred"], merged["prediction"])[0, 1] > 0.99
    assert np.sqrt(np.mean(diff**2)) < 0.05


def test_python_pipeline_matches_r_fixed_effects_numerically() -> None:
    baseline_beta = pd.read_csv("baselines/loc_BeiJ/mixed_solve_beta.csv")
    result = kin_blup_from_paths(GENOTYPE_PATH, PHENOTYPE_PATH, "loc_BeiJ")
    python_beta = pd.DataFrame({"beta": result.mixed_solve.beta})
    merged = baseline_beta.merge(
        python_beta.reset_index(names="row_id"),
        left_index=True,
        right_index=True,
        suffixes=("_r", "_py"),
    )
    diff = merged["beta_r"] - merged["beta_py"]
    assert np.sqrt(np.mean(diff**2)) < 0.05


def test_python_pipeline_matches_r_random_effects_numerically() -> None:
    baseline_u = pd.read_csv("baselines/loc_BeiJ/mixed_solve_u.csv")
    result = kin_blup_from_paths(GENOTYPE_PATH, PHENOTYPE_PATH, "loc_BeiJ")
    merged = baseline_u.merge(
        pd.DataFrame({"sample": result.predictions.index, "u": result.mixed_solve.u}),
        on="sample",
        suffixes=("_r", "_py"),
    )
    diff = merged["u_r"] - merged["u_py"]
    assert np.corrcoef(merged["u_r"], merged["u_py"])[0, 1] > 0.99
    assert np.sqrt(np.mean(diff**2)) < 0.05


def test_python_pipeline_handles_near_zero_vu_boundary_for_loc_jiangx() -> None:
    baseline_pred = pd.read_csv("baselines/loc_JiangX/kin_blup_predictions.csv")
    result = kin_blup_from_paths(GENOTYPE_PATH, PHENOTYPE_PATH, "loc_JiangX")
    merged = baseline_pred.merge(
        result.predictions.rename("prediction").rename_axis("sample").reset_index(),
        on="sample",
    )
    diff = merged["pred"] - merged["prediction"]
    assert np.sqrt(np.mean(diff**2)) < 1e-4
