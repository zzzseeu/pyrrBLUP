import pandas as pd
import pytest

from pyrrblup import kin_blup
from pyrrblup.kin_blup import kin_blup_from_frames, kin_blup_from_paths


def test_kin_blup_returns_predictions_for_aligned_samples() -> None:
    geno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "m1": [0, 0, 1, 1],
            "m2": [0, 1, 0, 1],
            "m3": [0, 1, 1, 0],
        }
    )
    pheno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "value": [1.0, 2.0, 3.0, 2.5],
            "loc_A": [1, 1, 1, 1],
        }
    )

    result = kin_blup_from_frames(geno, pheno, "loc_A")

    assert list(result.predictions.index) == ["s1", "s2", "s3", "s4"]
    assert result.predictions.name == "prediction"
    assert result.kinship.shape == (4, 4)


def test_kin_blup_dispatches_dataframe_inputs() -> None:
    geno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "m1": [0, 0, 1, 1],
            "m2": [0, 1, 0, 1],
            "m3": [0, 1, 1, 0],
        }
    )
    pheno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "value": [1.0, 2.0, 3.0, 2.5],
            "loc_A": [1, 1, 1, 1],
        }
    )

    direct = kin_blup_from_frames(geno, pheno, "loc_A")
    dispatched = kin_blup(geno, pheno, "loc_A")

    assert dispatched.predictions.equals(direct.predictions)
    assert dispatched.kinship.equals(direct.kinship)


def test_kin_blup_dispatches_path_inputs(tmp_path) -> None:
    geno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "m1": [0, 0, 1, 1],
            "m2": [0, 1, 0, 1],
            "m3": [0, 1, 1, 0],
        }
    )
    pheno = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "value": [1.0, 2.0, 3.0, 2.5],
            "loc_A": [1, 1, 1, 1],
        }
    )
    geno_path = tmp_path / "geno.csv"
    pheno_path = tmp_path / "pheno.csv"
    geno.to_csv(geno_path, index=False)
    pheno.to_csv(pheno_path, index=False)

    direct = kin_blup_from_paths(str(geno_path), str(pheno_path), "loc_A")
    dispatched = kin_blup(geno_path, pheno_path, "loc_A")

    assert dispatched.predictions.equals(direct.predictions)
    assert dispatched.kinship.equals(direct.kinship)


def test_kin_blup_rejects_mixed_input_types(tmp_path) -> None:
    geno = pd.DataFrame(
        {
            "sample": ["s1", "s2"],
            "m1": [0, 1],
        }
    )
    pheno = pd.DataFrame(
        {
            "sample": ["s1", "s2"],
            "value": [1.0, 2.0],
            "loc_A": [1, 1],
        }
    )
    pheno_path = tmp_path / "pheno.csv"
    pheno.to_csv(pheno_path, index=False)

    with pytest.raises(TypeError, match="both be paths or both be pandas.DataFrame"):
        kin_blup(geno, pheno_path, "loc_A")
