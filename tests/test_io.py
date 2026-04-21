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
