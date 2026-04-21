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
