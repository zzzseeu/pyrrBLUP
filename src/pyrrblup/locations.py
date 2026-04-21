from __future__ import annotations

import pandas as pd


def discover_locations(phenotype_df: pd.DataFrame) -> list[str]:
    locations = [column for column in phenotype_df.columns if column.startswith("loc_")]
    if not locations:
        raise ValueError("No phenotype location columns with prefix 'loc_' were found")
    return locations


def filter_locations(discovered: list[str], requested: list[str] | None) -> list[str]:
    if not requested:
        return discovered
    missing = sorted(set(requested).difference(discovered))
    if missing:
        raise ValueError(f"Unknown requested locations: {missing}")
    filtered = [location for location in discovered if location in requested]
    if not filtered:
        raise ValueError("No locations remain after applying the requested filter")
    return filtered
