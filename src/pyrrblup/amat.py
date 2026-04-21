from __future__ import annotations

import numpy as np
import pandas as pd

from pyrrblup.models import AMatResult


def amat(markers: pd.DataFrame) -> AMatResult:
    if markers.empty:
        raise ValueError("Marker matrix is empty")

    marker_matrix = markers.astype(float).to_numpy(copy=True)
    n_lines = marker_matrix.shape[0]

    frac_missing = np.isnan(marker_matrix).mean(axis=0)
    freq = np.nanmean(marker_matrix + 1.0, axis=0) / 2.0
    maf = np.minimum(freq, 1.0 - freq)
    min_maf = 1.0 / (2.0 * n_lines)
    max_missing = 1.0 - 1.0 / (2.0 * n_lines)
    keep_markers = (maf >= min_maf) & (frac_missing <= max_missing)
    if not np.any(keep_markers):
        raise ValueError("No markers remain after rrBLUP-style filtering")

    kept_freq = freq[keep_markers]
    w = marker_matrix[:, keep_markers] + 1.0 - 2.0 * kept_freq
    var_a = 2.0 * np.mean(kept_freq * (1.0 - kept_freq))
    if np.isclose(var_a, 0.0):
        raise ValueError("A.mat scaling factor is zero after rrBLUP-style filtering")

    kinship = (w @ w.T) / (var_a * w.shape[1])
    matrix = pd.DataFrame(kinship, index=markers.index.copy(), columns=markers.index.copy())
    return AMatResult(matrix=matrix)
