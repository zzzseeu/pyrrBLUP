import numpy as np
import pandas as pd

from pyrrblup.amat import amat


def test_amat_returns_symmetric_dataframe() -> None:
    markers = pd.DataFrame(
        {"m1": [0, 0, 1], "m2": [0, 1, 1]},
        index=pd.Index(["s1", "s2", "s3"], name="sample"),
    )

    result = amat(markers)

    assert list(result.matrix.index) == ["s1", "s2", "s3"]
    assert np.allclose(result.matrix.values, result.matrix.values.T)
