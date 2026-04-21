# pyrrBLUP

Python implementation of core `rrBLUP` workflows with numeric checks against the R package.

## Installation

```bash
pip install pyrrblup
```

## Quick Start

```python
import numpy as np
import pandas as pd

from pyrrblup import amat, kin_blup, mixed_solve

markers = pd.DataFrame(
    {
        "m1": [0, 0, 2, 2],
        "m2": [0, 2, 0, 2],
        "m3": [0, 1, 1, 2],
    },
    index=["s1", "s2", "s3", "s4"],
)

kinship = amat(markers).matrix
y = np.array([1.0, 2.0, 1.5, 2.5])
X = pd.DataFrame(
    {
        "Intercept": [1.0, 1.0, 1.0, 1.0],
        "covariate": [0.0, 1.0, 0.5, 1.5],
    }
)
Z = pd.DataFrame(
    {
        "line_A": [1.0, 0.0, 1.0, 0.0],
        "line_B": [0.0, 1.0, 0.0, 1.0],
    }
)
K = np.array(
    [
        [1.0, 0.2],
        [0.2, 1.0],
    ]
)

print(kinship.shape)
fit = mixed_solve(y, X=X, Z=Z, K=K)
print(fit.beta_terms, fit.beta)
print(fit.Vu, fit.Ve)
```

Run `kin_blup` from file paths:

```python
from pyrrblup import kin_blup

result = kin_blup("genotype.csv", "phenotype.csv", "loc_BeiJ")
print(result.predictions.head())
```

`mixed_solve` follows an rrBLUP-like matrix interface:

```python
result = mixed_solve(y, Z=None, K=None, X=None)
```

`kin_blup` accepts either:

- genotype and phenotype file paths
- genotype and phenotype `pandas.DataFrame` objects

Both inputs must be the same kind. Mixed input types are rejected.

## Real Data Examples

This repository also includes deterministic subsets cut from the real `loc_BeiJ` data.

- `examples/minimal_real_data/`
  - `32` samples
  - `64` markers
  - intended for quick README-scale runs

- `examples/medium_real_data/`
  - `256` samples
  - `512` markers
  - intended for more realistic local demonstrations

Build or refresh the example data:

```bash
python scripts/build_examples.py
```

Run the minimal real-data example:

```bash
python examples/minimal_real_data/run_example.py
```

More detail:

- `docs/mixed_solve.md`
- `examples/README.md`

## Development

Install from source:

```bash
pip install -e .
```

Run the test suite with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest -v
```
