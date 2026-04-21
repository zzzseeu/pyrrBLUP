# `mixed_solve` Usage

## Signature

```python
from pyrrblup.mixed_solve import mixed_solve

result = mixed_solve(y, *, Z=None, K=None, X=None)
```

This interface is intentionally close to `rrBLUP::mixed.solve`, but uses Python matrix objects directly.

## Accepted Inputs

- `y`: `numpy.ndarray` or `pandas.Series`
- `X`: `numpy.ndarray` or `pandas.DataFrame`
- `Z`: `numpy.ndarray` or `pandas.DataFrame`
- `K`: `numpy.ndarray` or `pandas.DataFrame`

The function does not align on index labels. Inputs are used in the order provided.

## Default Behavior

- `X=None`: use a single intercept column
- `Z=None`: use the identity map
- `K=None`: use the identity matrix in random-effect space

## Return Value

`mixed_solve` returns a `MixedSolveResult` with these fields:

- `beta`: fixed-effect estimates
- `u`: random-effect estimates
- `Vu`: random-effect variance estimate
- `Ve`: residual variance estimate
- `beta_terms`: names for `beta`
- `u_terms`: names for `u`

If `X` or `Z` is passed as a `DataFrame`, term names come from the columns. Otherwise default names such as `x0` and `u0` are used.

## Examples

### 1. Intercept-only model

```python
import numpy as np
from pyrrblup.mixed_solve import mixed_solve

y = np.array([1.2, 0.7, 1.5, 1.1])
result = mixed_solve(y)

print(result.beta)
print(result.Vu, result.Ve)
```

### 2. Fixed effects only

```python
import numpy as np
import pandas as pd
from pyrrblup.mixed_solve import mixed_solve

y = np.array([1.0, 2.0, 1.5, 2.5])
X = pd.DataFrame(
    {
        "Intercept": [1.0, 1.0, 1.0, 1.0],
        "covariate": [0.0, 1.0, 0.5, 1.5],
    }
)

result = mixed_solve(y, X=X)
print(result.beta_terms)
print(result.beta)
```

### 3. General mixed model

```python
import numpy as np
import pandas as pd
from pyrrblup.mixed_solve import mixed_solve

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

result = mixed_solve(y, X=X, Z=Z, K=K)
print(result.beta_terms, result.beta)
print(result.u_terms, result.u)
print(result.Vu, result.Ve)
```

## Error Rules

These cases raise `ValueError`:

- `X` and `y` row counts differ
- `Z` and `y` row counts differ
- `K` shape does not match the number of `Z` columns
- `X` or `Z` is not two-dimensional after normalization
