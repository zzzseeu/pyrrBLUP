from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from pyrrblup.models import MixedSolveResult


LOWER_LOG_LAMBDA = -8.0
INITIAL_UPPER_LOG_LAMBDA = 8.0
MAX_UPPER_LOG_LAMBDA = 32.0
BOUNDARY_TOLERANCE = 1e-3


def _normalize_design_matrix(
    matrix: np.ndarray | pd.DataFrame | None,
    *,
    rows: int,
    default: np.ndarray | None = None,
    default_terms: list[str] | None = None,
    prefix: str,
    name: str,
) -> tuple[np.ndarray, list[str]]:
    if matrix is None:
        if default is None or default_terms is None:
            raise ValueError(f"{name} is required")
        return default.astype(float), default_terms

    if isinstance(matrix, pd.DataFrame):
        values = matrix.to_numpy(dtype=float)
        terms = [str(column) for column in matrix.columns]
    else:
        values = np.asarray(matrix, dtype=float)
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        terms = [f"{prefix}{idx}" for idx in range(values.shape[1])]

    if values.ndim != 2:
        raise ValueError(f"{name} must be a 2D matrix")
    if values.shape[0] != rows:
        raise ValueError(f"{name} must have the same number of rows as y")
    return values, terms


def _normalize_covariance(K: np.ndarray | pd.DataFrame | None, cols: int) -> np.ndarray:
    if K is None:
        return np.eye(cols, dtype=float)
    if isinstance(K, pd.DataFrame):
        values = K.to_numpy(dtype=float)
    else:
        values = np.asarray(K, dtype=float)
    if values.shape != (cols, cols):
        raise ValueError(f"K must have shape ({cols}, {cols}) to match Z")
    return values


def _neg_reml(log_lambda: float, y: np.ndarray, x: np.ndarray, zkzt: np.ndarray) -> float:
    lam = float(np.exp(log_lambda))
    v = zkzt + lam * np.eye(zkzt.shape[0])
    sign_v, logdet_v = np.linalg.slogdet(v)
    if sign_v <= 0:
        return np.inf

    v_inv = np.linalg.inv(v)
    xt_vinv = x.T @ v_inv
    xt_vinv_x = xt_vinv @ x
    sign_x, logdet_x = np.linalg.slogdet(xt_vinv_x)
    if sign_x <= 0:
        return np.inf

    beta = np.linalg.solve(xt_vinv_x, xt_vinv @ y)
    resid = y - x @ beta
    quad = float(resid.T @ v_inv @ resid)
    if quad <= 0.0:
        return -np.inf
    n = y.shape[0]
    p = x.shape[1]
    return 0.5 * (logdet_v + logdet_x + (n - p) * np.log(quad))


def _optimize_log_lambda(y: np.ndarray, x: np.ndarray, zkzt: np.ndarray) -> float:
    upper = INITIAL_UPPER_LOG_LAMBDA
    while True:
        opt = minimize_scalar(_neg_reml, bounds=(LOWER_LOG_LAMBDA, upper), method="bounded", args=(y, x, zkzt))
        if not opt.success:
            raise RuntimeError(f"REML optimization failed: {opt.message}")
        if upper >= MAX_UPPER_LOG_LAMBDA:
            return float(opt.x)
        if opt.x < upper - BOUNDARY_TOLERANCE:
            return float(opt.x)
        upper = min(upper * 2.0, MAX_UPPER_LOG_LAMBDA)


def mixed_solve(
    y: np.ndarray | pd.Series,
    *,
    Z: np.ndarray | pd.DataFrame | None = None,
    K: np.ndarray | pd.DataFrame | None = None,
    X: np.ndarray | pd.DataFrame | None = None,
) -> MixedSolveResult:
    """Solve the single-random-effect mixed model using an rrBLUP-like matrix interface.

    Signature
    ---------
    mixed_solve(y, *, Z=None, K=None, X=None)

    Parameters
    ----------
    y
        Response vector with length ``n``. Accepts ``numpy.ndarray`` or ``pandas.Series``.
    Z
        Random-effect design matrix with shape ``(n, q)``. Accepts ``numpy.ndarray`` or
        ``pandas.DataFrame``. Defaults to the identity map.
    K
        Random-effect covariance matrix with shape ``(q, q)``. Accepts ``numpy.ndarray``
        or ``pandas.DataFrame``. Defaults to the identity matrix.
    X
        Fixed-effect design matrix with shape ``(n, p)``. Accepts ``numpy.ndarray`` or
        ``pandas.DataFrame``. Defaults to an intercept column.

    Notes
    -----
    - This function does not align inputs by index or column label.
    - Inputs are consumed in the order they are provided.
    - Dimension mismatches raise ``ValueError``.

    Returns
    -------
    MixedSolveResult
        Result object with fields ``beta``, ``u``, ``Vu``, ``Ve``, ``beta_terms``, and
        ``u_terms``.

    Examples
    --------
    Intercept-only model::

        result = mixed_solve(y)

    Fixed effects with an explicit design matrix::

        X = pd.DataFrame({"Intercept": 1.0, "covariate": [0.2, 0.4, 0.6]})
        result = mixed_solve(y, X=X)

    General mixed model::

        result = mixed_solve(y, X=X, Z=Z, K=K)
    """
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.shape[0]

    x_values, beta_terms = _normalize_design_matrix(
        X,
        rows=n,
        default=np.ones((n, 1), dtype=float),
        default_terms=["Intercept"],
        prefix="x",
        name="X",
    )
    z_values, u_terms = _normalize_design_matrix(
        Z,
        rows=n,
        default=np.eye(n, dtype=float),
        default_terms=[f"u{idx}" for idx in range(n)],
        prefix="u",
        name="Z",
    )
    k_values = _normalize_covariance(K, z_values.shape[1])
    zkzt = z_values @ k_values @ z_values.T

    log_lambda = _optimize_log_lambda(y, x_values, zkzt)
    lam = float(np.exp(log_lambda))
    v = zkzt + lam * np.eye(n)
    v_inv = np.linalg.inv(v)
    beta = np.linalg.solve(x_values.T @ v_inv @ x_values, x_values.T @ v_inv @ y)
    resid = y - x_values @ beta
    Vu = float((resid.T @ v_inv @ resid) / max(n - x_values.shape[1], 1))
    Ve = float(Vu * lam)
    u = k_values @ z_values.T @ v_inv @ resid
    return MixedSolveResult(beta=beta, u=u, Vu=Vu, Ve=Ve, beta_terms=beta_terms, u_terms=u_terms)
