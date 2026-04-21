import numpy as np
import pandas as pd
import pytest

from pyrrblup.mixed_solve import mixed_solve


def test_mixed_solve_returns_expected_shapes() -> None:
    y = np.array([1.0, 2.0, 3.0])
    k = np.eye(3)

    result = mixed_solve(y, K=k)

    assert result.beta.shape == (1,)
    assert result.u.shape == (3,)
    assert result.beta_terms == ["Intercept"]
    assert result.u_terms == ["u0", "u1", "u2"]
    assert result.Vu >= 0.0
    assert result.Ve >= 0.0


def test_mixed_solve_preserves_dataframe_term_names() -> None:
    y = np.array([1.0, 2.0, 1.5, 2.5])
    x = pd.DataFrame(
        {
            "Intercept": [1.0, 1.0, 1.0, 1.0],
            "covariate": [0.0, 1.0, 0.5, 1.5],
        }
    )
    z = pd.DataFrame(
        {
            "line_A": [1.0, 0.0, 1.0, 0.0],
            "line_B": [0.0, 1.0, 0.0, 1.0],
        }
    )
    k = np.eye(2)

    result = mixed_solve(y, X=x, Z=z, K=k)

    assert result.beta.shape == (2,)
    assert result.u.shape == (2,)
    assert result.beta_terms == ["Intercept", "covariate"]
    assert result.u_terms == ["line_A", "line_B"]


def test_mixed_solve_matches_between_dataframe_and_ndarray_inputs() -> None:
    y = np.array([1.0, 2.0, 1.5, 2.5])
    x_df = pd.DataFrame(
        {
            "Intercept": [1.0, 1.0, 1.0, 1.0],
            "covariate": [0.0, 1.0, 0.5, 1.5],
        }
    )
    z_df = pd.DataFrame(
        {
            "line_A": [1.0, 0.0, 1.0, 0.0],
            "line_B": [0.0, 1.0, 0.0, 1.0],
        }
    )
    k = np.array([[1.0, 0.2], [0.2, 1.0]])

    frame_result = mixed_solve(y, X=x_df, Z=z_df, K=k)
    array_result = mixed_solve(y, X=x_df.to_numpy(), Z=z_df.to_numpy(), K=k)

    assert np.allclose(frame_result.beta, array_result.beta)
    assert np.allclose(frame_result.u, array_result.u)
    assert frame_result.beta_terms == ["Intercept", "covariate"]
    assert array_result.beta_terms == ["x0", "x1"]
    assert frame_result.u_terms == ["line_A", "line_B"]
    assert array_result.u_terms == ["u0", "u1"]


def test_mixed_solve_defaults_to_identity_k_when_missing() -> None:
    y = np.array([1.0, 2.0, 1.5, 2.5])
    x = pd.DataFrame({"Intercept": [1.0, 1.0, 1.0, 1.0]})
    z = pd.DataFrame(
        {
            "line_A": [1.0, 0.0, 1.0, 0.0],
            "line_B": [0.0, 1.0, 0.0, 1.0],
        }
    )

    implicit_k = mixed_solve(y, X=x, Z=z, K=None)
    explicit_k = mixed_solve(y, X=x, Z=z, K=np.eye(2))

    assert np.allclose(implicit_k.beta, explicit_k.beta)
    assert np.allclose(implicit_k.u, explicit_k.u)


def test_mixed_solve_rejects_mismatched_x_rows() -> None:
    y = np.array([1.0, 2.0, 3.0])
    x = np.ones((2, 1))

    with pytest.raises(ValueError, match="same number of rows"):
        mixed_solve(y, X=x, K=np.eye(3))


def test_mixed_solve_rejects_mismatched_z_and_k_dimensions() -> None:
    y = np.array([1.0, 2.0, 1.5, 2.5])
    x = pd.DataFrame({"Intercept": [1.0, 1.0, 1.0, 1.0]})
    z = pd.DataFrame(
        {
            "line_A": [1.0, 0.0, 1.0, 0.0],
            "line_B": [0.0, 1.0, 0.0, 1.0],
        }
    )

    with pytest.raises(ValueError, match="K must have shape"):
        mixed_solve(y, X=x, Z=z, K=np.eye(3))
