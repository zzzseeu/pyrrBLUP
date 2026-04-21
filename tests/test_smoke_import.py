from pyrrblup import __version__, amat, kin_blup, mixed_solve


def test_package_imports() -> None:
    assert __version__
    assert callable(amat)
    assert callable(mixed_solve)
    assert callable(kin_blup)
