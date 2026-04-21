from pathlib import Path


def test_baseline_script_exists() -> None:
    assert Path("scripts/build_r_baseline.R").exists()
