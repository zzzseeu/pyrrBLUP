from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pyrrblup import kin_blup


def main() -> int:
    here = Path(__file__).resolve().parent
    result = kin_blup(here / "genotype_012.csv", here / "phenotype.csv", "loc_BeiJ")
    print("samples", len(result.predictions))
    print("Vu", result.mixed_solve.Vu)
    print("Ve", result.mixed_solve.Ve)
    print(result.predictions.head().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
