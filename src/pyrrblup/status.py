from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from pyrrblup.models import BatchRunStatus


def write_status_json(status: BatchRunStatus, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(asdict(status), indent=2), encoding="utf-8")
