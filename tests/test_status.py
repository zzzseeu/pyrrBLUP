from pathlib import Path

from pyrrblup.models import BatchRunStatus
from pyrrblup.status import write_status_json


def test_write_status_json_creates_expected_fields(tmp_path: Path) -> None:
    status = BatchRunStatus(
        location="loc_A",
        tool="python",
        success=True,
        error_message="",
        started_at="2026-04-16T10:00:00",
        finished_at="2026-04-16T10:00:03",
    )

    out_path = tmp_path / "status.json"
    write_status_json(status, out_path)

    payload = out_path.read_text(encoding="utf-8")
    assert '"location": "loc_A"' in payload
    assert '"tool": "python"' in payload
    assert '"success": true' in payload
