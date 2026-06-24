"""Write JSON artifacts for EyeFlow runs."""

import json
from pathlib import Path


def write_json_file(path: str | Path, payload) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return target
