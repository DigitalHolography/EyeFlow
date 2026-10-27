"""Write JSON sidecars owned by EyeFlow."""

import json
from pathlib import Path
from typing import Mapping


def write_json_file(path: str | Path, values: Mapping[str, object]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(values, indent=2, ensure_ascii=False), encoding="utf-8")
    return target
