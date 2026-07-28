#!/usr/bin/env python3
"""Build a SHA-256 manifest for new Stage E deliverables."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


STAGE_E_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = STAGE_E_ROOT / "manifest.json"
EXCLUDED = {".pycache", "__pycache__", ".matplotlib-cache"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    files = []
    for path in sorted(STAGE_E_ROOT.rglob("*")):
        if (
            not path.is_file()
            or path == OUTPUT
            or any(part in EXCLUDED for part in path.parts)
        ):
            continue
        files.append(
            {
                "path": path.relative_to(STAGE_E_ROOT).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    payload = {
        "schema": "stage_E_new_files_manifest_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": "stage_E_reversed_candidate1_repeated_scattering",
        "runtime_caches_excluded": True,
        "file_count": len(files),
        "total_size_bytes": sum(row["size_bytes"] for row in files),
        "files": files,
    }
    OUTPUT.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
