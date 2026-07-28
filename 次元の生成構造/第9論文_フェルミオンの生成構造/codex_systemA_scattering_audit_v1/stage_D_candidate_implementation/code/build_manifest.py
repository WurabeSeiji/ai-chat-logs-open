#!/usr/bin/env python3
"""Build the Stage D manifest from new files inside the dedicated directory."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


STAGE_D_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = STAGE_D_ROOT / "manifest.json"
EXCLUDED_PARTS = {".pycache", "__pycache__", ".matplotlib-cache"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    records = []
    for path in sorted(STAGE_D_ROOT.rglob("*")):
        if (
            not path.is_file()
            or path == OUTPUT
            or any(part in EXCLUDED_PARTS for part in path.parts)
        ):
            continue
        records.append(
            {
                "path": path.relative_to(STAGE_D_ROOT).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    payload = {
        "schema": "stage_D_new_files_manifest_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": "stage_D_candidate_implementation",
        "policy": (
            "Only newly created Stage D files under the dedicated output "
            "directory are listed. Runtime caches are excluded."
        ),
        "file_count": len(records),
        "total_size_bytes": sum(row["size_bytes"] for row in records),
        "files": records,
    }
    OUTPUT.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
