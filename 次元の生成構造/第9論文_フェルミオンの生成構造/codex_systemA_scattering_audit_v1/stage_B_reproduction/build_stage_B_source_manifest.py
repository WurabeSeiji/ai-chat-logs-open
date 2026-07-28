#!/usr/bin/env python3
"""Build and compare Stage B hashes without importing or executing originals."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


SEMANTIC_NOTICE = (
    "B_to_A_transfer is spectral cosine similarity of the A-channel state "
    "to the initial B spectrum; it is NOT a path-exchange norm."
)

AUDIT_ROOT = Path(__file__).resolve().parents[1]
TARGET_ROOT = AUDIT_ROOT.parent
MANIFEST_ROOT = AUDIT_ROOT / "manifests"

SOURCES = (
    TARGET_ROOT
    / "対照実験_波束収縮_実行環境_v1"
    / "20260715"
    / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    TARGET_ROOT
    / "対照実験_波束収縮_実行環境_v1"
    / "20260713"
    / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py",
    TARGET_ROOT
    / "対照実験_波束収縮_実行環境_v1"
    / "20260715"
    / "run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py",
)

COMPARE_FIELDS = ("exists", "size_bytes", "mtime_ns", "line_count", "sha256")


def record(path: Path) -> dict:
    relative = path.relative_to(TARGET_ROOT).as_posix()
    if not path.is_file():
        return {"path": relative, "exists": False}
    data = path.read_bytes()
    stat = path.stat()
    return {
        "path": relative,
        "exists": True,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "line_count": len(data.splitlines()),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def build(phase: str) -> Path:
    MANIFEST_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "system_A_stage_B_source_manifest_v1",
        "phase": phase,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": (
            "Audited System A/System B originals are read and hashed only; "
            "they are not imported, executed, modified, or used as output targets."
        ),
        "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
        "source_count": len(SOURCES),
        "sources": [record(path) for path in SOURCES],
    }
    output = MANIFEST_ROOT / f"stage_B_source_manifest_{phase}.json"
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return output


def compare() -> Path:
    before_path = MANIFEST_ROOT / "stage_B_source_manifest_before.json"
    after_path = MANIFEST_ROOT / "stage_B_source_manifest_after.json"
    before = json.loads(before_path.read_text(encoding="utf-8"))
    after = json.loads(after_path.read_text(encoding="utf-8"))
    old_by_path = {row["path"]: row for row in before["sources"]}
    new_by_path = {row["path"]: row for row in after["sources"]}
    paths = sorted(set(old_by_path) | set(new_by_path))
    rows = []
    changed_paths = []
    for path in paths:
        old = old_by_path.get(path)
        new = new_by_path.get(path)
        changes = {}
        for field in COMPARE_FIELDS:
            old_value = old.get(field) if old else None
            new_value = new.get(field) if new else None
            if old_value != new_value:
                changes[field] = {"before": old_value, "after": new_value}
        changed = old is None or new is None or bool(changes)
        rows.append({"path": path, "changed": changed, "field_changes": changes})
        if changed:
            changed_paths.append(path)
    payload = {
        "schema": "system_A_stage_B_source_manifest_comparison_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
        "source_count_before": len(old_by_path),
        "source_count_after": len(new_by_path),
        "unchanged": not changed_paths,
        "changed_path_count": len(changed_paths),
        "changed_paths": changed_paths,
        "comparisons": rows,
    }
    output = MANIFEST_ROOT / "stage_B_source_manifest_comparison.json"
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if changed_paths:
        raise SystemExit(f"Stage B source integrity failure: {changed_paths}")
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("before", "after", "compare"))
    args = parser.parse_args()
    output = compare() if args.action == "compare" else build(args.action)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
