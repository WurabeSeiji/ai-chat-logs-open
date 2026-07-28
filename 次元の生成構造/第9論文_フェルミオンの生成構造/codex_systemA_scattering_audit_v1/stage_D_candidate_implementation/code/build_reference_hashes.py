#!/usr/bin/env python3
"""Hash Stage D references without importing or executing them."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


STAGE_D_ROOT = Path(__file__).resolve().parents[1]
AUDIT_ROOT = STAGE_D_ROOT.parent
TARGET_ROOT = AUDIT_ROOT.parent
DATA_ROOT = STAGE_D_ROOT / "data"

REFERENCES = (
    TARGET_ROOT / "Codex向け_StageD_基準量取得とCandidate013独立実装比較指示.md",
    AUDIT_ROOT / "reports" / "Stage_C_mathematical_design_report.md",
    AUDIT_ROOT / "reports" / "06_existing_behavior_reproduction.md",
    AUDIT_ROOT / "data" / "stage_B" / "current_behavior_baseline.csv",
    AUDIT_ROOT / "data" / "stage_B" / "current_behavior_diagnostics.json",
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


def file_record(path: Path) -> dict:
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
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "stage_D_reference_hashes_v1",
        "phase": phase,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": (
            "References are read and SHA-256 hashed only. They are not imported, "
            "executed, modified, or used as output targets."
        ),
        "reference_count": len(REFERENCES),
        "references": [file_record(path) for path in REFERENCES],
    }
    output = DATA_ROOT / f"reference_hashes_{phase}.json"
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return output


def compare() -> Path:
    before_path = DATA_ROOT / "reference_hashes_before.json"
    after_path = DATA_ROOT / "reference_hashes_after.json"
    before = json.loads(before_path.read_text(encoding="utf-8"))
    after = json.loads(after_path.read_text(encoding="utf-8"))
    old_by_path = {row["path"]: row for row in before["references"]}
    new_by_path = {row["path"]: row for row in after["references"]}
    paths = sorted(set(old_by_path) | set(new_by_path))
    comparisons = []
    changed_paths = []
    for path in paths:
        old = old_by_path.get(path)
        new = new_by_path.get(path)
        field_changes = {}
        for field in COMPARE_FIELDS:
            old_value = old.get(field) if old else None
            new_value = new.get(field) if new else None
            if old_value != new_value:
                field_changes[field] = {"before": old_value, "after": new_value}
        changed = old is None or new is None or bool(field_changes)
        comparisons.append(
            {"path": path, "changed": changed, "field_changes": field_changes}
        )
        if changed:
            changed_paths.append(path)
    payload = {
        "schema": "stage_D_reference_hash_comparison_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference_count_before": len(old_by_path),
        "reference_count_after": len(new_by_path),
        "unchanged": not changed_paths,
        "changed_path_count": len(changed_paths),
        "changed_paths": changed_paths,
        "comparisons": comparisons,
    }
    output = DATA_ROOT / "reference_hash_comparison.json"
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if changed_paths:
        raise SystemExit(f"Stage D reference changed: {changed_paths}")
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
