#!/usr/bin/env python3
"""Hash Stage E references without importing, executing, or modifying them."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


STAGE_E_ROOT = Path(__file__).resolve().parents[1]
AUDIT_ROOT = STAGE_E_ROOT.parent
PAPER_ROOT = AUDIT_ROOT.parent
DATA_ROOT = STAGE_E_ROOT / "data"

REFERENCES = (
    PAPER_ROOT
    / "対照実験_波束収縮_実行環境_v1"
    / "20260715"
    / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    PAPER_ROOT
    / "対照実験_波束収縮_実行環境_v1"
    / "20260713"
    / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py",
    AUDIT_ROOT
    / "stage_B_reproduction"
    / "reproduce_current_system_A_stage_B.py",
    AUDIT_ROOT / "data" / "stage_B" / "current_behavior_baseline.csv",
    AUDIT_ROOT
    / "stage_D_candidate_implementation"
    / "code"
    / "scattering_api.py",
    AUDIT_ROOT
    / "stage_D_candidate_implementation"
    / "code"
    / "parity_metrics.py",
    AUDIT_ROOT
    / "stage_D_candidate_implementation"
    / "data"
    / "stage_D_summary.json",
    AUDIT_ROOT
    / "stage_D_candidate_implementation"
    / "reports"
    / "Stage_D_report.md",
)

COMPARE_FIELDS = ("exists", "size_bytes", "mtime_ns", "sha256")


def _record(path: Path) -> dict:
    relative = path.relative_to(PAPER_ROOT).as_posix()
    if not path.is_file():
        return {"path": relative, "exists": False}
    data = path.read_bytes()
    stat = path.stat()
    return {
        "path": relative,
        "exists": True,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def _build(phase: str) -> Path:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "stage_E_reference_hashes_v1",
        "phase": phase,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": (
            "References are read and SHA-256 hashed only; they are not "
            "imported, executed, modified, or used as output targets."
        ),
        "reference_count": len(REFERENCES),
        "references": [_record(path) for path in REFERENCES],
    }
    output = DATA_ROOT / f"reference_hashes_{phase}.json"
    output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return output


def _compare() -> Path:
    before = json.loads(
        (DATA_ROOT / "reference_hashes_before.json").read_text(
            encoding="utf-8"
        )
    )
    after = json.loads(
        (DATA_ROOT / "reference_hashes_after.json").read_text(
            encoding="utf-8"
        )
    )
    old = {row["path"]: row for row in before["references"]}
    new = {row["path"]: row for row in after["references"]}
    comparisons = []
    changed_paths = []
    for path in sorted(set(old) | set(new)):
        changes = {}
        for field in COMPARE_FIELDS:
            before_value = old.get(path, {}).get(field)
            after_value = new.get(path, {}).get(field)
            if before_value != after_value:
                changes[field] = {
                    "before": before_value,
                    "after": after_value,
                }
        changed = bool(changes)
        comparisons.append(
            {"path": path, "changed": changed, "field_changes": changes}
        )
        if changed:
            changed_paths.append(path)
    payload = {
        "schema": "stage_E_reference_hash_comparison_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference_count_before": len(old),
        "reference_count_after": len(new),
        "unchanged": not changed_paths,
        "changed_path_count": len(changed_paths),
        "changed_paths": changed_paths,
        "comparisons": comparisons,
    }
    output = DATA_ROOT / "reference_hash_comparison.json"
    output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if changed_paths:
        raise SystemExit(f"Stage E reference changed: {changed_paths}")
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("before", "after", "compare"))
    args = parser.parse_args()
    output = _compare() if args.action == "compare" else _build(args.action)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
