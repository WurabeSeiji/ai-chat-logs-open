#!/usr/bin/env python3
"""Compare Stage A before/after source manifests and fail on any source change."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


COMPARE_FIELDS = ("exists", "size_bytes", "mtime_ns", "sha256", "line_count")


def indexed(payload: dict) -> dict[str, dict]:
    return {str(record["path"]): record for record in payload["sources"]}


def main() -> int:
    audit_root = Path(__file__).resolve().parents[1]
    manifest_dir = audit_root / "manifests"
    before_path = manifest_dir / "source_manifest_before.json"
    after_path = manifest_dir / "source_manifest_after.json"
    output_path = manifest_dir / "source_manifest_comparison.json"

    before = json.loads(before_path.read_text(encoding="utf-8"))
    after = json.loads(after_path.read_text(encoding="utf-8"))
    before_by_path = indexed(before)
    after_by_path = indexed(after)
    all_paths = sorted(set(before_by_path) | set(after_by_path))
    comparisons = []
    changed_paths = []
    for path in all_paths:
        old = before_by_path.get(path)
        new = after_by_path.get(path)
        field_changes = {}
        for field in COMPARE_FIELDS:
            old_value = old.get(field) if old is not None else None
            new_value = new.get(field) if new is not None else None
            if old_value != new_value:
                field_changes[field] = {"before": old_value, "after": new_value}
        changed = old is None or new is None or bool(field_changes)
        comparisons.append(
            {
                "path": path,
                "scope_before": old.get("scope") if old else None,
                "scope_after": new.get("scope") if new else None,
                "changed": changed,
                "field_changes": field_changes,
            }
        )
        if changed:
            changed_paths.append(path)

    payload = {
        "schema": "codex_systemA_source_manifest_comparison_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_count_before": len(before_by_path),
        "source_count_after": len(after_by_path),
        "unchanged": not changed_paths,
        "changed_path_count": len(changed_paths),
        "changed_paths": changed_paths,
        "comparisons": comparisons,
    }
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(output_path)
    if changed_paths:
        raise SystemExit(f"audited source changed: {changed_paths}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
