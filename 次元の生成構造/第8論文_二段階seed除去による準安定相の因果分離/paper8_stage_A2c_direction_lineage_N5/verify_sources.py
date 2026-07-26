#!/usr/bin/env python3
"""A2c固定原本・Stage A0成果物・既存関数名を実験前に検証する。"""

from __future__ import annotations

import ast
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LOGS = HERE / "logs"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def defined_names(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        if isinstance(node, ast.FunctionDef) and node.name == "__init__":
            names.add("__init__")
    return names


def main() -> None:
    cfg = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
    expected = json.loads((HERE / "expected_hashes.json").read_text(encoding="utf-8"))
    if cfg["stage"] != "A2c" or cfg["n"] != 5 or cfg["dtype"] != "float64":
        raise SystemExit("SOURCE_MISMATCH: config Stage/N/dtype")
    if cfg["trajectory"]["max_step"] != 5000 or cfg["transverse"]["t0"] != 4167:
        raise SystemExit("SOURCE_MISMATCH: max_step/t0")

    checks = []
    ok = True
    for group in ("sources", "dependencies", "stage_a0"):
        for name, item in expected[group].items():
            path = REPO / item["path"]
            exists = path.is_file()
            actual = sha256(path) if exists else None
            hash_match = exists and actual == item["sha256"]
            required = item.get("required", [])
            available = defined_names(path) if exists and required else set()
            missing = sorted(set(required) - available)
            item_ok = hash_match and not missing
            checks.append({
                "group": group,
                "name": name,
                "absolute_path": str(path),
                "exists": exists,
                "expected_sha256": item["sha256"],
                "actual_sha256": actual,
                "hash_match": hash_match,
                "required_names": required,
                "missing_names": missing,
                "ok": item_ok,
            })
            ok = ok and item_ok

    report = HERE.parent / "paper7_N5_reproduction" / "reports" / "paper7_N5_reproduction_report.md"
    a0_exact = report.is_file() and "REPRODUCED_EXACTLY" in report.read_text(encoding="utf-8")
    ok = ok and a0_exact
    LOGS.mkdir(parents=True, exist_ok=True)
    result = {
        "stage": "A2c",
        "status": "VERIFIED" if ok else "SOURCE_MISMATCH",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "stage_a0_report_exact_gate": a0_exact,
        "checks": checks,
    }
    (LOGS / "source_verification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if not ok:
        raise SystemExit("SOURCE_MISMATCH: " + ", ".join(c["name"] for c in checks if not c["ok"]))
    print(f"VERIFIED: {len(checks)} files and all required existing definitions")


if __name__ == "__main__":
    main()
