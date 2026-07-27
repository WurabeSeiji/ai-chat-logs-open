#!/usr/bin/env python3
"""Stage A2a固定原本・Stage A0/A1b入力のSHA-256検証。実験は行わない。"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
HASHES = HERE / "expected_hashes.json"
CONFIG = HERE / "config_locked.json"
LOGS = HERE / "logs"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    expected = json.loads(HASHES.read_text(encoding="utf-8"))
    if cfg["stage"] != "A2a" or cfg["n"] != 5 or cfg["dtype"] != "float64":
        raise SystemExit("SOURCE_MISMATCH: config_locked.jsonのStage/N/dtypeが固定仕様と不一致")
    if cfg["run_ids"] != ["A2a_N5_seedless_f64_e1", "A2a_N5_seedless_f64_e2"]:
        raise SystemExit("SOURCE_MISMATCH: run_id不一致")

    checks = []
    ok = True
    for group in ("sources", "dependencies", "stage_a0_inputs", "stage_a1b_inputs"):
        for name, item in expected[group].items():
            path = REPO / item["path"]
            exists = path.is_file()
            actual = sha256(path) if exists else None
            match = exists and actual == item["sha256"]
            checks.append({
                "group": group,
                "name": name,
                "path": str(path),
                "exists": exists,
                "expected_sha256": item["sha256"],
                "actual_sha256": actual,
                "match": match,
            })
            ok = ok and match

    a0_report = HERE.parent / "paper7_N5_reproduction" / "reports" / "paper7_N5_reproduction_report.md"
    a1b_report = HERE.parent / "paper7_N5_transition_anatomy" / "reports" / "paper7_N5_transition_anatomy_report.md"
    report_checks = {
        "stage_a0_report_exists": a0_report.is_file(),
        "stage_a0_reproduced_exactly": a0_report.is_file() and "REPRODUCED_EXACTLY" in a0_report.read_text(encoding="utf-8"),
        "stage_a1b_report_exists": a1b_report.is_file(),
        "stage_a1b_complete": a1b_report.is_file() and "TRANSITION_ANATOMY_COMPLETE" in a1b_report.read_text(encoding="utf-8"),
    }
    ok = ok and all(report_checks.values())

    LOGS.mkdir(parents=True, exist_ok=True)
    result = {
        "stage": "A2a",
        "status": "VERIFIED" if ok else "SOURCE_MISMATCH",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "checks": checks,
        "report_checks": report_checks,
    }
    (LOGS / "source_verification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if not ok:
        failed = [c["path"] for c in checks if not c["match"]]
        raise SystemExit("SOURCE_MISMATCH: " + ", ".join(failed or ["Stage A0/A1b report gate"]))
    print(f"VERIFIED: {len(checks)} files; N=5 float64 only")


if __name__ == "__main__":
    main()
