#!/usr/bin/env python3
"""Stage A0 source and baseline hash gate. This script performs no experiment."""

from __future__ import annotations

import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
EXPECTED_PATH = PACKAGE_ROOT / "expected_hashes.json"
CONFIG_PATH = PACKAGE_ROOT / "config_locked.json"
RESULT_PATH = PACKAGE_ROOT / "comparison" / "source_verification.json"
LOG_PATH = PACKAGE_ROOT / "logs" / "verify_sources.log"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_path(relative_path: str, reproduced: bool = False) -> Path:
    if reproduced:
        return PACKAGE_ROOT / relative_path
    return REPO_ROOT / relative_path


def check_entry(group: str, name: str, item: dict) -> dict:
    path_key = "expected" if group.endswith("_baselines") else "path"
    path = resolve_path(item[path_key])
    result = {
        "group": group,
        "name": name,
        "path": str(path),
        "expected_sha256": item["sha256"],
        "exists": path.is_file(),
    }
    if result["exists"]:
        stat = path.stat()
        actual = sha256(path)
        result.update(
            {
                "actual_sha256": actual,
                "sha256_match": actual == item["sha256"],
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
    else:
        result.update(
            {
                "actual_sha256": None,
                "sha256_match": False,
                "size_bytes": None,
                "mtime_ns": None,
            }
        )
    return result


def main() -> int:
    started = time.perf_counter()
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    if config.get("locked_n") != 5 or config.get("allowed_n") != [5]:
        raise RuntimeError("config_locked.json が N=5 だけに固定されていない")

    checks = []
    checks.append(check_entry("locked_spec", "Stage_A_実行前固定仕様書.md", expected["locked_spec"]))
    for group in ("sources", "dependencies", "csv_baselines", "json_baselines", "png_baselines"):
        for name, item in expected[group].items():
            checks.append(check_entry(group, name, item))

    failed = [item for item in checks if not item["sha256_match"]]
    result = {
        "stage": "A0",
        "locked_n": 5,
        "success": not failed,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "repo_root": str(REPO_ROOT),
        "package_root": str(PACKAGE_ROOT),
        "checks": checks,
        "failures": failed,
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = []
    for item in checks:
        status = "OK" if item["sha256_match"] else "NG"
        lines.append(f"[{status}] {item['group']} {item['name']} {item['actual_sha256']}")
    lines.append("SUCCESS" if result["success"] else "STOP: SHA-256不一致または必須ファイル欠落")
    LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if result["success"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
