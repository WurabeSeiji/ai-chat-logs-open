#!/usr/bin/env python3
"""Build a read-only source manifest for the System A Stage A audit.

This script never imports or executes an audited source file.  It only reads
bytes and stat metadata, then writes one manifest below the dedicated audit
directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


TARGET_ROOT_REL = Path("次元の生成構造/第9論文_フェルミオンの生成構造")
AUDIT_ROOT_REL = TARGET_ROOT_REL / "codex_systemA_scattering_audit_v1"

SOURCES = (
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "20260715"
        / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    ),
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "20260713"
        / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py",
    ),
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "20260715"
        / "run_system_A_localization_exchange_R_sweep_instrumented_v1.py",
    ),
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "20260715"
        / "run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py",
    ),
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "parity_suite_v1"
        / "run_parity_suite_v1.py",
    ),
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "production_dump_v1"
        / "run_production_dump_v1.py",
    ),
    (
        "inside_target",
        TARGET_ROOT_REL
        / "対照実験_波束収縮_実行環境_v1"
        / "additional_fullkernel_wavelength_v1"
        / "run_fullkernel_wavelength_experiment_v1.py",
    ),
    (
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    ),
    (
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "run_minimal_system_B_gray_direct_check_v5.py",
    ),
    (
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "phase5_eigenphase_resonance_v2.py",
    ),
    (
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "run_two_physical_roots_multiprecision_v1.py",
    ),
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def source_record(repo_root: Path, scope: str, relative_path: Path) -> dict:
    path = repo_root / relative_path
    if not path.is_file():
        return {
            "path": relative_path.as_posix(),
            "scope": scope,
            "exists": False,
        }
    data = path.read_bytes()
    stat = path.stat()
    return {
        "path": relative_path.as_posix(),
        "scope": scope,
        "exists": True,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_bytes(data),
        "line_count": len(data.splitlines()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--phase", choices=("before", "after"), required=True)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    audit_root = repo_root / AUDIT_ROOT_REL
    manifest_dir = audit_root / "manifests"
    log_dir = audit_root / "logs"
    report_dir = audit_root / "reports"
    stage_dir = audit_root / "stage_A_audit"
    for directory in (manifest_dir, log_dir, report_dir, stage_dir):
        directory.mkdir(parents=True, exist_ok=True)

    records = [
        source_record(repo_root, scope, relative_path)
        for scope, relative_path in SOURCES
    ]
    payload = {
        "schema": "codex_systemA_scattering_source_manifest_v1",
        "phase": args.phase,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": {
            "inside_target": "read-only audited originals; no overwrite",
            "outside_readonly": "read and SHA-256 only; no execution, import, modification, or output",
        },
        "source_count": len(records),
        "sources": records,
    }
    output_path = manifest_dir / f"source_manifest_{args.phase}.json"
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
