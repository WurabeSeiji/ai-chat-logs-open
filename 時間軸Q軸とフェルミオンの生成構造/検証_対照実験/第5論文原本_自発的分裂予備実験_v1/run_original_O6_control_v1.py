#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""前シリーズ第5論文・実験O6を隔離して検証・再現する。

使い方:
    python3 run_original_O6_control_v1.py validate
    python3 run_original_O6_control_v1.py audit
    python3 run_original_O6_control_v1.py run 40
    python3 run_original_O6_control_v1.py run 300
    python3 run_original_O6_control_v1.py reference 1000
    python3 run_original_O6_control_v1.py plot
    python3 run_original_O6_control_v1.py all

原本三ファイルは変更せず、走行器と作図器の出力先だけを
O6_原本対照結果_v1/ へ切り替える。
"""

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import sys
from pathlib import Path

import matplotlib
import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
SOURCE_DIR = REPO_ROOT / "次元の生成構造" / "自発的分裂予備実験_v1"
REFERENCE_DIR = SOURCE_DIR / "largeN_splitting_result_v1"
RESULT_DIR = HERE / "O6_原本対照結果_v1"

CORE_HASHES = {
    "run_n_scaling_lowrank_v1.py":
        "ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d",
    "run_spontaneous_splitting_largeN_v1.py":
        "13baf6f5158c53ee92d1a08a0a5b60832424a83222b0601563a2676a392515ac",
    "make_largeN_figure_v1.py":
        "f4cbc4efe2af2ffbffa45eec8c71ecb0340a1a3c940b78a2d83847fb03ea5804",
}

N_VALUES = (40, 300, 1000)
DELTA = 1e-15
SEED = 0
CAP = 9000
AFTER = 1500
TOL = 1e-12
TIMING_KEYS = {"t_parent_sec", "t_run_sec"}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tag_for(n):
    return f"N{n:05d}_delta{DELTA:.0e}_seed{SEED}"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"モジュールを読み込めません: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_core_hashes():
    rows = []
    failed = []
    for filename, expected in CORE_HASHES.items():
        local_path = HERE / filename
        source_path = SOURCE_DIR / filename
        local_hash = sha256(local_path)
        source_hash = sha256(source_path)
        ok = local_hash == expected and source_hash == expected
        rows.append({
            "file": filename,
            "expected_sha256": expected,
            "local_sha256": local_hash,
            "source_sha256": source_hash,
            "ok": ok,
        })
        if not ok:
            failed.append(filename)
        print(f"[{'OK' if ok else 'NG'}] {filename} {local_hash}")
    if failed:
        raise RuntimeError("原本ハッシュ不一致: " + ", ".join(failed))
    return rows


def environment():
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
    }


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(value, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def validate_dense():
    hashes = verify_core_hashes()
    sys.path.insert(0, str(HERE))
    engine = load_module("o6_engine_control", HERE / "run_n_scaling_lowrank_v1.py")
    result = engine.validate_against_dense(12, 0, steps=300)
    limits = {
        "err_matvec": 1e-12,
        "err_G": 1e-12,
        "err_sigma_spectrum": 1e-12,
        "max_traj_dev": 1e-11,
    }
    passed = all(result[key] <= limit for key, limit in limits.items())
    report = {
        "test": "low-rank engine versus dense matrix",
        "n": 12,
        "seed": 0,
        "core_hashes": hashes,
        "environment": environment(),
        "limits": limits,
        "result": result,
        "passed": passed,
    }
    write_json(RESULT_DIR / "validation_dense_N00012_seed0.json", report)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not passed:
        raise RuntimeError("低ランク実装と密行列実装の一致検証に失敗しました")
    return report


def compare_one(n, result_dir=RESULT_DIR, report_dir=None):
    if report_dir is None:
        report_dir = result_dir
    tag = tag_for(n)
    csv_name = f"fcurve_{tag}.csv"
    json_name = f"summary_{tag}.json"
    actual_csv = result_dir / csv_name
    actual_json = result_dir / json_name
    reference_csv = REFERENCE_DIR / csv_name
    reference_json = REFERENCE_DIR / json_name

    if not actual_csv.exists() or not actual_json.exists():
        raise FileNotFoundError(f"N={n} の対照出力がありません: {result_dir}")

    actual = json.loads(actual_json.read_text(encoding="utf-8"))
    reference = json.loads(reference_json.read_text(encoding="utf-8"))
    stable_actual = {key: value for key, value in actual.items()
                     if key not in TIMING_KEYS}
    stable_reference = {key: value for key, value in reference.items()
                        if key not in TIMING_KEYS}
    csv_actual_hash = sha256(actual_csv)
    csv_reference_hash = sha256(reference_csv)
    report = {
        "n": n,
        "tag": tag,
        "csv_actual_sha256": csv_actual_hash,
        "csv_reference_sha256": csv_reference_hash,
        "csv_byte_identical": csv_actual_hash == csv_reference_hash,
        "summary_stable_fields_identical": stable_actual == stable_reference,
        "actual_summary": actual,
        "reference_summary": reference,
    }
    report["passed"] = (
        report["csv_byte_identical"]
        and report["summary_stable_fields_identical"]
    )
    write_json(report_dir / f"comparison_{tag}.json", report)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return report


def audit_existing_copy():
    verify_core_hashes()
    existing = HERE / "largeN_splitting_result_v1"
    reports = [
        compare_one(n, result_dir=existing, report_dir=RESULT_DIR)
        for n in N_VALUES
    ]
    audit = {
        "target": str(existing),
        "reports": reports,
        "passed": all(report["passed"] for report in reports),
    }
    write_json(RESULT_DIR / "audit_existing_copy.json", audit)
    if not audit["passed"]:
        print("[注意] 既存コピーの保存出力は原本と完全一致していません。")
    return audit


def run_one(n):
    if n not in N_VALUES:
        raise ValueError(f"実験O6のNは {N_VALUES} のいずれかです: {n}")
    verify_core_hashes()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(HERE))
    runner = load_module("o6_runner_control", HERE / "run_spontaneous_splitting_largeN_v1.py")
    runner.RESULT_DIR = str(RESULT_DIR)
    print(
        f"=== O6原本対照 N={n}, delta={DELTA}, seed={SEED}, "
        f"cap={CAP}, after={AFTER}, tol={TOL} ===",
        flush=True,
    )
    runner.run(n, DELTA, SEED, CAP, AFTER, tol=TOL)
    report = compare_one(n)
    if not report["passed"]:
        raise RuntimeError(f"N={n} は原本保存出力との厳密比較に失敗しました")
    return report


def import_reference_one(n):
    if n not in N_VALUES:
        raise ValueError(f"実験O6のNは {N_VALUES} のいずれかです: {n}")
    verify_core_hashes()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    tag = tag_for(n)
    names = (
        f"fcurve_{tag}.csv",
        f"summary_{tag}.json",
    )
    copied = []
    for name in names:
        source = REFERENCE_DIR / name
        destination = RESULT_DIR / name
        shutil.copyfile(source, destination)
        copied.append({
            "source": str(source),
            "destination": str(destination),
            "sha256": sha256(destination),
        })
    provenance = {
        "n": n,
        "status": "published reference imported; not re-executed",
        "reason": "N=1000 control rerun skipped by research decision on 2026-07-24",
        "files": copied,
    }
    write_json(RESULT_DIR / f"provenance_reference_{tag}.json", provenance)
    report = compare_one(n)
    if not report["passed"]:
        raise RuntimeError(f"N={n} の原論文参照値の導入に失敗しました")
    print(f"[OK] N={n} は再実行せず、原論文の保存系列を来歴付きで導入しました")
    return report


def make_plot():
    verify_core_hashes()
    for n in N_VALUES:
        report = compare_one(n)
        if not report["passed"]:
            raise RuntimeError(f"N={n} の対照一致前には集約図を生成しません")
    plotter = load_module("o6_plotter_control", HERE / "make_largeN_figure_v1.py")
    plotter.DIR = str(RESULT_DIR)
    plotter.main()
    output = RESULT_DIR / "dormant_growth_large_n_v1.png"
    if not output.exists() or output.stat().st_size == 0:
        raise RuntimeError("集約図が生成されませんでした")
    print(f"[OK] 原本3系列だけから集約図を生成: {output}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=("validate", "audit", "run", "reference", "plot", "all"),
    )
    parser.add_argument("n", nargs="?", type=int)
    args = parser.parse_args()
    if args.command in ("run", "reference") and args.n is None:
        parser.error(f"{args.command} には N=40, 300, 1000 のいずれかが必要です")
    if args.command not in ("run", "reference") and args.n is not None:
        parser.error("Nを指定できるのは run と reference だけです")
    return args


def main():
    args = parse_args()
    if args.command == "validate":
        validate_dense()
    elif args.command == "audit":
        audit_existing_copy()
    elif args.command == "run":
        run_one(args.n)
    elif args.command == "reference":
        import_reference_one(args.n)
    elif args.command == "plot":
        make_plot()
    elif args.command == "all":
        validate_dense()
        for n in (40, 300):
            run_one(n)
        import_reference_one(1000)
        make_plot()


if __name__ == "__main__":
    main()
