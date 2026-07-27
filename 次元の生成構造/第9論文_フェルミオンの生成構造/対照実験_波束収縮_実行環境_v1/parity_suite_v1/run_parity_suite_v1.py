#!/usr/bin/env python3
"""対照テストスイート v1: 原本実験のパラメータを完全復元してコピー環境で再実行し、
原本CSVとの完全一致を検証する。

- 原本パラメータは原本 result_v1.json の記録（r_values / cases / fixed_norm / run_id）から機械的に復元する。
- 再実行はミラーコピー環境（本フォルダの親）のランナーで行い、出力は本フォルダ配下に隔離保存する。
- 比較は CSV をバイト単位で行い、不一致時は行・列単位の差分を記録する。
- PNG は --no-plots で生成しない（データ比較に不要）。
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent          # parity_suite_v1
ENV_DIR = HERE.parent                            # 対照実験_波束収縮_実行環境_v1
REPO = ENV_DIR.parent.parent.parent              # リポジトリルート
RUNNER = ENV_DIR / "20260715" / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py"
ORIG_ROOT = REPO / "波の情報読出し" / "20260715" / "system_A_localization_exchange_R_sweep_result_v1"

CSV_KINDS = ("rows", "summary", "best", "collision_terrain")


def fmt(value: float) -> str:
    return repr(float(value))


def build_tests() -> list[dict]:
    tests = []

    # --- Test 1: odd_kernel N=1,2,3,5,15,63 (既定Rセット63点, run_id なし, 明示 file_stem) ---
    orig_json = ORIG_ROOT / "odd_kernel_N_1_2_3_5_15_63_result_v1.json"
    d = json.loads(orig_json.read_text())
    tests.append({
        "name": "01_odd_kernel_N_1_2_3_5_15_63",
        "orig_json": orig_json,
        "orig_dir": ORIG_ROOT,
        "stem": d["file_stem"],
        "args": [
            "--pairs", "1:1,1:2,1:3,1:5,1:15,1:63",
            "--file-stem", d["file_stem"],
        ],
        "expect_r_values": d["r_values"],
    })

    # --- Test 2: complete_femtofocus_R137_B12_fixed_global_norm (R137極細121点, 固定正規化) ---
    orig_dir = ORIG_ROOT / "complete_femtofocus_R137_B12_fixed_global_norm"
    orig_json = next(orig_dir.glob("*result_v1.json"))
    d = json.loads(orig_json.read_text())
    tests.append({
        "name": "02_complete_femtofocus_R137_B12_fixed_global_norm",
        "orig_json": orig_json,
        "orig_dir": orig_dir,
        "stem": d["file_stem"],
        "args": [
            "--packet-a", "1", "--packet-b", "1,2",
            "--r-values", ",".join(fmt(v) for v in d["r_values"]),
            "--fixed-l-norm", fmt(d["fixed_l_norm"]),
            "--fixed-n-norm", fmt(d["fixed_n_norm"]),
            "--run-id", d["run_id"],
        ],
        "expect_r_values": d["r_values"],
    })

    # --- Test 3: phase10_B12 (位相ずれ対照, 既定Rセット) ---
    orig_dir = ORIG_ROOT / "phase10_B12"
    orig_json = next(orig_dir.glob("*result_v1.json"))
    d = json.loads(orig_json.read_text())
    case = d["cases"][0]
    tests.append({
        "name": "03_phase10_B12",
        "orig_json": orig_json,
        "orig_dir": orig_dir,
        "stem": d["file_stem"],
        "args": [
            "--packet-a", "1", "--packet-b", "1,2",
            "--packet-b-phases", ",".join(fmt(v) for v in case["packet_b_phases"]),
            "--run-id", d["run_id"],
        ],
        "expect_r_values": d["r_values"],
    })

    return tests


def compare_csv(orig: Path, new: Path) -> dict:
    if orig.read_bytes() == new.read_bytes():
        return {"identical": True}
    # 不一致の場合はセル単位差分を数える
    with open(orig) as fo, open(new) as fn:
        ro = list(csv.reader(fo))
        rn = list(csv.reader(fn))
    detail: dict = {"identical": False, "orig_rows": len(ro), "new_rows": len(rn), "cell_diffs": []}
    for i, (a, b) in enumerate(zip(ro, rn)):
        for j, (x, y) in enumerate(zip(a, b)):
            if x != y:
                detail["cell_diffs"].append({"row": i, "col": ro[0][j] if j < len(ro[0]) else j, "orig": x, "new": y})
                if len(detail["cell_diffs"]) >= 20:
                    detail["truncated"] = True
                    return detail
    return detail


def main() -> None:
    report: dict = {"suite": "parity_suite_v1", "tests": []}
    for test in build_tests():
        name = test["name"]
        out_dir = HERE / name / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(RUNNER), *test["args"], "--output-dir", str(out_dir), "--no-plots"]
        print(f"=== {name} ===", flush=True)
        print(" ".join(cmd[:6]) + " ...", flush=True)
        log_path = HERE / name / "run_log.txt"
        with open(log_path, "w") as log:
            log.write(" ".join(cmd) + "\n\n")
            log.flush()
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        entry: dict = {"name": name, "returncode": proc.returncode, "command": cmd, "files": {}}
        if proc.returncode != 0:
            entry["error"] = "runner failed; see run_log.txt"
            report["tests"].append(entry)
            print(f"  RUNNER FAILED rc={proc.returncode}", flush=True)
            continue
        # R値集合の一致確認
        new_json = next(out_dir.glob("*result_v1.json"))
        nd = json.loads(new_json.read_text())
        entry["r_values_match"] = nd["r_values"] == test["expect_r_values"]
        # CSV比較
        all_ok = True
        for kind in CSV_KINDS:
            orig_csv = test["orig_dir"] / f"{test['stem']}_{kind}_v1.csv"
            new_csv = out_dir / f"{test['stem']}_{kind}_v1.csv"
            if not orig_csv.exists() or not new_csv.exists():
                entry["files"][kind] = {"missing": {"orig": not orig_csv.exists(), "new": not new_csv.exists()}}
                all_ok = False
                continue
            result = compare_csv(orig_csv, new_csv)
            entry["files"][kind] = result
            all_ok = all_ok and result["identical"]
        entry["verdict"] = "EXACT_MATCH" if (all_ok and entry["r_values_match"]) else "MISMATCH"
        print(f"  verdict: {entry['verdict']}", flush=True)
        report["tests"].append(entry)

    out = HERE / "parity_suite_result_v1.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(f"saved: {out}", flush=True)
    print(json.dumps([{t['name']: t.get('verdict', 'ERROR')} for t in report["tests"]], ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
