#!/usr/bin/env python3
"""本番ダンプ v1: 計測版ランナーで倍音別複素係数を系統的に取得する。

対象（波束収縮と交差項の関係性調査.md §3.1 の優先系列に対応）:
  01_femtofocus_R137_B12  波束収縮の中心窓。原本 femtofocus と同一の R121点・固定正規化を
                          原本JSONから復元して再実行（既存列は原本とバイト一致するはず＝内蔵parity）。
  02_B12_keyR             同ケースを代表9点のRで取得（対照 0/0.5/1.0、R137厳密値、R128、
                          femto窓中心、off-resonance 0.68/0.70/0.71）。
  03_oddN_B{n}_keyR       倍音次数比較 N=1,2,3,5,15,63（odd_kernel、ペア 1:n）を同じ代表9点で。
                          ダンプ倍音範囲 M は各ケースのパケット支持 n+2 に絞る
                          （線形発展は新周波数を作らない。npz 内 coverage で無切り捨てを監査）。

全実行で --dump-stride 1（全衝突記録。選択直前時系列の解析要件を満たすため）。
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent           # production_dump_v1
ENV_DIR = HERE.parent
REPO = ENV_DIR.parent.parent.parent
RUNNER = ENV_DIR / "20260715" / "run_system_A_localization_exchange_R_sweep_instrumented_v1.py"
ORIG_FEMTO = next((REPO / "波の情報読出し/20260715/system_A_localization_exchange_R_sweep_result_v1/complete_femtofocus_R137_B12_fixed_global_norm").glob("*result_v1.json"))

R137_EXACT = 0.6971778791282474
R128_EXACT = 0.686671465671125
KEY_R = [0.0, 0.5, 0.68, R128_EXACT, R137_EXACT, 0.697177927, 0.7, 0.71, 1.0]


def fmt(v: float) -> str:
    return repr(float(v))


def build_runs() -> list[dict]:
    femto = json.loads(ORIG_FEMTO.read_text())
    runs = [
        {
            "name": "01_femtofocus_R137_B12",
            "args": [
                "--packet-a", "1", "--packet-b", "1,2",
                "--r-values", ",".join(fmt(v) for v in femto["r_values"]),
                "--fixed-l-norm", fmt(femto["fixed_l_norm"]),
                "--fixed-n-norm", fmt(femto["fixed_n_norm"]),
                "--run-id", femto["run_id"],
                "--dump-max-n", "6",
            ],
        },
        {
            "name": "02_B12_keyR",
            "args": [
                "--packet-a", "1", "--packet-b", "1,2",
                "--r-values", ",".join(fmt(v) for v in KEY_R),
                "--dump-max-n", "6",
            ],
        },
    ]
    for n in (1, 2, 3, 5, 15, 63):
        runs.append({
            "name": f"03_oddN_B{n}_keyR",
            "args": [
                "--pairs", f"1:{n}",
                "--r-values", ",".join(fmt(v) for v in KEY_R),
                "--dump-max-n", str(n + 2),
            ],
        })
    return runs


def main() -> None:
    manifest = {"suite": "production_dump_v1", "runner": RUNNER.name, "runs": []}
    for run in build_runs():
        out_dir = HERE / run["name"] / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(RUNNER), *run["args"], "--output-dir", str(out_dir), "--no-plots", "--dump-stride", "1"]
        log_path = HERE / run["name"] / "run_log.txt"
        print(f"=== {run['name']} ===", flush=True)
        started = time.time()
        with open(log_path, "w") as log:
            log.write(" ".join(cmd) + "\n\n")
            log.flush()
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        elapsed = round(time.time() - started, 1)
        entry = {"name": run["name"], "returncode": proc.returncode, "elapsed_s": elapsed, "command": cmd}
        if proc.returncode == 0:
            dump_dir = out_dir / "harmonic_dump_v1"
            npzs = sorted(dump_dir.glob("*.npz"))
            entry["npz_count"] = len(npzs)
            entry["dump_bytes"] = sum(p.stat().st_size for p in npzs)
            dm = dump_dir / "harmonic_dump_manifest_v1.json"
            if dm.exists():
                dj = json.loads(dm.read_text())
                covs = [e.get("coverage_min", float("nan")) for e in dj.get("entries", [])]
                entry["coverage_min"] = min(covs) if covs else None
        print(f"  rc={proc.returncode} elapsed={elapsed}s npz={entry.get('npz_count')} cov_min={entry.get('coverage_min')}", flush=True)
        manifest["runs"].append(entry)
    (HERE / "production_dump_manifest_v1.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    ok = all(r["returncode"] == 0 for r in manifest["runs"])
    print(json.dumps({"all_ok": ok, "total_MB": round(sum(r.get("dump_bytes", 0) for r in manifest["runs"]) / 1e6, 1)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
