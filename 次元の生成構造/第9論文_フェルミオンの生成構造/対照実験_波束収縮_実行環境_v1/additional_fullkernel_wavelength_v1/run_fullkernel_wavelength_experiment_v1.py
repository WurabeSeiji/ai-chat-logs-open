#!/usr/bin/env python3
"""追加実験: フルカーネル（強局在B）× 波長ずれ 決定版対照 v1

目的（考察md §7 の唯一の未走行対照）:
  原本の波長ずれ系列は B={1,63} の2成分で B 自身の局在が弱かった。
  本実験は B をフル奇数カーネル {1,3,5,...,63}（32成分、λ=1 で odd_kernel B=63 と
  同一の強局在状態）とし、波長ずれだけを加えて
  「強く局在した B でも、整数波長整合が壊れると転送がボゾン型へ退化するか」を検証する。

条件:
  λ ∈ {1.0（整合対照）, 1.03, 1.1, 1.3}（原本波長系列と同一の刻み）
  変種 uniform : B 全成分の波長を一様に λ 倍（B内部の整数関係は保持、Aとの整合のみ破壊）
  変種 harmonic: B の基本波は λ=1 のまま、倍音成分のみ λ 倍（原本 {1,63} 系列の直接拡張）

走行:
  1) 既定Rセット63点の掃引（基底ランナー、原本波長系列と同一の観測量）
  2) R137厳密値での複素係数ダンプ（計測版ランナー、M=128、洩れは coverage で監査）

検証（内蔵parity）:
  λ=1.0 の結果が原本 odd_kernel|A=1|B=63 の挙動（L_B≈5.2e-3、N_eff_A→32、
  R*≈0.697177879）を再現すること。
"""

from __future__ import annotations

import csv
import glob
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENV = HERE.parent
BASE_RUNNER = ENV / "20260715" / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py"
INST_RUNNER = ENV / "20260715" / "run_system_A_localization_exchange_R_sweep_instrumented_v1.py"

ODD = [str(n) for n in range(1, 64, 2)]          # 1,3,...,63（32成分）
PACKET_B = ",".join(ODD)
LAMBDAS = [1.0, 1.03, 1.1, 1.3]
R137 = "0.6971778791282474"


def conditions() -> list[dict]:
    conds = []
    for lam in LAMBDAS:
        if lam == 1.0:
            conds.append({"name": "matched_lambda1", "scales": "1"})
            continue
        tag = f"{lam:.2f}".replace(".", "p")
        conds.append({"name": f"uniform_lambda{tag}", "scales": f"{lam!r}"})
        conds.append({"name": f"harmonic_lambda{tag}",
                      "scales": ",".join(["1"] + [repr(lam)] * (len(ODD) - 1))})
    return conds


def run(cmd: list[str], log_path: Path) -> int:
    with open(log_path, "w") as log:
        log.write(" ".join(cmd) + "\n\n")
        log.flush()
        return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode


def main() -> None:
    manifest = {"experiment": "fullkernel_wavelength_v1", "packet_b": PACKET_B, "runs": []}
    for cond in conditions():
        name = cond["name"]
        out_sweep = HERE / name / "sweep"
        out_dump = HERE / name / "dump_R137"
        out_sweep.mkdir(parents=True, exist_ok=True)
        out_dump.mkdir(parents=True, exist_ok=True)
        common = ["--packet-a", "1", "--packet-b", PACKET_B,
                  "--packet-b-wavelength-scales", cond["scales"], "--no-plots"]
        t0 = time.time()
        rc1 = run([sys.executable, str(BASE_RUNNER), *common,
                   "--run-id", f"fullkernel_{name}", "--output-dir", str(out_sweep)],
                  HERE / name / "sweep_log.txt")
        rc2 = run([sys.executable, str(INST_RUNNER), *common,
                   "--r-values", R137, "--dump-max-n", "128",
                   "--run-id", f"fullkernel_{name}_R137", "--output-dir", str(out_dump)],
                  HERE / name / "dump_log.txt")
        elapsed = round(time.time() - t0, 1)

        entry = {"name": name, "scales": cond["scales"], "rc_sweep": rc1, "rc_dump": rc2,
                 "elapsed_s": elapsed}
        # 掃引の要約値を吸い上げ
        if rc1 == 0:
            terrain = glob.glob(str(out_sweep / "*collision_terrain_v1.csv"))
            if terrain:
                mLA = mNA = mLB = 0.0
                with open(terrain[0]) as f:
                    for row in csv.DictReader(f):
                        mLA = max(mLA, float(row["L_A"]))
                        mNA = max(mNA, float(row["N_eff_A"]))
                        mLB = max(mLB, float(row["L_B"]))
                entry.update({"max_L_A": mLA, "max_N_eff_A": mNA, "max_L_B": mLB})
        if rc2 == 0:
            dm = glob.glob(str(out_dump / "harmonic_dump_v1" / "harmonic_dump_manifest_v1.json"))
            if dm:
                dj = json.loads(Path(dm[0]).read_text())
                covs = [e.get("coverage_min") for e in dj.get("entries", [])]
                entry["dump_coverage_min"] = min(covs) if covs else None
        manifest["runs"].append(entry)
        print(f"{name}: sweep rc={rc1} dump rc={rc2} {elapsed}s "
              f"L_B={entry.get('max_L_B')} N_eff_A={entry.get('max_N_eff_A')} "
              f"L_A={entry.get('max_L_A')} cov={entry.get('dump_coverage_min')}", flush=True)

    (HERE / "fullkernel_wavelength_manifest_v1.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print("saved: fullkernel_wavelength_manifest_v1.json", flush=True)


if __name__ == "__main__":
    main()
