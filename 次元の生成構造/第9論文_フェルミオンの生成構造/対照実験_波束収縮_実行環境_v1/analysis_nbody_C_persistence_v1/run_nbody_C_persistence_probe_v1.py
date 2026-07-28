#!/usr/bin/env python3
"""N体系 交差項持続性プローブ v1

問い: 閉鎖系 Σxₙ² = B² において A²−B² = 2C の交差項 C は「普通は時間平均で消える」
      とされる。第8論文の自然発生 N=300 準安定状態（t=110000 チェックポイント）を
      既存力学そのままで発展させ、C(τ) が消えるか残るかを直接測る。

方法: 既存 checkpoint npz（Z 複素 44850 成分, wp 付き）から、第7/8論文と同一の
      更新ループ（set_theta(angle(Z)) → sigma_max_power → cayley_step）で継続発展。
      力学・正規化・パラメータは一切変更しない。計測のみ。

観測量（毎ステップ）:
  A  = Σᵢ Zᵢ          （合成波: 全成分の単純和）
  B2 = Σᵢ Zᵢ²         （双線形閉鎖量: この系では ≈0 の零二乗和）
  C  = (A² − B2) / 2   （双線形交差項総和）
  f  = 休眠フラクション（文脈用、largeN スクリプトと同定義）

判定量:
  持続率 |⟨C⟩_T| / ⟨|C|⟩_T   （1に近い=位相が揃ったまま残る、0に近い=平均消去）
  区間別持続率（T/4 ごと）と |C| の減衰有無
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent.parent
ENGINE_PATH = REPO / "次元の生成構造/自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py"
CKPT = REPO / ("次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/"
               "paper7_seedless_natural_figures3_4_v1/outputs/long_horizon_110000/"
               "checkpoints/N00300_state_t110000.npz")
STEPS = 2000

spec = importlib.util.spec_from_file_location("lowrank_engine", ENGINE_PATH)
eng = importlib.util.module_from_spec(spec)
sys.modules["lowrank_engine"] = eng
spec.loader.exec_module(eng)


def main() -> None:
    ck = np.load(CKPT)
    n, step0 = int(ck["N"]), int(ck["step"])
    Z = ck["Z"].astype(np.complex128)
    wp = ck["wp"].astype(np.float64)
    sys_lr = eng.LowRankSystem(n)
    print(f"checkpoint: N={n} step={step0} m={len(Z)}", flush=True)

    # 親平面基底（f の文脈計測用、largeN スクリプトと同じ構成）
    p = Z.real / np.linalg.norm(Z.real)
    q = Z.imag - (Z.imag @ p) * p
    q = q / np.linalg.norm(q)

    rows = []
    t0 = time.time()
    for t in range(STEPS + 1):
        A = complex(Z.sum())
        B2 = complex(Z @ Z)
        C = (A * A - B2) / 2.0
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        htot = float(np.real(np.conj(Z) @ Z))
        f = float(np.real(np.conj(Zp) @ Zp)) / htot
        rows.append([step0 + t, A.real, A.imag, B2.real, B2.imag, C.real, C.imag, abs(C), f])
        if t >= STEPS:
            break
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
        if t % 500 == 0 and t > 0:
            print(f"  t={t} ({(time.time()-t0)/t*1000:.1f} ms/step)", flush=True)

    out_csv = HERE / f"nbody_C_persistence_N{n:05d}_from{step0}_steps{STEPS}_v1.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "A_re", "A_im", "B2_re", "B2_im", "C_re", "C_im", "C_abs", "f"])
        w.writerows(rows)

    arr = np.array(rows)
    Cc = arr[:, 5] + 1j * arr[:, 6]
    absC = arr[:, 7]
    persistence = abs(Cc.mean()) / max(absC.mean(), 1e-300)
    quarters = []
    for i in range(4):
        seg = Cc[i * len(Cc) // 4:(i + 1) * len(Cc) // 4]
        segabs = absC[i * len(absC) // 4:(i + 1) * len(absC) // 4]
        quarters.append(abs(seg.mean()) / max(segabs.mean(), 1e-300))
    summary = {
        "N": n, "start_step": step0, "steps": STEPS,
        "mean_absC": float(absC.mean()), "final_absC": float(absC[-1]),
        "mean_absB2": float(np.abs(arr[:, 3] + 1j * arr[:, 4]).mean()),
        "persistence_ratio_total": float(persistence),
        "persistence_ratio_quarters": [float(x) for x in quarters],
        "absC_first_quarter_mean": float(absC[: len(absC) // 4].mean()),
        "absC_last_quarter_mean": float(absC[-len(absC) // 4:].mean()),
        "f_range": [float(arr[:, 8].min()), float(arr[:, 8].max())],
        "csv": out_csv.name,
    }
    (HERE / f"nbody_C_persistence_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
