#!/usr/bin/env python3
"""N体系 交差項 C の条件別バッテリー v1

問い: ボゾン的条件（閉包しない状態）では、交差項 C は
      (a) そもそも発生しないのか、(b) 発生するが位相が揃わず相殺・揺動するのか。

条件:
  locked   既存チェックポイント（自然発生の閉包状態）: N300 t110000 / N40 t55000,t110000 / N5 t110000
  random   乱数複素状態（単位ノルム、非閉包の代表）を既存力学で発展 — 理論文書の条件B1相当
           ※力学・エンジンは一切変更しない。初期状態のみ乱数（rng固定で再現可能）

判定量（各条件、STEPS ステップ）:
  mean|C|            交差項総和の振幅水準 —「発生しているか」
  CV|C| = std/mean   振幅の変動係数 — 閉包なら ~1e-13（剛体回転）、非閉包なら O(0.1~1)（揺動）
  |⟨C⟩|/⟨|C|⟩        素朴時間平均の残存率 —「相殺するか」
  drift = |C|末端四分位平均 / 冒頭四分位平均 — 振幅の成長/減衰
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent.parent
ENGINE_PATH = REPO / "次元の生成構造/自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py"
CKPT_DIR = REPO / ("次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/"
                   "paper7_seedless_natural_figures3_4_v1/outputs/long_horizon_110000/checkpoints")
STEPS = 2000

spec = importlib.util.spec_from_file_location("lowrank_engine", ENGINE_PATH)
eng = importlib.util.module_from_spec(spec)
sys.modules["lowrank_engine"] = eng
spec.loader.exec_module(eng)


def evolve_and_measure(n: int, Z: np.ndarray, wp: np.ndarray, steps: int) -> dict:
    sys_lr = eng.LowRankSystem(n)
    Cs = np.zeros(steps + 1, dtype=complex)
    B2s = np.zeros(steps + 1, dtype=complex)
    for t in range(steps + 1):
        A = complex(Z.sum())
        B2 = complex(Z @ Z)
        Cs[t] = (A * A - B2) / 2.0
        B2s[t] = B2
        if t >= steps:
            break
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    absC = np.abs(Cs)
    q = len(absC) // 4
    return {
        "mean_absC": float(absC.mean()),
        "cv_absC": float(absC.std() / max(absC.mean(), 1e-300)),
        "naive_avg_ratio": float(abs(Cs.mean()) / max(absC.mean(), 1e-300)),
        "drift_last_over_first_quarter": float(absC[-q:].mean() / max(absC[:q].mean(), 1e-300)),
        "mean_absB2": float(np.abs(B2s).mean()),
        "max_absB2": float(np.abs(B2s).max()),
    }


def main() -> None:
    results = {}

    # --- locked: 既存チェックポイント ---
    for name in ["N00300_state_t110000", "N00040_state_t055000",
                 "N00040_state_t110000", "N00005_state_t110000"]:
        p = CKPT_DIR / f"{name}.npz"
        if not p.exists():
            continue
        ck = np.load(p)
        n = int(ck["N"])
        r = evolve_and_measure(n, ck["Z"].astype(complex), ck["wp"].astype(float), STEPS)
        results[f"locked_{name}"] = {"N": n, **r}
        print(f"locked  {name}: mean|C|={r['mean_absC']:.4f} CV={r['cv_absC']:.2e} "
              f"naive_avg={r['naive_avg_ratio']:.3f} drift={r['drift_last_over_first_quarter']:.4f}", flush=True)

    # --- random: 非閉包対照（既存力学無変更、初期状態のみ乱数） ---
    for n, seed in [(300, 1), (300, 2), (40, 1), (40, 2)]:
        m = n * (n - 1) // 2
        rng = np.random.default_rng(90260728 + seed)
        Z = rng.normal(size=m) + 1j * rng.normal(size=m)
        Z = Z / np.linalg.norm(Z)
        wp = rng.normal(size=m)
        r = evolve_and_measure(n, Z, wp, STEPS)
        results[f"random_N{n:05d}_seed{seed}"] = {"N": n, "seed": seed, **r}
        print(f"random  N={n} seed={seed}: mean|C|={r['mean_absC']:.4f} CV={r['cv_absC']:.2e} "
              f"naive_avg={r['naive_avg_ratio']:.3f} drift={r['drift_last_over_first_quarter']:.4f}", flush=True)

    (HERE / "nbody_C_condition_battery_v1.json").write_text(
        json.dumps({"steps": STEPS, "results": results}, ensure_ascii=False, indent=2), encoding="utf-8")
    print("saved: nbody_C_condition_battery_v1.json", flush=True)


if __name__ == "__main__":
    main()
