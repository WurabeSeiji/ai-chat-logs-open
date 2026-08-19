#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 毎步分解能の輸送比 ρ_n = ‖ΔB‖₂ / g_min（判定 M3 の二段目）

背景: run_dl23_matter_v1.py（10步間隔）では ρ 中央値 3.74・全フレーム域外＝追跡不能。
本プローブは同一の物質走行（N=16, δ=0.1, F=v1, 決定的）を再走行し、
物質相後期（15000步以降）で毎步の ΔB と g_min を測る。
Davis–Kahan（[F2] 命題・第2編 [D2]）: ρ_n < 1/2 の步で枠の連続輸送が定理保証される。

出力: result_dl3_transport_perstep_v1.json・dl3_perstep_series_v1.npz
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent.parent / "統一万能関数_v1"
N, DELTA = 16, 0.1
T_SKIP = 15000   # 物質相後期へ（本走行と同一の決定的軌道）
T_MEAS = 500     # 毎步測定の步数


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    u1 = _load("uni_dl3p", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)

    def B_and_gmin():
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        ones = np.ones(N) / np.sqrt(N)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        l = np.sort(lam[keep])[::-1]
        return B, float(min(l[0] - l[1], l[1] - l[2], l[2] - l[3]))

    for _ in range(T_SKIP):
        eng.step()

    rho, gmins, dB2s = [], [], []
    B_prev, g_prev = B_and_gmin()
    for _ in range(T_MEAS):
        eng.step()
        B, g = B_and_gmin()
        d = float(np.linalg.norm(B - B_prev, 2))
        rho.append(d / max(g_prev, 1e-300))
        gmins.append(g)
        dB2s.append(d)
        B_prev, g_prev = B, g

    rho = np.array(rho)
    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "t_skip": T_SKIP, "t_meas": T_MEAS,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行"},
        "rho_median": float(np.median(rho)),
        "rho_mean": float(rho.mean()),
        "frac_in_zone": float(np.mean(rho < 0.5)),
        "rho_p10": float(np.percentile(rho, 10)),
        "rho_p90": float(np.percentile(rho, 90)),
        "elapsed_sec": time.time() - t0,
    }
    np.savez_compressed(HERE / "dl3_perstep_series_v1.npz",
                        rho=rho, gmin=np.array(gmins), dB2=np.array(dB2s))
    (HERE / "result_dl3_transport_perstep_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"毎步輸送比（物質相後期 {T_MEAS}步）: 中央値={res['rho_median']:.3f} "
          f"平均={res['rho_mean']:.3f} 域内割合={res['frac_in_zone']:.3f} "
          f"p10/p90={res['rho_p10']:.3f}/{res['rho_p90']:.3f} "
          f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
