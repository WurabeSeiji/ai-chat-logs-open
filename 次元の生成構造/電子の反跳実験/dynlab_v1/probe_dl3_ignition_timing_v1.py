#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 点火の同時性: τ_g（ギャップ開通）と τ_−（虚方向初出）の毎步測定

背景: run_dl23_matter_v1（10步サンプル）の図 H5 で、g_min の跳躍と負固有値の出現が
τ≈9200 で同時に見えた。本プローブは同一の決定的軌道を毎步分解能で測り、
両時刻を事前登録した定義で確定する。

事前登録定義（実行前固定）:
  基線 g_b = median(g_min, τ∈[8000,8500])（点火前窓）
  τ_g = g_min > 10·g_b が50步連続で成立する最初の步
  τ_− = n_neg > 0 が50步連続で成立する最初の步（τ>8000）
  Δτ = τ_g − τ_−。|Δτ| ≤ 1 なら「同一相転移」、それ以外は両時刻を記録。

出力: result_dl3_ignition_timing_v1.json・dl3_ignition_series_v1.npz
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
T_SKIP, T_END = 8000, 10500


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def first_sustained(tau, mask, run=50):
    """mask が run 步連続で立つ最初の tau（なければ None）。"""
    cnt = 0
    for t, m in zip(tau, mask):
        cnt = cnt + 1 if m else 0
        if cnt >= run:
            return int(t - run + 1)
    return None


def main():
    t0 = time.time()
    u1 = _load("uni_dl3ig", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    for _ in range(T_SKIP):
        eng.step()

    taus, gmins, nnegs = [], [], []
    for tau in range(T_SKIP, T_END):
        eng.step()
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        l = np.sort(lam[keep])[::-1]
        gmins.append(float(min(l[0] - l[1], l[1] - l[2], l[2] - l[3])))
        nnegs.append(int(np.sum(l < -1e-12 * max(float(l.sum()), 1e-300))))
        taus.append(tau + 1)

    taus = np.array(taus)
    gmins = np.array(gmins)
    nnegs = np.array(nnegs)
    base = taus <= 8500
    g_b = float(np.median(gmins[base]))
    tau_g = first_sustained(taus, gmins > 10 * g_b)
    tau_m = first_sustained(taus, nnegs > 0)
    dtau = (tau_g - tau_m) if (tau_g is not None and tau_m is not None) else None

    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "window": [T_SKIP, T_END], "run_required": 50,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行"},
        "g_baseline_median": g_b,
        "tau_g": tau_g, "tau_minus": tau_m, "delta_tau": dtau,
        "same_transition": (dtau is not None and abs(dtau) <= 1),
        "elapsed_sec": time.time() - t0,
    }
    np.savez_compressed(HERE / "dl3_ignition_series_v1.npz",
                        tau=taus, gmin=gmins, n_neg=nnegs)
    (HERE / "result_dl3_ignition_timing_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"点火時刻: τ_g={tau_g}  τ_−={tau_m}  Δτ={dtau}  "
          f"同一相転移(|Δτ|≤1)={res['same_transition']}  "
          f"基線g_b={g_b:.2e} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
