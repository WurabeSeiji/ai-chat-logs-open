#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL1 真空読出し v1 — 保存則の定理（補題DL1-2・命題DL1-3）の統一エンジン検証（判定 K1〜K6）

導出: DL1_導出ノート.md。走行は DL0 と同一（決定的・再走行）。
  K1 B𝟙=0（恒等・全フレーム）
  K2 Σ|x_e|² の変動幅（命題DL1-3: 共有O＋単一スライスで厳密保存）
  K3 Σx_e² の変動幅（同・DL0-2）
  K4 奇数帯パワー（系: 真空に物質は湧かない。v1 で厳密0）
  K5 方向読出しの無意味性: ラグ120步の上位3重なり（DL0 J8 と同一測定・値の記録）
  K6 tr(B)=c/N（恒等）
出力: result_dl1_readout_v1.json
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
N, M, T, SAMPLE = 16, 120, 20000, 10


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    u1 = _load("uni_dl1", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, 0.0)
    J = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)

    c0 = None
    bil0 = None
    K = {"B1_max": 0.0, "c_drift": 0.0, "bil_drift": 0.0,
         "odd_max": 0.0, "trB_rel_max": 0.0}
    V3s = []
    for tau in range(T):
        eng.step()
        if tau % SAMPLE:
            continue
        C2 = eng.C2()
        K["odd_max"] = max(K["odd_max"],
                           float((np.abs(C2[:, eng.odd_k, :]) ** 2).sum()))
        x = C2.sum(axis=(1, 2))
        c = float(np.sum(np.abs(x) ** 2))
        bil = complex(np.sum(x * x))
        if c0 is None:
            c0, bil0 = c, bil
        K["c_drift"] = max(K["c_drift"], abs(c - c0))
        K["bil_drift"] = max(K["bil_drift"], abs(bil - bil0))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * J @ D2 @ J
        K["B1_max"] = max(K["B1_max"], float(np.max(np.abs(B @ np.ones(N)))))
        K["trB_rel_max"] = max(K["trB_rel_max"], abs(np.trace(B) * N / c - 1.0))
        lam, V = np.linalg.eigh(B)
        ones = np.ones(N) / np.sqrt(N)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        V3s.append(Vk[:, o][:, :3].copy())

    LAG = 12
    ov = np.array([np.sum((V3s[i].T @ V3s[i + LAG]) ** 2) / 3.0
                   for i in range(len(V3s) - LAG)])
    late = ov[int(len(ov) * 0.75):]

    res = {
        "K1_B1": {"max": K["B1_max"], "pass": bool(K["B1_max"] < 1e-13)},
        "K2_hermitian": {"drift": K["c_drift"], "pass": bool(K["c_drift"] < 1e-12)},
        "K3_bilinear": {"drift": K["bil_drift"], "pass": bool(K["bil_drift"] < 1e-12)},
        "K4_odd_band": {"max": K["odd_max"], "pass": bool(K["odd_max"] == 0.0)},
        "K5_overlap_lag120": {"late_mean": float(late.mean()),
                              "late_min": float(late.min()),
                              "late_max": float(late.max()),
                              "note": "枠非持続の記録（DL0 J8 と同一測定）。"
                                      "持続閾値0.5を下回ることの確認",
                              "pass": bool(late.max() < 0.5)},
        "K6_trB": {"rel_max": K["trB_rel_max"],
                   "pass": bool(K["trB_rel_max"] < 1e-6)},
        "elapsed_sec": time.time() - t0,
    }
    res["all_pass"] = bool(all(res[k]["pass"] for k in
                               ("K1_B1", "K2_hermitian", "K3_bilinear",
                                "K4_odd_band", "K5_overlap_lag120", "K6_trB")))
    (HERE / "result_dl1_readout_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k, v in res.items():
        if isinstance(v, dict):
            print(f"  {k}: {'PASS' if v['pass'] else 'FAIL'}  {v}")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
