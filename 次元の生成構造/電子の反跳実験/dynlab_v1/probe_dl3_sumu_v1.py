#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 判定 M7: Σ_A u_A = 0 の恒等性（命題 DL3-5）

run_dl23_matter_v1.py と同一の決定的物質走行（N=16, δ=0.1, F=v1）を再走行し、
復元配置 X3（逐次整列後）の列和と、有限差分速度 u_A = Δ_n x_A の総和を測る。
理論: B𝟙=0（DL1-1）により非自明固有ベクトルは定数ベクトルと直交
→ X3 は中心化済み → Σx_A = 0 と Σu_A = 0 は恒等（機械精度）。

出力: result_dl3_sumu_v1.json
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
N, DELTA, T, SAMPLE = 16, 0.1, 20000, 10


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    u1 = _load("uni_dl3u", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    X_prev = None
    max_sum_x = 0.0
    max_sum_u = 0.0
    V3_al_prev = None
    for tau in range(T):
        eng.step()
        if tau % SAMPLE:
            continue
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        lamk, Vk = lamk[o], Vk[:, o]
        V3 = Vk[:, :3]
        if V3_al_prev is not None:
            U_, _, Vt_ = np.linalg.svd(V3_al_prev.T @ V3)
            V3 = V3 @ (U_ @ Vt_).T
        V3_al_prev = V3
        X3 = V3 * np.sqrt(np.maximum(lamk[:3], 0.0))[None, :]
        max_sum_x = max(max_sum_x, float(np.max(np.abs(X3.sum(axis=0)))))
        if X_prev is not None:
            u = X3 - X_prev
            max_sum_u = max(max_sum_u, float(np.max(np.abs(u.sum(axis=0)))))
        X_prev = X3

    res = {
        "config": {"N": N, "delta": DELTA, "T": T, "SAMPLE": SAMPLE,
                   "engine": "unified_interaction_v1",
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行"},
        "M7_sum_x_max": max_sum_x,
        "M7_sum_u_max": max_sum_u,
        "pass": bool(max_sum_x < 1e-12 and max_sum_u < 1e-12),
        "elapsed_sec": time.time() - t0,
    }
    (HERE / "result_dl3_sumu_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"M7: {'PASS' if res['pass'] else 'FAIL'}  "
          f"max|Σx|={max_sum_x:.3e}  max|Σu|={max_sum_u:.3e} "
          f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
