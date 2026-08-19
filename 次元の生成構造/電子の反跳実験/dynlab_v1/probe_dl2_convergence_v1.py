#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL2 追走行 v1 — 保存ドリフト 10⁻⁹ の収束試験（再査読対応・実行前固定）

主張「残るドリフトは頂点部（RK4）の積分器誤差」を事実にする。
頂点積分の刻み上限 s2.H_MAX（正本値 0.02）を {H0, H0/2, H0/4} と振り、
同一走行（N=16, δ=0.1）の窓 τ∈[9000, 10500]（点火を含む）での
双線形レジスタのドリフト E(h) = max|Σx²(τ) − Σx²(9000)| を測る。
RK4 が支配なら E(h) ∝ h⁴（刻み半減で約1/16）。

注: H_MAX はモジュール定数であり、正本プログラムは変更しない（実行時に設定）。

出力: result_dl2_convergence_v1.json
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
T_A, T_B = 9000, 10500
H0 = 0.02


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def drift_for(hmax):
    u1 = _load(f"uni_cv_{hmax}", UNI / "unified_interaction_v1.py")
    u1.s2.H_MAX = hmax
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    for _ in range(T_A):
        eng.step()
    x = eng.C2().sum(axis=(1, 2))
    ref = complex(np.sum(x * x))
    d = 0.0
    for _ in range(T_B - T_A):
        eng.step()
        x = eng.C2().sum(axis=(1, 2))
        d = max(d, abs(complex(np.sum(x * x)) - ref))
    return d


def main():
    t0 = time.time()
    hs = [H0, H0 / 2, H0 / 4]
    ds = []
    for h in hs:
        d = drift_for(h)
        ds.append(d)
        print(f"  H_MAX={h}: E={d:.3e}")
    r1 = ds[0] / ds[1] if ds[1] > 0 else float("inf")
    r2 = ds[1] / ds[2] if ds[2] > 0 else float("inf")
    p1 = float(np.log2(r1)) if np.isfinite(r1) and r1 > 0 else None
    p2 = float(np.log2(r2)) if np.isfinite(r2) and r2 > 0 else None
    res = {
        "config": {"N": N, "delta": DELTA, "window": [T_A, T_B],
                   "H_MAX_values": hs,
                   "note": "s2.H_MAX を実行時設定（正本は不変更）"},
        "drift": dict(zip(map(str, hs), ds)),
        "ratio_h_to_h2": r1, "ratio_h2_to_h4": r2,
        "order_est_1": p1, "order_est_2": p2,
        "verdict_integrator_error": bool(p1 is not None and p1 > 2.5
                                         and p2 is not None and p2 > 2.5),
        "elapsed_sec": time.time() - t0,
    }
    (HERE / "result_dl2_convergence_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"収束比: E(h)/E(h/2)={r1:.1f}（次数≈{p1}）  "
          f"E(h/2)/E(h/4)={r2:.1f}（次数≈{p2}）")
    print(f"積分器誤差と決着: {res['verdict_integrator_error']} "
          f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
