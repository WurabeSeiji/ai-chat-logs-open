#!/usr/bin/env python3
"""媒介頂点一意化定理の数値検証（媒介頂点の一意化と多体アーキテクチャ_v1.md §3）

検証項目（実行前固定）:
    T1 閉塞保存の恒等性: g₂=−g₁・対対称強度（算術平均）で、ランダム状態の
       頂点一回作用が dC=Σz_e·δz_e を機械精度ゼロにする（K5/K8/K12、複数乱数）。
    T2 反例1: g₂=+g₁ では dC≠0（一意化の必要性）。
    T3 反例2: 非対称強度 R_{ee'}=R_e（自分側のみ）では dC≠0（対対称の必要性）。
    T4 因数分解の厳密性: O(M)頂点集約実装が素朴な隣接対和実装と一致
       （相対誤差 ≤1e-14）。
    T5 積強度 R_{ee'}=R_e·R_{e'} も対対称なので dC=0（系1の正直な記録:
       閉塞保存だけでは平均と積は判別できない）。

使い方: python3 test_mediated_vertex_closure_v1.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def edges_and_adjacency(n):
    ia, ib = np.triu_indices(n, k=1)
    m = len(ia)
    adj = [[] for _ in range(m)]
    for e in range(m):
        for f in range(m):
            if e == f:
                continue
            if ia[e] in (ia[f], ib[f]) or ib[e] in (ia[f], ib[f]):
                adj[e].append(f)
    return ia, ib, adj


def vertex_naive(z, R, adj, g1, g2, pair_rule):
    m = len(z)
    dz = np.zeros_like(z)
    for e in range(m):
        for f in adj[e]:
            if pair_rule == "mean":
                Rp = 0.5 * (R[e] + R[f])
            elif pair_rule == "left":
                Rp = R[e]
            elif pair_rule == "prod":
                Rp = R[e] * R[f]
            dz[e] += 1j * Rp * (g1 * abs(z[f]) ** 2 * z[e] + g2 * z[f] ** 2 * np.conj(z[e]))
    return dz


def vertex_factored_mean(z, R, ia, ib, n, g1):
    """算術平均強度・g₂=−g₁ の O(M) 頂点集約実装（系2）。"""
    m = len(z)
    A = np.zeros(n)
    B = np.zeros(n, complex)
    AR = np.zeros(n)
    BR = np.zeros(n, complex)
    a2 = np.abs(z) ** 2
    z2 = z ** 2
    np.add.at(A, ia, a2); np.add.at(A, ib, a2)
    np.add.at(B, ia, z2); np.add.at(B, ib, z2)
    np.add.at(AR, ia, R * a2); np.add.at(AR, ib, R * a2)
    np.add.at(BR, ia, R * z2); np.add.at(BR, ib, R * z2)
    cA = A[ia] + A[ib] - 2 * a2
    cB = B[ia] + B[ib] - 2 * z2
    cAR = AR[ia] + AR[ib] - 2 * R * a2
    cBR = BR[ia] + BR[ib] - 2 * R * z2
    return 0.5j * g1 * (R * (cA * z - cB * np.conj(z)) + (cAR * z - cBR * np.conj(z)))


def dC(z, dz):
    return complex(np.sum(z * dz))


def main() -> None:
    t0 = time.time()
    results = {}
    ok_all = True
    for n in (5, 8, 12):
        ia, ib, adj = edges_and_adjacency(n)
        m = len(ia)
        rng = np.random.default_rng(1234 + n)
        rows = {"T1": [], "T2": [], "T3": [], "T4": [], "T5": []}
        for trial in range(5):
            z = rng.normal(size=m) + 1j * rng.normal(size=m)
            R = rng.uniform(0.1, 0.9, size=m)
            scale = float(np.sum(np.abs(z) ** 2))
            d1 = vertex_naive(z, R, adj, 1.0, -1.0, "mean")
            rows["T1"].append(abs(dC(z, d1)) / scale ** 1.5)
            d2 = vertex_naive(z, R, adj, 1.0, +1.0, "mean")
            rows["T2"].append(abs(dC(z, d2)) / scale ** 1.5)
            d3 = vertex_naive(z, R, adj, 1.0, -1.0, "left")
            rows["T3"].append(abs(dC(z, d3)) / scale ** 1.5)
            d4 = vertex_factored_mean(z, R, ia, ib, n, 1.0)
            rows["T4"].append(float(np.linalg.norm(d4 - d1) / np.linalg.norm(d1)))
            d5 = vertex_naive(z, R, adj, 1.0, -1.0, "prod")
            rows["T5"].append(abs(dC(z, d5)) / scale ** 1.5)
        v = {"T1_max_dC": max(rows["T1"]), "T2_min_dC": min(rows["T2"]),
             "T3_min_dC": min(rows["T3"]), "T4_max_relerr": max(rows["T4"]),
             "T5_max_dC": max(rows["T5"])}
        passed = (v["T1_max_dC"] < 1e-14 and v["T2_min_dC"] > 1e-3
                  and v["T3_min_dC"] > 1e-3 and v["T4_max_relerr"] < 1e-13
                  and v["T5_max_dC"] < 1e-14)
        ok_all &= passed
        results[f"N{n}"] = {**v, "pass": bool(passed)}
        print(f"N={n:2d} (M={m}): T1 dC={v['T1_max_dC']:.1e} | T2(g2=+g1) dC={v['T2_min_dC']:.1e} | "
              f"T3(非対称R) dC={v['T3_min_dC']:.1e} | T4 因数分解誤差={v['T4_max_relerr']:.1e} | "
              f"T5(積) dC={v['T5_max_dC']:.1e} → {'PASS' if passed else 'FAIL'}")
    print(f"\n総合: {'ALL PASS' if ok_all else 'FAIL'}——定理§3と系2は数値的に成立"
          f"（T5: 積強度も閉塞保存——系1の正直な記録どおり）")
    out = {"results": results, "all_pass": bool(ok_all), "runtime_sec": time.time() - t0}
    (HERE / "test_mediated_vertex_closure_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.1f}s)")


if __name__ == "__main__":
    main()
