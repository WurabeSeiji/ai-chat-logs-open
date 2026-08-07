#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V0-v3: 真空の準備——インフレーション的発展による安定空間の最小N掃引

正しい問い（木原指定・v1/v2の反省後）: 空間は初期条件でなくインフレーションの
到達点。無seed自然軌道（論文7/8の正本）——親凝縮真空だけから、潜伏→急拡大
（crossing）→三方向準安定構造への凝縮——を N ごとに走らせ、
  (i) 急拡大が起きるか（crossing）
  (ii) 三方向準安定構造が形成されるか（rank_Q=4・q2..q4有限）
  (iii) それが持続するか（後期窓での rank_Q=4 持続率＝真空の安定性）
を正本の診断器（qsv4・rank_Q・閉鎖）で測る。

正本: paper7_seedless_natural_figures3_4_v1/（build_seedless・engine=
set_theta+σ_max_power+cayley_step・診断=gram_reduce/dominant_plane/qsv4）。
本スクリプトは正本 base を read-only import し、観測のみ追加（コード改変なし）。
乱数は正本と同一規則（40260722+1000N）。

判定（事前固定）:
 (V0a) 閉鎖保存: 全Nで max|Z·Z| < 1e-8 かつ ノルム誤差 < 1e-8。
 (V0b) 急拡大: T=30000 内に crossing（f>0.05）が検出される N を記録。
 (V0c) 三方向凝縮: crossing 後、rank_Q=4 が観測される。
 (V0d) 安定真空の最小N*: 後期窓（crossing+GUARD 以降）の rank_Q=4 持続率
       ≥0.9 かつ q2,q3,q4 が閾（Q_REL_TAU·q1）超を維持する最小N。
       作業記憶 N>6 と照合。
使い方: python3 run_v0_vacuum_genesis_sweep_v3.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
P7 = (HERE.parent / "第8論文_二段階seed除去による準安定相の因果分離"
      / "paper7_seedless_natural_figures3_4_v1")
sys.path.insert(0, str(P7))
import run_seedless_natural_figures3_4_v1 as base  # noqa: E402

N_LIST = (4, 5, 6, 7, 8, 9, 10, 12, 16)
T_MAX = 30000
SAMPLE = 50
F_CROSS = 0.05


def main():
    t0 = time.time()
    base.verify_sources()
    base.prepare_import_paths()
    from run_n300_dimension_saturation_v2 import dominant_plane, gram_reduce

    out = {"params": {"T": T_MAX, "sample": SAMPLE, "N_list": list(N_LIST)},
           "N": {}}
    print(f"{'N':>3} {'M':>4} {'crossing':>9} {'rank4初出':>9} "
          f"{'後期持続率':>9} {'q234保持':>8} {'閉鎖':>9}")
    for n in N_LIST:
        (sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp,
         parent_residual, parent_sigma, method) = base.build_seedless(n)
        crossing = None
        rank4_first = None
        samples = []
        max_zz = 0.0; max_ne = 0.0
        for t in range(T_MAX + 1):
            Z_perp = Z - p * (p @ Z) - q * (q @ Z)
            totZ = float(np.real(np.conj(Z) @ Z))
            f = float(np.real(np.conj(Z_perp) @ Z_perp)) / totZ
            if crossing is None and f > F_CROSS:
                crossing = t
            if t % SAMPLE == 0:
                gr = gram_reduce(sys_lr, Z)
                _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
                qs = base.qsv4(B0, Bdom)
                rank_q = int(np.sum(qs > base.Q_REL_TAU * qs[0]))
                if rank4_first is None and rank_q >= 4:
                    rank4_first = t
                samples.append((t, f, rank_q,
                                [float(x) for x in qs[:4]]))
                max_zz = max(max_zz, abs(complex(Z @ Z)))
                max_ne = max(max_ne, abs(totZ - 1.0))
            if t == T_MAX:
                break
            sys_lr.set_theta(np.angle(Z))
            sigma_estimate, wp = sys_lr.sigma_max_power(wp)
            Z = sys_lr.cayley_step(Z, sigma_estimate)
        # 後期窓の安定性
        guard = getattr(base, "GUARD", 2000)
        if crossing is not None:
            late = [s for s in samples if s[0] >= crossing + guard]
        else:
            late = []
        if late:
            pers = float(np.mean([1.0 if s[2] >= 4 else 0.0 for s in late]))
            q234 = float(np.mean([
                1.0 if all(qv > base.Q_REL_TAU * s[3][0] for qv in s[3][1:4])
                else 0.0 for s in late]))
        else:
            pers, q234 = float("nan"), float("nan")
        out["N"][n] = {"M": int(sys_lr.m), "crossing": crossing,
                       "rank4_first": rank4_first,
                       "late_rank4_persistence": pers,
                       "late_q234_retention": q234,
                       "max_zero_square": max_zz, "max_norm_err": max_ne,
                       "parent_residual": float(parent_residual)}
        print(f"{n:>3} {sys_lr.m:>4} {str(crossing):>9} {str(rank4_first):>9} "
              f"{pers if pers == pers else float('nan'):>9.3f} "
              f"{q234 if q234 == q234 else float('nan'):>8.3f} {max_zz:>9.2e}")

    v0a = all(out["N"][n]["max_zero_square"] < 1e-8
              and out["N"][n]["max_norm_err"] < 1e-8 for n in N_LIST)
    crossed = [n for n in N_LIST if out["N"][n]["crossing"] is not None]
    cond = [n for n in N_LIST if out["N"][n]["rank4_first"] is not None]
    stable = [n for n in N_LIST
              if out["N"][n]["late_rank4_persistence"] == out["N"][n]["late_rank4_persistence"]
              and out["N"][n]["late_rank4_persistence"] >= 0.9
              and out["N"][n]["late_q234_retention"] >= 0.9]
    nstar = None
    for n in sorted(stable):
        if all(m in stable for m in N_LIST if m >= n):
            nstar = n
            break
    print(f"\n(V0a) 閉鎖保存 <1e-8 全N: {'通過' if v0a else '不成立'}")
    print(f"(V0b) 急拡大が起きたN: {crossed}")
    print(f"(V0c) 三方向凝縮(rank_Q=4)が出たN: {cond}")
    print(f"(V0d) 安定真空（後期 rank4 持続≥0.9 かつ q234 保持≥0.9）のN: {stable}"
          f" → 最小N*={nstar}（作業記憶 N>6 と照合）")
    out.update({"V0a": bool(v0a), "crossed": crossed, "condensed": cond,
                "stable": stable, "N_star": nstar,
                "runtime_sec": time.time() - t0})
    (HERE / "result_v0_vacuum_genesis_sweep_v3.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
