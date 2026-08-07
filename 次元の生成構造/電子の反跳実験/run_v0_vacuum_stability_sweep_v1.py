#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V0: 真空の準備——安定した空間ができる最小Nの実験（電子の反跳実験・基盤）

問い（木原指定）: 粒子のない真空＝空間だけの状態を先に準備する。タネなしの
単一閉鎖波から、N がいくつ以上なら3方向が安定に発展し、場の読出しで
読み出せるか。

構成（正本 run_nbody_rank_saturation_preliminary_v1 の部品をコピー移植）:
  N体の全二体関係波 X_e（M=C(N,2)本）・閉鎖 ΣX_e²=R²・
  位相差正弦生成子 K・Cayley直交更新（集団時計144步/周）。
  真空状態＝一般位置の閉鎖波のみ（局在種・粒子なし）。

読出し（すべて状態からの関係量）:
  - 閉鎖誤差（保存則）
  - rank K と理論値 2·min(N,⌊M/2⌋) の一致（一般位置等号・正本定理2）
  - 空間3平面: K の上位3特異平面（回転面）を XYZ 候補と読む
  - 安定性（v2修正）: 自己無撞着力学——各ステップで現在の位相から K_t を
    読み直し、その K_t のCayley步で状態を回す（空間が自分を生成する力学）。
    後期窓の主要3平面が収束するか: F_late = 後期スナップショット対の
    部分空間忠実度の最小値（空間が形成され定常になるか）。
    （v1 は固定生成子で発展させており空間は構成的に不変＝測るべきものを
    測っていなかった——設計ミスとして記録）
  - 可読性: 特異値ギャップ g=(σ₆−σ₇)/σ₆（3平面と内部モードの分離）

判定（事前固定）:
 (V0a) 閉鎖保存: 全Nで閉鎖誤差 < 1e-10（正本許容）。
 (V0b) ランク則: rank K = 2·min(N,⌊M/2⌋) が全N・全試行で成立。
 (V0c) 三方向の存在: 回転平面数 r=rank/2 ≥ 3 となる最小N（理論では
       rank≥6 ⟺ min(N,⌊M/2⌋)≥3 ⟺ N≥4…だが安定性は別問題）を確認。
 (V0d) 安定空間の最小N: 忠実度の最小値 F_min(N) > 0.9 が試行中央値で
       成立する最小 N* を特定し、N>6 という作業記憶と比較（記録・判定は
       閾値でなく曲線の形の報告を主とする）。
使い方: python3 run_v0_vacuum_stability_sweep_v1.py
"""
from __future__ import annotations
import json, math, time
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TAU = 2.0 * math.pi
RADIUS2 = 1.0
SEED_AMP = 0.35
T_STEPS = 720
PERIOD = 144
SNAP = 36
N_LIST = tuple(range(3, 13))
TRIALS = 8
RNG0 = 20260807


def relation_pairs(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


def relation_adjacency(pairs):
    M = len(pairs)
    A = np.zeros((M, M))
    for a in range(M):
        for b in range(M):
            if a != b and len(set(pairs[a]) & set(pairs[b])) > 0:
                A[a, b] = 1.0
    return A


def initial_closed_state(M, rng):
    re = rng.normal(size=M); re /= np.linalg.norm(re)
    im = rng.normal(size=M); im -= float(np.dot(im, re)) * re
    im /= np.linalg.norm(im)
    return math.sqrt(RADIUS2 + SEED_AMP ** 2) * re + 1j * SEED_AMP * im


def gen_from_phases(ph, A):
    d = ph[None, :] - ph[:, None]
    raw = A * np.sin(d)
    return 0.5 * (raw - raw.T)


def cayley(K, period):
    nrm = float(np.linalg.norm(K, ord=2))
    if nrm <= 1e-14:
        return np.eye(K.shape[0])
    Kn = K / nrm
    c = math.tan((TAU / period) / 2.0)
    I = np.eye(K.shape[0])
    return np.linalg.solve(I - c * Kn, I + c * Kn)


def top_planes_projector(K, n_planes=3):
    _, s, Vh = np.linalg.svd(K)
    V = Vh[:2 * n_planes].T
    return V @ V.T, s


def main():
    t0 = time.time()
    out = {"N": {}, "params": {"T": T_STEPS, "period": PERIOD, "trials": TRIALS}}
    print(f"{'N':>3} {'M':>4} {'閉鎖err':>9} {'rank':>5} {'期待':>4} "
          f"{'r面':>3} {'gap':>6} {'F_min中央値':>10}")
    for N in N_LIST:
        pairs = relation_pairs(N)
        M = len(pairs)
        A = relation_adjacency(pairs)
        exp_rank = 2 * min(N, M // 2)
        cls, rks, gaps, fmins = [], [], [], []
        for tr in range(TRIALS):
            rng = np.random.default_rng(RNG0 + 1000 * N + tr)
            X0 = initial_closed_state(M, rng)
            K0 = gen_from_phases(np.angle(X0), A)
            P0, s0 = top_planes_projector(K0)
            scale = max(s0[0], 1e-300)
            rank = int(np.sum(s0 > 1e-10 * scale))
            g = float((s0[5] - s0[6]) / max(s0[5], 1e-300)) if len(s0) > 6 else 1.0
            # 自己無撞着発展（空間が自分を生成する力学）と後期収束
            X = X0.copy()
            worst_cl = 0.0
            projs = []
            for t in range(1, T_STEPS + 1):
                Kt = gen_from_phases(np.angle(X), A)
                X = cayley(Kt, PERIOD) @ X
                cr = np.sum(X.real ** 2 - X.imag ** 2) - RADIUS2
                ci = 2.0 * np.sum(X.real * X.imag)
                worst_cl = max(worst_cl, float(np.hypot(cr, ci)))
                if t % SNAP == 0 and t > T_STEPS // 2:
                    Kt2 = gen_from_phases(np.angle(X), A)
                    Pt, _ = top_planes_projector(Kt2)
                    projs.append(Pt)
            fmin = 1.0
            for i in range(len(projs) - 1):
                F = float(np.trace(projs[i] @ projs[-1])) / 6.0
                fmin = min(fmin, F)
            cls.append(worst_cl); rks.append(rank); gaps.append(g); fmins.append(fmin)
        med_f = float(np.median(fmins))
        out["N"][N] = {"M": M, "closure_worst": float(np.max(cls)),
                       "rank": rks, "expected_rank": exp_rank,
                       "gap_median": float(np.median(gaps)),
                       "Fmin_median": med_f,
                       "Fmin_all": [float(f) for f in fmins]}
        print(f"{N:>3} {M:>4} {np.max(cls):>9.2e} {int(np.median(rks)):>5} "
              f"{exp_rank:>4} {int(np.median(rks))//2:>3} "
              f"{np.median(gaps):>6.3f} {med_f:>10.4f}")

    # 判定
    v0a = all(out["N"][N]["closure_worst"] < 1e-10 for N in N_LIST)
    v0b = all(all(r == out["N"][N]["expected_rank"] for r in out["N"][N]["rank"])
              for N in N_LIST)
    n3 = min((N for N in N_LIST if out["N"][N]["expected_rank"] >= 6), default=None)
    stab = [(N, out["N"][N]["Fmin_median"]) for N in N_LIST]
    nstar = min((N for N, f in stab if f > 0.9
                 and all(f2 > 0.9 for N2, f2 in stab if N2 >= N)), default=None)
    print(f"\n(V0a) 閉鎖保存 <1e-10 全N: {'通過' if v0a else '不成立'}")
    print(f"(V0b) ランク則等号 全N全試行: {'通過' if v0b else '不成立'}")
    print(f"(V0c) 3回転平面の最小N: N={n3}")
    print(f"(V0d) 安定空間(後期収束) F_late>0.9(中央値・以後維持) の最小N*: {nstar}"
          f"（作業記憶 N>6 と比較）")
    out.update({"V0a": bool(v0a), "V0b": bool(v0b), "N_three_planes": n3,
                "N_star_stable": nstar, "runtime_sec": time.time() - t0})
    (HERE / "result_v0_vacuum_stability_sweep_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
