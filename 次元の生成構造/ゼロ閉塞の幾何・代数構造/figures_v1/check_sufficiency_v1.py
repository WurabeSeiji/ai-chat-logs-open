#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""主張16-b の導出 — 「N = 2^d」も「rank = d」も、単独では十分条件でない

主張4 は、零閉鎖する配置が平行多面体であり頂点数が N = 2^d に限られることを
述べる。ここで確かめるのは、その **逆は成り立たない** ことである。

測る量は本ノートの定義そのもの（§0 用語表の偏差 s）を、配置が実際に張る
次元 r（数値ランク）で書いたものである。配置 V（N×r、重心が原点）に対し

    T = Vᵀ V,    c = r / N,    s_i = sqrt( v_iᵀ T⁻¹ v_i / c )

全頂点が慣性楕円面上にあれば全ての i で s_i = 1、したがって max/min = 1。

----------------------------------------------------------------------
先に確認しておくべき自明な場合（この検証で分かったこと）
----------------------------------------------------------------------
r = N−1（重心を除いた最大次元）のとき、s_i は配置によらず恒等的に 1 である。
射影行列 H = V T⁻¹ Vᵀ は 1 の直交補空間への射影に一致するので

    H = I − (1/N) 𝟏𝟏ᵀ,   H_ii = (N−1)/N = c   （全ての i）

となり s_i = 1 が恒等的に成立する。すなわち **次元を落としていない配置では
「楕円面に乗る」という条件は内容を持たない。** 条件が意味を持つのは
r < N−1 のとき、つまり配置の次元が落ちているときだけである。
数値模型で意味があるのはまさにこの領域である（物質側は rank が N−1 から
落ちる）。以下の比較はこの点を踏まえて組んである。

----------------------------------------------------------------------
比較する4条件（いずれも d = 4、本ノートが採る次元）
----------------------------------------------------------------------
  (0) rank = N−1（次元を落とさない）        → 自明に 1（内容が空・上の証明）
  (1) N = 2^d かつ rank = d だが平行多面体でない → 乗らない
  (2) rank = d だけ（頂点数は 2 の冪でない）     → 乗らない
  (3) 平行多面体                              → 厳密に乗る（s ≡ 1）

(1) と (2) の対比が要点である。頂点数を 2^d に合わせても、次元を d に
合わせても、その両方を満たしても、楕円面には乗らない。乗るのは配置が
平行多面体そのものであるときだけである。

乱数配置は試行ごとに引き直し、TRIALS 回の中央値と最悪値を報告する。
seed は固定してあるので再現する。

使い方:
    python3 check_sufficiency_v1.py
出力:
    check_sufficiency_v1.json（数値そのもの）と標準出力の表
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEED = 20260812
TRIALS = 200
D = 4                      # 本ノートが採る平行多面体の次元
N_POW = 2 ** D             # = 16
N_NONPOW = 13              # 2 の冪でない頂点数（条件(2)用）


def deviation_ratio(V: np.ndarray) -> tuple[float, float]:
    """配置 V（N×r、重心は呼び出し側で除去済み）の偏差 s の max/min と CV。"""
    N = V.shape[0]
    # 数値ランクまで落とす（張っていない方向を T⁻¹ に入れない）
    U, sv, _ = np.linalg.svd(V, full_matrices=False)
    r = int((sv > 1e-10 * sv[0]).sum())
    W = U[:, :r] * sv[:r]                      # N×r の座標
    # macOS の BLAS は 16×15 のような形状で偽の FP フラグを立てることがある。
    # 結果は有限で正しい（下の assert で確認する）ので警告だけ止める。
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        T = W.T @ W
        s2 = np.einsum("ij,jk,ik->i", W, np.linalg.inv(T), W) / (r / N)
    assert np.isfinite(s2).all(), "偏差の計算が発散した"
    s = np.sqrt(np.maximum(s2, 0.0))
    return float(s.max() / s.min()), float(s.std() / s.mean())


def centered(X: np.ndarray) -> np.ndarray:
    return X - X.mean(axis=0)


def case_full_rank(rng) -> np.ndarray:
    """(0) 次元を落とさない一般配置（rank = N−1）。自明に乗る。"""
    return centered(rng.normal(size=(N_POW, N_POW - 1)))


def case_pow_and_rank(rng) -> np.ndarray:
    """(1) 頂点数 2^d かつ rank = d。ただし平行多面体ではない一般配置。"""
    return centered(rng.normal(size=(N_POW, D)))


def case_rank_only(rng) -> np.ndarray:
    """(2) 次元だけ d に合わせた一般配置。頂点数は 2 の冪でない。"""
    return centered(rng.normal(size=(N_NONPOW, D)))


def case_parallelotope(rng) -> np.ndarray:
    """(3) 平行多面体そのもの。v_s = (1/2) Σ_i s_i u_i、s ∈ {±1}^d。"""
    A = rng.normal(size=(D, D))
    while abs(np.linalg.det(A)) < 1e-3:
        A = rng.normal(size=(D, D))
    S = np.array([[1 if (k >> i) & 1 else -1 for i in range(D)]
                  for k in range(2 ** D)], dtype=float)
    return centered(S @ A.T / 2.0)


CASES = [
    ("(0) rank = N-1（次元を落とさない）", case_full_rank, N_POW, "N-1"),
    ("(1) N = 2^d かつ rank = d（平行多面体でない）", case_pow_and_rank,
     N_POW, str(D)),
    ("(2) rank = d だけ（頂点数が 2 の冪でない）", case_rank_only,
     N_NONPOW, str(D)),
    ("(3) 平行多面体（主張4 が選び出す配置）", case_parallelotope, N_POW, str(D)),
]


def sweep_N_at_fixed_rank(rng, rank: int, Ns) -> list:
    """rank を d に固定して N を振る。自明化する境界 N = rank+1 を示すため。

    rank = N−1 となる N（= rank+1、単体）でだけ max/min が恒等的に 1 になる。
    それより N が大きい領域では、頂点を増やすほど一つの楕円面に乗りにくくなる。
    「次元を d に落とせば楕円面に乗る」わけではないことがここで分かる。
    """
    rows = []
    for N in Ns:
        rs = [deviation_ratio(centered(rng.normal(size=(N, rank))))[0]
              for _ in range(TRIALS)]
        rows.append({"N": N, "rank": rank, "N_minus_1": N - 1,
                     "degenerate": bool(N - 1 == rank),
                     "maxmin_median": float(np.median(rs)),
                     "maxmin_worst": float(np.max(rs))})
    return rows


def sweep_parallelotope_d(rng, ds) -> list:
    """平行多面体 N=2^d が自明ケース（rank = N−1）に落ちるのは d=1 だけ。"""
    rows = []
    for d in ds:
        N = 2 ** d
        rs = []
        for _ in range(50):
            A = rng.normal(size=(d, d))
            while abs(np.linalg.det(A)) < 1e-3:
                A = rng.normal(size=(d, d))
            S = np.array([[1 if (k >> i) & 1 else -1 for i in range(d)]
                          for k in range(N)], dtype=float)
            rs.append(deviation_ratio(centered(S @ A.T / 2.0))[0])
        rows.append({"d": d, "N": N, "rank": d, "N_minus_1": N - 1,
                     "degenerate": bool(d == N - 1),
                     "maxmin_median": float(np.median(rs))})
    return rows


def main() -> None:
    rng = np.random.default_rng(SEED)
    print(f"d = {D} / 2^d = {N_POW} / 試行 {TRIALS} 回 / seed = {SEED}")
    print(f"{'条件':<46} {'N':>4} {'rank':>5} "
          f"{'max/min 中央値':>16} {'最悪':>10} {'CV 中央値':>11}")
    rows = []
    for label, gen, nn, rk in CASES:
        ratios, cvs = [], []
        for _ in range(TRIALS):
            r_, cv_ = deviation_ratio(gen(rng))
            ratios.append(r_); cvs.append(cv_)
        med = float(np.median(ratios)); worst = float(np.max(ratios))
        cvm = float(np.median(cvs))
        rows.append({"case": label, "N": nn, "rank": rk,
                     "maxmin_median": med, "maxmin_worst": worst,
                     "cv_median": cvm})
        print(f"{label:<46} {nn:>4} {rk:>5} {med:>16.10f} {worst:>10.3f} "
              f"{cvm:>11.2e}")
    print(f"\n[rank = {D} に固定して N を振る]  "
          f"自明化するのは N = rank+1 = {D + 1} のときだけである")
    print(f"{'N':>4} {'rank':>5} {'N-1':>5} {'max/min 中央値':>16} {'最悪':>10}  判定")
    sweep_N = sweep_N_at_fixed_rank(rng, D, [5, 6, 7, 8, 10, 13, 16, 20])
    for r in sweep_N:
        print(f"{r['N']:>4} {r['rank']:>5} {r['N_minus_1']:>5} "
              f"{r['maxmin_median']:>16.10f} {r['maxmin_worst']:>10.3f}  "
              + ("自明（rank = N-1・単体）" if r["degenerate"] else "内容あり"))

    print(f"\n[平行多面体 N = 2^d]  rank = d と N-1 = 2^d-1 が一致するのは d=1 のみ")
    print(f"{'d':>3} {'N=2^d':>6} {'rank':>5} {'N-1':>5} {'max/min 中央値':>16}  判定")
    sweep_d = sweep_parallelotope_d(rng, [1, 2, 3, 4, 5])
    for r in sweep_d:
        print(f"{r['d']:>3} {r['N']:>6} {r['rank']:>5} {r['N_minus_1']:>5} "
              f"{r['maxmin_median']:>16.10f}  "
              + ("自明" if r["degenerate"] else "内容あり"))

    out = {"d": D, "N_pow": N_POW, "N_nonpow": N_NONPOW,
           "sweep_N_at_fixed_rank": sweep_N,
           "sweep_parallelotope_d": sweep_d,
           "trials": TRIALS, "seed": SEED,
           "definition": "s_i = sqrt(v_i^T T^-1 v_i / (rank/N)); "
                         "max/min = 1 iff all vertices lie on the inertia ellipsoid",
           "note_full_rank": "rank = N-1 のとき H = I - (1/N)11^T より s_i = 1 が恒等的に成立。"
                             "楕円面条件は次元が落ちている配置でのみ内容を持つ。",
           "rows": rows}
    p = HERE / "check_sufficiency_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"-> {p.name}")


if __name__ == "__main__":
    main()
