#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""主張17 の導出 — 「全実数で厳密解は無い」と「N>3 なら何でもあり」の切り分け

木原の問い（2026-08-12）:
  (Q1) 線分の長さは全て実数なのに Σ x_n² = 0 は全実数では成立しないのでは？
  (Q2) 実質、N > 3 なら何でもありなのでは？

この2問は **符号の与え方** を固定するか否かで答えが逆になる。ここを分ける。

----------------------------------------------------------------------
検査A（Q1）: 符号なし（全項が正）の零閉鎖
----------------------------------------------------------------------
x_n が全て実数なら Σ x_n² = 0 ⟺ 全ての n で x_n = 0。
長さは非負なので、非自明な配置は存在しない。これは証明であって測定ではない。
→ 非自明な解を持たせるには、どこかを虚数と読むしかない。Q1 は正しい。

----------------------------------------------------------------------
検査B（Q2 の「何でもあり」側）: 符号を自由に選べる場合
----------------------------------------------------------------------
どの関係を虚数と読むかを自由に選べるなら、条件は
    Σ_e s_e d_e² = 0,  s_e ∈ {+1, −1}
という **スカラー1本の方程式** にすぎない。未知数 M = N(N−1)/2 本に対し
拘束は1本なので、一般の配置でも符号の選び方でほぼ常に満たせてしまう。
実際に部分和問題として解き、達成できる最小の |Σ s d²| を測る。
→ 符号が自由なら Q2 は正しい（条件が事実上効かない）。

----------------------------------------------------------------------
検査C（Q2 の「何でもありでない」側）: 符号を面次元で固定した場合
----------------------------------------------------------------------
主張4 の符号則（両端点を含む最小の面の次元 k で符号を決める）を課すと、
交代和は配置ごとに一意に決まり、選ぶ自由が無い。平行多面体では恒等的に
零になるが、一般の中心対称配置では零にならない。
→ 符号が幾何で固定されているなら Q2 は誤り。何でもありではない。

----------------------------------------------------------------------
検査D: 平行多面体における厳密な全実数恒等式
----------------------------------------------------------------------
平行多面体（頂点 v_s = (1/2) Σ s_i u_i, s ∈ {±1}^d）では

    Σ_{辺} d²  =  Σ_{主対角線} d²   （厳密、全次元 d で成立）

が成り立つ。左辺は 2^{d−1} Σ_i |u_i|²、右辺も 2^{d−1} Σ_i |u_i|² である。
d=2 はオイラーの四辺形定理そのもの。すなわち **主対角線だけを虚数と読めば**
辺と主対角線の範囲で Σ x² = 0 が厳密に成立する。木原の予想
「虚数になるのは対角線」は、この部分集合の上では正しい。
ただし d≥3 では辺と主対角線だけでは全 M 本を尽くさない（面対角線が残る）。
その残りをどう符号付けるかが主張4 の面次元則である。

使い方:
    python3 check_real_solutions_v1.py
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEED = 20260812
TRIALS = 50


def pair_lengths(V: np.ndarray) -> np.ndarray:
    """全 M = N(N-1)/2 対の長さ。"""
    N = V.shape[0]
    ia, ib = np.triu_indices(N, k=1)
    return np.linalg.norm(V[ia] - V[ib], axis=1)


def best_sign_split(w: np.ndarray, iters: int = 20000, rng=None) -> float:
    """min |Σ s_e w_e| を貪欲＋局所探索で求める（w = d²）。

    厳密な部分和最小化は NP 困難だが、ここで必要なのは「ほぼ零にできるか」
    であって最適性ではない。降順貪欲（Karmarkar–Karp 的な単純版）で初期解を
    作り、1点反転と2点交換で改善する。
    """
    rng = rng or np.random.default_rng(0)
    o = np.argsort(-w)
    s = np.zeros(len(w))
    tot = 0.0
    for i in o:                      # 降順貪欲: 現在の偏りを打ち消す側へ置く
        s[i] = -1.0 if tot > 0 else 1.0
        tot += s[i] * w[i]
    best = abs(tot)
    for _ in range(iters):           # 1点反転の局所探索
        i = rng.integers(len(w))
        cand = abs(tot - 2 * s[i] * w[i])
        if cand < best - 1e-18:
            tot -= 2 * s[i] * w[i]; s[i] = -s[i]; best = cand
    return float(best)


def parallelotope(rng, d: int):
    A = rng.normal(size=(d, d))
    while abs(np.linalg.det(A)) < 1e-3:
        A = rng.normal(size=(d, d))
    S = np.array(list(itertools.product([-1.0, 1.0], repeat=d)))
    return S @ A.T / 2.0, S, A


def check_D(rng, ds=(2, 3, 4, 5, 6)) -> list:
    """検査D: Σ_辺 d² − Σ_主対角線 d² が厳密に零か。"""
    rows = []
    for d in ds:
        rel = []
        for _ in range(TRIALS):
            V, S, A = parallelotope(rng, d)
            N = 2 ** d
            ia, ib = np.triu_indices(N, k=1)
            diff = np.abs(S[ia] - S[ib]).sum(axis=1) / 2   # 異なる座標の個数
            L2 = np.linalg.norm(V[ia] - V[ib], axis=1) ** 2
            edge = L2[diff == 1].sum()                     # 辺: 1座標だけ違う
            diag = L2[diff == d].sum()                     # 主対角線: 全座標違う
            rel.append(abs(edge - diag) / edge)
        rows.append({"d": d, "N": 2 ** d,
                     "n_edges": int(d * 2 ** (d - 1)),
                     "n_main_diagonals": int(2 ** (d - 1)),
                     "n_all_pairs": int(2 ** d * (2 ** d - 1) // 2),
                     "rel_residual_median": float(np.median(rel)),
                     "rel_residual_max": float(np.max(rel))})
    return rows


def check_B(rng, Ns=(4, 5, 6, 7, 8)) -> list:
    """検査B: 符号自由なら |Σ s d²| をどこまで零に近づけられるか。"""
    rows = []
    for N in Ns:
        r = []
        for _ in range(TRIALS):
            V = rng.normal(size=(N, 3))
            w = pair_lengths(V) ** 2
            r.append(best_sign_split(w, rng=rng) / w.sum())
        rows.append({"N": N, "M": N * (N - 1) // 2,
                     "rel_min_median": float(np.median(r)),
                     "rel_min_max": float(np.max(r))})
    return rows


def check_E(rng, ds=(2, 3, 4, 5)) -> list:
    """検査E: 面次元クラス別の Σd² と、どのクラスが虚になるかの内訳。

    クラス k = 両端点で符号が異なる座標の個数（= 両端点を含む最小の面の次元）。
    符号は (−1)^{k+1}（k 奇数が実、k 偶数が虚）。

    測ってみると S_k ∝ C(d−1, k−1) であり、したがって
        Σ_k (−1)^{k+1} S_k ∝ Σ_{j=0}^{d−1} (−1)^j C(d−1, j) = (1−1)^{d−1} = 0
    となる（d ≥ 2）。交代和が消えるのは二項定理そのものである。
    """
    rows = []
    for d in ds:
        A = rng.normal(size=(d, d))
        while abs(np.linalg.det(A)) < 1e-3:
            A = rng.normal(size=(d, d))
        S = np.array(list(itertools.product([-1.0, 1.0], repeat=d)))
        V = S @ A.T / 2.0
        N = 2 ** d
        ia, ib = np.triu_indices(N, k=1)
        kcls = (np.abs(S[ia] - S[ib]).sum(axis=1) / 2).astype(int)
        L2 = np.linalg.norm(V[ia] - V[ib], axis=1) ** 2
        cls = []
        for k in range(1, d + 1):
            cls.append({"k": k, "count": int((kcls == k).sum()),
                        "sum_d2": float(L2[kcls == k].sum()),
                        "sign": 1 if k % 2 else -1,
                        "role": "辺" if k == 1 else
                                ("主対角線" if k == d else f"{k}次元面の対角線")})
        alt = sum(c["sign"] * c["sum_d2"] for c in cls)
        n_real = sum(c["count"] for c in cls if c["sign"] > 0)
        n_imag = sum(c["count"] for c in cls if c["sign"] < 0)
        rows.append({"d": d, "N": N, "M": N * (N - 1) // 2,
                     "classes": cls, "alternating_sum": float(alt),
                     "n_real": n_real, "n_imag": n_imag})
    return rows


def classes_from_lengths(d_vec: np.ndarray, N: int, tol: float = 1e-7):
    """長さだけから面次元クラス k を復元する（主張4-a の手続きの実行）。

      長さ → 二重中心化 → 配置（合同・鏡映を除いて一意, 主張3）
           → 凸包 → 各対を含む最小の面の次元 k

    make_figs.py の face_dim_classes は配置 V を入力に取る。ここでは
    **長さしか与えない**ところから同じ分類に到達できるかを確かめる。
    鏡映は面構造を変えないので、主張3 の (1:2) の不定性は分類に効かない。
    """
    from scipy.spatial import ConvexHull
    ia, ib = np.triu_indices(N, k=1)
    D2 = np.zeros((N, N))
    D2[ia, ib] = d_vec ** 2
    D2[ib, ia] = d_vec ** 2
    J = np.eye(N) - np.ones((N, N)) / N
    B = -0.5 * J @ D2 @ J
    B = 0.5 * (B + B.T)
    lam, U = np.linalg.eigh(B)
    o = np.argsort(-lam); lam, U = lam[o], U[:, o]
    r = int((lam > 1e-9 * max(lam[0], 1e-300)).sum())
    V = U[:, :r] * np.sqrt(lam[:r])            # N×r の復元配置
    H = ConvexHull(V)
    A_, b_ = H.equations[:, :r], H.equations[:, r]
    on = np.abs(V @ A_.T + b_) < tol           # (N, facets)
    out = np.zeros(len(ia), dtype=int)
    for e, (i, j) in enumerate(zip(ia, ib)):
        F = np.flatnonzero(on[i] & on[j])
        out[e] = r if len(F) == 0 else r - np.linalg.matrix_rank(A_[F], tol=1e-8)
    return out, r


def check_F(rng, ds=(2, 3, 4)) -> list:
    """検査F: 長さだけから対角線（クラス k）を復元できるか。"""
    rows = []
    for d in ds:
        N = 2 ** d
        ok = 0
        for _ in range(10):
            V, S, A = parallelotope(rng, d)
            ia, ib = np.triu_indices(N, k=1)
            truth = (np.abs(S[ia] - S[ib]).sum(axis=1) / 2).astype(int)
            dv = np.linalg.norm(V[ia] - V[ib], axis=1)
            got, r = classes_from_lengths(dv, N)
            ok += int(r == d and np.array_equal(got, truth))
        rows.append({"d": d, "N": N, "trials": 10, "exact_recovery": ok})
    return rows


def main() -> None:
    rng = np.random.default_rng(SEED)
    print("=" * 72)
    print("検査A（符号なし・全実数）: Σ x² = 0 ⟺ 全 x = 0。非自明解は存在しない。")
    print("  → 証明であって測定ではない。虚数を入れない限り解は無い。")

    print("=" * 72)
    print("検査B（符号自由）: min |Σ s d²| / Σ d²  ＝ 事実上いくらでも零にできるか")
    print(f"{'N':>4} {'M':>5} {'相対最小 中央値':>18} {'最悪':>14}")
    rowsB = check_B(rng)
    for r in rowsB:
        print(f"{r['N']:>4} {r['M']:>5} {r['rel_min_median']:>18.3e} "
              f"{r['rel_min_max']:>14.3e}")

    print("=" * 72)
    print("検査D（平行多面体・全実数の厳密恒等式）:")
    print("  Σ_辺 d² − Σ_主対角線 d² = 0 か（相対残差）")
    print(f"{'d':>3} {'N':>5} {'辺':>5} {'主対角':>7} {'全対':>6} "
          f"{'相対残差 中央値':>18} {'最悪':>12}")
    rowsD = check_D(rng)
    for r in rowsD:
        print(f"{r['d']:>3} {r['N']:>5} {r['n_edges']:>5} "
              f"{r['n_main_diagonals']:>7} {r['n_all_pairs']:>6} "
              f"{r['rel_residual_median']:>18.3e} {r['rel_residual_max']:>12.3e}")

    print("=" * 72)
    print("検査E（面次元クラス別の内訳・どの関係が虚になるか）")
    rowsE = check_E(rng)
    for r in rowsE:
        print(f" d={r['d']} N={r['N']} M={r['M']}:")
        for c in r["classes"]:
            print(f"   k={c['k']} {c['role']:<14} 本数={c['count']:>4} "
                  f"Σd²={c['sum_d2']:>10.4f} 符号={'+' if c['sign']>0 else '-'} "
                  f"→ {'実' if c['sign']>0 else '虚'}")
        print(f"   交代和 = {r['alternating_sum']:.3e}   "
              f"実 {r['n_real']} 本 / 虚 {r['n_imag']} 本")

    print("=" * 72)
    print("検査F（長さだけから対角線を復元できるか・主張4-a の手続き）")
    rowsF = check_F(rng)
    for r in rowsF:
        print(f"  d={r['d']} N={r['N']}: 完全一致 {r['exact_recovery']}/{r['trials']} "
              "（長さ → 二重中心化 → 凸包 → 最小面次元）")

    out = {"seed": SEED, "trials": TRIALS,
           "F_class_recovery_from_lengths": rowsF,
           "A": "Σx²=0 with all real x has only the trivial solution x=0",
           "B_free_signs": rowsB, "D_parallelotope_identity": rowsD,
           "E_face_dimension_classes": rowsE}
    p = HERE / "check_real_solutions_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"-> {p.name}")


if __name__ == "__main__":
    main()
