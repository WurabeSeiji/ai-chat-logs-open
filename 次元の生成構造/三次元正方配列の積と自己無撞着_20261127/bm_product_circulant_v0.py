#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3: 巡回対称 (並進不変 T_{i+g, j+g, k+g} = T_{ijk} mod N) な 3 次元正方配列の
2 次元縮約方程式 F(u, v) = T_{0, u, v} に対する解析的初期値パターンの探索。

縮約方程式:
    F'(u, v) = sum_{s=0}^{N-1} F(s, v) * F(u, s) * F((u - s) mod N, (v - s) mod N)

乱数なし。以下を調べる:
  1. 縮約方程式と N^3 全成分 BM 積の厳密一致検算 (N = 4)。
  2. 解析的パターン (2 次位相 ω^{p u^2 + q v^2 + r u v}、1 次位相、対角/境界ゼロ型、
     反対称・エルミート型など) を網羅的に調べ、
     - 0 に消えないもの
     - 固定点 / 有限周期 / 準周期軌道 (回り続ける自己無撞着)
     - 結果としての Σ F(u, v)^2 = 0 (すなわち Σ T_{ijk}^2 = 0)
     が現れるかを特定する。
"""
import os
import numpy as np


def bm_circulant(F):
    """F'(u, v) = sum_s F(s, v) * F(u, s) * F((u-s)%N, (v-s)%N)"""
    N = F.shape[0]
    F_next = np.zeros_like(F)
    for u in range(N):
        for v in range(N):
            s_val = 0j
            for s in range(N):
                s_val += F[s, v] * F[u, s] * F[(u - s) % N, (v - s) % N]
            F_next[u, v] = s_val
    return F_next


def verify_circulant(N=4):
    """巡回対称な N^3 配列に対する全成分 BM 積と縮約式の厳密一致を整数で検算"""
    # 乱数を使わず、決定的な整数テーブル F[u, v] = 3*u - 5*v + 2*u*v + 1
    F = np.array([[(3 * u - 5 * v + 2 * u * v + 1) for v in range(N)] for u in range(N)], dtype=np.int64)
    T = np.zeros((N, N, N), dtype=np.int64)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                T[i, j, k] = F[(j - i) % N, (k - i) % N]
    P_full = np.einsum("itk,ijt,tjk->ijk", T, T, T)
    F_red = np.zeros((N, N), dtype=np.int64)
    for u in range(N):
        for v in range(N):
            s_val = 0
            for s in range(N):
                s_val += F[s, v] * F[u, s] * F[(u - s) % N, (v - s) % N]
            F_red[u, v] = s_val
    for i in range(N):
        for j in range(N):
            for k in range(N):
                assert P_full[i, j, k] == F_red[(j - i) % N, (k - i) % N]
    return True


def proj_dist(A, B):
    """スカラー倍 (位相含む) を除いた射影距離"""
    na = np.linalg.norm(A)
    nb = np.linalg.norm(B)
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - ov * ov)))


def sum_z2_ratio(F):
    denom = float(np.sum(np.abs(F) ** 2))
    if denom == 0:
        return 0.0
    return float(abs(np.sum(F ** 2)) / denom)


def classify_orbit(F0, max_steps=60):
    """
    初期配列 F0 から 3 次斉次スケーリング分離で最大 max_steps 回 BM 積を適用し、
    軌道を分類する:
      - zero: 0 に消滅
      - fixed: 固定点 (period 1)
      - periodic(p): 周期 p (p <= 20)
      - quasi_or_drift: 60 ステップ以内に厳密には戻らない (準周期 or ドリフト/カオス)
    """
    nrm0 = np.linalg.norm(F0)
    if nrm0 < 1e-15:
        return "zero", 0, 0.0, 0.0, 0.0
    U = F0 / nrm0
    hist = [U]
    for step in range(1, max_steps + 1):
        P = bm_circulant(U)
        nrm = np.linalg.norm(P)
        if nrm < 1e-25:
            return "zero", step, 0.0, 0.0, 0.0
        U_next = P / nrm
        # 過去の状態との一致判定
        for prev_i in range(max(0, len(hist) - 20), len(hist)):
            if proj_dist(hist[prev_i], U_next) < 1e-8:
                period = step - prev_i
                sz2 = sum_z2_ratio(U_next)
                d_init = proj_dist(F0, U_next)
                return f"periodic({period})" if period > 1 else "fixed", step, sz2, d_init, float(nrm)
        hist.append(U_next)
        U = U_next

    # 60 ステップ戻らなかった場合: 後半 20 ステップの Σz^2 と、直前ステップとの距離・最小再帰距離
    sz2_end = sum_z2_ratio(U)
    d_prev = proj_dist(hist[-2], U)
    min_rec = min(proj_dist(h, U) for h in hist[-25:-1])
    return "non_periodic", max_steps, sz2_end, min_rec, float(nrm)


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("Step 3: 巡回対称 (並進不変) 2 次元縮約方程式 F(u, v) の解析的パターン探索")
    out("=" * 72)

    assert verify_circulant(4)
    out("[検算合格] 巡回縮約方程式 F'(u,v) は N = 4 の全成分 BM 積と厳密一致。\n")

    # パターン 1: 2 次位相 F(u, v) = ω^{p u^2 + q v^2 + r u v} (乱数なし、全パラメータ網羅)
    out("--- 1. 2 次位相 F(u, v) = ω^{p u^2 + q v^2 + r u v} (N = 3, 4, 5, 8) ---")
    for N in (3, 4, 5, 8):
        w = np.exp(2j * np.pi / N)
        counts = {}
        examples_sz2_zero = []
        examples_nonper = []
        # N=8 は重いので p, q, r は対称性を活かして p,q in {0,1,2,4}, r in 0..N-1 に絞らず N<=5 は全数、N=8 は部分
        rng_pq = range(N) if N <= 5 else [0, 1, 2, 3, 4]
        for p in rng_pq:
            for q in rng_pq:
                for r in range(N):
                    F0 = np.array([[w ** ((p * u * u + q * v * v + r * u * v) % N) for v in range(N)] for u in range(N)], dtype=np.complex128)
                    kind, st, sz2, aux, nrm = classify_orbit(F0, max_steps=40)
                    counts[kind] = counts.get(kind, 0) + 1
                    if kind != "zero" and sz2 < 1e-6 and len(examples_sz2_zero) < 5:
                        examples_sz2_zero.append((p, q, r, kind, st, sz2, aux))
                    if kind == "non_periodic" and len(examples_nonper) < 5:
                        examples_nonper.append((p, q, r, sz2, aux))
        out(f"[N = {N}] 軌道分類の集計: {dict(sorted(counts.items()))}")
        if examples_sz2_zero:
            out("  Σz^2 = 0 を満たす生き残り例 (p, q, r, 種別, ステップ, |Σz^2|/Σ|z|^2, aux):")
            for ex in examples_sz2_zero:
                out(f"    p={ex[0]}, q={ex[1]}, r={ex[2]}: {ex[3]} (step={ex[4]}), Σz^2={ex[5]:.2e}, d_init={ex[6]:.3f}")
        if examples_nonper:
            out("  非周期 (40 step 以内に戻らない) 例:")
            for ex in examples_nonper:
                out(f"    p={ex[0]}, q={ex[1]}, r={ex[2]}: Σz^2={ex[3]:.3e}, min_rec={ex[4]:.3e}")

    # パターン 2: 対角ゼロ + 2 次位相 (u=0 or v=0 or u=v を 0 にしたもの、および反対称/エルミート型)
    out("\n--- 2. 対角/境界ゼロ・反対称を含む解析パターン (N = 4, 5, 8) ---")
    for N in (4, 5, 8):
        w = np.exp(2j * np.pi / N)
        patterns = []
        # (a) 3 つ異なる成分のみ 2 次位相 (u!=0, v!=0, u!=v)
        for r in range(1, N):
            F = np.zeros((N, N), dtype=np.complex128)
            for u in range(N):
                for v in range(N):
                    if u != 0 and v != 0 and u != v:
                        F[u, v] = w ** ((r * u * v) % N)
            patterns.append((f"R4のみ ω^({r}uv)", F))
        # (b) 反対称型 F(u,v) = ω^{r u v} - ω^{r v u} (r は 1 次含む: ω^{a u + b v} - ω^{a v + b u})
        for a in range(1, min(N, 4)):
            for b in range(a + 1, min(N, 5)):
                F = np.array([[w ** ((a * u + b * v) % N) - w ** ((a * v + b * u) % N) for v in range(N)] for u in range(N)], dtype=np.complex128)
                patterns.append((f"反対称 ω^({a}u+{b}v)-対称", F))
        # (c) 差の 1 変数 + 回転位相: F(u, v) = w^{u} * c_v (境界ゼロなど)
        for p in (1, 2):
            F = np.zeros((N, N), dtype=np.complex128)
            for u in range(N):
                for v in range(N):
                    if u != v:
                        F[u, v] = w ** ((p * u - v) % N)
            patterns.append((f"u!=v で ω^({p}u-v)", F))
        # (d) N 巡回の「実部+虚部ねじれ」: F(u,v) = cos(2πu/N) + i sin(2πv/N) (u!=0 or v!=0)
        F_cs = np.array([[np.cos(2 * np.pi * u / N) + 1j * np.sin(2 * np.pi * v / N) for v in range(N)] for u in range(N)], dtype=np.complex128)
        patterns.append(("cos(u)+i sin(v)", F_cs))
        F_cs2 = np.array([[np.cos(2 * np.pi * (u - v) / N) + 1j * np.sin(2 * np.pi * (u + v) / N) for v in range(N)] for u in range(N)], dtype=np.complex128)
        patterns.append(("cos(u-v)+i sin(u+v)", F_cs2))

        out(f"[N = {N}]")
        for pname, F0 in patterns:
            kind, st, sz2, aux, nrm = classify_orbit(F0, max_steps=60)
            if kind != "zero":
                out(f"  {pname:26s}: {kind:14s} (step={st:2d}), |Σz^2|/Σ|z|^2={sz2:.3e}, aux={aux:.3e}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_circulant_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
