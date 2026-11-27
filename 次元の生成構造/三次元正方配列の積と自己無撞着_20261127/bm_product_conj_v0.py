#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
共役を 1 つ含む三項積 (共役 BM 積):
    Prod_*(T)_{ijk} = sum_{t=0}^{N-1} T_{itk} * conj(T_{ijt}) * T_{tjk}

の解析的検討と厳密検証:
  1. 完全無名 (S_N 不変) 5 変数縮約方程式の検算と、a=0, b1=0 の厳密解 (ΣT^2=0 を含む) の導出。
  2. N = 4, 8 の共役 BM 積結晶解に 1 マス欠陥を入れたときの 200 ステップ発展
     (0 に全滅するか、それとも溶け残り・周期・準周期として生き残るか)。
  3. 巡回対称 (N = 4, 5, 8) の反対称波・2 次位相・エルミート型に対する共役 BM 積の振る舞い
     (生き残り・周期・準周期・ΣT^2 = 0)。

乱数なし。3 次斉次スケーリング分離を用い、||P|| < 1e-12 を厳密ゼロとして停止する。
"""
import math
import os
import numpy as np


def bm_prod_conj(T):
    return np.einsum("itk,ijt,tjk->ijk", T, np.conj(T), T)


def bm_circulant_conj(F):
    N = F.shape[0]
    F_next = np.zeros_like(F)
    for u in range(N):
        for v in range(N):
            s_val = 0j
            for s in range(N):
                s_val += F[s, v] * np.conj(F[u, s]) * F[(u - s) % N, (v - s) % N]
            F_next[u, v] = s_val
    return F_next


def proj_dist(A, B):
    na = np.linalg.norm(A)
    nb = np.linalg.norm(B)
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - ov * ov)))


def sum_z2_ratio(T):
    denom = float(np.sum(np.abs(T) ** 2))
    return 0.0 if denom == 0 else float(abs(np.sum(T ** 2)) / denom)


def build_sn_tensor(a, b1, b2, b3, c, N):
    T = np.full((N, N, N), c, dtype=np.complex128)
    for i in range(N):
        T[i, :, i] = b2
    for i in range(N):
        T[i, i, :] = b3
    for j in range(N):
        T[:, j, j] = b1
    for i in range(N):
        T[i, i, i] = a
    return T


def bm_5var_conj_analytic(a, b1, b2, b3, c, N):
    """5 変数 (a, b1, b2, b3, c) に対する共役 BM 積の厳密縮約方程式"""
    a_n = a * abs(a) ** 2 + (N - 1) * b1 * b2 * np.conj(b3)
    b1_n = a * abs(b1) ** 2 + b1 * np.conj(b2) * b3 + (N - 2) * b1 * abs(c) ** 2
    b2_n = a * abs(b2) ** 2 + np.conj(b1) * b2 * b3 + (N - 2) * b2 * abs(c) ** 2
    b3_n = np.conj(a) * (b3 ** 2) + b1 * b2 * np.conj(b3) + (N - 2) * np.conj(b3) * (c ** 2)
    c_n = b3 * np.conj(b2) * c + np.conj(b1) * b3 * c + b1 * b2 * np.conj(c) + (N - 3) * c * abs(c) ** 2
    return a_n, b1_n, b2_n, b3_n, c_n


def active_indices(U, tol=1e-5):
    N = U.shape[0]
    mask = np.abs(U) > tol
    return [m for m in range(N) if np.any(mask[m, :, :]) or np.any(mask[:, m, :]) or np.any(mask[:, :, m])]


def track_orbit(T0, fn, max_steps=200):
    nrm0 = float(np.linalg.norm(T0))
    if nrm0 < 1e-15:
        return "zero", 0, 0.0, 0.0, T0
    U = T0 / nrm0
    hist = [U]
    for step in range(1, max_steps + 1):
        P = fn(U)
        nrm = float(np.linalg.norm(P))
        if nrm < 1e-12:
            return f"zero at step {step}", step, 0.0, 0.0, U
        U_next = P / nrm
        for prev_i in range(len(hist) - 1, max(-1, len(hist) - 21), -1):
            if proj_dist(hist[prev_i], U_next) < 1e-9 and step >= 10:
                p = step - prev_i
                return f"periodic({p}) at step {step}", step, sum_z2_ratio(U_next), nrm, U_next
        hist.append(U_next)
        U = U_next
    min_rec = min(proj_dist(h, U) for h in hist[-51:-1])
    d_prev = proj_dist(hist[-2], U)
    return f"non_periodic (d_prev={d_prev:.2e}, min_rec50={min_rec:.2e})", max_steps, sum_z2_ratio(U), nrm, U


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("共役を 1 つ含む三項積 Prod_*(T)_{ijk} = sum_t T_{itk} conj(T_{ijt}) T_{tjk} の検証")
    out("=" * 72)

    # 1. 5 変数縮約方程式の厳密検算
    test_vals = (1 + 2j, -3 + 1j, 2 - 4j, -1 + 5j, 3 + 2j)
    for N in (4, 5):
        T = build_sn_tensor(*test_vals, N)
        R_full = bm_prod_conj(T)
        an = bm_5var_conj_analytic(*test_vals, N)
        R_exp = build_sn_tensor(*an, N)
        assert np.allclose(R_full, R_exp, atol=1e-10), f"N={N} 5変数検算失敗"
    out("[検算合格] 共役 BM 積の 5 変数縮約方程式は N = 4, 5 の全成分と厳密一致。\n")

    # 2. S_N 結晶解 (a=0, b1=0, b2=i*r, b3=i/r, c=1) の検証
    # a=0, b1=0 のとき:
    #   b2' = (N-2) b2 |c|^2 = (N-2) b2
    #   b3' = (N-2) conj(b3) c^2 = -(N-2) b3  (b3 が純虚数なので反転!)
    #   c'  = c (conj(b2)*b3 + (N-3)|c|^2)
    # 1 ステップ目: b2=i*r, b3=i/r => conj(b2)*b3 = (-i*r)(i/r) = +1 => c' = (N-2) c!
    # 2 ステップ目: b2=i*r, b3=-i/r => conj(b2)*b3 = (-i*r)(-i/r) = -1 => c' = (N-4) c!
    # おっ! N = 4 では 2 ステップ目で c' = 0 になるが、b2, b3 が実数なら固定点!
    # もし b2 = r * e^{iθ}, b3 = (1/r) * e^{iθ}, c = 1 とすると:
    #   conj(b2)*b3 = 1 -> c' = (N-2) c (実数!)
    #   b2' = (N-2) r e^{iθ}, b3' = (N-2) (1/r) e^{-iθ}
    # よって θ = 0 (実数) なら固定点だが ΣT^2 > 0。
    # では b2, b3, c すべてに位相 (あるいは c = e^{iφ}) を許すと:
    #   もし c = e^{iφ}, b2 = r e^{iα}, b3 = (1/r) e^{iβ}:
    #   conj(b2)*b3 = e^{i(β-α)}。
    #   1 ステップで b3 の位相は 2φ - β、b2 の位相は α、c の位相は φ + arg(e^{i(β-α)} + N - 3)。
    out("--- 1. 共役 BM 積での N = 4, 8 の対称初期値と 1 マス欠陥発展 ---")
    # まず N=4, 8 で (a=0, b1=0, b2=r e^{iα}, b3=(1/r) e^{iβ}, c=1) や
    # 実数結晶 (a=1, b1=b2=b3=1, c=-1) に 1 マス欠陥を入れたときの挙動を見る
    for N in (4, 8):
        out(f"\n[N = {N}]")
        # 結晶候補 A: 実数結晶 a=1, b1=b2=b3=1, c=-1 (純粋 BM 積でも共役 BM 積でも厳密固定点 λ=N)
        T_real = build_sn_tensor(1, 1, 1, 1, -1, N)
        assert np.allclose(bm_prod_conj(T_real), N * T_real, atol=1e-12)
        # 結晶候補 B: 複素固定点候補 (a=0, b1=0, b2=i*sqrt(-u), b3=i/sqrt(-u), c=1 -> ΣT^2=0!)
        u = 0.5 * (-(N - 2.0) + math.sqrt((N - 2.0) ** 2 - 4.0))
        T_sz0 = build_sn_tensor(0, 0, 1j * math.sqrt(-u), 1j / math.sqrt(-u), 1.0, N)
        st_b, _, sz2_b, _, _ = track_orbit(T_sz0, bm_prod_conj, max_steps=100)
        out(f"  欠陥なし T_sz0 (a=0,b1=0,b2=i√|u|,b3=i/√|u|,c=1, 初期ΣT^2=0): {st_b}, 最終ΣT^2={sz2_b:.3e}")

        for cname, Tbase in [("実数結晶 (1,1,1,1,-1)", T_real), ("複素結晶 T_sz0 (初期ΣT^2=0)", T_sz0)]:
            for dlabel, pos, eps in [
                ("対角(0,0,0)+iε", (0, 0, 0), 1e-4j),
                ("b1(1,0,0)+ε", (1, 0, 0), 1e-4),
                ("c(0,1,2)+iε", (0, 1, 2), 1e-4j),
            ]:
                Td = Tbase.copy()
                Td[pos] += eps
                st, _, sz2, _, Uend = track_orbit(Td, bm_prod_conj, max_steps=150)
                nz = int(np.sum(np.abs(Uend) > 1e-5))
                act = active_indices(Uend)
                out(f"  [{cname} + {dlabel:14s}] → {st}, Σz^2={sz2:.3e}, nz={nz}/{N**3}, 生存要素={len(act)}")

    # 3. 巡回対称 (N = 4, 5, 8) に対する共役 BM 積の振る舞い
    out("\n--- 2. 巡回対称パターンに対する共役 BM 積 (N = 4, 5, 8) ---")
    for N in (4, 5, 8):
        w = np.exp(2j * np.pi / N)
        patterns = []
        # 反対称波
        for a, b in [(1, 2), (1, 3), (2, 3)]:
            if b < N:
                F = np.array([[w ** ((a * u + b * v) % N) - w ** ((a * v + b * u) % N) for v in range(N)] for u in range(N)], dtype=np.complex128)
                patterns.append((f"反対称 ({a},{b})", F))
        # 2 次位相
        for p, q, r in [(0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 2, 1)]:
            F = np.array([[w ** ((p * u * u + q * v * v + r * u * v) % N) for v in range(N)] for u in range(N)], dtype=np.complex128)
            patterns.append((f"2次位相 ({p},{q},{r})", F))
        # R4 のみ ω^(uv)
        F_r4 = np.zeros((N, N), dtype=np.complex128)
        for u in range(N):
            for v in range(N):
                if u != 0 and v != 0 and u != v:
                    F_r4[u, v] = w ** ((u * v) % N)
        patterns.append(("R4のみ ω^(uv)", F_r4))
        # cos(u) + i sin(v)
        F_cs = np.array([[np.cos(2 * np.pi * u / N) + 1j * np.sin(2 * np.pi * v / N) for v in range(N)] for u in range(N)], dtype=np.complex128)
        patterns.append(("cos(u)+i sin(v)", F_cs))

        out(f"\n[N = {N}]")
        for pname, F0 in patterns:
            st, _, sz2, _, Uend = track_orbit(F0, bm_circulant_conj, max_steps=150)
            nz = int(np.sum(np.abs(Uend) > 1e-5))
            out(f"  {pname:18s} → {st}, ΣF^2={sz2:.3e}, nz={nz}/{N*N}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_conj_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
