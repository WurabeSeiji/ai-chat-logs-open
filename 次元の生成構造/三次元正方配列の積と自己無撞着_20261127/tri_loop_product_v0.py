#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三角形ループ積 (3 次元正方配列の 3 辺ループ縮約):
    Tri(T)_{ijk} = sum_{p, q, r = 0}^{N-1} T_{ipq} * T_{qjr} * T_{rpk}

を以下について厳密に検証する (乱数なし):
  1. 複素数 (2^3)・四元数 (4^3)・八元数 (8^3) の積表 f_{abc} そのものに対する
     厳密整数計算 (割り算なし・浮動小数なし) と、3 軸並べ替え 6 通りの比較。
  2. 完全無名 (S_N 不変) な 3 変数・5 変数での振る舞い確認。
  3. 巡回対称 (N = 4, 5, 8) の反対称波・2 次位相・R4 のみ・cos+i sin の 100 ステップ追跡
     (ゼロ消滅、固定点、周期軌道、準周期軌道の有無、ΣT^2 = 0)。
  4. 見つかった自己無撞着結晶に対する 1 マス欠陥発展 (ゼロ全滅するか、溶け残りが残るか)。
"""
import itertools
import os
import numpy as np

from bm_product_algebra_tables_v0 import make_mult_table, permute_axes, stats, check_proportional, check_sum_z2_pairs


def tri_prod_int(T):
    """
    厳密整数版 (Python 任意精度整数):
      A[i, p, j, r] = sum_q T[i][p][q] * T[q][j][r]
      R[i, j, k]    = sum_{p, r} A[i, p, j, r] * T[r][p][k]
    """
    N = len(T)
    # まず (i, p, j, r) を計算
    A = [[[[0] * N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for p in range(N):
            row = T[i][p]
            for q in range(N):
                v = row[q]
                if v == 0:
                    continue
                Tq = T[q]
                for j in range(N):
                    Tqj = Tq[j]
                    Aipj = A[i][p][j]
                    for r in range(N):
                        if Tqj[r] != 0:
                            Aipj[r] += v * Tqj[r]
    R = [[[0] * N for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            Rij = R[i][j]
            for p in range(N):
                Aipj = A[i][p][j]
                for r in range(N):
                    a_val = Aipj[r]
                    if a_val == 0:
                        continue
                    Trp = T[r][p]
                    for k in range(N):
                        if Trp[k] != 0:
                            Rij[k] += a_val * Trp[k]
    return R


def tri_prod_np(T):
    return np.einsum("ipq,qjr,rpk->ijk", T, T, T, optimize=True)


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


def track_orbit(T0, max_steps=100):
    nrm0 = float(np.linalg.norm(T0))
    if nrm0 < 1e-15:
        return "zero", 0, 0.0, 0.0, T0
    U = T0 / nrm0
    hist = [U]
    for step in range(1, max_steps + 1):
        P = tri_prod_np(U)
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
    min_rec = min(proj_dist(h, U) for h in hist[-41:-1])
    d_prev = proj_dist(hist[-2], U)
    return f"non_periodic (d_prev={d_prev:.2e}, min_rec40={min_rec:.2e})", max_steps, sum_z2_ratio(U), nrm, U


def build_circulant_3d(F):
    N = F.shape[0]
    idx = np.arange(N)
    u = (idx[None, :, None] - idx[:, None, None]) % N
    v = (idx[None, None, :] - idx[:, None, None]) % N
    return F[u, v]


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("三角形ループ積 Tri(T)_{ijk} = sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk} の厳密検証")
    out("=" * 72)

    # まず厳密整数版と numpy 版の一致検算 (N = 4)
    T4_chk = make_mult_table(4)
    R_int = tri_prod_int(T4_chk)
    R_np = tri_prod_np(np.array(T4_chk, dtype=np.float64))
    assert np.allclose(np.array(R_int, dtype=np.float64), R_np, atol=1e-12)
    out("[検算合格] 厳密整数版と NumPy 版が N = 4 で全成分一致。\n")

    # --- 1. 複素数・四元数・八元数の積表に対する厳密整数計算 ---
    out("--- 1. 複素数 (2^3)・四元数 (4^3)・八元数 (8^3) の積表 (厳密整数) ---")
    surviving_tables = []
    for name, N in [("複素数 (C)", 2), ("四元数 (H)", 4), ("八元数 (O)", 8)]:
        out(f"\n[{name}: N = {N}]")
        T0 = make_mult_table(N)
        for perm in itertools.permutations((0, 1, 2)):
            P0 = permute_axes(T0, perm)
            P1 = tri_prod_int(P0)
            nz1, max1, sum1, sq1, freq1 = stats(P1)
            prop1 = check_proportional(P0, P1)
            prop1_s = "なし" if prop1 is None else (f"{prop1[0]}" if prop1[1] == 1 else f"{prop1[0]}/{prop1[1]}")
            if nz1 == 0:
                out(f"  perm={perm}: step 1 で 0 に消滅")
                continue
            P2 = tri_prod_int(P1)
            nz2, max2, sum2, sq2, freq2 = stats(P2)
            prop2_0 = check_proportional(P0, P2)
            prop2_1 = check_proportional(P1, P2)
            rel2 = []
            if prop2_0 is not None:
                rel2.append(f"P2=({prop2_0[0]}/{prop2_0[1]})*P0")
            if prop2_1 is not None:
                rel2.append(f"P2=({prop2_1[0]}/{prop2_1[1]})*P1")
            rel2_s = ", ".join(rel2) if rel2 else "比例なし"
            tp, nzp, mre, mim = check_sum_z2_pairs(P1)
            out(f"  perm={perm}: step1 nz={nz1}/{N**3}, max={max1}, P1∝P0:{prop1_s}, 分布={freq1 if len(freq1)<=6 else f'{len(freq1)}種'}, 2成分Σz^2非ゼロ={nzp}/{tp}")
            out(f"             step2 nz={nz2}/{N**3}, max={max2}, {rel2_s}")
            if prop1 is not None or prop2_0 is not None or prop2_1 is not None:
                surviving_tables.append((name, N, perm, P0, P1))

    # --- 2. 巡回対称パターン (N = 4, 5, 8) の 100 ステップ追跡 ---
    out("\n--- 2. 巡回対称パターン (N = 4, 5, 8) の 100 ステップ追跡 ---")
    surviving_circ = []
    for N in (4, 5, 8):
        w = np.exp(2j * np.pi / N)
        patterns = []
        for a, b in [(1, 2), (1, 3), (2, 3)]:
            if b < N:
                F = np.array([[w ** ((a * u + b * v) % N) - w ** ((a * v + b * u) % N) for v in range(N)] for u in range(N)], dtype=np.complex128)
                patterns.append((f"反対称 ({a},{b})", F))
        for p, q, r in [(0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 2, 1)]:
            F = np.array([[w ** ((p * u * u + q * v * v + r * u * v) % N) for v in range(N)] for u in range(N)], dtype=np.complex128)
            patterns.append((f"2次位相 ({p},{q},{r})", F))
        F_r4 = np.zeros((N, N), dtype=np.complex128)
        for u in range(N):
            for v in range(N):
                if u != 0 and v != 0 and u != v:
                    F_r4[u, v] = w ** ((u * v) % N)
        patterns.append(("R4のみ ω^(uv)", F_r4))
        F_cs = np.array([[np.cos(2 * np.pi * u / N) + 1j * np.sin(2 * np.pi * v / N) for v in range(N)] for u in range(N)], dtype=np.complex128)
        patterns.append(("cos(u)+i sin(v)", F_cs))

        out(f"\n[N = {N}]")
        for pname, F0 in patterns:
            T0 = build_circulant_3d(F0)
            st, step_end, sz2, nrm, Uend = track_orbit(T0, max_steps=100)
            nz = int(np.sum(np.abs(Uend) > 1e-5))
            out(f"  {pname:18s} → {st}, ΣT^2={sz2:.3e}, nz={nz}/{N**3}")
            if not st.startswith("zero") and sz2 < 1e-6:
                surviving_circ.append((N, pname, Uend, st))

    # --- 3. 見つかった自己無撞着結晶 (積表で比例保存されたもの・ΣT^2=0 巡回結晶) の 1 マス欠陥発展 ---
    out("\n--- 3. 見つかった結晶に対する 1 マス欠陥発展 (ε = 1e-4, 100 ステップ) ---")
    # (a) 積表のうち P1 ∝ P0 となった代表例 (各 N で 1 つ)
    seen_N = set()
    for name, N, perm, P0, P1 in surviving_tables:
        if N in seen_N:
            continue
        seen_N.add(N)
        T0 = np.array(P0, dtype=np.complex128)
        T0_unit = T0 / np.linalg.norm(T0)
        out(f"\n[積表 {name} N={N}, perm={perm}]")
        for dlabel, pos, eps in [("対角(0,0,0)+ε", (0, 0, 0), 1e-4), ("対角(0,0,0)+iε", (0, 0, 0), 1e-4j), ("(0,1,0)+ε", (0, 1, 0), 1e-4)]:
            Td = T0.copy()
            Td[pos] += eps
            st, _, sz2, _, Uend = track_orbit(Td, max_steps=100)
            nz = int(np.sum(np.abs(Uend) > 1e-5))
            d0 = proj_dist(T0_unit, Uend)
            out(f"  {dlabel:16s} → {st}, dist(T0)={d0:.4f}, ΣT^2={sz2:.3e}, nz={nz}/{N**3}")

    # (b) ΣT^2 = 0 で生き残った巡回結晶 (最大 4 例)
    for N, pname, Ucry, orig_st in surviving_circ[:4]:
        out(f"\n[巡回結晶 N={N}, {pname} ({orig_st})]")
        for dlabel, pos, eps in [("対角(0,0,0)+ε", (0, 0, 0), 1e-4), ("(0,1,0)+iε", (0, 1, 0), 1e-4j)]:
            Td = Ucry.copy()
            Td[pos] += eps
            st, _, sz2, _, Uend = track_orbit(Td, max_steps=100)
            nz = int(np.sum(np.abs(Uend) > 1e-5))
            d0 = proj_dist(Ucry, Uend)
            out(f"  {dlabel:16s} → {st}, dist(T0)={d0:.4f}, ΣT^2={sz2:.3e}, nz={nz}/{N**3}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_product_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
