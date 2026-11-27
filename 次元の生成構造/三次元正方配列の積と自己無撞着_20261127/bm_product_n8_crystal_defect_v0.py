#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
実験 I: N = 8 厳密結晶解 T* (a=0, b1=0, b2=i(sqrt(2)-1), b3=-i(sqrt(2)+1), c=1)
に対する欠陥発展 (200 ステップ) と「溶け残り」の構造検証。

性質:
  - Prod(T*, T*, T*) = 6 * T*
  - Σ (T*)^2 = 0
  - 非ゼロ成分数: 448 / 512 (8x7x2 + 8x7x6)

乱数なし。3 次斉次スケーリング分離を用い、||P|| < 1e-12 を厳密ゼロとして停止する。
"""
import math
import os
import numpy as np


def build_n8_crystal():
    N = 8
    b2 = 1j * (math.sqrt(2.0) - 1.0)
    b3 = -1j * (math.sqrt(2.0) + 1.0)
    c = 1.0 + 0j
    T = np.zeros((N, N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i == j == k or j == k:
                    T[i, j, k] = 0.0
                elif i == k:
                    T[i, j, k] = b2
                elif i == j:
                    T[i, j, k] = b3
                else:
                    T[i, j, k] = c
    return T


def bm_prod(T):
    return np.einsum("itk,ijt,tjk->ijk", T, T, T)


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


def active_indices(U, tol=1e-5):
    """どの要素番号 (0..N-1) が非ゼロ成分に関与しているかを調べる"""
    N = U.shape[0]
    mask = np.abs(U) > tol
    act = []
    for m in range(N):
        if np.any(mask[m, :, :]) or np.any(mask[:, m, :]) or np.any(mask[:, :, m]):
            act.append(m)
    return act


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("実験 I: N = 8 厳密結晶解 T* (ΣT^2=0, 448/512 非ゼロ) の欠陥発展と溶け残り")
    out("=" * 72)

    T0 = build_n8_crystal()
    P0 = bm_prod(T0)
    assert np.allclose(P0, 6.0 * T0, atol=1e-13), "Prod(T*) = 6 T* 検算失敗"
    assert abs(np.sum(T0 ** 2)) < 1e-12, "Σ(T*)^2 = 0 検算失敗"
    out("[検算合格] Prod(T*) = 6 T*,  Σ(T*)^2 = 0,  非ゼロ = 448/512.\n")

    T0_unit = T0 / np.linalg.norm(T0)

    defects = [
        ("対角 (0,0,0) に +iε", [((0, 0, 0), 1e-4j)]),
        ("対角 (0,0,0) に -iε", [((0, 0, 0), -1e-4j)]),
        ("対角 (0,0,0) に +ε", [((0, 0, 0), 1e-4)]),
        ("2対角 (0,0,0),(1,1,1) に +iε", [((0, 0, 0), 1e-4j), ((1, 1, 1), 1e-4j)]),
        ("3対角 (0,0,0),(1,1,1),(2,2,2) に +iε", [((0, 0, 0), 1e-4j), ((1, 1, 1), 1e-4j), ((2, 2, 2), 1e-4j)]),
        ("b1 (1,0,0) に +ε", [((1, 0, 0), 1e-4)]),
        ("b1 (1,0,0) に +iε", [((1, 0, 0), 1e-4j)]),
        ("b2 (0,1,0) に +ε", [((0, 1, 0), 1e-4)]),
        ("b2 (0,1,0) に +iε", [((0, 1, 0), 1e-4j)]),
        ("c (0,1,2) に +ε", [((0, 1, 2), 1e-4)]),
        ("c (0,1,2) に +iε", [((0, 1, 2), 1e-4j)]),
    ]

    for label, dlist in defects:
        out(f"--- [{label}] ---")
        U = T0.copy()
        for pos, val in dlist:
            U[pos] += val
        U = U / np.linalg.norm(U)
        hist = [U]
        status = "running"
        for step in range(1, 201):
            P = bm_prod(U)
            nrm = float(np.linalg.norm(P))
            if nrm < 1e-12:
                out(f"  step {step:3d}: ||P|| = {nrm:.3e} → 0 に消滅")
                status = f"zero at step {step}"
                break
            U_next = P / nrm
            d_crys = proj_dist(T0_unit, U_next)
            d_prev = proj_dist(U, U_next)
            sz2 = sum_z2_ratio(U_next)
            if step in (10, 20, 40, 80, 150, 200):
                nz_cnt = int(np.sum(np.abs(U_next) > 1e-5))
                out(f"  step {step:3d}: dist(T*)={d_crys:.4e}, dist(prev)={d_prev:.4e}, |Σz^2|/Σ|z|^2={sz2:.4e}, nz={nz_cnt}/512")
            found_p = None
            for prev_i in range(max(0, len(hist) - 12), len(hist)):
                if proj_dist(hist[prev_i], U_next) < 1e-10:
                    found_p = step - prev_i
                    break
            hist.append(U_next)
            U = U_next
            if found_p is not None and step >= 15:
                nz_cnt = int(np.sum(np.abs(U) > 1e-5))
                act = active_indices(U)
                out(f"  step {step:3d}: 周期 {found_p} に自発停止 (誤差<1e-10), dist(T*)={d_crys:.4f}, |Σz^2|/Σ|z|^2={sz2:.4f}, nz={nz_cnt}/512, 生存要素={act} ({len(act)}個)")
                # もし k 要素の部分結晶 (a=0, b1=0) なら nz = k(k-1)*2 + k(k-1)(k-2) = k^2(k-1)
                k = len(act)
                if nz_cnt == k * k * (k - 1) and k >= 2:
                    out(f"    → ちょうど k={k} 要素のサブ結晶 (非ゼロ数 k^2(k-1) = {k*k*(k-1)}) と一致！")
                status = f"converged period={found_p} at step {step}"
                break
        if status == "running":
            min_rec = min(proj_dist(h, U) for h in hist[-51:-1])
            act = active_indices(U)
            nz_cnt = int(np.sum(np.abs(U) > 1e-5))
            out(f"  step 200 未収束: dist(T*)={proj_dist(T0_unit, U):.4f}, min_rec_50={min_rec:.3e}, nz={nz_cnt}/512, 生存要素={act}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_n8_crystal_defect_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
