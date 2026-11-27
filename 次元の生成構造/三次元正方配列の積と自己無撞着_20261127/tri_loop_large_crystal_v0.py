#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
方向 A (課題 2): 結合した大結晶 N = 32 (H4 ⊗ O8) および N = 64 (O8 ⊗ O8) の
三角形ループ積による短めテスト走行 (最大 80 ステップ)。

解析的根拠:
  三角形ループ積 Tri(T)_{ijk} = sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk} は
  テンソル積に対して厳密に Tri(A ⊗ B) = Tri(A) ⊗ Tri(B) を満たす。
  したがって H4 (4^3, 周期2) と O8 (8^3, 周期2) のテンソル積:
    - N = 32: H4 ⊗ O8 (32^3 = 32,768 成分, 初期非ゼロ 1,024, 厳密周期 2)
    - N = 64: O8 ⊗ O8 (64^3 = 262,144 成分, 初期非ゼロ 4,096, 厳密周期 2)
  は独立ブロックではなく全体が結合した厳密周期 2 結晶となる。

高速化 (厳密同値):
  BLAS matmul 2 回で計算し、N = 8 で einsum 版と誤差 < 1e-13 をアサート検算する。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table


def tri_prod_fast(T):
    """
    Tri(T)_{ijk} = sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}
      Step 1: X[i, p, j, r] = sum_q T[i, p, q] * T[q, j, r]
      Step 2: R[i, j, k]    = sum_{r, p} X[i, j, r, p] * T[r, p, k]
    """
    N = T.shape[0]
    X = (T.reshape(N * N, N) @ T.reshape(N, N * N)).reshape(N, N, N, N)
    # X[i, p, j, r] -> transpose(0, 2, 3, 1) -> [i, j, r, p]
    X_ijrp = np.ascontiguousarray(X.transpose(0, 2, 3, 1)).reshape(N * N, N * N)
    R = (X_ijrp @ T.reshape(N * N, N)).reshape(N, N, N)
    return R


def tensor_prod_3d(A, B):
    """(A ⊗ B)_{(i1 i2),(j1 j2),(k1 k2)} = A_{i1,j1,k1} * B_{i2,j2,k2}"""
    m = A.shape[0]
    n = B.shape[0]
    C = np.einsum("ijk,abc->iajbkc", A, B).reshape(m * n, m * n, m * n)
    return C


def proj_dist(A, B):
    na = float(np.linalg.norm(A))
    nb = float(np.linalg.norm(B))
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - min(1.0, ov * ov))))


def sum_z2_ratio(T):
    denom = float(np.sum(np.abs(T) ** 2))
    return 0.0 if denom == 0 else float(abs(np.sum(T ** 2)) / denom)


def support_summary(U, tol=1e-5):
    mask = np.abs(U) > tol
    nz = int(np.sum(mask))
    ni = int(np.sum(np.any(mask, axis=(1, 2))))
    nj = int(np.sum(np.any(mask, axis=(0, 2))))
    nk = int(np.sum(np.any(mask, axis=(0, 1))))
    return nz, (ni, nj, nk)


def run_short(label, T0, defects, max_steps, out, is_cplx=False):
    t_start = time.time()
    U = T0.astype(np.complex128) if is_cplx else T0.astype(np.float64).copy()
    for pos, val in defects:
        U[pos] += val
    U /= np.linalg.norm(U)
    U_prev1 = U.copy()
    N = U.shape[0]
    out(f"\n  [{label}]")
    checkpoints = {10, 20, 30, 40, 50, 60, 80}
    for step in range(1, max_steps + 1):
        P = tri_prod_fast(U)
        nrm = float(np.linalg.norm(P))
        if nrm < 1e-12:
            out(f"    step {step:3d}: ||P||={nrm:.3e} → 0 に消滅")
            break
        U_next = P / nrm
        d1 = proj_dist(U, U_next)
        d2 = proj_dist(U_prev1, U_next)
        if step in checkpoints or step == max_steps:
            nz, dims = support_summary(U_next)
            sz2 = sum_z2_ratio(U_next)
            out(f"    step {step:3d}: nz={nz:6d}/{N**3}, 生存軸={dims}, ||P||={nrm:.6f}, d1={d1:.2e}, d2={d2:.2e}, ΣU^2={sz2:.2e}")
        U_prev1 = U
        U = U_next
    out(f"    (経過 {time.time() - t_start:.2f} 秒)")


def run():
    t0 = time.time()
    lines = []
    def out(s=""):
        print(s, flush=True)
        lines.append(s)

    out("=" * 72)
    out("方向 A: 結合大結晶 N = 32 (H4⊗O8) および N = 64 (O8⊗O8) の短めテスト走行")
    out("=" * 72)

    H4 = np.array(make_mult_table(4), dtype=np.float64)
    O8 = np.array(make_mult_table(8), dtype=np.float64)

    # 1. 高速版と einsum 版の一致検算 (N = 8)
    P_fast = tri_prod_fast(O8)
    P_ein = np.einsum("ipq,qjr,rpk->ijk", O8, O8, O8, optimize=True)
    assert np.allclose(P_fast, P_ein, atol=1e-13), "tri_prod_fast 検算失敗"
    out("[検算合格] BLAS 高速版は N = 8 で einsum 版と厳密一致 (誤差 < 1e-13)。")

    # 2. N = 32 (H4 ⊗ O8)
    out("\n--- 1. N = 32: H4 ⊗ O8 (32^3 = 32,768 成分, 初期非ゼロ 1,024) ---")
    T32 = tensor_prod_3d(H4, O8)
    t_b = time.time()
    T32_1 = tri_prod_fast(T32)
    T32_2 = tri_prod_fast(T32_1)
    dt32 = (time.time() - t_b) / 2.0
    assert np.allclose(T32_2, (16.0 * 1296.0) * T32, atol=1e-6), "N=32 周期2検算失敗"
    out(f"[検算合格] N = 32 厳密周期 2 (P2 = 20736 * P0), 1 ステップ約 {dt32*1000:.1f} ms")

    run_short("実数 1マス欠陥 (0,1,0)+1e-4 [O8内]", T32, [((0, 1, 0), 1e-4)], 80, out, is_cplx=False)
    run_short("実数 1マス欠陥 (0,8,0)+1e-4 [H4間]", T32, [((0, 8, 0), 1e-4)], 80, out, is_cplx=False)
    run_short("実数 2マス微差欠陥 (0,1,0)+1e-4, (8,9,8)+1.01e-4", T32, [((0, 1, 0), 1e-4), ((8, 9, 8), 1.01e-4)], 80, out, is_cplx=False)
    run_short("虚数 1マス欠陥 (0,1,2)+1e-4j [O8内]", T32, [((0, 1, 2), 1e-4j)], 80, out, is_cplx=True)
    run_short("虚数 1マス欠陥 (0,8,16)+1e-4j [H4間]", T32, [((0, 8, 16), 1e-4j)], 80, out, is_cplx=True)

    # 3. N = 64 (O8 ⊗ O8)
    out("\n--- 2. N = 64: O8 ⊗ O8 (64^3 = 262,144 成分, 初期非ゼロ 4,096) ---")
    T64 = tensor_prod_3d(O8, O8)
    t_b = time.time()
    T64_1 = tri_prod_fast(T64)
    T64_2 = tri_prod_fast(T64_1)
    dt64 = (time.time() - t_b) / 2.0
    assert np.allclose(T64_2, (1296.0 * 1296.0) * T64, atol=1e-4), "N=64 周期2検算失敗"
    out(f"[検算合格] N = 64 厳密周期 2 (P2 = 1679616 * P0), 1 ステップ約 {dt64*1000:.1f} ms")

    run_short("実数 1マス欠陥 (0,1,0)+1e-4 [内側O8]", T64, [((0, 1, 0), 1e-4)], 80, out, is_cplx=False)
    run_short("実数 1マス欠陥 (0,8,0)+1e-4 [外側O8]", T64, [((0, 8, 0), 1e-4)], 80, out, is_cplx=False)
    run_short("実数 2マス微差欠陥 (0,1,0)+1e-4, (8,9,8)+1.01e-4", T64, [((0, 1, 0), 1e-4), ((8, 9, 8), 1.01e-4)], 80, out, is_cplx=False)
    run_short("虚数 1マス欠陥 (0,1,2)+1e-4j [内側O8]", T64, [((0, 1, 2), 1e-4j)], 80, out, is_cplx=True)
    run_short("虚数 1マス欠陥 (0,8,16)+1e-4j [外側O8]", T64, [((0, 8, 16), 1e-4j)], 80, out, is_cplx=True)

    out(f"\n総所要時間: {time.time() - t0:.2f} 秒")
    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_large_crystal_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
