#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
課題 1: 八元数積表 (8^3) の虚数欠陥 +iε による 32 マス準安定構造 (四元数部分代数) の
20,000 ステップ長期追跡。

確認項目:
  - step 20,000 までに飽和・準安定化するか、8 マス固定点へ崩壊するか、
    あるいは 32 マス (四元数部分空間) 上で回り続けるか。
  - 各観測ステップでの nz (非ゼロ数)、||P||、dist(prev1)、dist(prev2)、
    |ΣU^2|/Σ|U|^2、非ゼロ添字集合 (I, J, K)。
  - 比較のため 8 マス固定点 (0,0,1)+iε と実数 8 マス回転 (0,1,0)+ε も 20,000 ステップ確認。

乱数なし。||P|| < 1e-12 で厳密ゼロ停止。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table


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


def idx_sets(U, tol=1e-5):
    mask = np.abs(U) > tol
    N = U.shape[0]
    I = [i for i in range(N) if np.any(mask[i, :, :])]
    J = [j for j in range(N) if np.any(mask[:, j, :])]
    K = [k for k in range(N) if np.any(mask[:, :, k])]
    return I, J, K, int(np.sum(mask))


def run():
    t0 = time.time()
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("課題 1: 32 マス準安定構造 (四元数部分代数) 等の 20,000 ステップ長期追跡")
    out("=" * 72)

    T0_int = make_mult_table(8)
    T0_cplx = np.array(T0_int, dtype=np.complex128)
    T0_real = np.array(T0_int, dtype=np.float64)

    checkpoints = {10, 20, 30, 40, 50, 60, 80, 100, 200, 500, 1000, 2000, 5000, 10000, 20000}

    cases = [
        ("32マス型1: (0,1,2)+iε [i=0, all_diff]", T0_cplx, (0, 1, 2), 1e-4j),
        ("32マス型2: (1,2,4)+iε [all>0, all_diff]", T0_cplx, (1, 2, 4), 1e-4j),
        ("32マス型3: (1,1,2)+iε [all>0, i=j]", T0_cplx, (1, 1, 2), 1e-4j),
        ("比較A: 8マス固定点 (0,0,1)+iε", T0_cplx, (0, 0, 1), 1e-4j),
        ("比較B: 実数8マス回転 (0,1,0)+ε", T0_real, (0, 1, 0), 1e-4),
    ]

    for label, Tbase, pos, eps in cases:
        out(f"\n--- [{label}] ---")
        U = Tbase.copy()
        U[pos] += eps
        U /= np.linalg.norm(U)
        U_prev1 = U.copy()
        U_prev2 = U.copy()
        settled_step = None

        for step in range(1, 20001):
            P = tri_prod_np(U)
            nrm = float(np.linalg.norm(P))
            if nrm < 1e-12:
                out(f"  step {step:5d}: ||P||={nrm:.3e} → 0 に消滅")
                break
            U_next = P / nrm
            d1 = proj_dist(U, U_next)
            d2 = proj_dist(U_prev1, U_next)

            if settled_step is None and step >= 100 and (d1 < 1e-11 or d2 < 1e-11):
                settled_step = step

            if step in checkpoints:
                I, J, K, nz = idx_sets(U_next)
                sz2 = sum_z2_ratio(U_next)
                out(f"  step {step:5d}: nz={nz:3d}/512, ||P||={nrm:.6f}, d1={d1:.2e}, d2={d2:.2e}, ΣU^2={sz2:.2e}, I={I}, J={J}, K={K}")

            U_prev2 = U_prev1
            U_prev1 = U
            U = U_next

        # step 20,000 付近の直近 20 ステップでの最小再帰距離と、飽和判定
        if nrm >= 1e-12:
            tail = [U]
            cur = U
            for _ in range(20):
                P = tri_prod_np(cur)
                cur = P / np.linalg.norm(P)
                tail.append(cur)
            min_rec20 = min(proj_dist(tail[0], tail[m]) for m in range(1, 21))
            d1_end = proj_dist(tail[0], tail[1])
            d2_end = proj_dist(tail[0], tail[2])
            out(f"  [20,000 step 後要約] 初回収束(d1/d2<1e-11)={settled_step}, 終端 d1={d1_end:.2e}, d2={d2_end:.2e}, min_rec20={min_rec20:.2e}")

    out(f"\n所要時間: {time.time() - t0:.2f} 秒")
    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_32cell_long_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
