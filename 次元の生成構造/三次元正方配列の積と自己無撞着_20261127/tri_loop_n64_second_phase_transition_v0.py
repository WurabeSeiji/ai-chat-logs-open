#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
512 マス状態 H = O8 ⊗ M8 (step 200) は「内側だけが溶け、外側 O8 がまだ無傷の
不完全な結晶溶解（準安定・過冷却状態）」である可能性を検証する。

外側 O8 のゼロ領域に微小な種を入れ、第 2 相転移の全過程（step 0 -> 250）を追跡する:
  Case 1: 外側 O8 ゼロ領域 (0,1,0) に実数種 +1e-4
  Case 2: 実験 (E) と同等の SVD 再構成による大域丸め誤差種 (~1e-16)
  Case 3: 外側 O8 ゼロ領域 (0,1,0) に虚数種 +1e-4j
乱数なし。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d, proj_dist


def step_tri(U):
    P = tri_prod_fast(U)
    nrm = float(np.linalg.norm(P))
    if nrm < 1e-15:
        return np.zeros_like(U), 0.0
    return P / nrm, nrm


def svd_rank_info(U):
    U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    sv = np.linalg.svd(U6, compute_uv=False)
    return float(sv[0]), float(sv[1]), int(np.sum(sv > 1e-4))


def run_case(label, U_init, max_steps, log):
    log("\n" + "=" * 78)
    log(f"[{label}]")
    log("=" * 78)
    U = U_init / np.linalg.norm(U_init)
    checkpoints = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 250}
    U_prev = U.copy()
    for t in range(0, max_steps + 1):
        if t == 0:
            P0 = tri_prod_fast(U)
            nrm = float(np.linalg.norm(P0))
            d1 = 0.0
        else:
            U_next, nrm = step_tri(U)
            if nrm == 0.0:
                log(f"  t={t:3d}: ゼロ消滅")
                break
            d1 = proj_dist(U, U_next)
            U_prev = U
            U = U_next

        if t in checkpoints:
            absU = np.abs(U)
            mask = absU > 1e-5
            nz = int(np.sum(mask))
            di = int(np.sum(np.any(mask, axis=(1, 2))))
            dj = int(np.sum(np.any(mask, axis=(0, 2))))
            dk = int(np.sum(np.any(mask, axis=(0, 1))))
            ipr = float(np.sum(absU ** 4))
            n_eff = (1.0 / ipr) if ipr > 0 else 0.0
            sum_u2 = float(abs(np.sum(U ** 2)))
            s1, s2, rnk = svd_rank_info(U)
            log(f"  t={t:3d}: nz={nz:5d}, N_eff={n_eff:7.2f}, 軸=({di:2d},{dj:2d},{dk:2d}), "
                f"||P||={nrm:.6f}, d1={d1:.3f}, |ΣU²|={sum_u2:.3e}, SVD=(s1={s1:.4f}, s2={s2:.2e}, rank={rnk})")

    # 最終状態の構造詳細 (外側・内側の非ゼロ数と 2×2 構造)
    U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, Sv, Vtv = np.linalg.svd(U6, full_matrices=False)
    out_f = Uu[:, 0].reshape(8, 8, 8)
    in_f = Vtv[0, :].reshape(8, 8, 8)
    nz_out = int(np.sum(np.abs(out_f) > 1e-5))
    nz_in = int(np.sum(np.abs(in_f) > 1e-5))
    log(f"  [最終 t={max_steps} SVD第1成分] 外側 nz={nz_out}, 内側 nz={nz_in}, σ1={Sv[0]:.6f}, σ2={Sv[1]:.2e}")
    # 最終状態のエネルギー階層 (上位ユニーク値)
    a2 = np.abs(U) ** 2
    ne = (np.sum(a2, axis=(1, 2)) + np.sum(a2, axis=(0, 2)) + np.sum(a2, axis=(0, 1))) / 3.0
    nz_ne = np.sort(ne[ne > 1e-6])[::-1]
    uniq_e, cnt_e = np.unique(np.round(nz_ne, 5), return_counts=True)
    log(f"  [最終 生存要素数={len(nz_ne)}] エネルギー階層 (値:要素数): " +
        ", ".join(f"{v:.5f}:{c}" for v, c in zip(uniq_e[::-1], cnt_e[::-1])))


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()

    # 外側・内側抽出
    H6 = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, _, Vtv = np.linalg.svd(H6, full_matrices=False)
    outer = Uu[:, 0].reshape(8, 8, 8)
    inner = Vtv[0, :].reshape(8, 8, 8)

    # Case 1: 外側 O8 のゼロ領域 (0,1,0) に実数種 +1e-4
    out_r = outer.copy()
    out_r[0, 1, 0] += 1e-4
    U_c1 = tensor_prod_3d(out_r, inner)
    run_case("Case 1: 外側 O8 に実数種 (0,1,0)+1e-4 (厳密直積から開始)", U_c1, 250, log)

    # Case 2: 実験 (E) と同じ大域丸め誤差種 (~1e-16) の追跡
    U_c2 = tensor_prod_3d(outer, inner)
    run_case("Case 2: 大域丸め誤差種 (~1e-16, SVD再構成)", U_c2, 250, log)

    # Case 3: 外側 O8 のゼロ領域 (0,1,0) に虚数種 +1e-4j
    out_i = outer.astype(np.complex128).copy()
    out_i[0, 1, 0] += 1e-4j
    U_c3 = tensor_prod_3d(out_i, inner.astype(np.complex128))
    run_case("Case 3: 外側 O8 に虚数種 (0,1,0)+1e-4j (カスケード検証)", U_c3, 250, log)

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_second_phase_transition_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
