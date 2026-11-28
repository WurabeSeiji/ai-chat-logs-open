#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 3 層のその先（第 4 相転移・空間次元分離の兆候があるか）を検証するため、
Case 1 (64 マス最大エンタングル相, ||P||=0.4817) を t=1,000 まで、
Case 2 (1,554 マス・22 要素の多成分最大エンタングル相, ||P||=0.4310) を t=1,500 まで
現在の系のまま長ステップ走行する。

乱数なし。途中で相転移（||P|| や nz や SVD ランクの急変）が起きた場合は即座に詳細記録する。
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


def svd_top(U, k=4):
    U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    sv = np.linalg.svd(U6, compute_uv=False)
    return sv[:k], int(np.sum(sv > 1e-4))


def long_run(label, U_init, max_steps, log):
    log("\n" + "=" * 78)
    log(f"[{label}] (max_steps={max_steps})")
    log("=" * 78)
    U = U_init / np.linalg.norm(U_init)
    prev_nrm = None
    prev_nz = None
    t_start = time.time()

    for t in range(0, max_steps + 1):
        if t == 0:
            P0 = tri_prod_fast(U)
            nrm = float(np.linalg.norm(P0))
            d1 = 0.0
        else:
            U_next, nrm = step_tri(U)
            if nrm == 0.0:
                log(f"  t={t:4d}: ゼロ消滅")
                break
            d1 = proj_dist(U, U_next)
            U = U_next

        # 記録タイミング：定期 + 急変検知
        is_periodic = (t in (0, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500))
        absU = np.abs(U)
        mask = absU > 1e-5
        nz = int(np.sum(mask))

        is_jump = False
        if prev_nrm is not None and abs(nrm - prev_nrm) > 1e-3:
            is_jump = True
        if prev_nz is not None and abs(nz - prev_nz) > max(16, int(0.1 * prev_nz)):
            is_jump = True

        if is_periodic or is_jump:
            di = int(np.sum(np.any(mask, axis=(1, 2))))
            dj = int(np.sum(np.any(mask, axis=(0, 2))))
            dk = int(np.sum(np.any(mask, axis=(0, 1))))
            ipr = float(np.sum(absU ** 4))
            n_eff = (1.0 / ipr) if ipr > 0 else 0.0
            sv, rnk = svd_top(U, 4)
            sv_str = ",".join(f"{v:.4f}" for v in sv[:3])
            tag = " [★急変]" if is_jump else ""
            log(f"  t={t:4d}: nz={nz:5d}, N_eff={n_eff:7.2f}, 軸=({di:2d},{dj:2d},{dk:2d}), "
                f"||P||={nrm:.6f}, d1={d1:.3f}, sv=({sv_str}, r={rnk}){tag}")

        prev_nrm = nrm
        prev_nz = nz

    # 最終状態の構造解析（外側軸・内側軸・連結成分・エネルギー階層）
    U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, Sv, Vtv = np.linalg.svd(U6, full_matrices=False)
    out1 = Uu[:, 0].reshape(8, 8, 8)
    in1 = Vtv[0, :].reshape(8, 8, 8)
    out_nz = int(np.sum(np.abs(out1) > 1e-5))
    in_nz = int(np.sum(np.abs(in1) > 1e-5))
    # 外側 8 軸のエネルギー分布（空間次元の分離兆候を見る）
    eo_i = np.sum(out1 ** 2, axis=(1, 2))
    eo_j = np.sum(out1 ** 2, axis=(0, 2))
    eo_k = np.sum(out1 ** 2, axis=(0, 1))
    eo_avg = (eo_i + eo_j + eo_k) / 3.0
    log(f"  [最終 t={max_steps}] 所要={time.time()-t_start:.1f}s, SVD上位4={np.round(Sv[:4], 6)}")
    log(f"    外側(第1特異) nz={out_nz}, 8軸エネルギー={np.round(eo_avg, 5).tolist()}")
    log(f"    内側(第1特異) nz={in_nz}")
    return U


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)

    H6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, _, Vtv = np.linalg.svd(H6, full_matrices=False)
    outer = Uu[:, 0].reshape(8, 8, 8)
    inner = Vtv[0, :].reshape(8, 8, 8)

    # Case 1: 外側 O8 に実数種 (0,1,0)+1e-4 -> 1,000 ステップ
    out_r = outer.copy()
    out_r[0, 1, 0] += 1e-4
    U_c1 = tensor_prod_3d(out_r, inner)
    long_run("Case 1 長距離: 外側 O8 実数種 (0,1,0)+1e-4", U_c1, 1000, log)

    # Case 2: 大域丸め誤差種 (~1e-16) -> 1,500 ステップ
    U_c2 = tensor_prod_3d(outer, inner)
    long_run("Case 2 長距離: 大域丸め誤差種 (~1e-16, SVD再構成)", U_c2, 1500, log)

    log(f"\n総所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_long_run_layer3_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
