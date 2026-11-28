#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N = 64 (O8 ⊗ O8) 第 3 層・第 4 層のさらに先（終着点・追加相転移の有無）を探るため、
更新ステップ τ = 0 -> 10,000 までフル 64×64×64 厳密計算でスイープする。
（注: 反復回数は時間 t ではなく更新ステップ τ を使用）

引数:
  --case 1 : Case 1 (外側 O8 実数種 (0,1,0)+1e-4)
  --case 2 : Case 2 (大域丸め誤差種 ~1e-16, SVD 再構成)
  --max-tau 10000
乱数なし。
"""
import argparse
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


def svd_top(U, k=6):
    U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    sv = np.linalg.svd(U6, compute_uv=False)
    return sv[:k], int(np.sum(sv > 1e-4))


def run_sweep(case_id, max_tau):
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

    if case_id == 1:
        label = "Case 1 (外側 O8 実数種 (0,1,0)+1e-4)"
        out_r = outer.copy()
        out_r[0, 1, 0] += 1e-4
        U = tensor_prod_3d(out_r, inner)
    elif case_id == 2:
        label = "Case 2 (大域丸め誤差種 ~1e-16, SVD再構成)"
        U = tensor_prod_3d(outer, inner)
    else:
        raise ValueError("case_id must be 1 or 2")

    U /= np.linalg.norm(U)

    log("=" * 80)
    log(f"[{label}] フル 64×64×64 厳密スイープ: τ = 0 -> {max_tau}")
    log("=" * 80)

    checkpoints = {
        0, 50, 100, 200, 500, 800, 900, 950, 1000,
        1500, 2000, 2500, 3000, 4000, 5000,
        6000, 7000, 8000, 9000, 10000
    }

    prev_nrm = None
    last_jump_log = -999
    jump_count = 0

    for tau in range(0, max_tau + 1):
        if tau == 0:
            P0 = tri_prod_fast(U)
            nrm = float(np.linalg.norm(P0))
            d1 = 0.0
        else:
            U_next, nrm = step_tri(U)
            if nrm == 0.0:
                log(f"  τ={tau:5d}: ゼロ消滅")
                break
            d1 = proj_dist(U, U_next)
            U = U_next

        is_cp = (tau in checkpoints)
        is_jump = (prev_nrm is not None and abs(nrm - prev_nrm) > 2e-3)

        # 急変ログが連続しすぎないよう、同一遷移内でも重要変化のみ記録（最初の50回は全記録）
        if is_cp or (is_jump and (jump_count < 60 or tau - last_jump_log >= 10)):
            absU = np.abs(U)
            mask = absU > 1e-5
            nz = int(np.sum(mask))
            di = int(np.sum(np.any(mask, axis=(1, 2))))
            dj = int(np.sum(np.any(mask, axis=(0, 2))))
            dk = int(np.sum(np.any(mask, axis=(0, 1))))
            ipr = float(np.sum(absU ** 4))
            n_eff = (1.0 / ipr) if ipr > 0 else 0.0
            sv, rnk = svd_top(U, 6)
            sv_str = ",".join(f"{v:.4f}" for v in sv[:4])
            tag = " [★急変]" if is_jump else ""
            if is_jump:
                jump_count += 1
                last_jump_log = tau
            log(f"  τ={tau:5d}: nz={nz:5d}, N_eff={n_eff:7.2f}, 軸=({di:2d},{dj:2d},{dk:2d}), "
                f"||P||={nrm:.6f}, d1={d1:.3f}, sv=({sv_str}; r={rnk}){tag}")

        prev_nrm = nrm

    # 最終状態の構造解剖
    U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, Sv, Vtv = np.linalg.svd(U6, full_matrices=False)
    out1 = Uu[:, 0].reshape(8, 8, 8)
    in1 = Vtv[0, :].reshape(8, 8, 8)
    eo_avg = (np.sum(out1 ** 2, axis=(1, 2)) + np.sum(out1 ** 2, axis=(0, 2)) + np.sum(out1 ** 2, axis=(0, 1))) / 3.0
    ei_avg = (np.sum(in1 ** 2, axis=(1, 2)) + np.sum(in1 ** 2, axis=(0, 2)) + np.sum(in1 ** 2, axis=(0, 1))) / 3.0
    # 全体の外側 8 軸エネルギー（ランク > 1 でも正確な全エネルギー射影）
    U_abs2 = np.abs(U) ** 2
    U_abs2_6 = U_abs2.reshape(8, 8, 8, 8, 8, 8)
    eA_i = np.sum(U_abs2_6, axis=(1, 2, 3, 4, 5))
    eA_j = np.sum(U_abs2_6, axis=(0, 1, 3, 4, 5))
    eA_k = np.sum(U_abs2_6, axis=(0, 1, 2, 3, 5))
    eA_tot = (eA_i + eA_j + eA_k) / 3.0

    log("\n" + "-" * 80)
    log(f"[最終状態 τ={max_tau} 解剖]")
    log(f"  SVD 上位 8 特異値: {np.round(Sv[:8], 6).tolist()}")
    log(f"  全体エネルギーの外側 8 軸 (A=0..7) 分布: {np.round(eA_tot, 5).tolist()}")
    log(f"  第1特異ベクトル外側 8 軸分布:            {np.round(eo_avg, 5).tolist()}")
    log(f"  第1特異ベクトル内側 8 軸分布:            {np.round(ei_avg, 5).tolist()}")
    log(f"  総所要時間: {time.time() - t0:.2f} 秒")

    out_path = os.path.join("results", f"tri_loop_n64_tau{max_tau}_case{case_id}.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"  保存先: {out_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", type=int, required=True)
    p.add_argument("--max-tau", type=int, default=10000)
    args = p.parse_args()
    run_sweep(args.case, args.max_tau)
