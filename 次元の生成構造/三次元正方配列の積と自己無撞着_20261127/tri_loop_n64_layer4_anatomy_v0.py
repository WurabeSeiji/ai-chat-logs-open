#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Case 1 (外側 O8 実数種 (0,1,0)+1e-4) を τ = 0 -> 1,200 まで厳密計算し、
  1. 図 1 用の全ステップ時系列データを results/tri_loop_n64_case1_layer4.npz に保存
  2. 第 4 層 (τ = 1,000, ||P|| = 0.999533, sv = 0.5×4, d1 = 0) の
     4 つの特異ベクトル・64 マス成分・固定点構造を厳密解剖して
     results/tri_loop_n64_layer4_anatomy_v0.txt に保存する。
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

    out_r = outer.copy()
    out_r[0, 1, 0] += 1e-4
    U = tensor_prod_3d(out_r, inner)
    U /= np.linalg.norm(U)

    max_tau = 1200
    taus = []
    nzs = []
    neffs = []
    nrms = []
    d1s = []
    sv4s = []

    U_layer4 = None

    for tau in range(0, max_tau + 1):
        if tau == 0:
            P0 = tri_prod_fast(U)
            nrm = float(np.linalg.norm(P0))
            d1 = 0.0
        else:
            U_next, nrm = step_tri(U)
            d1 = proj_dist(U, U_next)
            U = U_next

        # 高密度サンプリング: 遷移期は毎ステップ、安定期は 5〜10 ステップおき
        keep = (
            tau <= 70
            or (850 <= tau <= 960)
            or (tau % 10 == 0)
        )
        if keep:
            absU = np.abs(U)
            nz = int(np.sum(absU > 1e-5))
            ipr = float(np.sum(absU ** 4))
            n_eff = (1.0 / ipr) if ipr > 0 else 0.0
            U6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
            sv = np.linalg.svd(U6, compute_uv=False)
            taus.append(tau)
            nzs.append(nz)
            neffs.append(n_eff)
            nrms.append(nrm)
            d1s.append(d1)
            sv4s.append(sv[:4].copy())

        if tau == 1000:
            U_layer4 = U.copy()

    np.savez_compressed(
        os.path.join("results", "tri_loop_n64_case1_layer4.npz"),
        tau=np.array(taus, dtype=np.int32),
        nz=np.array(nzs, dtype=np.int32),
        n_eff=np.array(neffs, dtype=np.float64),
        norm_P=np.array(nrms, dtype=np.float64),
        d1=np.array(d1s, dtype=np.float64),
        sv4=np.array(sv4s, dtype=np.float64),
    )

    # 第 4 層 (τ = 1000) の完全解剖
    log("=" * 78)
    log("第 4 層 (Case 1, τ = 1,000, ||P|| = 0.999533, d1 = 0) の厳密解剖")
    log("=" * 78)
    absU = np.abs(U_layer4)
    mask = absU > 1e-5
    idx_i = sorted(np.where(np.any(mask, axis=(1, 2)))[0].tolist())
    idx_j = sorted(np.where(np.any(mask, axis=(0, 2)))[0].tolist())
    idx_k = sorted(np.where(np.any(mask, axis=(0, 1)))[0].tolist())
    log(f"生存軸 i={idx_i}, j={idx_j}, k={idx_k}")

    # 4×4×4 の部分テンソルをそのまま取り出す
    sub = U_layer4[np.ix_(idx_i, idx_j, idx_k)]
    log(f"4×4×4 部分テンソルのノルム: {np.linalg.norm(sub):.10f}")
    # 振幅の大きい上位 16 成分を表示
    flat_idx = np.argsort(np.abs(sub).ravel())[::-1]
    log("\n4×4×4 部分テンソルの成分（絶対値降順・上位 20 個）:")
    for r in range(20):
        ii, jj, kk = np.unravel_index(flat_idx[r], (4, 4, 4))
        val = sub[ii, jj, kk]
        log(f"  [{r+1:2d}] sub[{ii},{jj},{kk}] (orig=({idx_i[ii]},{idx_j[jj]},{idx_k[kk]})) = {val:+.10f}  (|val|^2={val**2:.6f})")

    # 4×4×4 上での 3 軸対称性・直交性・クロネッカー構造チェック
    # j 固定スライス M_j[ii, kk] = sub[ii, jj, kk] の特異値と直交性を見る
    log("\nj 固定 4×4 スライス M_j (j=0..3) の特異値とフロベニウスノルム:")
    for jj in range(4):
        Mj = sub[:, jj, :]
        sv_j = np.linalg.svd(Mj, compute_uv=False)
        log(f"  j={jj} (orig {idx_j[jj]}): ||M_j||={np.linalg.norm(Mj):.6f}, sv={np.round(sv_j, 6).tolist()}")
        # M_j の行列そのものを表示
        for row in Mj:
            log("      " + " ".join(f"{x:+.5f}" for x in row))

    # SVD (外側 2×2×2 × 内側 2×2×2) の 4 つの特異ベクトルを確認
    sub6 = U_layer4.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    u4, s4, vt4 = np.linalg.svd(sub6, full_matrices=False)
    log(f"\nSVD 上位 8 特異値: {np.round(s4[:8], 6).tolist()}")

    # 4×4×4 テンソルをさらに基底変換（M_0..M_3 の対角化・特異値）で解析
    # 特に ||P|| = 0.999533 ≈ 1.0 は、クロネッカー単位テンソル (Tri(E)=E, ||E||=1) に近いかチェック！
    # 単位的 4×4 行列の直和・置換行列への近さ: 各 (i,j) に対する k の最大絶対値
    abs_sorted = np.sort(np.abs(sub).ravel())[::-1]
    log(f"\n上位 8 成分の二乗和: {np.sum(abs_sorted[:8]**2):.8f} (N_eff={1.0/np.sum(sub**4):.4f})")
    log(f"上位 8 成分の絶対値: {np.round(abs_sorted[:8], 6).tolist()}")
    log(f"9〜16 位の絶対値:    {np.round(abs_sorted[8:16], 6).tolist()}")

    # 4×4×4 上の Tri(sub) との一致度（厳密固定点チェック）
    P_sub = tri_prod_fast(sub)
    log(f"4×4×4 部分テンソル単体の ||Tri(sub)|| = {np.linalg.norm(P_sub):.10f}")
    log(f"proj_dist(sub, Tri(sub)) = {proj_dist(sub, P_sub):.3e}")
    # 固有基底（4×4 の直交変換 O で sub が対角・置換ブロック化されるか確認）
    # グラム行列 G_i = sum_{j,k} sub[i,j,k] sub[i',j,k]
    Gi = np.einsum("ijk,ljk->il", sub, sub)
    Gj = np.einsum("ijk,ilk->jl", sub, sub)
    Gk = np.einsum("ijk,ijl->kl", sub, sub)
    ei, Vi = np.linalg.eigh(Gi)
    ej, Vj = np.linalg.eigh(Gj)
    ek, Vk = np.linalg.eigh(Gk)
    log(f"\n4×4 グラム行列の固有値:")
    log(f"  G_i: {np.round(ei[::-1], 6).tolist()}")
    log(f"  G_j: {np.round(ej[::-1], 6).tolist()}")
    log(f"  G_k: {np.round(ek[::-1], 6).tolist()}")

    # Vj 基底で回転したときの各スライスを見る（直交変換 Q^T M Q）
    sub_rot = np.einsum("ia,jb,kc,ijk->abc", Vi, Vj, Vk, sub)
    sr_abs = np.sort(np.abs(sub_rot).ravel())[::-1]
    log(f"\n固有基底 (Vi, Vj, Vk) 変換後の成分絶対値 上位 12 個:")
    log(f"  {np.round(sr_abs[:12], 6).tolist()}")
    log(f"  上位 1 成分の二乗 = {sr_abs[0]**2:.8f}, 上位 2 成分の二乗和 = {np.sum(sr_abs[:2]**2):.8f}")
    log(f"  上位 4 成分の二乗和 = {np.sum(sr_abs[:4]**2):.8f}, 上位 8 成分の二乗和 = {np.sum(sr_abs[:8]**2):.8f}")

    # 同一の直交変換 Q (= Vi) で i, k を変換（i, k は同じ基底空間のはず）
    sub_ik = np.einsum("ia,jb,kc,ijk->abc", Vi, Vj, Vi, sub)
    sik_abs = np.sort(np.abs(sub_ik).ravel())[::-1]
    log(f"共通直交基底 (Vi, Vj, Vi) 変換後の上位 12 成分:")
    log(f"  {np.round(sik_abs[:12], 6).tolist()}")

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_layer4_anatomy_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
