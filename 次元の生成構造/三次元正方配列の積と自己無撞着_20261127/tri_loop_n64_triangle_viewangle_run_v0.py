#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N = 64 (O8 ⊗ O8) 実数 1 マス欠陥 (0,1,0)+1e-4 (step 0〜100) について、
三角形ループ Tri(U)_{ijk} = sum_{p,q,r} U_{ipq} U_{qjr} U_{rpk} を構成する
3 辺の強さ (a, b, c) = (|U_{ipq}|, |U_{qjr}|, |U_{rpk}|) の比 (視野角 ρ = L3 / L1、
指数部の差 Δln = ln(L1 / L3)) のスペクトルと、粒子の凝縮 (512 マスへのエネルギー集中)
を厳密に計測し、再現可能なデータ (.csv, .npz, .txt) として保存する。

乱数なし。
  - まず step 100 の最終状態から生き残り 512 マス集合 S_512 と、
    初期結晶 4,096 マスのうち溶け落ちた 3,584 マス集合 S_melt を厳密に特定する。
  - 粒子内部の代表出力位置 pos_part (S_512 内の最大振幅成分) と、
    溶け落ちた背景の代表出力位置 pos_melt (S_melt 内の初期非ゼロ成分 (0,0,0)) の
    それぞれについて、全 64^3 = 262,144 個の三角形ループ (p, q, r) のうち
    3 辺すべてが > 1e-300 となる三角形の 3 辺比 ρ = L3 / L1 を毎ステップ全数計測する。
"""
import csv
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d


def triangle_ratios_for_ijk(absU, i, j, k, eps_zero=1e-300):
    """
    出力成分 (i, j, k) に寄与する全 64^3 個の三角形 (p, q, r) の 3 辺
      a = |U[i, p, q]|,  b = |U[q, j, r]|,  c = |U[r, p, k]|
    から、3 辺とも > eps_zero の三角形について ρ = L3 / L1 (L1 >= L2 >= L3) を全数返す。
    """
    # shape: (p, q, r) = (64, 64, 64)
    A = absU[i, :, :][:, :, None]           # (p, q, 1)
    B = absU[:, j, :][None, :, :]           # (1, q, r)
    C = absU[:, :, k].T[:, None, :]         # (p, 1, r)  since absU[r, p, k] -> transpose -> (p, r)

    # 3 辺すべて非ゼロのマスク (ブロードキャスト)
    mask = (A > eps_zero) & (B > eps_zero) & (C > eps_zero)
    if not np.any(mask):
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    a_v = np.broadcast_to(A, mask.shape)[mask]
    b_v = np.broadcast_to(B, mask.shape)[mask]
    c_v = np.broadcast_to(C, mask.shape)[mask]

    L1 = np.maximum(np.maximum(a_v, b_v), c_v)
    L3 = np.minimum(np.minimum(a_v, b_v), c_v)
    rho = L3 / L1
    # 寄与の重み (3 辺の積 a*b*c) も返す
    weight = a_v * b_v * c_v
    return rho, weight


def summarize_rho(rho):
    if len(rho) == 0:
        return 0, np.nan, np.nan, np.nan, np.nan, np.nan
    log_r = np.log10(np.maximum(rho, 1e-300))
    dln = -np.log(np.maximum(rho, 1e-300))
    return (
        len(rho),
        float(np.median(rho)),
        float(np.percentile(rho, 10)),
        float(np.percentile(rho, 90)),
        float(np.median(dln)),
        float(np.max(dln)),
    )


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    mask_init = np.abs(T0) > 1e-5  # 4,096 マス

    # まず step 100 まで 1 回走らせて S_512 (生き残り 512 マス) と代表位置を厳密に特定
    U_tmp = T0.copy()
    U_tmp[0, 1, 0] += 1e-4
    U_tmp /= np.linalg.norm(U_tmp)
    for _ in range(100):
        P_tmp = tri_prod_fast(U_tmp)
        U_tmp = P_tmp / np.linalg.norm(P_tmp)
    mask_512 = np.abs(U_tmp) > 1e-5
    assert int(np.sum(mask_512)) == 512, f"step 100 の非ゼロ数が 512 でない: {int(np.sum(mask_512))}"
    mask_melt = mask_init & (~mask_512)
    # 注意: 512 マスのうち初期非ゼロと重なる個数も確認
    overlap_init_512 = int(np.sum(mask_init & mask_512))
    n_melt = int(np.sum(mask_melt))

    # 粒子代表: S_512 内で step 100 の振幅最大の位置
    pos_part = tuple(int(x) for x in np.unravel_index(np.argmax(np.abs(U_tmp)), U_tmp.shape))
    # 背景代表: 初期非ゼロ (0,0,0) が mask_melt にあることを確認
    assert mask_melt[0, 0, 0], "(0,0,0) が mask_melt に含まれない"
    pos_melt = (0, 0, 0)

    # step 0〜100 の本計測
    max_steps = 100
    snap_steps = [0, 20, 30, 35, 40, 60, 70, 80, 100]
    snap_set = set(snap_steps)

    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)

    rows = []
    snaps = {"snap_steps": np.array(snap_steps, dtype=np.int32)}

    for step in range(0, max_steps + 1):
        absU = np.abs(U)
        # 凝縮指標
        e_512 = float(np.sum(absU[mask_512] ** 2))
        e_melt = float(np.sum(absU[mask_melt] ** 2))
        mean_512 = float(np.mean(absU[mask_512]))
        mean_melt = float(np.mean(absU[mask_melt]))
        max_melt = float(np.max(absU[mask_melt]))
        contrast = mean_512 / max(mean_melt, 1e-300)
        dln_contrast = float(np.log(max(contrast, 1.0)))
        ipr = float(np.sum(absU ** 4))
        n_eff = 1.0 / ipr

        # 三角形視野角の全数計測 (粒子代表 & 背景代表)
        rho_p, w_p = triangle_ratios_for_ijk(absU, *pos_part)
        rho_m, w_m = triangle_ratios_for_ijk(absU, *pos_melt)

        # 振幅 > 1e-12 の実質的な三角形だけで見た統計も併記
        rho_p_sig, _ = triangle_ratios_for_ijk(absU, *pos_part, eps_zero=1e-12)
        rho_m_sig, _ = triangle_ratios_for_ijk(absU, *pos_melt, eps_zero=1e-12)

        np_cnt, np_med, np_p10, np_p90, np_dln_med, np_dln_max = summarize_rho(rho_p)
        nm_cnt, nm_med, nm_p10, nm_p90, nm_dln_med, nm_dln_max = summarize_rho(rho_m)
        _, np_sig_med, _, _, np_sig_dln, _ = summarize_rho(rho_p_sig)
        _, nm_sig_med, _, _, nm_sig_dln, _ = summarize_rho(rho_m_sig)

        rows.append({
            "step": step,
            "e_512": e_512,
            "e_melt": e_melt,
            "mean_512": mean_512,
            "mean_melt": mean_melt,
            "max_melt": max_melt,
            "contrast": contrast,
            "dln_contrast": dln_contrast,
            "n_eff": n_eff,
            "tri_part_cnt": np_cnt,
            "tri_part_rho_med": np_med,
            "tri_part_rho_p10": np_p10,
            "tri_part_rho_p90": np_p90,
            "tri_part_dln_med": np_dln_med,
            "tri_part_dln_max": np_dln_max,
            "tri_part_sig_rho_med": np_sig_med,
            "tri_part_sig_dln_med": np_sig_dln,
            "tri_melt_cnt": nm_cnt,
            "tri_melt_rho_med": nm_med,
            "tri_melt_rho_p10": nm_p10,
            "tri_melt_rho_p90": nm_p90,
            "tri_melt_dln_med": nm_dln_med,
            "tri_melt_dln_max": nm_dln_max,
            "tri_melt_sig_rho_med": nm_sig_med,
            "tri_melt_sig_dln_med": nm_sig_dln,
        })

        if step in snap_set:
            snaps[f"rho_part_{step}"] = rho_p
            snaps[f"rho_melt_{step}"] = rho_m

        if step < max_steps:
            P = tri_prod_fast(U)
            U = P / np.linalg.norm(P)

    # CSV 保存
    csv_path = os.path.join("results", "tri_loop_n64_viewangle_timeseries.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # NPZ 保存
    npz_path = os.path.join("results", "tri_loop_n64_viewangle_snapshots.npz")
    np.savez_compressed(npz_path, **snaps)

    # テキストサマリー保存
    lines = []
    lines.append("=" * 82)
    lines.append("N = 64 (O8 ⊗ O8) 三角形ループの視野角 (ρ = L3/L1) と粒子凝縮の厳密計測 (step 0〜100)")
    lines.append("=" * 82)
    lines.append(f"粒子代表位置 pos_part = {pos_part} (512マス内最大成分), 初期4096との重なり = {overlap_init_512}/512")
    lines.append(f"背景代表位置 pos_melt = {pos_melt} (初期4096から溶け落ちた {n_melt} マスの代表)")
    lines.append("")
    show = {0, 10, 20, 25, 30, 35, 40, 50, 60, 65, 70, 75, 80, 90, 100}
    for r in rows:
        if r["step"] in show:
            lines.append(
                f"step {r['step']:3d}: E512={r['e_512']:.4f}, ln(粒子/背景)={r['dln_contrast']:7.2f}, "
                f"粒子側Δln(中央/最大)={r['tri_part_dln_med']:6.2f}/{r['tri_part_dln_max']:7.2f} (n={r['tri_part_cnt']}), "
                f"背景側Δln(中央/最大)={r['tri_melt_dln_med']:6.2f}/{r['tri_melt_dln_max']:7.2f} (n={r['tri_melt_cnt']})"
            )
    lines.append(f"\n所要時間: {time.time() - t0:.2f} 秒")
    lines.append(f"保存ファイル: {csv_path}, {npz_path}")

    txt_path = os.path.join("results", "tri_loop_n64_triangle_viewangle_run_v0.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    run()
