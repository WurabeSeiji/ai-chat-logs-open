#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
「倍音読み（Harmonic Spectrum Analysis）」検証スクリプト:
  - 非局在な波（背景 H、大域複素波 P3a、スピン1ゲージ波）
  - 局在凝縮体（単一サイト塊+衣 P4、四元数塊+衣 P5、Case 1 第3層）
について、
  Part 1: 空間（外側 O8 添字空間 Z_2^3）上のウォルシュ＝アダマール変換（厳密フーリエ倍音次数 k=0,1,2,3）
          および 1D FFT 空間波数スペクトルの比較。
  Part 2: 更新ステップ τ 方向（256 ステップ時系列）の FFT 時間倍音スペクトル
          （基本波 vs 高調波ピーク数・スペクトルエントロピー・奇数倍音構造）の比較。
乱数なし。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d


def step_tri(U):
    P = tri_prod_fast(U)
    nrm = float(np.linalg.norm(P))
    if nrm < 1e-15:
        return np.zeros_like(U), 0.0
    return P / nrm, nrm


def walsh_hadamard_matrix(n=8):
    """サイズ 8 (= 2^3) の正規化ウォルシュ＝アダマール行列 H8 と各行のハミング重み（倍音次数 k=0,1,2,3）"""
    H = np.array([[1.0]])
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    H /= np.sqrt(n)
    orders = np.array([bin(i).count("1") for i in range(n)], dtype=int)
    return H, orders


def spatial_harmonic_profile(T_part, H8, orders):
    """
    64×64×64 テンソルのうち着目する成分 T_part を外側 (A,B,C) × 内側 (a,b,c) に分解し、
    外側 3 軸 (A,B,C) にウォルシュ＝アダマール変換 (Z_2^3 上の厳密フーリエ変換) を適用して、
    倍音次数 K = order(kA) + order(kB) + order(kC) (0..9) ごとのエネルギー分布を返す。
    また、1 軸射影エネルギー E(A) (A=0..7) の WHT 倍音次数 k=0,1,2,3 も返す。
    """
    nrm2 = float(np.sum(np.abs(T_part) ** 2))
    if nrm2 < 1e-24:
        return np.zeros(10), np.zeros(4), 0.0
    T6 = T_part.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)  # (A,B,C,a,b,c)
    # 外側 3 軸に WHT を適用
    W = np.einsum("pA,qB,rC,ABCabc->pqrabc", H8, H8, H8, T6)
    e_pqr = np.sum(np.abs(W) ** 2, axis=(3, 4, 5)) / nrm2  # (8,8,8), sum=1
    by_K = np.zeros(10)
    for p in range(8):
        for q in range(8):
            for r in range(8):
                K = orders[p] + orders[q] + orders[r]
                by_K[K] += e_pqr[p, q, r]

    # 1 軸エネルギー分布 E(A) の WHT (k=0..3)
    eA = (
        np.sum(np.abs(T6) ** 2, axis=(1, 2, 3, 4, 5))
        + np.sum(np.abs(T6) ** 2, axis=(0, 2, 3, 4, 5))
        + np.sum(np.abs(T6) ** 2, axis=(0, 1, 3, 4, 5))
    ) / (3.0 * nrm2)
    wA = H8 @ eA
    pA = (wA ** 2) / np.sum(wA ** 2)
    by_k1 = np.zeros(4)
    for p in range(8):
        by_k1[orders[p]] += pA[p]

    # スペクトルエントロピー (by_K 上)
    nz = by_K[by_K > 1e-12]
    ent = float(-np.sum(nz * np.log(nz)))
    return by_K, by_k1, ent


def temporal_fft_summary(series):
    """複素または実数時系列 (長さ T=256) の FFT パワースペクトル要約"""
    x = np.asarray(series)
    x = x - np.mean(x)
    if np.linalg.norm(x) < 1e-12:
        return 0.0, 0, []
    spec = np.abs(np.fft.rfft(x.real)) ** 2
    if np.iscomplexobj(x):
        spec = np.abs(np.fft.fft(x)) ** 2
        spec = spec[: len(spec) // 2 + 1]
    tot = float(np.sum(spec))
    if tot < 1e-24:
        return 0.0, 0, []
    p = spec / tot
    nz = p[p > 1e-12]
    ent = float(-np.sum(nz * np.log(nz)))
    n_peaks_1pct = int(np.sum(p > 0.01))
    top_idx = np.argsort(p)[::-1][:6]
    top_list = [(int(k), float(p[k])) for k in top_idx if p[k] > 0.005]
    return ent, n_peaks_1pct, top_list


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H4 = np.array(make_mult_table(4), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    H8, orders = walsh_hadamard_matrix(8)

    # 第 1 層 H = U^(200)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()

    H6 = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    Uu, _, Vtv = np.linalg.svd(H6, full_matrices=False)
    outer = Uu[:, 0].reshape(8, 8, 8)
    inner = Vtv[0, :].reshape(8, 8, 8)
    if np.sum(outer * O8) < 0:
        outer = -outer
        inner = -inner

    # プローブ構成 (b={2,3} チャネル)
    D_in_b = np.zeros_like(inner)
    for ib, b in enumerate([2, 3]):
        M = inner[np.ix_([4, 5], [ib], [4, 5])][:, 0, :]
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                D_in_b[a, b, c] = M[ia, ic]

    # b={2,3} チャネルマスク (プローブ本体 + 衣をまとめて抽出するため)
    mask_b23 = np.zeros((64, 64, 64), dtype=bool)
    for B in range(8):
        for b in (2, 3):
            mask_b23[:, 8 * B + b, :] = True

    G_lx = np.zeros((8, 8)); G_lx[2, 3] = -1; G_lx[3, 2] = 1
    dO_lx = np.einsum("ia,ajk->ijk", G_lx, outer) + np.einsum("jb,ibk->ijk", G_lx, outer) + np.einsum("kc,ijc->ijk", G_lx, outer)
    out_loc1 = np.zeros((8, 8, 8)); out_loc1[0, 0, 0] = 1.0
    out_H4 = np.zeros((8, 8, 8)); out_H4[:4, :4, :4] = H4

    cases = [
        ("1. 背景 H (O8 ⊗ M8, 非局在時間結晶)", H, None),
        ("2. 大域複素波 P3a (O8 ⊗ M8[b=2,3], 非局在)", H + 0.1 * (tensor_prod_3d(outer, D_in_b) / np.linalg.norm(tensor_prod_3d(outer, D_in_b))), mask_b23),
        ("3. ゲージ/ベクトル波 (Lx O8 ⊗ M8[b=2,3], 非局在)", H + 0.1 * (tensor_prod_3d(dO_lx, D_in_b) / np.linalg.norm(tensor_prod_3d(dO_lx, D_in_b))), mask_b23),
        ("4. 局在四元数凝縮体+衣 P5 (H4 ⊗ M8[b=2,3])", H + 0.1 * (tensor_prod_3d(out_H4, D_in_b) / np.linalg.norm(tensor_prod_3d(out_H4, D_in_b))), mask_b23),
        ("5. 局在単一サイト凝縮体+衣 P4 (δ000 ⊗ M8[b=2,3])", H + 0.1 * (tensor_prod_3d(out_loc1, D_in_b) / np.linalg.norm(tensor_prod_3d(out_loc1, D_in_b))), mask_b23),
    ]

    log("=" * 80)
    log("Part 1. 空間倍音読み (外側 O8 = Z_2^3 上のウォルシュ＝アダマール倍音次数 K=0..9)")
    log("  (τ = 30 安定後の状態について、非局在波 vs 局在凝縮体+衣 の空間倍音分布を比較)")
    log("=" * 80)

    ts_dict = {}
    for name, U_init, extract_mask in cases:
        Uc = U_init / np.linalg.norm(U_init)
        # まず τ=30 まで緩和させて安定状態（衣の形成後）を作る
        for _ in range(30):
            Uc, _ = step_tri(Uc)
        T_target = Uc * extract_mask if extract_mask is not None else Uc
        by_K, by_k1, ent_sp = spatial_harmonic_profile(T_target, H8, orders)
        # 基本波(K<=2) vs 中倍音(K=3..5) vs 高倍音(K>=6)
        low_K = float(np.sum(by_K[:3]))
        mid_K = float(np.sum(by_K[3:6]))
        hi_K = float(np.sum(by_K[6:]))
        # 1軸エネルギー分布の DC(k=0) vs 倍音(k=1,2,3)
        dc_1 = by_k1[0]
        harm_1 = float(np.sum(by_k1[1:]))
        log(f"\n[{name}]")
        log(f"  3軸倍音次数 K=0..9 分布: {np.round(by_K, 4).tolist()}")
        log(f"  低次(K<=2)={low_K:.4f}, 中次(K=3..5)={mid_K:.4f}, 高次(K>=6)={hi_K:.4f} | 空間倍音エントロピー H_sp={ent_sp:.4f}")
        log(f"  1軸射影 E(A) の WHT: 基本DC(k=0)={dc_1:.4f}, 空間倍音(k=1,2,3)={np.round(by_k1[1:], 4).tolist()} (倍音合計={harm_1:.4f})")

        # Part 2 用に τ=30 -> 286 (256 ステップ) の代表成分時系列を採取
        # 最大振幅成分の座標を特定
        coord = tuple(int(x) for x in np.unravel_index(np.argmax(np.abs(T_target)), T_target.shape))
        series = np.zeros(256)
        Ut = Uc.copy()
        for t in range(256):
            series[t] = Ut[coord]
            Ut, _ = step_tri(Ut)
        ts_dict[name] = (coord, series)

    log("\n" + "=" * 80)
    log("Part 2. 更新ステップ τ 方向の時間倍音スペクトル (256 ステップ FFT)")
    log("=" * 80)
    for name, (coord, series) in ts_dict.items():
        ent_t, npk, tops = temporal_fft_summary(series)
        top_str = ", ".join(f"f={k}:{v:.3f}" for k, v in tops)
        log(f"  {name:48s}: H_tau={ent_t:.3f}, 有意倍音(>1%)={npk:2d}本 | 上位: {top_str}")

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_harmonic_localization_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
