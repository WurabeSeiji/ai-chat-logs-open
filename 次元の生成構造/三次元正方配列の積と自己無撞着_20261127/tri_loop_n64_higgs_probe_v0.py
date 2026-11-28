#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
512 マス複合状態 H = O8 ⊗ M8 (step 200) について、
  Part 1: ヒッグス場的特徴の検証（スピン・転置/巡回対称性・パウリ分解・各モードの応答固有値）
  Part 2: 近傍プローブとの相互作用比較
          - 光子的な波（ヌル/カイラル波 ΣU^2=0、大域直交複素波）
          - 質量を持ちそうな凝縮体（局在 8 マス回転塊、局在 32 マス四元数塊）
          それぞれを H に弱結合で重ね、単体走行 vs H 共存走行で
          エネルギー推移・回転角速度（減速/ロック）・交差結合 ρ_cross の差を厳密計測する。
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


def pauli_decompose_2x2(M):
    """2×2 行列 M を I, σ_x, iσ_y, σ_z に分解: M = c0*I + cx*σ_x + cy*(iσ_y) + cz*σ_z"""
    c0 = 0.5 * (M[0, 0] + M[1, 1])
    cz = 0.5 * (M[0, 0] - M[1, 1])
    cx = 0.5 * (M[0, 1] + M[1, 0])
    cy = 0.5 * (M[0, 1] - M[1, 0])  # iσ_y = [[0, 1], [-1, 0]]
    return c0, cx, cy, cz


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H4 = np.array(make_mult_table(4), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    # 基準となる 512 マス状態 H = U200 を生成
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, nrm_H = step_tri(U)
    H = U.copy()

    # 外側 O8 と 内側 M8 を SVD 抽出
    H6 = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    U_svd, S_svd, Vt_svd = np.linalg.svd(H6, full_matrices=False)
    outer = U_svd[:, 0].reshape(8, 8, 8)
    inner = Vt_svd[0, :].reshape(8, 8, 8)
    # 符号を O8 と揃える
    if np.sum(outer * O8) < 0:
        outer = -outer
        inner = -inner

    log("=" * 78)
    log("Part 1. 512 マス格子 H = O8 ⊗ M8 のスピン・対称性・パウリ基底・モード応答")
    log("=" * 78)

    # 1.1 転置対称性 (i <-> k) と巡回置換対称性 (i,j,k) -> (j,k,i)
    def sym_ratio(X, X_perm):
        s_part = 0.5 * (X + X_perm)
        a_part = 0.5 * (X - X_perm)
        tot = float(np.sum(np.abs(X) ** 2))
        return float(np.sum(np.abs(s_part) ** 2) / tot), float(np.sum(np.abs(a_part) ** 2) / tot)

    s_ik, a_ik = sym_ratio(H, H.transpose(2, 1, 0))
    s_ij, a_ij = sym_ratio(H, H.transpose(1, 0, 2))
    s_jk, a_jk = sym_ratio(H, H.transpose(0, 2, 1))
    ov_cyc = float(abs(np.vdot(H, H.transpose(1, 2, 0))))
    log("[1.1] H 全体の添字置換対称性（スピン/テンソル階数指標）:")
    log(f"  i<->k 転置: 対称={s_ik:.6f}, 反対称={a_ik:.6f}")
    log(f"  i<->j 転置: 対称={s_ij:.6f}, 反対称={a_ij:.6f}")
    log(f"  j<->k 転置: 対称={s_jk:.6f}, 反対称={a_jk:.6f}")
    log(f"  巡回置換 (i,j,k)->(j,k,i) 重なり |<H, cyc(H)>| = {ov_cyc:.6f}")

    # 外側 O8 と内側 M8 それぞれの i<->k 転置対称性
    so_ik, ao_ik = sym_ratio(outer, outer.transpose(2, 1, 0))
    si_ik, ai_ik = sym_ratio(inner, inner.transpose(2, 1, 0))
    log(f"  [分解] 外側 O8 の i<->k: 対称={so_ik:.6f}, 反対称={ao_ik:.6f}")
    log(f"  [分解] 内側 M8 の a<->c: 対称={si_ik:.6f}, 反対称={ai_ik:.6f}")

    # 1.2 内側 M8 のパウリ基底分解 (b=0, 1) と時間変化
    log("\n[1.2] 内側 M8 (a,c in {4,5}) のパウリ分解 M = c0*I + cx*σ_x + cy*(iσ_y) + cz*σ_z:")
    for b in [0, 1]:
        M = inner[np.ix_([4, 5], [b], [4, 5])][:, 0, :]
        c0, cx, cy, cz = pauli_decompose_2x2(M)
        log(f"  b={b}: c0(I,スカラー)={c0:+.8f}, cy(iσ_y,回転/擬スカラー)={cy:+.8f}, "
            f"cx(σ_x)={cx:+.2e}, cz(σ_z)={cz:+.2e} | |z|=sqrt(c0^2+cy^2)={np.hypot(c0,cy):.8f}")

    # 1.3 線形応答解析（各スピン/対称性モードへの微小摂動 δ = 1e-5 に対する 1 ステップ増幅率）
    log("\n[1.3] 各対称性モードへの微小摂動 (ε=1e-5) の 1〜5 ステップ増幅率 |δ(t)|/|δ(0)|:")
    log("  (<1: 復元/質量ギャップあり, =1: 平坦/ゴールドストーン, >1: 不安定/成長モード)")

    # 基準軌道
    H_ref = [H.copy()]
    Ut = H.copy()
    for _ in range(10):
        Ut, _ = step_tri(Ut)
        H_ref.append(Ut.copy())

    def eval_mode(label, D, eps=1e-5):
        D = D / np.linalg.norm(D)
        # H に直交する成分に射影
        D = D - np.vdot(H, D) * H
        nd = np.linalg.norm(D)
        if nd < 1e-12:
            log(f"  {label:42s}: H と平行（純粋スカラー伸縮＝規格化で自明消去）")
            return
        D = D / nd
        Up = H + eps * D
        Up /= np.linalg.norm(Up)
        # 基準からの射影距離を t=1,2,3,5,10 で計測
        ratios = []
        for t in range(1, 11):
            Up, _ = step_tri(Up)
            if t in (1, 2, 3, 5, 10):
                dt = proj_dist(Up, H_ref[t])
                ratios.append((t, dt / eps))
        r_str = ", ".join(f"t={t}:{r:7.3f}" for t, r in ratios)
        log(f"  {label:42s}: {r_str}")

    # Mode 1: 内側 b=0 と b=1 の相対振幅比 (z1/z0 の絶対値) を変える（動径・対称性保持）
    D_rel_amp = inner.copy()
    D_rel_amp[:, 0, :] *= 1.0
    D_rel_amp[:, 1, :] *= -1.0
    eval_mode("M1: 内側 相対振幅 (|z1|/|z0| 変更)", tensor_prod_3d(outer, D_rel_amp))

    # Mode 2: 内側 b=0,1 共通位相回転 (U(1) 全体位相)
    R90 = np.array([[0.0, -1.0], [1.0, 0.0]])
    D_u1 = np.zeros_like(inner)
    for b in [0, 1]:
        M = inner[np.ix_([4, 5], [b], [4, 5])][:, 0, :]
        M90 = R90 @ M
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                D_u1[a, b, c] = M90[ia, ic]
    eval_mode("M2: 内側 共通U(1)位相回転 (θ0,θ1同方向)", tensor_prod_3d(outer, D_u1))

    # Mode 3: 内側 相対位相回転 (θ1 - θ0 変更)
    D_rel_ph = D_u1.copy()
    D_rel_ph[:, 1, :] *= -1.0
    eval_mode("M3: 内側 相対位相 (θ1-θ0 変更)", tensor_prod_3d(outer, D_rel_ph))

    # Mode 4: 内側 鏡映/四重極モード (σ_x: 対称トレースレス＝スピン2的変形)
    D_sx = np.zeros_like(inner)
    D_sx[4, 0, 5] = 1.0; D_sx[5, 0, 4] = 1.0
    eval_mode("M4: 内側 σ_x 変形 (鏡映/スピン2的)", tensor_prod_3d(outer, D_sx))

    # Mode 5: 内側 σ_z 変形 (対角異方性＝スピン2的変形)
    D_sz = np.zeros_like(inner)
    D_sz[4, 0, 4] = 1.0; D_sz[5, 0, 5] = -1.0
    eval_mode("M5: 内側 σ_z 変形 (対角異方性/スピン2的)", tensor_prod_3d(outer, D_sz))

    # Mode 6: 内側 直交チャネル a,c in {6,7}, b in {0,1} (同じ j 軸を共有する別の複素波!)
    D_orth_ac = np.zeros_like(inner)
    for b in [0, 1]:
        M = inner[np.ix_([4, 5], [b], [4, 5])][:, 0, :]
        for ia, a in enumerate([6, 7]):
            for ic, c in enumerate([6, 7]):
                D_orth_ac[a, b, c] = M[ia, ic]
    eval_mode("M6: 内側 直交チャネル a,c={6,7} (j共有)", tensor_prod_3d(outer, D_orth_ac))

    # Mode 7: 内側 直交チャネル a,c in {4,5}, b in {2,3} (i,k 軸共有、別 b チャネル)
    D_orth_b = np.zeros_like(inner)
    for ib, b in enumerate([2, 3]):
        M = inner[np.ix_([4, 5], [ib], [4, 5])][:, 0, :]
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                D_orth_b[a, b, c] = M[ia, ic]
    eval_mode("M7: 内側 直交チャネル b={2,3} (i,k共有)", tensor_prod_3d(outer, D_orth_b))

    # Mode 8: 外側 O8 の 1 マス振幅ゆらぎ（空間局所的な格子歪み）
    D_outer_loc = np.zeros_like(outer)
    D_outer_loc[0, 0, 0] = 1.0
    eval_mode("M8: 外側 O8 局所振幅ゆらぎ (A=B=C=0)", tensor_prod_3d(D_outer_loc, inner))

    # Mode 9: 外側 O8 の非ゼロ外 (0 領域) への 1 マス摂動
    D_outer_zero = np.zeros_like(outer)
    D_outer_zero[0, 1, 0] = 1.0
    eval_mode("M9: 外側 O8 ゼロ領域への局所摂動", tensor_prod_3d(D_outer_zero, inner))

    # Mode 10: 内側 O8 全体の結晶モード（512 マスから元の O8 結晶へ戻る方向）
    eval_mode("M10: 内側 O8 結晶復帰モード (O8 ⊗ O8方向)", T0)

    # =========================================================================
    # Part 2: 近傍プローブとの相互作用比較（光子的な波 vs 質量を持ちそうな凝縮体）
    # =========================================================================
    log("\n" + "=" * 78)
    log("Part 2. ヒッグス場候補 H (= O8 ⊗ M8) と近傍プローブの相互作用比較")
    log("=" * 78)

    # まず各プローブ単体（背景 H なし）と、H + eps * Probe（弱結合 eps = 0.1）を比較する
    # 複素数で回すプローブもあるため complex128 に統一
    H_c = H.astype(np.complex128)
    mask_H = np.abs(H_c) > 1e-8

    # プローブ生成ヘルパー
    # (P1) 光子的な波 1: ヌル・カイラル波 (ΣU^2 = 0, パウリ射影 I + i*σ_y 型) を外側 O8 全体に広げた波
    #      a,c in {6,7}, b in {0,1} に配置（H の j 軸 b={0,1} を共有し、i,k 軸は直交！）
    P1_in = np.zeros((8, 8, 8), dtype=np.complex128)
    # パウリ射影 (I + σ_y) = [[1, -i], [i, 1]] -> 2乗和 1^2 + (-i)^2 + i^2 + 1^2 = 0 (完全ヌル！)
    proj_null = np.array([[1.0, -1.0j], [1.0j, 1.0]], dtype=np.complex128)
    for b in [0, 1]:
        for ia, a in enumerate([6, 7]):
            for ic, c in enumerate([6, 7]):
                P1_in[a, b, c] = proj_null[ia, ic]
    P1 = tensor_prod_3d(outer.astype(np.complex128), P1_in)
    P1 /= np.linalg.norm(P1)

    # (P2) 光子的な波 2: 実数・大域的な直交複素回転波 (H と同じ 2×2 回転行列だが a,c in {6,7}, j 軸共有)
    P2 = tensor_prod_3d(outer.astype(np.complex128), D_orth_ac.astype(np.complex128))
    P2 /= np.linalg.norm(P2)

    # (P3) 光子的な波 3: 実数・大域的な直交複素回転波 (a,c in {4,5} 共有, b in {2,3} 直交)
    P3 = tensor_prod_3d(outer.astype(np.complex128), D_orth_b.astype(np.complex128))
    P3 /= np.linalg.norm(P3)

    # (P3b) i,k 共有・大域ヌル/カイラル波 (O8 ⊗ (I+σ_y)[b=2,3])
    P3b_in = np.zeros((8, 8, 8), dtype=np.complex128)
    for b in [2, 3]:
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                P3b_in[a, b, c] = proj_null[ia, ic]
    P3b = tensor_prod_3d(outer.astype(np.complex128), P3b_in)
    P3b /= np.linalg.norm(P3b)

    # (P3c) i,k 共有・大域スピン2(鏡映 σ_x)波 (O8 ⊗ σ_x[b=2,3])
    P3c_in = np.zeros((8, 8, 8), dtype=np.complex128)
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    for b in [2, 3]:
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                P3c_in[a, b, c] = sx[ia, ic]
    P3c = tensor_prod_3d(outer.astype(np.complex128), P3c_in)
    P3c /= np.linalg.norm(P3c)

    # (P4) 質量を持ちそうな凝縮体 1: 局在した 8 マス実数回転塊（単一サイト ⊗ M8[b=2,3], i,k共有）
    out_loc1 = np.zeros((8, 8, 8), dtype=np.complex128)
    out_loc1[0, 0, 0] = 1.0
    P4 = tensor_prod_3d(out_loc1, D_orth_b.astype(np.complex128))
    P4 /= np.linalg.norm(P4)

    # (P5) 質量を持ちそうな凝縮体 2: 局在した四元数凝縮体（H4 ⊗ M8[b=2,3], i,k共有）
    out_H4 = np.zeros((8, 8, 8), dtype=np.complex128)
    out_H4[:4, :4, :4] = H4
    P5 = tensor_prod_3d(out_H4, D_orth_b.astype(np.complex128))
    P5 /= np.linalg.norm(P5)

    # (P6) 質量を持ちそうな凝縮体 3: 局在した四元数ヌル塊（H4 ⊗ (I+σ_y)[b=2,3], i,k共有）
    P6 = tensor_prod_3d(out_H4, P3b_in)
    P6 /= np.linalg.norm(P6)

    probes = [
        ("P1: 大域ヌル波 (j共有, i,k直交)", P1),
        ("P2: 大域複素波 (j共有, i,k直交)", P2),
        ("P3a: 大域複素回転波 (O8 ⊗ M8[b=2,3], i,k共有)", P3),
        ("P3b: 大域ヌル/カイラル波 (O8 ⊗ (I+σ_y)[b=2,3], i,k共有)", P3b),
        ("P3c: 大域スピン2(σ_x)波 (O8 ⊗ σ_x[b=2,3], i,k共有)", P3c),
        ("P4: 局在単一サイト凝縮体 (δ000 ⊗ M8[b=2,3], i,k共有)", P4),
        ("P5: 局在四元数凝縮体 (H4 ⊗ M8[b=2,3], i,k共有)", P5),
        ("P6: 局在四元数ヌル塊 (H4 ⊗ (I+σ_y)[b=2,3], i,k共有)", P6),
    ]

    def measure_cross_rho(U_cur, mask_A, mask_B):
        """
        成分 (i,j,k) in mask_B へ寄与する三角形 (p,q,r) のうち、
        3辺の中に mask_A の成分と mask_B の成分が両方含まれる「交差三角形」の
        最大視野角 ρ_cross と寄与エネルギー比を計測する。
        """
        coords_B = np.argwhere(mask_B)
        if len(coords_B) == 0:
            return 0.0, 0.0
        # 代表 4 点で高速計測
        max_rho = 0.0
        cross_contrib = 0.0
        total_contrib = 0.0
        absU = np.abs(U_cur)
        nz_all = absU > 1e-12
        idx_j_all = np.where(np.any(nz_all, axis=(0, 2)))[0]
        idx_ik_all = np.where(np.any(nz_all, axis=(1, 2)) & np.any(nz_all, axis=(0, 1)))[0]
        if len(idx_j_all) == 0 or len(idx_ik_all) == 0:
            return 0.0, 0.0
        for (i, j, k) in coords_B[:4]:
            a1 = absU[i][np.ix_(idx_j_all, idx_ik_all)]
            a2 = absU[:, j, :][np.ix_(idx_ik_all, idx_ik_all)]
            a3 = absU[:, :, k][np.ix_(idx_ik_all, idx_j_all)]
            mA1 = mask_A[i][np.ix_(idx_j_all, idx_ik_all)]
            mA2 = mask_A[:, j, :][np.ix_(idx_ik_all, idx_ik_all)]
            mA3 = mask_A[:, :, k][np.ix_(idx_ik_all, idx_j_all)]

            e1 = a1[:, :, None]
            e2 = a2[None, :, :]
            e3 = a3.transpose(1, 0)[:, None, :]
            prod = e1 * e2 * e3
            valid = prod > 1e-24
            if not np.any(valid):
                continue
            has_A = (mA1[:, :, None] | mA2[None, :, :] | mA3.transpose(1, 0)[:, None, :])
            total_contrib += float(np.sum(prod[valid]))
            cross_v = valid & has_A
            if np.any(cross_v):
                cross_contrib += float(np.sum(prod[cross_v]))
                e1_b, e2_b, e3_b = np.broadcast_arrays(e1, e2, e3)
                tmin = np.minimum(np.minimum(e1_b[cross_v], e2_b[cross_v]), e3_b[cross_v])
                tmax = np.maximum(np.maximum(e1_b[cross_v], e2_b[cross_v]), e3_b[cross_v])
                mr = float(np.max(tmin / tmax))
                if mr > max_rho:
                    max_rho = mr
        ratio = (cross_contrib / total_contrib) if total_contrib > 0 else 0.0
        return max_rho, ratio

    eps_mix = 0.1  # プローブ初期振幅（エネルギー比 ~1%）
    for name, Pr in probes:
        mask_P = np.abs(Pr) > 1e-8
        # 1) 単体走行 (20 step)
        Us = Pr.copy()
        solo_nrms = []
        solo_d1s = []
        for t in range(1, 21):
            Us_next, nrm_s = step_tri(Us)
            if nrm_s == 0.0:
                solo_nrms.append(0.0)
                solo_d1s.append(0.0)
                Us = Us_next
                break
            solo_nrms.append(nrm_s)
            solo_d1s.append(proj_dist(Us, Us_next))
            Us = Us_next
        solo_status = f"生存(||P||={solo_nrms[-1]:.5f}, d1={solo_d1s[-1]:.3f})" if solo_nrms[-1] > 0 else f"step {len(solo_nrms)} 消滅"

        # 2) H との共存走行 (50 step)
        Um = H_c + eps_mix * Pr
        Um /= np.linalg.norm(Um)
        eP0 = float(np.sum(np.abs(Um[mask_P]) ** 2))
        eH0 = float(np.sum(np.abs(Um[mask_H]) ** 2))
        rho_c0, cr_r0 = measure_cross_rho(Um, mask_H, mask_P)
        log(f"\n[{name}]")
        log(f"  単体走行(20step): {solo_status}")
        log(f"  H共存 t= 0: E_H={eH0:.5f}, E_P={eP0:.5f}, E_other={1-eH0-eP0:.2e}, max_ρ_cross={rho_c0:.4f}, cross寄与比={cr_r0:.3f}")
        for t in range(1, 51):
            Um_next, nrm_m = step_tri(Um)
            d1_P = proj_dist(Um[mask_P], Um_next[mask_P]) if np.linalg.norm(Um_next[mask_P]) > 1e-12 else 0.0
            d1_H = proj_dist(Um[mask_H], Um_next[mask_H]) if np.linalg.norm(Um_next[mask_H]) > 1e-12 else 0.0
            Um = Um_next
            if t in (1, 2, 5, 10, 20, 50):
                eH = float(np.sum(np.abs(Um[mask_H]) ** 2))
                eP = float(np.sum(np.abs(Um[mask_P]) ** 2))
                eO = max(0.0, 1.0 - eH - eP)
                nz_m = int(np.sum(np.abs(Um) > 1e-5))
                rho_c, cr_r = measure_cross_rho(Um, mask_H, mask_P)
                sum_u2_P = float(abs(np.sum(Um[mask_P] ** 2)) / eP) if eP > 1e-12 else 0.0
                log(f"  H共存 t={t:2d}: nz={nz_m:4d}, ||P||={nrm_m:.5f}, E_H={eH:.5f}, E_P={eP:.5f}, E_oth={eO:.2e}, "
                    f"d1_P={d1_P:.3f}, |ΣP²|/E_P={sum_u2_P:.3f}, ρ_cross={rho_c:.4f}(寄与{cr_r:.2f})")

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_higgs_probe_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
