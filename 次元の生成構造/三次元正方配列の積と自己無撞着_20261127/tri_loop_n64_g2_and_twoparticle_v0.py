#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 3 速報のダメ押し検証（所要 ~25 秒）:
  Part 1: 外側 O8 の全 28 次元 SO(8) 回転生成子（虚数 7 次元 SO(7) = 21 次元 + 実数混合 7 次元）
          を第 1 層 H = O8 ⊗ M8 に作用させ、
          - 厳密に H を不変にする次元数（= G2 自己同型群 14 次元！）
          - H に直交し |λ| = 1.000（無質量ゴールドストーン/ゲージモード）となる次元数
          - さらに 1 つの虚数軸 (例 A=7) を固定する G2 の部分群（= SU(3) 8 グルーオン！）の次元数
          を厳密に SVD カウントする。
  Part 2: 局在単一サイト凝縮体 (P4) を 2 個、
          (a) 同一サイト (A=0) に同相/逆相で重ねた場合
          (b) 別サイト (A=0 と A=1) に同相/逆相で置いた場合
          (c) 異なるチャネル (b={2,3} と b={6,7}) に 2 個置いた場合
          の 30 ステップ相互作用（共存・反発・融合・衣の共有）を検証する。
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


def apply_outer_op(T8, G):
    return (
        np.einsum("ia,ajk->ijk", G, T8)
        + np.einsum("jb,ibk->ijk", G, T8)
        + np.einsum("kc,ijc->ijk", G, T8)
    )


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

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

    # =========================================================================
    # Part 1: SO(8) 28 次元 / SO(7) 21 次元 / SO(6) 15 次元の厳密カーネル次元カウント
    #   G2 = Ker(SO(7) -> dO8) = 14 次元
    #   SU(3) = Ker(SO(6) -> dO8) (1虚数軸固定) = 8 次元 (8 グルーオン!)
    # =========================================================================
    log("=" * 78)
    log("Part 1. 外側 O8 結晶の厳密ゲージ対称性次元カウント (G2 = 14, SU(3) = 8)")
    log("=" * 78)

    def count_kernel(label, pairs):
        cols = []
        for (a, b) in pairs:
            G = np.zeros((8, 8))
            G[a, b] = -1.0
            G[b, a] = +1.0
            dO = apply_outer_op(outer, G)
            cols.append(dO.ravel())
        M = np.column_stack(cols)  # (512, len(pairs))
        sv = np.linalg.svd(M, compute_uv=False)
        n_tot = len(pairs)
        n_nonzero = int(np.sum(sv > 1e-8))
        n_ker = n_tot - n_nonzero
        log(f"  {label:46s}: 総生成子={n_tot:2d}, 厳密不変(ゲージ)={n_ker:2d}, 破れ(無質量|λ|=1)={n_nonzero:2d}")
        return n_tot, n_ker, n_nonzero

    pairs_so8 = [(a, b) for a in range(8) for b in range(a + 1, 8)]
    pairs_so7 = [(a, b) for a in range(1, 8) for b in range(a + 1, 8)]
    pairs_so6 = [(a, b) for a in range(1, 7) for b in range(a + 1, 7)]  # A=0, 7 固定
    pairs_so3 = [(1, 2), (2, 3), (3, 1)]  # 四元数3空間

    count_kernel("SO(8) 全回転 (0..7)", pairs_so8)
    count_kernel("SO(7) 虚数7空間回転 (1..7) -> G2 検出", pairs_so7)
    count_kernel("SO(6) 1虚数軸(A=7)固定回転 (1..6) -> SU(3) 検出", pairs_so6)
    count_kernel("SO(3) 四元数3空間回転 (1..3)", pairs_so3)

    # 破れモード（SO(7) の非ゼロ 7 モード）がすべて |λ|=1.000 無質量であることを1本検算
    cols7 = []
    for (a, b) in pairs_so7:
        G = np.zeros((8, 8)); G[a, b] = -1; G[b, a] = 1
        cols7.append(apply_outer_op(outer, G).ravel())
    M7 = np.column_stack(cols7)
    U7, S7, _ = np.linalg.svd(M7, full_matrices=False)
    # 7 本の非ゼロ基底すべてで τ=1,5 の |λ| を確認
    H1, _ = step_tri(H)
    lams = []
    for k in range(7):
        D = tensor_prod_3d(U7[:, k].reshape(8, 8, 8), inner)
        D = D - np.vdot(H, D) * H
        D /= np.linalg.norm(D)
        Up = (H + 1e-5 * D) / np.linalg.norm(H + 1e-5 * D)
        Up1, _ = step_tri(Up)
        lams.append(proj_dist(Up1, H1) / 1e-5)
    log(f"  SO(7)/G2 の 7 本の破れモードの τ=1 応答 |λ|: {[round(x, 6) for x in lams]}")

    # =========================================================================
    # Part 2: 局在単一サイト凝縮体 (P4) 2 個の相互作用テスト (30 ステップ)
    # =========================================================================
    log("\n" + "=" * 78)
    log("Part 2. 局在単一サイト凝縮体 (P4) 2 個の相互作用テスト (τ = 0 -> 30)")
    log("=" * 78)

    def make_in_ch(b_pair):
        D = np.zeros_like(inner)
        for ib, b in enumerate(b_pair):
            M = inner[np.ix_([4, 5], [ib], [4, 5])][:, 0, :]
            for ia, a in enumerate([4, 5]):
                for ic, c in enumerate([4, 5]):
                    D[a, b, c] = M[ia, ic]
        return D

    in_b23 = make_in_ch([2, 3])
    in_b67 = make_in_ch([6, 7])

    def loc_probe(site_A, in_part):
        out_s = np.zeros((8, 8, 8))
        out_s[site_A, site_A, site_A] = 1.0
        P = tensor_prod_3d(out_s, in_part)
        return P / np.linalg.norm(P)

    P_A0_23 = loc_probe(0, in_b23)
    P_A1_23 = loc_probe(1, in_b23)
    P_A0_67 = loc_probe(0, in_b67)
    P_A1_67 = loc_probe(1, in_b67)

    mask_H = np.abs(H) > 1e-8
    mask_A0_23 = np.abs(P_A0_23) > 1e-8
    mask_A1_23 = np.abs(P_A1_23) > 1e-8
    mask_A0_67 = np.abs(P_A0_67) > 1e-8
    mask_A1_67 = np.abs(P_A1_67) > 1e-8

    def test_2p(label, P1, P2, m1, m2, sign=+1.0):
        eps = 0.08
        Um = H + eps * P1 + sign * eps * P2
        Um /= np.linalg.norm(Um)
        e1_0 = float(np.sum(Um[m1] ** 2))
        e2_0 = float(np.sum(Um[m2] ** 2))
        for _ in range(30):
            Um, nrm = step_tri(Um)
        eH = float(np.sum(Um[mask_H] ** 2))
        e1 = float(np.sum(Um[m1] ** 2))
        e2 = float(np.sum(Um[m2] ** 2))
        # m1 と m2 が重なる場合に対応
        m_union = mask_H | m1 | m2
        e_cloth = float(np.sum(Um[~m_union] ** 2))
        nz = int(np.sum(np.abs(Um) > 1e-5))
        log(f"  {label:46s}: nz={nz:4d}, ||P||={nrm:.5f}, E1={e1:.6f}, E2={e2:.6f}, 衣E_oth={e_cloth:.6f}")

    test_2p("(a1) 同一サイト A=0 同チャネル 同相 (+)", P_A0_23, P_A0_23, mask_A0_23, mask_A0_23, +1.0)
    test_2p("(a2) 同一サイト A=0 同チャネル 逆相 (-)", P_A0_23, P_A0_23, mask_A0_23, mask_A0_23, -1.0)
    test_2p("(b1) 別サイト A=0,1 同チャネル 同相 (+)", P_A0_23, P_A1_23, mask_A0_23, mask_A1_23, +1.0)
    test_2p("(b2) 別サイト A=0,1 同チャネル 逆相 (-)", P_A0_23, P_A1_23, mask_A0_23, mask_A1_23, -1.0)
    test_2p("(c1) 同一サイト A=0 別チャネル (23 & 67)", P_A0_23, P_A0_67, mask_A0_23, mask_A0_67, +1.0)
    test_2p("(c2) 別サイト A=0,1 別チャネル (23 & 67)", P_A0_23, P_A1_67, mask_A0_23, mask_A1_67, +1.0)

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_g2_and_twoparticle_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
