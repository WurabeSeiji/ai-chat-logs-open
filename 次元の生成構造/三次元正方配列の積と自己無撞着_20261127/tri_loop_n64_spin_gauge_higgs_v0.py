#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
空間回転 SO(3) の表現（スピン 0 = ヒッグス、スピン 1 = ゲージ/光子、スピン 2 = グラビトン）の
観点から、三角形ループ積の励起モードの応答（無質量 |λ|=1 vs 質量ギャップ |λ|<1 vs 不安定 |λ|>1）
を厳密に検証する。

Part 1: 第 1 層 H = O8 ⊗ M8 の外側 O8 内の 3 次元空間 {1,2,3} に対する
        SO(3) スピン 0 (1成分) / スピン 1 (3成分: 回転生成子 & ベクトル場) / スピン 2 (5成分: 計量変形)
        の微小励起 (ε = 1e-5) の 1〜10 ステップ応答。
Part 2: 第 4 層 (Case 1, τ = 1,000 完全静止固定点 U = a ⊗ b ⊗ c) の
        4×4×4 = 64 次元接空間における厳密ヤコビアン行列 J (64×64) の全固有値スペクトルと
        各固有モード（|λ|=1 無質量ゲージ/ゴールドストーン、0<|λ|<1 復元、|λ|=0 遮断）の対称性分類。

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


def apply_outer_linear_op(T8, G):
    """8×8×8 テンソル T8 の 3 軸すべてに 8×8 生成子 G を微分作用させる: (G⊗I⊗I + I⊗G⊗I + I⊗I⊗G) T8"""
    d1 = np.einsum("ia,ajk->ijk", G, T8)
    d2 = np.einsum("jb,ibk->ijk", G, T8)
    d3 = np.einsum("kc,ijc->ijk", G, T8)
    return d1 + d2 + d3


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    # 第 1 層 H = U^(200) を生成
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
    # Part 1: 第 1 層 H = O8 ⊗ M8 上での SO(3) スピン 0 / 1 / 2 励起応答
    # =========================================================================
    log("=" * 78)
    log("Part 1. 第 1 層 H = O8 ⊗ M8 上での SO(3) スピン 0 / 1 / 2 励起テスト")
    log("  (外側 O8 の四元数 3 次元空間 {1, 2, 3} に対する表現別応答, ε = 1e-5)")
    log("=" * 78)

    H_ref = [H.copy()]
    Ut = H.copy()
    for _ in range(10):
        Ut, _ = step_tri(Ut)
        H_ref.append(Ut.copy())

    def test_spin_mode(label, D, eps=1e-5):
        D = D - np.vdot(H, D) * H
        nd = float(np.linalg.norm(D))
        if nd < 1e-12:
            log(f"  {label:44s}: H と平行 (恒等変換 / 自明)")
            return
        D = D / nd
        Up = H + eps * D
        Up /= np.linalg.norm(Up)
        res = []
        for tau in range(1, 11):
            Up, _ = step_tri(Up)
            if tau in (1, 2, 3, 5, 10):
                dt = proj_dist(Up, H_ref[tau])
                res.append((tau, dt / eps))
        s = ", ".join(f"τ={k}:{v:7.3f}" for k, v in res)
        log(f"  {label:44s}: {s}")

    # --- スピン 0 (スカラー・1重項): 3空間軸 {1,2,3} の等方伸縮 ---
    G_s0 = np.zeros((8, 8))
    for a in (1, 2, 3):
        G_s0[a, a] = 1.0
    test_spin_mode("[Spin 0] 3空間 {1,2,3} 等方スカラー伸縮", tensor_prod_3d(apply_outer_linear_op(outer, G_s0), inner))

    # --- スピン 0 (スカラー): 残り4軸 {4,5,6,7} vs {1,2,3} 相対伸縮 ---
    G_s0b = np.zeros((8, 8))
    for a in (1, 2, 3):
        G_s0b[a, a] = 1.0
    for a in (4, 5, 6, 7):
        G_s0b[a, a] = -0.75
    test_spin_mode("[Spin 0] 3空間 vs 内部4空間 相対伸縮", tensor_prod_3d(apply_outer_linear_op(outer, G_s0b), inner))

    # --- スピン 1 (SO(3) 空間回転生成子 Lx, Ly, Lz: 3重項) ---
    # Lz: (1,2) 回転, Lx: (2,3) 回転, Ly: (3,1) 回転
    for name, (a, b) in [("Lx (2-3回転)", (2, 3)), ("Ly (3-1回転)", (3, 1)), ("Lz (1-2回転)", (1, 2))]:
        G_rot = np.zeros((8, 8))
        G_rot[a, b] = -1.0
        G_rot[b, a] = +1.0
        test_spin_mode(f"[Spin 1] SO(3) 空間回転 {name}", tensor_prod_3d(apply_outer_linear_op(outer, G_rot), inner))

    # --- スピン 1 (G2 自己同型: O8 全体の厳密ゲージ回転 = 1-2回転 + 5-6回転の組など) ---
    # O8 の標準自己同型: (1,2) を θ 回転するとき (4,7) や (5,6) も連動して積表を保つ生成子を数値探索で厳密チェック
    # まず (1,2) + (5,6) と (1,2) - (5,6) で積表不変性を確認
    for sgn, nm in [(+1, "(1,2)+(5,6)"), (-1, "(1,2)-(5,6)"), (+1, "(1,2)+(4,7)"), (-1, "(1,2)-(4,7)")]:
        G = np.zeros((8, 8))
        G[1, 2] = -1; G[2, 1] = 1
        if "5,6" in nm:
            G[5, 6] = -sgn; G[6, 5] = sgn
        else:
            G[4, 7] = -sgn; G[7, 4] = sgn
        dO = apply_outer_linear_op(outer, G)
        test_spin_mode(f"[Spin 1/G2?] 連動回転 {nm}", tensor_prod_3d(dO, inner))

    # --- スピン 2 (3空間 {1,2,3} のトレースレス対称計量変形: 5重項) ---
    # (1) diag(1, -1, 0)  (2) diag(1, 1, -2)  (3) sym(1,2)  (4) sym(2,3)  (5) sym(3,1)
    s2_ops = [
        ("h_xx - h_yy (diag 1,-1,0)", [(1, 1, 1.0), (2, 2, -1.0)]),
        ("h_xx + h_yy - 2h_zz", [(1, 1, 1.0), (2, 2, 1.0), (3, 3, -2.0)]),
        ("h_xy (sym 1-2)", [(1, 2, 1.0), (2, 1, 1.0)]),
        ("h_yz (sym 2-3)", [(2, 3, 1.0), (3, 2, 1.0)]),
        ("h_zx (sym 3-1)", [(3, 1, 1.0), (1, 3, 1.0)]),
    ]
    for nm, entries in s2_ops:
        G = np.zeros((8, 8))
        for r, c, val in entries:
            G[r, c] = val
        test_spin_mode(f"[Spin 2] 3空間計量変形 {nm}", tensor_prod_3d(apply_outer_linear_op(outer, G), inner))

    # --- i,k 共有チャネル b={2,3} に乗せたスピン 0 / 1 / 2 場（プローブ型ゲージ・ヒッグス・グラビトン）---
    log("\n--- (補足) i,k 共有 b={2,3} チャネル上のスピン 0 / 1 / 2 励起 ---")
    D_in_b = np.zeros_like(inner)
    for ib, b in enumerate([2, 3]):
        M = inner[np.ix_([4, 5], [ib], [4, 5])][:, 0, :]
        for ia, a in enumerate([4, 5]):
            for ic, c in enumerate([4, 5]):
                D_in_b[a, b, c] = M[ia, ic]
    # スカラー (O8 そのまま)、ベクトル (Lx 作用)、テンソル (h_xy 作用)
    test_spin_mode("[b=2,3] スカラー場 (O8 ⊗ M8)", tensor_prod_3d(outer, D_in_b))
    G_lx = np.zeros((8, 8)); G_lx[2, 3] = -1; G_lx[3, 2] = 1
    test_spin_mode("[b=2,3] ベクトル場 (Lx O8 ⊗ M8)", tensor_prod_3d(apply_outer_linear_op(outer, G_lx), D_in_b))
    G_hxy = np.zeros((8, 8)); G_hxy[1, 2] = 1; G_hxy[2, 1] = 1
    test_spin_mode("[b=2,3] テンソル場 (h_xy O8 ⊗ M8)", tensor_prod_3d(apply_outer_linear_op(outer, G_hxy), D_in_b))

    # =========================================================================
    # Part 2: 第 4 層 (Case 1, τ = 1,000 完全静止固定点) の厳密 64×64 ヤコビアン全固有値
    # =========================================================================
    log("\n" + "=" * 78)
    log("Part 2. 第 4 層 (Case 1, τ = 1,000 完全静止固定点) の厳密ヤコビアン固有値スペクトル")
    log("=" * 78)

    out_r = outer.copy()
    out_r[0, 1, 0] += 1e-4
    U4 = tensor_prod_3d(out_r, inner)
    U4 /= np.linalg.norm(U4)
    for _ in range(1000):
        U4, nrm4 = step_tri(U4)

    mask4 = np.abs(U4) > 1e-5
    idx_i = sorted(np.where(np.any(mask4, axis=(1, 2)))[0].tolist())
    idx_j = sorted(np.where(np.any(mask4, axis=(0, 2)))[0].tolist())
    idx_k = sorted(np.where(np.any(mask4, axis=(0, 1)))[0].tolist())
    sub = U4[np.ix_(idx_i, idx_j, idx_k)]
    log(f"第 4 層 4×4×4 部分テンソル: ||P||={nrm4:.8f}, 軸 i={idx_i}, j={idx_j}, k={idx_k}")

    # 4×4×4 (64次元) での写像 f(x) = Tri(x)/||Tri(x)|| の厳密ヤコビアン J (64×64) を中心差分で計算
    def f_map(x_vec):
        X = x_vec.reshape(4, 4, 4)
        P = tri_prod_fast(X)
        return (P / np.linalg.norm(P)).ravel()

    x0 = sub.ravel()
    eps_j = 1e-7
    J = np.zeros((64, 64), dtype=np.float64)
    for col in range(64):
        dx = np.zeros(64)
        dx[col] = eps_j
        J[:, col] = (f_map(x0 + dx) - f_map(x0 - dx)) / (2.0 * eps_j)

    evals, evecs = np.linalg.eig(J)
    order = np.argsort(np.abs(evals))[::-1]
    evals = evals[order]
    evecs = evecs[:, order]

    # 固有値の分類集計
    n_gt1 = int(np.sum(np.abs(evals) > 1.001))
    n_eq1 = int(np.sum((np.abs(evals) >= 0.999) & (np.abs(evals) <= 1.001)))
    n_mid = int(np.sum((np.abs(evals) > 1e-3) & (np.abs(evals) < 0.999)))
    n_zero = int(np.sum(np.abs(evals) <= 1e-3))

    log(f"\n64×64 厳密ヤコビアン固有値の分類:")
    log(f"  |λ| > 1 (不安定モード):          {n_gt1:2d} 個")
    log(f"  |λ| = 1 (無質量/平坦モード):     {n_eq1:2d} 個")
    log(f"  0 < |λ| < 1 (復元/質量ギャップ): {n_mid:2d} 個")
    log(f"  |λ| = 0 (即時遮断/超重):         {n_zero:2d} 個")

    log("\n非ゼロ固有値 (|λ| > 1e-4) 全一覧と、各固有モードの転置対称性 (i<->k):")
    for idx in range(64):
        lam = evals[idx]
        if abs(lam) <= 1e-4:
            continue
        v = evecs[:, idx].reshape(4, 4, 4)
        # i<->k 対称度
        sym_p = 0.5 * (v + v.transpose(2, 1, 0))
        s_ratio = float(np.sum(np.abs(sym_p) ** 2) / np.sum(np.abs(v) ** 2))
        # x0 との内積（接空間直交性チェック）
        ov0 = float(abs(np.vdot(x0, v.ravel())))
        log(f"  [{idx+1:2d}] λ = {lam.real:+.6f}{lam.imag:+.6f}i  (|λ|={abs(lam):.6f}) | i<->k対称={s_ratio:.3f}, |<v,U0>|={ov0:.3e}")

    log(f"\n所要時間: {time.time() - t0:.2f} 秒")
    out_path = os.path.join("results", "tri_loop_n64_spin_gauge_higgs_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"保存先: {out_path}")


if __name__ == "__main__":
    run()
