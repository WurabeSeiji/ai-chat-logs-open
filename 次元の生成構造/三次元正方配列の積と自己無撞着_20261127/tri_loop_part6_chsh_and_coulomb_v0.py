#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part VI 精密検証 [SC-1b, SC-2b, SC-4b, SC-5b]：
  - [SC-1b] U(1) ゲージリンク媒介による同種電荷 (+,+) vs 異種電荷 (+,-) の相互作用エネルギー・計量分裂
  - [SC-2b] 空間分離 2 粒子間のスピン 2 (σ_x) 動的散乱（因果遅延 τ=1 で 0 -> τ>=2 で位相散乱）の振幅スケーリング
  - [SC-4b] 正しい 2 量子ビット指標 (Alice=(A, a in {4,5}), Bob=(B, c in {4,5})) での
            空間分離 2 サイト CHSH ベル不等式の動的破れ（真空基底・積状態・ベル状態・Tri(T) 40 ステップ発展）
  - [SC-5b] 振幅掃引による CHSH 破れの動的安定性と地平面崩壊時のエンタングルメント転移
乱数なし。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d
from tri_loop_part6_scattering_chsh_v0 import (
    step_tri,
    build_H_and_M8,
    build_site_scalar,
    build_site_spin2,
    build_gauge_packet,
    site_observables,
    chsh_metrics,
)


def extract_ab_rho_13(U, site_A, site_B, ch="all"):
    """
    M_8 の物理二重項構造に従い、第 1 指標 (Alice = site_A, a in {4,5}) と
    第 3 指標 (Bob = site_B, c in {4,5}) から 4x4 密度行列 ρ_AB = Tr_Env(|Ψ><Ψ|) を抽出。
    """
    U6 = U.reshape(8, 8, 8, 8, 8, 8)  # (A,a, M,m, B,c)
    if ch == "all":
        amp = U6[site_A, 4:6, :, :, site_B, 4:6].transpose(0, 3, 1, 2).reshape(4, 64)
    elif ch == "bg":
        amp = U6[site_A, 4:6, :, 0:2, site_B, 4:6].transpose(0, 3, 1, 2).reshape(4, 16)
    elif ch == "pr":
        amp = U6[site_A, 4:6, :, 2:4, site_B, 4:6].transpose(0, 3, 1, 2).reshape(4, 16)
    else:
        amp = U6[site_A, 4:6, :, :, site_B, 4:6].transpose(0, 3, 1, 2).reshape(4, 64)
    rho = amp @ amp.T
    tr = float(np.trace(rho))
    if tr < 1e-20:
        return np.eye(4) / 4.0, 0.0
    return rho / tr, tr


def build_link_bell_perturbation(site_A, site_B, mid_C, state_type="psi_minus", ch_b=(0, 1)):
    """
    O_8[A, C, B] != 0 となる空間リンク (A, C, B) の第 1・第 3 指標 {4,5} x {4,5} に
    2 量子ビット状態（product / phi_plus / psi_minus）を埋め込む。
    """
    P = np.zeros((64, 64, 64), dtype=np.float64)
    mat = np.zeros((2, 2), dtype=np.float64)
    if state_type == "product_00":
        mat[0, 0] = 1.0
    elif state_type == "phi_plus":
        mat[0, 0] = 1.0 / np.sqrt(2.0)
        mat[1, 1] = 1.0 / np.sqrt(2.0)
    elif state_type == "psi_minus":
        mat[0, 1] = 1.0 / np.sqrt(2.0)
        mat[1, 0] = -1.0 / np.sqrt(2.0)
    for m in ch_b:
        for ia, a in enumerate([4, 5]):
            for ib, c in enumerate([4, 5]):
                P[site_A * 8 + a, mid_C * 8 + m, site_B * 8 + c] = mat[ia, ib]
    return P / np.linalg.norm(P)


def evolve_final(U0, steps=40):
    U = U0 / np.linalg.norm(U0)
    nrm = 0.0
    for _ in range(steps):
        U, nrm = step_tri(U)
    return U, nrm


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H, M8, normP_H = build_H_and_M8(O8)

    log("=" * 88)
    log("Part VI 精密検証：ゲージ媒介電荷分離・遠隔スピン 2 散乱・空間分離 CHSH ベル不等式")
    log("=" * 88)

    # ==================================================================
    # [SC-1b] U(1) ゲージリンク媒介による同種電荷 (+,+) vs 異種電荷 (+,-) の分裂
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-1b] U(1) ゲージリンク媒介 2 粒子相互作用：同種 (+,+) vs 異種 (+,-) (A=0, B=1, τ=20)")
    log("-" * 88)
    P0 = build_site_scalar(0, M8)
    P1 = build_site_scalar(1, M8)
    G01 = build_gauge_packet(0, 1, M8)  # 反対称リンク (0->1: +, 1->0: -)
    # 荷電粒子対＋媒介ゲージ場：
    #   同種 (+,+): P0 + P1 + q_A*G01 + q_B*(-G01) ではなく、A->B の向きに依存する結合として
    #   粒子 A=0 の局所電荷結合 (+P0 + g G01) と 粒子 B=1 の局所電荷結合 (+P1 ± g G01)
    #   -> 同種 (q_A=q_B=+1): G01 が強め合う (+2g G01) vs 異種 (q_A=+1, q_B=-1): 相殺 (0) または反転
    #   さらに非自明な検証として、A=0 に +G01、B=1 に +G12 (共有しない別リンク) vs +G02 を比較
    G12 = build_gauge_packet(1, 2, M8)
    G02 = build_gauge_packet(0, 2, M8)
    # 共有中間サイト C=2 を介した 2 荷電粒子（A=0 と B=1 が中間サイト 2 のゲージ場に結合）
    # q_A = +1: P0 + g*G02, q_B = ±1: P1 ± g*G12
    _, n_vac20 = evolve_final(H, 20)
    log(f"  {'配置 (中間サイト 2 媒介)':<30} | {'g':>6} | {'||P||(20)':>10} | {'Δ||P||_int':>11} | {'E_bg(0)/E_bg(7)':>15} | {'E_bg(2)/E_bg(7)':>15}")
    for g in [0.002, 0.005, 0.010, 0.020]:
        eps = 0.08
        _, nA = evolve_final(H + eps * P0 + g * G02, 20)
        _, nBp = evolve_final(H + eps * P1 + g * G12, 20)
        _, nBm = evolve_final(H + eps * P1 - g * G12, 20)
        for label, U0, nb in [
            ("同種電荷結合 q=(+1,+1)", H + eps * P0 + eps * P1 + g * G02 + g * G12, nBp),
            ("異種電荷結合 q=(+1,-1)", H + eps * P0 + eps * P1 + g * G02 - g * G12, nBm),
        ]:
            Uf, nf = evolve_final(U0, 20)
            d_int = nf - nA - nb + n_vac20
            e0, _, _ = site_observables(Uf, 0)
            e2, _, _ = site_observables(Uf, 2)
            e7, _, _ = site_observables(Uf, 7)
            r07 = e0 / e7 if e7 > 1e-30 else float("inf")
            r27 = e2 / e7 if e7 > 1e-30 else float("inf")
            log(f"  {label:<30} | {g:6.3f} | {nf:10.6f} | {d_int:+11.4e} | {r07:15.6f} | {r27:15.6f}")

    # ==================================================================
    # [SC-2b] 空間分離 2 粒子間のスピン 2 (σ_x) 動的散乱スケーリング (A=0 標的, B=1 源)
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-2b] 空間分離 2 粒子間のスピン 2 (σ_x) 動的散乱：標的 A=0 位相シフト δθ(τ) の振幅スケーリング")
    log("-" * 88)
    S1_x = build_site_spin2(1)
    eps_t = 0.08
    U_f = (H + eps_t * P0) / np.linalg.norm(H + eps_t * P0)
    th_free = {}
    for st in range(1, 11):
        U_f, _ = step_tri(U_f)
        th_free[st] = site_observables(U_f, 0)[2]

    def wrap(x): return float((x + np.pi) % (2 * np.pi) - np.pi)
    log(f"  {'ε_s2':>6} | {'δθ(τ=1)':>12} | {'δθ(τ=2)':>12} | {'δθ(τ=4)':>12} | {'δθ(τ=8)':>12} | {'δθ(τ=2)/ε_s2^2':>15}")
    for es2 in [0.005, 0.01, 0.02, 0.04, 0.08]:
        U = (H + eps_t * P0 + es2 * S1_x) / np.linalg.norm(H + eps_t * P0 + es2 * S1_x)
        dths = {}
        for st in range(1, 11):
            U, _ = step_tri(U)
            if st in (1, 2, 4, 8):
                dths[st] = wrap(site_observables(U, 0)[2] - th_free[st])
        log(f"  {es2:6.3f} | {dths[1]:+12.4e} | {dths[2]:+12.4e} | {dths[4]:+12.4e} | {dths[8]:+12.4e} | {dths[2]/(es2**2):+15.4f}")

    # ==================================================================
    # [SC-4b] 空間分離 2 サイト (Alice=1, Bob=2) の CHSH ベル不等式（第 1・第 3 指標 {4,5}）
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-4b] 空間分離 2 サイト (Alice=1, Bob=2, 環境=全 64) の CHSH ベル不等式動的検証")
    log("-" * 88)
    # O_8[1, 3, 2] = -1 != 0 なので (A=1, C=3, B=2) は厳密非ゼロリンク！
    # まず真空 H そのものの (A=1, B=2) と (A=0, B=1) の 2 サイト量子相関を監査
    for (sA, sB) in [(1, 2), (0, 1), (1, 4)]:
        rho, tr = extract_ab_rho_13(H, sA, sB, ch="all")
        pur, se, smax, sstd, Tmat = chsh_metrics(rho)
        log(f"  [真空 H 基底 (A={sA}, B={sB})] Tr={tr:.4e}, purity={pur:.4f}, S_E={se:.4f} bit, S_max={smax:.6f}, S_std={sstd:.6f}")

    # 空間分離 2 サイト (A=1, B=2) に初期状態（積状態 vs ベル状態）を投入し Tri(T) 発展
    # 1) 真空背景への摂動 (ε=0.10, ch="all")
    for stype, label in [
        ("product_00", "積状態 |↑_1 ↑_2> (ch=bg)"),
        ("phi_plus",   "ベル状態 |Φ^+> (ch=bg)"),
        ("psi_minus",  "ベル状態 |Ψ^-> (ch=bg)"),
    ]:
        P_b = build_link_bell_perturbation(1, 2, 3, stype, ch_b=(0, 1))
        U = (H + 0.10 * P_b) / np.linalg.norm(H + 0.10 * P_b)
        log(f"\n  --- {label} (ε=0.10, A=1, B=2) ---")
        log(f"  {'τ':>3} | {'||P||':>9} | {'Tr(ρ_AB)':>10} | {'purity':>7} | {'S_E [bit]':>9} | {'S_max (Horodecki)':>17} | {'S_std (angle)':>13}")
        for st in range(0, 41):
            if st > 0:
                U, nrm = step_tri(U)
            else:
                nrm = normP_H
            if st in (0, 1, 2, 3, 5, 10, 20, 40):
                rho, tr = extract_ab_rho_13(U, 1, 2, ch="all")
                pur, se, smax, sstd, _ = chsh_metrics(rho)
                log(f"  {st:3d} | {nrm:9.6f} | {tr:10.3e} | {pur:7.4f} | {se:9.4f} | {smax:17.6f} | {sstd:13.6f}")

    # ==================================================================
    # [SC-5b] 振幅掃引：弱結合共存相 vs 臨界崩壊相（地平面形成）での空間分離 CHSH 破れ
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-5b] 振幅掃引（A=0 スカラー質量投入時の遠隔 2 サイト (1,2) CHSH 相関, τ=40）")
    log("-" * 88)
    log(f"  {'ε_mass':>7} | {'||P||(40)':>10} | {'相':<10} | {'purity':>7} | {'S_E [bit]':>9} | {'S_max':>9} | {'S_std':>9} | {'CHSH > 2?':>9}")
    for eps in [0.00, 0.05, 0.15, 0.25, 0.27, 0.28, 0.35]:
        Uf, nf = evolve_final(H + eps * P0, 40)
        rho, tr = extract_ab_rho_13(Uf, 1, 2, ch="all")
        pur, se, smax, sstd, _ = chsh_metrics(rho)
        phase = "Collapsed" if nf > 0.25 else "Coexist"
        flag = "VIOLATED" if smax > 2.0 + 1e-6 else "CLASSICAL"
        log(f"  {eps:7.2f} | {nf:10.6f} | {phase:<10} | {pur:7.4f} | {se:9.4f} | {smax:9.4f} | {sstd:9.4f} | {flag:>9}")

    log(f"\nDone in {time.time() - t0:.2f}s.")
    out_path = "results/tri_loop_part6_chsh_and_coulomb_v0.txt"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    run()
