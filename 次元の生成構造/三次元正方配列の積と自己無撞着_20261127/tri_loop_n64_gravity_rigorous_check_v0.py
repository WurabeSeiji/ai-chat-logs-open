#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
本丸（重力時間遅延・臨界崩壊・計量・曲率）の厳密反証チェック・精密検証スクリプト
  Check 1: 複素位相角 θ(A, τ) = atan2(c_J, c_0) の直接追跡
           -> A=0 の時計は本当に「周波数 ω(0) < ω(far)」で遅れ続けるのか、
              それとも「定常位相差（位相ズレのみ）」なのかを τ=1..120 で完全判別！
  Check 2: 弱重力極限 ε ∈ [0.02, 0.25] の精密スキャン (τ=30..100 長時間平均)
           -> 遅延率 δphi = 1 - phi(0)/phi(far) と計量窪み δg = 1 - E_bg(0)/E_bg(far) の
              比 δphi/δg が弱重力極限で 1/2 (= sqrt(g_00)) に収束するか、別指数かを精密測定。
  Check 3: 臨界崩壊 ε ∈ [0.25, 0.30] の精密スキャン (0.005 刻み) と崩壊後状態の完全解剖
           -> 単なるチャネル乗っ取りか、A=0 共存＋遠方 A=1..7 消滅か、nz・SVD・エネルギー内訳を確認。
  Check 4: 虚数サイト A=1 に置いた質量 P4_A1 での対照実験（実数軸 A=0 の特殊性でないことの確認）。
  Check 5: 2 質量 (A=0, A=1) 同時配置での時間遅延・重力井戸の重ね合わせ確認。
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


def extract_bg_site_complex_phases(U):
    """
    背景セクター b=0 において、各外側サイト A=0..7 の対角成分 (A, B=0, C=A) の
    内側 2x2 ブロック (a,c in {4,5}) から、
      - 複素位相 θ(A) = atan2(c_J, c_0)
      - 複素振幅 r(A) = hypot(c_0, c_J)
      - 異方性漏れ aniso(A) = hypot(c_x, c_z) / r(A)
    を返す。
    """
    T6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)
    theta = np.zeros(8)
    rad = np.zeros(8)
    aniso = np.zeros(8)
    for A in range(8):
        # O8[A, 0, A] は A=0 で +1、A=1..7 でも +1 (e_A * 1 = e_A)
        B = T6[A, 0, A, 4:6, 0, 4:6]
        c0 = 0.5 * (B[0, 0] + B[1, 1])
        cz = 0.5 * (B[0, 0] - B[1, 1])
        cx = 0.5 * (B[0, 1] + B[1, 0])
        cJ = 0.5 * (B[1, 0] - B[0, 1])
        r = float(np.hypot(c0, cJ))
        theta[A] = float(np.arctan2(cJ, c0))
        rad[A] = r
        aniso[A] = float(np.hypot(cx, cz) / r) if r > 1e-15 else 0.0
    return theta, rad, aniso


def site_metrics(U_prev, U_curr, b_slice):
    T_p = U_prev.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, b_slice, :]
    T_c = U_curr.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, b_slice, :]
    E = np.zeros(8)
    phi = np.zeros(8)
    # また、サイト A の中で b=0 と b=1 のエネルギー比（内部構造の揺れ vs 純位相回転の分離）
    b0_frac = np.zeros(8)
    for A in range(8):
        vp = T_p[A].ravel()
        vc = T_c[A].ravel()
        np_v = float(np.linalg.norm(vp))
        nc_v = float(np.linalg.norm(vc))
        E[A] = nc_v ** 2
        if np_v > 1e-14 and nc_v > 1e-14:
            cval = float(np.clip(np.dot(vp, vc) / (np_v * nc_v), -1.0, 1.0))
            phi[A] = float(np.arccos(abs(cval)))
            e_b0 = float(np.sum(T_c[A, :, :, :, 0, :] ** 2))
            b0_frac[A] = e_b0 / E[A]
    return E, phi, b0_frac


def wrap_pi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    log("=" * 78)
    log("本丸・厳密反証チェック（複素位相直接追跡・弱重力べき則・臨界崩壊解剖・対称性）")
    log("=" * 78)

    # 第 1 層 H = U^(200)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt_svd = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt_svd[0].reshape(8, 8, 8)

    def make_P4_at(A_site):
        P = np.zeros((64, 64, 64), dtype=np.float64)
        for a in [4, 5]:
            for b in [2, 3]:
                for c in [4, 5]:
                    P[A_site * 8 + a, 0 * 8 + b, A_site * 8 + c] = M8[a, b - 2, c]
        return P / np.linalg.norm(P)

    P4_A0 = make_P4_at(0)
    P4_A1 = make_P4_at(1)

    # ------------------------------------------------------------------
    # Check 1: 複素位相角 θ(A, τ) の直接追跡 (eps=0.16, τ=1..100)
    # ------------------------------------------------------------------
    log("\n" + "-" * 78)
    log("[Check 1] 複素位相角 θ(A, τ) の直接追跡 (P4 at A=0, eps=0.16, τ=0..100)")
    log("  問い: θ(0, τ) - θ(1, τ) は時間とともに単調に開く（本当の周波数差）か、定常位相ズレか？")
    log("-" * 78)
    Ut = (H + 0.16 * P4_A0)
    Ut /= np.linalg.norm(Ut)
    th_prev, _, _ = extract_bg_site_complex_phases(Ut)
    dtheta_diff_list = []
    phase_gap_list = []
    for st in range(1, 101):
        Un, _ = step_tri(Ut)
        th_curr, rad_curr, aniso_curr = extract_bg_site_complex_phases(Un)
        # 3倍角写像からのズレ: delta_3 = wrap(th_curr - 3*th_prev)
        d3 = wrap_pi(th_curr - 3.0 * th_prev)
        gap_01 = wrap_pi(th_curr[0] - th_curr[1])
        gap_17 = wrap_pi(th_curr[1] - th_curr[7])
        E_bg, phi_bg, b0f = site_metrics(Ut, Un, slice(0, 2))
        if st >= 20:
            dtheta_diff_list.append(wrap_pi(th_curr[0] - th_prev[0]) - wrap_pi(th_curr[1] - th_prev[1]))
            phase_gap_list.append(gap_01)
        if st in [1, 2, 5, 10, 20, 30, 50, 75, 100]:
            log(
                f"  tau={st:3d}: th(0)-th(1)={gap_01:+.6e} rad, th(1)-th(7)={gap_17:+.2e} | "
                f"d3(0)={d3[0]:+.4e}, d3(1)={d3[1]:+.4e} | "
                f"aniso(0)={aniso_curr[0]:.2e} | "
                f"b0_frac(0)={b0f[0]:.6f}, b0_frac(1)={b0f[1]:.6f} | "
                f"phi(0)/phi(1)={phi_bg[0]/phi_bg[1]:.6f}"
            )
        th_prev = th_curr
        Ut = Un
    log(
        f"  [Summary τ=20..100] mean|th(0)-th(1)|={np.mean(np.abs(phase_gap_list)):.6e}, "
        f"max|th(0)-th(1)|={np.max(np.abs(phase_gap_list)):.6e}, "
        f"std(th(0)-th(1))={np.std(phase_gap_list):.6e}"
    )

    # ------------------------------------------------------------------
    # Check 2: 弱重力極限 ε ∈ [0.02..0.24] の精密スキャン (τ=30..80 平均)
    # ------------------------------------------------------------------
    log("\n" + "-" * 78)
    log("[Check 2] 弱重力べき則の精密スキャン (τ=30..80 平均): δphi/δg, δb0_frac, δr")
    log("-" * 78)
    log(
        f"{'eps':>6s} | {'dg=1-E0/Ef':>11s} {'dphi=1-p0/pf':>12s} {'dphi/dg':>9s} | "
        f"{'mean|th0-th1|':>13s} {'|b0f(0)-b0f(1)|':>15s} {'E_P_tot':>9s}"
    )
    for eps in [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.24]:
        Ut = (H + eps * P4_A0)
        Ut /= np.linalg.norm(Ut)
        dg_l, dp_l, gap_l, b0d_l, ep_l = [], [], [], [], []
        for st in range(1, 81):
            Un, _ = step_tri(Ut)
            if st >= 30:
                E_bg, phi_bg, b0f = site_metrics(Ut, Un, slice(0, 2))
                E_pr, _, _ = site_metrics(Ut, Un, slice(2, 4))
                th_c, _, _ = extract_bg_site_complex_phases(Un)
                dg_l.append(1.0 - E_bg[0] / np.mean(E_bg[1:]))
                dp_l.append(1.0 - phi_bg[0] / np.mean(phi_bg[1:]))
                gap_l.append(abs(wrap_pi(th_c[0] - th_c[1])))
                b0d_l.append(abs(b0f[0] - np.mean(b0f[1:])))
                ep_l.append(float(np.sum(E_pr)))
            Ut = Un
        dg = float(np.mean(dg_l))
        dp = float(np.mean(dp_l))
        ratio = dp / dg if abs(dg) > 1e-12 else 0.0
        log(
            f"{eps:6.3f} | {dg:11.7f} {dp:12.7f} {ratio:9.4f} | "
            f"{np.mean(gap_l):13.6e} {np.mean(b0d_l):15.6e} {np.mean(ep_l):9.6f}"
        )

    # ------------------------------------------------------------------
    # Check 3: 臨界崩壊 ε ∈ [0.25..0.30] の精密スキャン & 崩壊後状態の解剖
    # ------------------------------------------------------------------
    log("\n" + "-" * 78)
    log("[Check 3] 臨界崩壊 ε ∈ [0.250..0.300] 精密スキャン (τ=60 時点)")
    log("-" * 78)
    log(f"{'eps':>6s} | {'nz':>5s} {'||P||':>8s} | {'E_bg(0)':>9s} {'E_bg(far)':>9s} {'E_P(0)':>9s} {'E_P(far)':>9s} | {'s1':>7s} {'s2':>7s}")
    for eps in [0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.300]:
        Ut = (H + eps * P4_A0)
        Ut /= np.linalg.norm(Ut)
        nP = 0.0
        for _ in range(60):
            Ut, nP = step_tri(Ut)
        nz = int(np.sum(np.abs(Ut) > 1e-12))
        T6 = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)
        e_bg0 = float(np.sum(T6[0, :, :, :, 0:2, :] ** 2))
        e_bgf = float(np.sum(T6[1:, :, :, :, 0:2, :] ** 2))
        e_p0 = float(np.sum(T6[0, :, :, :, 2:4, :] ** 2))
        e_pf = float(np.sum(T6[1:, :, :, :, 2:4, :] ** 2))
        _, sv, _ = np.linalg.svd(T6.reshape(512, 512), full_matrices=False)
        log(
            f"{eps:6.3f} | {nz:5d} {nP:8.5f} | {e_bg0:9.6f} {e_bgf:9.6f} {e_p0:9.6f} {e_pf:9.6f} | "
            f"{sv[0]:7.4f} {sv[1]:7.4f}"
        )

    # ------------------------------------------------------------------
    # Check 4: 虚数サイト A=1 に置いた P4_A1 での対照実験 (eps=0.16, 0.30)
    # ------------------------------------------------------------------
    log("\n" + "-" * 78)
    log("[Check 4] 対照実験：虚数サイト A=1 に質量 P4_A1 を置いた場合 (eps=0.16, 0.30)")
    log("-" * 78)
    for eps in [0.16, 0.30]:
        Ut = (H + eps * P4_A1)
        Ut /= np.linalg.norm(Ut)
        dg_l, dp_l = [], []
        for st in range(1, 61):
            Un, nP = step_tri(Ut)
            if st >= 30:
                E_bg, phi_bg, _ = site_metrics(Ut, Un, slice(0, 2))
                other_idx = [a for a in range(8) if a != 1]
                ef = float(np.mean(E_bg[other_idx]))
                pf = float(np.mean(phi_bg[other_idx]))
                if ef > 1e-12 and pf > 1e-12:
                    dg_l.append(1.0 - E_bg[1] / ef)
                    dp_l.append(1.0 - phi_bg[1] / pf)
            Ut = Un
        nz = int(np.sum(np.abs(Ut) > 1e-12))
        T6 = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)
        e_bg1 = float(np.sum(T6[1, :, :, :, 0:2, :] ** 2))
        e_bg_oth = float(np.sum(T6[:, :, :, :, 0:2, :] ** 2)) - e_bg1
        e_p1 = float(np.sum(T6[1, :, :, :, 2:4, :] ** 2))
        log(
            f"  P4_A1 (eps={eps:.2f}, tau=60): nz={nz}, ||P||={nP:.5f} | "
            f"E_bg(1)={e_bg1:.6f}, E_bg(oth)={e_bg_oth:.6f}, E_P(1)={e_p1:.6f} | "
            f"dg={np.mean(dg_l) if dg_l else np.nan:.6f}, dphi={np.mean(dp_l) if dp_l else np.nan:.6f}"
        )

    # ------------------------------------------------------------------
    # Check 5: 2 質量 (A=0, A=1) 同時配置での重ね合わせ確認 (各 eps=0.12)
    # ------------------------------------------------------------------
    log("\n" + "-" * 78)
    log("[Check 5] 2 質量 (A=0, A=1, 各 eps=0.12) 同時配置での時間遅延・重力井戸")
    log("-" * 78)
    Ut = (H + 0.12 * P4_A0 + 0.12 * P4_A1)
    Ut /= np.linalg.norm(Ut)
    E_bg_acc = np.zeros(8)
    phi_bg_acc = np.zeros(8)
    cnt = 0
    for st in range(1, 61):
        Un, _ = step_tri(Ut)
        if st >= 30:
            E_bg, phi_bg, _ = site_metrics(Ut, Un, slice(0, 2))
            E_bg_acc += E_bg
            phi_bg_acc += phi_bg
            cnt += 1
        Ut = Un
    E_bg_m = E_bg_acc / cnt
    phi_bg_m = phi_bg_acc / cnt
    ef = float(np.mean(E_bg_m[2:]))
    pf = float(np.mean(phi_bg_m[2:]))
    log(f"  Site A=0 (mass): E_bg/E_far={E_bg_m[0]/ef:.6f} (dg={1-E_bg_m[0]/ef:.6f}), phi/phi_far={phi_bg_m[0]/pf:.6f} (dphi={1-phi_bg_m[0]/pf:.6f})")
    log(f"  Site A=1 (mass): E_bg/E_far={E_bg_m[1]/ef:.6f} (dg={1-E_bg_m[1]/ef:.6f}), phi/phi_far={phi_bg_m[1]/pf:.6f} (dphi={1-phi_bg_m[1]/pf:.6f})")
    log(f"  Sites A=2..7   : E_bg std={np.std(E_bg_m[2:]):.2e}, phi std={np.std(phi_bg_m[2:]):.2e}")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s.")

    with open("results/tri_loop_n64_gravity_rigorous_check_v0.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
