#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part VI（2 粒子動的散乱・重力合算/ゲージ電荷分離・遠隔スピン 2 散乱・
空間分離 CHSH ベル不等式の動的生成と地平面切断転移）用の
図 2 種（日英計 4 枚）生成スクリプト。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def setup_font(lang):
    if lang == "ja":
        plt.rcParams["font.sans-serif"] = ["Hiragino Sans", "Hiragino Maru Gothic Pro", "YuGothic", "DejaVu Sans"]
    else:
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["axes.unicode_minus"] = False


def plot_fig1(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), dpi=200)

    # (a) 空間分離 2 粒子間のスピン 2 (σ_x) 動的散乱：標的 A=0 位相シフト |δθ(τ=2)| ∝ ε_s2^2
    ax = axes[0]
    es = np.array([0.005, 0.01, 0.02, 0.04, 0.08])
    dth2 = np.array([3.7350e-5, 1.4940e-4, 5.9750e-4, 2.3886e-3, 9.5315e-3])
    dth8 = np.array([9.3240e-3, 3.7278e-2, 1.4883e-1, 5.9073e-1, 2.2890e0])
    fit2 = 1.4940 * (es ** 2)

    ax.loglog(es, dth2, "o-", color="#1f77b4", lw=2.3, ms=7.8,
              label=r"$\tau=2$ 到着時 $|\delta\theta(2)|$（実測）" if ja else r"Arrival $\tau=2$: $|\delta\theta(2)|$ (measured)")
    ax.loglog(es, fit2, "--", color="#1f77b4", lw=1.5, alpha=0.75,
              label=r"厳密二乗則 $1.494\,\epsilon_{\mathrm{s2}}^2$（4桁一致）" if ja else r"Exact quadratic $1.494\,\epsilon_{\mathrm{s2}}^2$ (4-digit fit)")
    ax.loglog(es, dth8, "s-", color="#d62728", lw=2.3, ms=7.5,
              label=r"$\tau=8$ 累積散乱 $|\delta\theta(8)|$（最大 2.29 rad）" if ja else r"Accumulated $\tau=8$: $|\delta\theta(8)|$ (up to 2.29 rad)")
    ax.set_xlabel(r"遠隔サイト $B=1$ スピン 2 ($\sigma_x$) 源振幅 $\epsilon_{\mathrm{s2}}$" if ja else r"Remote site $B=1$ spin-2 ($\sigma_x$) amplitude $\epsilon_{\mathrm{s2}}$", fontsize=10.0)
    ax.set_ylabel(r"標的 $A=0$ 散乱位相シフト $|\delta\theta(\tau)|$ [rad]" if ja else r"Target $A=0$ scattering phase shift $|\delta\theta(\tau)|$ [rad]", fontsize=10.0)
    ax.set_title(
        r"(a) 空間分離 2 粒子間のスピン 2 ($\sigma_x$) 動的散乱スケーリング"
        if ja else
        r"(a) Remote Two-Particle Spin-2 ($\sigma_x$) Dynamical Scattering Scaling",
        fontsize=10.3, fontweight="bold"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.0, framealpha=0.94)
    ax.text(
        0.04, 0.74,
        ("因果律テスト ($\\tau=1$):\n"
         "全振幅で $\\delta\\theta(\\tau=1) = 0.0000$\n"
         "$\\tau=2$ で厳密 $\\propto \\epsilon_{\\mathrm{s2}}^2$ 到着")
        if ja else
        ("Causality test ($\\tau=1$):\n"
         "$\\delta\\theta(\\tau=1) = 0.0000$ at all $\\epsilon_{\\mathrm{s2}}$\n"
         "Arrives at $\\tau=2$ with exact $\\propto \\epsilon_{\\mathrm{s2}}^2$"),
        transform=ax.transAxes, fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) ゲージリンク媒介による同種電荷 (+,+) vs 異種電荷 (+,-) の相互作用分裂 & 2 粒子重力合算
    ax = axes[1]
    gs = np.array([0.002, 0.005, 0.010, 0.020])
    d_same = np.array([0.23111, 0.049178, 0.093327, 0.0048884])
    d_opp = np.array([0.077669, 0.020825, 0.031906, 0.010441])
    x = np.arange(len(gs))
    w = 0.35
    ax.bar(x - w / 2, d_same, width=w, color="#d62728", edgecolor="#222",
           label=r"同種電荷結合 $q=(+1,+1)$" if ja else r"Like charges $q=(+1,+1)$")
    ax.bar(x + w / 2, d_opp, width=w, color="#1f77b4", edgecolor="#222",
           label=r"異種電荷結合 $q=(+1,-1)$" if ja else r"Opposite charges $q=(+1,-1)$")
    ax.set_xticks(x)
    ax.set_xticklabels([f"g={g:.3f}" for g in gs], fontsize=9.0)
    ax.set_xlabel(r"共有中間サイト $C=2$ ゲージリンク結合 $g$" if ja else r"Shared intermediate site $C=2$ gauge link coupling $g$", fontsize=10.0)
    ax.set_ylabel(r"非線形相互作用強度 $|\Delta\|P\|_{\mathrm{int}}|$ ($\tau=20$)" if ja else r"Nonlinear interaction strength $|\Delta\|P\|_{\mathrm{int}}|$ ($\tau=20$)", fontsize=10.0)
    ax.set_title(
        r"(b) $U(1)$ ゲージリンク媒介による同種 vs 異種電荷分裂"
        if ja else
        r"(b) Like vs. Opposite Charge Splitting via $U(1)$ Gauge Links",
        fontsize=10.3, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.4, framealpha=0.94)
    ax.text(
        0.36, 0.50,
        ("ゲージ場なし対照 [SC-1, SC-3]:\n"
         "・局所 $U(1)$ 位相差 = 厳密 `0.000`\n"
         "・2 粒子重力崩壊: $\\epsilon_c^{(2)}/\\epsilon_c^{(1)} = 0.9842$\n"
         "  （符号非依存の普遍引力合算）")
        if ja else
        ("No-gauge control [SC-1, SC-3]:\n"
         "• Local $U(1)$ phase diff = exact `0.000`\n"
         "• 2-body horizon: $\\epsilon_c^{(2)}/\\epsilon_c^{(1)} = 0.9842$\n"
         "  (sign-independent universal gravity)"),
        transform=ax.transAxes, fontsize=7.9,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#d62728", alpha=0.94)
    )

    fig.suptitle(
        r"Part VI 検証 (1)：空間分離 2 粒子間のスピン 2 動的散乱と $U(1)$ ゲージ媒介電荷分離"
        if ja else
        r"Part VI Verification (1): Remote Spin-2 Dynamical Scattering & $U(1)$ Gauge-Mediated Charge Splitting",
        fontsize=11.6, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p6_two_particle_scattering_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), dpi=200)

    # (a) 空間分離 2 サイト (Alice=1, Bob=2) の CHSH ベル不等式動的生成・維持
    ax = axes[0]
    taus = np.array([0, 1, 2, 3, 5, 10, 20, 40])
    s_prod = np.array([2.146230, 2.697259, 2.787219, 2.825880, 2.807192, 2.804504, 2.823910, 2.811652])
    s_phi = np.array([2.793107, 2.827395, 2.825559, 2.827394, 2.827391, 2.827390, 2.827391, 2.827391])
    s_psi = np.array([2.827301, 2.827521, 2.828041, 2.828138, 2.828183, 2.828188, 2.828188, 2.828188])

    ax.plot(taus, s_prod, "o-", color="#1f77b4", lw=2.3, ms=7.2,
            label=r"積状態 $|\uparrow_1\uparrow_2\rangle$ 投入（動的生成）" if ja else r"Product $|\uparrow_1\uparrow_2\rangle$ (dynamical generation)")
    ax.plot(taus, s_phi, "s--", color="#2ca02c", lw=2.1, ms=6.8,
            label=r"ベル状態 $|\Phi^+\rangle$ 投入（動的維持）" if ja else r"Bell $|\Phi^+\rangle$ (dynamical preservation)")
    ax.plot(taus, s_psi, "^-", color="#d62728", lw=2.1, ms=6.8,
            label=r"ベル状態 $|\Psi^-\rangle$ 投入（$2.8282$ 維持）" if ja else r"Bell $|\Psi^-\rangle$ ($2.8282$ preserved)")
    ax.axhline(2.0, color="#555", ls=":", lw=1.6, label=r"局所実在論（古典）限界 $S=2$" if ja else r"Classical Bell bound $S=2$")
    ax.axhline(2 * np.sqrt(2), color="#9467bd", ls="-.", lw=1.4, label=r"チレルソン量子限界 $2\sqrt{2}\approx 2.8284$" if ja else r"Tsirelson bound $2\sqrt{2}\approx 2.8284$")
    ax.set_xlabel(r"反復ステップ $\tau$" if ja else r"Iteration step $\tau$", fontsize=10.2)
    ax.set_ylabel(r"空間分離 CHSH 相関値 $S_{\mathrm{CHSH}}(\tau)$" if ja else r"Spatially separated CHSH value $S_{\mathrm{CHSH}}(\tau)$", fontsize=10.2)
    ax.set_ylim(1.85, 2.95)
    ax.set_title(
        r"(a) 空間分離 2 サイト (Alice=1, Bob=2) の CHSH 動的生成・維持"
        if ja else
        r"(a) Dynamical Generation & Preservation of CHSH Violation (Alice=1, Bob=2)",
        fontsize=10.1, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=7.8, framealpha=0.94)

    # (b) 事象の地平面崩壊に伴う CHSH ベル不等式破れの不連続破壊（ブラックホール・エンタングルメント転移）
    ax = axes[1]
    eps_m = np.array([0.00, 0.05, 0.15, 0.25, 0.27, 0.28, 0.35])
    s_vals = np.array([2.8284, 2.8284, 2.8284, 2.8284, 2.8284, 0.0000, 0.0000])
    pur_vals = np.array([1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.2500, 0.2500])

    ax.plot(eps_m, s_vals, "o-", color="#d62728", lw=2.5, ms=8.0,
            label=r"遠隔 2 サイト CHSH 値 $S_{\mathrm{CHSH}}(\tau=40)$" if ja else r"Remote CHSH value $S_{\mathrm{CHSH}}(\tau=40)$")
    ax.plot(eps_m, pur_vals * 2.0, "s--", color="#1f77b4", lw=2.0, ms=6.8,
            label=r"純粋度 $2\times\mathrm{Tr}(\rho_{AB}^2)$（$1.0\to 0.25$）" if ja else r"Purity $2\times\mathrm{Tr}(\rho_{AB}^2)$ ($1.0\to 0.25$)")
    ax.axvline(0.27358, color="#333", ls=":", lw=1.8)
    ax.axhline(2.0, color="#555", ls=":", lw=1.4)
    ax.text(0.276, 1.55, r"地平面臨界点" + "\n" + r"$\epsilon_c \approx 0.2736$" if ja else r"Horizon threshold" + "\n" + r"$\epsilon_c \approx 0.2736$",
            fontsize=8.2, fontweight="bold", color="#333")
    ax.text(0.03, 2.48, "弱結合共存相:\n$S = 2.8284$ (VIOLATED)" if ja else "Weak coexistence:\n$S = 2.8284$ (VIOLATED)",
            fontsize=8.2, fontweight="bold", color="#b22222",
            bbox=dict(boxstyle="round,pad=0.25", fc="#fff8f8", ec="#d62728", alpha=0.92))
    ax.text(0.285, 0.16, "地平面崩壊相:\n$S = 0.0000$\n(完全混合 $\\mathbb{I}/4$)" if ja else "Collapsed phase:\n$S = 0.0000$\n(Max. mixed $\\mathbb{I}/4$)",
            fontsize=8.0, fontweight="bold", color="#1f77b4",
            bbox=dict(boxstyle="round,pad=0.25", fc="#f4faff", ec="#1f77b4", alpha=0.92))
    ax.set_xlabel(r"中心 $A=0$ スカラー質量振幅 $\epsilon_{\mathrm{mass}}$" if ja else r"Center $A=0$ scalar mass amplitude $\epsilon_{\mathrm{mass}}$", fontsize=10.2)
    ax.set_ylabel(r"CHSH 相関値 $S_{\mathrm{CHSH}}$ / スケール純粋度" if ja else r"CHSH value $S_{\mathrm{CHSH}}$ / Scaled purity", fontsize=10.2)
    ax.set_ylim(-0.15, 3.05)
    ax.set_xlim(-0.01, 0.37)
    ax.set_title(
        r"(b) 事象の地平面崩壊に伴う CHSH 破れの不連続破壊（$2\sqrt{2}\to 0$）"
        if ja else
        r"(b) Discontinuous Entanglement Collapse at the Horizon ($2\sqrt{2}\to 0$)",
        fontsize=10.1, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", fontsize=8.0, framealpha=0.94)

    fig.suptitle(
        r"Part VI 検証 (2)：空間分離 2 サイト CHSH ベル不等式の動的破れと事象の地平面エンタングルメント転移"
        if ja else
        r"Part VI Verification (2): Spatially Separated CHSH Violation & Horizon Entanglement Collapse",
        fontsize=11.6, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p6_chsh_and_horizon_transition_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    for lang in ("ja", "en"):
        plot_fig1(lang)
        plot_fig2(lang)
