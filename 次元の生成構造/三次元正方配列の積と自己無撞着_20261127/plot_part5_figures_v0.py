#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part V（多殻大格子 O_8^{⊗K} ⊗ M_8 の厳密閉包エンジン・臨界崩壊スケーリング・
弱結合/強結合多殻重力・光円錐因果伝播・独立時計 H^{⊗2} g_00 監査）用の
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

    # (a) 臨界崩壊振幅 ε_c(K) の殻数 K=1..3 厳密スケーリング & K=4 外挿
    ax = axes[0]
    Ks = np.array([1, 2, 3])
    ecs = np.array([0.273619, 0.013474, 0.001997])
    ax.semilogy(Ks, ecs, "o-", color="#d62728", lw=2.5, ms=8.5,
                label=r"実測臨界崩壊点 $\epsilon_c(K)$（厳密閉包 $\tau=40$）" if ja else r"Measured critical collapse $\epsilon_c(K)$ (exact, $\tau=40$)")
    ax.semilogy([3, 4], [0.001997, 2.96e-4], "s--", color="#ff7f0e", lw=2.0, ms=7.2,
                label=r"$K=4$ 外挿予測 ($r_{23}=0.1482 \Rightarrow 2.96\times 10^{-4}$)" if ja else r"$K=4$ extrapolation ($r_{23}=0.1482 \Rightarrow 2.96\times 10^{-4}$)")
    for k, e, nstr in [(1, 0.273619, "N=64"), (2, 0.013474, "N=512"), (3, 0.001997, "N=4,096")]:
        ax.text(k, e * 1.38, f"{e:.4f}\n({nstr})", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color="#b22222")
    ax.text(4, 2.96e-4 * 1.45, "2.96e-4\n(N=32,768)", ha="center", va="bottom", fontsize=8.2, color="#d95f02", fontweight="bold")
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["$K=1$\n(8 sites)", "$K=2$\n(64 sites)", "$K=3$\n(512 sites)", "$K=4$\n(4,096 sites)"], fontsize=8.8)
    ax.set_ylim(1e-4, 1.5)
    ax.set_xlim(0.6, 4.4)
    ax.set_ylabel(r"事象の地平面臨界振幅 $\epsilon_c(K)$（対数軸）" if ja else r"Critical horizon amplitude $\epsilon_c(K)$ (log scale)", fontsize=10.2)
    ax.set_title(
        r"(a) 多殻大格子 $O_8^{\otimes K}\otimes M_8$ の臨界崩壊スケーリング"
        if ja else
        r"(a) Critical Horizon Scaling on Multi-Shell $O_8^{\otimes K}\otimes M_8$",
        fontsize=10.5, fontweight="bold"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.0, framealpha=0.94)
    ax.text(
        0.52, 0.72,
        ("厳密閉包次元: 外側 5D × 内側 4D\n"
         "$K=1$ 密配列対照誤差: $8.89\\times 10^{-14}$\n"
         "$K=4$ (280TB相当) を 9.5s で厳密計算")
        if ja else
        ("Exact closure: Outer 5D × Inner 4D\n"
         "$K=1$ dense control err: $8.89\\times 10^{-14}$\n"
         "$K=4$ (280TB eq.) exact in 9.5s"),
        transform=ax.transAxes, fontsize=7.9,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) 強結合地平面相 (ε=0.16, τ=40) の殻深さ d=1..K 線形勾配 Δφ(d) vs 弱結合 K=4 の 1/5 減衰
    ax = axes[1]
    d4 = np.array([1, 2, 3, 4])
    phi4_str = np.array([11.4299, 7.8055, 4.1808, 0.0000])
    d3 = np.array([1, 2, 3])
    phi3_str = np.array([7.8490, 4.0299, 0.0000])
    d2 = np.array([1, 2])
    phi2_str = np.array([4.0935, 0.0000])

    ax.plot(d4, phi4_str, "D-", color="#d62728", lw=2.4, ms=7.8,
            label=r"強結合 $K=4$ ($N=32{,}768$, 地平面外側 $d\geq 1$)" if ja else r"Strong $K=4$ ($N=32{,}768$, outside horizon $d\geq 1$)")
    ax.plot(d3, phi3_str, "s-", color="#1f77b4", lw=2.2, ms=7.2,
            label=r"強結合 $K=3$ ($N=4{,}096$)" if ja else r"Strong $K=3$ ($N=4{,}096$)")
    ax.plot(d2, phi2_str, "o-", color="#2ca02c", lw=2.1, ms=7.0,
            label=r"強結合 $K=2$ ($N=512$)" if ja else r"Strong $K=2$ ($N=512$)")

    ax.set_xticks([1, 2, 3, 4])
    ax.set_xlabel(r"中心からの殻深さ（ハミング距離）$d$" if ja else r"Shell depth (Hamming distance) $d$", fontsize=10.2)
    ax.set_ylabel(r"対数計量ポテンシャル $\Delta\phi(d) = \frac{1}{2}\ln(E_{\mathrm{bg}}(d)/E_{\mathrm{bg}}(K))$" if ja else r"Log metric potential $\Delta\phi(d) = \frac{1}{2}\ln(E_{\mathrm{bg}}(d)/E_{\mathrm{bg}}(K))$", fontsize=9.8)
    ax.set_ylim(-0.8, 13.0)
    ax.set_title(
        r"(b) 地平面外側の多殻勾配：等間隔 $\Delta\phi(d) \approx 3.8(K-d)$"
        if ja else
        r"(b) Multi-Shell Horizon Exterior: Linear $\Delta\phi(d) \approx 3.8(K-d)$",
        fontsize=10.4, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.1, framealpha=0.94)
    ax.text(
        0.04, 0.08,
        ("弱結合共存相 ($K=4, \\epsilon=0.6\\epsilon_c$):\n"
         "$d=0\\to 1\\to 2$ で減衰比 $\\approx 0.20\\;(1/5)$\n"
         "→ 対数変換 $r=5^d$ で $\\Delta\\phi(r) \\propto 1/r$ 型")
        if ja else
        ("Weak-field phase ($K=4, \\epsilon=0.6\\epsilon_c$):\n"
         "Decay ratio $\\approx 0.20\\;(1/5)$ for $d=0\\to 1\\to 2$\n"
         "→ Under log map $r=5^d$, $\\Delta\\phi(r) \\propto 1/r$"),
        transform=ax.transAxes, fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#d62728", alpha=0.94)
    )

    fig.suptitle(
        r"Part V 検証 (1)：多殻大格子 $O_8^{\otimes K}\otimes M_8$ の臨界崩壊スケーリングと多殻計量ポテンシャル"
        if ja else
        r"Part V Verification (1): Critical Horizon Scaling & Multi-Shell Metric Potential on $O_8^{\otimes K}\otimes M_8$",
        fontsize=11.8, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p5_multishell_scaling_and_potential_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), dpi=200)

    # (a) 弱結合 K=3 (N=4,096) の因果伝播フロント (τ=0..8)
    ax = axes[0]
    taus = np.arange(0, 9)
    d0 = np.abs(np.array([0.0, -0.001744, 0.000032, -0.000056, 0.000002, -0.000014, -0.000005, -0.000011, -0.000009]))
    d1 = np.abs(np.array([0.0, 0.0, 0.0, -0.000007, 0.000002, -0.000002, 0.0, -0.000002, 0.0]))
    d2 = np.abs(np.array([0.0, 0.0, 0.000016, -0.000003, 0.000002, -0.000001, 0.000001, 0.0, 0.000001]))

    ax.plot(taus, d0 * 1e4, "o-", color="#d62728", lw=2.3, ms=7.0,
            label=r"中心 $d=0$（$\tau=1$ で即時応答）" if ja else r"Center $d=0$ (responds at $\tau=1$)")
    ax.plot(taus, d2 * 1e4, "^-", color="#1f77b4", lw=2.1, ms=7.0,
            label=r"外殻 $d=2$（$\tau=2$ で到達）" if ja else r"Shell $d=2$ (arrives at $\tau=2$)")
    ax.plot(taus, d1 * 1e4, "s--", color="#2ca02c", lw=2.1, ms=6.8,
            label=r"中間殻 $d=1$（$\tau=3$ で到達）" if ja else r"Shell $d=1$ (arrives at $\tau=3$)")
    ax.set_xlabel(r"反復ステップ $\tau$" if ja else r"Iteration step $\tau$", fontsize=10.2)
    ax.set_ylabel(r"計量応答振幅 $|\Delta\phi(d)| \times 10^4$" if ja else r"Metric response amplitude $|\Delta\phi(d)| \times 10^4$", fontsize=10.2)
    ax.set_yscale("symlog", linthresh=0.02)
    ax.set_title(
        r"(a) 弱結合 $K=3$ ($N=4{,}096$) の有限速度因果伝播（光円錐フロント）"
        if ja else
        r"(a) Finite-Speed Causal Propagation on Weak-Field $K=3$ ($N=4{,}096$)",
        fontsize=10.2, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.2, framealpha=0.94)

    # (b) g_00(d) 厳密監査：スカラー弱結合ロック vs スピン2 (σ_x) 弱結合12桁解禁 vs 強結合崩壊解禁
    ax = axes[1]
    cases = [
        "共通 $K=4$\nスカラー弱\n($0.6\\epsilon_c$)" if ja else "Shared $K=4$\nScalar weak\n($0.6\\epsilon_c$)",
        "独立 $H^{\\otimes 2}$\nスカラー弱\n($\\epsilon=0.020$)" if ja else "Indep. $H^{\\otimes 2}$\nScalar weak\n($\\epsilon=0.020$)",
        "多殻 $K=2$\nスピン2 $\\sigma_x$弱\n($\\epsilon=0.005$)" if ja else "Multi $K=2$\nSpin-2 $\\sigma_x$\n($\\epsilon=0.005$)",
        "多殻 $K=2$\nスピン2 $\\sigma_x$弱\n($\\epsilon=0.015$)" if ja else "Multi $K=2$\nSpin-2 $\\sigma_x$\n($\\epsilon=0.015$)",
        "独立 $H^{\\otimes 2}$\n強結合崩壊\n($\\epsilon=0.050$)" if ja else "Indep. $H^{\\otimes 2}$\nCollapsed\n($\\epsilon=0.050$)",
        "共通 $K=4$\n強結合崩壊\n($\\epsilon=0.160$)" if ja else "Shared $K=4$\nCollapsed\n($\\epsilon=0.160$)",
    ]
    dth_vals = [1.732e-14, 1.346e-13, 0.1411, 0.8784, 2.759, 3.049]
    cols = ["#1f77b4", "#1f77b4", "#2ca02c", "#2ca02c", "#d62728", "#d62728"]
    x = np.arange(len(cases))
    ax.bar(x, dth_vals, width=0.58, color=cols, edgecolor="#222")
    ax.set_yscale("log")
    ax.set_ylim(1e-16, 60.0)
    ax.axhline(1e-9, color="#555", ls=":", lw=1.5)
    ax.text(-0.42, 3.0e-9, "単一時計ロック上限\n($10^{-9}$ rad)" if ja else "Single-clock lock\nceiling ($10^{-9}$ rad)", fontsize=7.6, color="#333", va="bottom")
    for xi, v, c in zip(x, dth_vals, cols):
        lbl = f"{v:.1e}" if v < 1e-3 else f"{v:.2f} rad"
        tcol = "#1f77b4" if c == "#1f77b4" else ("#1b7a1b" if c == "#2ca02c" else "#b22222")
        ax.text(xi, v * 1.9, lbl, ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=tcol)
    ax.set_xticks(x)
    ax.set_xticklabels(cases, fontsize=7.6)
    ax.set_ylabel(r"最大殻間時計位相差 $\max_d |\Delta\theta(d)|$ [rad]（対数軸）" if ja else r"Max inter-shell clock phase shift $|\Delta\theta(d)|$ [rad] (log)", fontsize=9.8)
    ax.set_title(
        r"(b) $g_{00}(d)$ 厳密監査：スカラー弱結合ロック vs スピン 2 ($\sigma_x$) 弱結合 12 桁解禁"
        if ja else
        r"(b) $g_{00}(d)$ Audit: Scalar Weak Lock vs. Spin-2 ($\sigma_x$) Weak-Field 12-Decade Unlock",
        fontsize=9.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        r"Part V 検証 (2)：多殻光円錐因果伝播と $g_{00}(d)$ 厳密監査（独立時計ロック vs スピン 2 弱結合解禁）"
        if ja else
        r"Part V Verification (2): Multi-Shell Causal Propagation & $g_{00}(d)$ Audit (Lock vs. Spin-2 Unlock)",
        fontsize=11.6, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p5_causal_and_clock_audit_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    for lang in ("ja", "en"):
        plot_fig1(lang)
        plot_fig2(lang)
