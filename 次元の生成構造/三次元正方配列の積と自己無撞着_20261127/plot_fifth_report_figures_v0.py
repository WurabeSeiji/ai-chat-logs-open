#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 5 速報（球対称 4 クラス空間計量歪み・離散曲率・臨界特異点崩壊・スピン 2 時計引きずり）用の
図 2 種（日本語版・英語版、計 4 枚）を生成する専用スクリプト。
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
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=200)

    # (a) 球対称 4 クラス空間計量歪み (E / E_ref - 1 [%])
    ax = axes[0]
    classes = [
        ("中心自己\n(0,0,0)\n[1本]" if ja else "Center self\n(0,0,0)\n[1 edge]"),
        ("遠方接線\n(m,n,p)\n[42本]" if ja else "Far tangent\n(m,n,p)\n[42 edges]"),
        ("動径接続\n(0,m,m)\n[14本]" if ja else "Radial link\n(0,m,m)\n[14 edges]"),
        ("遠方自己\n(m,0,m)\n[7本, 基準]" if ja else "Far self\n(m,0,m)\n[7 edges, ref]"),
    ]
    e08 = [11.73232, 2.57321, 0.54810, 0.0]
    e16 = [32.48375, 7.18855, 1.22646, 0.0]
    x = np.arange(len(classes))
    w = 0.35
    ax.bar(x - w / 2, e08, width=w, color="#1f77b4", edgecolor="#222", label=r"$\epsilon = 0.08$ ($\tau=40$)")
    ax.bar(x + w / 2, e16, width=w, color="#d62728", edgecolor="#222", label=r"$\epsilon = 0.16$ ($\tau=40$)")
    for xi, v1, v2 in zip(x, e08, e16):
        if v1 > 0:
            ax.text(xi - w / 2, v1 + 0.7, f"+{v1:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#1f77b4")
            ax.text(xi + w / 2, v2 + 0.7, f"+{v2:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#d62728")
        else:
            ax.text(xi, 0.7, "0.0%\n(ref)", ha="center", va="bottom", fontsize=8.5, color="#444")
    ax.set_xticks(x)
    ax.set_xticklabels(classes, fontsize=9.2)
    ax.set_ylim(0, 39.5)
    ax.set_ylabel(r"空間エッジエネルギー比の増分 $(E/E_{\mathrm{ref}} - 1)$ [%]" if ja else r"Edge energy excess $(E/E_{\mathrm{ref}} - 1)$ [%]", fontsize=10.2)
    ax.set_title(
        "(a) 球対称 4 クラス空間計量歪み（7 方向 7 桁完全等方）"
        if ja else
        "(a) Spherically Symmetric 4-Class Metric Distortion (7-Digit Isotropic)",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=9.0, framealpha=0.92)
    ax.text(
        0.34, 0.56,
        ("動径エッジ (0,m,m) が最も細く抑えられ\n"
         "中心 (0,0,0) へ質量集中 (+32.5%)\n"
         "全エッジ複素位相ズレ < $10^{-9}$ rad")
        if ja else
        ("Radial links (0,m,m) are suppressed\n"
         "while center (0,0,0) concentrates (+32.5%)\n"
         "All edge phase diffs < $10^{-9}$ rad"),
        transform=ax.transAxes, fontsize=8.4,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#d62728", alpha=0.94)
    )

    # (b) 離散曲率 (重み付き視野角 w_mean_rho) の真空平坦 vs 質量直上 A=0 vs 遠方 A=1,7
    ax = axes[1]
    labels = [
        "真空 $H$\n(平坦: $A=0,1,7$)\n$n_{\\mathrm{tri}}=64$" if ja else "Vacuum $H$\n(Flat: $A=0,1,7$)\n$n_{\\mathrm{tri}}=64$",
        "質量あり: 遠方\n($A=1,7$, $\\epsilon=0.20$)\n$n_{\\mathrm{tri}}=72$" if ja else "With Mass: Far\n($A=1,7$, $\\epsilon=0.20$)\n$n_{\\mathrm{tri}}=72$",
        "質量あり: 直上\n($A=0$, $\\epsilon=0.20$)\n$n_{\\mathrm{tri}}=72$" if ja else "With Mass: Center\n($A=0$, $\\epsilon=0.20$)\n$n_{\\mathrm{tri}}=72$",
    ]
    w_rhos = [0.576699, 0.787796, 0.796889]
    m_rhos = [0.621227, 0.709985, 0.716065]
    xx = np.arange(len(labels))
    ax.bar(xx - 0.18, w_rhos, width=0.34, color=["#7f7f7f", "#ff7f0e", "#d62728"], edgecolor="#222",
           label=r"重み付き視野角 $\langle\rho\rangle_w$" if ja else r"Weighted view-angle $\langle\rho\rangle_w$")
    ax.bar(xx + 0.18, m_rhos, width=0.34, color=["#c7c7c7", "#ffbb78", "#ff9896"], edgecolor="#222",
           label=r"算術平均視野角 $\bar{\rho}$" if ja else r"Mean view-angle $\bar{\rho}$")
    for xi, wr, mr in zip(xx, w_rhos, m_rhos):
        ax.text(xi - 0.18, wr + 0.012, f"{wr:.4f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
        ax.text(xi + 0.18, mr + 0.012, f"{mr:.4f}", ha="center", va="bottom", fontsize=8.2, color="#444")
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, fontsize=9.2)
    ax.set_ylim(0.45, 0.90)
    ax.set_ylabel(r"三角形ループ視野角 $\rho = L_3 / L_1$" if ja else r"Triangle view-angle $\rho = L_3 / L_1$", fontsize=10.5)
    ax.set_title(
        r"(b) 離散曲率（視野角 $\rho$）の平坦性破れと球対称勾配"
        if ja else
        r"(b) Discrete Curvature ($\rho$) Gradient: Flat Vacuum vs Mass",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.6, framealpha=0.92)
    ax.annotate(
        r"曲率差 $+0.0091$" if ja else r"Curvature gap $+0.0091$",
        xy=(1.82, 0.797), xytext=(1.15, 0.855),
        arrowprops=dict(arrowstyle="->", color="#d62728"), fontsize=8.8, fontweight="bold", color="#d62728"
    )

    fig.suptitle(
        r"局在質量 P4 が引き起こす球対称 4 クラス空間計量歪みと離散曲率勾配（第 1 層 $H = O_8 \otimes M_8$）"
        if ja else
        r"Spherically Symmetric 4-Class Spatial Metric Distortion and Discrete Curvature Gradient on $H = O_8 \otimes M_8$",
        fontsize=12.3, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r5_metric_and_curvature_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=200)

    # (a) 臨界質量 ε_c ∈ (0.270, 0.275) での遠方空間切断・1サイト特異点崩壊
    ax = axes[0]
    eps = np.array([0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.300])
    e_bg_far = np.array([0.875300, 0.875428, 0.875458, 0.875232, 0.873854, 0.0, 0.0, 0.0, 0.0, 0.0])
    e_bg_0   = np.array([0.124510, 0.124530, 0.124530, 0.124476, 0.124188, 0.276774, 0.145363, 0.094612, 0.075447, 0.061260])
    e_p_0    = np.array([0.000031, 0.000007, 0.000002, 0.000048, 0.000321, 0.723226, 0.854637, 0.905388, 0.924553, 0.938740])
    normP    = np.array([0.04682, 0.04682, 0.04682, 0.04682, 0.04682, 0.49998, 0.49998, 0.49998, 0.49998, 0.49998])

    ax.plot(eps, e_bg_far, "o-", color="#1f77b4", lw=2.3, ms=6, label=r"遠方空間 $E_{\mathrm{bg}}(\mathrm{far})$ ($A=1..7$)" if ja else r"Far space $E_{\mathrm{bg}}(\mathrm{far})$ ($A=1..7$)")
    ax.plot(eps, e_p_0, "s-", color="#d62728", lw=2.3, ms=6, label=r"中心粒子 $E_P(0)$" if ja else r"Center particle $E_P(0)$")
    ax.plot(eps, e_bg_0, "^--", color="#2ca02c", lw=1.8, ms=6, label=r"中心背景 $E_{\mathrm{bg}}(0)$" if ja else r"Center bg $E_{\mathrm{bg}}(0)$")
    ax.plot(eps, normP, "d:", color="#9467bd", lw=2.0, ms=5.5, label=r"自己無撞着ノルム $\|P\|$ ($0.047 \to 1/2$)" if ja else r"Norm $\|P\|$ ($0.047 \to 1/2$)")
    ax.axvspan(0.270, 0.275, color="#ffcccc", alpha=0.65, label=r"臨界点 $\epsilon_c \in (0.270, 0.275)$" if ja else r"Critical $\epsilon_c \in (0.270, 0.275)$")
    ax.set_xlim(0.248, 0.302)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel(r"局在質量 P4 の初期振幅 $\epsilon$" if ja else r"Initial amplitude $\epsilon$ of localized mass P4", fontsize=10.5)
    ax.set_ylabel("エネルギー占有率 / ノルム" if ja else "Energy fraction / Norm", fontsize=10.5)
    ax.set_title(
        r"(a) 臨界質量 $\epsilon_c$ での遠方空間切断・特異点崩壊"
        if ja else
        r"(a) Critical Mass Collapse & Far-Space Disconnection at $\epsilon_c$",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(0.01, 0.50), fontsize=8.0, framealpha=0.92)
    ax.annotate(
        "nz: 576 → 16\n$E_{\\mathrm{bg}}(\\mathrm{far})=0$\n($A=0,1$ 6桁一致)" if ja else "nz: 576 → 16\n$E_{\\mathrm{bg}}(\\mathrm{far})=0$\n($A=0,1$ 6-digit match)",
        xy=(0.275, 0.50), xytext=(0.279, 0.26),
        arrowprops=dict(arrowstyle="->", color="#9467bd"), fontsize=8.4,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#9467bd", alpha=0.92)
    )

    # (b) スピン 2 (σ_x) グラビトン候補の普遍 1/sqrt(2) ロックと背景時計引きずり
    ax = axes[1]
    conds = [
        "真空 $H$\n(基準)" if ja else "Vacuum $H$\n(ref)",
        r"$\sigma_x$" + "\n(0.05)",
        r"$\sigma_x$" + "\n(0.10)",
        r"$\sigma_x$" + "\n(0.20)",
        "P4\n(0.10)",
        r"P4+$\sigma_x$" + "\n(0.10+0.10)",
    ]
    phi_bg = [1.00000, 1.06761, 0.90386, 0.65788, 0.93453, 0.86635]
    phi_pr = [np.nan,  0.65778, 0.65778, 0.65778, 0.70088, 0.70996]
    xx = np.arange(len(conds))
    w = 0.34
    ax.bar(xx - w / 2, phi_bg, width=w, color="#1f77b4", edgecolor="#222",
           label=r"背景時空時計 $\phi_{\mathrm{bg}}/\phi_0$" if ja else r"Background clock $\phi_{\mathrm{bg}}/\phi_0$")
    ax.bar(xx[1:] + w / 2, phi_pr[1:], width=w, color="#ff7f0e", edgecolor="#222",
           label=r"励起セクター時計 $\phi_{b23}/\phi_0$" if ja else r"Excitation clock $\phi_{b23}/\phi_0$")
    ax.axhline(0.65778, color="#d62728", ls="--", lw=1.4,
               label=r"$\sigma_x$ 普遍ロック水準 $0.6578$ ($d_1=1/\sqrt{2}$)" if ja else r"$\sigma_x$ universal lock $0.6578$ ($d_1=1/\sqrt{2}$)")
    for xi, vb, vp in zip(xx, phi_bg, phi_pr):
        ax.text(xi - w / 2, vb + 0.015, f"{vb:.3f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color="#1f77b4")
        if not np.isnan(vp):
            ax.text(xi + w / 2, vp + 0.015, f"{vp:.3f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color="#d95f02")
    ax.set_xticks(xx)
    ax.set_xticklabels(conds, fontsize=8.8)
    ax.set_ylim(0, 1.24)
    ax.set_ylabel(r"相対時計進行速度 $\phi / \phi_0$（$\tau=20..40$ 平均）" if ja else r"Relative clock rate $\phi / \phi_0$ ($\tau=20..40$ avg)", fontsize=10.2)
    ax.set_title(
        r"(b) スピン 2（$\sigma_x$）の普遍ロックと時空時計引きずり"
        if ja else
        r"(b) Spin-2 ($\sigma_x$) Universal Lock & Background Clock Dragging",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

    fig.suptitle(
        r"臨界質量での離散重力崩壊（$\epsilon_c \in (0.270, 0.275)$）とスピン 2（$\sigma_x$）による時空時計引きずり"
        if ja else
        r"Discrete Gravitational Collapse at $\epsilon_c \in (0.270, 0.275)$ and Spin-2 ($\sigma_x$) Clock Dragging",
        fontsize=12.3, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r5_collapse_and_spin2_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    os.makedirs("results", exist_ok=True)
    for lang in ["ja", "en"]:
        plot_fig1(lang)
        plot_fig2(lang)


if __name__ == "__main__":
    main()
