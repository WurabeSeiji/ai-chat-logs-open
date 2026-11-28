#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part VII（例外ジョルダン代数 J_3(O) の Tri(T) 結晶性・F_4 -> SO(8) トライアリティ 3 世代構造・
3 世代非線形階層増幅・素の結合定数幾何統一比）用の図 2 種（日英計 4 枚）生成スクリプト。
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

    # (a) J_3(O) 対称群連鎖 F_4(52D) -> SO(8)(28D) -> SO(7)_1(21D) -> G_2(14D) -> SU(3)_C(8D)
    ax = axes[0]
    groups = [
        r"$\mathrm{Der}(J_3(\mathbb{O}))$" + "\n" + r"$F_4$",
        r"対角固定" + "\n" + r"$SO(8)$" if ja else "Diag-fix\n$SO(8)$",
        r"$O_1$ の $e_0$ 固定" + "\n" + r"$SO(7)_1$" if ja else "$O_1$ $e_0$-fix\n$SO(7)_1$",
        r"3世代 $e_0$ 固定" + "\n" + r"$G_2$" if ja else "All $e_0$-fix\n$G_2$",
        r"さらに $e_1$ 固定" + "\n" + r"$SU(3)_C$" if ja else "Fix $e_1$\n$SU(3)_C$",
    ]
    dims = [52, 28, 21, 14, 8]
    cols = ["#9467bd", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    x = np.arange(len(groups))
    ax.bar(x, dims, width=0.58, color=cols, edgecolor="#222")
    for xi, d in zip(x, dims):
        ax.text(xi, d + 1.2, f"{d}D", ha="center", va="bottom", fontsize=9.2, fontweight="bold", color="#222")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=8.4)
    ax.set_ylim(0, 62)
    ax.set_ylabel("SVD 核空間から抽出された厳密次元" if ja else "Exact dimension from SVD null space", fontsize=10.0)
    ax.set_title(
        r"(a) $J_3(\mathbb{O})$ (27D) の厳密対称群連鎖とトライアリティ 3 世代"
        if ja else
        r"(a) Exact Symmetry Chain of $J_3(\mathbb{O})$ (27D) & Triality 3 Generations",
        fontsize=10.1, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(
        0.34, 0.58,
        ("厳密検証結果 (誤差 $\\leq 2\\times 10^{-15}$):\n"
         "• $J_{\\mathrm{off}}$ (24D): $\\mathrm{Tri}(J_{\\mathrm{off}}) = -\\frac{1}{64}J_{\\mathrm{off}}$\n"
         "  （不動点残差 `0.000e+00`）\n"
         "• $SO(7)_1$ 固定次元 $(O_1,O_2,O_3) = (1,0,0)$\n"
         "  $\\Rightarrow 8_v \\oplus 8_s \\oplus 8_c$ トライアリティ\n"
         "• $SU(3)_C$ カシミール: 全 3 世代\n"
         "  $[0,0, 8/9\\times 6]$ (レプトン2+クォーク6)")
        if ja else
        ("Exact verification (err $\\leq 2\\times 10^{-15}$):\n"
         "• $J_{\\mathrm{off}}$ (24D): $\\mathrm{Tri}(J_{\\mathrm{off}}) = -\\frac{1}{64}J_{\\mathrm{off}}$\n"
         "  (fixed-point residual 0.000e+00)\n"
         "• $SO(7)_1$ fixed dims $(O_1,O_2,O_3)=(1,0,0)$\n"
         "  $\\Rightarrow 8_v \\oplus 8_s \\oplus 8_c$ triality\n"
         "• $SU(3)_C$ Casimir: all 3 generations\n"
         "  $[0,0, 8/9\\times 6]$ (2 leptons + 6 quarks)"),
        transform=ax.transAxes, fontsize=7.6,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) 3 世代結晶 J_off 上での 2 段階対称性破れによる指数的 3 世代階層増幅 (対数軸)
    ax = axes[1]
    labels = [
        "対称結晶\n(0, 0)" if ja else "Symmetric\n(0, 0)",
        "微小破れ\n(0.01, 0.005)" if ja else "Tiny break\n(0.01, 0.005)",
        "2段階破れ\n(0.14, 0.08)" if ja else "2-step break\n(0.14, 0.08)",
        "強破れ\n(0.18, 0.12)" if ja else "Strong break\n(0.18, 0.12)",
    ]
    e1 = [0.333333, 0.775423, 0.964830, 0.992123]
    e2 = [0.333333, 0.221966, 0.033242, 0.007413]
    e3 = [0.333333, 0.002610, 0.001928, 0.000464]
    x = np.arange(len(labels))
    w = 0.25
    ax.bar(x - w, e1, width=w, color="#d62728", edgecolor="#222", label=r"第 1 世代 $E(O_1)$" if ja else r"Gen 1 $E(O_1)$")
    ax.bar(x,     e2, width=w, color="#ff7f0e", edgecolor="#222", label=r"第 2 世代 $E(O_2)$" if ja else r"Gen 2 $E(O_2)$")
    ax.bar(x + w, e3, width=w, color="#1f77b4", edgecolor="#222", label=r"第 3 世代 $E(O_3)$" if ja else r"Gen 3 $E(O_3)$")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 3.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.4)
    ax.set_ylabel(r"各世代エネルギー占有率 $E(O_g)$（$\tau=40$, 対数軸）" if ja else r"Generation energy share $E(O_g)$ ($\tau=40$, log)", fontsize=9.8)
    ax.set_title(
        r"(b) $J_{\mathrm{off}}$ 結晶上の $\mathrm{Tri}(T)$ 非線形 3 世代階層増幅"
        if ja else
        r"(b) Nonlinear 3-Generation Hierarchy Amplification on $J_{\mathrm{off}}$",
        fontsize=10.1, fontweight="bold"
    )
    ax.grid(True, which="both", axis="y", alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.1, framealpha=0.94)
    ax.text(
        0.46, 0.08,
        ("比 $E_1:E_2:E_3$ の自発増幅:\n"
         "(0.14, 0.08) $\\Rightarrow 500 : 17.2 : 1$\n"
         "(0.18, 0.12) $\\Rightarrow 2139 : 16.0 : 1$\n"
         "※ 小出比 $Q$ は初期比に依存\n"
         "  （普遍定数 $2/3$ へは非ロック）")
        if ja else
        ("Spontaneous hierarchy ratio:\n"
         "(0.14, 0.08) $\\Rightarrow 500 : 17.2 : 1$\n"
         "(0.18, 0.12) $\\Rightarrow 2139 : 16.0 : 1$\n"
         "* Koide Q depends on initial ratio\n"
         "  (no universal lock to 2/3)"),
        transform=ax.transAxes, fontsize=7.7,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#d62728", alpha=0.94)
    )

    fig.suptitle(
        r"Part VII 検証 (1)：例外ジョルダン代数 $J_3(\mathbb{O})$ の 3 世代トライアリティと $\mathrm{Tri}(T)$ 階層増幅"
        if ja else
        r"Part VII Verification (1): $J_3(\mathbb{O})$ 3-Generation Triality & Nonlinear Hierarchy Amplification",
        fontsize=11.6, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p7_j3o_triality_and_hierarchy_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), dpi=200)

    # (a) 64D 1 世代 & 27D 3 世代の素の結合定数幾何統一比 (無パラメータ)
    ax = axes[0]
    cats = [
        r"$SU(3)_C$" + "\n" + r"$\mathrm{Tr}_{64}(T_3^{C\,2})$",
        r"$SU(2)_L$" + "\n" + r"$\mathrm{Tr}_{64}(T_3^{L\,2})$",
        r"$U(1)_Y$ (GUT)" + "\n" + r"$\mathrm{Tr}_{64}(Y_1^2)$",
        r"$U(1)_Y$ (素)" + "\n" + r"$\mathrm{Tr}_{64}(Y^2)$" if ja else r"$U(1)_Y$ (raw)" + "\n" + r"$\mathrm{Tr}_{64}(Y^2)$",
        r"$U(1)_{\mathrm{em}}$" + "\n" + r"$\mathrm{Tr}_{64}(Q_{\mathrm{em}}^2)$",
    ]
    vals = [8.0, 8.0, 8.0, 40.0 / 3.0, 64.0 / 3.0]
    cols = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd"]
    x = np.arange(len(cats))
    ax.bar(x, vals, width=0.56, color=cols, edgecolor="#222")
    ax.axhline(8.0, color="#333", ls="--", lw=1.4, label=r"厳密ゲージ統一ライン $\mathrm{Tr}=8$ ($\alpha_3:\alpha_2:\alpha_1 = 1:1:1$)" if ja else r"Exact unification line $\mathrm{Tr}=8$ ($\alpha_3:\alpha_2:\alpha_1 = 1:1:1$)")
    for xi, v, s in zip(x, vals, ["8", "8", "8", "40/3", "64/3"]):
        ax.text(xi, v + 0.55, s, ha="center", va="bottom", fontsize=9.0, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=8.3)
    ax.set_ylim(0, 26)
    ax.set_ylabel("64D 表現上の無パラメータ幾何トレース" if ja else "Parameter-free geometric trace on 64D", fontsize=9.8)
    ax.set_title(
        r"(a) 64D 表現上の素のゲージ結合定数幾何統一（$1:1:1$）"
        if ja else
        r"(a) Parameter-Free Bare Gauge Unification on 64D ($1:1:1$)",
        fontsize=10.2, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.0, framealpha=0.94)

    # (b) 素の幾何比 (高エネルギー) vs 低エネルギー実測値 (M_Z) の境界比較
    ax = axes[1]
    obs = [
        r"$\sin^2\theta_W$",
        r"$\alpha_{\mathrm{em}} / \alpha_s$",
        r"$\alpha_2 / \alpha_s$",
    ]
    bare_vals = [0.3750, 0.3750, 1.0000]
    exp_vals = [0.2312, 0.0663, 0.2869]  # α_em(M_Z)=1/127.9, α_2=α_em/sin^2θ_W=0.03383, α_s=0.1179 -> α_2/α_s=0.2869
    x = np.arange(len(obs))
    w = 0.34
    ax.bar(x - w / 2, bare_vals, width=w, color="#1f77b4", edgecolor="#222",
           label=r"本モデルの素の幾何値（無パラメータ）" if ja else r"Bare geometric value (zero free params)")
    ax.bar(x + w / 2, exp_vals, width=w, color="#2ca02c", edgecolor="#222",
           label=r"低エネルギー実測値（$M_Z$ スケール）" if ja else r"Low-energy experiment ($M_Z$ scale)")
    for xi, vb, ve in zip(x, bare_vals, exp_vals):
        ax.text(xi - w / 2, vb + 0.025, f"{vb:.4f}", ha="center", va="bottom", fontsize=8.4, fontweight="bold", color="#1f77b4")
        ax.text(xi + w / 2, ve + 0.025, f"{ve:.4f}", ha="center", va="bottom", fontsize=8.4, fontweight="bold", color="#1b7a1b")
    ax.set_xticks(x)
    ax.set_xticklabels(obs, fontsize=10.0)
    ax.set_ylim(0, 1.22)
    ax.set_ylabel("無次元結合定数比" if ja else "Dimensionless coupling ratio", fontsize=10.0)
    ax.set_title(
        r"(b) 素の幾何統一値 vs $M_Z$ 実測値（RG 走行境界の明示）"
        if ja else
        r"(b) Bare Unification Values vs. $M_Z$ Experiment (RG Boundary)",
        fontsize=10.2, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.2, framealpha=0.94)
    ax.text(
        0.04, 0.56,
        ("反証・境界の明示:\n"
         "• 素の幾何値: $\\alpha_3:\\alpha_2:\\alpha_1 = 1:1:1$,\n"
         "  $\\sin^2\\theta_W^{\\mathrm{bare}} = \\alpha_{\\mathrm{em}}/\\alpha_s = 3/8$\n"
         "• 低エネルギー実測 ($\\alpha \\approx 1/137.036$,\n"
         "  $\\alpha_s \\approx 0.1179$) との差は RG 走行が必須")
        if ja else
        ("Explicit falsification boundary:\n"
         "• Bare geometry: $\\alpha_3:\\alpha_2:\\alpha_1 = 1:1:1$,\n"
         "  $\\sin^2\\theta_W^{\\mathrm{bare}} = \\alpha_{\\mathrm{em}}/\\alpha_s = 3/8$\n"
         "• Gap to low-energy data ($\\alpha\\approx 1/137.036$,\n"
         "  $\\alpha_s\\approx 0.1179$) requires RG running"),
        transform=ax.transAxes, fontsize=7.8,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#d62728", alpha=0.94)
    )

    fig.suptitle(
        r"Part VII 検証 (2)：素の結合定数幾何統一比（$\alpha_3:\alpha_2:\alpha_1=1:1:1$）と低エネルギー RG 境界"
        if ja else
        r"Part VII Verification (2): Bare Gauge Coupling Unification ($1:1:1$) & Low-Energy RG Boundary",
        fontsize=11.6, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p7_coupling_unification_and_rg_boundary_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    for lang in ("ja", "en"):
        plot_fig1(lang)
        plot_fig2(lang)
