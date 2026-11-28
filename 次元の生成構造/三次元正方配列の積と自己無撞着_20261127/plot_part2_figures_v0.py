#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part II（統合仮説論文：重力・電磁気・QCD・ディラックの統一接続）用の図 2 種（日英計 4 枚）生成スクリプト。
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

    # (a) 三つ組閉じ残り: Maxwell (可換 U(1)_Q) vs Yang-Mills (非可換 O_8 & SU(3)) の O(ε^4) 収束
    ax = axes[0]
    eps = np.array([0.200, 0.100, 0.050, 0.025])
    ym_o8 = np.array([0.982373, 0.995565, 0.998889, 0.999722])
    ym_su3 = np.array([0.995540, 0.998888, 0.999722, 0.999931])
    ax.plot(eps, ym_su3, "o-", color="#d62728", lw=2.3, ms=7.5,
            label=r"非可換 $SU(3)$ 8 グルーオン（QCD ヤン＝ミルズ）" if ja else r"Non-Abelian $SU(3)$ 8 gluons (QCD Yang–Mills)")
    ax.plot(eps, ym_o8, "s--", color="#1f77b4", lw=2.2, ms=7.0,
            label=r"非可換八元数三つ組 $\mathrm{Re}(o_A o_B o_C)$" if ja else r"Non-Abelian octonion triple $\mathrm{Re}(o_A o_B o_C)$")
    ax.axhline(1.0, color="#2ca02c", ls=":", lw=1.8,
               label=r"厳密 $O(\varepsilon^4)$ ヤン＝ミルズ理論極限 $= 1.000$" if ja else r"Exact $O(\varepsilon^4)$ Yang–Mills limit $= 1.000$")
    ax.invert_xaxis()
    ax.set_ylim(0.978, 1.004)
    ax.set_xlabel(r"ループ結合スケール $\varepsilon$（右ほど連続極限 $\varepsilon \to 0$）" if ja else r"Loop scale $\varepsilon$ (continuum limit $\varepsilon \to 0$ to right)", fontsize=10.2)
    ax.set_ylabel(r"実測閉じ残り / 理論 $O(\varepsilon^4)$ ヤン＝ミルズ作用" if ja else r"Measured deficit / Exact $O(\varepsilon^4)$ YM action", fontsize=10.2)
    ax.set_title(
        r"(a) 同一三つ組からの Maxwell（可換）vs Yang–Mills（非可換）"
        if ja else
        r"(a) Abelian Maxwell vs. Non-Abelian Yang–Mills from One Triple",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.3, framealpha=0.92)
    ax.text(
        0.04, 0.25,
        ("可換 $U(1)_Q$ 平面 ($e_7$): 非可換閉じ残り $= 0.000\\times 10^0$\n"
         "$SU(3)$ ビアンキ恒等式 (全56組): 誤差 $0.000\\times 10^0$\n"
         "$[G_{\\mathrm{ext}}(SO(3))\\otimes I,\\, I\\otimes G_{\\mathrm{int}}]=0$: 誤差 $0.000\\times 10^0$")
        if ja else
        ("Abelian $U(1)_Q$ ($e_7$): non-Abelian deficit $= 0.000\\times 10^0$\n"
         "$SU(3)$ Bianchi identity (56 triples): err $0.000\\times 10^0$\n"
         "$[G_{\\mathrm{ext}}(SO(3))\\otimes I,\\, I\\otimes G_{\\mathrm{int}}]=0$: err $0.000\\times 10^0$"),
        transform=ax.transAxes, fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) 1 世代 8 状態の選択的ゲージ結合 (EM 光子 vs QCD グルーオン)
    ax = axes[1]
    sectors = [
        "$\\nu$ ($Q=0$)\n1重項" if ja else "$\\nu$ ($Q=0$)\nSinglet",
        "$\\bar{d}$ ($Q=+1/3$)\n3重項" if ja else "$\\bar{d}$ ($Q=+1/3$)\nTriplet",
        "$u$ ($Q=+2/3$)\n3重項" if ja else "$u$ ($Q=+2/3$)\nTriplet",
        "$e^+$ ($Q=+1$)\n1重項" if ja else "$e^+$ ($Q=+1$)\nSinglet",
    ]
    em_vals = [0.000000, 0.242733, 0.481878, 0.713899]
    qcd_vals = [0.000000, 0.504491, 0.580425, 0.000000]
    x = np.arange(4)
    w = 0.35
    ax.bar(x - w / 2, em_vals, width=w, color="#1f77b4", edgecolor="#222",
           label=r"電磁結合 $U(1)_Q$ 光子 ($\theta=0.73$)" if ja else r"EM coupling $U(1)_Q$ photon ($\theta=0.73$)")
    ax.bar(x + w / 2, qcd_vals, width=w, color="#d62728", edgecolor="#222",
           label=r"強い力 $SU(3)$ グルーオン ($\theta_g=0.51$)" if ja else r"Strong $SU(3)$ gluons ($\theta_g=0.51$)")
    for xi, ev, qv in zip(x, em_vals, qcd_vals):
        ax.text(xi - w / 2, ev + 0.015, f"{ev:.3f}", ha="center", va="bottom", fontsize=8.2, fontweight="bold", color="#1f77b4")
        ax.text(xi + w / 2, qv + 0.015, f"{qv:.3f}", ha="center", va="bottom", fontsize=8.2, fontweight="bold", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(sectors, fontsize=9.2)
    ax.set_ylim(0, 0.90)
    ax.set_ylabel(r"ゲージ変換応答ノルム $\|U\psi - \psi\|$" if ja else r"Gauge-transformation response $\|U\psi - \psi\|$", fontsize=10.5)
    ax.set_title(
        r"(b) 1 世代 8 状態の標準模型そのままの選択的ゲージ結合"
        if ja else
        r"(b) Exact Standard-Model Selective Gauge Couplings (1 Generation)",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.3, framealpha=0.92)
    ax.text(
        0.04, 0.68,
        ("ニュートリノ $\\nu$: EM=0, QCD=0\n"
         "陽電子 $e^+$: EM 最大, QCD=0\n"
         "クォーク $\\bar{d}, u$: EM・QCD 両結合")
        if ja else
        ("Neutrino $\\nu$: EM=0, QCD=0\n"
         "Positron $e^+$: max EM, QCD=0\n"
         "Quarks $\\bar{d}, u$: both EM & QCD"),
        transform=ax.transAxes, fontsize=8.1,
        bbox=dict(boxstyle="round,pad=0.28", fc="#fff8f0", ec="#d62728", alpha=0.93)
    )

    fig.suptitle(
        r"Part II 機械精度検証：Maxwell / Yang–Mills 読み分けと 1 世代フェルミオンの選択的ゲージ結合"
        if ja else
        r"Part II Machine-Precision Verification: Maxwell vs. Yang–Mills Splitting and Selective Fermion Couplings",
        fontsize=12.2, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p2_ym_maxwell_and_selective_coupling_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, ax = plt.subplots(figsize=(13.0, 6.2), dpi=200)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc, ec, title_c="#111"):
        rect = plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, lw=2.0, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h - 0.28, title, ha="center", va="top", fontsize=10.0, fontweight="bold", color=title_c, zorder=3)
        ax.text(x + w / 2, y + h - 0.62, body, ha="center", va="top", fontsize=8.4, color="#222", zorder=3, linespacing=1.35)

    # Top axiom box
    box(2.2, 4.95, 5.6, 1.05,
        "単一公理：実 3 階テンソル三角形ループ積  Tri(T)_{ijk} = Σ T_{ipq} T_{qjr} T_{rpk}" if ja else "Single Axiom: Real 3-Tensor Triangle-Loop Product  Tri(T)_{ijk} = Σ T_{ipq} T_{qjr} T_{rpk}",
        "直積準同型 Tri(A $\\otimes$ B) = Tri(A) $\\otimes$ Tri(B) (誤差 2.86e-16)  →  $O_8^{\\mathrm{ext}} \\otimes O_8^{\\mathrm{int}}$ 分離構造（Coleman–Mandula 回避）" if ja else "Homomorphism Tri(A $\\otimes$ B) = Tri(A) $\\otimes$ Tri(B) (err 2.86e-16)  →  $O_8^{\\mathrm{ext}} \\otimes O_8^{\\mathrm{int}}$ Separation (Coleman–Mandula safe)",
        "#f0f4ff", "#1f77b4", "#0b3d91")

    # Left: External Spacetime & Gravity
    box(0.3, 2.15, 4.3, 2.25,
        "外側セクター O_8^ext（外部時空・重力・物差し I）" if ja else "Outer Sector O_8^ext (External Spacetime & Gravity)",
        ("・四元数 3 空間 {e1,e2,e3}: SO(3) スピン 0/1/2 選別\n"
         "・球対称 4 クラス空間計量歪み (+32.5% / +7.2% / +1.2%)\n"
         "・離散曲率 (視野角 ρ) 勾配・臨界崩壊 ε_c ∈ (0.270, 0.275)\n"
         "・ディラック空間微分 D = Σ L_a p_a + m (D² = -|p|² + m²)\n"
         "・[K7] 物差し仮説 1: (1+3)+(1+3) 空間骨格")
        if ja else
        ("• Quaternionic 3-space {e1,e2,e3}: SO(3) spin 0/1/2 split\n"
         "• 7-digit isotropic 4-class metric distortion (+32.5% / +7.2% / +1.2%)\n"
         "• Discrete curvature (ρ) gradient & critical collapse ε_c\n"
         "• Dirac spatial operator D = Σ L_a p_a + m (D² = -|p|² + m²)\n"
         "• [K7] Ruler Hypothesis 1: (1+3)+(1+3) spatial skeleton"),
        "#fff6f0", "#d62728", "#8b0000")

    # Right: Internal Gauge & Fermions
    box(5.4, 2.15, 4.3, 2.25,
        "内側セクター O_8^int / M_8（内部ゲージ・電磁・QCD・時計）" if ja else "Inner Sector O_8^int / M_8 (Gauge, EM, QCD & Clock)",
        ("・e_7 固定: SU(3) 8 グルーオン × U(1)_Q 電荷 (可換誤差 0)\n"
         "・可換 e_7 平面 → 電磁気 (Maxwell, 閉じ残り 0)\n"
         "・非可換 SU(3) 方向 → 強い力 (QCD ヤン＝ミルズ, O(ε⁴))\n"
         "・CAR パウリ排他則・1 世代 8 状態 (ν, d_bar, u, e+) 選択結合\n"
         "・[K7] 物差し仮説 2: M_8 複素スカラー時計 (スピン 2 引きずり)")
        if ja else
        ("• Fixed e_7: SU(3) 8 gluons × U(1)_Q charge (comm err 0)\n"
         "• Abelian e_7 plane → Electromagnetism (Maxwell, deficit 0)\n"
         "• Non-Abelian SU(3) → Strong force (QCD Yang–Mills, O(ε⁴))\n"
         "• CAR Pauli exclusion & selective 1-gen couplings (ν, d_bar, u, e+)\n"
         "• [K7] Ruler Hypothesis 2: M_8 scalar clock (spin-2 dragging)"),
        "#f2fbf2", "#2ca02c", "#145a32")

    # Bottom: Coupling via Layer 3 & Triality
    box(1.3, 0.25, 7.4, 1.35,
        "両セクターの統一結合機構（第 3 層単一時計ロック × Triality × 階層入れ子 $O_8^{\\otimes K}$）" if ja else "Unified Coupling Mechanism (Layer-3 Single-Clock Lock × Triality × Nesting $O_8^{\\otimes K}$)",
        ("・$[G_{\\mathrm{ext}} \\otimes I,\\; I \\otimes G_{\\mathrm{int}}] = 0$ (誤差 0.000e+00) で独立対称性を保持しつつ、第 3 層単一複素構造ロックで単一時計を共有\n"
         "・同一三つ組 $\\mathrm{Re}(o_A o_B o_C)$ の読み分け: 3つともつなぎ → Maxwell/Yang–Mills、両端が物質 → Dirac ($\\bar{o}_A^T L(o_B) o_C$)\n"
         "・大格子 $O_8^{\\otimes K}$ 展開 → $1/r^2$ クーロン力・$1/r$ 重力・光円錐速度 $c$・連続一般相対論／標準模型への次段階検証")
        if ja else
        ("• Preserves $[G_{\\mathrm{ext}} \\otimes I,\\; I \\otimes G_{\\mathrm{int}}] = 0$ (err 0.000e+00) while sharing one global clock via Layer-3 single-complex lock\n"
         "• One triple $\\mathrm{Re}(o_A o_B o_C)$: 3 links → Maxwell/Yang–Mills; 2 matter + 1 link → Dirac ($\\bar{o}_A^T L(o_B) o_C$)\n"
         "• Large-lattice $O_8^{\\otimes K}$ expansion → next-stage test of $1/r^2$ Coulomb, $1/r$ gravity, light-cone $c$, and continuum GR/SM"),
        "#f9f4ff", "#9467bd", "#4a235a")

    # Arrows
    ax.annotate("", xy=(2.45, 4.42), xytext=(4.2, 4.93), arrowprops=dict(arrowstyle="->", lw=2.0, color="#1f77b4"))
    ax.annotate("", xy=(7.55, 4.42), xytext=(5.8, 4.93), arrowprops=dict(arrowstyle="->", lw=2.0, color="#1f77b4"))
    ax.annotate("", xy=(3.5, 1.62), xytext=(2.45, 2.13), arrowprops=dict(arrowstyle="->", lw=2.0, color="#9467bd"))
    ax.annotate("", xy=(6.5, 1.62), xytext=(7.55, 2.13), arrowprops=dict(arrowstyle="->", lw=2.0, color="#9467bd"))
    ax.annotate("", xy=(5.38, 3.28), xytext=(4.62, 3.28), arrowprops=dict(arrowstyle="<->", lw=1.8, color="#555", ls="--"))
    ax.text(5.0, 3.45, "[G_ext,G_int]=0", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color="#333")

    fig.suptitle(
        "Part II 統一見取図：外部時空（重力）× 内部ゲージ（電磁・QCD）× 物質（ディラック）の分離と結合"
        if ja else
        "Part II Unified Architecture: Separation and Coupling of Gravity, Electromagnetism, QCD, and Dirac Matter",
        fontsize=12.0, fontweight="bold", y=0.98
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = f"results/fig_p2_unified_architecture_{lang}.png"
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
