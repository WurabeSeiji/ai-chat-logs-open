#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part III（時空・内部空間読出し接続仮説：対数変換・過去光円錐・分解能ゲージ・1+3+4）用の
図 2 種（日英計 4 枚）生成スクリプト。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


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

    # (a) 対数変換による局所曲率半径 R_loc/R_0 = exp(-Δφ) の縮小と臨界事象の地平面
    ax = axes[0]
    eps = np.array([0.000, 0.050, 0.100, 0.160, 0.220, 0.270, 0.275])
    r_loc = np.array([1.000000, 0.977051, 0.922917, 0.868798, 0.845041, 0.832794, 0.000000])
    dphi = np.array([0.000000, 0.023217, 0.080216, 0.140645, 0.168370, 0.182969, np.nan])

    ax.plot(eps[:-1], r_loc[:-1], "o-", color="#1f77b4", lw=2.4, ms=7.5,
            label=r"局所曲率半径比 $R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$" if ja else r"Local curvature radius $R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$")
    ax.plot([0.270, 0.275], [0.832794, 0.0], "--", color="#d62728", lw=2.0)
    ax.plot([0.275], [0.0], "X", color="#d62728", ms=10.5,
            label=r"臨界地平面切断 ($\epsilon_c \in (0.270, 0.275)$)" if ja else r"Critical horizon cut ($\epsilon_c \in (0.270, 0.275)$)")
    ax.axvline(0.2725, color="#d62728", ls=":", lw=1.6, alpha=0.8)

    ax2 = ax.twinx()
    ax2.plot(eps[:-1], dphi[:-1], "s--", color="#ff7f0e", lw=2.0, ms=6.2,
             label=r"対数距離ズレ $\Delta\phi = \frac{1}{2}\ln(E_{00}/E_{\mathrm{far}})$" if ja else r"Log distance shift $\Delta\phi = \frac{1}{2}\ln(E_{00}/E_{\mathrm{far}})$")
    ax2.set_ylabel(r"対数変換ズレ $\Delta\phi$" if ja else r"Log shift $\Delta\phi$", fontsize=10.2, color="#d95f02")
    ax2.tick_params(axis="y", labelcolor="#d95f02")
    ax2.set_ylim(-0.02, 0.24)

    ax.set_ylim(-0.06, 1.08)
    ax.set_xlim(-0.01, 0.295)
    ax.set_xlabel(r"局在質量 P4 の振幅 $\epsilon$ ($\tau=40$)" if ja else r"Localized mass P4 amplitude $\epsilon$ ($\tau=40$)", fontsize=10.2)
    ax.set_ylabel(r"局所曲率半径比 $R_{\mathrm{loc}}/R_0$" if ja else r"Local curvature radius ratio $R_{\mathrm{loc}}/R_0$", fontsize=10.2, color="#1f77b4")
    ax.tick_params(axis="y", labelcolor="#1f77b4")
    ax.set_title(
        r"(a) 対数変換 $O=e^{\ln A+n\theta}$ による局所曲率半径縮小と地平面"
        if ja else
        r"(a) Local Curvature Shrinking & Horizon Cut via Log Map",
        fontsize=10.6, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center left", bbox_to_anchor=(0.01, 0.54), fontsize=7.9, framealpha=0.94)

    # (b) 自分を原点とする過去光円錐（扇型座標系）における真空等方性 vs 遠方観測者の潮汐分裂
    ax = axes[1]
    labels = [
        "真空 $H$\n(観測者 $A=0$)\n全7方向" if ja else "Vacuum $H$\n(Obs $A=0$)\nAll 7 dirs",
        "真空 $H$\n(観測者 $A=1$)\n全7方向" if ja else "Vacuum $H$\n(Obs $A=1$)\nAll 7 dirs",
        "質量直上\n(観測者 $A=0$)\n全7方向" if ja else "On mass\n(Obs $A=0$)\nAll 7 dirs",
        "遠方観測者\n($A=1\\to$質量$B=0$)\n動径方向" if ja else "Far obs ($A=1$)\nTo mass ($B=0$)\nRadial dir",
        "遠方観測者\n($A=1\\to B=2..7$)\n接線6方向" if ja else "Far obs ($A=1$)\nTo $B=2..7$\n6 tang. dirs",
    ]
    vals = [0.000000, 0.000000, 0.134551, -0.006094, -0.034709]
    colors = ["#7f7f7f", "#7f7f7f", "#2ca02c", "#d62728", "#1f77b4"]
    x = np.arange(len(labels))
    bars = ax.bar(x, vals, width=0.56, color=colors, edgecolor="#222")
    ax.axhline(0.0, color="#333", lw=1.1)
    for xi, v in zip(x, vals):
        va = "bottom" if v >= 0 else "top"
        off = 0.005 if v >= 0 else -0.006
        ax.text(xi, v + off, f"{v:+.4f}", ha="center", va=va, fontsize=8.5, fontweight="bold")

    ax.annotate(
        r"潮汐分裂 $\Delta D_{\mathrm{tidal}} = +0.0286$" if ja else r"Tidal split $\Delta D_{\mathrm{tidal}} = +0.0286$",
        xy=(3.5, -0.020), xytext=(2.65, -0.062),
        fontsize=8.8, fontweight="bold", color="#b22222",
        bbox=dict(boxstyle="round,pad=0.25", fc="#fff4f4", ec="#d62728"),
        arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.4)
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.3)
    ax.set_ylim(-0.082, 0.172)
    ax.set_ylabel(r"対数過去距離 $D_{\log}(A_{\mathrm{obs}}\to B)$" if ja else r"Log past distance $D_{\log}(A_{\mathrm{obs}}\to B)$", fontsize=10.2)
    ax.set_title(
        r"(b) 自分中心の過去光円錐：真空等方性（等価原理）vs 遠方潮汐分裂"
        if ja else
        r"(b) Observer-Centered Past Light Cone: Vacuum Isotropy vs. Tidal Split",
        fontsize=10.6, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        r"Part III 検証 (1)：対数変換による局所曲率半径縮小・地平面と、過去光円錐上の潮汐分裂"
        if ja else
        r"Part III Verification (1): Curvature Shrinking via Log Map and Tidal Splitting on Past Light Cone",
        fontsize=12.0, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p3_log_curvature_and_tidal_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.5), dpi=200)

    # (a) 分解能ゲージ：Z_2^3 ローパス窓 K_cut と Hurwitz 4 段階粒子の可視化境界
    ax = axes[0]
    k_cuts = np.array([0, 1, 2, 3])
    o8_c = np.array([1.000000, 1.000000, 1.000000, 1.000000])
    h4_c = np.array([0.900000, 1.000000, 1.000000, 1.000000])
    c2_c = np.array([0.749250, 0.916417, 1.000000, 1.000000])
    r1_c = np.array([0.556705, 0.746689, 0.936672, 1.000000])

    ax.plot(k_cuts, o8_c, "o-", color="#7f7f7f", lw=2.2, ms=7.0,
            label=r"$\mathbb{O}_8$ ($m=0$: 背景, $k_{\max}=0$)" if ja else r"$\mathbb{O}_8$ ($m=0$: bg, $k_{\max}=0$)")
    ax.plot(k_cuts, h4_c, "s-", color="#2ca02c", lw=2.3, ms=7.2,
            label=r"$\mathbb{H}_4$ ($m=1$: P5, $k_{\max}=1$)" if ja else r"$\mathbb{H}_4$ ($m=1$: P5, $k_{\max}=1$)")
    ax.plot(k_cuts, c2_c, "^-", color="#1f77b4", lw=2.3, ms=7.5,
            label=r"$\mathbb{C}_2$ ($m=2$, $k_{\max}=2$)" if ja else r"$\mathbb{C}_2$ ($m=2$, $k_{\max}=2$)")
    ax.plot(k_cuts, r1_c, "D-", color="#d62728", lw=2.4, ms=7.2,
            label=r"$\mathbb{R}_1$ ($m=3$: P4, $k_{\max}=3$)" if ja else r"$\mathbb{R}_1$ ($m=3$: P4, $k_{\max}=3$)")

    ax.set_xticks(k_cuts)
    ax.set_xticklabels(
        ["$K_{\\mathrm{cut}}=0$\n(DCのみ)", "$K_{\\mathrm{cut}}\\leq 1$", "$K_{\\mathrm{cut}}\\leq 2$", "$K_{\\mathrm{cut}}\\leq 3$\n(全可視)"]
        if ja else
        ["$K_{\\mathrm{cut}}=0$\n(DC only)", "$K_{\\mathrm{cut}}\\leq 1$", "$K_{\\mathrm{cut}}\\leq 2$", "$K_{\\mathrm{cut}}\\leq 3$\n(Full)"]
    )
    ax.set_ylim(0.32, 1.06)
    ax.set_xlabel(r"観測者の $\mathbb{Z}_2^3$ 倍音ローパス分解能窓 $K_{\mathrm{cut}}$" if ja else r"Observer $\mathbb{Z}_2^3$ harmonic low-pass resolution window $K_{\mathrm{cut}}$", fontsize=10.0)
    ax.set_ylabel(r"可視エネルギー比（累積ローパス通過率）" if ja else r"Visible energy fraction (cumulative low-pass)", fontsize=10.2)
    ax.set_title(
        r"(a) 分解能窓 $K_{\mathrm{cut}}=m$ による Hurwitz 粒子の可視化境界"
        if ja else
        r"(a) Hurwitz Particle Visibility Boundary at $K_{\mathrm{cut}} = m$",
        fontsize=10.6, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.3, framealpha=0.92)
    ax.text(
        0.03, 0.04,
        ("辺長閾値 $\\delta_{\\mathrm{res}} = 10^{-2} \\sim 10^{-20}$（18桁掃引）:\n"
         "・見える三角形数 $= 64$ 個（不変）\n"
         "・視野角 $\\langle\\rho\\rangle_w = 0.7767719701$（10桁不変）\n"
         "→ 分解能はゲージ、視野角はゲージ不変量")
        if ja else
        ("Edge threshold $\\delta_{\\mathrm{res}} = 10^{-2}..10^{-20}$ (18 orders):\n"
         "• Visible triangles $= 64$ (invariant)\n"
         "• View angle $\\langle\\rho\\rangle_w = 0.7767719701$ (10-digit inv.)\n"
         "→ Resolution = gauge, view angle = invariant"),
        transform=ax.transAxes, fontsize=7.9,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) O_8^ext (1+3+4) x O_8^int (2+6) -> (x,y,z,t,R,Q1,Q2,Q3) 対応アーキテクチャ
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def box(x0, y0, w, h, text, fc, ec, fs=8.3, fw="normal"):
        p = FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.22", fc=fc, ec=ec, lw=1.7)
        ax.add_patch(p)
        ax.text(x0 + w / 2, y0 + h / 2, text, ha="center", va="center", fontsize=fs, fontweight=fw)

    box(0.4, 8.2, 9.2, 1.35,
        ("単一公理 $\\mathrm{Tri}(T)$・状態空間 $O_8^{\\mathrm{ext}} \\otimes O_8^{\\mathrm{int}}$（64実自由度）\n"
         "真空 $H$: 全7種 $(1+3)+(1+3)$ 分割で $|a|^2=|b|^2=0.5$（誤差 $3.3\\times 10^{-16}$）")
        if ja else
        ("Single Axiom $\\mathrm{Tri}(T)$: State Space $O_8^{\\mathrm{ext}} \\otimes O_8^{\\mathrm{int}}$ (64 real DOF)\n"
         "Vacuum $H$: $|a|^2=|b|^2=0.5$ for all 7 $(1+3)+(1+3)$ splits (err $3.3\\times 10^{-16}$)"),
        "#f0f4ff", "#1f4e79", fs=8.5, fw="bold")

    box(0.4, 3.8, 4.35, 3.7,
        ("外側 $O_8^{\\mathrm{ext}}$：1 + 3 + 4 分割\n\n"
         "・1 ($e_0$): $\\ln A$・曲率 $R$・時間 $t$\n"
         "・3 ($e_{1..3}$): 過去光円錐の方向 $n$\n"
         "   $SO(3)$・重力潮汐 ($+0.0286$)\n"
         "・4 ($e_{4..7}$): 純粋 $SU(2)$ (3D)\n"
         "   3空間作用誤差 $8.3\\times 10^{-17}$\n"
         "   スピン $1/2$ 誤差 $5.8\\times 10^{-15}$\n"
         "   → 弱アイソスピン / ヒッグス4成分")
        if ja else
        ("Outer $O_8^{\\mathrm{ext}}$: 1 + 3 + 4 Split\n\n"
         "• 1 ($e_0$): $\\ln A$, curvature $R$, time $t$\n"
         "• 3 ($e_{1..3}$): Past light-cone dir $n$\n"
         "   $SO(3)$, tidal split ($+0.0286$)\n"
         "• 4 ($e_{4..7}$): Pure $SU(2)$ (3D)\n"
         "   3-space err $8.3\\times 10^{-17}$\n"
         "   Spin-1/2 err $5.8\\times 10^{-15}$\n"
         "   → Weak isospin / Higgs 4-comp."),
        "#eef9f0", "#2ca02c", fs=8.0)

    box(5.25, 3.8, 4.35, 3.7,
        ("内側 $O_8^{\\mathrm{int}}$：2 + 6 分割 ($e_7$ 固定)\n\n"
         "・$SU(3)_C$ (8グルーオン) × $U(1)_Q$\n"
         "・CAR 1世代8状態 ($1:3:3:1$)\n"
         "   色・電荷 $(Q_1, Q_2, Q_3)$\n"
         "・凝縮相 $M_8$: 複素時計 $J$\n"
         "   カイラル非対称分裂 ($-1.41\\%$)\n"
         "・外側 $SO(3)$ と厳密可換\n"
         "   全27組誤差 $0.000\\times 10^0$")
        if ja else
        ("Inner $O_8^{\\mathrm{int}}$: 2 + 6 Split ($e_7$ fixed)\n\n"
         "• $SU(3)_C$ (8 gluons) × $U(1)_Q$\n"
         "• CAR 1-gen 8 states ($1:3:3:1$)\n"
         "   Color & charge $(Q_1, Q_2, Q_3)$\n"
         "• Condensate $M_8$: complex clock $J$\n"
         "   Chiral split ($-1.41\\%$)\n"
         "• Strictly commutes with $SO(3)$\n"
         "   All 27 pairs err $0.000\\times 10^0$"),
        "#fff5ee", "#d95f02", fs=8.0)

    box(0.4, 0.45, 9.2, 2.65,
        ("読出し写像：対数変換 $O = e^{\\ln A + n\\theta}$ ＋ 分解能ゲージ（$\\delta_{\\mathrm{res}}, K_{\\mathrm{cut}}$）\n"
         "積の世界 ($\\mathrm{Tri}$) $\\leftrightarrow$ 加法の時空・共役双対 ($\\cosh\\phi$ ローレンツ vs $\\cos\\theta$ ボルン)\n"
         "・見える側（分解能以上・過去光円錐）: 時空 $(x, y, z, t)$・潮汐・粒子軌道\n"
         "・見えない側（分解能以下・投影軸 $R_0$）: 曲率 $R$・内部量子数 $(Q_1, Q_2, Q_3)$・$SU(2)_L$")
        if ja else
        ("Readout Map: Log Transform $O = e^{\\ln A + n\\theta}$ + Resolution Gauge ($\\delta_{\\mathrm{res}}, K_{\\mathrm{cut}}$)\n"
         "Multiplicative ($\\mathrm{Tri}$) $\\leftrightarrow$ Additive Spacetime & Duality ($\\cosh\\phi$ Lorentz vs. $\\cos\\theta$ Born)\n"
         "• Visible (above resolution, past light cone): Spacetime $(x, y, z, t)$, tides, particles\n"
         "• Hidden (below resolution, proj. axis $R_0$): Curvature $R$, internal $(Q_1, Q_2, Q_3)$, $SU(2)_L$"),
        "#fdf2f8", "#9b2265", fs=8.1, fw="bold")

    ax.annotate("", xy=(2.55, 7.52), xytext=(3.8, 8.18), arrowprops=dict(arrowstyle="->", lw=1.7, color="#333"))
    ax.annotate("", xy=(7.45, 7.52), xytext=(6.2, 8.18), arrowprops=dict(arrowstyle="->", lw=1.7, color="#333"))
    ax.annotate("", xy=(3.4, 3.12), xytext=(2.55, 3.78), arrowprops=dict(arrowstyle="->", lw=1.7, color="#333"))
    ax.annotate("", xy=(6.6, 3.12), xytext=(7.45, 3.78), arrowprops=dict(arrowstyle="->", lw=1.7, color="#333"))

    ax.set_title(
        r"(b) 64 次元状態空間から $(x,y,z,t,R,Q_1,Q_2,Q_3)$＋$SU(2)_L$ への対応"
        if ja else
        r"(b) Complete 64D Mapping to $(x,y,z,t,R,Q_1,Q_2,Q_3)$ + $SU(2)_L$",
        fontsize=10.6, fontweight="bold"
    )

    fig.suptitle(
        r"Part III 検証 (2)：分解能ゲージによる可視化境界と 64 次元時空・内部空間対応図"
        if ja else
        r"Part III Verification (2): Resolution-Gauge Visibility Boundary & 64D Spacetime–Internal Mapping",
        fontsize=12.0, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p3_resolution_gauge_and_mapping_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    for lang in ("ja", "en"):
        plot_fig1(lang)
        plot_fig2(lang)
