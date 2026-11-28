#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part IV（電弱対称性の動的破れ・パリティ最大破れ・W/Z 質量獲得・sin^2 θ_W = 3/8）用の
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

    # (a) G_2 内の 3+4 分割における自己双対(SD) vs 反自己双対(ASD) 100%:0% 厳密分離（パリティ最大破れ）
    ax = axes[0]
    groups = [
        "$SU(2)_{\\mathrm{pure4}}$\n(3空間不変・内部弱)\n$G_2$ 微分残差 0"
        if ja else
        "$SU(2)_{\\mathrm{pure4}}$\n(3-space inv., weak)\n$G_2$ res. 0",
        "$SU(2)_{\\mathrm{space}}$\n(3空間回転ロック)\n$SO(3)$ 随伴"
        if ja else
        "$SU(2)_{\\mathrm{space}}$\n(Locked to 3-space)\n$SO(3)$ adjoint",
    ]
    sd_vals = [1.0, 0.0]
    asd_vals = [0.0, 1.0]
    x = np.arange(2)
    w = 0.34
    ax.bar(x - w / 2, sd_vals, width=w, color="#1f77b4", edgecolor="#222",
           label=r"自己双対 SD（+カイラル / 左巻き $SU(2)_L$）" if ja else r"Self-dual SD (+chiral / left-handed $SU(2)_L$)")
    ax.bar(x + w / 2, asd_vals, width=w, color="#d62728", edgecolor="#222",
           label=r"反自己双対 ASD（−カイラル / 空間ロック）" if ja else r"Anti-self-dual ASD (-chiral / space-locked)")
    for xi, sv, av in zip(x, sd_vals, asd_vals):
        ax.text(xi - w / 2, sv + 0.025, f"{sv*100:.0f}%", ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#1f77b4")
        ax.text(xi + w / 2, av + 0.025, f"{av*100:.0f}%", ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=9.0)
    ax.set_ylim(0, 1.24)
    ax.set_ylabel(r"4 空間 $\{e_4,e_5,e_6,e_7\}$ のカイラル射影占有率" if ja else r"Chiral projection fraction on $\{e_4,e_5,e_6,e_7\}$", fontsize=10.0)
    ax.set_title(
        r"(a) $G_2 \supset SO(4)$ の 100%:0% カイラル分離（パリティ最大破れ）"
        if ja else
        r"(a) Exact 100%:0% Chiral Split in $G_2 \supset SO(4)$ (Maximal Parity Violation)",
        fontsize=10.4, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper center", fontsize=8.3, framealpha=0.94)
    ax.text(
        0.04, 0.50,
        ("左巻きダブレット ($e_4$): $W_1$ 応答 $= 0.3236$\n"
         "右巻きシングレット ($e_0,e_1$): $0.000\\times 10^0$\n"
         "$[SO(3),\\,\\mathrm{SM}_{12}]=0$ (全36組誤差 $0.000$)")
        if ja else
        ("Left doublet ($e_4$): $W_1$ resp. $= 0.3236$\n"
         "Right singlets ($e_0,e_1$): $0.000\\times 10^0$\n"
         "$[SO(3),\\,\\mathrm{SM}_{12}]=0$ (36 pairs err $0.000$)"),
        transform=ax.transAxes, fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) Tri(T) 実動力学によるヒッグス凝縮 P4_4 の W ボソン質量生成 vs 残余 8D 無質量
    ax = axes[1]
    eps_h = np.array([0.00, 0.05, 0.10, 0.16, 0.22])
    w_mass = np.array([0.000000, 0.012593, 0.025105, 0.035323, 0.035517])
    rem_8d = np.array([3.802e-14, 3.821e-14, 3.821e-14, 3.828e-14, 3.821e-14])

    ax.plot(eps_h, w_mass, "o-", color="#d62728", lw=2.5, ms=8.0,
            label=r"$W_1, W_2, W_3$ ($SU(2)_L$ 3成分 6桁完全縮退)" if ja else r"$W_1, W_2, W_3$ ($SU(2)_L$ 3-comp. exact degeneracy)")
    ax.plot(eps_h, rem_8d, "s--", color="#2ca02c", lw=2.2, ms=7.0,
            label=r"残余 8D 対称群 ($e_4$ 固定, $\leq 3.8\times 10^{-14}$ 無質量)" if ja else r"Residual 8D symmetry ($e_4$ fixed, $\leq 3.8\times 10^{-14}$ massless)")
    for eh, wm in zip(eps_h, w_mass):
        if eh > 0:
            ax.text(eh, wm + 0.0014, f"{wm:.4f}", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color="#b22222")
    ax.set_xlim(-0.01, 0.24)
    ax.set_ylim(-0.003, 0.043)
    ax.set_xlabel(r"ヒッグス凝縮体 P4$_4$（site $A=4$）振幅 $\epsilon_H$ ($\tau=40$)" if ja else r"Higgs condensate P4$_4$ (site $A=4$) amplitude $\epsilon_H$ ($\tau=40$)", fontsize=10.0)
    ax.set_ylabel(r"対称性破れ秩序パラメータ $\|G(U_{\mathrm{higgs}})\|_F$" if ja else r"Symmetry-breaking order parameter $\|G(U_{\mathrm{higgs}})\|_F$", fontsize=10.0)
    ax.set_title(
        r"(b) $\mathrm{Tri}(T)$ 実動力学：$W_1,W_2,W_3$ 質量生成 vs 残余無質量"
        if ja else
        r"(b) $\mathrm{Tri}(T)$ Dynamics: $W_1,W_2,W_3$ Mass Generation vs. Massless Residual",
        fontsize=10.4, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center right", bbox_to_anchor=(0.99, 0.36), fontsize=8.1, framealpha=0.94)

    fig.suptitle(
        r"Part IV 検証 (1)：$G_2$ によるパリティ最大破れの代数的証明と $\mathrm{Tri}(T)$ 動的ヒッグス質量生成"
        if ja else
        r"Part IV Verification (1): Algebraic Maximal Parity Violation in $G_2$ and Dynamical Higgs Mass in $\mathrm{Tri}(T)$",
        fontsize=11.8, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p4_parity_and_higgs_mass_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), dpi=200)

    # (a) 4x4 電弱質量行列の固有値 [m_γ^2, m_W1^2, m_W2^2, m_Z^2]
    ax = axes[0]
    bosons = [
        "$\\gamma$ (光子)\n厳密無質量" if ja else "$\\gamma$ (Photon)\nExact massless",
        "$W^+$ ボソン\n(縮退)" if ja else "$W^+$ Boson\n(Degenerate)",
        "$W^-$ ボソン\n(縮退)" if ja else "$W^-$ Boson\n(Degenerate)",
        "$Z^0$ ボソン\n($m_W^2/m_Z^2=5/8$)" if ja else "$Z^0$ Boson\n($m_W^2/m_Z^2=5/8$)",
    ]
    evals = [0.0, 0.25, 0.25, 0.40]
    colors = ["#2ca02c", "#1f77b4", "#1f77b4", "#d62728"]
    x = np.arange(4)
    ax.bar(x, evals, width=0.52, color=colors, edgecolor="#222")
    ax.text(0, 0.015, "$-2.1\\times 10^{-17}$\n(= 0)", ha="center", va="bottom", fontsize=8.4, fontweight="bold", color="#2ca02c")
    ax.text(1, 0.262, "0.2500", ha="center", va="bottom", fontsize=9.0, fontweight="bold", color="#1f77b4")
    ax.text(2, 0.262, "0.2500", ha="center", va="bottom", fontsize=9.0, fontweight="bold", color="#1f77b4")
    ax.text(3, 0.412, "0.4000", ha="center", va="bottom", fontsize=9.0, fontweight="bold", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(bosons, fontsize=8.8)
    ax.set_ylim(0, 0.50)
    ax.set_ylabel(r"電弱質量行列 $M^2$ の固有値（無次元規格化）" if ja else r"Electroweak mass-squared matrix $M^2$ eigenvalues", fontsize=10.0)
    ax.set_title(
        r"(a) 4×4 電弱質量スペクトル（$g'/g = \sqrt{3/5}$、調整パラメータ 0）"
        if ja else
        r"(a) 4×4 Electroweak Mass Spectrum ($g'/g = \sqrt{3/5}$, Zero Free Params)",
        fontsize=10.3, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(
        0.04, 0.72,
        ("複素ヒッグス二重項 ($T_3=-1/2, Y=+1/2 \\Rightarrow Q=0$):\n"
         "・$m_\\gamma^2 = 0$（機械ゼロ）\n"
         "・$m_{W^+}^2 = m_{W^-}^2 = 0.2500$\n"
         "・$m_W^2 / m_Z^2 = 0.6250 = 5/8 = \\cos^2\\theta_W^{\\mathrm{bare}}$")
        if ja else
        ("Complex Higgs doublet ($T_3=-1/2, Y=+1/2 \\Rightarrow Q=0$):\n"
         "• $m_\\gamma^2 = 0$ (machine zero)\n"
         "• $m_{W^+}^2 = m_{W^-}^2 = 0.2500$\n"
         "• $m_W^2 / m_Z^2 = 0.6250 = 5/8 = \\cos^2\\theta_W^{\\mathrm{bare}}$"),
        transform=ax.transAxes, fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff8f0", ec="#d62728", alpha=0.94)
    )

    # (b) 64D トレース比による bare ワインバーグ角 sin^2(θ_W) = 3/8 = 0.375
    ax = axes[1]
    tr_labels = [
        "$\\mathrm{Tr}(T_3^2)$\n(弱 $SU(2)_L$)" if ja else "$\\mathrm{Tr}(T_3^2)$\n(Weak $SU(2)_L$)",
        "$\\mathrm{Tr}(Y^2)$\n(ハイパーチャージ)" if ja else "$\\mathrm{Tr}(Y^2)$\n(Hypercharge)",
        "$\\mathrm{Tr}(Q_{\\mathrm{em}}^2)$\n(電磁電荷)" if ja else "$\\mathrm{Tr}(Q_{\\mathrm{em}}^2)$\n(EM Charge)",
    ]
    tr_vals = [8.0, 40.0 / 3.0, 64.0 / 3.0]
    tr_cols = ["#1f77b4", "#ff7f0e", "#9467bd"]
    x2 = np.arange(3)
    ax.bar(x2, tr_vals, width=0.50, color=tr_cols, edgecolor="#222")
    ax.text(0, 8.5, "8.000\n(= 24/3)", ha="center", va="bottom", fontsize=8.8, fontweight="bold", color="#1f77b4")
    ax.text(1, 13.8, "13.333\n(= 40/3)", ha="center", va="bottom", fontsize=8.8, fontweight="bold", color="#d95f02")
    ax.text(2, 21.8, "21.333\n(= 64/3)", ha="center", va="bottom", fontsize=8.8, fontweight="bold", color="#6a3d9a")
    ax.set_xticks(x2)
    ax.set_xticklabels(tr_labels, fontsize=9.0)
    ax.set_ylim(0, 26.5)
    ax.set_ylabel(r"64 次元空間 $O_8^{\mathrm{ext}}\otimes O_8^{\mathrm{int}}$ 上のトレース二乗和" if ja else r"Trace-squared over 64D $O_8^{\mathrm{ext}}\otimes O_8^{\mathrm{int}}$", fontsize=10.0)
    ax.set_title(
        r"(b) 素のワインバーグ角 $\sin^2\theta_W^{\mathrm{bare}} = \mathrm{Tr}(T_3^2)/\mathrm{Tr}(Q^2) = 3/8$"
        if ja else
        r"(b) Bare Weinberg Angle $\sin^2\theta_W^{\mathrm{bare}} = \mathrm{Tr}(T_3^2)/\mathrm{Tr}(Q^2) = 3/8$",
        fontsize=10.3, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(
        0.04, 0.72,
        ("$SU(5)/SO(10)$ を仮定せず 64D トレース比だけで:\n"
         "$\\sin^2\\theta_W^{\\mathrm{bare}} = \\frac{8}{64/3} = \\frac{3}{8} = 0.375000$\n"
         "（大統一理論の理論値 $3/8$ と 6 桁厳密一致）\n"
         "※ 低エネルギー値 $\\approx 0.231$ への RG 走行は残課題")
        if ja else
        ("Without assuming $SU(5)/SO(10)$, 64D trace ratio yields:\n"
         "$\\sin^2\\theta_W^{\\mathrm{bare}} = \\frac{8}{64/3} = \\frac{3}{8} = 0.375000$\n"
         "(Exact 6-digit match to canonical GUT value $3/8$)\n"
         "* RG running to low-energy $\\approx 0.231$ is open"),
        transform=ax.transAxes, fontsize=7.9,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    fig.suptitle(
        r"Part IV 検証 (2)：4×4 電弱質量スペクトル（$\gamma, W^\pm, Z^0$）と素のワインバーグ角 $\sin^2\theta_W^{\mathrm{bare}} = 3/8$"
        if ja else
        r"Part IV Verification (2): 4×4 Electroweak Mass Spectrum ($\gamma, W^\pm, Z^0$) and Bare Weinberg Angle $\sin^2\theta_W^{\mathrm{bare}} = 3/8$",
        fontsize=11.8, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_p4_ew_spectrum_and_weinberg_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    for lang in ("ja", "en"):
        plot_fig1(lang)
        plot_fig2(lang)
