#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 7 速報（八元数 Triality スピン 1/2・CAR フェルミオン電荷 1:3:3:1・カイラル分裂）用の
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
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=200)

    # (a) Triality 線形系 (512×84) の 84 次元分解と厳密機械精度検証バー
    ax = axes[0]
    labels = [
        "最小非ゼロ特異値\n(rank=56)" if ja else "Min nonzero SV\n(rank=56)",
        "最大ゼロ特異値\n(Spin(8): 28次元)" if ja else "Max zero SV\n(Spin(8): dim 28)",
        "有限回転誤差\n($\\theta=0.7$)" if ja else "Finite rot. err\n($\\theta=0.7$)",
        "$||R_v(2\\pi)-I||$\n(ベクトル)" if ja else "$||R_v(2\\pi)-I||$\n(Vector)",
        "$||R_s(2\\pi)+I||$\n(スピノル$-I$)" if ja else "$||R_s(2\\pi)+I||$\n(Spinor $-I$)",
    ]
    vals = [4.898979, 4.489e-15, 8.934e-16, 7.622e-16, 1.584e-14]
    cols = ["#7f7f7f", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    x = np.arange(len(labels))
    ax.bar(x, vals, width=0.52, color=cols, edgecolor="#222", alpha=0.88)
    ax.set_yscale("log")
    ax.set_ylim(1e-17, 100.0)
    ax.axhline(1e-13, color="#888", ls="--", lw=1.1, label="倍精度機械ゼロ水準" if ja else "Machine-zero level")
    for xi, v in zip(x, vals):
        lbl = f"{v:.2f}" if v > 1 else f"{v:.1e}"
        ax.text(xi, v * 2.3, lbl, ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.4)
    ax.set_ylabel("特異値 / 残差ノルム（対数）" if ja else "Singular value / residual norm (log)", fontsize=10.2)
    ax.set_title(
        r"(a) $Spin(8)$ Triality（28 次元）の厳密不変性と機械精度"
        if ja else
        r"(a) Exact 28D $Spin(8)$ Triality Invariance (Machine Precision)",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.92)

    # (b) 空間回転角 θ ∈ [0, 4π] に対するベクトル (周期 2π) vs スピノル (周期 4π, 2π で -1)
    ax = axes[1]
    th = np.linspace(0, 4 * np.pi, 300)
    # (e1, e2) 回転: ベクトルは 2 次元だけ cos(θ), 残り 6 次元は +1 -> 回転平面成分 cos(θ)
    # スピノルは全 8 次元が ±θ/2 で回転 -> 固有値 e^{±iθ/2} -> 実部 cos(θ/2)
    vec_rot = np.cos(th)
    spin_rot = np.cos(0.5 * th)
    ax.plot(th / np.pi, vec_rot, "--", color="#1f77b4", lw=2.3,
            label=r"ベクトル $8_v$: 回転平面 $\cos\theta$（$2\pi$ で $+1$）" if ja else r"Vector $8_v$: plane $\cos\theta$ ($+1$ at $2\pi$)")
    ax.plot(th / np.pi, spin_rot, "-", color="#d62728", lw=2.6,
            label=r"スピノル $8_s, 8_c$: $\frac{1}{8}\mathrm{Tr}\,R_{s,c}(\theta) = \cos(\theta/2)$（$2\pi$ で $-1$）" if ja else r"Spinors $8_s, 8_c$: $\frac{1}{8}\mathrm{Tr}\,R_{s,c}(\theta) = \cos(\theta/2)$ ($-1$ at $2\pi$)")
    ax.scatter([2.0, 4.0], [1.0, 1.0], s=85, color="#1f77b4", edgecolor="#111", zorder=5)
    ax.scatter([2.0, 4.0], [-1.0, 1.0], s=100, color="#d62728", edgecolor="#111", zorder=6)
    ax.axhline(0, color="#aaa", lw=0.9)
    ax.axvline(2.0, color="#d62728", ls=":", lw=1.4, alpha=0.7)
    ax.annotate(
        "$\\theta = 2\\pi$: $R_{s,c} = -I$\n(スピン $1/2$ 符号反転！)" if ja else "$\\theta = 2\\pi$: $R_{s,c} = -I$\n(Spin-1/2 sign flip!)",
        xy=(2.0, -1.0), xytext=(2.25, -0.78),
        arrowprops=dict(arrowstyle="->", color="#d62728"), fontsize=8.8, fontweight="bold", color="#d62728",
        bbox=dict(boxstyle="round,pad=0.25", fc="#fff5f5", ec="#d62728", alpha=0.92)
    )
    ax.set_xlim(0, 4.15)
    ax.set_ylim(-1.25, 1.32)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(["0", r"$\pi$", r"$2\pi$", r"$3\pi$", r"$4\pi$"], fontsize=10)
    ax.set_xlabel(r"空間 $(e_1, e_2)$ 回転角 $\theta$" if ja else r"Spatial $(e_1, e_2)$ rotation angle $\theta$", fontsize=10.5)
    ax.set_ylabel("回転表現の位相符号" if ja else "Representation phase sign", fontsize=10.5)
    ax.set_title(
        r"(b) $2\pi$ 回転でのスピノル符号反転（$R_{s,c}(2\pi) = -I$, スピン $1/2$）"
        if ja else
        r"(b) Spin-1/2 Sign Flip Under $2\pi$ Rotation ($R_{s,c}(2\pi) = -I$)",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

    fig.suptitle(
        r"外側八元数結晶 $O_8$ の 3 添字 $Spin(8)$ Triality とスピン $1/2$ フェルミオン表現の厳密実証"
        if ja else
        r"Exact 3-Index $Spin(8)$ Triality and Spin-1/2 Fermionic Representation on Outer Octonion Crystal $O_8$",
        fontsize=12.3, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r7_triality_spin_half_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=200)

    # (a) e_7 固定 CAR フェルミオン代数による標準模型 1 世代電荷スペクトル (1:3:3:1)
    ax = axes[0]
    q_labels = [
        "$Q = 0$\n($N=0$)\n$\\nu$ 型" if ja else "$Q = 0$\n($N=0$)\n$\\nu$-type",
        "$Q = +1/3$\n($N=1$)\n$\\bar{d}$ 型" if ja else "$Q = +1/3$\n($N=1$)\n$\\bar{d}$-type",
        "$Q = +2/3$\n($N=2$)\n$u$ 型" if ja else "$Q = +2/3$\n($N=2$)\n$u$-type",
        "$Q = +1$\n($N=3$)\n$e^+$ 型" if ja else "$Q = +1$\n($N=3$)\n$e^+$-type",
    ]
    degs = [1, 3, 3, 1]
    cols = ["#1f77b4", "#d62728", "#ff7f0e", "#1f77b4"]
    x = np.arange(4)
    bars = ax.bar(x, degs, width=0.52, color=cols, edgecolor="#222", alpha=0.88)
    c2_lbls = [
        "1重項\n($C_2=0$)" if ja else "Singlet\n($C_2=0$)",
        "3重項\n($C_2=16/3$)" if ja else "Triplet\n($C_2=16/3$)",
        "反3重項\n($C_2=16/3$)" if ja else "Anti-triplet\n($C_2=16/3$)",
        "1重項\n($C_2=0$)" if ja else "Singlet\n($C_2=0$)",
    ]
    for xi, d, cl in zip(x, degs, c2_lbls):
        ax.text(xi, d + 0.12, f"{d}\n{cl}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(q_labels, fontsize=9.2)
    ax.set_ylim(0, 5.25)
    ax.set_ylabel("状態数（縮退度：1 : 3 : 3 : 1）" if ja else "State multiplicity (1 : 3 : 3 : 1)", fontsize=10.5)
    ax.set_title(
        r"(a) CAR フェルミオン代数と標準模型 1 世代電荷・カラー構造"
        if ja else
        r"(a) CAR Fermion Algebra & One-Generation Charge/Color Spectrum",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(
        0.03, 0.78,
        ("パウリ排他則 $\\{\\alpha_m,\\alpha_n\\}=0$: 誤差 $0.000\\times 10^0$\n"
         "$\\{\\alpha_m,\\alpha_n^\\dagger\\}=\\delta_{mn}I$: 誤差 $0.000\\times 10^0$\n"
         "$[\\Lambda_A, Q]=0$ ($SU(3)$ 8本): 誤差 $0.000\\times 10^0$")
        if ja else
        ("Pauli exclusion $\\{\\alpha_m,\\alpha_n\\}=0$: err $0.000\\times 10^0$\n"
         "$\\{\\alpha_m,\\alpha_n^\\dagger\\}=\\delta_{mn}I$: err $0.000\\times 10^0$\n"
         "$[\\Lambda_A, Q]=0$ (8 $SU(3)$ gens): err $0.000\\times 10^0$"),
        transform=ax.transAxes, fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faff", ec="#1f77b4", alpha=0.94)
    )

    # (b) 時間結晶背景 H = O_8 ⊗ M_8 上のカイラル・パリティ・CP 応答 (||P|| 分裂)
    ax = axes[1]
    conds = [
        "$O_8 \\otimes M_8$\n(同カイラル R,R)" if ja else "$O_8 \\otimes M_8$\n(Co-chiral R,R)",
        "$O_8 \\otimes M_8^T$\n(逆カイラル R,L)" if ja else "$O_8 \\otimes M_8^T$\n(Anti-chiral R,L)",
        "$O_8^T \\otimes M_8$\n(空間反転 L,R)" if ja else "$O_8^T \\otimes M_8$\n(Parity L,R)",
        "$O_8^T \\otimes M_8^T$\n(CP 反転 L,L)" if ja else "$O_8^T \\otimes M_8^T$\n(CP-flip L,L)",
    ]
    nP_vals = [0.046873, 0.046212, 0.046873, 0.046212]
    bcols = ["#1f77b4", "#d62728", "#1f77b4", "#d62728"]
    xx = np.arange(4)
    ax.bar(xx, nP_vals, width=0.50, color=bcols, edgecolor="#222", alpha=0.88)
    for xi, v in zip(xx, nP_vals):
        ax.text(xi, v + 0.00008, f"{v:.6f}", ha="center", va="bottom", fontsize=8.8, fontweight="bold")
    ax.set_xticks(xx)
    ax.set_xticklabels(conds, fontsize=8.8)
    ax.set_ylim(0.0455, 0.0475)
    ax.set_ylabel(r"$\tau = 30$ 自己無撞着ノルム $\|P\|$" if ja else r"Self-consistent norm $\|P\|$ at $\tau = 30$", fontsize=10.5)
    ax.set_title(
        r"(b) 時間結晶内側 $M_8$ によるカイラル非対称分裂（$\Delta\|P\| = -1.41\%$）"
        if ja else
        r"(b) Chiral Splitting Driven by Inner Time-Crystal $M_8$ ($\Delta\|P\| = -1.41\%$)",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.annotate(
        "内側回転向き ($M_8$ vs $M_8^T$) のみが\n左右カイラリティを差別！" if ja else "Only inner rotation ($M_8$ vs $M_8^T$)\nsplits left/right chirality!",
        xy=(1, 0.046212), xytext=(1.55, 0.04718),
        arrowprops=dict(arrowstyle="->", color="#d62728"), fontsize=8.4,
        bbox=dict(boxstyle="round,pad=0.25", fc="#fff8f0", ec="#d62728", alpha=0.92)
    )

    fig.suptitle(
        r"電荷の代数的起源（1:3:3:1＝空間倍音と同一の $\mathbb{Z}_2^3$ 構造）と時間結晶カイラル分裂"
        if ja else
        r"Algebraic Origin of Electric Charge (1:3:3:1 Matching $\mathbb{Z}_2^3$ Harmonics) and Time-Crystal Chiral Splitting",
        fontsize=12.2, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r7_fermion_charge_and_chiral_{lang}.png"
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
