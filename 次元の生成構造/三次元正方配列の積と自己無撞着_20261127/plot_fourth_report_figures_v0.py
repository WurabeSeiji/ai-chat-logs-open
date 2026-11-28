#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 4 速報（量子構造の創発）用の図 2 種（日本語版・英語版、計 4 枚）を生成する専用スクリプト。
実測データソース:
  - results/tri_loop_n64_quantum_foundations_verify_v0.txt / .npz
出力:
  - results/fig_r4_complex_and_born_ja.png / _en.png
  - results/fig_r4_uncertainty_and_cutoff_ja.png / _en.png
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
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), dpi=200)

    # (a) 実配列 M_8 のパウリ成分 (対数スケール): 複素構造 I, J の自発選択
    ax = axes[0]
    labels = [
        r"$|c_0|$ ($I$)" + ("\n実部" if ja else "\nReal"),
        r"$|c_J|$ ($i\sigma_y$)" + ("\n虚部 $J$" if ja else "\nImag $J$"),
        r"$|c_x|$ ($\sigma_x$)" + ("\n異方性" if ja else "\nAniso."),
        r"$|c_z|$ ($\sigma_z$)" + ("\n異方性" if ja else "\nAniso."),
    ]
    b0_vals = [0.27343360, 0.49200681, 1.943e-16, 6.661e-16]
    b1_vals = [0.20405223, 0.37619949, 5.551e-17, 5.274e-16]
    x = np.arange(len(labels))
    w = 0.34
    ax.bar(x - w / 2, b0_vals, width=w, color="#1f77b4", edgecolor="#222", label=r"$b=0$ ブロック" if ja else r"$b=0$ block")
    ax.bar(x + w / 2, b1_vals, width=w, color="#ff7f0e", edgecolor="#222", label=r"$b=1$ ブロック" if ja else r"$b=1$ block")
    ax.set_yscale("log")
    ax.set_ylim(1e-18, 5.0)
    ax.axhline(1e-15, color="#888", ls="--", lw=1.1, label="倍精度機械ゼロ水準 ($10^{-15}$)" if ja else "Machine-zero level ($10^{-15}$)")
    for xi, v0, v1 in zip(x, b0_vals, b1_vals):
        if v0 > 1e-3:
            ax.text(xi - w / 2, v0 * 1.4, f"{v0:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#1f77b4")
            ax.text(xi + w / 2, v1 * 1.4, f"{v1:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#d95f02")
        else:
            ax.text(xi, max(v0, v1) * 3.5, r"$\sim 10^{-16}$", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#555")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel("パウリ基底成分の大きさ（対数）" if ja else "Pauli-basis component magnitude (log)", fontsize=10.5)
    ax.set_title(
        r"(a) 実配列 $M_8$ からの純粋複素構造 $\mathbb{C}$ の自発創発"
        if ja else
        r"(a) Spontaneous Complex Structure $\mathbb{C}$ from Real $M_8$",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.92)
    ax.text(
        0.03, 0.06,
        ("第3層では $\\mathbb{C}_{\\mathrm{out}}\\otimes_{\\mathbb{R}}\\mathbb{C}_{\\mathrm{in}} \\to \\mathbb{C}_{\\mathrm{single}}$ にロック\n"
         "$\\sigma_1=\\sigma_2=1/\\sqrt{2}$ (誤差 $3.2\\times 10^{-7}$), $\\|P\\|\\times 2$")
        if ja else
        ("Layer 3 locks $\\mathbb{C}_{\\mathrm{out}}\\otimes_{\\mathbb{R}}\\mathbb{C}_{\\mathrm{in}} \\to \\mathbb{C}_{\\mathrm{single}}$\n"
         "$\\sigma_1=\\sigma_2=1/\\sqrt{2}$ (err $3.2\\times 10^{-7}$), $\\|P\\|\\times 2$"),
        transform=ax.transAxes, fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f7fbff", ec="#1f77b4", alpha=0.93)
    )

    # (b) O(64) 回転共変性と L^p (p=1,2,4) の回転変化率比較 (Born 則の一意性)
    ax = axes[1]
    cats = [
        r"$L^1$ ノルム" if ja else r"$L^1$ norm",
        r"$L^2$ ノルム" + ("\n(ヒルベルト / Born)" if ja else "\n(Hilbert / Born)"),
        r"$L^4$ ノルム" if ja else r"$L^4$ norm",
        r"$\mathrm{Tri}$ 共変誤差" if ja else r"$\mathrm{Tri}$ cov. error",
    ]
    vals = [3.518, 2.220e-16, 0.4202, 8.843e-16]
    cols = ["#d62728", "#2ca02c", "#d62728", "#1f77b4"]
    xx = np.arange(len(cats))
    bars = ax.bar(xx, vals, width=0.52, color=cols, edgecolor="#222", alpha=0.88)
    ax.set_yscale("log")
    ax.set_ylim(1e-18, 90.0)
    ax.axhline(1e-15, color="#888", ls="--", lw=1.1)
    for xi, v in zip(xx, vals):
        lbl = f"{v*100:.1f}%" if v > 1e-3 else f"{v:.2e}"
        ax.text(xi, v * 2.3, lbl, ha="center", va="bottom", fontsize=9.0, fontweight="bold")
    ax.set_xticks(xx)
    ax.set_xticklabels(cats, fontsize=9.5)
    ax.set_ylabel(r"$O(64)$ 回転に対する相対変化率（対数）" if ja else r"Relative change under $O(64)$ rotation (log)", fontsize=10.5)
    ax.set_title(
        r"(b) $O(64)$ 共変性と $L^2$ ノルムの一意性（Born 則 $P=|\psi|^2$）"
        if ja else
        r"(b) $O(64)$ Covariance & Uniqueness of $L^2$ Norm (Born Rule)",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(
        0.26, 0.36,
        ("直交回転で不変なのは $p=2$ のみ\n"
         "＋複素構造 $\\mathbb{C}$・次元 $\\geq 3$\n"
         "$\\Rightarrow$ Gleason の定理で $P=|\\psi|^2$ 一意決定")
        if ja else
        ("Only $p=2$ is $O(N)$-invariant\n"
         "+ complex structure $\\mathbb{C}$, dim $\\geq 3$\n"
         "$\\Rightarrow$ Gleason's theorem fixes $P=|\\psi|^2$"),
        transform=ax.transAxes, fontsize=8.3,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4fbf4", ec="#2ca02c", alpha=0.95)
    )

    fig.suptitle(
        r"定理 Q2・Q3：実 3 階テンソルからの単一複素構造 $\mathbb{C}$ の自発創発と Born 則（$L^2$ 一意性）"
        if ja else
        r"Theorems Q2 & Q3: Spontaneous Single Complex Structure $\mathbb{C}$ and Uniqueness of the $L^2$ Born Rule",
        fontsize=12.5, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r4_complex_and_born_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.5), dpi=200)
    ln8 = np.log(8.0)

    # (a) (H_space, H_harmonic) 平面上のエントロピー不確定性境界と裸／ドレスド状態
    ax = axes[0]
    hx_line = np.linspace(0, ln8, 100)
    hk_line = ln8 - hx_line
    ax.fill_between(hx_line, 0, hk_line, color="#ffdddd", alpha=0.55,
                    label=r"禁止領域 $H_x + H_k < \ln 8$" if ja else r"Forbidden $H_x + H_k < \ln 8$")
    ax.plot(hx_line, hk_line, "-", color="#d62728", lw=2.3,
            label=r"厳密不確定性下限 $H_x + H_k = \ln 8 \approx 2.079$" if ja else r"Exact lower bound $H_x + H_k = \ln 8 \approx 2.079$")

    # 裸の芯 3 点 (厳密飽和)
    bare_x = [ln8, np.log(4.0), 0.0]
    bare_k = [0.0, np.log(2.0), ln8]
    ax.scatter(bare_x, bare_k, s=110, color="#1f77b4", edgecolor="#111", zorder=5, marker="o",
               label="裸の芯 (非局在 / $H_4$ / $\\delta_{000}$): 下限を厳密飽和" if ja else "Bare cores (Nonlocal / $H_4$ / $\\delta_{000}$): exact saturation")

    # ドレスド安定状態 (τ=30)
    dr_x = [2.079442, 2.058650, 2.022809, 1.793664]
    dr_k = [0.000000, 0.072019, 0.129825, 0.598750]
    ax.scatter(dr_x[1:], dr_k[1:], s=110, color="#ff7f0e", edgecolor="#111", zorder=6, marker="s",
               label=r"$\tau=30$ ドレスド状態 ($L_x$, P5, P4): 衣の分だけ上方シフト" if ja else r"$\tau=30$ dressed states ($L_x$, P5, P4): shifted above bound")

    # 注釈
    ax.annotate("単一サイト芯 $\\delta_{000}$\n(最小分解能 $\\Delta x=1$)" if ja else "1-site core $\\delta_{000}$\n(min resolution $\\Delta x=1$)",
                xy=(0.0, ln8), xytext=(0.25, 1.96), arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=8.5)
    ax.annotate("四元数芯 $H_4$\n(4サイト × 2モード)" if ja else "Quaternion core $H_4$\n(4 sites × 2 modes)",
                xy=(np.log(4), np.log(2)), xytext=(1.02, 1.24), arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=8.5)
    ax.annotate("非局在 $H$, P3a\n(8サイト × $k=0$のみ)" if ja else "Nonlocal $H$, P3a\n(8 sites × $k=0$ only)",
                xy=(ln8, 0.0), xytext=(1.02, 0.06), arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=8.5)
    ax.annotate("P4 (1サイト+衣)\n+0.313" if ja else "P4 (1-site+cloud)\n+0.313",
                xy=(1.7937, 0.5988), xytext=(1.55, 1.00), arrowprops=dict(arrowstyle="->", color="#d95f02"), fontsize=8.5, fontweight="bold")
    ax.annotate("P5 ($H_4$+衣), $L_x$" if ja else "P5 ($H_4$+cloud), $L_x$",
                xy=(2.0228, 0.1298), xytext=(1.55, 0.36), arrowprops=dict(arrowstyle="->", color="#d95f02"), fontsize=8.5)

    ax.set_xlim(-0.12, 2.30)
    ax.set_ylim(-0.12, 2.30)
    ax.set_xlabel(r"空間エントロピー $H_{\mathrm{space}}(f)$" if ja else r"Spatial entropy $H_{\mathrm{space}}(f)$", fontsize=10.5)
    ax.set_ylabel(r"空間倍音エントロピー $H_{\mathrm{harmonic}}(\hat{f})$" if ja else r"Harmonic entropy $H_{\mathrm{harmonic}}(\hat{f})$", fontsize=10.5)
    ax.set_title(
        r"(a) $\mathbb{Z}_2^3$ エントロピー不確定性平面（裸の飽和 vs ドレスド）"
        if ja else
        r"(a) $\mathbb{Z}_2^3$ Entropic Uncertainty Plane (Bare Saturation vs Dressed)",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

    # (b) Donoho-Stark 積 |supp(x)| * |supp(k)| = 8 と最小分解能 (プランクスケール原型)
    ax = axes[1]
    sx = np.linspace(0.85, 8.5, 200)
    ax.fill_between(sx, 0, 8.0 / sx, color="#ffdddd", alpha=0.50,
                    label=r"禁止領域 $|\mathrm{supp}(x)|\cdot|\mathrm{supp}(k)| < 8$" if ja else r"Forbidden $|\mathrm{supp}(x)|\cdot|\mathrm{supp}(k)| < 8$")
    ax.plot(sx, 8.0 / sx, "-", color="#d62728", lw=2.3,
            label=r"Donoho–Stark 飽和曲線 $xy = 8$" if ja else r"Donoho–Stark saturation $xy = 8$")
    pts_x = [1, 4, 8]
    pts_k = [8, 2, 1]
    ax.scatter(pts_x, pts_k, s=120, color="#1f77b4", edgecolor="#111", zorder=5)
    ann_texts = [
        ("1サイト局在 ($\\delta_{000}$)\n$1\\times 8 = 8$\n(UV上限 $k_{\\max}=3$)" if ja else "1-site ($\\delta_{000}$)\n$1\\times 8 = 8$\n(UV cutoff $k_{\\max}=3$)", (1, 8), (2.2, 7.2)),
        ("4サイト局在 ($H_4$)\n$4\\times 2 = 8$" if ja else "4-site ($H_4$)\n$4\\times 2 = 8$", (4, 2), (4.6, 3.3)),
        ("非局在波 ($H$, P3a)\n$8\\times 1 = 8$" if ja else "Nonlocal ($H$, P3a)\n$8\\times 1 = 8$", (8, 1), (5.5, 1.55)),
    ]
    for txt, xy, xytext in ann_texts:
        ax.annotate(txt, xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=8.8)

    ax.axvline(1.0, color="#2ca02c", ls="--", lw=1.6, label=r"最小空間分解能 $\Delta x_{\min}=1$ ($\ell_P$ 原型)" if ja else r"Min spatial resolution $\Delta x_{\min}=1$ ($\ell_P$ prototype)")
    ax.axhline(8.0, color="#9467bd", ls=":", lw=1.6, label=r"最大波数モード数 $=8$ ($k_{\max}=3$)" if ja else r"Max harmonic modes $=8$ ($k_{\max}=3$)")
    ax.set_xlim(0.3, 9.0)
    ax.set_ylim(0.0, 9.5)
    ax.set_xlabel(r"空間サポート $|\mathrm{supp}(x)|$（サイト数）" if ja else r"Spatial support $|\mathrm{supp}(x)|$ (sites)", fontsize=10.5)
    ax.set_ylabel(r"空間倍音サポート $|\mathrm{supp}(k)|$（モード数）" if ja else r"Harmonic support $|\mathrm{supp}(k)|$ (modes)", fontsize=10.5)
    ax.set_title(
        r"(b) Donoho–Stark 飽和と最小分解能（$\ell_P, t_P, \hbar$ 原型）"
        if ja else
        r"(b) Donoho–Stark Saturation & Minimum Resolution ($\ell_P, t_P, \hbar$)",
        fontsize=11.0, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

    fig.suptitle(
        r"定理 Q4：$\mathbb{Z}_2^3$ 離散不確定性原理の厳密飽和と最小分解能（プランクスケール原型）"
        if ja else
        r"Theorem Q4: Exact Saturation of $\mathbb{Z}_2^3$ Discrete Uncertainty and Minimum Resolution",
        fontsize=12.5, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r4_uncertainty_and_cutoff_{lang}.png"
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
