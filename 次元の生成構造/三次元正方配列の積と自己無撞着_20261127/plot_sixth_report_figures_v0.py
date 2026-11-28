#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 6 速報（Hurwitz 4 段階質量階層と直積準同型指数則）用の図 2 種（日英計 4 枚）生成スクリプト。
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

    # (a) 衣/芯エネルギー比 E_cloud / E_core vs 厳密メルセンヌ理論 (8-d)/d = 2^m - 1
    ax = axes[0]
    m_vals = np.array([0, 1, 2, 3])
    theo = np.array([0.0, 1.0, 3.0, 7.0])
    e05 = np.array([0.00000, 1.00000, 2.98804, 6.68202])
    e10 = np.array([0.00000, 1.00000, 2.95264, 5.94756])
    m_cont = np.linspace(0, 3.1, 100)
    ax.plot(m_cont, 2.0 ** m_cont - 1.0, "--", color="#2ca02c", lw=2.2,
            label=r"厳密理論曲線 $(8-d)/d = 2^m - 1$ (0, 1, 3, 7)" if ja else r"Exact law $(8-d)/d = 2^m - 1$ (0, 1, 3, 7)")
    w = 0.26
    ax.bar(m_vals - w / 2, e05, width=w, color="#1f77b4", edgecolor="#222", label=r"実測 $\epsilon=0.05$ ($\tau=30$)" if ja else r"Measured $\epsilon=0.05$ ($\tau=30$)")
    ax.bar(m_vals + w / 2, e10, width=w, color="#ff7f0e", edgecolor="#222", label=r"実測 $\epsilon=0.10$ ($\tau=30$)" if ja else r"Measured $\epsilon=0.10$ ($\tau=30$)")
    for mi, v05, th in zip(m_vals, e05, theo):
        ax.text(mi - w / 2, v05 + 0.16, f"{v05:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#1f77b4")
        ax.scatter([mi], [th], s=70, color="#2ca02c", edgecolor="#111", zorder=6)
    ax.set_xticks(m_vals)
    ax.set_xticklabels([
        "$\\mathbb{O}_8$ ($d=8$)\n$m=0$",
        "$\\mathbb{H}_4$ ($d=4$)\n$m=1$",
        "$\\mathbb{C}_2$ ($d=2$)\n$m=2$",
        "$\\mathbb{R}_1$ ($d=1$)\n$m=3$",
    ], fontsize=9.5)
    ax.set_ylim(0, 8.2)
    ax.set_ylabel(r"衣／芯エネルギー比 $E_{\mathrm{cloud}} / E_{\mathrm{core}}$" if ja else r"Cloud-to-core ratio $E_{\mathrm{cloud}} / E_{\mathrm{core}}$", fontsize=10.5)
    ax.set_title(
        r"(a) Hurwitz 4 段階の衣／芯比：メルセンヌ則 $2^m - 1$"
        if ja else
        r"(a) Hurwitz 4-Stage Cloud/Core Ratio: Mersenne Law $2^m - 1$",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.92)

    # (b) 空間倍音スペクトル k=0,1,2,3 (パスカルの三角形・二項則)
    ax = axes[1]
    k0 = np.array([100.0, 90.0, 74.93, 55.67])
    k1 = np.array([0.0,   10.0, 16.72, 19.00])
    k2 = np.array([0.0,    0.0,  8.36, 19.00])
    k3 = np.array([0.0,    0.0,  0.00,  6.33])
    xx = np.arange(4)
    ax.bar(xx, k1, width=0.52, color="#1f77b4", edgecolor="#222", label=r"$k=1$ (1次倍音)" if ja else r"$k=1$ (1st harmonic)")
    ax.bar(xx, k2, width=0.52, bottom=k1, color="#ff7f0e", edgecolor="#222", label=r"$k=2$ (2次倍音)" if ja else r"$k=2$ (2nd harmonic)")
    ax.bar(xx, k3, width=0.52, bottom=k1 + k2, color="#d62728", edgecolor="#222", label=r"$k=3$ (3次倍音)" if ja else r"$k=3$ (3rd harmonic)")
    ann = ["0.00%\n(DCのみ)" if ja else "0.00%\n(DC only)", "10.00%\n[1]", "25.07%\n[2 : 1]", "44.33%\n[3 : 3 : 1]"]
    tot_h = k1 + k2 + k3
    for xi, th, a in zip(xx, tot_h, ann):
        ax.text(xi, th + 1.3, a, ha="center", va="bottom", fontsize=8.8, fontweight="bold")
    ax.set_xticks(xx)
    ax.set_xticklabels([
        "$\\mathbb{O}_8$ ($m=0$)",
        "$\\mathbb{H}_4$ ($m=1$)",
        "$\\mathbb{C}_2$ ($m=2$)",
        "$\\mathbb{R}_1$ ($m=3$)",
    ], fontsize=9.5)
    ax.set_ylim(0, 56)
    ax.set_ylabel("空間倍音エネルギー占有率 (%)" if ja else "Spatial harmonic energy fraction (%)", fontsize=10.5)
    ax.set_title(
        r"(b) 空間倍音のパスカル二項則 $\binom{m}{k}$ ($1,\; 2:1,\; 3:3:1$)"
        if ja else
        r"(b) Pascal Binomial Harmonic Law $\binom{m}{k}$ ($1,\; 2:1,\; 3:3:1$)",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.92)

    fig.suptitle(
        r"Hurwitz 4 段階局在系列（$\mathbb{O}_8 \to \mathbb{H}_4 \to \mathbb{C}_2 \to \mathbb{R}_1$）の厳密衣比・倍音法則"
        if ja else
        r"Exact Dressing and Spatial Harmonic Laws Across the Hurwitz 4-Stage Series ($\mathbb{O}_8 \to \mathbb{H}_4 \to \mathbb{C}_2 \to \mathbb{R}_1$)",
        fontsize=12.3, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r6_hurwitz_dressing_and_harmonics_{lang}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=200)

    # (a) 総エネルギー d/8 則 & コアエネルギー (d/8)^2 則 (1 : 1/4 : 1/16 : 1/64)
    ax = axes[0]
    d_arr = np.array([8.0, 4.0, 2.0, 1.0])
    e_tot_rel = np.array([0.002494, 0.001251, 0.000623, 0.000294]) / 0.002494
    e_cor_rel = np.array([0.002494, 0.000625, 0.000156, 0.000038]) / 0.002494
    d_cont = np.linspace(0.85, 9.0, 100)
    ax.plot(d_cont, d_cont / 8.0, "--", color="#1f77b4", lw=2.0,
            label=r"理論線 $d/8$ (総エネルギー $E_{\mathrm{tot}}$)" if ja else r"Exact $d/8$ (total energy $E_{\mathrm{tot}}$)")
    ax.plot(d_cont, (d_cont / 8.0) ** 2, "--", color="#d62728", lw=2.0,
            label=r"理論線 $(d/8)^2$ (コア $E_{\mathrm{core}}$: $1, 1/4, 1/16, 1/64$)" if ja else r"Exact $(d/8)^2$ (core $E_{\mathrm{core}}$: $1, 1/4, 1/16, 1/64$)")
    ax.scatter(d_arr, e_tot_rel, s=100, color="#1f77b4", edgecolor="#111", zorder=5, marker="o",
               label=r"実測 $E_{\mathrm{tot}}(d)/E_{\mathrm{tot}}(8)$" if ja else r"Measured $E_{\mathrm{tot}}(d)/E_{\mathrm{tot}}(8)$")
    ax.scatter(d_arr, e_cor_rel, s=100, color="#d62728", edgecolor="#111", zorder=6, marker="s",
               label=r"実測 $E_{\mathrm{core}}(d)/E_{\mathrm{core}}(8)$" if ja else r"Measured $E_{\mathrm{core}}(d)/E_{\mathrm{core}}(8)$")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks([1, 2, 4, 8])
    ax.set_xticklabels(["1 ($\\mathbb{R}_1$)", "2 ($\\mathbb{C}_2$)", "4 ($\\mathbb{H}_4$)", "8 ($\\mathbb{O}_8$)"], fontsize=9.5)
    ax.set_yticks([1/64, 1/16, 1/4, 1])
    ax.set_yticklabels(["1/64", "1/16", "1/4", "1"], fontsize=9.5)
    ax.set_xlim(0.8, 10.0)
    ax.set_ylim(1/90, 1.8)
    ax.set_xlabel(r"Hurwitz 部分代数次元 $d$（対数）" if ja else r"Hurwitz subalgebra dimension $d$ (log)", fontsize=10.5)
    ax.set_ylabel(r"正規化エネルギー比（対数）" if ja else r"Normalized energy ratio (log)", fontsize=10.5)
    ax.set_title(
        r"(a) 1 層内の厳密べき則：$E_{\mathrm{tot}} \propto d/8$, $E_{\mathrm{core}} \propto (d/8)^2$"
        if ja else
        r"(a) Single-Layer Power Laws: $E_{\mathrm{tot}} \propto d/8$, $E_{\mathrm{core}} \propto (d/8)^2$",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.0, framealpha=0.92)

    # (b) 直積準同型定理 (誤差 2.86e-16) による K 段指数関数的質量階層 r^{-K}
    ax = axes[1]
    K = np.arange(0, 36)
    ax.plot(K, 64.0 ** K, "-", color="#d62728", lw=2.4,
            label=r"コアエネルギー比 $64^K$ ($K=3$: $10^{5.4}$, $K=9$: $10^{16.2}$)" if ja else r"Core energy ratio $64^K$ ($K=3$: $10^{5.4}$, $K=9$: $10^{16.2}$)")
    ax.plot(K, 8.0 ** K, "-", color="#ff7f0e", lw=2.4,
            label=r"コア占有率比 $8^K$ ($K=6$: $10^{5.4}$, $K=18$: $10^{16.3}$)" if ja else r"Core fraction ratio $8^K$ ($K=6$: $10^{5.4}$, $K=18$: $10^{16.3}$)")
    ax.plot(K, 3.0 ** K, "-", color="#1f77b4", lw=2.4,
            label=r"スピン0/2 ギャップ比 $3^K$ ($K=12$: $10^{5.7}$, $K=34$: $10^{16.2}$)" if ja else r"Spin-0/2 gap ratio $3^K$ ($K=12$: $10^{5.7}$, $K=34$: $10^{16.2}$)")
    ax.axhline(3.4e5, color="#555", ls="--", lw=1.3)
    ax.axhline(1.0e16, color="#222", ls=":", lw=1.5)
    ax.text(20.5, 3.4e5 * 2.2, "電子 / トップクォーク比 ($3.4\\times 10^5$)" if ja else "Electron / top quark ratio ($3.4\\times 10^5$)", fontsize=8.2, color="#444")
    ax.text(19.0, 1.0e16 * 2.2, "電弱 / プランク比 ($10^{16}$)" if ja else "Electroweak / Planck ratio ($10^{16}$)", fontsize=8.2, color="#111", fontweight="bold")
    ax.set_yscale("log")
    ax.set_xlim(0, 35)
    ax.set_ylim(1, 1e19)
    ax.set_xlabel(r"テンソル積の入れ子段数 $K$（$O_8^{\otimes K}$）" if ja else r"Tensor-product nesting depth $K$ ($O_8^{\otimes K}$)", fontsize=10.5)
    ax.set_ylabel(r"階層倍率 $r^{-K}$（対数）" if ja else r"Hierarchy magnification $r^{-K}$ (log)", fontsize=10.5)
    ax.set_title(
        r"(b) 直積準同型定理（誤差 $2.86\times 10^{-16}$）による指数階層 $r^{-K}$"
        if ja else
        r"(b) Exponential Hierarchy $r^{-K}$ via Tensor-Product Homomorphism",
        fontsize=10.8, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=7.8, framealpha=0.93)

    fig.suptitle(
        r"1 層内の $(d/8)^2$ 離散階層から $K$ 段入れ子の $10^5 \sim 10^{16}$ 指数階層へ"
        if ja else
        r"From Single-Layer $(d/8)^2$ Discrete Scaling to $10^5$–$10^{16}$ Exponential Hierarchies via Nesting",
        fontsize=12.3, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"results/fig_r6_core_scaling_and_exponential_{lang}.png"
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
