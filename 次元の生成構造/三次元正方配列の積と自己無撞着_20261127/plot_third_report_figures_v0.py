#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 3 速報用の図 2 種（日本語版・英語版、計 4 枚）を生成する専用スクリプト（計算なし）。
実測データソース:
  - results/tri_loop_n64_spin_gauge_higgs_v0.txt
  - results/tri_loop_n64_harmonic_localization_v0.txt
  - results/tri_loop_n64_g2_and_twoparticle_v0.txt
出力:
  - results/fig_r3_spin_and_gauge_ja.png / _en.png
  - results/fig_r3_harmonic_and_2p_ja.png / _en.png
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
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), dpi=200)

    # (a) SO(3) スピン 0 / 1 / 2 励起応答 (τ = 0..10)
    ax = axes[0, 0]
    taus = np.array([0, 1, 2, 3, 5, 10])
    s1_vals = np.ones_like(taus, dtype=float)
    s0_vals = np.array([1.0, 0.333, 0.111, 0.037, 0.004, 1.0 / (3 ** 10)])
    t_line = np.linspace(0, 10, 100)
    ax.plot(taus, s1_vals, "o-", color="#1f77b4", lw=2.4, ms=7,
            label="スピン1: SO(3)回転 Lx,Ly,Lz (3重項, |λ|=1.000 無質量)" if ja else "Spin-1: SO(3) Lx,Ly,Lz (triplet, |λ|=1.000 massless)")
    ax.plot(t_line, (1.0 / 3.0) ** t_line, "--", color="#d62728", lw=1.6,
            label=r"理論線 $(1/3)^\tau$ (質量ギャップ $\ln 3$)" if ja else r"Exact $(1/3)^\tau$ (mass gap $\ln 3$)")
    ax.plot(taus[:5], s0_vals[:5], "s", color="#d62728", ms=8,
            label="スピン0 (2種) & スピン2 (5重項): 全7モード実測一致" if ja else "Spin-0 (2) & Spin-2 (5): all 7 modes match")
    ax.set_yscale("log")
    ax.set_xlim(-0.3, 10.3)
    ax.set_ylim(1e-4, 2.5)
    ax.set_xlabel(r"更新ステップ $\tau$" if ja else r"Update step $\tau$", fontsize=11)
    ax.set_ylabel(r"相対振幅 $|\delta(\tau)|/|\delta(0)|$ (対数)" if ja else r"Relative amplitude $|\delta(\tau)|/|\delta(0)|$ (log)", fontsize=11)
    ax.set_title("(a) $SO(3)$ スピン 0 / 1 / 2 励起の厳密選別則" if ja else "(a) Exact $SO(3)$ Spin 0 / 1 / 2 Selection", fontsize=11.0, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.5, framealpha=0.92)

    # (b) 別チャネル b={2,3} 上のスピン s=0,1,2 の 1/sqrt(s+1) 法則
    ax = axes[0, 1]
    s_arr = np.array([0, 1, 2])
    meas = np.array([1.000, 0.707, 0.577])
    theo = 1.0 / np.sqrt(s_arr + 1.0)
    s_cont = np.linspace(0, 2.2, 100)
    ax.plot(s_cont, 1.0 / np.sqrt(s_cont + 1.0), "--", color="#2ca02c", lw=2.0,
            label=r"理論曲線 $1/\sqrt{s+1}$" if ja else r"Exact curve $1/\sqrt{s+1}$")
    ax.bar(s_arr, meas, width=0.38, color=["#1f77b4", "#ff7f0e", "#9467bd"], edgecolor="#333", alpha=0.85)
    for s_i, m_i in zip(s_arr, meas):
        ax.text(s_i, m_i + 0.03, f"{m_i:.3f}\n($1/\\sqrt{{{s_i+1}}}$)", ha="center", va="bottom", fontsize=9.5, fontweight="bold")
    ax.set_xticks(s_arr)
    ax.set_xticklabels(["Spin 0\n(スカラー)" if ja else "Spin 0\n(Scalar)",
                        "Spin 1\n(ベクトル)" if ja else "Spin 1\n(Vector)",
                        "Spin 2\n(テンソル)" if ja else "Spin 2\n(Tensor)"], fontsize=10)
    ax.set_ylim(0, 1.22)
    ax.set_ylabel(r"$\tau=10$ 安定振幅比" if ja else r"Stable amplitude at $\tau=10$", fontsize=11)
    ax.set_title(r"(b) 別チャネル ($b=\{2,3\}$) の離散 $1/\sqrt{s+1}$ 法則" if ja else r"(b) Exact $1/\sqrt{s+1}$ Law in $b=\{2,3\}$ Channel", fontsize=11.0, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.92)

    # (c) 外側 O8 の厳密ゲージ対称性次元カウント (G2=14, SU(3)=8)
    ax = axes[1, 0]
    groups = ["SO(8)\n(28)", "SO(7)\n(21)", "SO(6)\n[A=7固定](15)" if ja else "SO(6)\n[A=7 fixed](15)", "SO(3)\n(3)"]
    inv_dims = np.array([14, 14, 8, 0])
    brk_dims = np.array([14, 7, 7, 3])
    y = np.arange(len(groups))
    ax.barh(y, inv_dims, color="#1f77b4", edgecolor="#333", label="厳密不変ゲージ対称性 (G₂=14, SU(3)=8)" if ja else "Exact gauge symmetry (G₂=14, SU(3)=8)")
    ax.barh(y, brk_dims, left=inv_dims, color="#ffbb78", edgecolor="#333", label="破れ・無質量モード (|λ|=1.000)" if ja else "Broken massless modes (|λ|=1.000)")
    for yi, iv, bv in zip(y, inv_dims, brk_dims):
        if iv > 0:
            ax.text(iv / 2, yi, f"{iv}", ha="center", va="center", color="white", fontweight="bold", fontsize=10.5)
        ax.text(iv + bv / 2, yi, f"{bv}", ha="center", va="center", color="#222", fontweight="bold", fontsize=10.5)
    ax.set_yticks(y)
    ax.set_yticklabels(groups, fontsize=9.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 31)
    ax.set_xlabel("生成子の次元数" if ja else "Number of generators", fontsize=11)
    ax.set_title(r"(c) 外側 $O_8$ の厳密ゲージ本数 ($G_2=14$, $SU(3)=8$)" if ja else r"(c) Exact Gauge Counts on Outer $O_8$ ($G_2=14, SU(3)=8$)", fontsize=11.0, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.92)

    # (d) 第 4 層 (Case 1, τ=1,000) 64×64 厳密ヤコビアン固有値分類
    ax = axes[1, 1]
    cats = ["不安定\n(|λ|>1)" if ja else "Unstable\n(|λ|>1)",
            "無質量平坦\n(|λ|=1.000)" if ja else "Massless flat\n(|λ|=1.000)",
            "中間復元\n(0<|λ|<1)" if ja else "Intermediate\n(0<|λ|<1)",
            "即時遮断\n(|λ|=0)" if ja else "Expelled\n(|λ|=0)"]
    cnts = [0, 9, 0, 55]
    cols = ["#d62728", "#1f77b4", "#ff7f0e", "#7f7f7f"]
    bars = ax.bar(cats, cnts, color=cols, edgecolor="#333", width=0.52)
    for b, c in zip(bars, cnts):
        ax.text(b.get_x() + b.get_width() / 2, c + 1.2, f"{c}", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_ylim(0, 64)
    ax.set_ylabel("固有モード数（全 64 次元）" if ja else "Eigenmode count (total 64)", fontsize=11)
    ax.set_title(r"(d) 第 4 層静止固定点の 64×64 厳密ヤコビアン" if ja else r"(d) Exact 64×64 Jacobian of Layer 4 Fixed Point", fontsize=11.0, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.annotate("S³×S³×S³ 回転\n(3+3+3 = 9 本)" if ja else "S³×S³×S³ rotations\n(3+3+3 = 9 modes)",
                xy=(1, 9), xytext=(1.35, 28),
                arrowprops=dict(arrowstyle="->", color="#1f77b4"), fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#1f77b4", alpha=0.9))

    fig.suptitle(
        r"三角形ループ積における $SO(3)$ スピン選別、$G_2/SU(3)$ ゲージ対称性、および第 4 層厳密ヤコビアン"
        if ja else
        r"$SO(3)$ Spin Selection, $G_2/SU(3)$ Gauge Symmetries, and Exact Layer-4 Jacobian in the Triangle-Loop Product",
        fontsize=13.0, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join("results", f"fig_r3_spin_and_gauge_{lang}.png")
    plt.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


def plot_fig2(lang):
    setup_font(lang)
    ja = (lang == "ja")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=200)

    # (a) 空間倍音読み (WHT k=0,1,2,3)
    ax = axes[0]
    labels = [
        "背景 H\n(非局在)" if ja else "Background H\n(Nonlocal)",
        "大域波 P3a\n(非局在)" if ja else "Global P3a\n(Nonlocal)",
        "ゲージ波 Lx\n(ほぼ非局在)" if ja else "Gauge Lx\n(Near-nonlocal)",
        "四元数塊 P5\n(4サイト局在)" if ja else "H4 Lump P5\n(4-site local)",
        "単一サイト P4\n(1サイト極局在)" if ja else "1-Site P4\n(Point-local)",
    ]
    k0 = np.array([100.0, 100.0, 96.34, 90.00, 54.15])
    k1 = np.array([0.0, 0.0, 2.44, 10.00, 19.65])
    k2 = np.array([0.0, 0.0, 1.22, 0.00, 19.65])
    k3 = np.array([0.0, 0.0, 0.00, 0.00, 6.55])
    x = np.arange(len(labels))
    w = 0.55
    ax.bar(x, k0, width=w, color="#aec7e8", edgecolor="#333", label="基本 DC (k=0)" if ja else "Fundamental DC (k=0)")
    ax.bar(x, k1, bottom=k0, width=w, color="#ff7f0e", edgecolor="#333", label="1次空間倍音 (k=1)" if ja else "1st harmonic (k=1)")
    ax.bar(x, k2, bottom=k0 + k1, width=w, color="#d62728", edgecolor="#333", label="2次空間倍音 (k=2)" if ja else "2nd harmonic (k=2)")
    ax.bar(x, k3, bottom=k0 + k1 + k2, width=w, color="#9467bd", edgecolor="#333", label="3次空間倍音 (k=3)" if ja else "3rd harmonic (k=3)")
    tot_h = k1 + k2 + k3
    for xi, th in zip(x, tot_h):
        ax.text(xi, 102.0, f"倍音 {th:.1f}%" if ja else f"Harm {th:.1f}%", ha="center", va="bottom", fontsize=9, fontweight="bold",
                color="#d62728" if th > 9 else "#444")
    ax.annotate("厳密 3:3:1\n二項比!" if ja else "Exact 3:3:1\nbinomial!",
                xy=(4, 77), xytext=(2.65, 68),
                arrowprops=dict(arrowstyle="->", color="#d62728"), fontsize=9, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#d62728", alpha=0.92))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.8)
    ax.set_ylim(0, 114)
    ax.set_ylabel("1 軸エネルギーの空間倍音占有率 (%)" if ja else "Spatial harmonic energy share (%)", fontsize=10.5)
    ax.set_title(r"(a) $\mathbb{Z}_2^3$ 空間倍音読み：非局在波 vs 局在凝縮体" if ja else r"(a) $\mathbb{Z}_2^3$ Spatial Harmonics: Nonlocal vs Localized", fontsize=11.0, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.3, framealpha=0.92)

    # (b) 2 局在粒子 (P4 × 2) の共存・対消滅比較 (τ=30)
    ax = axes[1]
    c_labels = [
        "(a1) 同位置\n同ch 同相(+)" if ja else "(a1) Same site\nSame ch (+)",
        "(a2) 同位置\n同ch 逆相(-)" if ja else "(a2) Same site\nSame ch (-)",
        "(b1) 別位置\n同ch 同相(+)" if ja else "(b1) Diff sites\nSame ch (+)",
        "(b2) 別位置\n同ch 逆相(-)" if ja else "(b2) Diff sites\nSame ch (-)",
        "(c1) 同位置\n別ch (23&67)" if ja else "(c1) Same site\nDiff ch",
        "(c2) 別位置\n別ch (23&67)" if ja else "(c2) Diff sites\nDiff ch",
    ]
    e1 = np.array([0.243, 0.0, 0.094, 0.094, 0.084, 0.094])
    e2 = np.array([0.243, 0.0, 0.094, 0.094, 0.084, 0.094])
    ec = np.array([1.266, 0.0, 1.196, 1.196, 0.970, 1.196])
    x2 = np.arange(len(c_labels))
    w2 = 0.25
    ax.bar(x2 - w2, e1, width=w2, color="#1f77b4", edgecolor="#333", label=r"粒子1 本体 $E_1$" if ja else r"Particle 1 core $E_1$")
    ax.bar(x2, e2, width=w2, color="#2ca02c", edgecolor="#333", label=r"粒子2 本体 $E_2$" if ja else r"Particle 2 core $E_2$")
    ax.bar(x2 + w2, ec, width=w2, color="#ffbb78", edgecolor="#333", label=r"共有の衣 $E_{\mathrm{oth}}$" if ja else r"Shared cloud $E_{\mathrm{oth}}$")
    ax.text(1, 0.08, "対消滅\n(nz=512)" if ja else "Annihilation\n(nz=512)", ha="center", va="bottom", fontsize=8.5, color="#d62728", fontweight="bold")
    for xi in (2, 3, 4, 5):
        ax.text(xi, ec[xi] + 0.04, "1:1共存\n(nz=640)" if ja else "1:1 Coexist\n(nz=640)", ha="center", va="bottom", fontsize=8.0, color="#1f77b4", fontweight="bold")
    ax.set_xticks(x2)
    ax.set_xticklabels(c_labels, fontsize=8.2)
    ax.set_ylim(0, 1.72)
    ax.set_ylabel(r"$\tau=30$ エネルギー ($\times 10^{-3}$)" if ja else r"Energy at $\tau=30$ ($\times 10^{-3}$)", fontsize=10.5)
    ax.set_title(r"(b) 背景 $H$ 上での 2 局在粒子 (P4×2) の完全共存と対消滅" if ja else r"(b) Two Localized Particles (P4×2) on $H$: Coexistence & Annihilation", fontsize=11.0, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper center", ncol=3, fontsize=8.2, framealpha=0.92)

    fig.suptitle(
        r"空間倍音による局在粒子の形成（3:3:1 二項比）と 2 局在粒子の完全共存"
        if ja else
        r"Spatial Harmonic Localization (Exact 3:3:1 Binomial Ratio) and Two-Particle Coexistence on Background $H$",
        fontsize=13.0, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join("results", f"fig_r3_harmonic_and_2p_{lang}.png")
    plt.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    for lg in ("ja", "en"):
        plot_fig1(lg)
        plot_fig2(lg)
