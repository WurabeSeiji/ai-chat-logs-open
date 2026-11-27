#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
図化専用プログラム (計算は一切行わない):
tri_loop_n64_triangle_viewangle_run_v0.py が保存したデータ
  - results/tri_loop_n64_viewangle_timeseries.csv
  - results/tri_loop_n64_viewangle_snapshots.npz
を読み込み、以下の 2 枚の図を生成する。

  1. results/fig_n64_viewangle_vs_condensation_v0.png
       粒子凝縮 (E_512, N_eff) と、粒子/背景の指数部差 ln(粒子/背景)、
       三角形ループの視野角 ρ = L3/L1 および引き伸ばし Δln = ln(L1/L3) の対比 (4 パネル)
  2. results/fig_n64_viewangle_spectrum_v0.png
       代表時点 (step 0, 30, 40, 60, 70, 80) における三角形視野角 log10(L3/L1) の
       スペクトル比較 (粒子内部の三角形 vs 背景の三角形)
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Hiragino Sans", "Hiragino Kaku Gothic ProN", "Yu Gothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def load_csv(path):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            for k, v in row.items():
                data.setdefault(k, []).append(v)
    out = {}
    for k, vals in data.items():
        if k in ("step", "tri_part_cnt", "tri_melt_cnt"):
            out[k] = np.array([int(x) for x in vals], dtype=np.int64)
        else:
            out[k] = np.array([float(x) for x in vals], dtype=np.float64)
    return out


def plot_vs_condensation(d, out_path):
    st = d["step"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) 粒子の凝縮と三角形視野角 ($\rho = L_3/L_1$)・指数部の引き伸ばし (step 0〜100)", fontsize=13, fontweight="bold")

    # Panel 1: 粒子凝縮 (E_512 と N_eff)
    ax = axes[0, 0]
    ax.plot(st, d["e_512"], label="512マスへのエネルギー集中率 E_512", color="#d62728", lw=2.0)
    ax.set_ylabel("エネルギー集中率 E_512", color="#d62728")
    ax.tick_params(axis="y", labelcolor="#d62728")
    ax.set_ylim(-0.05, 1.08)
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(st, d["n_eff"], label="有効成分数 N_eff", color="#1f77b4", lw=1.5, ls="--")
    ax2.set_yscale("log")
    ax2.set_ylabel("有効成分数 N_eff (対数)", color="#1f77b4")
    ax2.tick_params(axis="y", labelcolor="#1f77b4")
    ax.set_title("1. 粒子の凝縮 (512マス集中率と有効成分数)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center left", fontsize=8.5)

    # Panel 2: 粒子/背景の指数部の差 ln(<|U_512|> / <|U_melt|>)
    ax = axes[0, 1]
    ax.plot(st, d["dln_contrast"], label="ln(粒子平均振幅 / 融解背景平均振幅)", color="#2ca02c", lw=2.0)
    ax.set_ylabel("指数部の差 Δln (粒子 vs 背景)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_title("2. 粒子と溶け落ちた背景の指数部コントラスト (急拡大と飽和)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    # Panel 3: 三角形の視野角 ρ = L3 / L1 (対数: 1 が正三角形、小さいほど細長い)
    ax = axes[1, 0]
    ax.plot(st, d["tri_part_rho_med"], label="粒子代表の三角形 (中央値)", color="#d62728", lw=1.8)
    ax.fill_between(st, np.maximum(d["tri_part_rho_p10"], 1e-30), d["tri_part_rho_p90"], color="#d62728", alpha=0.15, label="粒子側 10%〜90% 帯")
    ax.plot(st, d["tri_melt_rho_med"], label="背景代表の三角形 (中央値)", color="#1f77b4", lw=1.8, ls="--")
    ax.fill_between(st, np.maximum(d["tri_melt_rho_p10"], 1e-30), d["tri_melt_rho_p90"], color="#1f77b4", alpha=0.15, label="背景側 10%〜90% 帯")
    ax.set_yscale("log")
    ax.set_ylim(1e-25, 2.0)
    ax.set_xlabel("ステップ数")
    ax.set_ylabel("視野角指標 ρ = L3 / L1 (対数)")
    ax.set_title("3. 三角形の視野角 ρ = L3/L1 (1=正三角形, 急減=遠ざかり)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)

    # Panel 4: 実質的三角形 (3辺 > 1e-12) の指数部引き伸ばし Δln = ln(L1/L3)
    ax = axes[1, 1]
    ax.plot(st, d["tri_part_sig_dln_med"], label="粒子代表 (3辺>1e-12 中央値)", color="#d62728", lw=1.8)
    ax.plot(st, d["tri_melt_sig_dln_med"], label="背景代表 (3辺>1e-12 中央値)", color="#1f77b4", lw=1.8, ls="--")
    ax.plot(st, d["tri_part_dln_med"], label="粒子代表 (全非ゼロ中央値)", color="#ff7f0e", lw=1.2, ls=":")
    ax.plot(st, d["tri_melt_dln_med"], label="背景代表 (全非ゼロ中央値)", color="#17becf", lw=1.2, ls=":")
    ax.set_yscale("symlog", linthresh=0.5)
    ax.set_xlabel("ステップ数")
    ax.set_ylabel("指数部の引き伸ばし Δln = ln(L1/L3)")
    ax.set_title("4. 指数部での引き伸ばし Δln (正三角形=0)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_spectrum(npz, out_path):
    chosen = [0, 30, 40, 60, 70, 80]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) 三角形視野角 $\log_{10}(L_3 / L_1)$ スペクトルの変遷 (0 = 正三角形)", fontsize=13, fontweight="bold")

    bins = np.linspace(-25.0, 0.0, 60)
    colors = ["#7f7f7f", "#ff7f0e", "#2ca02c", "#1f77b4", "#9467bd", "#d62728"]

    for ax, prefix, title in [
        (axes[0], "rho_part", "1. 粒子代表位置の三角形スペクトル"),
        (axes[1], "rho_melt", "2. 溶け落ちた背景代表位置の三角形スペクトル"),
    ]:
        for st, col in zip(chosen, colors):
            vals = npz[f"{prefix}_{st}"]
            if len(vals) > 0:
                logv = np.log10(np.maximum(vals, 1e-25))
                ax.hist(logv, bins=bins, histtype="step", lw=1.7, label=f"step {st} (n={len(vals)})", color=col)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\log_{10}(\rho) = \log_{10}(L_3 / L_1)$   (0: 正三角形, 左ほど細長い)")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=8.5)
    axes[0].set_ylabel("三角形の個数 (対数)")

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run():
    d = load_csv(os.path.join("results", "tri_loop_n64_viewangle_timeseries.csv"))
    npz = np.load(os.path.join("results", "tri_loop_n64_viewangle_snapshots.npz"))

    p1 = os.path.join("results", "fig_n64_viewangle_vs_condensation_v0.png")
    p2 = os.path.join("results", "fig_n64_viewangle_spectrum_v0.png")
    plot_vs_condensation(d, p1)
    plot_spectrum(npz, p2)
    print("図を生成しました:")
    for p in (p1, p2):
        print(" ", p)


if __name__ == "__main__":
    run()
