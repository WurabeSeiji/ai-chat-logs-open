#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
図化専用プログラム (計算は一切行わない):
tri_loop_n64_512melt_run_v0.py が保存したデータ
  - results/tri_loop_n64_512melt_timeseries.csv
  - results/tri_loop_n64_512melt_snapshots.npz
を読み込み、以下の 3 枚の図を生成する。

  1. results/fig_n64_512melt_evolution_v0.png
       全体推移 4 パネル (有効成分数 N_eff と非ゼロ数 nz、成長率 ||P|| と結晶崩壊度、
       生存軸数と連結成分数、ステップ間変化 d1/d2)
  2. results/fig_n64_512melt_snapshots_v0.png
       代表 6 時点の 64x64 射影強度マップ sqrt(sum_j |U[i,j,k]|^2)
  3. results/fig_n64_512melt_spectrum_v0.png
       代表時点の振幅スペクトル (|U_{ijk}| 分布) と 64 要素エネルギー分布
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# macOS 日本語フォント設定 (無ければ DejaVu Sans)
plt.rcParams["font.sans-serif"] = ["Hiragino Sans", "Hiragino Kaku Gothic ProN", "Yu Gothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def load_csv(path):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k, v in row.items():
                data.setdefault(k, []).append(v)
    out = {}
    for k, vals in data.items():
        if k == "comp_sizes":
            out[k] = vals
        elif k in ("step", "nz_1e5", "nz_1e3", "nz_1e1", "dim_i", "dim_j", "dim_k", "n_comp", "max_comp"):
            out[k] = np.array([int(x) for x in vals], dtype=np.int64)
        else:
            out[k] = np.array([float(x) for x in vals], dtype=np.float64)
    return out


def plot_evolution(d, out_path):
    steps = d["step"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) 実数 1 マス欠陥 (0,1,0)+1e-4: 巨大結晶の融解と再結合 (step 0〜500)", fontsize=13, fontweight="bold")

    # Panel 1: 非ゼロ数と有効成分数 N_eff
    ax = axes[0, 0]
    ax.plot(steps, d["nz_1e5"], label="非ゼロ数 (|U| > 1e-5)", color="#1f77b4", lw=1.6)
    ax.plot(steps, d["n_eff"], label="有効成分数 N_eff = 1/Σ|U|^4", color="#d62728", lw=1.6, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel("成分数 (対数)")
    ax.set_title("1. 結晶崩壊と塊サイズ (非ゼロ数・有効成分数)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    # Panel 2: 成長率 ||P|| と結晶からの距離 dist(T0)
    ax = axes[0, 1]
    ax.plot(steps, d["norm_P"], label="成長率 ||P||", color="#2ca02c", lw=1.6)
    ax.set_ylabel("成長率 ||P||", color="#2ca02c")
    ax.tick_params(axis="y", labelcolor="#2ca02c")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(steps, d["dist_crys"], label="初期結晶からの距離 dist(T0)", color="#9467bd", lw=1.4, alpha=0.85)
    ax2.set_ylabel("結晶からの距離 dist(T0)", color="#9467bd")
    ax2.tick_params(axis="y", labelcolor="#9467bd")
    ax2.set_ylim(-0.05, 1.08)
    ax.set_title("2. 成長率 ||P|| の段差と初期結晶からの離脱")
    # 凡例統合
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center right", fontsize=8.5)

    # Panel 3: 生存軸数と連結成分 (塊) 数
    ax = axes[1, 0]
    ax.plot(steps, d["dim_i"], label="生存軸数 dim_i (= dim_k)", color="#ff7f0e", lw=1.6)
    ax.plot(steps, d["dim_j"], label="生存軸数 dim_j", color="#17becf", lw=1.3, ls="--")
    ax.plot(steps, d["max_comp"], label="最大連結塊の頂点数", color="#8c564b", lw=1.3, ls=":")
    ax.set_ylabel("軸数 / 頂点数")
    ax.set_xlabel("ステップ数")
    ax.set_title("3. 生存軸数と最大連結塊サイズ")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8.5)

    # Panel 4: 状態変化 d1, d2 と連結塊の個数
    ax = axes[1, 1]
    ax.plot(steps[1:], np.maximum(d["d1"][1:], 1e-16), label="d1 (1ステップ前との距離)", color="#1f77b4", lw=1.1, alpha=0.85)
    ax.plot(steps[2:], np.maximum(d["d2"][2:], 1e-16), label="d2 (2ステップ前との距離)", color="#e377c2", lw=1.1, alpha=0.85)
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 2.0)
    ax.set_ylabel("射影距離 (対数)")
    ax.set_xlabel("ステップ数")
    ax.set_title("4. ステップ間変化 d1, d2 (回転・準安定の推移)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_snapshots(npz, d, out_path):
    chosen = [0, 30, 40, 60, 80, 500]
    avail = set(int(x) for x in npz["snap_steps"])
    chosen = [s for s in chosen if s in avail]
    step_to_idx = {int(s): idx for idx, s in enumerate(d["step"])}

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5))
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) 結晶融解・再結合の 64×64 射影強度マップ  $\sqrt{\Sigma_j |U[i,j,k]|^2}$", fontsize=13, fontweight="bold")

    for ax, st in zip(axes.flat, chosen):
        M = npz[f"map_ik_{st}"]
        # 対数スケールで見やすく (0 は背景黒/濃紺)
        M_log = np.log10(np.maximum(M, 1e-6))
        im = ax.imshow(M_log, origin="upper", cmap="inferno", vmin=-5.5, vmax=0.0, aspect="equal")
        idx = step_to_idx[st]
        nz = d["nz_1e5"][idx]
        neff = d["n_eff"][idx]
        di, dj, dk = d["dim_i"][idx], d["dim_j"][idx], d["dim_k"][idx]
        ax.set_title(f"step {st}: nz={nz}, N_eff={neff:.1f}, 軸=({di},{dj},{dk})", fontsize=10)
        ax.set_xlabel("k (0..63)", fontsize=8.5)
        ax.set_ylabel("i (0..63)", fontsize=8.5)
        ax.set_xticks([0, 16, 32, 48, 63])
        ax.set_yticks([0, 16, 32, 48, 63])

    cbar_ax = fig.add_axes([0.92, 0.12, 0.018, 0.76])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("log10 sqrt(Σ_j |U[i,j,k]|²)", fontsize=9.5)
    plt.tight_layout(rect=[0, 0, 0.90, 0.94])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_spectrum(npz, d, out_path):
    chosen = [0, 30, 40, 80, 500]
    avail = set(int(x) for x in npz["snap_steps"])
    chosen = [s for s in chosen if s in avail]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    fig.suptitle("N = 64 融解・再結合における振幅スペクトルと 64 要素エネルギー分布", fontsize=13, fontweight="bold")

    # Panel 1: 非ゼロ振幅 |U_{ijk}| の分布 (対数ヒストグラム)
    ax = axes[0]
    bins = np.linspace(-5.5, 0.0, 56)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for st, col in zip(chosen, colors):
        vals = npz[f"abs_{st}"]
        if len(vals) > 0:
            ax.hist(np.log10(np.maximum(vals, 1e-6)), bins=bins, histtype="step", lw=1.8, label=f"step {st} (nz={len(vals)})", color=col)
    ax.set_yscale("log")
    ax.set_xlabel("log10 |U_{ijk}|")
    ax.set_ylabel("成分の個数 (対数)")
    ax.set_title("1. 成分振幅スペクトルの変遷")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5)

    # Panel 2: 最終時点 (step 500) と step 80 の 64 要素ごとのエネルギー分布
    ax = axes[1]
    x = np.arange(64)
    e0 = npz["node_e_0"]
    e80 = npz["node_e_80"]
    e500 = npz["node_e_500"]
    ax.plot(x, e0, label="step 0 (初期結晶: 全64要素均等)", color="#7f7f7f", lw=1.2, ls=":")
    ax.bar(x - 0.18, e80, width=0.36, label="step 80", color="#d62728", alpha=0.7)
    ax.bar(x + 0.18, e500, width=0.36, label="step 500", color="#1f77b4", alpha=0.7)
    ax.set_xlabel("要素番号 (0..63)")
    ax.set_ylabel("要素あたりエネルギー (3軸平均)")
    ax.set_title("2. 64 要素へのエネルギー集中 (step 80 vs step 500)")
    ax.set_xticks(np.arange(0, 65, 8))
    ax.set_ylim(0, 0.052)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", ncol=3, fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run():
    d = load_csv(os.path.join("results", "tri_loop_n64_512melt_timeseries.csv"))
    npz = np.load(os.path.join("results", "tri_loop_n64_512melt_snapshots.npz"))

    p1 = os.path.join("results", "fig_n64_512melt_evolution_v0.png")
    p2 = os.path.join("results", "fig_n64_512melt_snapshots_v0.png")
    p3 = os.path.join("results", "fig_n64_512melt_spectrum_v0.png")

    plot_evolution(d, p1)
    plot_snapshots(npz, d, p2)
    plot_spectrum(npz, d, p3)
    print("図を生成しました:")
    for p in (p1, p2, p3):
        print(" ", p)


if __name__ == "__main__":
    run()
