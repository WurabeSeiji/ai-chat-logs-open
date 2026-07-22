#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""note記事「シリーズ続編（第3〜第6論文まとめ）」用の説明図5点を生成する。

図1: 関係の数は爆発、回転はゆっくり、読める方向は三つ（第3論文）
図2: あふれた関係波は回転平面（レジスタ）に整理される（第4論文）
図3: 眠っていた成分の幾何級数的増幅＝自発的分裂（第5論文、実データ）
図4: 分裂の行き先は三つ——体の数と保存量が決める（第5論文）
図5: 存在は無限に開き、読出しだけが飽和する（第6論文、実データ）
"""

import csv
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle

BASE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(BASE, "..", "自発的分裂予備実験_v1")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams["font.family"] = "Hiragino Sans"
plt.rcParams["axes.unicode_minus"] = False

C_BLUE = "#1f77b4"
C_GREEN = "#2ca02c"
C_RED = "#d62728"
C_ORANGE = "#ff7f0e"
C_GRAY = "#888888"


# ---------------------------------------------------------------- 図1
def fig1():
    ns = np.arange(3, 13)
    m = ns * (ns - 1) // 2
    rank = 2 * np.minimum(ns, m // 2)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ns, m, "o-", color=C_BLUE, lw=2, label="関係波の数（体の数の2乗で爆発）")
    ax.plot(ns, rank, "s-", color=C_GREEN, lw=2, label="回転の自由度（線形にしか増えない）")
    ax.axhline(3, color="k", ls="--", lw=1.5, label="一意に読める空間方向（三つで固定）")
    ax.set_xlabel("体の数 N")
    ax.set_ylabel("個数")
    ax.set_title("関係は爆発し、回転はゆっくり増え、読める方向は三つで止まる")
    ax.set_xticks(ns)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left")
    ax.annotate("12体で関係波は66本", xy=(12, 66), xytext=(9.0, 60),
                arrowprops=dict(arrowstyle="->", color=C_BLUE), color=C_BLUE)
    ax.annotate("それでも三方向", xy=(11, 3), xytext=(10.2, 12),
                arrowprops=dict(arrowstyle="->", color="k"))
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_relations_explode_directions_three_v1.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- 図2
def fig2():
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("あふれた関係波は、少数の回転平面に整理される")

    # 左：15本の関係波（N=6）
    rng = np.random.default_rng(7)
    t = np.linspace(0, 2 * np.pi, 120)
    for i in range(15):
        y0 = 0.9 + i * 0.55
        ph = rng.uniform(0, 2 * np.pi)
        ax.plot(0.5 + 2.2 * t / t[-1], y0 + 0.18 * np.sin(3 * t + ph),
                color=C_BLUE, lw=1.0, alpha=0.8)
    ax.text(1.6, 9.5, "関係波 15本（6体の場合）", ha="center", fontsize=11, color=C_BLUE)

    # 矢印
    ax.add_patch(FancyArrowPatch((3.2, 5.0), (5.0, 5.0),
                                 arrowstyle="-|>", mutation_scale=25, color="k"))
    ax.text(4.1, 5.35, "整理", ha="center", fontsize=11)

    # 右：3つの回転平面
    labels = ["回転平面 1", "回転平面 2", "回転平面 3"]
    for k, lab in enumerate(labels):
        cy = 7.6 - k * 2.3
        e = Ellipse((6.9, cy), 3.0, 1.5, facecolor="#e8f0fe",
                    edgecolor=C_GREEN, lw=2)
        ax.add_patch(e)
        th = np.linspace(0, 2 * np.pi, 60)
        ax.plot(6.9 + 1.05 * np.cos(th), cy + 0.5 * np.sin(th),
                color=C_GREEN, lw=1, ls=":")
        ax.add_patch(FancyArrowPatch((6.9 + 1.05, cy), (6.9 + 1.02, cy + 0.09),
                                     arrowstyle="-|>", mutation_scale=14, color=C_GREEN))
        ax.text(6.9, cy, lab, ha="center", va="center", fontsize=11)
    ax.text(6.9, 1.0, "各平面の中身は、二つの波の関係（AB）と同じかたち\n＝ひとまとまりの「粒子のような」単位",
            ha="center", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_planes_as_registers_v1.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- 図3
def _load_largeN_curves(wanted=(300, 1000)):
    """largeN_splitting_result_v1 から指定Nの f 曲線を読む（存在するもののみ、N昇順）。"""
    import glob
    d = os.path.join(EXP, "largeN_splitting_result_v1")
    found = {}
    for path in glob.glob(os.path.join(d, "fcurve_N*.csv")):
        n = int(os.path.basename(path).split("_")[1][1:])
        if n in wanted:
            found[n] = path
    out = []
    for n in sorted(found):
        taus, fs = [], []
        with open(found[n]) as fh:
            for row in csv.DictReader(fh):
                taus.append(float(row["tau"]))
                fs.append(float(row["f"]))
        out.append((n, np.array(taus), np.array(fs)))
    return out


def fig3():
    # 小N（3体）の実データ
    path = os.path.join(EXP, "spontaneous_splitting_result_v1", "dormant_fraction_curves_v1.csv")
    taus3, f3 = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            taus3.append(float(row["tau"]))
            f3.append(float(row["f_median_delta_1e-05"]))
    taus3 = np.array(taus3)
    f3 = np.array(f3)

    # 大Nの実データ（300体・1000体）
    curves = _load_largeN_curves()

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    ax.semilogy(taus3, f3, color=C_BLUE, lw=1.6,
                label="3体（関係波3本、種はエネルギー比で100億分の1）")
    big_colors = {300: C_ORANGE, 1000: C_RED}
    for n_big, taus_b, f_b in curves:
        m_big = n_big * (n_big - 1) // 2
        ax.semilogy(taus_b, f_b, color=big_colors.get(n_big, C_GREEN), lw=1.8,
                    label=f"{n_big}体（関係波{m_big:,}本、種は10のマイナス30乗）")
    ax.set_xlabel("時間（ステップ）")
    ax.set_ylabel("新しい波に移った割合（対数目盛）")
    ax.set_title("ほぼゼロだった成分が、幾何級数で立ち上がる（実データ）")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=10)
    ax.set_ylim(1e-31, 2.0)
    if curves:
        n_big, taus_b, f_b = curves[-1]
        i_mid = int(np.argmin(np.abs(f_b - 1e-15)))
        ax.annotate("20桁以上を同じ倍率で\nまっすぐ駆け上がる",
                    xy=(taus_b[i_mid], f_b[i_mid]), xytext=(taus_b[i_mid] + 1100, 1e-19),
                    arrowprops=dict(arrowstyle="->", color="k"), fontsize=10.5)
        ax.annotate("出発点は10のマイナス30乗\n（測定にかからない「眠った」成分）",
                    xy=(taus_b[0], f_b[0]), xytext=(500, 3e-28),
                    arrowprops=dict(arrowstyle="->", color="k"), fontsize=10.5)
    ax.annotate("一定の割合で安定\n＝新しい波が定着",
                xy=(taus3[-80], f3[-80]), xytext=(700, 4e-5),
                arrowprops=dict(arrowstyle="->", color="k"), fontsize=10.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_dormant_growth_v1.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- 図4
def fig4():
    # 実験O4（閉鎖ゼロの層）と実験O5（保存量が大きい層）の帰結クラス
    rows = ["閉鎖の保存量が\nゼロの層", "閉鎖の保存量が\n大きい層"]
    ns = [3, 4, 5, 6]
    classes = [
        ["有界混合", "有界混合", "拡大", "拡大"],
        ["回帰", "拡大", "拡大", "拡大"],
    ]
    color = {"拡大": "#f4b6b6", "有界混合": "#d9d9d9", "回帰": "#b6c8f4"}

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.set_xlim(0, len(ns))
    ax.set_ylim(0, len(rows))
    for i, row in enumerate(classes):
        for j, cls in enumerate(row):
            y = len(rows) - 1 - i
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor=color[cls],
                                   edgecolor="white", lw=3))
            ax.text(j + 0.5, y + 0.5, cls, ha="center", va="center", fontsize=12)
    ax.set_xticks([j + 0.5 for j in range(len(ns))])
    ax.set_xticklabels([f"{n} 体" for n in ns], fontsize=11)
    ax.set_yticks([len(rows) - 1 - i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels(rows, fontsize=10.5)
    ax.set_title("分裂の行き先は三つ——体の数と、保存量の値が決める")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.text(0.5, 0.015,
             "拡大＝新しい波が定着して広がる ／ 有界混合＝広がるが定着しない ／ 回帰＝元の波の近くへ戻り続ける",
             ha="center", fontsize=9.5)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(os.path.join(FIG, "fig4_three_fates_v1.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- 図5
def fig5():
    with open(os.path.join(EXP, "counting_ceiling_result_v1", "summary_v1.json")) as fh:
        d = json.load(fh)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # 左：存在の層——同じ閉鎖は何項にでも書き換えられる
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.axis("off")
    ax1.set_title("存在の層：波の数に天井はない")
    levels = [(8.6, "1 項", 1), (6.8, "2 項", 2), (5.0, "4 項", 4), (3.2, "8 項", 8)]
    t = np.linspace(0, 2 * np.pi, 100)
    for y0, lab, k in levels:
        ax1.text(0.6, y0, lab, fontsize=11, va="center")
        width = 7.0 / k
        for i in range(k):
            x0 = 2.2 + i * width
            ax1.plot(x0 + width * 0.85 * t / t[-1],
                     y0 + 0.35 * np.sin(k * t / 2 + i), color=C_BLUE, lw=1.4)
    for ya, yb in [(8.2, 7.4), (6.4, 5.6), (4.6, 3.8)]:
        ax1.add_patch(FancyArrowPatch((1.2, ya), (1.2, yb),
                                      arrowstyle="-|>", mutation_scale=14, color=C_GRAY))
    ax1.text(5.5, 1.9, "…… 何項にでも書き換えられる（無限に続く）", fontsize=11, ha="center")
    ax1.text(5.5, 0.7, "どれも同じ一つの閉鎖の正当な読み方", fontsize=10.5, ha="center", color=C_GRAY)

    # 右：読出しの層——上限があり、実際はその1〜4割で止まる（実データ）
    u3 = d["u3"]
    ns = sorted(int(k) for k in u3)
    filled = {n: {fam: u3[str(n)][fam] / u3[str(n)]["ceiling"]
                  for fam in ("thermal", "frozen", "parent")} for n in ns}
    fams = [("thermal", "波が混ざり合った状態", C_BLUE),
            ("frozen", "凍結した状態", C_ORANGE),
            ("parent", "単一の波", C_GREEN)]
    x = np.arange(len(ns))
    w = 0.26
    for k, (fam, lab, col) in enumerate(fams):
        vals = [filled[n].get(fam, 0.0) for n in ns]
        ax2.bar(x + (k - 1) * w, vals, w, color=col, label=lab)
    ax2.axhline(1.0, color=C_RED, ls="--", lw=1.8)
    ax2.text(len(ns) - 0.55, 1.03, "規約が許す上限（天井）", color=C_RED, fontsize=10.5, ha="right")
    ax2.set_ylim(0, 1.15)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{n} 体" for n in ns])
    ax2.set_ylabel("読出し数 ÷ 上限")
    ax2.set_title("読出しの層：上限のさらに1〜4割で止まる（実データ）")
    ax2.legend(loc="upper left", fontsize=9.5)
    ax2.grid(axis="y", alpha=0.3)

    fig.suptitle("存在は無限に開き、読出しだけが飽和する", fontsize=14, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(FIG, "fig5_existence_vs_readout_v1.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    print("figures written to", FIG)
