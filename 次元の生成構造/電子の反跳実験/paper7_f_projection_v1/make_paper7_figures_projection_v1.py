#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 図1/2/3（5色）・比較図・横安定性図・λ vs N。PNG+SVG。CSVから生成。解釈なし。

共通横軸: 絶対step 0..55000, 目盛り5000刻み, crossing不動（各N位置に点線）。
図2=5色stackplot(linear)+黒線f(log右), 図3=5色占有(log y, 線)+黒線f。
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
# ★変更1: 原本の paper7_longtime/ を明示（コピー先では CODE.parent が別物になるため）
P7 = (CODE / ".." / ".." / ".." / "時間軸Q軸とフェルミオンの生成構造" / "検証_対照実験"
      / "第5論文原本_自発的分裂予備実験_v1" / "exact_lowN_eigenspectrum_v2"
      / "paper7_longtime").resolve()
# ★変更2: 出力先を本フォルダへ（公開図を上書きしないため）
FIG = CODE / "figures_projection"; FIG.mkdir(parents=True, exist_ok=True)
NS = [5, 40, 300]
XMAX = 55000
COLORS = ["#4C78A8", "#E45756", "#F58518", "#B0B0B0", "#54A24B"]  # P1, d3, d4, other残, 核
LABELS = ["P1 (dominant plane)", "direction 3", "direction 4", "remaining other-rotation", "kernel"]


def load5(n):
    # ★変更3: 両形式CSV（既存16列は公開CSVと一致検証済み）を読む
    f = CODE / f"dual_f_timeseries_N{n:05d}.csv"
    if not f.exists():
        return None
    r = list(csv.DictReader(open(f)))
    t = np.array([float(x["time"]) for x in r])
    P1 = np.array([float(x["direction_1_occupation"]) + float(x["direction_2_occupation"]) for x in r])
    d3 = np.array([float(x["direction_3_occupation"]) for x in r])
    d4 = np.array([float(x["direction_4_occupation"]) for x in r])
    oth = np.array([float(x["other_rotating_occupation"]) for x in r])
    ker = np.array([float(x["kernel_occupation"]) for x in r])
    f_ = np.array([float(x["splitting_fraction_projection"]) for x in r])  # ★変更4: 射影形
    return {"t": t, "bands": [P1, d3, d4, oth, ker], "f": f_}


def crossing(n):
    return json.load(open(P7 / "summary" / f"N{n:05d}_5color_meta.json"))["crossing"]


def setx(ax):
    ax.set_xlim(0, XMAX); ax.set_xticks(np.arange(0, XMAX + 1, 5000))


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=130)
    fig.savefig(FIG / f"{name}.svg")
    plt.close(fig)


def fig1():
    avail = [n for n in NS if load5(n) and (P7 / "summary" / f"N{n:05d}_5color_meta.json").exists()]
    for n in avail:
        d = load5(n); fig, ax = plt.subplots(figsize=(10, 5))
        ax.semilogy(d["t"], d["f"], "k-", lw=1.2)  # ★変更5: 線形→対数
        ax.axvline(crossing(n), color="k", ls=":", lw=0.8)
        setx(ax); ax.set_xlabel("step (absolute)"); ax.set_ylabel("splitting fraction f (projection, log)")
        ax.set_title(f"N={n} Figure1: splitting fraction (common axis 0-55000)")
        save(fig, f"figure1_N{n:05d}")
    # 比較3段
    fig, axes = plt.subplots(len(avail), 1, figsize=(11, 3 * len(avail)), sharex=True, squeeze=False)
    for ax, n in zip(axes[:, 0], avail):
        d = load5(n); ax.semilogy(d["t"], d["f"], "k-", lw=1.1); ax.axvline(crossing(n), color="r", ls=":", lw=0.8)
        setx(ax); ax.set_ylabel(f"N={n}\nf")  # ★変更6: semilogy 化・set_ylim(0,1) 削除
    axes[-1, 0].set_xlabel("step (absolute)")
    fig.suptitle("Figure1 compare: splitting fraction f (common axis)")
    save(fig, "figure1_compare_N5_N40_N300")


def fig23():
    avail = [n for n in NS if load5(n) and (P7 / "summary" / f"N{n:05d}_5color_meta.json").exists()]
    # 図2: stack linear + f log右
    for n in avail:
        d = load5(n); fig, ax = plt.subplots(figsize=(10, 5.5))
        ax.stackplot(d["t"], *d["bands"], labels=LABELS, colors=COLORS, alpha=0.9)
        ax.axvline(crossing(n), color="k", ls=":", lw=0.8)
        setx(ax); ax.set_ylim(0, 1); ax.set_ylabel("occupation fraction")
        ax.set_xlabel("step (absolute)"); ax.legend(loc="center right", fontsize=7)
        ax2 = ax.twinx(); ax2.semilogy(d["t"], d["f"], "k-", lw=1.2)  # ★変更7: clip 削除
        ax2.set_ylabel("f (projection, log)")
        ax.set_title(f"N={n} Figure2 (5-color): P1 / dir3 / dir4 / remaining-other / kernel")
        save(fig, f"figure2_N{n:05d}_5color")
    # 図3: 5色占有 log y（線）+ 黒線f
    for n in avail:
        d = load5(n); fig, ax = plt.subplots(figsize=(10, 5.5))
        for band, lab, c in zip(d["bands"], LABELS, COLORS):
            ax.semilogy(d["t"], np.clip(band, 1e-6, None), lw=0.9, color=c, label=lab)
        ax.semilogy(d["t"], d["f"], "k-", lw=1.0, label="f")  # ★変更8: clip 削除
        ax.axvline(crossing(n), color="k", ls=":", lw=0.8)
        setx(ax); ax.set_xlabel("step (absolute)"); ax.set_ylabel("occupation fraction (log)")
        ax.legend(loc="center right", fontsize=7)
        ax.set_title(f"N={n} Figure3 (5-color, log): occupation per direction")
        save(fig, f"figure3_N{n:05d}_5color")
    # 比較3段（図2 stack）
    for figtag, logy in [("figure2_compare_N5_N40_N300", False), ("figure3_compare_N5_N40_N300", True)]:
        fig, axes = plt.subplots(len(avail), 1, figsize=(11, 3.2 * len(avail)), sharex=True, squeeze=False)
        for ax, n in zip(axes[:, 0], avail):
            d = load5(n)
            if logy:
                for band, c in zip(d["bands"], COLORS):
                    ax.semilogy(d["t"], np.clip(band, 1e-6, None), lw=0.8, color=c)
                ax.semilogy(d["t"], d["f"], "k-", lw=0.9)  # ★変更9: clip 削除
            else:
                ax.stackplot(d["t"], *d["bands"], colors=COLORS, alpha=0.9); ax.set_ylim(0, 1)
            ax.axvline(crossing(n), color="k", ls=":", lw=0.8); setx(ax); ax.set_ylabel(f"N={n}")
        axes[-1, 0].set_xlabel("step (absolute)")
        fig.suptitle(("Figure3 compare (5-color, log)" if logy else "Figure2 compare (5-color stack)")
                     + " — common axis")
        if not logy:
            axes[0, 0].legend(LABELS, fontsize=6, loc="center right")
        save(fig, figtag)


def load_trans(n):
    f = P7 / "raw" / f"N{n:05d}" / "transverse_stability_timeseries.csv"
    if not f.exists():
        return None
    r = list(csv.DictReader(open(f)))
    return r


def fig_transverse():
    avail = [n for n in NS if load_trans(n)]
    # 9.1 log10 A_perp(t)
    for n in avail:
        r = load_trans(n); fig, ax = plt.subplots(figsize=(10, 5.5))
        keys = sorted(set((x["seed"], x["epsilon"]) for x in r))
        for (sd, ep) in keys:
            sub = [x for x in r if x["seed"] == sd and x["epsilon"] == ep and x["local_transverse_growth_rate"] == "nan"]
            if not sub:
                continue
            t = np.array([float(x["time"]) for x in sub]); A = np.array([float(x["normalized_transverse_amplification"]) for x in sub])
            ax.plot(t, np.log10(np.clip(A, 1e-30, None)), lw=0.5)
        setx(ax); ax.set_xlabel("step (absolute)"); ax.set_ylabel("log10 A_perp(t)")
        ax.set_title(f"N={n} transverse amplification (all seed/eps)")
        save(fig, f"transverse_growth_N{n:05d}")
    if avail:
        fig, axes = plt.subplots(len(avail), 1, figsize=(11, 3.2 * len(avail)), sharex=True, squeeze=False)
        for ax, n in zip(axes[:, 0], avail):
            r = load_trans(n)
            for (sd, ep) in sorted(set((x["seed"], x["epsilon"]) for x in r)):
                sub = [x for x in r if x["seed"] == sd and x["epsilon"] == ep and x["local_transverse_growth_rate"] == "nan"]
                if not sub:
                    continue
                t = np.array([float(x["time"]) for x in sub]); A = np.array([float(x["normalized_transverse_amplification"]) for x in sub])
                ax.plot(t, np.log10(np.clip(A, 1e-30, None)), lw=0.4)
            setx(ax); ax.set_ylabel(f"N={n}\nlog10 A_perp")
        axes[-1, 0].set_xlabel("step (absolute)")
        fig.suptitle("Transverse amplification compare (common axis)")
        save(fig, "transverse_growth_compare_N5_N40_N300")


def fig_lambda_vs_N():
    pts = []
    for n in NS:
        f = P7 / "summary" / f"N{n:05d}_transverse_meta.json"
        if f.exists():
            d = json.load(open(f)); pts.append((n, d["lambda_max_for_N"], d.get("lambda_max_normalized")))
    if not pts:
        return
    Nv = np.array([p[0] for p in pts], float); lam = np.array([p[1] for p in pts]); lamn = np.array([p[2] for p in pts], float)
    for xv, xl, name in [(Nv, "N", "lambda_vs_N"), (1 / Nv, "1/N", "lambda_vs_inverse_N"),
                         (1 / np.log(Nv), "1/log N", "lambda_vs_inverse_logN")]:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(xv, lam, "o-", label="lambda_max")
        ax.plot(xv, lamn, "s--", label="lambda_max / sigma1")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel(xl); ax.set_ylabel("lambda_transverse_max"); ax.legend()
        ax.set_title(f"lambda_max vs {xl} (N=5,40,300; no asymptotic law claimed)")
        save(fig, name)


if __name__ == "__main__":
    fig1(); fig23(); fig_transverse(); fig_lambda_vs_N()
    print(f"[図] {FIG} に出力（PNG+SVG）")
