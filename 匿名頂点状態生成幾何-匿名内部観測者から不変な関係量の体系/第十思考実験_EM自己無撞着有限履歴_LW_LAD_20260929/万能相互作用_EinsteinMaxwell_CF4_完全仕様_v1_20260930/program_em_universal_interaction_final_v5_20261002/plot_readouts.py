#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_readouts.py — 節点写像の走行結果（CSV + summary JSON）を仕様 §7 の読み出しに沿って図にする。
仕様: ../coefficients_relational_full_list_ja_final_v6_20261002.md §2.2・§7

図（SVG と PNG を保存）
  F1 ロック      |S₁|、|S₂|（放射の源。0 に収束すればロック）
  F2 エネルギー   Δε（力学的エネルギーの変化）と ∫P_rad dθ（Larmor の放射損失）— 同じ単位、同じ軸
  F3 倍音分布    〈m²〉 − 1 の時間変化、最終 |ψ_m|²
  F4 呼吸        q = A ζ_c（呼吸の振幅比）
  F5 時計の比    γ_k − 1、(v/c)/α − 1、dθ_a/dθ_Y（a のサイクル/周）— 節点時計での比。時計因子の有無は仕様 §2.2
  F6 位相        arg B − θ_X（呼吸位相のずれ）、ψ₁ の位相の進み／周（π 単位）
  F7 倍音の重み  |ψ₀|²、|ψ₁|²、|ψ₂|²（〈m〉・〈m²〉から復元）と S₁ の二つの成分 |ψ₀ψ₁|、|ψ₁ψ₂|
横軸はすべて共通時計：呼吸の周（節点/62）。

表示規約（dataviz）：一つの軸に一つの単位（二重軸なし）、細い線、控えめな格子、2 系列以上は凡例、
文字は文字色（系列色を文字に使わない）、色の順序は固定（青・橙・水色）。
記録点が多いとき（> 5000）は窓ごとの平均線と最小–最大の帯で描く。
"""
import argparse
import csv
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
import gzip  # noqa: E402

# ---- 配色（dataviz 参照パレット、light。検証済み：adjacent CVD ΔE ≥ 8、normal ≥ 15） ----
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]      # 1 青、2 橙、3 水色（固定順）
INK = "#0b0b0b"
INK2 = "#52514e"
INK3 = "#8a8984"
GRID = "#e6e5e2"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
SEQ_BLUE = "#2a78d6"
ALPHA_FS = 7.2973525643e-3


def setup_fonts():
    names = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ["Hiragino Sans", "Hiragino Kaku Gothic ProN", "Noto Sans CJK JP", "YuGothic", "Osaka"]:
        if cand in names:
            plt.rcParams["font.family"] = cand
            break
    plt.rcParams.update({
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.titlecolor": INK, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.grid": True, "axes.grid.axis": "y", "axes.spines.top": False, "axes.spines.right": False,
        "axes.facecolor": SURFACE, "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "lines.linewidth": 1.6, "font.size": 10, "axes.titlesize": 11, "legend.frameon": False,
        "axes.unicode_minus": False, "svg.fonttype": "none",
    })


def load_run(csv_path):
    opener = gzip.open if csv_path.endswith(".gz") else open
    with opener(csv_path, "rt") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) > 3:
        ks = [int(float(r["k"])) for r in rows]
        step = ks[1] - ks[0]
        if step > 1 and (ks[-1] - ks[-2]) != step:
            rows = rows[:-1]                       # 最終行は記録間隔から外れた節点（呼吸位相が違う）なので落とす
    data = {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}
    summ_path = os.path.splitext(csv_path[:-3] if csv_path.endswith(".gz") else csv_path)[0] + "_summary.json"
    summ = json.load(open(summ_path)) if os.path.exists(summ_path) else {}
    return data, summ


def binned(x, y, nbins=2000):
    """記録点が多いときの窓集約：窓の中心 x、平均 y、最小 y、最大 y。"""
    n = len(x)
    if n <= 5000:
        return x, y, None, None
    edges = np.linspace(0, n, nbins + 1).astype(int)
    xc = np.array([x[a:b].mean() for a, b in zip(edges[:-1], edges[1:]) if b > a])
    ym = np.array([y[a:b].mean() for a, b in zip(edges[:-1], edges[1:]) if b > a])
    ylo = np.array([y[a:b].min() for a, b in zip(edges[:-1], edges[1:]) if b > a])
    yhi = np.array([y[a:b].max() for a, b in zip(edges[:-1], edges[1:]) if b > a])
    return xc, ym, ylo, yhi


def draw_series(ax, x, y, color, label=None):
    xc, ym, ylo, yhi = binned(x, y)
    if ylo is not None:
        ax.fill_between(xc, ylo, yhi, color=color, alpha=0.18, linewidth=0)
    ax.plot(xc, ym, color=color, label=label)
    return xc, ym


def end_label(ax, x, y, text):
    ax.annotate(text, xy=(x[-1], y[-1]), xytext=(4, 0), textcoords="offset points", va="center",
                fontsize=9, color=INK2)


def finish(fig, out_base):
    for a in fig.get_axes():
        a.xaxis.set_major_locator(MaxNLocator(nbins=6))
    fig.tight_layout()
    fig.savefig(out_base + ".svg")
    fig.savefig(out_base + ".png", dpi=150)
    plt.close(fig)
    return [out_base + ".svg", out_base + ".png"]


def run_label(summ):
    ini = summ.get("init", {})
    return f"δ = {ini.get('delta', '?')}, f = {ini.get('f', '?')}, スピン {ini.get('sa', '?')}/{ini.get('sb', '?')}, 節点 {ini.get('nodes', '?')}"


def main():
    ap = argparse.ArgumentParser(description="走行結果の図化（仕様 §7 の読み出し）")
    ap.add_argument("csv", nargs="+")
    ap.add_argument("--out", default="results/figures")
    args = ap.parse_args()
    setup_fonts()
    os.makedirs(args.out, exist_ok=True)
    written = []
    for path in args.csv:
        d, summ = load_run(path)
        name = os.path.splitext(os.path.basename(path[:-3] if path.endswith(".gz") else path))[0]
        base = os.path.join(args.out, name)
        lab = run_label(summ)
        cyc = d["k"] / 62.0                     # 共通時計：呼吸の周
        xl = "呼吸の周（節点 / 62）"

        # ---------- F1 ロック ----------
        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        draw_series(ax, cyc, d["S1_abs"], SERIES[0], "|S₁|（双極子）")
        draw_series(ax, cyc, d["S2_abs"], SERIES[1], "|S₂|（四重極）")
        ax.set_xlabel(xl)
        ax.set_ylabel("二次の関係の大きさ")
        ax.set_title(f"F1 ロックの指標：放射の源 |S₁|、|S₂|（0 に収束すればロック）— {lab}", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        written += finish(fig, base + "_F1_lock")

        # ---------- F2 エネルギー ----------
        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        de = d["ebar"] - d["ebar"][0]
        dk = np.diff(d["k"])
        step = float(np.median(dk)) if len(dk) else 1.0
        dtheta = d["dth"] * step                     # 記録間隔の座標時間幅（記録が 1 周ごとなら 62 節点ぶん）
        prad = np.concatenate([[0.0], np.cumsum(0.5 * (d["eloss_em"][1:] + d["eloss_em"][:-1]) * dtheta[1:])])
        draw_series(ax, cyc, de, SERIES[0], "Δε（分布平均の力学的エネルギーの変化）")
        draw_series(ax, cyc, prad, SERIES[1], "∫P_rad dθ（Larmor、差分 D″）")
        ax.set_xlabel(xl)
        ax.set_ylabel("エネルギー（μ で割った値）")
        note = "" if step <= 1 else f"（P_rad は {int(step)} 節点ごとの標本の台形積分。拍動の折り返しを含む）"
        ax.set_title(f"F2 エネルギー収支：Δε と放射損失の積分 {note}— {lab}", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        written += finish(fig, base + "_F2_energy")

        # ---------- F3 倍音分布 ----------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [2, 1]})
        xc, ym = draw_series(ax1, cyc, d["m2"] - 1.0, SERIES[0])
        ax1.set_xlabel(xl)
        ax1.set_ylabel("〈m²〉 − 1")
        ax1.set_title("F3a 倍音分布の二次モーメント 〈m²〉 − 1（0 で単一倍音）")
        w = summ.get("psi_w_final")
        if w:
            M = (len(w) - 1) // 2
            ms = np.arange(-M, M + 1)
            sel = np.abs(ms) <= 6
            ax2.bar(ms[sel], np.array(w)[sel], width=0.6, color=SEQ_BLUE, linewidth=0)
            ax2.set_yscale("log")
            ax2.set_ylim(1e-12, 2)
            ax2.set_xlabel("倍音 m")
            ax2.set_ylabel("|ψ_m|²（対数）")
            ax2.set_title("F3b 最終の倍音分布（|m| ≦ 6）")
            ax2.grid(True, axis="y")
        fig.suptitle(lab, color=INK2, fontsize=9, y=0.995)
        written += finish(fig, base + "_F3_harmonics")

        # ---------- F4 呼吸 ----------
        fig, ax = plt.subplots(figsize=(9, 4.2))
        xc, ym = draw_series(ax, cyc, d["q"], SERIES[0])
        ax.set_xlabel(xl)
        ax.set_ylabel("q = A ζ_c")
        ax.set_title(f"F4 呼吸の振幅比 q = A/ξ̄_c（放射で減れば円軌道化）— {lab}")
        written += finish(fig, base + "_F4_breathing")

        # ---------- F5 時計の比 ----------
        fig, axs = plt.subplots(1, 3, figsize=(12, 4.0))
        draw_series(axs[0], cyc, d["gamma"] - 1.0, SERIES[0])
        axs[0].set_title("F5a γ_k − 1 = Δθ_Y/Δθ_X − 1（近点での値）")
        axs[0].set_ylabel("γ_k − 1")
        draw_series(axs[1], cyc, d["v_over_c"] / ALPHA_FS - 1.0, SERIES[0])
        axs[1].set_title("F5b (v/c)/α − 1（近点での値）")
        axs[1].set_ylabel("(v/c)/α − 1")
        ratio = np.where(d["dthY"] != 0, d["dtha"] / np.where(d["dthY"] != 0, d["dthY"], 1.0), np.nan)
        draw_series(axs[2], cyc, ratio, SERIES[0])
        axs[2].set_title("F5c dθ_a/dθ_Y（a のサイクル／周）")
        axs[2].set_ylabel("dθ_a/dθ_Y")
        for a in axs:
            a.set_xlabel(xl)
        fig.suptitle(f"F5 時計の比（節点時計での読み出し）— {lab}", color=INK2, fontsize=9, y=0.995)
        written += finish(fig, base + "_F5_clocks")

        # ---------- F6 位相 ----------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.0))
        wrapped = np.angle(np.exp(1j * d["argB_minus_thX"]))      # (−π, π] に折り返す
        draw_series(ax1, cyc, wrapped, SERIES[0])
        ax1.set_xlabel(xl)
        ax1.set_ylabel("arg B − θ_X [rad]")
        ax1.set_title("F6a 呼吸位相のずれ arg B − θ_X（読み出し）")
        dph = np.diff(d["cum_phase1"]) / math.pi
        dorb = np.diff(d["cum_thY_over_2pi"])
        ok = dorb > 0
        r = np.full(len(dph), np.nan)
        r[ok] = dph[ok] / dorb[ok]
        draw_series(ax2, cyc[1:], r, SERIES[0])
        ax2.set_xlabel(xl)
        ax2.set_ylabel("ψ₁ の位相の進み／周 [π]")
        ax2.set_title("F6b ψ₁ の位相の進み／周（単一倍音の Bohr 値は 1）")
        fig.suptitle(lab, color=INK2, fontsize=9, y=0.995)
        written += finish(fig, base + "_F6_phase")

        # ---------- F7 倍音の重みと双極子の成分 ----------
        # 台 {0,1,2} のとき |ψ₂|² = (〈m²〉−〈m〉)/2、|ψ₁|² = 2〈m〉−〈m²〉、|ψ₀|² = 1 − |ψ₁|² − |ψ₂|²（|m| ≧ 3 と m = −1 が無視できる間だけ正しい）
        w2 = (d["m2"] - d["m1"]) / 2.0
        w1 = 2.0 * d["m1"] - d["m2"]
        w0 = d["norm"] - w1 - w2
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
        draw_series(ax1, cyc, w0, SERIES[0], "|ψ₀|²")
        draw_series(ax1, cyc, w1, SERIES[1], "|ψ₁|²")
        draw_series(ax1, cyc, w2, SERIES[2], "|ψ₂|²")
        ax1.set_yscale("log")
        ax1.set_xlabel(xl)
        ax1.set_ylabel("|ψ_m|²（対数）")
        ax1.set_title("F7a 倍音の重み（〈m〉・〈m²〉から復元。台 {0,1,2} の間だけ正しい）", loc="left")
        ax1.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
        a_ = np.sqrt(np.clip(w0 * w1, 0, None))
        b_ = np.sqrt(np.clip(w1 * w2, 0, None))
        draw_series(ax2, cyc, a_, SERIES[0], "|ψ₀||ψ₁|（0↔1 の成分）")
        draw_series(ax2, cyc, b_, SERIES[1], "|ψ₁||ψ₂|（1↔2 の成分）")
        ax2.set_xlabel(xl)
        ax2.set_ylabel("双極子 S₁ の成分の大きさ")
        ax2.set_title("F7b S₁ = ψ₀ψ₁* + ψ₁ψ₂* の二つの成分（拍動の最大 = 和、最小 = 差）", loc="left")
        ax2.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.suptitle(lab, color=INK2, fontsize=9, y=0.995)
        written += finish(fig, base + "_F7_weights")
    for p in written:
        print(p)


if __name__ == "__main__":
    main()
