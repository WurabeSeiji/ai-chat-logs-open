#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_readouts_v7.py — v7（トーラス波 Ψ_{n,m}）の走行結果を仕様 §7 の読み出しに沿って図にする。SVG と PNG。
  F1 ロック      |D|（双極子）、|Q|（四重極）
  F2 エネルギー   Δε（分布平均の力学的エネルギー）と ∫P_rad dθ（Larmor、差分 D″）
  F3 倍音の重み  w(0,1)、w(0,2)、w(1,1)、その他（対数）
  F4 呼吸        〈ξ〉（X̂ の期待値。隣り合う n の拍動）
  F5 時計の比    γ_k − 1（拍動の比 ρ̄/ω_X）、(v/c)/α − 1、dθ_a/dθ_Y
  F6 位相        Ψ_{0,1} の位相の進み／呼吸（π 単位）、D の位相
横軸は共通時計：呼吸の拍動の周（節点/62）。表示規約は plot_readouts.py と同じ（dataviz）。
"""
import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from plot_readouts import (setup_fonts, load_run, draw_series, finish, SERIES, SEQ_BLUE, INK2, ALPHA_FS)  # noqa: E402


def run_label(summ):
    ini = summ.get("init", {})
    return f"v7 δ = {ini.get('delta', '?')}, f = {ini.get('f', '?')}, スピン {ini.get('sa', '?')}/{ini.get('sb', '?')}, 節点 {ini.get('nodes', '?')}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="+")
    ap.add_argument("--out", default="results_v7/figures")
    args = ap.parse_args()
    setup_fonts()
    os.makedirs(args.out, exist_ok=True)
    written = []
    for path in args.csv:
        d, summ = load_run(path)
        name = os.path.splitext(os.path.basename(path[:-3] if path.endswith(".gz") else path))[0]
        base = os.path.join(args.out, name)
        lab = run_label(summ)
        cyc = d["k"] / 62.0
        xl = "呼吸の拍動の周（節点 / 62）"

        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        draw_series(ax, cyc, d["D_abs"], SERIES[0], "|D| = |〈r e^{iφ}〉|（双極子）")
        draw_series(ax, cyc, d["Q_abs"] / 300.0, SERIES[1], "|Q|/300（四重極、目盛り合わせ）")
        ax.set_xlabel(xl)
        ax.set_ylabel("大きさ（a_a 単位）")
        ax.set_title(f"F1 ロックの指標：放射の源 |D|、|Q|（0 に収束すればロック）— {lab}", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        written += finish(fig, base + "_F1_lock")

        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        de = d["ebar"] - d["ebar"][0]
        dk = np.diff(d["k"])
        step = float(np.median(dk)) if len(dk) else 1.0
        dtheta = d["dth"] * step
        prad = np.concatenate([[0.0], np.cumsum(0.5 * (d["eloss_em"][1:] + d["eloss_em"][:-1]) * dtheta[1:])])
        draw_series(ax, cyc, de, SERIES[0], "Δε（分布平均の力学的エネルギーの変化）")
        draw_series(ax, cyc, prad, SERIES[1], "∫P_rad dθ（Larmor、差分 D″）")
        ax.set_xlabel(xl)
        ax.set_ylabel("エネルギー（μ で割った値）")
        note = "" if step <= 1 else f"（P_rad は {int(step)} 節点ごとの標本の台形積分）"
        ax.set_title(f"F2 エネルギー収支 {note}— {lab}", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        written += finish(fig, base + "_F2_energy")

        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        draw_series(ax, cyc, d["w01"], SERIES[0], "w(0,1)（円軌道、床）")
        draw_series(ax, cyc, d["w02"], SERIES[1], "w(0,2)")
        draw_series(ax, cyc, d["w11"], SERIES[2], "w(1,1)")
        rest = np.clip(d["w_rest"], 1e-18, None)
        ax.plot(*((cyc, rest) if len(cyc) <= 5000 else (cyc[::max(1, len(cyc) // 2000)], rest[::max(1, len(cyc) // 2000)])),
                color=INK2, linewidth=1.0, label="その他の (n,m)")
        ax.set_yscale("log")
        ax.set_ylim(1e-16, 2)
        ax.set_xlabel(xl)
        ax.set_ylabel("|Ψ_{n,m}|²（対数）")
        ax.set_title(f"F3 倍音の重み — {lab}", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        written += finish(fig, base + "_F3_weights")

        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        draw_series(ax, cyc, d["xi_mean"], SERIES[0])
        ax.set_xlabel(xl)
        ax.set_ylabel("〈ξ〉（a_a 単位）")
        ax.set_title(f"F4 呼吸：〈ξ〉 = 〈Ψ|X̂|Ψ〉（隣り合う n の拍動。帯は窓内の最小–最大）— {lab}", loc="left")
        written += finish(fig, base + "_F4_breathing")

        fig, axs = plt.subplots(1, 3, figsize=(12, 4.0))
        draw_series(axs[0], cyc, d["gamma"] - 1.0, SERIES[0])
        axs[0].set_title("F5a γ_k − 1 = ρ̄/ω_X − 1（拍動の比）")
        axs[0].set_ylabel("γ_k − 1")
        draw_series(axs[1], cyc, d["v_over_c"] / ALPHA_FS - 1.0, SERIES[0])
        axs[1].set_title("F5b (v/c)/α − 1（〈2gm/ξ〉）")
        axs[1].set_ylabel("(v/c)/α − 1")
        ratio = np.where(d["dthY"] != 0, d["dtha"] / np.where(d["dthY"] != 0, d["dthY"], 1.0), np.nan)
        draw_series(axs[2], cyc, ratio, SERIES[0])
        axs[2].set_title("F5c dθ_a/dθ_Y（a のサイクル／回転の拍動 2π）")
        axs[2].set_ylabel("dθ_a/dθ_Y")
        for a in axs:
            a.set_xlabel(xl)
        fig.suptitle(f"F5 時計の比 — {lab}", color=INK2, fontsize=9, y=0.995)
        written += finish(fig, base + "_F5_clocks")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.0))
        dph = np.diff(d["cum_phase01"]) / math.pi
        dcyc = np.diff(cyc)
        ok = dcyc > 0
        r = np.full(len(dph), np.nan)
        r[ok] = dph[ok] / dcyc[ok]
        draw_series(ax1, cyc[1:], r, SERIES[0])
        ax1.set_xlabel(xl)
        ax1.set_ylabel("Ψ(0,1) の位相の進み／呼吸 [π]")
        ax1.set_title("F6a Ψ(0,1) の位相の進み／呼吸（Newton では 8/3）")
        wrapped = np.angle(np.exp(1j * d["D_arg"]))
        draw_series(ax2, cyc, wrapped, SERIES[0])
        ax2.set_xlabel(xl)
        ax2.set_ylabel("arg D [rad]")
        ax2.set_title("F6b 双極子の位相 arg D（回転の拍動）")
        fig.suptitle(lab, color=INK2, fontsize=9, y=0.995)
        written += finish(fig, base + "_F6_phase")
    for p in written:
        print(p)


if __name__ == "__main__":
    main()
