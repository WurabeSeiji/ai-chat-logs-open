#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数値実験02 の生成データを、そのまま図にする。

出力（既定: このスクリプトと同じフォルダの figures/）
    obs_Dk_ja / _en    観測データそのもの（観測1の天球上の位置、観測2の重力波の歪み）。
                       解析側から見える情報だけを載せる（データセットID以外の条件は書かない）。
    truth_Dk_ja / _en  正解の時系列（M 単位の相対軌道、相対距離、2体の固有時間の差）。
    それぞれ PNG と SVG。日本語フォントが見つからない場合は英語版のみ出力する。

使い方
    python plot_binary_pn_observations.py
    python plot_binary_pn_observations.py --data データフォルダ --out 図の出力先
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                 # noqa: E402
from matplotlib import font_manager             # noqa: E402

JA_FONT_CANDIDATES = [
    "Yu Gothic", "YuGothic", "Meiryo", "BIZ UDGothic", "MS Gothic",
    "Noto Sans CJK JP", "Noto Sans JP", "IPAexGothic", "IPAGothic",
    "Hiragino Sans", "Hiragino Kaku Gothic ProN", "TakaoGothic", "VL Gothic",
]

TEXT = {
    "en": {
        "obs_title": "Experiment 02 observation data {did} (as observed by c)",
        "sky": "(a) Observation 1: sky positions",
        "sky_x": "x [nrad]", "sky_y": "y [nrad]",
        "body1": "body 1", "body2": "body 2",
        "full": "(b) Observation 2: strain h+ (whole record)",
        "start": "(c) Observation 2: beginning",
        "end": "(d) Observation 2: end",
        "t_s": "observer time t [s]", "strain": "strain",
        "truth_title": "Experiment 02 ground truth {did}: a = {a:g}, b = {b:g}, m_A/m_B = {q:g} (nu = {nu:.4f})",
        "orbit": "(a) Relative orbit (units of M)",
        "sep": "(b) Separation r (units of M)",
        "tau": "(c) Proper-time difference tau_A - tau_B",
        "t_M": "t [M]", "x_M": "x [M]", "y_M": "y [M]", "r_M": "r [M]", "dtau": "tau_A - tau_B [M]",
        "time": "time (normalized)",
    },
    "ja": {
        "obs_title": "数値実験02　観測データ {did}（観測者 c が得るもの）",
        "sky": "(a) 観測1：天球上の位置",
        "sky_x": "x [nrad]", "sky_y": "y [nrad]",
        "body1": "天体1", "body2": "天体2",
        "full": "(b) 観測2：重力波の歪み h+（全体）",
        "start": "(c) 観測2：はじめの部分",
        "end": "(d) 観測2：終わりの部分",
        "t_s": "観測者の時刻 t [s]", "strain": "歪み",
        "truth_title": "数値実験02　正解 {did}：a = {a:g}, b = {b:g}, m_A/m_B = {q:g}（ν = {nu:.4f}）",
        "orbit": "(a) 相対軌道（M 単位）",
        "sep": "(b) 相対距離 r（M 単位）",
        "tau": "(c) 固有時間の差 τ_A − τ_B",
        "t_M": "t [M]", "x_M": "x [M]", "y_M": "y [M]", "r_M": "r [M]", "dtau": "τ_A − τ_B [M]",
        "time": "時刻（規格化）",
    },
}


def find_ja_font():
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in JA_FONT_CANDIDATES:
        if name in available:
            return name
    return None


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f))
    header, body = rows[0], rows[1:]
    return {h: np.array([float(r[i]) for r in body]) for i, h in enumerate(header)}


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def plot_obs(did, o1, o2, T):
    fig, axs = plt.subplots(2, 2, figsize=(13.0, 9.5))
    ax = axs[0, 0]
    ax.scatter(o1["body1_x_nrad"], o1["body1_y_nrad"], s=1.5, color="C0", linewidths=0, label=T["body1"])
    ax.scatter(o1["body2_x_nrad"], o1["body2_y_nrad"], s=1.5, color="C1", linewidths=0, label=T["body2"])
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(T["sky_x"])
    ax.set_ylabel(T["sky_y"])
    ax.set_title(T["sky"], fontsize=10)
    ax.legend(fontsize=8, markerscale=6, loc="upper right")

    t = o2["t_s"]
    ax = axs[0, 1]
    ax.plot(t, o2["h_plus"], lw=0.3, color="C2")
    ax.set_xlabel(T["t_s"])
    ax.set_ylabel(T["strain"])
    ax.set_title(T["full"], fontsize=10)

    span = t[-1] - t[0]
    win = span * 0.04
    for axw, (lo, hi), title in ((axs[1, 0], (t[0], t[0] + win), T["start"]),
                                 (axs[1, 1], (t[-1] - win, t[-1]), T["end"])):
        m = (t >= lo) & (t <= hi)
        axw.plot(t[m], o2["h_plus"][m], lw=0.8, color="C2", label="h+")
        axw.plot(t[m], o2["h_cross"][m], lw=0.8, color="C3", label="h×")
        axw.set_xlabel(T["t_s"])
        axw.set_ylabel(T["strain"])
        axw.set_title(title, fontsize=10)
        axw.legend(fontsize=8, loc="upper right")
    for a_ in axs.ravel():
        a_.tick_params(labelsize=7)
    fig.suptitle(T["obs_title"].format(did=did), fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def plot_truth(did, meta, ts, T):
    sp = meta["orbit_shape_parameters"]
    ms = meta["masses"]
    fig, axs = plt.subplots(1, 3, figsize=(16.0, 5.2))
    tn = ts["t_M"] / ts["t_M"][-1]
    ax = axs[0]
    sc = ax.scatter(ts["x_rel"], ts["y_rel"], c=tn, s=0.6, cmap="viridis", linewidths=0)
    ax.plot([0.0], [0.0], marker="+", color="k", ms=8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(T["x_M"])
    ax.set_ylabel(T["y_M"])
    ax.set_title(T["orbit"], fontsize=10)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(T["time"], fontsize=8)
    cb.ax.tick_params(labelsize=7)

    axs[1].plot(ts["t_M"], ts["r"], lw=0.4, color="C0")
    axs[1].set_xlabel(T["t_M"])
    axs[1].set_ylabel(T["r_M"])
    axs[1].set_title(T["sep"], fontsize=10)

    axs[2].plot(ts["t_M"], ts["tau_A"] - ts["tau_B"], lw=0.8, color="C4")
    axs[2].set_xlabel(T["t_M"])
    axs[2].set_ylabel(T["dtau"])
    axs[2].set_title(T["tau"], fontsize=10)
    for a_ in axs:
        a_.tick_params(labelsize=7)
    fig.suptitle(T["truth_title"].format(did=did, a=sp["a"], b=sp["b"], q=ms["q_mA_over_mB"], nu=ms["nu"]),
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def save(fig, out, stem):
    fig.savefig(out / f"{stem}.png", dpi=150)
    fig.savefig(out / f"{stem}.svg")
    plt.close(fig)
    print(f"wrote {stem}.png / {stem}.svg")


def main():
    parser = argparse.ArgumentParser(description="数値実験02 の生成データを図にする")
    parser.add_argument("--data", default=None, help="データフォルダ（既定: このスクリプトと同じフォルダの data）")
    parser.add_argument("--out", default=None, help="図の出力先（既定: このスクリプトと同じフォルダの figures）")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    data = Path(args.data) if args.data else here / "data"
    out = Path(args.out) if args.out else here / "figures"
    out.mkdir(parents=True, exist_ok=True)

    manifest = load_json(data / "truth" / "manifest.json")
    langs = ["en"]
    ja_font = find_ja_font()
    if ja_font:
        langs.insert(0, "ja")
    else:
        print("Japanese font not found: only English figures are written.")

    for d in manifest["datasets"]:
        did = d["dataset_id"]
        o1 = read_csv(data / d["obs1_file"])
        o2 = read_csv(data / d["obs2_file"])
        meta = load_json(data / d["meta_file"])
        ts = read_csv(data / d["timeseries_file"])
        for lang in langs:
            rc = {"font.family": ja_font, "axes.unicode_minus": False} if lang == "ja" else {}
            with plt.rc_context(rc):
                save(plot_obs(did, o1, o2, TEXT[lang]), out, f"obs_{did}_{lang}")
                save(plot_truth(did, meta, ts, TEXT[lang]), out, f"truth_{did}_{lang}")


if __name__ == "__main__":
    main()
