#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数値実験01 の生成データを図化する。

出力（既定: このスクリプトと同じフォルダの figures/）
    overview_input_ja / _en   解く側から見える入力（形と並び順だけ。データセットID以外は載せない）
    overview_truth_ja / _en   正解の重ね描き（引力中心の焦点・近点方向・進行方向・時刻）
    それぞれ PNG と SVG。日本語フォントが見つからない場合は英語版のみ出力する。

使い方
    python plot_conic_orbit_samples.py
    python plot_conic_orbit_samples.py --data データフォルダ --out 図の出力先
"""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                 # noqa: E402
from matplotlib import font_manager             # noqa: E402
from matplotlib.lines import Line2D             # noqa: E402

JA_FONT_CANDIDATES = [
    "Yu Gothic", "YuGothic", "Meiryo", "BIZ UDGothic", "MS Gothic",
    "Noto Sans CJK JP", "Noto Sans JP", "IPAexGothic", "IPAGothic",
    "Hiragino Sans", "Hiragino Kaku Gothic ProN", "TakaoGothic", "VL Gothic",
]

TEXT = {
    "en": {
        "input_title": "Experiment 01 - solver input: ordered point sets (shape and order only, no time stamps)",
        "truth_title": "Experiment 01 - ground truth: attracting focus, periapsis direction, motion direction and time",
        "pattern": {"circle": "circle", "ellipse": "ellipse", "parabola": "parabola", "hyperbola": "hyperbola"},
        "mode": {"A": "A: jitter", "B": "B: clustered"},
        "first": "first input point",
        "curve": "generating curve",
        "focus": "attracting focus",
        "empty": "empty focus",
        "peri": "periapsis direction",
        "motion": "motion direction",
        "time": "time along the motion (normalized, 0 to 1)",
    },
    "ja": {
        "input_title": "数値実験01　解く側に渡す入力：順序付き点列（形と並び順だけ・時刻なし）",
        "truth_title": "数値実験01　正解：引力中心の焦点・近点方向・進行方向・時刻",
        "pattern": {"circle": "円", "ellipse": "楕円", "parabola": "放物線", "hyperbola": "双曲線"},
        "mode": {"A": "A：一様ゆらぎ", "B": "B：疎密変調"},
        "first": "入力の先頭点",
        "curve": "生成曲線",
        "focus": "引力中心（焦点）",
        "empty": "もう一方の焦点",
        "peri": "近点方向",
        "motion": "進行方向",
        "time": "進行方向に沿った時刻（規格化、0〜1）",
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


def draw_input_panel(ax, did, x, y):
    ax.plot(x, y, "-", lw=0.4, color="0.65", zorder=1)
    ax.scatter(x, y, s=2.0, color="C0", linewidths=0, zorder=2)
    ax.plot(x[0], y[0], marker="^", color="C3", ms=5, linestyle="none", zorder=3)
    ax.set_title(did, fontsize=10)


def draw_truth_panel(ax, did, meta, tru, x, y, T):
    g = meta["derived_geometry"]
    p = meta["shape_parameters"]["p"]
    e = meta["shape_parameters"]["e"]
    O = meta["pose"]["focus_offset_O_xy"]
    phi = meta["pose"]["rotation_phi_rad"]
    closed = meta["closed_curve"]

    th0, th1 = g["theta_range_rad"]
    th = np.linspace(th0, th1, 3000)
    r = p / (1.0 + e * np.cos(th))
    xc, yc = r * np.cos(th), r * np.sin(th)
    c, s = math.cos(phi), math.sin(phi)
    ax.plot(O[0] + c * xc - s * yc, O[1] + s * xc + c * yc, "-", lw=0.8, color="0.8", zorder=1)

    te = tru["tau_elapsed"]
    total = g["period_tau"] if closed else te.max()
    sc = ax.scatter(x, y, c=te / total, cmap="viridis", vmin=0.0, vmax=1.0, s=3.0, linewidths=0, zorder=2)

    ax.plot(O[0], O[1], marker="*", color="C3", ms=10, linestyle="none", zorder=4)
    if g["pattern"] == "ellipse" and closed:
        ef = g["empty_focus_xy"]
        ax.plot(ef[0], ef[1], marker="o", mfc="none", mec="C3", ms=6, linestyle="none", zorder=4)
    if g["pattern"] != "circle":
        ax.annotate("", xy=g["periapsis_point_xy"], xytext=O,
                    arrowprops=dict(arrowstyle="->", color="C3", lw=1.0), zorder=5)

    order = np.argsort(tru["time_rank"])
    i0 = order[0]
    i1 = order[min(len(order) - 1, max(3, len(order) // 40))]
    ax.annotate("", xy=(x[i1], y[i1]), xytext=(x[i0], y[i0]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2), zorder=6)

    ax.set_title(f'{did}  {T["pattern"][g["pattern"]]}  e={e:.2f}  {T["mode"][meta["sampling"]["mode"]]}',
                 fontsize=9)
    return sc


def plot_overview(items, lang, truth):
    T = TEXT[lang]
    n = len(items)
    ncol = 4
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 4.0 * nrow + 1.2))
    axes = np.atleast_1d(axes).ravel()
    sc = None
    for ax, (d, inp, meta, tru) in zip(axes, items):
        x, y = inp["x"], inp["y"]
        if truth:
            sc = draw_truth_panel(ax, d["dataset_id"], meta, tru, x, y, T)
        else:
            draw_input_panel(ax, d["dataset_id"], x, y)
        ax.set_aspect("equal", adjustable="datalim")
        ax.tick_params(labelsize=6)
    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle(T["truth_title"] if truth else T["input_title"], fontsize=13)
    if truth:
        handles = [
            Line2D([], [], color="0.8", lw=1.2, label=T["curve"]),
            Line2D([], [], marker="*", color="C3", linestyle="none", ms=10, label=T["focus"]),
            Line2D([], [], marker="o", mfc="none", mec="C3", linestyle="none", ms=6, label=T["empty"]),
            Line2D([], [], color="C3", lw=1.0, label=T["peri"]),
            Line2D([], [], color="k", lw=1.2, label=T["motion"]),
        ]
    else:
        handles = [Line2D([], [], marker="^", color="C3", linestyle="none", ms=6, label=T["first"])]
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 0.97))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    if truth and sc is not None:
        cax = fig.add_axes([0.3, 0.045, 0.4, 0.008])
        cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cb.set_label(T["time"], fontsize=8)
        cb.ax.tick_params(labelsize=7)
    return fig


def save(fig, out, stem):
    fig.savefig(out / f"{stem}.png", dpi=150)
    fig.savefig(out / f"{stem}.svg")
    plt.close(fig)
    print(f"wrote {stem}.png / {stem}.svg")


def main():
    parser = argparse.ArgumentParser(description="数値実験01 の生成データを図化する")
    parser.add_argument("--data", default=None, help="データフォルダ（既定: このスクリプトと同じフォルダの data）")
    parser.add_argument("--out", default=None, help="図の出力先（既定: このスクリプトと同じフォルダの figures）")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    data = Path(args.data) if args.data else here / "data"
    out = Path(args.out) if args.out else here / "figures"
    out.mkdir(parents=True, exist_ok=True)

    manifest = load_json(data / "truth" / "manifest.json")
    items = []
    for d in manifest["datasets"]:
        items.append((d,
                      read_csv(data / d["input_file"]),
                      load_json(data / d["meta_file"]),
                      read_csv(data / d["truth_file"])))

    langs = ["en"]
    ja_font = find_ja_font()
    if ja_font:
        langs.insert(0, "ja")
    else:
        print("Japanese font not found: only English figures are written.")

    for lang in langs:
        rc = {"font.family": ja_font, "axes.unicode_minus": False} if lang == "ja" else {}
        with plt.rc_context(rc):
            save(plot_overview(items, lang, truth=False), out, f"overview_input_{lang}")
            save(plot_overview(items, lang, truth=True), out, f"overview_truth_{lang}")


if __name__ == "__main__":
    main()
