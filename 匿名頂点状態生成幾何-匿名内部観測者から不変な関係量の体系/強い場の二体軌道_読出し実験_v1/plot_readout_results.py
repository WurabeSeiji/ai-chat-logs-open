#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_readout_results.py  (v1.0.0)
readout_binary_pn.py の結果（results/）を図にする。

出力（既定: figures/）
  readout_summary_ja / _en : 読み出した量の相対誤差（正解との差）と、読出しが見積もった 1σ
  readout_details_ja / _en : データセットごとの読出しの中身
                              (a) 近点の方向の進み（近点移動）  (b) 動径周期の減少
                              (c) 相対距離 r(t) の相対誤差        (d) 2天体の固有時間の差 τ_A − τ_B
  それぞれ SVG と PNG（--formats で選択）。日本語フォントが無い環境では英語版だけを書く。

使い方
  python plot_readout_results.py [--results results] [--out figures] [--formats svg,png]
  入力は results/summary.csv, features.csv, reconstruction.csv の3つ。
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
from matplotlib import font_manager          # noqa: E402

HERE = Path(__file__).resolve().parent
JA_FONT_CANDIDATES = ["Yu Gothic", "YuGothic", "Meiryo", "BIZ UDGothic", "MS Gothic", "Noto Sans CJK JP",
                      "Noto Sans JP", "IPAexGothic", "IPAGothic", "Hiragino Sans", "Hiragino Kaku Gothic ProN",
                      "TakaoGothic", "VL Gothic", "Noto Sans CJK TC", "Noto Sans CJK KR"]
COLORS = ["#1f77b4", "#d62728", "#2ca02c"]
MARKERS = ["o", "s", "^"]

CATEGORIES = [  # (quantity, method, ja, en)
    ("q", "obs1 centre of mass", "q\n重心(観測1)", "q\nCoM (obs 1)"),
    ("q", "obs2 GW (from nu)", "q\n重力波(ν から)", "q\nGW (from ν)"),
    ("nu", "obs1 centre of mass", "ν\n重心(観測1)", "ν\nCoM (obs 1)"),
    ("nu", "obs1 orbit decay", "ν\n周期の減少(観測1)", "ν\nperiod decay (obs 1)"),
    ("nu", "obs2 GW", "ν\n重力波", "ν\nGW"),
    ("a", "obs1", "a\n観測1", "a\nobs 1"), ("a", "obs2 GW", "a\n重力波", "a\nGW"),
    ("b", "obs1", "b\n観測1", "b\nobs 1"), ("b", "obs2 GW", "b\n重力波", "b\nGW"),
    ("T_M", "obs1", "T_M\n観測1", "T_M\nobs 1"), ("T_M", "obs2 GW", "T_M\n重力波", "T_M\nGW"),
    ("D", "obs1", "D\n観測1", "D\nobs 1"), ("D", "obs2 GW", "D\n重力波", "D\nGW"),
    ("tau_A - tau_B end", "obs1", "τ_A−τ_B\n観測1", "τ_A−τ_B\nobs 1"),
    ("tau_A - tau_B end", "obs2 GW", "τ_A−τ_B\n重力波", "τ_A−τ_B\nGW"),
    ("radar d_A - d_B mean", "obs1", "d_A−d_B\n観測1", "d_A−d_B\nobs 1"),
    ("radar d_A - d_B mean", "obs2 GW", "d_A−d_B\n重力波", "d_A−d_B\nGW"),
]

TEXT = {
    "ja": {
        "sum_title": "数値実験02 読出し：観測1・観測2から読み出した量の相対誤差（正解と照合）",
        "sum_y": "|読出し − 正解| / |正解|",
        "hollow": "白抜きの印: 正解が 0（等質量）の差なので、τ_A または d_A で割った誤差",
        "sigma": "読出しが見積もった 1σ", "res": "1回の観測の分解能 (1/1000)",
        "case": "{d}（正解 a={a:g}, b={b:g}, q={q:g}）",
        "det_title": "数値実験02 読出しの中身（点: 観測から取り出した特徴量、線: 当てはめたモデル）",
        "prec": "(a) 近点の方向の進み", "prec_y": "累積の近点移動 [rad]", "n": "近点通過の番号 n",
        "per": "(b) 動径周期の減少", "per_y": "動径周期 [s]",
        "rerr": "(c) 相対距離 r(t) の相対誤差", "rerr_y": "(r読出し − r正解) / r正解", "tM": "t [M]（正解）",
        "dtau": "(d) 固有時間の差 τ_A − τ_B", "dtau_y": "τ_A − τ_B [M]",
        "o1": "観測1", "o2": "重力波（観測2）", "fit": "モデル", "truth": "正解",
        "o1d": "観測1の角度 ÷ 読み出した尺度", "o1m": "観測1の読出しで再構成", "o2m": "重力波の読出しで再構成",
    },
    "en": {
        "sum_title": "Experiment 02 readout: relative errors of the quantities read from observations 1 and 2",
        "sum_y": "|readout − truth| / |truth|",
        "hollow": "open symbols: the true difference is 0 (equal masses); error divided by τ_A or d_A",
        "sigma": "1σ estimated by the readout", "res": "single-observation resolution (1/1000)",
        "case": "{d} (truth a={a:g}, b={b:g}, q={q:g})",
        "det_title": "Experiment 02 readout details (points: features from the data, lines: fitted model)",
        "prec": "(a) Advance of the periapsis direction", "prec_y": "accumulated advance [rad]",
        "n": "periapsis passage n",
        "per": "(b) Decay of the radial period", "per_y": "radial period [s]",
        "rerr": "(c) Relative error of r(t)", "rerr_y": "(r_read − r_true) / r_true", "tM": "t [M] (truth)",
        "dtau": "(d) Proper-time difference τ_A − τ_B", "dtau_y": "τ_A − τ_B [M]",
        "o1": "obs 1", "o2": "GW (obs 2)", "fit": "model", "truth": "truth",
        "o1d": "obs-1 angle ÷ read scale", "o1m": "reconstructed (obs-1 readout)", "o2m": "reconstructed (GW readout)",
    },
}


def find_ja_font():
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in JA_FONT_CANDIDATES:
        if name in available:
            return name
    return None


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0].keys():
        vals = []
        for r in rows:
            try:
                vals.append(float(r[k]) if r[k] != "" else math.nan)
            except ValueError:
                vals.append(r[k])
        out[k] = np.array(vals, dtype=object if isinstance(vals[0], str) else float)
    return out


def load_evals(res_dir):
    """results/summary.csv（読出しと正解の対照表）をデータセットごとにまとめる。"""
    evals = {}
    with open(res_dir / "summary.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            row = {"quantity": r["quantity"], "method": r["method"], "normalized_by": r["normalized_by"]}
            for k in ("readout", "sigma", "truth", "relative_error"):
                row[k] = float(r[k]) if r[k] != "" else None
            evals.setdefault(r["dataset"], []).append(row)
    return [{"dataset": d, "rows": rows} for d, rows in evals.items()]


def case_labels(evals, T):
    labels = {}
    for ds in evals:
        truth = {(r["quantity"], r["method"]): r["truth"] for r in ds["rows"]}
        labels[ds["dataset"]] = T["case"].format(d=ds["dataset"], a=truth[("a", "obs1")], b=truth[("b", "obs1")],
                                                 q=truth[("q", "obs1 centre of mass")])
    return labels


def plot_summary(evals, T, lang):
    fig, ax = plt.subplots(figsize=(15, 6.2))
    labels = case_labels(evals, T)
    x = np.arange(len(CATEGORIES))
    for k, ds in enumerate(evals):
        rows = {(r["quantity"], r["method"]): r for r in ds["rows"]}
        err, sig, hol = [], [], []
        for q, m, _, _ in CATEGORIES:
            r = rows[(q, m)]
            err.append(abs(r["relative_error"]))
            hol.append(r.get("normalized_by") == "reference")
            sig.append(r["sigma"] / abs(r["truth"]) if r.get("sigma") else math.nan)
        xs = x + (k - 1) * 0.22
        e_arr, h_arr = np.maximum(err, 1e-7), np.array(hol)
        ax.scatter(xs[~h_arr], e_arr[~h_arr], color=COLORS[k], marker=MARKERS[k], s=36, zorder=3,
                   label=labels[ds["dataset"]])
        if h_arr.any():
            ax.scatter(xs[h_arr], e_arr[h_arr], facecolors="none", edgecolors=COLORS[k], marker=MARKERS[k],
                       s=36, zorder=3)
        ax.scatter(xs, sig, color=COLORS[k], marker="_", s=160, linewidths=1.6, zorder=2,
                   label=T["sigma"] if k == 0 else None)
    ax.axhline(1e-3, color="0.4", ls="--", lw=1.0, label=T["res"])
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 1e-1)
    ax.set_xticks(x)
    ax.set_xticklabels([c[2] if lang == "ja" else c[3] for c in CATEGORIES], fontsize=8)
    for xv in np.arange(len(CATEGORIES) - 1) + 0.5:
        ax.axvline(xv, color="0.9", lw=0.6, zorder=0)
    ax.set_ylabel(T["sum_y"])
    ax.grid(axis="y", which="major", color="0.85", lw=0.6)
    ax.legend(fontsize=8, loc="upper left", ncol=2)
    ax.set_title(T["sum_title"], fontsize=12)
    ax.text(0.995, 0.02, T["hollow"], transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="0.3")
    fig.tight_layout()
    return fig


def plot_details(res_dir, evals, T):
    names = [ds["dataset"] for ds in evals]
    feat, rec = read_csv(res_dir / "features.csv"), read_csv(res_dir / "reconstruction.csv")
    labels = case_labels(evals, T)
    fig, axs = plt.subplots(len(names), 4, figsize=(17, 3.6 * len(names)), squeeze=False)
    for i, name in enumerate(names):
        f = {k: v[feat["dataset"] == name] for k, v in feat.items()}
        rc = {k: v[rec["dataset"] == name] for k, v in rec.items()}
        for chan, color, mk, lab, turn in (("obs1", COLORS[0], "o", T["o1"], 1.0), ("obs2", COLORS[1], "x", T["o2"], 2.0)):
            s = f["channel"] == chan
            n = f["n"][s].astype(float)
            ref = f["angle_fit"][s][0] - 2.0 * math.pi * turn * n[0]
            y_obs = (f["angle_rad"][s] - ref) / turn - 2.0 * math.pi * n
            y_fit = (f["angle_fit"][s] - ref) / turn - 2.0 * math.pi * n
            axs[i, 0].plot(n, y_obs, mk, color=color, ms=4, mfc="none", label=lab)
            axs[i, 0].plot(n, y_fit, "-", color=color, lw=0.9, alpha=0.7, label=T["fit"] if chan == "obs1" else None)
            p_obs, p_fit = np.diff(f["t_s"][s]), np.diff(f["t_fit"][s])
            axs[i, 1].plot(n[:-1], p_obs, mk, color=color, ms=4, mfc="none", label=lab)
            axs[i, 1].plot(n[:-1], p_fit, "-", color=color, lw=0.9, alpha=0.7, label=T["fit"] if chan == "obs1" else None)
        t = rc["t_M_true"]
        axs[i, 2].plot(t, (rc["r_obs1_data"] - rc["r_true"]) / rc["r_true"], ".", color="0.6", ms=2.5, label=T["o1d"])
        axs[i, 2].plot(t, (rc["r_obs1_model"] - rc["r_true"]) / rc["r_true"], "-", color=COLORS[0], lw=1.0, label=T["o1m"])
        axs[i, 2].plot(t, (rc["r_obs2_model"] - rc["r_true"]) / rc["r_true"], "-", color=COLORS[1], lw=1.0, label=T["o2m"])
        axs[i, 3].plot(t, rc["dtau_true"], "-", color="k", lw=2.2, alpha=0.35, label=T["truth"])
        axs[i, 3].plot(t, rc["dtau_obs1"], "--", color=COLORS[0], lw=1.1, label=T["o1m"])
        axs[i, 3].plot(t, rc["dtau_obs2"], ":", color=COLORS[1], lw=1.4, label=T["o2m"])
        for j, (title, xl, yl) in enumerate(((T["prec"], T["n"], T["prec_y"]), (T["per"], T["n"], T["per_y"]),
                                             (T["rerr"], T["tM"], T["rerr_y"]), (T["dtau"], T["tM"], T["dtau_y"]))):
            ax = axs[i, j]
            ax.set_title(f"{labels[name]}\n{title}" if j == 0 else title, fontsize=9)
            ax.set_xlabel(xl, fontsize=8)
            ax.set_ylabel(yl, fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(color="0.9", lw=0.5)
            ax.legend(fontsize=6.5, loc="best")
    fig.suptitle(T["det_title"], fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=str(HERE / "results"))
    ap.add_argument("--out", default=str(HERE / "figures"))
    ap.add_argument("--formats", default="svg,png")
    args = ap.parse_args()
    res_dir, out = Path(args.results), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    evals = load_evals(res_dir)
    ja_font = find_ja_font()
    langs = ["ja", "en"] if ja_font else ["en"]
    if not ja_font:
        print("日本語フォントが見つからないため、英語版だけを書きます。")
    for lang in langs:
        rc = {"svg.fonttype": "none", "axes.unicode_minus": False}
        if lang == "ja":
            rc["font.family"] = ja_font
        with plt.rc_context(rc):
            for stem, fig in ((f"readout_summary_{lang}", plot_summary(evals, TEXT[lang], lang)),
                              (f"readout_details_{lang}", plot_details(res_dir, evals, TEXT[lang]))):
                for fmt in args.formats.split(","):
                    fig.savefig(out / f"{stem}.{fmt}", dpi=150)
                plt.close(fig)
                print("書き出し:", stem)


if __name__ == "__main__":
    main()
