#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
条件A/B の縦積み拡大図。ZOOMS は paper7_f_projection_v1 / paper8_a2a_seedless_v1 と同一仕様。

  python3 make_ablation_figures_v1.py 40

入力: run_ablation_everystep_v1.py が出した ablation_N{n}_cond{A,B}_everystep.csv
      （既存27列が公開CSVと一致検証済み。追加列 f_projection のみが新規）

図は3種:
  figureZ_stacked_zoom_condA_N{n}.png   無シード
  figureZ_stacked_zoom_condB_N{n}.png   δ=1e-15
  figureZ_compare_AB_N{n}.png           両条件を同一軸で重ねる
縦軸は f_projection（射影形）。引き算形は 1e-16 で床に張り付くため参考として点線で併記する。
"""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

ZOOMS = [(5000, 500, "full  0-5000", False),
         (250, 25, "x20  0-250", True),
         (25, 5, "x200  0-25  (per-step)", True)]
LABEL = {"A": "condition A: no seed (Z0=v)",
         "B": "condition B: seeded (delta=1e-15)"}


def load(n, cond):
    p = os.path.join(HERE, f"ablation_N{n:05d}_cond{cond}_everystep.csv")
    r = list(csv.DictReader(open(p, encoding="utf-8")))
    return {"t": np.array([float(x["step"]) for x in r]),
            "fp": np.array([float(x["f_projection"]) for x in r]),
            "fs": np.array([float(x["f_outside_parent"]) for x in r])}


def meta(n, cond):
    p = os.path.join(HERE, f"ablation_N{n:05d}_cond{cond}_everystep_meta.json")
    return json.load(open(p, encoding="utf-8"))


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{stem}.{ext}"), dpi=150 if ext == "png" else None,
                    bbox_inches="tight")
    plt.close(fig)


def stacked(n, cond):
    d = load(n, cond); m = meta(n, cond)
    fig, axes = plt.subplots(len(ZOOMS), 1, figsize=(11, 4.0 * len(ZOOMS)), squeeze=False)
    for ax, (xm, st, tag, fit) in zip(axes[:, 0], ZOOMS):
        k = d["t"] <= xm
        ax.semilogy(d["t"][k], d["fp"][k], "k-", lw=1.2,
                    marker=("o" if xm <= 25 else None), ms=3, label="f (projection)")
        ax.semilogy(d["t"][k], np.abs(d["fs"][k]), color="#b0b0b0", ls=":", lw=1.0,
                    label="|f| (subtraction, recorded)")
        if m["crossing"] is not None and m["crossing"] <= xm:
            ax.axvline(m["crossing"], color="r", ls=":", lw=0.9)
        ax.set_xlim(0, xm); ax.set_xticks(np.arange(0, xm + 1, st))
        ax.grid(alpha=0.25, which="both"); ax.set_ylabel(f"f\n{tag}")
        if fit:
            ax.set_ylim(d["fp"][k].min() / 3.0, d["fp"][k].max() * 3.0)
    axes[0, 0].legend(fontsize=8, loc="lower right")
    axes[-1, 0].set_xlabel("absolute step")
    fig.suptitle(f"N={n}  {LABEL[cond]}:  f stacked zooms  "
                 f"(f(0)={m['f_projection_first']:.6e}, crossing={m['crossing']})")
    save(fig, f"figureZ_stacked_zoom_cond{cond}_N{n:05d}")


def compare(n):
    a, b = load(n, "A"), load(n, "B")
    ma, mb = meta(n, "A"), meta(n, "B")
    fig, axes = plt.subplots(len(ZOOMS), 1, figsize=(11, 4.0 * len(ZOOMS)), squeeze=False)
    for ax, (xm, st, tag, fit) in zip(axes[:, 0], ZOOMS):
        ka, kb = a["t"] <= xm, b["t"] <= xm
        ax.semilogy(a["t"][ka], a["fp"][ka], "-", color="#d95f02", lw=1.3,
                    marker=("o" if xm <= 25 else None), ms=3, label="A: no seed")
        ax.semilogy(b["t"][kb], b["fp"][kb], "-", color="#1b6ca8", lw=1.3,
                    marker=("s" if xm <= 25 else None), ms=3, label="B: delta=1e-15")
        if ma["crossing"] is not None and ma["crossing"] <= xm:
            ax.axvline(ma["crossing"], color="r", ls=":", lw=0.9)
        ax.set_xlim(0, xm); ax.set_xticks(np.arange(0, xm + 1, st))
        ax.grid(alpha=0.25, which="both"); ax.set_ylabel(f"f (projection)\n{tag}")
        if fit:
            lo = min(a["fp"][ka].min(), b["fp"][kb].min())
            hi = max(a["fp"][ka].max(), b["fp"][kb].max())
            ax.set_ylim(lo / 3.0, hi * 3.0)
    axes[0, 0].legend(fontsize=9, loc="lower right")
    axes[-1, 0].set_xlabel("absolute step")
    fig.suptitle(f"N={n}  seeded vs seedless:  f(0)={ma['f_projection_first']:.4e} (A) / "
                 f"{mb['f_projection_first']:.4e} (B),  crossing={ma['crossing']} / {mb['crossing']}")
    save(fig, f"figureZ_compare_AB_N{n:05d}")


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    for cond in ("A", "B"):
        stacked(n, cond)
    compare(n)
    print(f"[図] {FIG} に出力（PNG+SVG）")


if __name__ == "__main__":
    main()
