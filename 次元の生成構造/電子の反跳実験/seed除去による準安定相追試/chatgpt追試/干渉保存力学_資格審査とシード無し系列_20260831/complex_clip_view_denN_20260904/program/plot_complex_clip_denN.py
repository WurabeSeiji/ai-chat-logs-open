#!/usr/bin/env python3
"""den=N・step2000 状態の複素平面クリップ図（抑圧max 正規化）。

各 N (3..33) の den=N ケースについて step2000 の状態 Z（全成分・複素数）を読み、
抑圧max で正規化した複素平面散布図を 31 パネルで描く。

- 生存/抑圧の判定: |Z| >= 0.1 * |Z|max を生存（凝縮成分）、それ未満を抑圧成分とする。
- 抑圧max = 抑圧成分の最大振幅。
- 各パネルの表示範囲は実部・虚部とも ±1（= ±抑圧max）に固定。
  これを超える成分（凝縮成分）は視野外にトリムされる。
- 赤破線は |z| = 抑圧max（半径1）の参照円。

副産物として、全分母（den=N-2..N+2, 124）の抑圧成分統計
suppressed_summary_N3_N33_step2000.csv も出力する（チャット報告値の正本）。
"""
import argparse
import csv
import glob
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPLIT = 0.1  # 生存/抑圧の判定しきい値（|Z|max に対する比）


def final_state(path):
    z = np.load(path)
    return z["Z"][-1]


def split_state(v):
    """(生存成分, 抑圧成分, 抑圧max or None) を返す。"""
    a = np.abs(v)
    mask = a >= SPLIT * a.max()
    sup = a[~mask]
    smax = sup.max() if len(sup) else None
    return mask, sup, smax


def write_suppressed_summary(src, out_csv):
    rows = []
    for f in sorted(glob.glob(os.path.join(src, "hm_N*_den_*_states_2000.npz"))):
        m = re.search(r"hm_N(\d+)_den_(\d+)_states_2000", f)
        N, den = int(m.group(1)), int(m.group(2))
        v = final_state(f)
        mask, sup, smax = split_state(v)
        a2 = np.abs(v) ** 2
        pr = a2.sum() ** 2 / (a2 ** 2).sum()
        rows.append(
            [N, den, len(v), int(mask.sum()), f"{pr:.6f}",
             f"{smax:.6e}" if smax is not None else "",
             f"{sup.min():.6e}" if len(sup) else "",
             len(sup)]
        )
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_csv, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["N", "den", "components", "survivors", "PR",
                    "sup_max", "sup_min", "sup_count"])
        w.writerows(rows)
    print("saved:", out_csv)


def plot_complex_clip(src, out_png, out_csv):
    plt.rcParams["font.family"] = ["Hiragino Sans", "Arial"]
    Ns = list(range(3, 34))
    nrow, ncol = 6, 6
    fig, axes = plt.subplots(nrow, ncol, figsize=(18, 18))
    th = np.linspace(0, 2 * np.pi, 361)
    rows = []
    for k, N in enumerate(Ns):
        ax = axes[k // ncol][k % ncol]
        v = final_state(os.path.join(src, f"hm_N{N}_den_{N}_states_2000.npz"))
        mask, sup, smax = split_state(v)
        if smax is None:
            ax.text(0, 0, "凝縮なし\n(抑圧max 未定義)",
                    ha="center", va="center", fontsize=10, color="gray")
            rows.append([N, N, len(v), int(mask.sum()), "", 0])
        else:
            w = v / smax
            ax.scatter(w.real, w.imag, s=12, alpha=0.65,
                       color="#1f6feb", edgecolors="none")
            ax.plot(np.cos(th), np.sin(th), ls="--", lw=0.8, color="#c0392b")
            trimmed = int(np.sum((np.abs(w.real) > 1) | (np.abs(w.imag) > 1)))
            rows.append([N, N, len(v), int(mask.sum()),
                         f"{smax:.6e}", trimmed])
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_aspect("equal")
        ax.axhline(0, lw=0.4, color="gray")
        ax.axvline(0, lw=0.4, color="gray")
        ax.set_title(f"N={N} (den={N})", fontsize=10)
        ax.tick_params(labelsize=7)
    for k in range(len(Ns), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(
        "den=N・step2000 状態の複素平面図"
        "（抑圧max 正規化、表示範囲 ±1 でトリム、赤破線 = |z| = 抑圧max）",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_png, dpi=150)
    print("saved:", out_png)

    with open(out_csv, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["N", "den", "components", "survivors", "sup_max",
                    "trimmed_points"])
        w.writerows(rows)
    print("saved:", out_csv)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--src",
        default="../ChatGPT_denominator_controls_N3_N33_legacyparent_20260903"
                "/results_2000steps",
        help="状態 npz（hm_N*_den_*_states_2000.npz）のあるフォルダ",
    )
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    plot_complex_clip(
        args.src,
        os.path.join(args.out, "fig_complex_clip_denN_N3_N33_step2000.png"),
        os.path.join(args.out, "supmax_denN_step2000.csv"),
    )
    write_suppressed_summary(
        args.src,
        os.path.join(args.out, "suppressed_summary_N3_N33_step2000.csv"),
    )


if __name__ == "__main__":
    main()
