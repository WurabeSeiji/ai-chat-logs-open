#!/usr/bin/env python3
"""E-M9 図化 v1：f(t)（親平面外占有＝拡大座標）と rank_Q(t)（方向数）の時系列

各 N について、対照（既存 build_init 初期値）＋倍音海の各段を再走行し、
f(t) を毎 step、rank_Q(t) を SAMPLE 間隔で記録して2段組の図にする。
色: 対照=黒、破れ族=青系、等振幅族(N−1)=橙系。crossing を縦破線で表示。
測定は駆動スクリプトと同一（abl read-only import・既存ファイル無改変）。

使い方: python3 plot_paper8_em9_v1.py 5   （N=5 → plot_em9_N5.png）
        python3 plot_paper8_em9_v1.py 40  （N=40 → plot_em9_N40.png）
        python3 plot_paper8_em9_v1.py 300 （N=300 → plot_em9_N300.png）
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_m9p", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("mph_m9p", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

XMAX = 12000
CFG = {5: {"H": 8, "seed": 40260801},
       40: {"H": 4, "seed": 40260802},
       300: {"H": 4, "seed": 40260803}}


def trajectory(n, v0, wp):
    """駆動スクリプトと同一の測定（f は毎step、rank_Q は SAMPLE 間隔）。"""
    sample_ev = abl.SAMPLE[n]
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)
    gr0 = abl.gram_reduce(sys_lr, v0)
    _, B0, _, _, _ = abl.dominant_plane(sys_lr, gr0)

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    Z = v0.copy()
    fs = np.zeros(XMAX + 1)
    rk_t, rk_v = [], []
    crossing = None
    for t in range(XMAX + 1):
        fs[t] = fval(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
        if t % sample_ev == 0 or t == XMAX:
            gr = abl.gram_reduce(sys_lr, Z)
            _, Bdom, _, _, _ = abl.dominant_plane(sys_lr, gr)
            qs = abl.qsv4(B0, Bdom)
            rk_t.append(t); rk_v.append(int(np.sum(qs > abl.Q_REL_TAU * qs[0])))
        if t < XMAX:
            Z, wp = abl.evolve(sys_lr, Z, wp)
    return fs, np.array(rk_t), np.array(rk_v), crossing


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cfg = CFG[n]
    H = cfg["H"]

    print(f"N={n}: 対照を走行中…", flush=True)
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
    runs = [("対照（既存make_parent）", "control", trajectory(n, Z0, wp0.copy()))]

    print(f"N={n}: 倍音海 {H} 段を走行中…", flush=True)
    Z, info = mph.make_parent_harmonic(n, H, cfg["seed"], iters=2000, restarts=10, tol=1e-12)
    for h in range(H):
        lv = info["levels"][h]
        fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
        v0 = Z[:, h] * np.sqrt(H)
        wp = np.random.default_rng(90000 + h).normal(size=len(v0))
        runs.append((f"段n={h+1} σ₁={lv['sigma1']:.4f}", fam, trajectory(n, v0, wp)))
        print(f"  段n={h+1}（{fam}）完了", flush=True)

    blues = plt.get_cmap("winter")
    oranges = plt.get_cmap("autumn")
    nb = sum(1 for _, f, _ in runs if f == "broken")
    no = sum(1 for _, f, _ in runs if f == "N-1")
    ib = io = 0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True,
                                    height_ratios=[2, 1])
    for label, fam, (fs, rt, rv, crossing) in runs:
        if fam == "control":
            color, lw, z = "black", 2.0, 5
        elif fam == "broken":
            color, lw, z = blues(0.15 + 0.7 * ib / max(nb, 1)), 1.4, 4; ib += 1
        else:
            color, lw, z = oranges(0.1 + 0.6 * io / max(no, 1)), 1.4, 3; io += 1
        ax1.plot(np.arange(XMAX + 1), fs, color=color, lw=lw, zorder=z,
                 label=f"{label}" + (f"  crossing={crossing}" if crossing else "  crossingなし"))
        if crossing:
            ax1.axvline(crossing, color=color, ls="--", lw=0.6, alpha=0.5, zorder=1)
        ax2.step(rt, rv, where="post", color=color, lw=lw, zorder=z)

    ax1.axhline(0.05, color="red", ls=":", lw=0.8, label="crossing 閾値 f=0.05")
    ax1.set_ylabel("f（親平面外占有）")
    ax1.set_title(f"E-M9  N={n}（M={n*(n-1)//2}）種なし条件A・XMAX={XMAX}\n"
                  "黒=対照 / 青系=破れ族（σ₁<N−1）/ 橙系=等振幅族（σ₁=N−1）")
    ax1.legend(fontsize=7, loc="center right")
    ax1.set_ylim(-0.03, 1.02)
    ax2.set_ylabel("rank_Q（方向数）")
    ax2.set_xlabel("step")
    ax2.set_yticks(range(0, 6))
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out = HERE / f"plot_em9_N{n}.png"
    fig.savefig(out, dpi=150)
    print(f"saved: {out.name}")


if __name__ == "__main__":
    main()
