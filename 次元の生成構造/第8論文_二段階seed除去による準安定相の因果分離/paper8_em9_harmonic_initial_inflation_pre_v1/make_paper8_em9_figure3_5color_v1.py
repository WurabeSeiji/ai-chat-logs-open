#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E-M9 Figure3（5色, log）——note公開図と同一形式・同一観測量。解釈なし。

複製元（無改変 read-only）:
    観測量: run_paper7_5color_timeseries.py run() 112-142行
        P1 = direction_1+direction_2（親平面2軸）, direction_3/4（align_2d 済み新方向）,
        remaining other-rotation, kernel, 黒線 f = 1 − E_P1
    図形式: make_paper7_figures.py 100-116行 figure3_compare（5色 semilogy, clip 1e-6,
        lw0.8/0.9, crossing点線, 共通横軸 0..55000・5000刻み, COLORS/suptitle 同一）

行構成: 各 N について 対照＋倍音海全段（公開図の N 段を系列段に置き換え）。
使い方: python3 make_paper8_em9_figure3_5color_v1.py <N>
出力: fig3_5color_em9_N{n:05d}.png/.svg
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_f3", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
color = importlib.import_module("run_paper7_5color_timeseries")
spec2 = importlib.util.spec_from_file_location("mph_f3", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

XMAX = 55000
CFG = {5: {"H": 8, "seed": 40260801},
       40: {"H": 4, "seed": 40260802},
       300: {"H": 4, "seed": 40260803}}
COLORS = ["#4C78A8", "#E45756", "#F58518", "#B0B0B0", "#54A24B"]  # P1, d3, d4, other残, 核
LABELS = ["P1 (dominant plane)", "direction 3", "direction 4",
          "remaining other-rotation", "kernel"]


def five_color_run(n, v0, wp):
    """run_paper7_5color_timeseries.run() 112-142行の複製（初期状態を引数化）。"""
    se = abl.SAMPLE[n]
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    if n <= 40:
        _, B_p1, B_rot, _ = abl.parent_plane_split_exact(sys_lr, v0)
    else:
        _, B_p1, B_rot, _, _ = abl.parent_plane_split_approx(sys_lr, v0, abl.SIG_REL)
    gr0 = abl.gram_reduce(sys_lr, v0)
    _, B0, _, _, _ = abl.dominant_plane(sys_lr, gr0)
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    Zr = v0.copy()
    f_prev = None
    ts, bands, fs = [], [[], [], [], [], []], []
    crossing = None
    t = 0
    while True:
        fcur = fval(Zr)
        if crossing is None and fcur > 0.05:
            crossing = t
        if t % se == 0 or t == XMAX:
            totZ = float(np.real(np.conj(Zr) @ Zr))
            E_P1 = abl.occ(B_p1, Zr)
            E_other = abl.occ(B_rot, Zr)
            E_ker = totZ - E_P1 - E_other
            gr = abl.gram_reduce(sys_lr, Zr)
            _, Bdom, _, _, _ = abl.dominant_plane(sys_lr, gr)
            e34 = color.s4_new_dirs(B0, Bdom)
            proj = B_rot @ (B_rot.T @ e34)
            fq, _ = np.linalg.qr(proj)
            f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
            f34 = color.align_2d(f_prev, f34); f_prev = f34
            E_d3 = abl.occ(f34[:, [0]], Zr)
            E_d4 = abl.occ(f34[:, [1]], Zr) if f34.shape[1] > 1 else 0.0
            E_rem = max(0.0, E_other - E_d3 - E_d4)
            ts.append(t)
            for arr, val in zip(bands, (E_P1, E_d3, E_d4, E_rem, E_ker)):
                arr.append(val / totZ)
            fs.append(1.0 - E_P1 / totZ)
        if t >= XMAX:
            break
        Zr, wp = abl.evolve(sys_lr, Zr, wp); t += 1
    return (np.array(ts), [np.array(b) for b in bands], np.array(fs), crossing)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cfg = CFG[n]; H = cfg["H"]
    series = []
    print(f"N={n}: 対照走行…", flush=True)
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
    series.append(("control", five_color_run(n, Z0, wp0.copy())))
    print(f"  control 完了 crossing={series[0][1][3]}", flush=True)
    Zh, info = mph.make_parent_harmonic(n, H, cfg["seed"], iters=2000, restarts=10, tol=1e-12)
    for h in range(1, H + 1):
        lv = info["levels"][h - 1]
        fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
        v0 = Zh[:, h - 1] * np.sqrt(H)
        wp = np.random.default_rng(90000 + (h - 1)).normal(size=len(v0))
        d = five_color_run(n, v0, wp)
        series.append((f"n={h} ({fam})", d))
        print(f"  段n={h}（{fam}）完了 crossing={d[3]}", flush=True)

    rows = len(series)
    fig, axes = plt.subplots(rows, 1, figsize=(11, 2.6 * rows), sharex=True, squeeze=False)
    for ax, (label, (t, bands, f, cr)) in zip(axes[:, 0], series):
        for band, c in zip(bands, COLORS):
            ax.semilogy(t, np.clip(band, 1e-6, None), lw=0.8, color=c)
        ax.semilogy(t, np.clip(f, 1e-6, None), "k-", lw=0.9)
        if cr is not None:
            ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlim(0, XMAX); ax.set_xticks(np.arange(0, XMAX + 1, 5000))
        ax.set_ylabel(label, fontsize=8)
    axes[0, 0].legend(LABELS + ["f"], fontsize=6, loc="center right")
    axes[-1, 0].set_xlabel("step (absolute)")
    fig.suptitle(f"Figure3 compare (5-color, log) — common axis (E-M9 N={n}: "
                 "control + harmonic-sea levels)")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(HERE / f"fig3_5color_em9_N{n:05d}.png", dpi=130)
    fig.savefig(HERE / f"fig3_5color_em9_N{n:05d}.svg")
    plt.close(fig)
    print(f"saved: fig3_5color_em9_N{n:05d}.png/.svg", flush=True)


if __name__ == "__main__":
    main()
