#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E-M9 方向3・4占有図（Stage A2c figure03 と同一の観測量・同一形式）

観測量の複製元（無改変 read-only）:
    paper8_stage_A2c_direction_lineage_N5/replay_and_extract_bases.py 150-173行
      B0   = dominant_plane(gram_reduce(v0))     初期支配平面（方向1・2）
      D34  = s4_new_dirs(B0, Bdom(t))            新方向（方向3・4）
      disp = QR(B_rot(B_rotᵀD34)) を align_2d で前step と連続化
      direction_3_occupation = occ(disp[:,0], Z)/‖Z‖²   （方向4 も同様）
    図形式: make_figures.py 110-119行（semilogy, lw=0.9, 表示床 1e-34,
      xlabel "absolute step", ylabel "occupation (display floor 1e-34)"）

系列: 対照（条件A初期値）＋倍音海の各段。step 0..5000（A2c と同一範囲）、
      サンプル間隔 SAMPLE[n]。各系列を1パネル（figure03 形式）とし格子に並べる。

判定の読み方（論文と同一）: 新方向3・4の占有が crossing 後に立ち上がり
      rank_Q=4（S4=方向3本＋完結軸1本の4軸）が完成すれば「三方向の誕生」。
      立ち上がらなければ方向は生まれていない。
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

spec = importlib.util.spec_from_file_location("abl_d34", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
color = importlib.import_module("run_paper7_5color_timeseries")  # abl が sys.path 設定済み
spec2 = importlib.util.spec_from_file_location("mph_d34", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

TMAX = 5000                    # A2c と同一範囲
CFG = {5: {"H": 8, "seed": 40260801},
       40: {"H": 4, "seed": 40260802},
       300: {"H": 4, "seed": 40260803}}


def direction34_series(n, v0, wp):
    """A2c replay 150-173行と同一の方向3・4占有時系列（条件A・注入なし）。"""
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

    Z = v0.copy()
    prev_disp = None
    ts, d3s, d4s = [], [], []
    crossing = None
    for t in range(TMAX + 1):
        f = fval(Z)
        if crossing is None and f > 0.05:
            crossing = t
        if t % se == 0 or t == TMAX:
            gr = abl.gram_reduce(sys_lr, Z)
            _, Bdom, _, _, _ = abl.dominant_plane(sys_lr, gr)
            D34 = color.s4_new_dirs(B0, Bdom)
            projected = B_rot @ (B_rot.T @ D34)
            disp, _ = np.linalg.qr(projected)
            disp = color.align_2d(prev_disp, disp[:, :2])
            prev_disp = disp
            total = float(np.real(np.vdot(Z, Z)))
            d3 = abl.occ(disp[:, [0]], Z) / total
            d4 = abl.occ(disp[:, [1]], Z) / total
            ts.append(t); d3s.append(d3); d4s.append(d4)
        if t < TMAX:
            Z, wp = abl.evolve(sys_lr, Z, wp)
    return np.array(ts), np.array(d3s), np.array(d4s), crossing


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cfg = CFG[n]; H = cfg["H"]

    series = []
    print(f"N={n}: 対照走行…", flush=True)
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
    series.append(("control (cond A init)", direction34_series(n, Z0, wp0.copy())))
    Zh, info = mph.make_parent_harmonic(n, H, cfg["seed"], iters=2000, restarts=10, tol=1e-12)
    for h in range(H):
        lv = info["levels"][h]
        fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
        v0 = Zh[:, h] * np.sqrt(H)
        wp = np.random.default_rng(90000 + h).normal(size=len(v0))
        series.append((f"harmonic n={h+1} ({fam}, σ₁={lv['sigma1']:.4f})",
                       direction34_series(n, v0, wp)))
        print(f"  段n={h+1}（{fam}）完了", flush=True)

    ncols = 3
    nrows = int(np.ceil(len(series) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.6 * nrows),
                              sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    for i, (label, (ts, d3, d4, cr)) in enumerate(series):
        ax = axes[i // ncols][i % ncols]
        ax.semilogy(ts, np.maximum(d3, 1e-34), lw=0.9, label="direction 3 occupation")
        ax.semilogy(ts, np.maximum(d4, 1e-34), lw=0.9, label="direction 4 occupation")
        if cr is not None:
            ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_title(label + (f"  crossing={cr}" if cr is not None else "  crossingなし"),
                     fontsize=9)
        ax.set_ylim(1e-34, 2.0)
        if i // ncols == nrows - 1:
            ax.set_xlabel("absolute step")
        if i % ncols == 0:
            ax.set_ylabel("occupation (display floor 1e-34)")
        if i == 0:
            ax.legend(fontsize=7)
    for j in range(len(series), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(f"N={n} E-M9: existing aligned direction 3/4 occupations "
                 f"(Stage A2c figure03 形式・step 0..{TMAX})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = HERE / f"fig_em9_direction34_N{n:05d}.png"
    fig.savefig(out, dpi=130)
    print(f"saved: {out.name}", flush=True)


if __name__ == "__main__":
    main()
