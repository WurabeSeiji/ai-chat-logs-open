#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E-M9 図（第8論文 図生成コードと同一形式）。解釈なし。指示外の図は作らない。

形式の複製元: code/make_preliminary_seed_ablation_figures_v1.py（fig01〜03、無改変 read-only）
    fig01: f_outside_parent
    fig02: q3, q4
    fig03: rank_Q
    共通横軸 絶対step 0..55000, 5000刻み, crossing 不動。lw=0.9, dpi=130,
    ylim/yticks/ラベル/タイトル書式は原本と同一。
系列: 対照（条件A＝既存 make_parent 初期値）＝原本条件Aと同色 #7F7F7F。
    倍音海の各段: 破れ族=青系, 等振幅族(N−1)=橙系（凡例に σ₁ と族を明記）。
観測: 原本 run() と同一（f/q3/q4/rank_Q を SAMPLE[n] 間隔で記録、XMAX=55000）。

使い方: python3 make_paper8_em9_figures_v1.py <N>
出力: fig01_f_em9_N{n:05d}.png, fig02_q3q4_em9_N{n:05d}.png, fig03_rankQ_em9_N{n:05d}.png
      ＋ 走行データ em9_traj_N{n:05d}.npz（再図化用）
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_m9f", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("mph_m9f", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

XMAX = 55000                       # 原本と同一（共通横軸）
CFG = {5: {"H": 8, "seed": 40260801},
       40: {"H": 4, "seed": 40260802},
       300: {"H": 4, "seed": 40260803}}
COL_CONTROL = "#7F7F7F"            # 原本 条件A と同色
BLUES = ["#1F3F8F", "#2E6DB4", "#4C9FD0", "#7BC4E0"]
ORANGES = ["#F58518", "#FDB45C", "#E4762A", "#FFD08A"]


def sample_run(n, v0, wp):
    """原本 run()（96-173行）と同一の観測量を SAMPLE[n] 間隔で記録（条件A・注入なし）。"""
    se = abl.SAMPLE[n]
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    gr0 = abl.gram_reduce(sys_lr, v0)
    _, B0, _, _, _ = abl.dominant_plane(sys_lr, gr0)
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    Z = v0.copy()
    ts, fs, q3s, q4s, rks = [], [], [], [], []
    crossing = None
    t = 0
    while True:
        f = fval(Z)
        if crossing is None and f > 0.05:
            crossing = t
        if t % se == 0 or t == XMAX:
            gr = abl.gram_reduce(sys_lr, Z)
            _, Bdom, _, _, _ = abl.dominant_plane(sys_lr, gr)
            qs = abl.qsv4(B0, Bdom)
            ts.append(t); fs.append(f)
            q3s.append(qs[2]); q4s.append(qs[3])
            rks.append(int(np.sum(qs > abl.Q_REL_TAU * qs[0])))
        if t >= XMAX:
            break
        Z, wp = abl.evolve(sys_lr, Z, wp); t += 1
    return (np.array(ts), np.array(fs), np.array(q3s), np.array(q4s),
            np.array(rks), crossing)


def setx(ax):
    ax.set_xlim(0, XMAX); ax.set_xticks(np.arange(0, XMAX + 1, 5000))


def main():
    n = int(sys.argv[1])
    cfg = CFG[n]; H = cfg["H"]
    series = []          # (label, color, data)

    print(f"N={n}: 対照走行（XMAX={XMAX}）…", flush=True)
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
    series.append(("control (cond A init)", COL_CONTROL, sample_run(n, Z0, wp0.copy())))
    print(f"  control 完了 crossing={series[0][2][5]}", flush=True)

    Zh, info = mph.make_parent_harmonic(n, H, cfg["seed"], iters=2000, restarts=10, tol=1e-12)
    ib = io = 0
    for h in range(H):
        lv = info["levels"][h]
        fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
        if fam == "broken":
            color = BLUES[ib % len(BLUES)]; ib += 1
        else:
            color = ORANGES[io % len(ORANGES)]; io += 1
        v0 = Zh[:, h] * np.sqrt(H)
        wp = np.random.default_rng(90000 + h).normal(size=len(v0))
        d = sample_run(n, v0, wp)
        series.append((f"harmonic n={h+1} ({fam}, σ₁={lv['sigma1']:.4f})", color, d))
        print(f"  段n={h+1}（{fam}）完了 crossing={d[5]}", flush=True)

    np.savez(HERE / f"em9_traj_N{n:05d}.npz",
             **{f"s{i}_{k}": arr for i, (_, _, d) in enumerate(series)
                for k, arr in zip(("t", "f", "q3", "q4", "rank"), d[:5])},
             labels=[s[0] for s in series], colors=[s[1] for s in series],
             crossings=[(-1 if s[2][5] is None else s[2][5]) for s in series])

    cr_ctrl = series[0][2][5]

    # fig01 f（原本 64-71行と同一形式）
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, color, (t, f, q3, q4, rk, cr) in series:
        ax.plot(t, f, lw=0.9, color=color, label=label)
    if cr_ctrl is not None:
        ax.axvline(cr_ctrl, color="k", ls=":", lw=0.8)
    setx(ax); ax.set_ylim(0, 1)
    ax.set_xlabel("step (absolute)"); ax.set_ylabel("f_outside_parent = 1 - E_P1")
    ax.set_title(f"N={n} fig01: f_outside_parent (E-M9 harmonic-sea initial states)")
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(HERE / f"fig01_f_em9_N{n:05d}.png", dpi=130); plt.close(fig)

    # fig02 q3,q4（原本 73-82行と同一形式）
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for label, color, (t, f, q3, q4, rk, cr) in series:
        a1.plot(t, q3, lw=0.9, color=color, label=label)
        a2.plot(t, q4, lw=0.9, color=color)
    for a in (a1, a2):
        if cr_ctrl is not None:
            a.axvline(cr_ctrl, color="k", ls=":", lw=0.8)
        setx(a)
    a1.set_ylabel("q3"); a2.set_ylabel("q4"); a2.set_xlabel("step (absolute)")
    a1.legend(fontsize=7)
    a1.set_title(f"N={n} fig02: q3, q4 (E-M9 harmonic-sea initial states)")
    fig.tight_layout(); fig.savefig(HERE / f"fig02_q3q4_em9_N{n:05d}.png", dpi=130); plt.close(fig)

    # fig03 rank_Q（原本 84-90行と同一形式）
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, color, (t, f, q3, q4, rk, cr) in series:
        ax.plot(t, rk, lw=0.9, color=color, label=label)
    if cr_ctrl is not None:
        ax.axvline(cr_ctrl, color="k", ls=":", lw=0.8)
    setx(ax); ax.set_ylim(-0.2, 4.5); ax.set_yticks([0, 2, 4])
    ax.set_xlabel("step (absolute)"); ax.set_ylabel("rank_Q")
    ax.set_title(f"N={n} fig03: rank_Q (E-M9 harmonic-sea initial states)")
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(HERE / f"fig03_rankQ_em9_N{n:05d}.png", dpi=130); plt.close(fig)

    print(f"saved: fig01/02/03_em9_N{n:05d}.png", flush=True)


if __name__ == "__main__":
    main()
