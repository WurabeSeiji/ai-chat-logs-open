#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""平面流入の直接計測：状態エネルギーを親生成子の回転平面へ分解して追う。

親生成子 K(arg v)（実反対称 M×M）を σ値でグループ化した回転平面
（+ σ=0 核）へ状態 Z のエネルギーを完全射影し、時間発展で
「支配平面 P₁ から第2平面・核へ振幅が流れ込む」ことを数値化する。

力学は原本 run_spontaneous_splitting_largeN_v1.run() の忠実複製（不変更）。
対照テスト：全ステップで f == 1 − (P₁エネルギー比) を検査（両者は恒等のはず）。
さらに記録した f を metastable_series の正本 fcurve とも突き合わせる。

出力（plane_flow_result_v1/）:
    planeflow_N#####.csv       tau, f, frac_σ... （各回転平面群のエネルギー比）
    planeflow_N#####.json      平面群の (σ:次元)、対照テスト結果
    planeflow_N#####.png       積み上げ面図（平面間のエネルギー移動）

使い方: python3 run_plane_flow_v1.py 5 --after=20000 --cap=30000 --record-every=20
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import (
    LowRankSystem, make_parent, zero_closure_kernel_seed,
)

RESULT_DIR = HERE / "plane_flow_result_v1"
REFERENCE_DIR = HERE / "metastable_series_result_v1"


def parent_plane_bases(sys_lr, v, tol=1e-9):
    """親生成子 K(arg v) の回転平面群（σ値ごと）+ σ=0 核の正規直交実基底。"""
    M = sys_lr.m
    sys_lr.set_theta(np.angle(v))
    K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
    w, V = np.linalg.eig(K)
    sigs = w.imag
    raw = {}
    for i in np.where(sigs > tol)[0]:
        key = round(float(sigs[i]), 4)
        raw.setdefault(key, []).extend([np.real(V[:, i]), np.imag(V[:, i])])
    bases = {}
    for k, vs in raw.items():
        Q, _ = np.linalg.qr(np.column_stack(vs))
        bases[k] = Q
    U, S, _ = np.linalg.svd(K)
    ker = U[:, S < tol]
    if ker.shape[1] > 0:
        bases[0.0] = ker
    return bases, K


def run(n, delta, seed, cap, after, record_every, tol=1e-12):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + seed)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=tol)
    bases, K = parent_plane_bases(sys_lr, v)
    order = sorted(bases, reverse=True)
    dims = {k: bases[k].shape[1] for k in order}

    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    def energy_fracs(Z):
        a, b = Z.real, Z.imag
        tot = a @ a + b @ b
        return {k: (np.sum((B.T @ a) ** 2) + np.sum((B.T @ b) ** 2)) / tot
                for k, B in bases.items()}

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    rows = []
    max_id_dev = 0.0        # |f - (1 - P1比)| の対照
    crossed = None
    wp = rng.normal(size=sys_lr.m)
    for t in range(cap + 1):
        f = fval(Z)
        if crossed is None and f > 0.05:
            crossed = t
        if crossed is not None and t >= crossed + after:
            break
        if t % record_every == 0:
            e = energy_fracs(Z)
            max_id_dev = max(max_id_dev, abs(f - (1.0 - e[order[0]])))
            rows.append((t, f, *[e[k] for k in order]))
        sys_lr.set_theta(np.angle(Z))
        se, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, se)

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"N{n:05d}"
    with open(RESULT_DIR / f"planeflow_{tag}.csv", "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["tau", "f"] + [f"frac_sigma_{k}" for k in order])
        wtr.writerows(rows)

    # 対照2: f が正本 fcurve と一致
    ref_csv = REFERENCE_DIR / f"fcurve_{tag}_delta{delta:.0e}_seed{seed}.csv"
    ref_dev = None
    if ref_csv.exists():
        ref = {}
        with open(ref_csv) as fh:
            for r in csv.DictReader(fh):
                ref[int(r["tau"])] = float(r["f"])
        ref_dev = max(abs(row[1] - ref[row[0]]) for row in rows if row[0] in ref)

    control = {"identity_f_eq_1_minus_P1_maxdev": float(max_id_dev),
               "identity_passed": bool(max_id_dev < 1e-12),
               "f_vs_reference_maxdev": (None if ref_dev is None else float(ref_dev)),
               "f_vs_reference_passed": bool(ref_dev is None or ref_dev < 1e-12)}
    meta = {"n": n, "m": sys_lr.m, "delta": delta, "seed": seed,
            "crossing_tau": crossed, "parent_residual": residual,
            "plane_groups_sigma_dim": {str(k): dims[k] for k in order},
            "control_test": control}
    with open(RESULT_DIR / f"planeflow_{tag}.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)

    # 図：積み上げ面（平面群ごとのエネルギー比）
    arr = np.array(rows)
    taus = arr[:, 0]
    fracs = arr[:, 2:]
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    labels = [f"σ={k:.3f} plane ({dims[k]}D)" if k > 0 else f"σ=0 kernel ({dims[k]}D)"
              for k in order]
    ax.stackplot(taus, *fracs.T, labels=labels, alpha=0.85)
    if crossed is not None:
        ax.axvline(crossed, color="k", ls=":", lw=1.0)
        ax.text(crossed, 1.01, f" crossing τ={crossed}", fontsize=8)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("energy fraction per generator rotation plane")
    ax.set_ylim(0, 1)
    ax.set_title(f"N={n}: amplitude flows out of dominant plane P1 into the "
                 f"2nd rotation plane and kernel")
    ax.legend(loc="center right", fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / f"planeflow_{tag}.png", dpi=150)
    plt.close(fig)

    print(json.dumps({"n": n, "crossing_tau": crossed,
                      "plane_groups": {str(k): dims[k] for k in order},
                      "control_test": control}, ensure_ascii=False))
    if not control["identity_passed"]:
        raise SystemExit(f"恒等対照失敗: {max_id_dev:.2e}")
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--delta", type=float, default=1e-15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cap", type=int, default=30000)
    ap.add_argument("--after", type=int, default=20000)
    ap.add_argument("--record-every", type=int, default=20)
    args = ap.parse_args()
    run(args.n, args.delta, args.seed, args.cap, args.after, args.record_every)


if __name__ == "__main__":
    main()
