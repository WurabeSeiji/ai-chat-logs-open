#!/usr/bin/env python3
"""親（自己無撞着円偏波固有モード）の星座図と数値表を生成する。

原本力学 run_n_scaling_lowrank_v1.py の make_parent /
zero_closure_kernel_seed をそのまま import して使う（再実装禁止規約）。

出力（lowN_metastable_result_v1/ 内）:
  parent_constellation_N{n:05d}_seed{s:03d}.png   複素平面上の星座図
  parent_constellation_N{n:05d}_seed{s:03d}.json  15 辺の複素初期値・振幅位相グループ・検算値

使い方:
  python3 make_parent_constellation_v1.py 6 --seed=0 --delta=1e-15
"""

import argparse
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_n_scaling_lowrank_v1 import (
    LowRankSystem,
    make_parent,
    zero_closure_kernel_seed,
    build_edges,
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "lowN_metastable_result_v1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--delta", type=float, default=1e-15)
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    sys_lr = LowRankSystem(args.n)
    parent, residual, sigmas = make_parent(
        sys_lr, rng, iters=1200, tol=args.tol, restarts=8
    )
    seed_vec = zero_closure_kernel_seed(sys_lr, rng)
    z0 = parent + args.delta * seed_vec
    z0 = z0 / np.linalg.norm(z0)

    ea, eb = build_edges(args.n)
    groups = {}
    for e in range(sys_lr.m):
        key = (round(abs(parent[e]), 6), round(float(np.angle(parent[e])), 6))
        groups.setdefault(key, []).append((int(ea[e]), int(eb[e])))

    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    for r in sorted({k[0] for k in groups}):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, ls="--", color="gray", lw=0.8))
    for (r, th), edges in groups.items():
        x, y = r * np.cos(th), r * np.sin(th)
        ax.plot([0, x], [0, y], color="tab:blue", lw=1.2)
        ax.plot(x, y, "o", ms=11, color="tab:red")
        ax.annotate(
            f"x{len(edges)}", (x, y),
            textcoords="offset points", xytext=(10, 8), fontsize=13,
        )
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    lim = 1.45 * max(abs(parent))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_title(
        f"parent constellation N={args.n} seed={args.seed} "
        f"(sigma_max={sigmas[0]:.4f})"
    )

    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = f"N{args.n:05d}_seed{args.seed:03d}"
    png_path = os.path.join(RESULT_DIR, f"parent_constellation_{tag}.png")
    json_path = os.path.join(RESULT_DIR, f"parent_constellation_{tag}.json")
    fig.savefig(png_path, dpi=110, bbox_inches="tight")
    plt.close(fig)

    payload = {
        "n": args.n,
        "m": sys_lr.m,
        "seed": args.seed,
        "delta": args.delta,
        "parent_residual": float(residual),
        "sigma_max": float(sigmas[0]),
        "norm2_z0": float(np.real(np.vdot(z0, z0))),
        "abs_ztz_z0": float(abs(complex(z0 @ z0))),
        "edges": [[int(a), int(b)] for a, b in zip(ea, eb)],
        "z0_real": [float(x) for x in z0.real],
        "z0_imag": [float(x) for x in z0.imag],
        "amplitude_phase_groups": [
            {
                "amplitude": r,
                "phase_rad": th,
                "count": len(edges),
                "edges": edges,
            }
            for (r, th), edges in sorted(groups.items())
        ],
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    print(f"親残差 {residual:.3e}  σ_max {sigmas[0]:.6f}")
    print(f"振幅位相グループ数 {len(groups)}")
    print(f"保存: {png_path}")
    print(f"保存: {json_path}")


if __name__ == "__main__":
    main()
