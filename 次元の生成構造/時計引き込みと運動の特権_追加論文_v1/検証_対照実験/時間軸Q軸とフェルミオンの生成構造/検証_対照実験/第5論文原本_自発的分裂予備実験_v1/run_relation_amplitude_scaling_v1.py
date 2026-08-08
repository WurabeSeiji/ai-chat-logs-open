#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=8〜300 の関係波準安定振幅を個別CSVへ測定する。

大N側は前論文O6と同じ低ランクエンジン、親構成、零閉鎖核種、
warm-start sigma_1 推定を使う。N<=12では正当性対照のため厳密sigma_1も選べる。

使用例:
    python3 run_relation_amplitude_scaling_v1.py 40 1e-15 --seed=0
    python3 run_relation_amplitude_scaling_v1.py 300 1e-15 --seed=0 --sub=20
"""

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
from run_n_scaling_lowrank_v1 import (  # noqa: E402
    LowRankSystem,
    make_parent,
    progress,
    zero_closure_kernel_seed,
)

RESULT_DIR = os.path.join(BASE_DIR, "relation_amplitude_scaling_result_v1")
CORE_SEED_BASE = 40260722


def participation_ratio(z):
    weights = np.abs(z) ** 2
    return float(np.sum(weights) ** 2 / np.sum(weights * weights))


def normalized_entropy(z):
    weights = np.abs(z) ** 2
    weights = weights / np.sum(weights)
    nz = weights[weights > 0.0]
    return float(-np.sum(nz * np.log(nz)) / np.log(len(weights)))


def complement_fraction(z, p, q):
    zp = z - p * (p @ z) - q * (q @ z)
    return float(np.real(np.vdot(zp, zp)) / np.real(np.vdot(z, z)))


def amplitude_metrics(z):
    m = len(z)
    root_m = math.sqrt(m)
    absolute = np.abs(z)
    q05, median, q95 = np.quantile(absolute, [0.05, 0.5, 0.95])
    pr_fraction = participation_ratio(z) / m
    return {
        "pr_over_m": pr_fraction,
        "entropy_normalized": normalized_entropy(z),
        "relation_abs_mean": float(np.mean(absolute)),
        "relation_abs_median": float(median),
        "relation_abs_q05": float(q05),
        "relation_abs_q95": float(q95),
        "scaled_mean_sqrt_m": float(np.mean(absolute) * root_m),
        "scaled_median_sqrt_m": float(median * root_m),
        "scaled_q05_sqrt_m": float(q05 * root_m),
        "scaled_q95_sqrt_m": float(q95 * root_m),
        "scaled_width_q90_sqrt_m": float((q95 - q05) * root_m),
    }


def tail_stats(rows, key, count):
    selected = rows[-min(count, len(rows)):]
    values = np.array([float(row[key]) for row in selected])
    taus = np.array([float(row["tau"]) for row in selected])
    slope = (
        float(np.polyfit(taus, values, 1)[0])
        if len(values) >= 3 and np.ptp(taus) > 0
        else 0.0
    )
    return {
        "n": len(values),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "q05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "q95": float(np.quantile(values, 0.95)),
        "slope_per_step": slope,
    }


def run(args):
    if args.n < 8:
        raise ValueError("この軽量走行器はN>=8用です。N=3〜7はlowN走行器を使います")
    sys_lr = LowRankSystem(args.n)
    rng = np.random.default_rng(CORE_SEED_BASE + 1000 * args.n + args.seed)
    start_parent = time.time()
    parent, residual, parent_sigmas = make_parent(
        sys_lr,
        rng,
        iters=1200,
        tol=args.tol,
        restarts=8,
    )
    parent_runtime = time.time() - start_parent
    if residual > args.parent_residual_max:
        raise RuntimeError(
            f"親構成残差 {residual:.3e} が上限 "
            f"{args.parent_residual_max:.3e} を超えました"
        )

    seed_vector = zero_closure_kernel_seed(sys_lr, rng)
    z = parent + args.delta * seed_vector
    z = z / np.linalg.norm(z)
    c0 = complex(z @ z)
    norm0 = float(np.real(np.vdot(z, z)))

    p = parent.real / np.linalg.norm(parent.real)
    q = parent.imag - (parent.imag @ p) * p
    q = q / np.linalg.norm(q)

    normalization = (
        "exact_sigma"
        if args.n <= args.exact_sigma_until
        else "warm_start_power_sigma"
    )
    wp = rng.normal(size=sys_lr.m)
    crossing_tau = None
    plateau_tau = None
    stable_streak = 0
    rows = []
    max_closure_deviation = 0.0
    max_norm_deviation = 0.0
    start_run = time.time()
    final_tau = args.cap

    for tau in range(args.cap + 1):
        f = complement_fraction(z, p, q)
        if crossing_tau is None and f > args.crossing:
            crossing_tau = tau
            progress(f"amplitude N={args.n} 閾値交差 tau={tau}")

        max_closure_deviation = max(
            max_closure_deviation, abs(complex(z @ z) - c0)
        )
        max_norm_deviation = max(
            max_norm_deviation,
            abs(float(np.real(np.vdot(z, z))) - norm0),
        )

        sys_lr.set_theta(np.angle(z))
        if tau % args.sub == 0 or tau == args.cap:
            metrics = amplitude_metrics(z)
            sigmas = sys_lr.sigma_spectrum()
            sigma_ratio = (
                float(sigmas[1] / sigmas[0])
                if len(sigmas) >= 2 and sigmas[0] > 0.0
                else ""
            )
            row = {
                "tau": tau,
                "f_initial_plane": f,
                **metrics,
                "sigma_1_exact": float(sigmas[0]) if len(sigmas) else "",
                "sigma_ratio_2_1": sigma_ratio,
                "abs_ztz": abs(complex(z @ z)),
                "norm2": float(np.real(np.vdot(z, z))),
            }
            rows.append(row)

            is_equalized = (
                metrics["pr_over_m"] >= 1.0 - args.pr_tolerance
                and abs(metrics["scaled_median_sqrt_m"] - 1.0)
                    <= args.median_tolerance
                and metrics["scaled_width_q90_sqrt_m"]
                    <= args.width_tolerance
            )
            stable_streak = stable_streak + 1 if is_equalized else 0
            enough_after_crossing = (
                crossing_tau is not None
                and tau >= crossing_tau + args.min_after_crossing
            )
            if (
                enough_after_crossing
                and stable_streak >= args.stable_records
            ):
                plateau_tau = tau - (args.stable_records - 1) * args.sub
                final_tau = tau
                progress(
                    f"amplitude N={args.n} 等振幅準安定判定 "
                    f"tau={plateau_tau}〜{tau}"
                )
                break

        if tau < args.cap:
            if normalization == "exact_sigma":
                sigmas_step = sys_lr.sigma_spectrum()
                sigma_step = float(sigmas_step[0])
            else:
                sigma_step, wp = sys_lr.sigma_max_power(wp)
            z = sys_lr.cayley_step(z, sigma_step)

        if tau > 0 and tau % args.progress_every == 0:
            latest = rows[-1] if rows else {}
            progress(
                f"amplitude N={args.n} tau={tau} f={f:.3e} "
                f"PR/M={latest.get('pr_over_m', float('nan')):.6f} "
                f"Amed*sqrt(M)="
                f"{latest.get('scaled_median_sqrt_m', float('nan')):.6f}"
            )

    run_runtime = time.time() - start_run
    os.makedirs(args.result_dir, exist_ok=True)
    tag = f"N{args.n:05d}_delta{args.delta:.0e}_seed{args.seed:03d}"
    csv_path = os.path.join(args.result_dir, f"amplitude_{tag}.csv")
    json_path = os.path.join(args.result_dir, f"summary_{tag}.json")
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    if plateau_tau is not None:
        # 判定成立時は、判定を満たした連続点だけを準安定窓とする。
        # sub が大きい高Nで固定個数の末尾を取ると、遷移前まで
        # 遡ってしまうためである。
        tail_count = min(args.stable_records, len(rows))
        tail_window = "equalization_streak"
    else:
        tail_count = min(args.tail_records, len(rows))
        tail_window = "terminal_window_without_equalization_detection"
    tail = {
        key: tail_stats(rows, key, tail_count)
        for key in (
            "f_initial_plane",
            "pr_over_m",
            "entropy_normalized",
            "relation_abs_mean",
            "relation_abs_median",
            "relation_abs_q05",
            "relation_abs_q95",
            "scaled_mean_sqrt_m",
            "scaled_median_sqrt_m",
            "scaled_width_q90_sqrt_m",
            "sigma_ratio_2_1",
        )
    }
    measured = tail["relation_abs_median"]["mean"]
    exact_prediction = 1.0 / math.sqrt(sys_lr.m)
    summary = {
        "experiment": "relation_amplitude_scaling_v1",
        "n": args.n,
        "m": sys_lr.m,
        "delta": args.delta,
        "seed": args.seed,
        "normalization": normalization,
        "cap": args.cap,
        "steps_run": final_tau,
        "sub": args.sub,
        "tol": args.tol,
        "parent_residual": residual,
        "parent_rank_planes": int(len(parent_sigmas)),
        "parent_runtime_sec": parent_runtime,
        "abs_ztz_initial": abs(c0),
        "max_closure_deviation": max_closure_deviation,
        "max_norm_deviation": max_norm_deviation,
        "crossing_tau": crossing_tau,
        "plateau_tau": plateau_tau,
        "equalization_criteria": {
            "pr_tolerance": args.pr_tolerance,
            "median_tolerance": args.median_tolerance,
            "width_tolerance": args.width_tolerance,
            "stable_records": args.stable_records,
            "min_after_crossing": args.min_after_crossing,
        },
        "tail_window": tail_window,
        "tail_records": tail_count,
        "tail": tail,
        "duality": {
            "measured_relation_abs_median": measured,
            "prediction_1_over_sqrt_m": exact_prediction,
            "measured_times_sqrt_m": measured * math.sqrt(sys_lr.m),
            "relative_error": measured / exact_prediction - 1.0,
            "n_times_measured": args.n * measured,
        },
        "run_runtime_sec": run_runtime,
        "trajectory_csv": os.path.basename(csv_path),
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return summary


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("delta", nargs="?", type=float, default=1e-15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--cap", type=int, default=12000)
    parser.add_argument("--sub", type=int, default=10)
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--parent-residual-max", type=float, default=1e-10)
    parser.add_argument("--exact-sigma-until", type=int, default=12)
    parser.add_argument("--crossing", type=float, default=0.05)
    parser.add_argument("--min-after-crossing", type=int, default=1500)
    parser.add_argument("--pr-tolerance", type=float, default=1e-5)
    parser.add_argument("--median-tolerance", type=float, default=2e-3)
    parser.add_argument("--width-tolerance", type=float, default=1e-2)
    parser.add_argument("--stable-records", type=int, default=20)
    parser.add_argument("--tail-records", type=int, default=200)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument("--result-dir", default=RESULT_DIR)
    return parser.parse_args()


def main():
    run(parse_args())


if __name__ == "__main__":
    main()
