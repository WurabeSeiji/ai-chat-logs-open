#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""原本力学を変えずに N=3〜7 の準安定域を長時間測定する。

各 N・delta・seed を独立に実行し、時系列 CSV と要約 JSON を生成する。

使用例:
    python3 run_lowN_metastable_scan_v1.py 3 1e-15 --seed=0 --cap=100000
    python3 run_lowN_metastable_scan_v1.py 7 1e-15 --seed=0 --cap=100000

原本との関係:
  - 位相差正弦生成子、スペクトルノルム正規化、Cayley 更新は
    run_n_scaling_lowrank_v1.py を import して使う。
  - 小 N では毎ステップ厳密 sigma_1 を使う。warm-start 冪反復は使わない。
  - 追加するのは読み取り専用の観測量と N=3 の厳密零閉鎖近傍種だけである。
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
)

RESULT_DIR = os.path.join(BASE_DIR, "lowN_metastable_result_v1")
CORE_SEED_BASE = 40260722


def parent_plane(v):
    """複素親ベクトルから実正規直交平面 p,q と円偏波表現を作る。"""
    p = v.real.copy()
    p = p / np.linalg.norm(p)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    vc = (p + 1j * q) / math.sqrt(2.0)
    # 元の v と反対のカイラリティを選ばないようにする。
    if abs(np.vdot(v, np.conj(vc))) > abs(np.vdot(v, vc)):
        q = -q
        vc = (p + 1j * q) / math.sqrt(2.0)
    return p, q, vc


def eigenmode_residual(sys_lr, v):
    sys_lr.set_theta(np.angle(v))
    kv = sys_lr.kmatvec(v)
    mu = float(np.real(np.conj(v) @ (1j * kv)))
    return float(np.linalg.norm(1j * kv - mu * v))


def n3_exact_zero_closure_state(v, p, q, rng, delta):
    """N=3 の1実次元補空間に厳密零閉鎖近傍種を置く。

    補方向 u に delta u を加えるだけでは Z^T Z=delta^2 となる。
    v^T conjugate(v)=1 を使い、-(delta^2/2) conjugate(v) を加えて相殺する。
    """
    u = rng.normal(size=len(v))
    u = u - (u @ p) * p - (u @ q) * q
    nu = np.linalg.norm(u)
    if nu < 1e-12:
        raise RuntimeError("N=3 補方向の構成に失敗しました")
    u = u / nu
    z = v + delta * u - 0.5 * delta * delta * np.conj(v)
    z = z / np.linalg.norm(z)
    return z, "n3_exact_transverse_with_conjugate_counterterm"


def parent_complement_zero_closure_state(v, p, q, rng, delta):
    """親平面の実直交補に零閉鎖正規直交対を置く（N>=4）。"""
    def project_out(x):
        return x - (x @ p) * p - (x @ q) * q

    u = project_out(rng.normal(size=len(v)))
    nu = np.linalg.norm(u)
    if nu < 1e-12:
        raise RuntimeError("親平面補方向 u の構成に失敗しました")
    u = u / nu
    w = project_out(rng.normal(size=len(v)))
    w = w - (w @ u) * u
    nw = np.linalg.norm(w)
    if nw < 1e-12:
        raise RuntimeError("親平面補方向 w の構成に失敗しました")
    w = w / nw
    g = (u + 1j * w) / math.sqrt(2.0)
    z = v + delta * g
    z = z / np.linalg.norm(z)
    return z, "zero_closure_parent_complement_pair"


def make_initial_state(n, delta, seed, tol, parent_residual_max):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(CORE_SEED_BASE + 1000 * n + seed)
    v_raw, residual_raw, _ = make_parent(
        sys_lr, rng, iters=1200, tol=tol, restarts=8
    )
    p, q, v = parent_plane(v_raw)
    residual = eigenmode_residual(sys_lr, v)
    if residual > parent_residual_max:
        raise RuntimeError(
            f"親構成残差 {residual:.3e} が上限 "
            f"{parent_residual_max:.3e} を超えました"
        )

    sys_lr.set_theta(np.angle(v))
    if n == 3:
        z, seed_kind = n3_exact_zero_closure_state(v, p, q, rng, delta)
    else:
        z, seed_kind = parent_complement_zero_closure_state(
            v, p, q, rng, delta
        )

    return {
        "system": sys_lr,
        "rng": rng,
        "z": z,
        "parent": v,
        "p0": p,
        "q0": q,
        "parent_residual_raw": residual_raw,
        "parent_residual": residual,
        "seed_kind": seed_kind,
    }


def participation_ratio(z):
    weights = np.abs(z) ** 2
    return float(np.sum(weights) ** 2 / np.sum(weights * weights))


def normalized_entropy(z):
    weights = np.abs(z) ** 2
    weights = weights / np.sum(weights)
    nz = weights[weights > 0.0]
    if len(weights) <= 1:
        return 0.0
    return float(-np.sum(nz * np.log(nz)) / np.log(len(weights)))


def dense_generator(sys_lr):
    """現在 set_theta 済みの低ランク表現から小N用密行列を再構成する。"""
    eye = np.eye(sys_lr.m)
    return np.column_stack(
        [sys_lr.kmatvec(eye[:, j]) for j in range(sys_lr.m)]
    )


def instantaneous_planes(sys_lr, z, max_planes):
    """瞬時生成子の全 sigma と各実回転平面への状態エネルギーを返す。"""
    k_dense = dense_generator(sys_lr)
    eigenvalues, eigenvectors = np.linalg.eigh(1j * k_dense)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    indices = np.where(eigenvalues > 1e-11 * scale)[0][::-1]
    sigmas = []
    energies = []
    for idx in indices[:max_planes]:
        u = eigenvectors[:, idx]
        p = math.sqrt(2.0) * u.real
        q = math.sqrt(2.0) * u.imag
        p = p / np.linalg.norm(p)
        q = q - (q @ p) * p
        q = q / np.linalg.norm(q)
        energy = abs(p @ z) ** 2 + abs(q @ z) ** 2
        sigmas.append(float(eigenvalues[idx]))
        energies.append(float(np.real(energy)))
    n_planes = len(sigmas)
    while len(sigmas) < max_planes:
        sigmas.append(0.0)
        energies.append(0.0)
    plane_sum = float(sum(energies))
    norm2 = float(np.real(np.vdot(z, z)))
    kernel_energy = max(0.0, norm2 - plane_sum)
    return sigmas, energies, kernel_energy, n_planes


def pack_complex(z):
    return np.concatenate([z.real, z.imag])


def unpack_complex(x):
    m = len(x) // 2
    return x[:m] + 1j * x[m:]


def exact_map(sys_lr, z):
    """厳密 sigma_1 正規化による原本1ステップ写像。"""
    sys_lr.set_theta(np.angle(z))
    sigmas = sys_lr.sigma_spectrum()
    if len(sigmas) == 0 or sigmas[0] <= 0.0:
        raise RuntimeError("正の sigma_1 が得られません")
    return sys_lr.cayley_step(z, float(sigmas[0]))


def constrained_tangent_basis(z):
    """ノルム、Re/Im(Z^T Z)、共通位相方向を除いた実接空間基底。"""
    x = z.real
    y = z.imag
    rows = np.stack([
        np.concatenate([x, y]),       # Hermitian norm
        np.concatenate([x, -y]),      # Re(Z^T Z)
        np.concatenate([y, x]),       # Im(Z^T Z)
        np.concatenate([-y, x]),      # common phase neutral direction
    ])
    _, singular, vh = np.linalg.svd(rows, full_matrices=True)
    threshold = max(rows.shape) * np.finfo(float).eps * max(singular[0], 1.0)
    rank = int(np.sum(singular > threshold))
    return vh[rank:].T


def local_transverse_multiplier(sys_lr, z, eps):
    """制約接空間上の1ステップ最大特異値を中心差分で測る。"""
    x0 = pack_complex(z)
    dim = len(x0)
    jac = np.empty((dim, dim))
    for j in range(dim):
        dx = np.zeros(dim)
        dx[j] = eps
        zp = unpack_complex(x0 + dx)
        zm = unpack_complex(x0 - dx)
        fp = pack_complex(exact_map(sys_lr, zp))
        fm = pack_complex(exact_map(sys_lr, zm))
        jac[:, j] = (fp - fm) / (2.0 * eps)

    z1 = exact_map(sys_lr, z)
    q0 = constrained_tangent_basis(z)
    q1 = constrained_tangent_basis(z1)
    restricted = q1.T @ jac @ q0
    singular = np.linalg.svd(restricted, compute_uv=False)
    smax = float(singular[0]) if len(singular) else 0.0
    return smax, math.log(smax) if smax > 0.0 else None


def quick_f(z, p0, q0):
    zp = z - p0 * (p0 @ z) - q0 * (q0 @ z)
    norm2 = float(np.real(np.vdot(z, z)))
    return float(np.real(np.vdot(zp, zp))) / norm2


def measure_row(sys_lr, z, p0, q0, tau, max_planes,
                previous_tau, previous_f, jacobian_every, jacobian_eps):
    sys_lr.set_theta(np.angle(z))
    f = quick_f(z, p0, q0)
    a_complement = math.sqrt(max(f, 0.0))
    norm2 = float(np.real(np.vdot(z, z)))
    ztz = complex(z @ z)
    sigmas, energies, kernel_energy, n_planes = instantaneous_planes(
        sys_lr, z, max_planes
    )
    has_second_plane = n_planes >= 2
    sigma_ratio = (
        sigmas[1] / sigmas[0]
        if has_second_plane and sigmas[0] > 0.0
        else ""
    )
    epsilon_half = 0.5 - sigma_ratio if has_second_plane else ""
    a_plane2 = (
        math.sqrt(max(energies[1], 0.0))
        if has_second_plane
        else ""
    )
    abs_z = np.abs(z)
    log_growth = ""
    if (
        previous_tau is not None
        and previous_f is not None
        and f > 0.0
        and previous_f > 0.0
    ):
        log_growth = (
            math.log(f) - math.log(previous_f)
        ) / (tau - previous_tau)

    jac_smax = ""
    jac_lambda = ""
    if jacobian_every > 0 and tau % jacobian_every == 0:
        jac_smax, jac_lambda = local_transverse_multiplier(
            sys_lr, z, jacobian_eps
        )

    row = {
        "tau": tau,
        "f_initial_plane": f,
        "a_complement": a_complement,
        "log_growth_f_per_step": log_growth,
        "norm2": norm2,
        "abs_ztz": abs(ztz),
        "pr": participation_ratio(z),
        "entropy_normalized": normalized_entropy(z),
        "relation_abs_mean": float(np.mean(abs_z)),
        "relation_abs_median": float(np.median(abs_z)),
        "relation_abs_q05": float(np.quantile(abs_z, 0.05)),
        "relation_abs_q95": float(np.quantile(abs_z, 0.95)),
        "n_instantaneous_planes": n_planes,
        "sigma_ratio_2_1": sigma_ratio,
        "epsilon_half": epsilon_half,
        "a_plane2": a_plane2,
        "kernel_energy": kernel_energy,
        "jacobian_smax_transverse": jac_smax,
        "jacobian_lambda_transverse": jac_lambda,
    }
    for i in range(max_planes):
        row[f"sigma_{i + 1}"] = sigmas[i]
        row[f"h_plane_{i + 1}"] = energies[i]
    return row


def scalar_stats(rows, key):
    pairs = [
        (float(row["tau"]), float(row[key]))
        for row in rows
        if row[key] not in ("", None)
    ]
    if not pairs:
        return {
            "n": 0,
            "mean": None,
            "std": None,
            "q05": None,
            "median": None,
            "q95": None,
            "oscillation_half_q90": None,
            "slope_per_step": None,
        }
    taus = np.array([pair[0] for pair in pairs], dtype=float)
    values = np.array([pair[1] for pair in pairs], dtype=float)
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
        "oscillation_half_q90": float(
            0.5 * (np.quantile(values, 0.95) - np.quantile(values, 0.05))
        ),
        "slope_per_step": slope,
    }


def tail_summary(rows, tail_window):
    last_tau = int(rows[-1]["tau"])
    tail = [row for row in rows if row["tau"] >= last_tau - tail_window]
    keys = (
        "f_initial_plane",
        "a_complement",
        "a_plane2",
        "sigma_ratio_2_1",
        "epsilon_half",
        "pr",
        "entropy_normalized",
        "relation_abs_median",
        "kernel_energy",
    )
    result = {
        "start_tau": int(tail[0]["tau"]),
        "end_tau": last_tau,
        "n_samples": len(tail),
    }
    for key in keys:
        result[key] = scalar_stats(tail, key)

    half = len(tail) // 2
    if half >= 2:
        first = tail[:half]
        second = tail[half:]
        result["adjacent_window_mean_change"] = {
            key: mean_change(first, second, key)
            for key in ("a_complement", "a_plane2", "epsilon_half")
        }
    return result


def mean_change(first, second, key):
    mean_first = scalar_stats(first, key)["mean"]
    mean_second = scalar_stats(second, key)["mean"]
    if mean_first is None or mean_second is None:
        return None
    return mean_second - mean_first


def run(args):
    if args.n < 3:
        raise ValueError("Nは3以上でなければなりません")
    if args.record_every < 1:
        raise ValueError("record_everyは1以上でなければなりません")
    if args.cap < args.record_every:
        raise ValueError("capはrecord_every以上でなければなりません")

    initial = make_initial_state(
        args.n,
        args.delta,
        args.seed,
        args.tol,
        args.parent_residual_max,
    )
    sys_lr = initial["system"]
    z = initial["z"]
    p0 = initial["p0"]
    q0 = initial["q0"]
    c0 = complex(z @ z)
    norm0 = float(np.real(np.vdot(z, z)))
    f0 = quick_f(z, p0, q0)
    progress(
        f"low-N N={args.n} 親残差={initial['parent_residual']:.3e} "
        f"|Z^T Z|={abs(c0):.3e} f0={f0:.3e}"
    )

    rows = []
    crossing_tau = None
    max_closure_deviation = 0.0
    max_norm_deviation = 0.0
    previous_tau = None
    previous_f = None
    start = time.time()
    for tau in range(args.cap + 1):
        f = quick_f(z, p0, q0)
        if crossing_tau is None and f > args.crossing:
            crossing_tau = tau
            progress(f"low-N N={args.n} 閾値交差 tau={tau}")

        closure_deviation = abs(complex(z @ z) - c0)
        norm_deviation = abs(float(np.real(np.vdot(z, z))) - norm0)
        max_closure_deviation = max(max_closure_deviation, closure_deviation)
        max_norm_deviation = max(max_norm_deviation, norm_deviation)

        if tau % args.record_every == 0 or tau == args.cap:
            row = measure_row(
                sys_lr,
                z,
                p0,
                q0,
                tau,
                args.n,
                previous_tau,
                previous_f,
                args.jacobian_every,
                args.jacobian_eps,
            )
            rows.append(row)
            previous_tau = tau
            previous_f = row["f_initial_plane"]

        if tau < args.cap:
            z = exact_map(sys_lr, z)
        if tau > 0 and tau % args.progress_every == 0:
            progress(
                f"low-N N={args.n} tau={tau} f={f:.4e} "
                f"経過={time.time() - start:.1f}s"
            )

    runtime = time.time() - start
    os.makedirs(args.result_dir, exist_ok=True)
    tag = f"N{args.n:05d}_delta{args.delta:.0e}_seed{args.seed:03d}"
    csv_path = os.path.join(args.result_dir, f"trajectory_{tag}.csv")
    json_path = os.path.join(args.result_dir, f"summary_{tag}.json")

    fieldnames = list(rows[0].keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    tail = tail_summary(rows, args.tail_window)
    relation_amplitude_measured = tail["relation_abs_median"]["mean"]
    relation_amplitude_equal = 1.0 / math.sqrt(sys_lr.m)
    summary = {
        "experiment": "lowN_metastable_scan_v1",
        "n": args.n,
        "m": sys_lr.m,
        "delta": args.delta,
        "seed": args.seed,
        "cap": args.cap,
        "record_every": args.record_every,
        "tol": args.tol,
        "parent_residual_max": args.parent_residual_max,
        "parent_residual_raw": initial["parent_residual_raw"],
        "parent_residual": initial["parent_residual"],
        "seed_kind": initial["seed_kind"],
        "abs_ztz_initial": abs(c0),
        "f_initial": f0,
        "crossing_threshold": args.crossing,
        "crossing_tau": crossing_tau,
        "max_closure_deviation": max_closure_deviation,
        "max_norm_deviation": max_norm_deviation,
        "jacobian_every": args.jacobian_every,
        "jacobian_eps": args.jacobian_eps,
        "tail_window": args.tail_window,
        "tail": tail,
        "relation_amplitude_duality": {
            "definition": "tail mean of median_e |Z_e|",
            "measured": relation_amplitude_measured,
            "equal_amplitude_prediction_1_over_sqrt_m":
                relation_amplitude_equal,
            "measured_times_sqrt_m":
                relation_amplitude_measured * math.sqrt(sys_lr.m),
            "asymptotic_n_times_measured":
                args.n * relation_amplitude_measured,
        },
        "runtime_sec": runtime,
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
    parser.add_argument("--cap", type=int, default=100000)
    parser.add_argument("--record-every", type=int, default=5)
    parser.add_argument("--tail-window", type=int, default=20000)
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--parent-residual-max", type=float, default=1e-10)
    parser.add_argument("--crossing", type=float, default=0.05)
    parser.add_argument("--jacobian-every", type=int, default=500)
    parser.add_argument("--jacobian-eps", type=float, default=1e-7)
    parser.add_argument("--progress-every", type=int, default=5000)
    parser.add_argument("--result-dir", default=RESULT_DIR)
    return parser.parse_args()


def main():
    run(parse_args())


if __name__ == "__main__":
    main()
