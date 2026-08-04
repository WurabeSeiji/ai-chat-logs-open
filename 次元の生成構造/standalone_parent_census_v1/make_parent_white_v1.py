#!/usr/bin/env python3
"""Standalone white-null parent generator.

Theory choice v1:
    resolution enters only through gamma_N = tan(pi / resolution).
    No harmonic count, finite-order projection, or particle label is supplied.

The program intentionally imports no research module from the repository.  Its
output is the versioned ``closed_wave_trajectory_v1`` file contract described
in README.md next to this file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np


SCHEMA_VERSION = "closed_wave_trajectory_v1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def build_edges(n_body: int) -> tuple[np.ndarray, np.ndarray]:
    edge_a, edge_b = np.triu_indices(n_body, k=1)
    return edge_a.astype(np.int64), edge_b.astype(np.int64)


def participation_ratio(z: np.ndarray) -> float:
    power = np.abs(z) ** 2
    denom = float(np.sum(power * power))
    if denom == 0.0:
        return 0.0
    return float(np.sum(power) ** 2 / denom)


class LowRankClosedWaveSystem:
    """Low-rank form of the real antisymmetric relation-wave generator."""

    def __init__(self, n_body: int, resolution: int):
        self.n_body = n_body
        self.resolution = resolution
        self.gamma = math.tan(math.pi / resolution)
        self.edge_a, self.edge_b = build_edges(n_body)
        self.m_relations = len(self.edge_a)

        self.J = np.zeros((2 * n_body, 2 * n_body), dtype=float)
        self.J[:n_body, n_body:] = np.eye(n_body)
        self.J[n_body:, :n_body] = -np.eye(n_body)
        self.G = np.zeros_like(self.J)
        self.cos_theta = np.empty(self.m_relations)
        self.sin_theta = np.empty(self.m_relations)

    def vertex_sum(self, values: np.ndarray) -> np.ndarray:
        n = self.n_body
        if np.iscomplexobj(values):
            real = (
                np.bincount(self.edge_a, weights=values.real, minlength=n)
                + np.bincount(self.edge_b, weights=values.real, minlength=n)
            )
            imag = (
                np.bincount(self.edge_a, weights=values.imag, minlength=n)
                + np.bincount(self.edge_b, weights=values.imag, minlength=n)
            )
            return real + 1j * imag
        return (
            np.bincount(self.edge_a, weights=values, minlength=n)
            + np.bincount(self.edge_b, weights=values, minlength=n)
        )

    def set_theta(self, theta: np.ndarray) -> None:
        n = self.n_body
        self.cos_theta = np.cos(theta)
        self.sin_theta = np.sin(theta)

        theta_matrix = np.zeros((n, n), dtype=float)
        theta_matrix[self.edge_a, self.edge_b] = theta
        theta_matrix[self.edge_b, self.edge_a] = theta
        ct = np.cos(theta_matrix)
        st = np.sin(theta_matrix)
        np.fill_diagonal(ct, 0.0)
        np.fill_diagonal(st, 0.0)

        gcc = ct * ct
        gcs = ct * st
        gss = st * st
        np.fill_diagonal(gcc, gcc.sum(axis=1))
        np.fill_diagonal(gcs, gcs.sum(axis=1))
        np.fill_diagonal(gss, gss.sum(axis=1))

        self.G[:n, :n] = gcc
        self.G[:n, n:] = gcs
        self.G[n:, :n] = gcs
        self.G[n:, n:] = gss

    def wt(self, z: np.ndarray) -> np.ndarray:
        return np.concatenate(
            [
                self.vertex_sum(self.cos_theta * z),
                self.vertex_sum(self.sin_theta * z),
            ]
        )

    def w(self, y: np.ndarray) -> np.ndarray:
        n = self.n_body
        yc, ys = y[:n], y[n:]
        return self.cos_theta * (yc[self.edge_a] + yc[self.edge_b]) + self.sin_theta * (
            ys[self.edge_a] + ys[self.edge_b]
        )

    def kmatvec(self, z: np.ndarray) -> np.ndarray:
        vs = self.vertex_sum(self.sin_theta * z)
        vc = self.vertex_sum(self.cos_theta * z)
        return self.cos_theta * (vs[self.edge_a] + vs[self.edge_b]) - self.sin_theta * (
            vc[self.edge_a] + vc[self.edge_b]
        )

    def sigma_spectrum(self) -> np.ndarray:
        eigenvalues = np.linalg.eigvals(self.J @ self.G)
        positive = eigenvalues.imag[eigenvalues.imag > 1e-12]
        return np.sort(positive)[::-1]

    def sigma_max_power(
        self,
        warm_vector: np.ndarray,
        tolerance: float,
        max_iterations: int,
    ) -> tuple[float, np.ndarray, int, float]:
        """Largest singular value of K from power iteration on -K^2.

        The residual is dimensionless:
        ||-K^2 w - sigma^2 w|| / sigma^2.
        """

        w = np.asarray(warm_vector, dtype=float).copy()
        norm_w = float(np.linalg.norm(w))
        if norm_w == 0.0 or not np.isfinite(norm_w):
            w = np.ones(self.m_relations, dtype=float)
            norm_w = float(np.linalg.norm(w))
        w /= norm_w

        residual = math.inf
        sigma = 0.0
        for iteration in range(1, max_iterations + 1):
            minus_k2_w = -self.kmatvec(self.kmatvec(w))
            norm_k2 = float(np.linalg.norm(minus_k2_w))
            if norm_k2 == 0.0 or not np.isfinite(norm_k2):
                return 0.0, w, iteration, math.inf
            w = minus_k2_w / norm_k2
            kw = self.kmatvec(w)
            sigma = float(np.linalg.norm(kw))
            if sigma == 0.0 or not np.isfinite(sigma):
                return 0.0, w, iteration, math.inf
            check = -self.kmatvec(kw)
            residual = float(np.linalg.norm(check - sigma * sigma * w) / (sigma * sigma))
            if residual <= tolerance:
                return sigma, w, iteration, residual
        return sigma, w, max_iterations, residual

    def cayley_step(self, z: np.ndarray, sigma_max: float) -> np.ndarray:
        if sigma_max <= 0.0 or not np.isfinite(sigma_max):
            raise FloatingPointError(f"invalid sigma_max: {sigma_max}")
        gamma_over_sigma = self.gamma / sigma_max
        right = z + gamma_over_sigma * self.kmatvec(z)
        reduced = (sigma_max / self.gamma) * self.J + self.G
        lifted = np.linalg.solve(reduced, self.wt(right))
        return right - self.w(lifted)


def make_white_null_state(rng: np.random.Generator, m_relations: int) -> np.ndarray:
    """Haar-random oriented real 2-plane represented as one complex null vector."""

    gaussian = rng.normal(size=(m_relations, 2))
    q, r = np.linalg.qr(gaussian, mode="reduced")
    # Fix the QR sign convention so a seed has a stable, explicit orientation.
    signs = np.where(np.diag(r) < 0.0, -1.0, 1.0)
    q = q * signs[np.newaxis, :]
    return (q[:, 0] + 1j * q[:, 1]) / math.sqrt(2.0)


def eigenmode_residual(system: LowRankClosedWaveSystem, vector: np.ndarray) -> tuple[float, float]:
    kv = system.kmatvec(vector)
    mu = float(np.real(np.vdot(vector, 1j * kv)))
    residual = float(np.linalg.norm(1j * kv - mu * vector) / max(abs(mu), 1.0))
    return mu, residual


def make_single_mode_parent(
    system: LowRankClosedWaveSystem,
    rng: np.random.Generator,
    iterations: int,
    beta: float,
    tolerance: float,
    restarts: int,
) -> tuple[np.ndarray, dict[str, object]]:
    """Original single self-consistent circular-mode construction, as a control."""

    best_vector: np.ndarray | None = None
    best_residual = math.inf
    best_mu = math.nan
    best_iteration = 0
    best_restart = 0

    for restart in range(1, restarts + 1):
        theta = rng.uniform(0.0, 2.0 * math.pi, system.m_relations)
        vector: np.ndarray | None = None
        for iteration in range(1, iterations + 1):
            system.set_theta(theta)
            eigenvalues, eigenvectors = np.linalg.eig(system.J @ system.G)
            index = int(np.argmin(eigenvalues.imag))
            vector = system.w(eigenvectors[:, index].astype(complex))
            vector_norm = float(np.linalg.norm(vector))
            if vector_norm == 0.0 or not np.isfinite(vector_norm):
                break
            vector /= vector_norm
            proposed = np.angle(vector)
            mixed = (1.0 - beta) * np.exp(1j * theta) + beta * np.exp(1j * proposed)
            theta = np.angle(mixed)
            if iteration % 10 == 0:
                system.set_theta(np.angle(vector))
                mu, residual = eigenmode_residual(system, vector)
                if residual <= tolerance:
                    break
        if vector is None:
            continue
        system.set_theta(np.angle(vector))
        mu, residual = eigenmode_residual(system, vector)
        if residual < best_residual:
            best_vector = vector.copy()
            best_residual = residual
            best_mu = mu
            best_iteration = iteration
            best_restart = restart
        if residual <= tolerance:
            break

    if best_vector is None:
        raise RuntimeError("single-mode parent construction produced no finite vector")
    spectrum = system.sigma_spectrum()
    return best_vector, {
        "construction": "single_self_consistent_circular_mode_control",
        "residual": best_residual,
        "mu": best_mu,
        "iterations_used": best_iteration,
        "restart_used": best_restart,
        "sigma_spectrum": spectrum.tolist(),
    }


def state_metrics(z: np.ndarray) -> tuple[float, float, float]:
    norm = float(np.real(np.vdot(z, z)))
    closure_abs = abs(complex(z @ z))
    return norm, closure_abs, participation_ratio(z)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate one closed N-body relation-wave trajectory from white null noise."
    )
    parser.add_argument("--n-body", type=int, required=True, help="N-body model size N (N>=3)")
    parser.add_argument(
        "--resolution", type=int, required=True, help="system resolution mathcal N (>=4)"
    )
    parser.add_argument("--seed", type=int, required=True, help="reproducibility seed")
    parser.add_argument("--steps", type=int, required=True, help="fixed numerical observation steps")
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    parser.add_argument(
        "--source",
        choices=("white_null", "single_mode"),
        default="white_null",
        help="white-null main condition or original single-mode control",
    )
    parser.add_argument("--record-stride", type=int, default=1)
    parser.add_argument(
        "--storage-dtype", choices=("complex64", "complex128"), default="complex128"
    )
    parser.add_argument("--sigma-tolerance", type=float, default=1e-12)
    parser.add_argument("--sigma-max-iterations", type=int, default=100)
    parser.add_argument("--parent-iterations", type=int, default=1200)
    parser.add_argument("--parent-beta", type=float, default=0.5)
    parser.add_argument("--parent-tolerance", type=float, default=1e-10)
    parser.add_argument("--parent-restarts", type=int, default=5)
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.n_body < 3:
        raise SystemExit("--n-body must be >= 3 so the complex null cone is nontrivial")
    if args.resolution < 4:
        raise SystemExit("--resolution must be >= 4")
    if args.steps < 1:
        raise SystemExit("--steps must be >= 1")
    if args.record_stride < 1:
        raise SystemExit("--record-stride must be >= 1")
    if args.sigma_tolerance <= 0.0 or args.sigma_max_iterations < 1:
        raise SystemExit("invalid sigma power-iteration controls")
    if args.output.exists():
        raise SystemExit(f"output already exists; refusing to overwrite: {args.output}")


def main() -> None:
    args = parse_args()
    validate_args(args)
    started = time.time()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    system = LowRankClosedWaveSystem(args.n_body, args.resolution)
    rng = np.random.default_rng(args.seed)

    if args.source == "white_null":
        z = make_white_null_state(rng, system.m_relations)
        source_info: dict[str, object] = {
            "construction": "gaussian_Mx2_QR_then_(q1+i*q2)/sqrt(2)",
            "preassigned_harmonics": 0,
            "finite_order_projection": False,
        }
    else:
        z, source_info = make_single_mode_parent(
            system,
            rng,
            iterations=args.parent_iterations,
            beta=args.parent_beta,
            tolerance=args.parent_tolerance,
            restarts=args.parent_restarts,
        )

    z = np.asarray(z, dtype=np.complex128)
    initial_norm, initial_closure, _ = state_metrics(z)

    record_steps = np.arange(0, args.steps + 1, args.record_stride, dtype=np.int64)
    if record_steps[-1] != args.steps:
        record_steps = np.append(record_steps, args.steps)
    np.save(output / "steps.npy", record_steps)
    np.save(output / "edges.npy", np.column_stack([system.edge_a, system.edge_b]))

    storage_dtype = np.complex64 if args.storage_dtype == "complex64" else np.complex128
    trajectory = np.lib.format.open_memmap(
        output / "trajectory.npy",
        mode="w+",
        dtype=storage_dtype,
        shape=(len(record_steps), system.m_relations),
    )
    norms = np.empty(len(record_steps), dtype=float)
    closures = np.empty(len(record_steps), dtype=float)
    prs = np.empty(len(record_steps), dtype=float)
    sigmas = np.empty(len(record_steps), dtype=float)
    sigma_residuals = np.empty(len(record_steps), dtype=float)
    sigma_iterations = np.empty(len(record_steps), dtype=np.int64)

    warm = rng.normal(size=system.m_relations)
    record_index = 0
    progress_every = max(1, args.steps // 20)

    for step in range(args.steps + 1):
        system.set_theta(np.angle(z))
        sigma, warm, n_iter, sigma_residual = system.sigma_max_power(
            warm,
            tolerance=args.sigma_tolerance,
            max_iterations=args.sigma_max_iterations,
        )

        if record_index < len(record_steps) and step == int(record_steps[record_index]):
            norm, closure_abs, pr = state_metrics(z)
            trajectory[record_index] = z.astype(storage_dtype, copy=False)
            norms[record_index] = norm
            closures[record_index] = closure_abs
            prs[record_index] = pr
            sigmas[record_index] = sigma
            sigma_residuals[record_index] = sigma_residual
            sigma_iterations[record_index] = n_iter
            record_index += 1

        if step < args.steps:
            z = system.cayley_step(z, sigma)
            if not np.all(np.isfinite(z)):
                raise FloatingPointError(f"non-finite state at step {step + 1}")

        if step % progress_every == 0 or step == args.steps:
            print(
                f"[make_parent] step={step}/{args.steps} "
                f"|Z^T Z|={abs(complex(z @ z)):.3e} sigma_res={sigma_residual:.3e}",
                file=sys.stderr,
                flush=True,
            )

    trajectory.flush()
    del trajectory
    np.savez(
        output / "diagnostics.npz",
        steps=record_steps,
        norm=norms,
        closure_abs=closures,
        participation_ratio=prs,
        sigma_max=sigmas,
        sigma_power_residual=sigma_residuals,
        sigma_power_iterations=sigma_iterations,
    )

    script_path = Path(__file__).resolve()
    manifest = {
        "schema": SCHEMA_VERSION,
        "created_at_unix": time.time(),
        "generator": script_path.name,
        "generator_sha256": sha256_file(script_path),
        "source": args.source,
        "source_info": source_info,
        "theory_inputs": {
            "resolution_mathcal_N": args.resolution,
            "lambda0": 2.0 * math.pi / args.resolution,
            "gamma": system.gamma,
            "n_body": args.n_body,
            "m_relations": system.m_relations,
        },
        "reproducibility": {"seed": args.seed},
        "numerical_observation": {
            "steps": args.steps,
            "record_stride": args.record_stride,
            "storage_dtype": args.storage_dtype,
            "sigma_tolerance": args.sigma_tolerance,
            "sigma_max_iterations": args.sigma_max_iterations,
        },
        "explicit_absences": {
            "harmonic_count_H": None,
            "preassigned_harmonic_weights": False,
            "finite_order_projection": False,
            "U_resolution_equals_I_enforced": False,
            "particle_labels_in_generator": False,
            "particle_completion_stop_rule": False,
        },
        "arrays": {
            "edges": {"file": "edges.npy", "shape": [system.m_relations, 2]},
            "steps": {"file": "steps.npy", "shape": [len(record_steps)]},
            "trajectory": {
                "file": "trajectory.npy",
                "shape": [len(record_steps), system.m_relations],
                "dtype": args.storage_dtype,
            },
            "diagnostics": {"file": "diagnostics.npz"},
        },
        "invariant_audit": {
            "initial_norm": initial_norm,
            "initial_closure_abs": initial_closure,
            "max_norm_drift": float(np.max(np.abs(norms - norms[0]))),
            "max_closure_abs": float(np.max(closures)),
            "max_sigma_power_residual": float(np.max(sigma_residuals)),
        },
        "runtime_seconds": time.time() - started,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(json.dumps(manifest["invariant_audit"], ensure_ascii=False, indent=2))
    print(f"saved: {output}")


if __name__ == "__main__":
    main()
