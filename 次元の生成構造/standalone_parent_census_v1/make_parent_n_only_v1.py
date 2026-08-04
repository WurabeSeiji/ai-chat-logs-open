#!/usr/bin/env python3
"""唯一の理論入力Nから、白色雑音起源の多周波数親状態を作る。

公開関数:
    make_parent(N, seed=None, max_retries=3) -> ParentResult

seed指定時はその系列で1回だけ試し、未収束なら例外で終了する。
seed未指定時はOSから128 bit seedを取得し、初回+max_retries回まで試す。
未収束の最良値を成功として返さない。

周波数候補数を別引数では与えない。分解能Nの各セル n=1..N に対して
白色複素Gaussian係数を生成し、既存make_parentと同じ自己無撞着反復を行う。
各列は一つの全M関係成分にまたがる零閉塞波である。
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import secrets
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "n_only_multifrequency_parent_v1"
SOLVER_ITERATIONS = 2000
SOLVER_BETA = 0.5
SOLVER_TOLERANCE = 1e-12
CHECK_EVERY = 10


def build_edges(n: int) -> tuple[np.ndarray, np.ndarray]:
    a, b = np.triu_indices(n, k=1)
    return a.astype(np.int64), b.astype(np.int64)


class LowRankSystem:
    """既存make_parentと同じ実反対称関係波生成子の低rank表現。"""

    def __init__(self, n: int):
        self.n = n
        self.edge_a, self.edge_b = build_edges(n)
        self.m = len(self.edge_a)
        self.J = np.zeros((2 * n, 2 * n), dtype=float)
        self.J[:n, n:] = np.eye(n)
        self.J[n:, :n] = -np.eye(n)
        self.G = np.zeros_like(self.J)
        self.ct = np.empty(self.m)
        self.st = np.empty(self.m)

    def vertex_sum(self, values: np.ndarray) -> np.ndarray:
        if np.iscomplexobj(values):
            real = np.bincount(self.edge_a, weights=values.real, minlength=self.n)
            real += np.bincount(self.edge_b, weights=values.real, minlength=self.n)
            imag = np.bincount(self.edge_a, weights=values.imag, minlength=self.n)
            imag += np.bincount(self.edge_b, weights=values.imag, minlength=self.n)
            return real + 1j * imag
        result = np.bincount(self.edge_a, weights=values, minlength=self.n)
        result += np.bincount(self.edge_b, weights=values, minlength=self.n)
        return result

    def set_theta(self, theta: np.ndarray) -> None:
        self.ct = np.cos(theta)
        self.st = np.sin(theta)
        matrix = np.zeros((self.n, self.n), dtype=float)
        matrix[self.edge_a, self.edge_b] = theta
        matrix[self.edge_b, self.edge_a] = theta
        ct = np.cos(matrix)
        st = np.sin(matrix)
        np.fill_diagonal(ct, 0.0)
        np.fill_diagonal(st, 0.0)
        gcc = ct * ct
        gcs = ct * st
        gss = st * st
        np.fill_diagonal(gcc, gcc.sum(axis=1))
        np.fill_diagonal(gcs, gcs.sum(axis=1))
        np.fill_diagonal(gss, gss.sum(axis=1))
        self.G[: self.n, : self.n] = gcc
        self.G[: self.n, self.n :] = gcs
        self.G[self.n :, : self.n] = gcs
        self.G[self.n :, self.n :] = gss

    def w(self, y: np.ndarray) -> np.ndarray:
        yc, ys = y[: self.n], y[self.n :]
        return self.ct * (yc[self.edge_a] + yc[self.edge_b]) + self.st * (
            ys[self.edge_a] + ys[self.edge_b]
        )

    def kmatvec(self, z: np.ndarray) -> np.ndarray:
        vs = self.vertex_sum(self.st * z)
        vc = self.vertex_sum(self.ct * z)
        return self.ct * (vs[self.edge_a] + vs[self.edge_b]) - self.st * (
            vc[self.edge_a] + vc[self.edge_b]
        )

    def sigma_spectrum(self) -> np.ndarray:
        eigenvalues = np.linalg.eigvals(self.J @ self.G)
        positive = eigenvalues.imag[eigenvalues.imag > 1e-12]
        return np.sort(positive)[::-1]


def eigenmode_residual(system: LowRankSystem, vector: np.ndarray) -> tuple[float, float]:
    kv = system.kmatvec(vector)
    mu = float(np.real(np.vdot(vector, 1j * kv)))
    residual = float(np.linalg.norm(1j * kv - mu * vector))
    return mu, residual


@dataclasses.dataclass
class ModeAttempt:
    frequency_multiple: int
    converged: bool
    iterations: int
    residual: float
    mu: float
    sigma_max: float
    closure_abs: float
    failure_reason: str | None
    vector: np.ndarray | None = dataclasses.field(repr=False)

    def audit_dict(self) -> dict[str, Any]:
        return {
            "frequency_multiple": self.frequency_multiple,
            "converged": self.converged,
            "iterations": self.iterations,
            "residual": self.residual,
            "mu": self.mu,
            "sigma_max": self.sigma_max,
            "closure_abs": self.closure_abs,
            "failure_reason": self.failure_reason,
        }


@dataclasses.dataclass
class SeedAttempt:
    attempt_number: int
    seed: int
    success: bool
    failure_reason: str | None
    modes: list[ModeAttempt]

    def audit_dict(self) -> dict[str, Any]:
        return {
            "attempt_number": self.attempt_number,
            "seed": str(self.seed),
            "success": self.success,
            "failure_reason": self.failure_reason,
            "modes": [mode.audit_dict() for mode in self.modes],
        }


@dataclasses.dataclass
class ParentResult:
    n: int
    m: int
    seed: int
    seed_was_explicit: bool
    white_input: np.ndarray
    parent_modes: np.ndarray
    edges: np.ndarray
    mode_weights: np.ndarray
    attempts: list[SeedAttempt]
    runtime_seconds: float


class ParentConstructionError(RuntimeError):
    def __init__(self, message: str, attempts: list[SeedAttempt]):
        super().__init__(message)
        self.attempts = attempts


def solve_self_consistent_mode(
    system: LowRankSystem,
    theta_initial: np.ndarray,
    frequency_multiple: int,
) -> ModeAttempt:
    """既存make_parentのargmin固有枝を、一つの白色初期位相から1回解く。"""

    theta = np.asarray(theta_initial, dtype=float).copy()
    vector: np.ndarray | None = None
    mu = math.nan
    residual = math.inf
    iteration = 0
    for iteration in range(1, SOLVER_ITERATIONS + 1):
        system.set_theta(theta)
        try:
            eigenvalues, eigenvectors = np.linalg.eig(system.J @ system.G)
        except np.linalg.LinAlgError:
            return ModeAttempt(
                frequency_multiple,
                False,
                iteration,
                math.inf,
                math.nan,
                math.nan,
                math.inf,
                "eigendecomposition_failed",
                None,
            )
        index = int(np.argmin(eigenvalues.imag))
        vector = system.w(eigenvectors[:, index].astype(complex))
        norm = float(np.linalg.norm(vector))
        if norm == 0.0 or not np.isfinite(norm):
            vector = None
            break
        vector /= norm
        proposed = np.angle(vector)
        mixed = (1.0 - SOLVER_BETA) * np.exp(1j * theta)
        mixed += SOLVER_BETA * np.exp(1j * proposed)
        if np.any(np.abs(mixed) == 0.0):
            vector = None
            break
        theta = np.angle(mixed)
        if iteration % CHECK_EVERY == 0:
            system.set_theta(np.angle(vector))
            mu, residual = eigenmode_residual(system, vector)
            if residual < SOLVER_TOLERANCE:
                break

    if vector is None:
        return ModeAttempt(
            frequency_multiple,
            False,
            iteration,
            math.inf,
            math.nan,
            math.nan,
            math.inf,
            "nonfinite_or_zero_eigenvector",
            None,
        )
    system.set_theta(np.angle(vector))
    mu, residual = eigenmode_residual(system, vector)
    sigmas = system.sigma_spectrum()
    sigma_max = float(sigmas[0]) if len(sigmas) else 0.0
    closure = float(abs(complex(vector @ vector)))
    converged = bool(
        np.isfinite(residual)
        and residual < SOLVER_TOLERANCE
        and closure < 100.0 * np.finfo(float).eps * system.m
    )
    return ModeAttempt(
        frequency_multiple,
        converged,
        iteration,
        residual,
        mu,
        sigma_max,
        closure,
        None if converged else "self_consistency_tolerance_not_reached",
        vector.copy(),
    )


def _attempt_seed(n: int, seed: int, attempt_number: int) -> tuple[SeedAttempt, Any]:
    m = n * (n - 1) // 2
    rng = np.random.default_rng(seed)
    system = LowRankSystem(n)
    white = np.empty((m, n), dtype=np.complex128)
    modes: list[ModeAttempt] = []
    for index in range(n):
        # CN(0,1) の極座標表示。位相を先に消費することで、既存make_parentの
        # rng.uniform(0,2π,M) 初期化契約を各周波数候補でそのまま保つ。
        theta = rng.uniform(0.0, 2.0 * math.pi, m)
        radius = np.sqrt(rng.exponential(scale=1.0, size=m))
        white[:, index] = radius * np.exp(1j * theta)
        mode = solve_self_consistent_mode(system, theta, index + 1)
        modes.append(mode)
        if not mode.converged:
            attempt = SeedAttempt(
                attempt_number,
                seed,
                False,
                f"frequency_multiple_{index + 1}_did_not_converge",
                modes,
            )
            return attempt, None
    column_norms = np.linalg.norm(white, axis=0)
    if np.any(column_norms == 0.0) or not np.all(np.isfinite(column_norms)):
        attempt = SeedAttempt(attempt_number, seed, False, "white_input_nonfinite", modes)
        return attempt, None
    weights = column_norms / np.linalg.norm(column_norms)
    parent = np.column_stack([mode.vector for mode in modes]) * weights[np.newaxis, :]
    attempt = SeedAttempt(attempt_number, seed, True, None, modes)
    return attempt, (white, parent, weights, system)


def make_parent(N: int, seed: int | None = None, max_retries: int = 3) -> ParentResult:
    """白色雑音からN個の位相分解された自己無撞着零閉塞波を構成する。

    N以外に理論パラメータは受け取らない。seedとmax_retriesは実験の
    再現・失敗処理条件である。seed明示時はmax_retriesを使用しない。
    """

    if isinstance(N, bool) or not isinstance(N, (int, np.integer)):
        raise TypeError("N must be an integer")
    if N < 3:
        raise ValueError("N must be >= 3")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, (int, np.integer))):
        raise TypeError("seed must be an integer or None")
    if seed is not None and seed < 0:
        raise ValueError("seed must be >= 0")
    if isinstance(max_retries, bool) or not isinstance(max_retries, (int, np.integer)):
        raise TypeError("max_retries must be an integer")
    if max_retries < 0:
        raise ValueError("max_retries must be >= 0")

    started = time.time()
    seed_was_explicit = seed is not None
    maximum_attempts = 1 if seed_was_explicit else 1 + max_retries
    attempts: list[SeedAttempt] = []
    for attempt_number in range(1, maximum_attempts + 1):
        current_seed = int(seed) if seed_was_explicit else secrets.randbits(128)
        attempt, payload = _attempt_seed(N, current_seed, attempt_number)
        attempts.append(attempt)
        if attempt.success:
            white, parent, weights, system = payload
            return ParentResult(
                n=N,
                m=system.m,
                seed=current_seed,
                seed_was_explicit=seed_was_explicit,
                white_input=white,
                parent_modes=parent,
                edges=np.column_stack([system.edge_a, system.edge_b]),
                mode_weights=weights,
                attempts=attempts,
                runtime_seconds=time.time() - started,
            )
        if seed_was_explicit:
            break
    raise ParentConstructionError(
        f"make_parent aborted after {len(attempts)} unsuccessful attempt(s)", attempts
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def base_manifest(N: int, seed: int | None, max_retries: int) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "generator": Path(__file__).name,
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "numpy_version": np.__version__,
        "prng": "numpy.random.PCG64 via default_rng",
        "os_seed_source": "secrets.randbits(128)",
        "function_contract": {
            "N": N,
            "seed_argument": None if seed is None else str(seed),
            "max_retries_argument": max_retries,
            "explicit_seed_means_single_attempt": True,
            "os_seed_total_attempt_limit": 1 + max_retries,
        },
        "derived": {
            "M": N * (N - 1) // 2,
            "lambda0": 2.0 * math.pi / N,
            "frequency_multiples": list(range(1, N + 1)),
        },
        "fixed_solver_contract": {
            "iterations": SOLVER_ITERATIONS,
            "beta": SOLVER_BETA,
            "tolerance": SOLVER_TOLERANCE,
            "check_every": CHECK_EVERY,
            "branch": "argmin(Im eigenvalue of JG)",
            "unconverged_best_value_is_never_returned": True,
        },
    }


def write_success(output: Path, result: ParentResult, manifest: dict[str, Any]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    np.save(output / "white_input.npy", result.white_input)
    np.save(output / "parent_modes.npy", result.parent_modes)
    np.save(output / "edges.npy", result.edges)
    np.save(output / "mode_weights.npy", result.mode_weights)
    layer_closures = np.abs(np.sum(result.parent_modes * result.parent_modes, axis=0))
    total_closure = float(abs(complex(np.sum(result.parent_modes * result.parent_modes))))
    manifest.update(
        {
            "status": "success",
            "accepted_seed": str(result.seed),
            "seed_was_explicit": result.seed_was_explicit,
            "attempts": [attempt.audit_dict() for attempt in result.attempts],
            "arrays": {
                "white_input": {"file": "white_input.npy", "shape": list(result.white_input.shape)},
                "parent_modes": {"file": "parent_modes.npy", "shape": list(result.parent_modes.shape)},
                "edges": {"file": "edges.npy", "shape": list(result.edges.shape)},
                "mode_weights": {"file": "mode_weights.npy", "shape": list(result.mode_weights.shape)},
            },
            "invariant_audit": {
                "all_modes_converged": True,
                "max_mode_residual": max(
                    mode.residual for mode in result.attempts[-1].modes
                ),
                "max_mode_closure_abs": float(np.max(layer_closures)),
                "total_closure_abs": total_closure,
                "frobenius_norm": float(np.linalg.norm(result.parent_modes)),
                "all_mode_weights_nonzero": bool(np.all(result.mode_weights > 0.0)),
            },
            "runtime_seconds": result.runtime_seconds,
        }
    )
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def write_failure(output: Path, manifest: dict[str, Any], error: ParentConstructionError) -> None:
    output.mkdir(parents=True, exist_ok=False)
    manifest.update(
        {
            "status": "aborted",
            "failure_reason": str(error),
            "attempts": [attempt.audit_dict() for attempt in error.attempts],
        }
    )
    (output / "failure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nだけを理論入力とする多周波数make_parent")
    parser.add_argument("N", type=int)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"output already exists; refusing to overwrite: {output}")
    manifest = base_manifest(args.N, args.seed, args.max_retries)
    try:
        result = make_parent(args.N, seed=args.seed, max_retries=args.max_retries)
    except ParentConstructionError as error:
        write_failure(output, manifest, error)
        print(str(error), file=sys.stderr)
        print(f"failure audit: {output / 'failure_manifest.json'}", file=sys.stderr)
        raise SystemExit(2)
    write_success(output, result, manifest)
    print(
        json.dumps(
            {
                "status": "success",
                "N": result.n,
                "M": result.m,
                "accepted_seed": str(result.seed),
                "attempt_count": len(result.attempts),
                "parent_modes_shape": list(result.parent_modes.shape),
                "output": str(output),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
