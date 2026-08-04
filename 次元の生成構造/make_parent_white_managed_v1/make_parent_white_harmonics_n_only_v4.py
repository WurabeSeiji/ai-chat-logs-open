#!/usr/bin/env python3
"""Nだけを理論入力とする白色雑音・零閉塞・倍音事後読出し親生成器 v4（純進行波・対積負債ゼロ版）。

v3 からの変更（2026-08-04 木原氏裁定・±k対積負債の解消）:
    定理: 同一関係波内の ±k 均衡（定常対）は零閉塞と両立しない
    （Σ_n(c₊e^{+ikn}+c₋e^{−ikn})² = 2N·c₊c₋）。よって各関係波は
    各 |k| につき一方向のみを占有する純進行波構成でなければならない。
    v4 は各関係波を許容モード（2k≢0 mod N）の符号付き集合の上で直接
    スペクトル構成する:
      - 各 |k|=1..⌊(N−1)/2⌋ について方向 s∈{+,−} を乱数硬貨で選択
      - 選択モードに複素 Gaussian 振幅、逆方向は恒等 0
      - 行閉塞 Σc_k c_{N−k} = 0 が項ごとに恒等成立（QR も射影も不要）
      - 方向選択は関係ごと独立 → 海全体では両方向が統計的に均衡
        （定常性・C対称性は関係横断のアンサンブルで実現）
    DC（k=0）・Nyquist（k=N/2）は v3 と同じく不在。
    乱数消費順: 親ソルバが先（v2/v3 と同一）→ 同一シードで parent_vector は
    bitwise 一致。白色部の raw 配列は構成が異なるため v3 と非互換
    （raw_white_noise には (符号ビット, 複素振幅) を保存する）。


v2 からの変更（2026-08-04 木原氏合意・自己対セクター欠陥の修正）:
    ガウス白色雑音の各列から、QR の前に
      - DC 方向（k=0、全1ベクトル）を射影除去
      - 偶数 N ではナイキスト方向（k=N/2、交代±1ベクトル）も射影除去
    する。根拠は仮定の追加ではなく既存公理の定理:
      (1) 閉塞=存在資格——2k≡0 (mod N) の自己対モードは単独で零閉塞できず存在資格がない
      (2) 分解能の倍音しか存在してはならない——k=0 は倍音ではない
      (3) 関係主義——一定オフセットは差を作らないゲージ
    射影はガウス性を保つため、これは「許容モード上の白色雑音」の厳密な実装である。
    乱数消費順は v2 と同一（同一シードで parent_vector と raw_white_noise は bitwise 一致）。


公開関数:
    make_parent(N, seed=None, max_retries=3) -> ParentResult

定義を固定する。

* N体の完全二体関係数は M=N(N-1)/2。
* M本の各関係波 w_m は N 個の標本 w[m,n] を持つ。
* 各関係波は sum_n w[m,n]^2=0 を満たす。
* 全体系も sum_m sum_n w[m,n]^2=0 を満たす。
* 倍音は生成器で配置しない。保存後の読出し器がN点DFTから読む。

既存 make_parent の自己無撞着親 v in C^M を、各関係波の振幅・位相階層
として保持する。各関係波のN点波形は、N x 2 実Gaussian白色雑音から
QRで得る直交二成分 (q1+i q2)/sqrt(2) とする。したがって閉塞は補正値を
足して作るのではなく、等ノルム直交二成分から恒等的に成立する。
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


SCHEMA = "n_only_white_closed_harmonic_parent_v4_pure_traveling"
SOLVER_ITERATIONS = 2000
SOLVER_BETA = 0.5
SOLVER_TOLERANCE = 1e-12
CHECK_EVERY = 10


def build_edges(n: int) -> tuple[np.ndarray, np.ndarray]:
    a, b = np.triu_indices(n, k=1)
    return a.astype(np.int64), b.astype(np.int64)


class LowRankSystem:
    """既存make_parentと同じ実反対称関係波生成子。"""

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

    def operator_matrix(self) -> np.ndarray:
        """JGを作り、警告の有無とは独立に全要素の有限性を検査する。"""

        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            operator = self.J @ self.G
        if not np.all(np.isfinite(operator)):
            raise FloatingPointError("nonfinite JG operator")
        return operator

    def sigma_spectrum(self) -> np.ndarray:
        operator = self.operator_matrix()
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            eigenvalues = np.linalg.eigvals(operator)
        if not np.all(np.isfinite(eigenvalues)):
            raise FloatingPointError("nonfinite JG eigenvalues")
        positive = eigenvalues.imag[eigenvalues.imag > 1e-12]
        return np.sort(positive)[::-1]


def eigenmode_residual(system: LowRankSystem, vector: np.ndarray) -> tuple[float, float]:
    kv = system.kmatvec(vector)
    mu = float(np.real(np.vdot(vector, 1j * kv)))
    residual = float(np.linalg.norm(1j * kv - mu * vector))
    return mu, residual


@dataclasses.dataclass
class SeedAttempt:
    attempt_number: int
    seed: int
    converged: bool
    iterations: int
    residual: float
    mu: float
    sigma_max: float
    parent_closure_abs: float
    failure_reason: str | None
    parent_vector: np.ndarray | None = dataclasses.field(repr=False)

    def audit_dict(self) -> dict[str, Any]:
        return {
            "attempt_number": self.attempt_number,
            "seed": str(self.seed),
            "converged": self.converged,
            "iterations": self.iterations,
            "residual": self.residual,
            "mu": self.mu,
            "sigma_max": self.sigma_max,
            "parent_closure_abs": self.parent_closure_abs,
            "failure_reason": self.failure_reason,
        }


@dataclasses.dataclass
class ParentResult:
    n: int
    m: int
    seed: int
    seed_was_explicit: bool
    raw_white_noise: np.ndarray
    parent_vector: np.ndarray
    relation_waves: np.ndarray
    edges: np.ndarray
    attempts: list[SeedAttempt]
    runtime_seconds: float


class ParentConstructionError(RuntimeError):
    def __init__(self, message: str, attempts: list[SeedAttempt]):
        super().__init__(message)
        self.attempts = attempts


def solve_parent(system: LowRankSystem, rng: np.random.Generator) -> SeedAttempt:
    """既存make_parentのargmin固有枝を一つのseed系列で一回解く。"""

    theta = rng.uniform(0.0, 2.0 * math.pi, system.m)
    vector: np.ndarray | None = None
    iteration = 0
    for iteration in range(1, SOLVER_ITERATIONS + 1):
        system.set_theta(theta)
        try:
            operator = system.operator_matrix()
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                eigenvalues, eigenvectors = np.linalg.eig(operator)
            if not (
                np.all(np.isfinite(eigenvalues))
                and np.all(np.isfinite(eigenvectors))
            ):
                raise FloatingPointError("nonfinite eigendecomposition output")
        except (np.linalg.LinAlgError, FloatingPointError):
            return SeedAttempt(
                0, 0, False, iteration, math.inf, math.nan, math.nan, math.inf,
                "eigendecomposition_failed", None,
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
            _, residual = eigenmode_residual(system, vector)
            if residual < SOLVER_TOLERANCE:
                break

    if vector is None:
        return SeedAttempt(
            0, 0, False, iteration, math.inf, math.nan, math.nan, math.inf,
            "nonfinite_or_zero_eigenvector", None,
        )
    system.set_theta(np.angle(vector))
    mu, residual = eigenmode_residual(system, vector)
    try:
        sigmas = system.sigma_spectrum()
    except (np.linalg.LinAlgError, FloatingPointError):
        return SeedAttempt(
            0, 0, False, iteration, residual, mu, math.nan, closure,
            "sigma_spectrum_failed", vector.copy(),
        )
    sigma_max = float(sigmas[0]) if len(sigmas) else 0.0
    closure = float(abs(complex(vector @ vector)))
    converged = bool(
        np.isfinite(residual)
        and residual < SOLVER_TOLERANCE
        and closure < 100.0 * np.finfo(float).eps * system.m
    )
    return SeedAttempt(
        0,
        0,
        converged,
        iteration,
        residual,
        mu,
        sigma_max,
        closure,
        None if converged else "self_consistency_tolerance_not_reached",
        vector.copy(),
    )


def white_closed_waveforms(
    rng: np.random.Generator, m: int, n: int
) -> tuple[np.ndarray, np.ndarray]:
    """M本のN点白色波形を、各行が零閉塞する形で構成する。"""

    n_half = (n - 1) // 2                             # 許容 |k| = 1..n_half（自己対は除外）
    signs = rng.integers(0, 2, size=(m, n_half))      # 0→+k, 1→−k（関係×|k| の硬貨）
    amps = rng.normal(size=(m, n_half, 2))            # 複素 Gaussian 振幅（Re, Im）
    raw = np.concatenate([signs[..., np.newaxis].astype(float), amps], axis=2)
    waves = np.zeros((m, n), dtype=np.complex128)
    grid = np.arange(n)
    for relation in range(m):
        for idx in range(n_half):
            k = idx + 1
            sk = k if signs[relation, idx] == 0 else (n - k)   # −k は bin n−k
            c = (amps[relation, idx, 0] + 1j * amps[relation, idx, 1])
            waves[relation] += c * np.exp(2j * np.pi * sk * grid / n)
        norm = np.linalg.norm(waves[relation])
        if norm > 0:
            waves[relation] /= norm
    return raw, waves


def _attempt_seed(n: int, seed: int, attempt_number: int) -> tuple[SeedAttempt, Any]:
    system = LowRankSystem(n)
    rng = np.random.default_rng(seed)
    attempt = solve_parent(system, rng)
    attempt.attempt_number = attempt_number
    attempt.seed = seed
    if not attempt.converged or attempt.parent_vector is None:
        return attempt, None

    raw, unit_waves = white_closed_waveforms(rng, system.m, n)
    relation_waves = attempt.parent_vector[:, np.newaxis] * unit_waves
    row_power = np.sum(np.abs(relation_waves) ** 2, axis=1)
    row_closure = np.sum(relation_waves * relation_waves, axis=1)
    if not np.all(np.isfinite(relation_waves)):
        attempt.converged = False
        attempt.failure_reason = "nonfinite_relation_wave"
        return attempt, None
    if float(np.max(np.abs(row_closure) / np.maximum(row_power, np.finfo(float).tiny))) >= 1e-12:
        attempt.converged = False
        attempt.failure_reason = "relation_wave_zero_closure_failed"
        return attempt, None
    edges = np.column_stack([system.edge_a, system.edge_b])
    return attempt, (raw, relation_waves, edges)


def make_parent(N: int, seed: int | None = None, max_retries: int = 3) -> ParentResult:
    """Nだけを理論入力としてM本のN点白色零閉塞関係波を作る。"""

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
    explicit = seed is not None
    limit = 1 if explicit else 1 + max_retries
    attempts: list[SeedAttempt] = []
    for attempt_number in range(1, limit + 1):
        current_seed = int(seed) if explicit else secrets.randbits(128)
        attempt, payload = _attempt_seed(N, current_seed, attempt_number)
        attempts.append(attempt)
        if attempt.converged:
            raw, waves, edges = payload
            return ParentResult(
                n=N,
                m=N * (N - 1) // 2,
                seed=current_seed,
                seed_was_explicit=explicit,
                raw_white_noise=raw,
                parent_vector=attempt.parent_vector.copy(),
                relation_waves=waves,
                edges=edges,
                attempts=attempts,
                runtime_seconds=time.time() - started,
            )
        if explicit:
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
        "function_contract": {
            "N": N,
            "seed_argument": None if seed is None else str(seed),
            "max_retries_argument": max_retries,
            "explicit_seed_means_single_attempt": True,
        },
        "derived_only_from_N": {
            "M": N * (N - 1) // 2,
            "lambda0": 2.0 * math.pi / N,
            "sample_count_per_relation_wave": N,
        },
        "fixed_solver_contract": {
            "iterations": SOLVER_ITERATIONS,
            "beta": SOLVER_BETA,
            "tolerance": SOLVER_TOLERANCE,
            "check_every": CHECK_EVERY,
            "branch": "argmin(Im eigenvalue of JG)",
        },
        "explicit_absences": {
            "harmonic_count_H": None,
            "preassigned_harmonic_orders": False,
            "preassigned_harmonic_amplitudes": False,
            "preassigned_harmonic_phases": False,
            "dc_component_k0": "projected out before QR (v3)",
            "nyquist_component_kN2": "projected out before QR for even N (v3)",
        },
    }


def write_success(
    output: Path,
    result: ParentResult,
    manifest: dict[str, Any],
    seed_search: list[dict[str, Any]] | None = None,
) -> None:
    output.mkdir(parents=True, exist_ok=False)
    np.save(output / "raw_white_noise.npy", result.raw_white_noise)
    np.save(output / "parent_vector.npy", result.parent_vector)
    np.save(output / "relation_waves.npy", result.relation_waves)
    np.save(output / "edges.npy", result.edges)
    row_power = np.sum(np.abs(result.relation_waves) ** 2, axis=1)
    row_closure = np.abs(np.sum(result.relation_waves**2, axis=1))
    manifest.update(
        {
            "status": "success",
            "accepted_seed": str(result.seed),
            "seed_was_explicit": result.seed_was_explicit,
            "attempts": [attempt.audit_dict() for attempt in result.attempts],
            "sequential_seed_search": seed_search,
            "arrays": {
                "raw_white_noise": {
                    "file": "raw_white_noise.npy",
                    "shape": list(result.raw_white_noise.shape),
                },
                "parent_vector": {
                    "file": "parent_vector.npy",
                    "shape": list(result.parent_vector.shape),
                },
                "relation_waves": {
                    "file": "relation_waves.npy",
                    "shape": list(result.relation_waves.shape),
                },
                "edges": {"file": "edges.npy", "shape": list(result.edges.shape)},
            },
            "invariant_audit": {
                "parent_self_consistency_residual": result.attempts[-1].residual,
                "parent_zero_closure_abs": result.attempts[-1].parent_closure_abs,
                "all_M_relation_waves_zero_closed": bool(
                    np.all(row_closure / np.maximum(row_power, np.finfo(float).tiny) < 1e-12)
                ),
                "max_relation_wave_closure_abs": float(np.max(row_closure)),
                "max_relation_wave_closure_relative": float(
                    np.max(row_closure / np.maximum(row_power, np.finfo(float).tiny))
                ),
                "nested_total_closure_abs": float(abs(complex(np.sum(result.relation_waves**2)))),
                "frobenius_norm": float(np.linalg.norm(result.relation_waves)),
            },
            "runtime_seconds": result.runtime_seconds,
        }
    )
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def write_failure(
    output: Path,
    manifest: dict[str, Any],
    attempts: list[SeedAttempt],
    message: str,
) -> None:
    """失敗したseed・残差を隠さず、配列を作らずに監査記録だけ保存する。"""

    output.mkdir(parents=True, exist_ok=False)
    manifest.update(
        {
            "status": "failure",
            "failure_message": message,
            "attempts": [attempt.audit_dict() for attempt in attempts],
            "arrays": None,
        }
    )
    (output / "failure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="NだけからM本の白色零閉塞関係波を生成")
    parser.add_argument("N", type=int)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--seed", type=int, default=None)
    group.add_argument("--search-seed-from", type=int, default=None)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"output already exists; refusing to overwrite: {output}")

    if args.search_seed_from is None:
        manifest = base_manifest(args.N, args.seed, args.max_retries)
        try:
            result = make_parent(args.N, seed=args.seed, max_retries=args.max_retries)
        except ParentConstructionError as error:
            write_failure(output, manifest, error.attempts, str(error))
            print(str(error), file=sys.stderr)
            raise SystemExit(2)
        write_success(output, result, manifest)
    else:
        if args.search_seed_from < 0:
            raise SystemExit("--search-seed-from must be >= 0")
        search_started = time.time()
        search: list[dict[str, Any]] = []
        seed = args.search_seed_from
        while True:
            attempt, payload = _attempt_seed(args.N, seed, len(search) + 1)
            search.append(attempt.audit_dict())
            print(
                f"seed={seed} converged={attempt.converged} residual={attempt.residual:.3e}",
                file=sys.stderr,
                flush=True,
            )
            if attempt.converged:
                raw, waves, edges = payload
                result = ParentResult(
                    n=args.N,
                    m=args.N * (args.N - 1) // 2,
                    seed=seed,
                    seed_was_explicit=True,
                    raw_white_noise=raw,
                    parent_vector=attempt.parent_vector.copy(),
                    relation_waves=waves,
                    edges=edges,
                    attempts=[attempt],
                    runtime_seconds=time.time() - search_started,
                )
                manifest = base_manifest(args.N, seed, args.max_retries)
                manifest["seed_selection_protocol"] = (
                    f"integer seeds tested sequentially from {args.search_seed_from}; "
                    "first converged seed accepted"
                )
                write_success(output, result, manifest, seed_search=search)
                break
            seed += 1

    print(
        json.dumps(
            {
                "status": "success",
                "N": result.n,
                "M": result.m,
                "accepted_seed": str(result.seed),
                "relation_waves_shape": list(result.relation_waves.shape),
                "output": str(output),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
