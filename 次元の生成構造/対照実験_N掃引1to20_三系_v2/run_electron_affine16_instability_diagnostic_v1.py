#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""electron affine16→64 instability diagnostic v1.

F v1 と既存 runner は変更せず、初期微小注入または full-step 後の射影を
外側から与える専用診断器。引数なしでは本走行しない。

  python3 run_electron_affine16_instability_diagnostic_v1.py --list
  python3 run_electron_affine16_instability_diagnostic_v1.py --smoke
  python3 run_electron_affine16_instability_diagnostic_v1.py --all
  python3 run_electron_affine16_instability_diagnostic_v1.py --arm A_standard
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_nsweep_three_series_v2.py"
BASE = HERE / "run_tb_nsweep_1to20_v1.py"
PREREG = HERE / "事前登録_electron_affine16_16to64診断_v1.md"
REFERENCE = HERE / "nsweep_electron_T42000_d0.1_N12_v2.npz"
MANIFEST = HERE / "manifest_electron_affine16_instability_v1.json"
UF = HERE.parent / "統一万能関数_v1"

N = 12
NN = 16
NETA = 8
T_FULL = 42000
DELTA = 0.1
EPSILONS = (1e-15, 1e-10, 1e-5)
ALLOWED_CELL = (1, 1)
FORBIDDEN_CELL = (1, 0)
THRESHOLDS = (1e-30, 1e-24, 1e-20, 1e-16, 1e-12,
              1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
RAW_EVERY = 500
RAW_FIXED = (0, 1, 2, 10, 100, 1000)

EXPECTED_HASHES = {
    RUNNER: "c4b3b33a82657232734223f089cff79e103ce7a18fb5a064fc71a058d1d006fe",
    BASE: "1615b8769cf765a444a8e057718ff4a760f040376335f434f9903b8cf4f9dba4",
    UF / "unified_interaction_v1.py":
        "b4db659b9bf958246d192c852a48fd97f5793f9bf2f891d5c75ce22ae5a5ed63",
    UF / "unified_dimension_v1.py":
        "cfc9cf5d835e39ec6959eeeadff102e5dd23b93deed9dfeff52778e3ee48830a",
    UF / "unified_readout_v3.py":
        "4dec07c9811c3c06b3c34941ade85593ac3feecbeb82e6902cecad1263958953",
    UF / "selection_v1.py":
        "90d546e61e20ce58b4e9a962bdcd761ef96d46a9061a60decca2586ad324eed2",
}
EXPECTED_REFERENCE_SHA256 = (
    "a42375dd17c6ba698de537fc70fd307f5b6380fdaa20bfd5a5c4187edecefde4"
)


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_frozen_sources(expected: dict[Path, str] = EXPECTED_HASHES) -> dict[str, str]:
    actual: dict[str, str] = {}
    bad: list[str] = []
    for path, wanted in expected.items():
        got = sha256(path) if path.is_file() else None
        actual[path.name] = got or "MISSING"
        if got != wanted:
            bad.append(f"{path.name}: expected {wanted}, got {got}")
    if bad:
        raise RuntimeError("frozen source SHA-256 mismatch: " + " | ".join(bad))
    got_ref = sha256(REFERENCE) if REFERENCE.is_file() else None
    if got_ref != EXPECTED_REFERENCE_SHA256:
        raise RuntimeError(
            f"reference SHA-256 mismatch: expected {EXPECTED_REFERENCE_SHA256}, got {got_ref}"
        )
    return actual


def load_existing_runner():
    """electron, N=12, T=42000, delta=.1 として既存runnerをread-only import。"""
    cache = Path(tempfile.gettempdir()) / "electron_affine16_diag_mpl_v1"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    old_argv = sys.argv[:]
    try:
        sys.argv = [str(RUNNER), "electron", "12", "12", "42000", "0.1"]
        spec = importlib.util.spec_from_file_location("electron_affine16_base_v1", RUNNER)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot load {RUNNER}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        sys.argv = old_argv


BASEMOD = load_existing_runner()


def ideal_eta(k: int) -> int:
    return int((-3 * (k - 2)) % NETA)


IDEAL_MASK = np.zeros((NN, NETA), dtype=bool)
for _k in range(NN):
    IDEAL_MASK[_k, ideal_eta(_k)] = True
PARITY_MASK = np.fromfunction(
    lambda k, e: (k.astype(int) % 2) == (e.astype(int) % 2),
    (NN, NETA), dtype=int,
).astype(bool)
ALLOWED_OFF_MASK = PARITY_MASK & ~IDEAL_MASK
FORBIDDEN_MASK = ~PARITY_MASK
assert int(IDEAL_MASK.sum()) == 16
assert int(ALLOWED_OFF_MASK.sum()) == 48
assert int(FORBIDDEN_MASK.sum()) == 64
assert np.all((IDEAL_MASK.astype(int) + ALLOWED_OFF_MASK.astype(int)
               + FORBIDDEN_MASK.astype(int)) == 1)
assert ideal_eta(1) == 3 and ideal_eta(2) == 0 and ideal_eta(3) == 5
assert ALLOWED_OFF_MASK[ALLOWED_CELL] and not IDEAL_MASK[ALLOWED_CELL]
assert FORBIDDEN_MASK[FORBIDDEN_CELL] and not PARITY_MASK[FORBIDDEN_CELL]


@dataclass(frozen=True)
class Arm:
    arm_id: str
    kind: str
    epsilon: float = 0.0
    inject_cell: tuple[int, int] | None = None

    @property
    def tag(self) -> str:
        if self.epsilon:
            return f"{self.arm_id}_e{self.epsilon:g}"
        return self.arm_id

    @property
    def npz_path(self) -> Path:
        return HERE / f"electron_affine16_{self.tag}_T{T_FULL}_v1.npz"

    @property
    def json_path(self) -> Path:
        return HERE / f"electron_affine16_{self.tag}_T{T_FULL}_v1.json"


ARMS = (
    Arm("A_standard", "standard"),
    Arm("B_project16", "project"),
    *(Arm("C_allowed", "inject_allowed", e, ALLOWED_CELL) for e in EPSILONS),
    *(Arm("D_forbidden", "inject_forbidden", e, FORBIDDEN_CELL) for e in EPSILONS),
)
ARM_BY_TAG = {a.tag: a for a in ARMS}


def total_power(C2: np.ndarray) -> float:
    return float(np.sum(np.abs(C2) ** 2, dtype=np.float64))


def state_norm(C2: np.ndarray) -> float:
    return float(np.linalg.norm(C2.reshape(-1)))


def normalize_to_norm(C2: np.ndarray, target_norm: float) -> tuple[float, float]:
    """全状態を一様再正規化し、保存する norm 読出しをbit exactに合わせる。"""
    before = state_norm(C2)
    if before == 0.0:
        raise FloatingPointError("cannot normalize zero state")
    original = C2.copy()
    center = target_norm / before
    candidates = [center]
    lo = hi = center
    for _ in range(32):
        lo = float(np.nextafter(lo, -np.inf))
        hi = float(np.nextafter(hi, np.inf))
        candidates.extend((lo, hi))
    best: tuple[float, float, np.ndarray] | None = None
    for scale in candidates:
        trial = original * scale
        error = state_norm(trial) - target_norm
        if best is None or abs(error) < abs(best[1]):
            best = (scale, error, trial)
        if error == 0.0:
            C2[...] = trial
            return float(scale), 0.0
    assert best is not None
    C2[...] = best[2]
    raise AssertionError(
        f"bit-exact norm matching unavailable: target={target_norm:.17g}, "
        f"best_error={best[1]:.3e}, scale={best[0]:.17g}"
    )


def sector_vector(C2: np.ndarray) -> np.ndarray:
    P = np.sum(np.abs(C2) ** 2, axis=0)
    ideal = float(P[IDEAL_MASK].sum())
    allowed = float(P[ALLOWED_OFF_MASK].sum())
    forbidden = float(P[FORBIDDEN_MASK].sum())
    total = float(P.sum())
    odd = float(P[np.arange(NN) % 2 == 1].sum())
    even_nonpump = float(P[(np.arange(NN) % 2 == 0) & (np.arange(NN) != 2)].sum())
    return np.array((ideal, allowed, forbidden, total, odd, even_nonpump), dtype=float)


SECTOR_NAMES = (
    "ideal16_power", "allowed_off_power", "forbidden_power", "total_power",
    "odd_power", "even_nonpump_power",
)


def initial_power_ledger(C2: np.ndarray, injected: tuple[int, int] | None) -> dict[str, float]:
    P = np.sum(np.abs(C2) ** 2, axis=0)
    odd = float(P[np.arange(NN) % 2 == 1].sum())
    even_nonpump = float(P[(np.arange(NN) % 2 == 0) & (np.arange(NN) != 2)].sum())
    return {
        "PF": odd,
        "PB": even_nonpump,
        "Ptotal": float(P.sum()),
        "norm": state_norm(C2),
        "pump_power": float(P[2, 0]),
        "primary_seed_power": float(P[1, 3]),
        "injected_cell_power": float(P[injected]) if injected is not None else 0.0,
    }


CURRENT_ARM: Arm = ARMS[0]
CURRENT_T = T_FULL
DIAG_ENGINES: list["DiagnosticEngine"] = []


def raw_steps(T: int) -> set[int]:
    return ({x for x in RAW_FIXED if x <= T}
            | set(range(RAW_EVERY, T + 1, RAW_EVERY)) | {T})


class DiagnosticEngine(BASEMOD.RecordingEngine):
    """F v1 full-stepをそのまま実行し、その外側で診断・任意射影だけを行う。"""

    def __init__(self, *args, arm: Arm, is_vacuum: bool, initial_ledger: dict, **kwargs):
        super().__init__(*args, **kwargs)
        self.diag_arm = arm
        self.is_vacuum = is_vacuum
        self.initial_ledger = initial_ledger
        self.sector_pre: list[np.ndarray] = []
        self.sector_post: list[np.ndarray] = []
        self.projection_removed: list[float] = []
        self.projection_scale: list[float] = []
        self.norm_pre_step: list[float] = []
        self.norm_post: list[float] = []
        self.norm_error: list[float] = []
        self.raw_post_steps: list[int] = [0] if not is_vacuum else []
        self.raw_post: list[np.ndarray] = [self.C2().copy()] if not is_vacuum else []
        self.raw_pre_steps: list[int] = []
        self.raw_pre: list[np.ndarray] = []
        self._diag_step = 0
        self._raw_wanted = raw_steps(CURRENT_T)

    def step(self, *args, **kwargs):
        norm_before = state_norm(self.C2())
        out = super().step(*args, **kwargs)  # F v1 full-step + existing recorder
        self._diag_step += 1
        C2 = self.C2()
        pre = sector_vector(C2)
        removed = 0.0
        scale = 1.0
        norm_error = state_norm(C2) - norm_before
        if self.diag_arm.kind == "project" and not self.is_vacuum:
            if self._diag_step in self._raw_wanted:
                self.raw_pre_steps.append(self._diag_step)
                self.raw_pre.append(C2.copy())
            removed = float(np.sum(np.abs(C2[:, ~IDEAL_MASK]) ** 2))
            C2[:, ~IDEAL_MASK] = 0.0
            scale, norm_error = normalize_to_norm(C2, norm_before)
        post = sector_vector(C2)
        if self.diag_arm.kind == "project" and not self.is_vacuum:
            if post[1] != 0.0 or post[2] != 0.0:
                raise AssertionError("projection arm retains off-orbit power")
        self.sector_pre.append(pre)
        self.sector_post.append(post)
        self.projection_removed.append(removed)
        self.projection_scale.append(scale)
        self.norm_pre_step.append(norm_before)
        self.norm_post.append(state_norm(C2))
        self.norm_error.append(norm_error)
        if not self.is_vacuum and self._diag_step in self._raw_wanted:
            self.raw_post_steps.append(self._diag_step)
            self.raw_post.append(C2.copy())
        return out


def build_diag_universe(n: int, delta: float, Nn: int = 5, Neta: int = 8, seed: int = 2):
    """既存electron buildと同じ。delta>0のときだけ宣言した診断腕を加える。"""
    m = n * (n - 1) // 2
    _, _v, _, _, _, _, _, Z0c, wp0 = BASEMOD.F1.abl.build_init(n, False)
    r2 = BASEMOD.F1.gen3.make_parent(n, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    C2_0[:, 2, 0] = Z0c
    is_vacuum = not (delta > 0)
    if not is_vacuum:
        C2_0[:, 1, 3] += DELTA * seed_state
        standard_norm = state_norm(C2_0)
        if CURRENT_ARM.inject_cell is not None:
            kk, ee = CURRENT_ARM.inject_cell
            C2_0[:, kk, ee] += CURRENT_ARM.epsilon * seed_state
            normalize_to_norm(C2_0, standard_norm)
    initial = initial_power_ledger(C2_0, CURRENT_ARM.inject_cell if not is_vacuum else None)
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    with np.errstate(divide="ignore", invalid="ignore"):
        q2 = q2 / np.linalg.norm(q2)
    engine = DiagnosticEngine(
        n, C2_0, wp0, arm=CURRENT_ARM, is_vacuum=is_vacuum,
        initial_ledger=initial,
    )
    DIAG_ENGINES.append(engine)
    return engine, p2, q2


def array_bytes(a: np.ndarray) -> bytes:
    return np.ascontiguousarray(a).tobytes(order="C")


def base_payload(Hm, Rm, Am, Ccm, Csm, Hv, Av,
                 hm: dict[str, np.ndarray], hv: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    payload: dict[str, np.ndarray] = {
        **{f"m_{k}": Hm[k] for k in BASEMOD.ns.KEYS},
        **{f"v_{k}": Hv[k] for k in BASEMOD.ns.KEYS},
        "m_resid": Rm, "m_acq": Am, "v_acq": Av,
        "m_cond_closure": Ccm, "m_seed_closure": Csm,
        **{f"rec_m_{k}": v for k, v in hm.items()},
        **{f"rec_v_{k}": v for k, v in hv.items()},
        "rec_m_conc_partner": hm["conc_k3"],
        "rec_m_partner_dom_m": hm["dom_m"],
        "rec_m_partner_q_hat": hm["q_hat"],
        "rec_m_partner_readable": hm["readable"],
        "rec_m_partner_power": hm["k3_power"],
        "rec_v_conc_partner": hv["conc_k3"],
        "rec_v_partner_dom_m": hv["dom_m"],
        "rec_v_partner_q_hat": hv["q_hat"],
        "rec_v_partner_readable": hv["readable"],
        "rec_v_partner_power": hv["k3_power"],
        "seed_cells_index": np.array([(1, 3)], dtype=int),
        "seed_cells_delta": np.array([DELTA], dtype=float),
        "targets_index_2d": np.array([(3, 5)], dtype=int),
        "targets_index": np.array([(3, 5)], dtype=float),
    }
    return payload


def first_ge(series: np.ndarray, threshold: float) -> int | None:
    idx = np.flatnonzero(np.asarray(series) >= threshold)
    return int(idx[0]) + 1 if idx.size else None


def threshold_summary(sectors: np.ndarray) -> dict[str, dict[str, int | None]]:
    total = sectors[:, 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        allowed = sectors[:, 1] / total
        forbidden = sectors[:, 2] / total
        off = (sectors[:, 1] + sectors[:, 2]) / total
    return {
        name: {f"{q:.0e}": first_ge(values, q) for q in THRESHOLDS}
        for name, values in (
            ("allowed_off_fraction", allowed),
            ("forbidden_fraction", forbidden),
            ("off_total_fraction", off),
        )
    }


def diagnostic_payload(engine: DiagnosticEngine) -> dict[str, np.ndarray]:
    pre = np.asarray(engine.sector_pre, dtype=float)
    post = np.asarray(engine.sector_post, dtype=float)
    raw_post = np.asarray(engine.raw_post, dtype=complex)
    raw_pre = np.asarray(engine.raw_pre, dtype=complex)
    if raw_pre.size == 0:
        raw_pre = np.empty((0, engine.m, NN, NETA), dtype=complex)
    return {
        "diag_sector_names": np.asarray(SECTOR_NAMES, dtype="U32"),
        "diag_sector_pre": pre,
        "diag_sector_post": post,
        "diag_projection_removed_power": np.asarray(engine.projection_removed, dtype=float),
        "diag_projection_scale": np.asarray(engine.projection_scale, dtype=float),
        "diag_norm_pre_step": np.asarray(engine.norm_pre_step, dtype=float),
        "diag_norm_post": np.asarray(engine.norm_post, dtype=float),
        "diag_norm_error": np.asarray(engine.norm_error, dtype=float),
        "diag_raw_post_steps": np.asarray(engine.raw_post_steps, dtype=int),
        "diag_raw_post_C": raw_post,
        "diag_raw_pre_steps": np.asarray(engine.raw_pre_steps, dtype=int),
        "diag_raw_pre_C": raw_pre,
        "diag_ideal_mask": IDEAL_MASK,
        "diag_allowed_off_mask": ALLOWED_OFF_MASK,
        "diag_forbidden_mask": FORBIDDEN_MASK,
    }


def compare_standard(candidate: dict[str, np.ndarray]) -> dict[str, Any]:
    checked = 0
    with np.load(REFERENCE, allow_pickle=False) as ref:
        missing = [k for k in ref.files if k not in candidate]
        if missing:
            raise AssertionError(f"standard candidate missing reference keys: {missing}")
        for key in ref.files:
            x = np.asarray(ref[key])
            y = np.asarray(candidate[key])
            if x.dtype != y.dtype or x.shape != y.shape:
                raise AssertionError(
                    f"standard {key}: {x.dtype}/{x.shape} != {y.dtype}/{y.shape}"
                )
            if array_bytes(x) != array_bytes(y):
                raise AssertionError(f"standard {key}: bytes differ")
            checked += 1
    return {
        "reference": REFERENCE.name,
        "arrays_checked": checked,
        "byte_exact": True,
        "max_abs_difference": 0.0,
    }


def compare_standard_prefix(candidate: dict[str, np.ndarray], steps: int) -> dict[str, Any]:
    """smoke軌道を既存T42000正本の同じ先頭へ、全99キーで照合する。"""
    checked = 0
    with np.load(REFERENCE, allow_pickle=False) as ref:
        missing = [k for k in ref.files if k not in candidate]
        if missing:
            raise AssertionError(f"standard prefix missing reference keys: {missing}")
        for key in ref.files:
            x = np.asarray(ref[key])
            y = np.asarray(candidate[key])
            if x.shape != y.shape:
                if (x.ndim < 1 or y.ndim != x.ndim or
                        x.shape[1:] != y.shape[1:] or x.shape[0] < y.shape[0]):
                    raise AssertionError(
                        f"standard prefix {key}: incompatible {x.shape} vs {y.shape}"
                    )
                x = x[:y.shape[0]]
            if x.dtype != y.dtype or x.shape != y.shape or array_bytes(x) != array_bytes(y):
                raise AssertionError(f"standard prefix {key}: bytes differ")
            checked += 1
    return {
        "reference": REFERENCE.name,
        "prefix_steps": steps,
        "arrays_checked": checked,
        "byte_exact": True,
        "max_abs_difference": 0.0,
    }


def compare_vacuum(candidate: dict[str, np.ndarray], steps: int) -> dict[str, Any]:
    """各腕の内部真空46配列を標準electron正本とbyte照合する。"""
    checked = 0
    prefix = steps != T_FULL
    with np.load(REFERENCE, allow_pickle=False) as ref:
        keys = [k for k in ref.files if k.startswith("v_") or k.startswith("rec_v_")]
        missing = [k for k in keys if k not in candidate]
        if missing:
            raise AssertionError(f"vacuum candidate missing reference keys: {missing}")
        for key in keys:
            x = np.asarray(ref[key])
            y = np.asarray(candidate[key])
            if prefix and x.shape != y.shape:
                if (x.ndim < 1 or y.ndim != x.ndim or
                        x.shape[1:] != y.shape[1:] or x.shape[0] < y.shape[0]):
                    raise AssertionError(
                        f"vacuum prefix {key}: incompatible {x.shape} vs {y.shape}"
                    )
                x = x[:y.shape[0]]
            if x.dtype != y.dtype or x.shape != y.shape or array_bytes(x) != array_bytes(y):
                raise AssertionError(f"vacuum {key}: bytes differ")
            checked += 1
    return {
        "reference": REFERENCE.name,
        "comparison": "prefix" if prefix else "full",
        "prefix_steps": steps if prefix else None,
        "arrays_checked": checked,
        "byte_exact": True,
        "max_abs_difference": 0.0,
    }


def run_arm(arm: Arm, T: int, save: bool) -> dict[str, Any]:
    global CURRENT_ARM, CURRENT_T
    CURRENT_ARM, CURRENT_T = arm, T
    BASEMOD.ns.T = T
    BASEMOD.ns.WIN = (T // 2, T)
    BASEMOD.ns.DELTA = DELTA
    BASEMOD.ns.F.build_standard_universe = build_diag_universe
    DIAG_ENGINES.clear()
    t0 = time.time()
    Hm, Rm, Am, Ccm, Csm = BASEMOD.ns.run_one(N, DELTA)
    matter = DIAG_ENGINES[-1]
    Hv, _Rv, Av, _Ccv, _Csv = BASEMOD.ns.run_one(N, 0.0)
    vacuum = DIAG_ENGINES[-1]
    hm, hv = BASEMOD.arrays(matter), BASEMOD.arrays(vacuum)
    payload = base_payload(Hm, Rm, Am, Ccm, Csm, Hv, Av, hm, hv)
    payload.update(diagnostic_payload(matter))
    sectors_pre = np.asarray(matter.sector_pre, dtype=float)
    sectors_post = np.asarray(matter.sector_post, dtype=float)
    partition_error = np.max(np.abs(
        sectors_post[:, 0] + sectors_post[:, 1] + sectors_post[:, 2]
        - sectors_post[:, 3]
    )) if T else 0.0
    summary: dict[str, Any] = {
        "arm": arm.tag,
        "kind": arm.kind,
        "epsilon": arm.epsilon,
        "inject_cell": list(arm.inject_cell) if arm.inject_cell else None,
        "conditions": {
            "N": N, "M": N * (N - 1) // 2, "Nn": NN, "Neta": NETA,
            "T": T, "delta": DELTA, "seed": 2, "cell": [2, 0],
            "order": 6, "window": [T // 2, T],
            "injection_vector": "Csec[:,1] / ||Csec[:,1]|| (same complex phase as primary seed)",
        },
        "ideal16": [[k, ideal_eta(k)] for k in range(NN)],
        "sector_counts": {"ideal16": 16, "allowed_off": 48, "forbidden": 64},
        "initial_power": matter.initial_ledger,
        "partition_max_abs_error": float(partition_error),
        "threshold_times_pre": threshold_summary(sectors_pre),
        "threshold_times_post": threshold_summary(sectors_post),
        "projection": {
            "removed_power_max": float(np.max(matter.projection_removed)),
            "post_allowed_off_max": float(np.max(sectors_post[:, 1])),
            "post_forbidden_max": float(np.max(sectors_post[:, 2])),
            "norm_error_max_abs": float(np.max(np.abs(matter.norm_error))),
        },
        "raw_state_steps_post": matter.raw_post_steps,
        "raw_state_steps_pre": matter.raw_pre_steps,
        "recording_semantics": {
            "legacy_rec_m": (
                "F-v1 full-step state before the external B projection; unchanged RecordingEngine"
            ),
            "m_H_D_G_S": "state after the external B projection and exact norm restoration",
            "diag_sector_pre": "F-v1 full-step state immediately before B projection",
            "diag_sector_post": "state after B projection and exact norm restoration",
        },
        "runtime_sec": time.time() - t0,
    }
    summary["vacuum_exact_control"] = compare_vacuum(payload, T)
    if arm.kind == "standard" and T == T_FULL:
        summary["standard_exact_control"] = compare_standard(payload)
    elif arm.kind == "standard":
        summary["standard_prefix_control"] = compare_standard_prefix(payload, T)
    if save:
        for path in (arm.npz_path, arm.json_path):
            if path.exists():
                raise FileExistsError(f"refusing overwrite: {path.name}")
        np.savez_compressed(arm.npz_path, **payload)
        summary["npz_sha256"] = sha256(arm.npz_path)
        arm.json_path.write_text(
            json.dumps(summary, ensure_ascii=False, indent=2, default=float) + "\n",
            encoding="utf-8",
        )
        summary["json_sha256"] = sha256(arm.json_path)
    return summary


def load_manifest() -> dict[str, Any]:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text(encoding="utf-8"))
    return {
        "schema": "electron-affine16-instability-v1",
        "created_at": now_iso(),
        "conditions": {
            "N": N, "T": T_FULL, "delta": DELTA,
            "epsilons": list(EPSILONS),
            "allowed_cell": list(ALLOWED_CELL),
            "forbidden_cell": list(FORBIDDEN_CELL),
        },
        "runs": {},
    }


def save_manifest(data: dict[str, Any]) -> None:
    data["updated_at"] = now_iso()
    MANIFEST.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, default=float) + "\n",
        encoding="utf-8",
    )


def control_passed(manifest: dict[str, Any]) -> bool:
    entry = manifest.get("runs", {}).get("A_standard", {})
    return bool(entry.get("summary", {}).get("standard_exact_control", {}).get("byte_exact"))


def execute_full_arm(arm: Arm, manifest: dict[str, Any]) -> None:
    if arm.kind != "standard" and not control_passed(manifest):
        raise RuntimeError("A_standard exact control must pass before B/C/D")
    verify_frozen_sources()
    if arm.tag in manifest.get("runs", {}):
        raise FileExistsError(f"manifest already contains arm {arm.tag}")
    manifest["runs"][arm.tag] = {"status": "running", "started_at": now_iso()}
    save_manifest(manifest)
    try:
        summary = run_arm(arm, T_FULL, save=True)
    except Exception as exc:
        manifest["runs"][arm.tag].update({
            "status": "failed", "finished_at": now_iso(),
            "error": f"{type(exc).__name__}: {exc}",
        })
        save_manifest(manifest)
        raise
    manifest["runs"][arm.tag].update({
        "status": "completed", "finished_at": now_iso(), "summary": summary,
        "artifacts": [arm.npz_path.name, arm.json_path.name],
    })
    save_manifest(manifest)


def smoke() -> dict[str, Any]:
    """短いin-memory構造検査。成果物・manifestは作らない。"""
    reports = []
    standard_norm: float | None = None
    for arm in ARMS:
        result = run_arm(arm, T=4, save=False)
        initial = result["initial_power"]
        if standard_norm is None:
            standard_norm = initial["norm"]
        ok = (result["partition_max_abs_error"] <= 1e-14
              and initial["norm"] == standard_norm
              and result["vacuum_exact_control"]["byte_exact"])
        if arm.kind == "project":
            p = result["projection"]
            ok = ok and p["post_allowed_off_max"] == 0.0 and p["post_forbidden_max"] == 0.0
        if arm.inject_cell is not None:
            expected_power = (arm.epsilon ** 2 * standard_norm ** 2
                              / (standard_norm ** 2 + arm.epsilon ** 2))
            ok = ok and math.isclose(
                initial["injected_cell_power"], expected_power,
                rel_tol=5e-14, abs_tol=0.0,
            )
        else:
            ok = ok and initial["injected_cell_power"] == 0.0
        if arm.kind == "standard":
            ok = ok and result["standard_prefix_control"]["byte_exact"]
        reports.append({
            "arm": arm.tag,
            "ok": bool(ok),
            "initial_power": result["initial_power"],
            "projection": result["projection"],
            "vacuum_exact_control": result["vacuum_exact_control"],
            "standard_prefix_control": result.get("standard_prefix_control"),
        })
    if not all(r["ok"] for r in reports):
        raise AssertionError(f"smoke failed: {reports}")
    return {
        "smoke": "PASS", "T": 4,
        "mask_counts": [int(IDEAL_MASK.sum()), int(ALLOWED_OFF_MASK.sum()),
                        int(FORBIDDEN_MASK.sum())],
        "arms": reports,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--list", action="store_true")
    group.add_argument("--smoke", action="store_true")
    group.add_argument("--all", action="store_true")
    group.add_argument("--arm", choices=tuple(ARM_BY_TAG))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.list:
        print(f"full_runs={len(ARMS)} N={N} T={T_FULL} delta={DELTA}")
        for i, arm in enumerate(ARMS, 1):
            print(i, arm.tag, arm.kind, arm.inject_cell, arm.epsilon)
        return 0
    if args.smoke:
        print(json.dumps(smoke(), ensure_ascii=False, indent=2, default=float))
        return 0

    sources_before = verify_frozen_sources()
    manifest = load_manifest()
    manifest["runner_sha256"] = sha256(Path(__file__).resolve())
    manifest["preregistration_sha256"] = sha256(PREREG)
    manifest["source_sha256_before"] = sources_before
    manifest["reference_sha256_before"] = sha256(REFERENCE)
    save_manifest(manifest)
    selected = ARMS if args.all else (ARM_BY_TAG[args.arm],)
    for arm in selected:
        execute_full_arm(arm, manifest)
    manifest["source_sha256_after"] = verify_frozen_sources()
    manifest["reference_sha256_after"] = sha256(REFERENCE)
    manifest["status"] = "complete" if args.all else "partial_complete"
    manifest["finished_at"] = now_iso()
    save_manifest(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
