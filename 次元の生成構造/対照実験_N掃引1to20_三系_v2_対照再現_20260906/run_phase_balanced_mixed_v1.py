#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-balanced mixed comparator for the fixed N=12, T=42000 experiment.

This is an upper-layer initial-condition wrapper.  It does not edit or fork the
F v1 dynamics.  The five k=1 fermion-family cells remain coherent.  In the
``phase-balanced`` arm only, the three k=6 B3 cells receive phases

    (k, eta) = (6, 0), (6, 3), (6, 5)
    phi      = 0,      2*pi/3, 4*pi/3,

so their complex vector sum is zero to floating-point tolerance while
``PF=5*delta**2``, ``PB=3*delta**2``, and the occupied addresses are unchanged.

The scientific run is deliberately locked to the existing comparator:
N=12, T=42000, Nn=16, Neta=8, delta=0.04357, seed=2, cell=(2,0), order=6,
and the latter half-window.  A separate explicit flag is required before a
T=42000 run can start.  Use ``smoke`` first; it writes no experiment outputs.

Examples (do not start the science runs before implementation review)::

    python3 run_phase_balanced_mixed_v1.py smoke
    python3 run_phase_balanced_mixed_v1.py run --arm control \
        --replicate ctl1 --allow-science-run
    python3 run_phase_balanced_mixed_v1.py run --arm phase-balanced \
        --replicate r1 --allow-science-run
    python3 run_phase_balanced_mixed_v1.py run --arm phase-balanced \
        --replicate r2 --compare-replicate r1 --allow-science-run

The control is checked array-by-array against the existing untagged mixed
delta=0.04357 NPZ.  A second phase-balanced replicate can likewise be required
to be byte-identical array-by-array to the first.  Every run has a distinct
``_rep-...`` tag and a SHA-256 manifest; existing files are never overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BASE_RUNNER = HERE / "run_nsweep_three_series_v2.py"
UF = HERE.parent / "統一万能関数_v1"

FIXED_N = 12
FIXED_T = 42000
FIXED_DELTA = 0.04357
FIXED_NN = 16
FIXED_NETA = 8
SMOKE_STEPS = 4

ARM_CONTROL = "control"
ARM_BALANCED = "phase-balanced"
ARMS = (ARM_CONTROL, ARM_BALANCED)

B3_PHASES = {
    (6, 0): 0.0,
    (6, 3): 2.0 * math.pi / 3.0,
    (6, 5): 4.0 * math.pi / 3.0,
}

SOURCE_PATHS = {
    "phase_runner": Path(__file__).resolve(),
    "base_runner": BASE_RUNNER,
    "mother_nsweep": HERE / "run_tb_nsweep_1to20_v1.py",
    "F_v1": UF / "unified_interaction_v1.py",
    "D_v1": UF / "unified_dimension_v1.py",
    "G_v3": UF / "unified_readout_v3.py",
    "S_v1": UF / "selection_v1.py",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(p) for p in SOURCE_PATHS.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("required source missing: " + ", ".join(missing))
    return {name: sha256(path) for name, path in SOURCE_PATHS.items()}


def safe_replicate(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,39}", value):
        raise argparse.ArgumentTypeError(
            "replicate must be 1-40 ASCII letters, digits, '_' or '-', "
            "starting with a letter or digit"
        )
    return value


def output_suffix(arm: str, replicate: str) -> str:
    prefix = "cohmixv1" if arm == ARM_CONTROL else "pbmixv1"
    return f"{prefix}-{replicate}"


def base_argv(t_steps: int, delta: float, suffix: str) -> list[str]:
    return [
        str(BASE_RUNNER), "mixed", str(FIXED_N), str(FIXED_N),
        str(t_steps), f"{delta:.17g}", suffix,
    ]


def load_base(t_steps: int, delta: float, suffix: str):
    """Load the existing runner without calling its main()."""
    sys.argv = base_argv(t_steps, delta, suffix)
    name = f"phase_balanced_base_{os.getpid()}_{suffix}"
    spec = importlib.util.spec_from_file_location(name, BASE_RUNNER)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {BASE_RUNNER}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    if mod.MODE != "mixed":
        raise AssertionError(f"base mode changed unexpectedly: {mod.MODE}")
    if (mod.NN, mod.NETA) != (FIXED_NN, FIXED_NETA):
        raise AssertionError(
            f"register mismatch: {(mod.NN, mod.NETA)} != "
            f"{(FIXED_NN, FIXED_NETA)}"
        )
    return mod


def cell_phases(base, arm: str) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for k, eta, amplitude in base.CELLS:
        phi = B3_PHASES[(k, eta % base.NETA)] if (
            arm == ARM_BALANCED and k == 6
        ) else 0.0
        phasor = float(amplitude) * complex(math.cos(phi), math.sin(phi))
        rows.append({
            "k": int(k),
            "eta_index": int(eta % base.NETA),
            "amplitude": float(amplitude),
            "phase_rad": float(phi),
            "phasor_real": float(phasor.real),
            "phasor_imag": float(phasor.imag),
        })
    return rows


def phase_summary(base, arm: str) -> dict[str, object]:
    rows = cell_phases(base, arm)
    total = sum(complex(r["phasor_real"], r["phasor_imag"]) for r in rows)
    fsum = sum(
        complex(r["phasor_real"], r["phasor_imag"])
        for r in rows if int(r["k"]) % 2 == 1
    )
    bsum = total - fsum
    amplitude_sum = sum(float(r["amplitude"]) for r in rows)
    tol = 64.0 * np.finfo(float).eps * max(
        amplitude_sum, np.finfo(float).tiny
    )
    if arm == ARM_BALANCED:
        if abs(bsum) > tol:
            raise AssertionError(
                f"B3 coherent-sum residual {abs(bsum):.17g} > {tol:.17g}"
            )
        expected = 5.0 * FIXED_DELTA
        if abs(abs(total) - expected) > tol:
            raise AssertionError(
                f"total coherent amplitude {abs(total):.17g} != {expected:.17g}"
            )
    return {
        "arm": arm,
        "cells": rows,
        "B3_phase_order": [[6, 0], [6, 3], [6, 5]],
        "B3_phases_rad": [B3_PHASES[(6, e)] for e in (0, 3, 5)],
        "coherent_sum": {"real": total.real, "imag": total.imag, "abs": abs(total)},
        "fermion_coherent_sum": {
            "real": fsum.real, "imag": fsum.imag, "abs": abs(fsum)
        },
        "B3_coherent_sum": {
            "real": bsum.real, "imag": bsum.imag, "abs": abs(bsum)
        },
        "B3_residual_tolerance": tol,
        "PF": 5.0 * FIXED_DELTA ** 2,
        "PB": 3.0 * FIXED_DELTA ** 2,
        "Pseed": 8.0 * FIXED_DELTA ** 2,
        "nominal_amplitude_sum": amplitude_sum,
        "addresses_unchanged_from_mixed": True,
    }


def phase_factor(phi: float, amplitude: float) -> complex | float:
    # Keeping phi=0 as a real multiply is necessary for an exact control.
    if phi == 0.0:
        return amplitude
    return amplitude * complex(math.cos(phi), math.sin(phi))


def build_phase_universe(base, arm: str, n: int, delta: float,
                         Nn: int = 5, Neta: int = 8, seed: int = 2):
    """Copy only the base runner's initial-state recipe; use F v1 unchanged."""
    if (Nn, Neta, seed) != (FIXED_NN, FIXED_NETA, 2):
        raise AssertionError(
            f"fixed recipe violated: Nn={Nn}, Neta={Neta}, seed={seed}"
        )
    m = n * (n - 1) // 2
    _, _v, _, _, _, _, _, z0c, wp0 = base.F1.abl.build_init(n, False)
    parent = base.F1.gen3.make_parent(n, seed=seed)
    csec = np.fft.fft(parent.relation_waves, axis=1) / n
    seed_state = csec[:, 1] / np.linalg.norm(csec[:, 1])
    c2_0 = np.zeros((m, Nn, Neta), complex)
    pump_eta = base.M_PUMP % Neta
    c2_0[:, base.K_PUMP, pump_eta] = z0c
    if delta > 0:
        for k, eta, amplitude in base.CELLS:
            phi = B3_PHASES[(k, eta % Neta)] if (
                arm == ARM_BALANCED and k == 6
            ) else 0.0
            c2_0[:, k, eta % Neta] += phase_factor(phi, amplitude) * seed_state
    p2 = c2_0[:, base.K_PUMP, pump_eta].real
    p2 = p2 / np.linalg.norm(p2)
    q2 = c2_0[:, base.K_PUMP, pump_eta].imag
    q2 = q2 - (q2 @ p2) * p2
    with np.errstate(divide="ignore", invalid="ignore"):
        q2 = q2 / np.linalg.norm(q2)
    engine = base.RecordingEngine(n, c2_0, wp0)
    base._ENGINES.append(engine)
    return engine, p2, q2


def assert_bytes_equal(a: np.ndarray, b: np.ndarray, label: str) -> None:
    if a.dtype != b.dtype or a.shape != b.shape or a.tobytes() != b.tobytes():
        max_abs = None
        if a.shape == b.shape and np.issubdtype(a.dtype, np.number):
            with np.errstate(invalid="ignore"):
                d = np.abs(a - b)
                max_abs = float(np.nanmax(d)) if d.size else 0.0
        raise AssertionError(
            f"byte mismatch {label}: dtype {a.dtype}/{b.dtype}, "
            f"shape {a.shape}/{b.shape}, max_abs={max_abs}"
        )


def compare_npz(left: Path, right: Path) -> dict[str, object]:
    if not left.is_file() or not right.is_file():
        raise FileNotFoundError(f"NPZ comparison missing: {left} / {right}")
    with np.load(left, allow_pickle=False) as a, np.load(right, allow_pickle=False) as b:
        if a.files != b.files:
            raise AssertionError(
                f"NPZ key/order mismatch: {a.files!r} != {b.files!r}"
            )
        for key in a.files:
            assert_bytes_equal(a[key], b[key], key)
        return {
            "left": left.name,
            "right": right.name,
            "array_count": len(a.files),
            "all_arrays_byte_identical": True,
            "max_abs_difference": 0.0,
        }


def compare_npz_reference_subset(reference: Path, candidate: Path) -> dict[str, object]:
    """Require every legacy reference array to be byte-identical in candidate."""
    if not reference.is_file() or not candidate.is_file():
        raise FileNotFoundError(
            f"NPZ subset comparison missing: {reference} / {candidate}"
        )
    with np.load(reference, allow_pickle=False) as ref, \
            np.load(candidate, allow_pickle=False) as cand:
        missing = [key for key in ref.files if key not in cand.files]
        if missing:
            raise AssertionError(f"candidate missing reference keys: {missing!r}")
        for key in ref.files:
            assert_bytes_equal(ref[key], cand[key], key)
        extras = [key for key in cand.files if key not in ref.files]
        return {
            "reference": reference.name,
            "candidate": candidate.name,
            "reference_array_count": len(ref.files),
            "candidate_array_count": len(cand.files),
            "all_reference_arrays_byte_identical": True,
            "max_abs_difference": 0.0,
            "candidate_extra_arrays": extras,
        }


def tagged_paths(arm: str, replicate: str) -> dict[str, Path]:
    suffix = output_suffix(arm, replicate)
    tag = f"_T{FIXED_T}_d{FIXED_DELTA:g}_rep-{suffix}"
    return {
        "result_json": HERE / f"result_nsweep_mixed{tag}_v2.json",
        "npz": HERE / f"nsweep_mixed{tag}_N{FIXED_N}_v2.npz",
        "fig_4panel": HERE / f"fig_mixed{tag}_4panel_N{FIXED_N}_v2.png",
        "fig_mix": HERE / f"fig_mixed{tag}_mix_N{FIXED_N}_v2.png",
        "fig_ledger": HERE / f"fig_mixed{tag}_ledger_N{FIXED_N}_v2.png",
        "fig_summary": HERE / f"fig_mixed{tag}_summary_v2.png",
        "fig_birth_matrix": HERE / f"fig_mixed{tag}_birth_matrix_v2.png",
    }


def canonical_control_npz() -> Path:
    return HERE / (
        f"nsweep_mixed_T{FIXED_T}_d{FIXED_DELTA:g}_N{FIXED_N}_v2.npz"
    )


def manifest_path(arm: str, replicate: str) -> Path:
    arm_tag = "control" if arm == ARM_CONTROL else "phase-balanced"
    return HERE / (
        f"manifest_{arm_tag}_mixed_T{FIXED_T}_d{FIXED_DELTA:g}_"
        f"N{FIXED_N}_rep-{replicate}_v1.json"
    )


def augment_result(path: Path, phase: dict[str, object],
                   hashes: dict[str, str]) -> None:
    result = json.loads(path.read_text(encoding="utf-8"))
    env = result.setdefault("env", {})
    env["phase_control"] = phase
    strength = env.setdefault("seed_strength", {})
    strength["Acoh"] = phase["coherent_sum"]["abs"]
    strength["Acoh_complex"] = phase["coherent_sum"]
    definitions = strength.setdefault("definitions", {})
    definitions["Acoh"] = "abs(sum_j delta_j * exp(i*phi_j))"
    env["source_sha256_phase_run"] = hashes
    path.write_text(
        json.dumps(result, ensure_ascii=False, indent=1, default=float),
        encoding="utf-8",
    )


def acquire_lock() -> tuple[int, Path]:
    lock = HERE / ".phase_balanced_mixed_v1.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(
            f"another phase comparator may be running: {lock.name}"
        ) from exc
    os.write(fd, f"pid={os.getpid()}\n".encode("ascii"))
    return fd, lock


def release_lock(fd: int, lock: Path) -> None:
    os.close(fd)
    try:
        lock.unlink()
    except FileNotFoundError:
        pass


def smoke() -> dict[str, object]:
    """Four deterministic steps, no files and no T=42000 experiment."""
    before = source_hashes()
    base = load_base(SMOKE_STEPS, FIXED_DELTA, "pbmixv1-smoke")
    phase = phase_summary(base, ARM_BALANCED)

    base._ENGINES.clear()
    ref, p_ref, q_ref = base.build_universe(
        FIXED_N, FIXED_DELTA, Nn=FIXED_NN, Neta=FIXED_NETA, seed=2
    )
    ctl, p_ctl, q_ctl = build_phase_universe(
        base, ARM_CONTROL, FIXED_N, FIXED_DELTA,
        Nn=FIXED_NN, Neta=FIXED_NETA, seed=2,
    )
    assert_bytes_equal(ref.C, ctl.C, "control initial C")
    assert_bytes_equal(p_ref, p_ctl, "control p2")
    assert_bytes_equal(q_ref, q_ctl, "control q2")
    for step in range(1, SMOKE_STEPS + 1):
        ref.step()
        ctl.step()
        assert_bytes_equal(ref.C, ctl.C, f"control C step {step}")
    assert_bytes_equal(
        np.asarray(ref._rec), np.asarray(ctl._rec), "control recorder"
    )

    pb1, p1, q1 = build_phase_universe(
        base, ARM_BALANCED, FIXED_N, FIXED_DELTA,
        Nn=FIXED_NN, Neta=FIXED_NETA, seed=2,
    )
    pb2, p2, q2 = build_phase_universe(
        base, ARM_BALANCED, FIXED_N, FIXED_DELTA,
        Nn=FIXED_NN, Neta=FIXED_NETA, seed=2,
    )
    assert_bytes_equal(pb1.C, pb2.C, "phase replicate initial C")
    assert_bytes_equal(p1, p2, "phase replicate p2")
    assert_bytes_equal(q1, q2, "phase replicate q2")
    for step in range(1, SMOKE_STEPS + 1):
        pb1.step()
        pb2.step()
        assert_bytes_equal(pb1.C, pb2.C, f"phase replicate C step {step}")
    assert_bytes_equal(
        np.asarray(pb1._rec), np.asarray(pb2._rec), "phase replicate recorder"
    )

    after = source_hashes()
    if before != after:
        raise AssertionError("source hashes changed during smoke")
    return {
        "status": "PASS",
        "science_run_started": False,
        "N": FIXED_N,
        "steps": SMOKE_STEPS,
        "Nn": FIXED_NN,
        "Neta": FIXED_NETA,
        "delta": FIXED_DELTA,
        "control_matches_base_bytewise": True,
        "phase_replicates_bytewise": True,
        "B3_coherent_sum_abs": phase["B3_coherent_sum"]["abs"],
        "B3_residual_tolerance": phase["B3_residual_tolerance"],
        "source_hashes_unchanged": True,
        "source_sha256": before,
    }


def run_science(arm: str, replicate: str,
                compare_replicate: str | None) -> dict[str, object]:
    before = source_hashes()
    suffix = output_suffix(arm, replicate)
    expected = tagged_paths(arm, replicate)
    occupied = [p.name for p in [*expected.values(), manifest_path(arm, replicate)]
                if p.exists()]
    if occupied:
        raise FileExistsError("refusing to overwrite: " + ", ".join(occupied))

    fd, lock = acquire_lock()
    try:
        base = load_base(FIXED_T, FIXED_DELTA, suffix)
        if (base.ns.T, base.ns.WIN) != (FIXED_T, (FIXED_T // 2, FIXED_T)):
            raise AssertionError(f"time condition changed: {base.ns.T}, {base.ns.WIN}")
        if (base.ns.SEED, tuple(base.ns.CELL), base.ns.ORDER) != (2, (2, 0), 6):
            raise AssertionError("seed/cell/order condition changed")
        phase = phase_summary(base, arm)
        if arm == ARM_BALANCED:
            def builder(n, delta, Nn=5, Neta=8, seed=2):
                return build_phase_universe(
                    base, ARM_BALANCED, n, delta, Nn=Nn, Neta=Neta, seed=seed
                )
            base.ns.F.build_standard_universe = builder
        # For ARM_CONTROL, the existing base builder is intentionally untouched.
        sys.argv = base_argv(FIXED_T, FIXED_DELTA, suffix)
        base.main()
    finally:
        release_lock(fd, lock)

    missing = [p.name for p in expected.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("run did not produce: " + ", ".join(missing))
    after = source_hashes()
    if before != after:
        raise AssertionError("source hashes changed during science run")

    augment_result(expected["result_json"], phase, before)
    comparisons: list[dict[str, object]] = []
    if arm == ARM_CONTROL:
        comparisons.append(compare_npz_reference_subset(
            canonical_control_npz(), expected["npz"]
        ))
    if compare_replicate is not None:
        other = tagged_paths(arm, compare_replicate)["npz"]
        comparisons.append(compare_npz(expected["npz"], other))

    artifact_hashes = {name: sha256(path) for name, path in expected.items()}
    manifest = {
        "schema": "phase-balanced-mixed-manifest-v1",
        "status": "complete",
        "arm": arm,
        "replicate": replicate,
        "conditions": {
            "N": FIXED_N,
            "T": FIXED_T,
            "Nn": FIXED_NN,
            "Neta": FIXED_NETA,
            "delta_per_cell": FIXED_DELTA,
            "seed": 2,
            "cell": [2, 0],
            "order": 6,
            "window": [FIXED_T // 2, FIXED_T],
            "mode": "mixed",
            "F_dynamics": "unified_interaction_v1.py unchanged",
        },
        "phase_control": phase,
        "comparisons": comparisons,
        "source_sha256_before_and_after": before,
        "source_hashes_unchanged": True,
        "artifacts_sha256": artifact_hashes,
        "software": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }
    mpath = manifest_path(arm, replicate)
    mpath.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, default=float),
        encoding="utf-8",
    )
    manifest["manifest_file"] = mpath.name
    manifest["manifest_sha256"] = sha256(mpath)
    return manifest


def validate_existing(arm: str, replicate: str,
                      compare_replicate: str | None) -> dict[str, object]:
    """Validate completed artifacts and create only their missing manifest."""
    expected = tagged_paths(arm, replicate)
    missing = [p.name for p in expected.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("existing run is incomplete: " + ", ".join(missing))
    mpath = manifest_path(arm, replicate)
    if mpath.exists():
        raise FileExistsError(f"refusing to overwrite manifest: {mpath.name}")

    result = json.loads(expected["result_json"].read_text(encoding="utf-8"))
    env = result.get("env", {})
    phase = env.get("phase_control")
    run_hashes = env.get("source_sha256_phase_run")
    if not isinstance(phase, dict) or phase.get("arm") != arm:
        raise AssertionError("existing result lacks matching phase metadata")
    if not isinstance(run_hashes, dict):
        raise AssertionError("existing result lacks run-time source hashes")

    validator_hashes = source_hashes()
    comparisons: list[dict[str, object]] = []
    if arm == ARM_CONTROL:
        comparisons.append(compare_npz_reference_subset(
            canonical_control_npz(), expected["npz"]
        ))
    if compare_replicate is not None:
        comparisons.append(compare_npz(
            expected["npz"], tagged_paths(arm, compare_replicate)["npz"]
        ))

    manifest = {
        "schema": "phase-balanced-mixed-manifest-v1",
        "status": "complete",
        "recovered_by_validation_only": True,
        "trajectory_rerun": False,
        "arm": arm,
        "replicate": replicate,
        "conditions": {
            "N": FIXED_N, "T": FIXED_T, "Nn": FIXED_NN,
            "Neta": FIXED_NETA, "delta_per_cell": FIXED_DELTA,
            "seed": 2, "cell": [2, 0], "order": 6,
            "window": [FIXED_T // 2, FIXED_T], "mode": "mixed",
            "F_dynamics": "unified_interaction_v1.py unchanged",
        },
        "phase_control": phase,
        "comparisons": comparisons,
        "science_run_source_sha256": run_hashes,
        "validation_source_sha256": validator_hashes,
        "artifacts_sha256": {
            name: sha256(path) for name, path in expected.items()
        },
        "software": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }
    mpath.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, default=float),
        encoding="utf-8",
    )
    manifest["manifest_file"] = mpath.name
    manifest["manifest_sha256"] = sha256(mpath)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("smoke", help="4-step read-only-in-memory determinism smoke")
    run = sub.add_parser("run", help="run one fixed T=42000 tagged arm")
    run.add_argument("--arm", choices=ARMS, required=True)
    run.add_argument("--replicate", type=safe_replicate, required=True)
    run.add_argument("--compare-replicate", type=safe_replicate)
    run.add_argument(
        "--allow-science-run",
        action="store_true",
        help="required acknowledgement for the fixed T=42000 run",
    )
    validate = sub.add_parser(
        "validate-existing",
        help="validate seven completed artifacts and create only their manifest",
    )
    validate.add_argument("--arm", choices=ARMS, required=True)
    validate.add_argument("--replicate", type=safe_replicate, required=True)
    validate.add_argument("--compare-replicate", type=safe_replicate)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "smoke":
        out = smoke()
    elif args.command == "validate-existing":
        out = validate_existing(args.arm, args.replicate, args.compare_replicate)
    else:
        if not args.allow_science_run:
            raise SystemExit(
                "T=42000 science run is locked pending implementation review; "
                "pass --allow-science-run only after approval"
            )
        out = run_science(args.arm, args.replicate, args.compare_replicate)
    print(json.dumps(out, ensure_ascii=False, indent=2, default=float))


if __name__ == "__main__":
    main()
