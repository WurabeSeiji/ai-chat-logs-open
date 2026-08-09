#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第3段階・有限位数 root passive probe v1。

母体 ``run_nsweep_three_series_v2.py`` と F v1 は変更しない。派生 subclass が、
各 step の共有線形部を通過した後、非線形頂点が実際に使用する直前の
関係辺別 ``R_e = scale*sin^2(theta_e)``（本条件で ``scale=1``）を
コピーする。通常の母体成果物には probe
配列を混ぜず、別 NPZ/JSON に保存する。

重要な型区別:
  * ``R_e``: ``M=N(N-1)/2`` 本の関係辺に添字を持つ、状態依存の頂点係数。
  * ``r_nopump``: 全関係辺・帯を集約し pump 帯を分母から除いた受動診断量。
  * ``rho``: 有限位数論文の外生的な二チャネル散乱パラメータ。

この runner は ``R_e == rho`` や ``U_rho^n == I`` を仮定しない。Jacobian、
二チャネル射影、monodromy も構成しない。
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
BASE_RUNNER = HERE / "run_nsweep_three_series_v2.py"
F_V1 = HERE.parent / "統一万能関数_v1" / "unified_interaction_v1.py"
DEFAULT_OUTPUT_DIR = HERE / "finite_order_root_probe_stage3_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def delta_tag(value: float) -> str:
    return f"{value:g}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Record the exact pre-vertex relation-edge R_e used by F v1 without "
            "changing the mother trajectory"
        )
    )
    parser.add_argument("--mode", required=True, choices=(
        "mixed", "neutral", "electron", "fermion_family", "boson_family"
    ))
    parser.add_argument("--delta", required=True, type=float)
    parser.add_argument("--N", type=int, default=12)
    parser.add_argument("--T", type=int, default=42000)
    parser.add_argument(
        "--output-suffix",
        required=True,
        help="ASCII identifier passed to the mother runner; existing outputs are refused",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--reference-npz",
        type=Path,
        help=(
            "Optional uninstrumented NPZ. Every common array must match exactly; "
            "a shorter probe is compared with the reference prefix."
        ),
    )
    parser.add_argument(
        "--require-identical-key-set",
        action="store_true",
        help=(
            "Require candidate and reference NPZ key sets to be identical. "
            "Without this flag, an older reference may be a complete subset, "
            "but every reference key must still exist and match."
        ),
    )
    parser.add_argument(
        "--describe-only",
        action="store_true",
        help="Print fixed constants and planned paths without importing/running the mother",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.N < 2:
        raise ValueError("--N must be at least 2")
    if args.T < 1:
        raise ValueError("--T must be positive")
    if not math.isfinite(args.delta) or args.delta < 0.0:
        raise ValueError("--delta must be a finite non-negative number")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.output_suffix):
        raise ValueError("--output-suffix must be 1..64 ASCII letters/digits/_/-")


def base_argv(args: argparse.Namespace) -> list[str]:
    return [
        str(BASE_RUNNER),
        args.mode,
        str(args.N),
        str(args.N),
        str(args.T),
        repr(float(args.delta)),
        args.output_suffix,
    ]


def load_base(argv: list[str]) -> Any:
    spec = importlib.util.spec_from_file_location("stage3_root_probe_mother_v1", BASE_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load mother runner: {BASE_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    old_argv = sys.argv
    try:
        sys.argv = argv
        spec.loader.exec_module(module)
    finally:
        sys.argv = old_argv
    return module


def install_passive_probe(base: Any) -> type:
    """Replace only the derived runner's class binding, never F v1 itself."""

    class PassivePrevertexRecordingEngine(base.RecordingEngine):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self._probe_capture_prevertex = False
            self._probe_R_prevertex: list[np.ndarray] = []
            self._probe_r_nopump_prevertex: list[float] = []
            self._probe_odd_power_prevertex: list[float] = []
            self._probe_even_power_nopump_prevertex: list[float] = []
            self._probe_readout_calls_per_step: list[int] = []
            self._probe_calls_this_step = 0

        def _readout(self):  # type: ignore[no-untyped-def]
            values = super()._readout()
            if self._probe_capture_prevertex:
                copied = np.array(values, dtype=float, copy=True)
                self._probe_R_prevertex.append(copied)
                # Synchronized global diagnostic from the exact same state as
                # copied R_e.  This remains separate from the coefficient used
                # by the vertex and from the ordinary post-step r_nopump.
                power_by_band = np.abs(self.C2()) ** 2
                power_by_band = power_by_band.sum(axis=(0, 2))
                odd_power = float(power_by_band[base.ODD_K].sum())
                even_nonpump_indices = [
                    k for k in base.EVEN_K if k != base.K_PUMP
                ]
                even_nonpump_power = float(
                    power_by_band[even_nonpump_indices].sum()
                )
                denominator = odd_power + even_nonpump_power
                r_nopump = (
                    odd_power / denominator if denominator > 0.0 else float("nan")
                )
                self._probe_r_nopump_prevertex.append(r_nopump)
                self._probe_odd_power_prevertex.append(odd_power)
                self._probe_even_power_nopump_prevertex.append(even_nonpump_power)
                self._probe_capture_prevertex = False
                self._probe_calls_this_step += 1
            return values

        def _nonlinear(self):  # type: ignore[no-untyped-def]
            if self._probe_capture_prevertex:
                raise AssertionError("pre-vertex capture flag leaked from prior step")
            self._probe_capture_prevertex = True
            result = super()._nonlinear()
            if self._probe_capture_prevertex:
                raise AssertionError("F v1 _nonlinear did not call _readout exactly once")
            return result

        def step(self, *args: Any, **kwargs: Any):  # type: ignore[no-untyped-def]
            self._probe_calls_this_step = 0
            result = super().step(*args, **kwargs)
            self._probe_readout_calls_per_step.append(self._probe_calls_this_step)
            if self._probe_calls_this_step != 1:
                raise AssertionError(
                    f"expected one pre-vertex readout, got {self._probe_calls_this_step}"
                )
            return result

    base.RecordingEngine = PassivePrevertexRecordingEngine
    return PassivePrevertexRecordingEngine


def is_declared_time_series_key(key: str) -> bool:
    """Keys whose first dimension is the mother runner's time/sample axis."""

    return key.startswith(("m_", "v_", "rec_m_", "rec_v_"))


def compare_npz_exact_or_prefix(
    candidate: Path,
    reference: Path,
    require_identical_key_set: bool = False,
) -> dict[str, Any]:
    """Compare reference arrays by dtype, shape/prefix, NaN bits, and bytes.

    Prefix comparison is restricted to declared trajectory/ledger keys.  Static
    arrays are never silently shortened.  An older reference may omit arrays
    newly emitted by the current mother runner, but it may not contain a key
    missing from the candidate.  Paired same-source controls should request an
    identical key set.
    """

    mismatches: list[dict[str, Any]] = []
    compared_arrays = 0
    compared_bytes = 0
    with np.load(candidate, allow_pickle=False) as cand, np.load(reference, allow_pickle=False) as ref:
        candidate_keys = set(cand.files)
        reference_keys = set(ref.files)
        common = sorted(candidate_keys & reference_keys)
        missing_reference_keys = sorted(reference_keys - candidate_keys)
        extra_candidate_keys = sorted(candidate_keys - reference_keys)
        for key in common:
            a = np.ascontiguousarray(cand[key])
            b_full = ref[key]
            prefix = False
            if a.shape == b_full.shape:
                b = np.ascontiguousarray(b_full)
            elif (
                a.ndim == b_full.ndim
                and a.ndim >= 1
                and is_declared_time_series_key(key)
                and a.shape[1:] == b_full.shape[1:]
                and a.shape[0] <= b_full.shape[0]
            ):
                b = np.ascontiguousarray(b_full[: a.shape[0]])
                prefix = True
            else:
                mismatches.append({
                    "key": key,
                    "reason": "shape",
                    "candidate_shape": list(a.shape),
                    "reference_shape": list(b_full.shape),
                })
                continue
            if a.dtype != b.dtype:
                mismatches.append({
                    "key": key,
                    "reason": "dtype",
                    "candidate_dtype": str(a.dtype),
                    "reference_dtype": str(b.dtype),
                })
                continue
            compared_arrays += 1
            compared_bytes += int(a.nbytes)
            if a.tobytes(order="C") != b.tobytes(order="C"):
                mismatches.append({
                    "key": key,
                    "reason": "bytes",
                    "prefix_comparison": prefix,
                    "shape": list(a.shape),
                    "nan_layout_equal": bool(
                        np.array_equal(np.isnan(a), np.isnan(b))
                        if np.issubdtype(a.dtype, np.inexact)
                        else True
                    ),
                })
        key_set_ok = not missing_reference_keys and (
            not require_identical_key_set or not extra_candidate_keys
        )
        return {
            "candidate": str(candidate),
            "reference": str(reference),
            "candidate_array_count": len(cand.files),
            "reference_array_count": len(ref.files),
            "common_array_count": len(common),
            "missing_reference_keys": missing_reference_keys,
            "extra_candidate_keys": extra_candidate_keys,
            "key_set_policy": (
                "identical" if require_identical_key_set else "reference_complete_subset"
            ),
            "key_set_ok": key_set_ok,
            "compared_array_count": compared_arrays,
            "compared_bytes": compared_bytes,
            "mismatch_count": len(mismatches),
            "mismatches": mismatches,
            "verdict": "PASS" if common and key_set_ok and not mismatches else "FAIL",
        }


def roots_payload() -> dict[str, Any]:
    rho_124 = math.cos(23.0 * math.pi / 124.0) ** 2
    rho_620 = math.cos(117.0 * math.pi / 620.0) ** 2
    rho_mz = 0.687822933884774
    return {
        "rho_124_23": {
            "value": rho_124,
            "kind": "finite_order_root",
            "n": 124,
            "m": 23,
        },
        "rho_620_117": {
            "value": rho_620,
            "kind": "finite_order_root",
            "n": 620,
            "m": 117,
        },
        "rho_MZ_physical_correspondence": {
            "value": rho_mz,
            "kind": "physical_correspondence_value_not_finite_root",
            "n": None,
            "m": None,
        },
        "separation_rho_620_minus_MZ": rho_620 - rho_mz,
    }


def main() -> None:
    args = parse_args()
    validate_args(args)
    output_dir = args.output_dir.resolve()
    planned = {
        "mode": args.mode,
        "delta": args.delta,
        "N": args.N,
        "T": args.T,
        "output_suffix": args.output_suffix,
        "output_dir": str(output_dir),
        "mother_runner": str(BASE_RUNNER),
        "mother_runner_sha256": sha256(BASE_RUNNER),
        "F_v1": str(F_V1),
        "F_v1_sha256": sha256(F_V1),
        "references": roots_payload(),
    }
    if args.describe_only:
        print(json.dumps(planned, ensure_ascii=False, indent=2))
        return

    reference = None
    if args.reference_npz is not None:
        reference = args.reference_npz.resolve()
        # Check before the expensive mother run.  This also prevents a planned
        # output path from becoming its own "reference" only after generation.
        if not reference.exists():
            raise FileNotFoundError(reference)

    output_dir.mkdir(parents=True, exist_ok=True)
    probe_stem = (
        f"finite_order_root_probe_{args.mode}_N{args.N}_T{args.T}_"
        f"d{delta_tag(args.delta)}_rep-{args.output_suffix}_v1"
    )
    probe_npz = output_dir / f"{probe_stem}.npz"
    probe_json = output_dir / f"{probe_stem}.json"
    if probe_npz.exists() or probe_json.exists():
        raise FileExistsError(f"probe output already exists: {probe_stem}")
    argv = base_argv(args)
    base = load_base(argv)
    install_passive_probe(base)
    # Redirect outputs only after all source imports have resolved from the mother folder.
    base.HERE = output_dir
    base.ns.HERE = output_dir
    ordinary_npz = output_dir / f"nsweep_{args.mode}{base.TAG}_N{args.N}_v2.npz"
    if reference is not None and reference == ordinary_npz.resolve():
        raise ValueError("--reference-npz must not be the candidate mother NPZ itself")

    started = time.time()
    old_argv = sys.argv
    try:
        sys.argv = argv
        base.main()
    finally:
        sys.argv = old_argv

    if len(base._ENGINES) != 2:
        raise AssertionError(f"expected matter+vacuum engines, got {len(base._ENGINES)}")
    matter, vacuum = base._ENGINES
    r_m = np.asarray(matter._probe_R_prevertex, dtype=float)
    r_v = np.asarray(vacuum._probe_R_prevertex, dtype=float)
    rnp_m = np.asarray(matter._probe_r_nopump_prevertex, dtype=float)
    rnp_v = np.asarray(vacuum._probe_r_nopump_prevertex, dtype=float)
    expected_shape = (args.T, args.N * (args.N - 1) // 2)
    if r_m.shape != expected_shape or r_v.shape != expected_shape:
        raise AssertionError(
            f"probe shape mismatch matter={r_m.shape} vacuum={r_v.shape} "
            f"expected={expected_shape}"
        )
    if rnp_m.shape != (args.T,) or rnp_v.shape != (args.T,):
        raise AssertionError(
            f"synchronized r_nopump shape mismatch matter={rnp_m.shape} "
            f"vacuum={rnp_v.shape}"
        )
    if not np.all(np.asarray(matter._probe_readout_calls_per_step) == 1):
        raise AssertionError("matter pre-vertex readout count was not exactly one per step")
    if not np.all(np.asarray(vacuum._probe_readout_calls_per_step) == 1):
        raise AssertionError("vacuum pre-vertex readout count was not exactly one per step")
    if not (np.all(np.isfinite(r_m)) and np.all(np.isfinite(r_v))):
        raise AssertionError("pre-vertex R_e contains non-finite values")
    if not (np.all((0.0 <= r_m) & (r_m <= 1.0)) and np.all((0.0 <= r_v) & (r_v <= 1.0))):
        raise AssertionError("pre-vertex R_e escaped [0,1]")

    if not ordinary_npz.exists():
        raise FileNotFoundError(f"mother NPZ missing after run: {ordinary_npz}")

    parity = None
    if reference is not None:
        parity = compare_npz_exact_or_prefix(
            ordinary_npz,
            reference,
            require_identical_key_set=args.require_identical_key_set,
        )
        if parity["verdict"] != "PASS":
            raise AssertionError(f"non-invasiveness parity failed: {parity['mismatches'][:3]}")

    source_hash_after = {
        "mother_runner_sha256": sha256(BASE_RUNNER),
        "F_v1_sha256": sha256(F_V1),
    }
    if source_hash_after["mother_runner_sha256"] != planned["mother_runner_sha256"]:
        raise AssertionError("mother runner hash changed during probe")
    if source_hash_after["F_v1_sha256"] != planned["F_v1_sha256"]:
        raise AssertionError("F v1 hash changed during probe")

    # Write probe-only artifacts only after all non-invasiveness gates pass.
    np.savez_compressed(
        probe_npz,
        step=np.arange(1, args.T + 1, dtype=np.int64),
        edge_index=np.arange(expected_shape[1], dtype=np.int64),
        edge_ia=np.asarray(matter.ia, dtype=np.int64),
        edge_ib=np.asarray(matter.ib, dtype=np.int64),
        R_prevertex_matter=r_m,
        R_prevertex_vacuum=r_v,
        r_nopump_prevertex_matter=rnp_m,
        r_nopump_prevertex_vacuum=rnp_v,
        odd_power_prevertex_matter=np.asarray(
            matter._probe_odd_power_prevertex, dtype=float
        ),
        odd_power_prevertex_vacuum=np.asarray(
            vacuum._probe_odd_power_prevertex, dtype=float
        ),
        even_power_nopump_prevertex_matter=np.asarray(
            matter._probe_even_power_nopump_prevertex, dtype=float
        ),
        even_power_nopump_prevertex_vacuum=np.asarray(
            vacuum._probe_even_power_nopump_prevertex, dtype=float
        ),
        readout_calls_per_step_matter=np.asarray(
            matter._probe_readout_calls_per_step, dtype=np.int8
        ),
        readout_calls_per_step_vacuum=np.asarray(
            vacuum._probe_readout_calls_per_step, dtype=np.int8
        ),
    )

    payload = {
        "schema": "finite_order_root_probe_stage3_v1",
        "status_boundary": {
            "instrument": (
                "passive pre-vertex relation-edge R_e recorder plus synchronized "
                "global r_nopump diagnostic"
            ),
            "does_not_identify": [
                "R_e with finite-order-paper rho",
                "r_nopump with rho",
                "instantaneous root proximity with U_rho^n=I",
            ],
            "jacobian_or_monodromy": "not constructed",
        },
        "conditions": planned,
        "source_hash_after": source_hash_after,
        "source_hash_unchanged": True,
        "timing_definition": (
            "R_e copied inside F v1 _nonlinear, after _linear and immediately when "
            "the parent _readout returns the coefficient consumed by _vertex_rate; "
            "r_nopump_prevertex is calculated from that same unchanged state"
        ),
        "index_definition": {
            "R_e_index": "relation edge e, not graph node",
            "M": expected_shape[1],
            "edge_endpoints": "edge_ia[e], edge_ib[e]",
            "readout_aggregation": (
                "edge power plus powers of other edges adjacent to both endpoints; "
                "the returned array remains indexed by relation edge"
            ),
        },
        "probe_npz": probe_npz.name,
        "ordinary_mother_npz": ordinary_npz.name,
        "probe_shapes": {
            "R_prevertex_matter": list(r_m.shape),
            "R_prevertex_vacuum": list(r_v.shape),
            "r_nopump_prevertex_matter": list(rnp_m.shape),
            "r_nopump_prevertex_vacuum": list(rnp_v.shape),
        },
        "probe_ranges": {
            "matter_min": float(r_m.min()),
            "matter_max": float(r_m.max()),
            "vacuum_min": float(r_v.min()),
            "vacuum_max": float(r_v.max()),
        },
        "one_prevertex_readout_per_step": True,
        "noninvasiveness_parity": parity,
        "runtime_sec": time.time() - started,
    }
    probe_json.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"saved passive probe: {probe_npz}")
    print(f"saved metadata: {probe_json}")
    if parity is not None:
        print(
            "non-invasiveness: "
            f"{parity['verdict']} / {parity['compared_array_count']} arrays / "
            f"{parity['compared_bytes']} bytes"
        )


if __name__ == "__main__":
    main()
