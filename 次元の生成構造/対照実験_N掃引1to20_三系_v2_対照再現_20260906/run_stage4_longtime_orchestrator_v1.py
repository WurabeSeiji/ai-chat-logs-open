#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第4段階長時間走行を、再現対照・prefix byte照合のゲート下で逐次実行する。

既存 run_nsweep_three_series_v2.py / F v1 は変更しない。全runに
固有output_suffixを付け、既存・部分成果物を上書きしない。プロセスは常に
1本ずつ。完成成果物は再開時に再検証してskipする。

実行順:
  Stage N4
    1. mixed / neutral / electron の T=4000 exact control
    2. fermion_family (F5) の新規 T=4000 control
    3. 同4 modeの T=42000。各first4000を上記controlとbyte照合
  Stage N12
    4. neutral / electron / fermion_family / mixed の delta=.01, T=300000
    5. 各既存T42000正本、および存在する独立T4000正本とprefix照合

N12 F5に独立T4000正本は存在しない。既存T42000正本と42000-step
prefix一致を必須とし、その内側first4000を独立証拠として二重計上しない。

安全のため引数なしでは実行しない。

  python3 run_stage4_longtime_orchestrator_v1.py --list
  python3 run_stage4_longtime_orchestrator_v1.py --smoke
  python3 run_stage4_longtime_orchestrator_v1.py --run-n4
  python3 run_stage4_longtime_orchestrator_v1.py --run-n12
  python3 run_stage4_longtime_orchestrator_v1.py --run-all
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterator

import numpy as np


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_nsweep_three_series_v2.py"
BASE = HERE / "run_tb_nsweep_1to20_v1.py"
UF = HERE.parent / "統一万能関数_v1"
MANIFEST = HERE / "manifest_stage4_longtime_v1.json"
LOCK = HERE / ".stage4_longtime_v1.lock"
MPL_CACHE = Path(tempfile.gettempdir()) / "stage4_longtime_mplconfig_v1"

DELTA = 1e-2
NN = 16
NETA = 8
PARENT_SEED = 2
CELL = (2, 0)
ORDER = 6
LEDGER_EVERY = 50
EXPECTED_ARRAY_COUNT = 99
DISK_SAFETY_FACTOR = 3.0
DISK_HEADROOM_BYTES = 2 * 1024 ** 3
RAM_SAFETY_FACTOR = 1.25
RAM_PHYSICAL_FACTOR = 1.50
ARTIFACT_OVERHEAD_BYTES = 12 * 1024 ** 2

EXPECTED_SOURCE_SHA256 = {
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

N4_CANONICAL = {
    "mixed": HERE / "nsweep_mixed_N4_v2.npz",
    "neutral": HERE / "nsweep_neutral_N4_v2.npz",
    "electron": HERE / "nsweep_electron_N4_v2.npz",
}
N4_VACUUM = HERE / "nsweep_vacuum_N4_v2.npz"
N12_SHORT_REFERENCES = {
    "mixed": HERE / "nsweep_mixed_N12_v2.npz",
    "neutral": HERE / "nsweep_neutral_N12_v2.npz",
    "electron": HERE / "nsweep_electron_N12_v2.npz",
}
N12_LONG_REFERENCES = {
    "mixed": HERE / "nsweep_mixed_T42000_N12_v2.npz",
    "neutral": HERE / "nsweep_neutral_T42000_N12_v2.npz",
    "electron": HERE / "nsweep_electron_T42000_d0.01_N12_v2.npz",
    "fermion_family": HERE / "nsweep_fermion_family_T42000_d0.01_N12_v2.npz",
}

REFERENCE_JSON_CONDITIONS = {
    HERE / "result_nsweep_mixed_v2.json": ("mixed", 4000, (4, 12)),
    HERE / "result_nsweep_neutral_v2.json": ("neutral", 4000, (4, 12)),
    HERE / "result_nsweep_electron_v2.json": ("electron", 4000, (4, 12)),
    HERE / "result_nsweep_mixed_T42000_v2.json": ("mixed", 42000, (12,)),
    HERE / "result_nsweep_neutral_T42000_v2.json": ("neutral", 42000, (12,)),
    HERE / "result_nsweep_electron_T42000_d0.01_v2.json":
        ("electron", 42000, (12,)),
    HERE / "result_nsweep_fermion_family_T42000_d0.01_v2.json":
        ("fermion_family", 42000, (12,)),
}

EXPECTED_REFERENCE_SHA256 = {
    N4_CANONICAL["mixed"]:
        "f52ce29ffbd2a550572e33b247c185dec4b5c0826dd8b810e17ec96832a465b9",
    N4_CANONICAL["neutral"]:
        "85f8d650ca10000e53e191c48c3dab957ebcb1055afab0e25cd4dd6830d26c5e",
    N4_CANONICAL["electron"]:
        "ba769fcf2dacdcd4ac26959c32c91ef4530a6f0556a3087b648555f1e547b416",
    N4_VACUUM:
        "77773a774f4a2a69333daf90f9a5b6f27c2497998c5bf06b3a4715ab2b3fbf25",
    N12_SHORT_REFERENCES["mixed"]:
        "0b114f10aacb892f586dc9174d772849ba3e35d48ceda2b00dfe7f1f41e5b05e",
    N12_SHORT_REFERENCES["neutral"]:
        "bc04e8beb99a796ed5316f3323b62739260fcc5683ef2705a449f5bd94949cb5",
    N12_SHORT_REFERENCES["electron"]:
        "e2fd0a9c37024e0e6af54223a4b9df088094fede59f95070e8359aea9d2d9747",
    N12_LONG_REFERENCES["mixed"]:
        "0de730a69d303fa7af0a88f09d768bea6f07159348a5c5b4b097cf60c1ad01ff",
    N12_LONG_REFERENCES["neutral"]:
        "6c78a39c37e50243ef78383d7ec8d94e7f741c207d5700d29bc55508f3be33c0",
    N12_LONG_REFERENCES["electron"]:
        "4d558b3caa07a3943f7a6418d2c515d9d7fded0d3c1f06940caa0dc7d1beb73a",
    N12_LONG_REFERENCES["fermion_family"]:
        "752fad29af40f896973acd844e5082fd7ff06df7cab6436f476dff4560aa7a9a",
    HERE / "result_nsweep_mixed_v2.json":
        "9092b2f262440915d7aff32702203e5ad5be2b05074cc0e95dcff4ebcbaed5b6",
    HERE / "result_nsweep_neutral_v2.json":
        "7c27fd7a849fc63649418c02f6fbaa21b60f86752a7f074bd55477e5c8883a70",
    HERE / "result_nsweep_electron_v2.json":
        "de50d6bcfaa5c01d125b172dd3af99b631423fbfc2d49835efedded42bfee6ec",
    HERE / "result_nsweep_mixed_T42000_v2.json":
        "af756ff7992e2c8342bd884d110fcd5d6c6cabcfb65e0f6b05b3095777934864",
    HERE / "result_nsweep_neutral_T42000_v2.json":
        "225ff0ec279090a6321cb3f2e05bd8c6386209354a0566ce4c44fb8941150c66",
    HERE / "result_nsweep_electron_T42000_d0.01_v2.json":
        "90ef2f194b6b0f88973f08c5994e77d28bbe2e87284a4a1b513a236126093e95",
    HERE / "result_nsweep_fermion_family_T42000_d0.01_v2.json":
        "1e3b284f7fccf42c7139eb7c4ad72c8df6157bd45bc7da4f9f02bf8c1446136c",
}


def mode_label(mode: str) -> str:
    return "F5" if mode == "fermion_family" else mode


@dataclass(frozen=True)
class RunSpec:
    run_id: str
    stage: str
    role: str
    mode: str
    n: int
    steps: int
    suffix: str

    @property
    def time_tag(self) -> str:
        return "" if self.steps == 4000 else f"_T{self.steps}"

    @property
    def tag(self) -> str:
        return f"{self.time_tag}_d{DELTA:g}_rep-{self.suffix}"

    @property
    def command(self) -> list[str]:
        return [
            sys.executable, "-u", str(RUNNER), self.mode,
            str(self.n), str(self.n), str(self.steps), repr(DELTA), self.suffix,
        ]

    @property
    def npz_path(self) -> Path:
        return HERE / f"nsweep_{self.mode}{self.tag}_N{self.n}_v2.npz"

    @property
    def json_path(self) -> Path:
        return HERE / f"result_nsweep_{self.mode}{self.tag}_v2.json"

    @property
    def figure_paths(self) -> tuple[Path, ...]:
        stem = f"fig_{self.mode}{self.tag}"
        return (
            HERE / f"{stem}_4panel_N{self.n}_v2.png",
            HERE / f"{stem}_mix_N{self.n}_v2.png",
            HERE / f"{stem}_ledger_N{self.n}_v2.png",
            HERE / f"{stem}_summary_v2.png",
            HERE / f"{stem}_birth_matrix_v2.png",
        )

    @property
    def artifacts(self) -> tuple[Path, ...]:
        return (self.json_path, self.npz_path, *self.figure_paths)

    @property
    def log_path(self) -> Path:
        return HERE / f"stage4_log_{self.run_id}_v1.txt"

    @property
    def all_outputs(self) -> tuple[Path, ...]:
        return (*self.artifacts, self.log_path)


N4_CONTROLS = tuple(
    RunSpec(
        f"n4_{mode_label(mode).lower()}_t4000_control", "N4", "control",
        mode, 4, 4000, f"s4-n4-{mode_label(mode).lower()}-t4000-c1",
    )
    for mode in ("mixed", "neutral", "electron", "fermion_family")
)
N4_LONG = tuple(
    RunSpec(
        f"n4_{mode_label(mode).lower()}_t42000_long", "N4", "long",
        mode, 4, 42000, f"s4-n4-{mode_label(mode).lower()}-t42000-l1",
    )
    for mode in ("mixed", "neutral", "electron", "fermion_family")
)
N12_RUNS = tuple(
    RunSpec(
        f"n12_{mode_label(mode).lower()}_t300000_long", "N12", "long",
        mode, 12, 300000, f"s4-n12-{mode_label(mode).lower()}-t300000-l1",
    )
    for mode in ("neutral", "electron", "fermion_family", "mixed")
)
ALL_SPECS = (*N4_CONTROLS, *N4_LONG, *N12_RUNS)
N4_CONTROL_BY_MODE = {s.mode: s for s in N4_CONTROLS}


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(HERE))
    except ValueError:
        return str(path)


def verify_hashes(expected: dict[Path, str], label: str) -> dict[str, str]:
    actual: dict[str, str] = {}
    errors: list[str] = []
    for path, wanted in expected.items():
        if not path.is_file():
            errors.append(f"missing {display_path(path)}")
            continue
        got = sha256(path)
        actual[display_path(path)] = got
        if got != wanted:
            errors.append(
                f"{display_path(path)}: expected {wanted}, got {got}"
            )
    if errors:
        raise RuntimeError(f"{label} SHA-256 mismatch: " + " | ".join(errors))
    return actual


def array_bytes(a: np.ndarray) -> bytes:
    return np.ascontiguousarray(a).tobytes(order="C")


def array_class(key: str, shape: tuple[int, ...], reference_steps: int) -> str:
    ledger_count = reference_steps // LEDGER_EVERY + 1
    if "ledger" in key and shape and shape[0] == ledger_count:
        return "ledger"
    if shape and shape[0] == reference_steps:
        return "time_series"
    return "static"


def compare_npz(
    reference: Path,
    candidate: Path,
    *,
    reference_steps: int,
    prefix: bool,
    keys: list[str] | None = None,
) -> dict[str, Any]:
    """referenceの全指定keyをcandidateとdtype/shape/byte exact照合。"""
    counts = {"time_series": 0, "ledger": 0, "static": 0}
    with np.load(reference, allow_pickle=False) as ref, \
            np.load(candidate, allow_pickle=False) as cand:
        names = list(ref.files) if keys is None else list(keys)
        missing = [key for key in names if key not in cand.files]
        if missing:
            raise AssertionError(f"{candidate.name}: missing keys {missing}")
        for key in names:
            x = np.asarray(ref[key])
            y = np.asarray(cand[key])
            kind = array_class(key, x.shape, reference_steps)
            if prefix and kind in {"time_series", "ledger"} and x.shape != y.shape:
                if (
                    x.ndim < 1 or y.ndim != x.ndim
                    or x.shape[1:] != y.shape[1:]
                    or y.shape[0] < x.shape[0]
                ):
                    raise AssertionError(
                        f"{key}: prefix-incompatible {x.shape} vs {y.shape}"
                    )
                y = y[:x.shape[0]]
            if x.dtype != y.dtype or x.shape != y.shape:
                raise AssertionError(
                    f"{key}: dtype/shape {x.dtype}/{x.shape} != {y.dtype}/{y.shape}"
                )
            if array_bytes(x) != array_bytes(y):
                raise AssertionError(f"{key}: bytes differ")
            counts[kind] += 1
    return {
        "reference": reference.name,
        "candidate": candidate.name,
        "prefix": prefix,
        "reference_steps": reference_steps,
        "arrays_checked": sum(counts.values()),
        "array_classes": counts,
        "byte_exact": True,
        "max_abs_difference": 0.0,
    }


def vacuum_state_keys(reference: Path) -> list[str]:
    state_rec = {
        "r_mean", "r_med", "r_min", "r_max", "r_raw",
        "dist_alpha", "absdist_alpha", "r_nopump", "dist_alpha_nopump",
        "even_power_nopump", "odd_power", "odd_amp_max", "even_power",
        "total_power", "bands", "ledger", "ledger_t",
    }
    with np.load(reference, allow_pickle=False) as z:
        return [
            key for key in z.files
            if key.startswith("v_")
            or (
                key.startswith("rec_v_")
                and key.removeprefix("rec_v_") in state_rec
            )
        ]


def validate_reference_jsons() -> dict[str, Any]:
    checked: dict[str, Any] = {}
    for path, (mode, steps, ns) in REFERENCE_JSON_CONDITIONS.items():
        data = json.loads(path.read_text(encoding="utf-8"))
        env = data.get("env", {})
        expected = {
            "mode": mode, "Nn": NN, "Neta": NETA, "T": steps,
            "delta": DELTA, "seed": PARENT_SEED, "cell": list(CELL),
            "order": ORDER, "window": [steps // 2, steps],
            "output_suffix": None,
        }
        for key, wanted in expected.items():
            if env.get(key) != wanted:
                raise AssertionError(
                    f"{path.name}: env.{key}={env.get(key)!r}, expected {wanted!r}"
                )
        for n in ns:
            rec = data.get("N", {}).get(str(n))
            if not rec or not rec.get("built"):
                raise AssertionError(f"{path.name}: N={n} is not built")
        checked[path.name] = {
            "mode": mode, "T": steps, "built_N": list(ns), "env_fixed": True,
        }
    return checked


def validate_reference_npz_shapes() -> dict[str, Any]:
    specs: list[tuple[Path, int]] = [
        *((path, 4000) for path in N4_CANONICAL.values()),
        (N4_VACUUM, 4000),
        *((path, 4000) for path in N12_SHORT_REFERENCES.values()),
        *((path, 42000) for path in N12_LONG_REFERENCES.values()),
    ]
    out: dict[str, Any] = {}
    for path, steps in specs:
        with np.load(path, allow_pickle=False) as z:
            for key in ("m_f2", "v_f2", "rec_m_r_mean", "rec_v_r_mean"):
                if key not in z.files or z[key].shape[0] != steps:
                    raise AssertionError(f"{path.name}: {key} length != {steps}")
            out[path.name] = {"arrays": len(z.files), "time_steps": steps}
    return out


def validate_existing_reference_prefixes() -> dict[str, Any]:
    """N12独立T4000正本は既存T42000軌道のprefixと一致するか。"""
    out: dict[str, Any] = {}
    for mode, short in N12_SHORT_REFERENCES.items():
        out[mode] = compare_npz(
            short, N12_LONG_REFERENCES[mode], reference_steps=4000, prefix=True,
        )
    out["fermion_family"] = {
        "independent_T4000_reference": None,
        "counted_as_independent_evidence": False,
        "reason": (
            "N12 F5 has no independent T4000 artifact; first4000 inside the "
            "T42000 reference is not counted twice"
        ),
    }
    return out


def validate_json(spec: RunSpec) -> dict[str, Any]:
    data = json.loads(spec.json_path.read_text(encoding="utf-8"))
    env = data.get("env", {})
    expected = {
        "mode": spec.mode, "Nn": NN, "Neta": NETA, "T": spec.steps,
        "delta": DELTA, "seed": PARENT_SEED, "cell": list(CELL),
        "order": ORDER, "window": [spec.steps // 2, spec.steps],
        "output_suffix": spec.suffix,
    }
    for key, wanted in expected.items():
        if env.get(key) != wanted:
            raise AssertionError(
                f"{spec.json_path.name}: env.{key}={env.get(key)!r}, expected {wanted!r}"
            )
    rec = data.get("N", {}).get(str(spec.n))
    if not rec or not rec.get("built"):
        raise AssertionError(f"{spec.json_path.name}: N={spec.n} is not built")
    if data.get("failed_N") != []:
        raise AssertionError(f"{spec.json_path.name}: failed_N={data.get('failed_N')}")
    return {"env_fixed": True, "built_N": spec.n, "failed_N": []}


def validate_npz(spec: RunSpec) -> dict[str, Any]:
    with np.load(spec.npz_path, allow_pickle=False) as z:
        if len(z.files) != EXPECTED_ARRAY_COUNT:
            raise AssertionError(
                f"{spec.npz_path.name}: arrays={len(z.files)}, expected {EXPECTED_ARRAY_COUNT}"
            )
        for key in ("m_f2", "v_f2", "rec_m_r_mean", "rec_v_r_mean"):
            if key not in z.files or z[key].shape[0] != spec.steps:
                raise AssertionError(
                    f"{spec.npz_path.name}: {key} length != {spec.steps}"
                )
        expected_t = np.asarray(
            [1, *range(LEDGER_EVERY, spec.steps + 1, LEDGER_EVERY)], dtype=float,
        )
        for key in ("rec_m_ledger_t", "rec_v_ledger_t"):
            got = np.asarray(z[key])
            if got.dtype != expected_t.dtype or array_bytes(got) != array_bytes(expected_t):
                raise AssertionError(f"{spec.npz_path.name}: {key} cadence mismatch")
        return {
            "arrays": len(z.files), "time_steps": spec.steps,
            "ledger_snapshots": int(expected_t.size),
        }


def artifact_hashes(spec: RunSpec) -> dict[str, dict[str, Any]]:
    return {
        path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
        for path in spec.all_outputs
    }


def assert_stored_artifact_hashes(
    spec: RunSpec,
    previous_entry: dict[str, Any],
    current_validation: dict[str, Any],
) -> None:
    """完成後の成果物変更をresume時に新しい正本として再基準化しない。"""
    stored = previous_entry.get("validation", {}).get("artifact_hashes")
    if stored is None:
        return
    current = current_validation.get("artifact_hashes")
    if stored != current:
        names = sorted(set(stored) | set(current or {}))
        changed = [
            name for name in names
            if stored.get(name) != (current or {}).get(name)
        ]
        raise AssertionError(
            f"{spec.run_id}: completed artifact hash drift: {changed}"
        )


def validate_spec(spec: RunSpec) -> dict[str, Any]:
    missing = [path.name for path in spec.all_outputs if not path.is_file()]
    if missing:
        raise AssertionError(f"{spec.run_id}: missing outputs {missing}")
    validation: dict[str, Any] = {
        "json": validate_json(spec),
        "npz": validate_npz(spec),
    }
    if spec.stage == "N4" and spec.role == "control":
        if spec.mode in N4_CANONICAL:
            validation["existing_T4000_exact_control"] = compare_npz(
                N4_CANONICAL[spec.mode], spec.npz_path,
                reference_steps=4000, prefix=False,
            )
        else:
            validation["new_T4000_control"] = {
                "mode": "fermion_family", "independent_existing_reference": None,
                "baseline_for_N4_T42000": True,
            }
        validation["internal_vacuum_exact"] = compare_npz(
            N4_VACUUM, spec.npz_path, reference_steps=4000, prefix=False,
            keys=vacuum_state_keys(N4_VACUUM),
        )
    elif spec.stage == "N4" and spec.role == "long":
        control = N4_CONTROL_BY_MODE[spec.mode]
        validation["first4000_vs_stage_control"] = compare_npz(
            control.npz_path, spec.npz_path, reference_steps=4000, prefix=True,
        )
    elif spec.stage == "N12":
        validation["first42000_vs_existing"] = compare_npz(
            N12_LONG_REFERENCES[spec.mode], spec.npz_path,
            reference_steps=42000, prefix=True,
        )
        if spec.mode in N12_SHORT_REFERENCES:
            validation["first4000_vs_independent_existing"] = compare_npz(
                N12_SHORT_REFERENCES[spec.mode], spec.npz_path,
                reference_steps=4000, prefix=True,
            )
        else:
            validation["first4000_independent_evidence"] = {
                "available": False,
                "counted": False,
                "reason": (
                    "No independent N12 F5 T4000 artifact. The first4000 already "
                    "contained in the required T42000 prefix is not double-counted."
                ),
            }
    else:
        raise AssertionError(f"unknown stage/role: {spec.stage}/{spec.role}")
    validation["artifact_hashes"] = artifact_hashes(spec)
    return validation


def output_state(spec: RunSpec) -> str:
    present = [path for path in spec.all_outputs if path.exists()]
    if not present:
        return "absent"
    if len(present) == len(spec.all_outputs) and all(path.is_file() for path in present):
        return "complete"
    raise FileExistsError(
        f"{spec.run_id}: partial outputs; refusing overwrite: "
        + ", ".join(path.name for path in present)
    )


def estimate_array_bytes(reference: Path, target_steps: int) -> tuple[int, int, float]:
    reference_steps = 42000
    reference_ledger = reference_steps // LEDGER_EVERY + 1
    target_ledger = target_steps // LEDGER_EVERY + 1
    with np.load(reference, allow_pickle=False) as z:
        source_raw = sum(int(z[key].nbytes) for key in z.files)
        target_raw = 0
        for key in z.files:
            a = z[key]
            if a.ndim and a.shape[0] == reference_steps:
                target_raw += math.ceil(a.nbytes * target_steps / reference_steps)
            elif a.ndim and "ledger" in key and a.shape[0] == reference_ledger:
                target_raw += math.ceil(a.nbytes * target_ledger / reference_ledger)
            else:
                target_raw += int(a.nbytes)
    compression_ratio = reference.stat().st_size / source_raw
    # 旧schema(80 keys)の欠落分とPython一時配列を過小評価しないよう、
    # 99-key正本かどうかに関わらずrawに25%のschema headroomを持たせる。
    estimated_raw = math.ceil(target_raw * 1.25)
    estimated_npz = math.ceil(estimated_raw * compression_ratio)
    return estimated_raw, estimated_npz, compression_ratio


def runtime_baseline_seconds(mode: str) -> float:
    path_by_mode = {
        "mixed": HERE / "result_nsweep_mixed_T42000_v2.json",
        "neutral": HERE / "result_nsweep_neutral_T42000_v2.json",
        "electron": HERE / "result_nsweep_electron_T42000_d0.01_v2.json",
        "fermion_family": HERE / "result_nsweep_fermion_family_T42000_d0.01_v2.json",
    }
    value = json.loads(path_by_mode[mode].read_text(encoding="utf-8")).get("runtime_sec")
    return float(value) if value is not None else 500.0


def estimate_spec(spec: RunSpec) -> dict[str, Any]:
    raw, npz_bytes, ratio = estimate_array_bytes(
        N12_LONG_REFERENCES[spec.mode], spec.steps
    )
    output_bytes = npz_bytes + ARTIFACT_OVERHEAD_BYTES
    peak_ram = raw * 8 + 1024 ** 3
    runtime = runtime_baseline_seconds(spec.mode) * spec.steps / 42000
    return {
        "run_id": spec.run_id,
        "raw_npz_arrays_bytes": raw,
        "estimated_compressed_npz_bytes": npz_bytes,
        "reference_compression_ratio": ratio,
        "estimated_all_output_bytes": output_bytes,
        "estimated_peak_ram_bytes": peak_ram,
        "estimated_runtime_sec": runtime,
        "conservative_notes": [
            "compressed estimate includes 25% schema margin",
            "peak RAM = raw NPZ estimate x8 + 1 GiB for Python lists and temporaries",
            "N4 estimate intentionally uses N12 reference and is conservative",
        ],
    }


def physical_memory_bytes() -> int | None:
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (ValueError, OSError, TypeError):
        return None


def available_memory_bytes() -> int | None:
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        match = re.search(r"^MemAvailable:\s+(\d+)\s+kB$", meminfo.read_text(), re.M)
        if match:
            return int(match.group(1)) * 1024
    try:
        proc = subprocess.run(
            ["vm_stat"], check=True, capture_output=True, text=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    size_match = re.search(r"page size of (\d+) bytes", proc.stdout)
    if not size_match:
        return None
    page_size = int(size_match.group(1))
    pages = 0
    for label in ("Pages free", "Pages inactive", "Pages speculative", "Pages purgeable"):
        match = re.search(rf"^{re.escape(label)}:\s+(\d+)\.", proc.stdout, re.M)
        if match:
            pages += int(match.group(1))
    return pages * page_size if pages else None


def capacity_report(specs: tuple[RunSpec, ...] | list[RunSpec]) -> dict[str, Any]:
    estimates = [estimate_spec(spec) for spec in specs]
    remaining = []
    for spec, estimate in zip(specs, estimates):
        state = output_state(spec)
        if state == "absent":
            remaining.append(estimate)
    predicted = sum(x["estimated_all_output_bytes"] for x in remaining)
    peak = max((x["estimated_peak_ram_bytes"] for x in remaining), default=0)
    disk = shutil.disk_usage(HERE)
    required_disk = math.ceil(predicted * DISK_SAFETY_FACTOR) + DISK_HEADROOM_BYTES
    physical = physical_memory_bytes()
    available = available_memory_bytes()
    disk_pass = disk.free >= required_disk
    physical_pass = physical is None or physical >= math.ceil(peak * RAM_PHYSICAL_FACTOR)
    available_pass = available is None or available >= math.ceil(peak * RAM_SAFETY_FACTOR)
    return {
        "spec_estimates": estimates,
        "remaining_run_count": len(remaining),
        "estimated_remaining_output_bytes": predicted,
        "estimated_peak_ram_bytes": peak,
        "disk_free_bytes": disk.free,
        "disk_required_bytes": required_disk,
        "disk_safety_factor": DISK_SAFETY_FACTOR,
        "disk_headroom_bytes": DISK_HEADROOM_BYTES,
        "physical_memory_bytes": physical,
        "available_memory_bytes": available,
        "ram_available_safety_factor": RAM_SAFETY_FACTOR,
        "ram_physical_safety_factor": RAM_PHYSICAL_FACTOR,
        "disk_gate_pass": disk_pass,
        "physical_ram_gate_pass": physical_pass,
        "available_ram_gate_pass": available_pass,
        "gate_pass": disk_pass and physical_pass and available_pass,
        "estimated_runtime_sec": sum(x["estimated_runtime_sec"] for x in remaining),
    }


def enforce_capacity(report: dict[str, Any]) -> None:
    if not report["disk_gate_pass"]:
        raise RuntimeError(
            f"disk capacity gate failed: free={report['disk_free_bytes']}, "
            f"required={report['disk_required_bytes']}"
        )
    if not report["physical_ram_gate_pass"]:
        raise RuntimeError(
            f"physical RAM gate failed: physical={report['physical_memory_bytes']}, "
            f"estimated_peak={report['estimated_peak_ram_bytes']}"
        )
    if not report["available_ram_gate_pass"]:
        raise RuntimeError(
            f"available RAM gate failed: available={report['available_memory_bytes']}, "
            f"estimated_peak={report['estimated_peak_ram_bytes']}"
        )


def new_manifest() -> dict[str, Any]:
    return {
        "schema": "stage4-longtime-orchestrator-v1",
        "created_at": now_iso(),
        "conditions": {
            "delta": DELTA, "Nn": NN, "Neta": NETA,
            "parent_seed": PARENT_SEED, "cell": list(CELL), "order": ORDER,
            "sequential": True, "all_runs_have_unique_output_suffix": True,
            "N4": {
                "controls": 4, "long_runs": 4,
                "T_control": 4000, "T_long": 42000,
                "modes": ["mixed", "neutral", "electron", "fermion_family"],
            },
            "N12": {
                "long_runs": 4, "T": 300000,
                "modes": ["neutral", "electron", "fermion_family", "mixed"],
            },
        },
        "evidence_policy": {
            "N4": (
                "existing T4000 exact controls for mixed/neutral/electron; new F5 "
                "T4000 control; every T42000 first4000 matches its stage control"
            ),
            "N12": (
                "every T300000 first42000 matches the existing T42000 reference; "
                "mixed/neutral/electron also match independent existing T4000 references"
            ),
            "N12_F5_T4000": {
                "independent_reference_exists": False,
                "double_count_first4000_inside_T42000": False,
                "required_evidence": "existing F5 T42000 full 42000-step prefix",
            },
        },
        "plan": [
            {
                "run_id": spec.run_id,
                "stage": spec.stage,
                "role": spec.role,
                "mode": spec.mode,
                "N": spec.n,
                "T": spec.steps,
                "delta": DELTA,
                "suffix": spec.suffix,
                "command": spec.command,
                "artifacts": [path.name for path in spec.artifacts],
                "log": spec.log_path.name,
            }
            for spec in ALL_SPECS
        ],
        "runs": {},
    }


def load_manifest() -> dict[str, Any]:
    if not MANIFEST.is_file():
        return new_manifest()
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if data.get("schema") != "stage4-longtime-orchestrator-v1":
        raise RuntimeError(f"unexpected manifest schema: {data.get('schema')}")
    current_plan = new_manifest()["plan"]
    if data.get("plan") != current_plan:
        raise RuntimeError(
            "stored Stage 4 plan differs from the current interpreter/path/schedule"
        )
    return data


def save_manifest(manifest: dict[str, Any]) -> None:
    manifest["updated_at"] = now_iso()
    fd, temp_name = tempfile.mkstemp(
        prefix=".manifest_stage4_longtime_v1.", suffix=".tmp", dir=HERE,
    )
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=1, ensure_ascii=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(temp_path, MANIFEST)
    finally:
        if temp_path.exists():
            temp_path.unlink()


@contextmanager
def campaign_lock() -> Iterator[None]:
    with LOCK.open("a+", encoding="utf-8") as lock_file:
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("another Stage 4 orchestrator is running") from exc
        lock_file.seek(0)
        lock_file.truncate()
        lock_file.write(f"pid={os.getpid()} started={now_iso()}\n")
        lock_file.flush()
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def execute_subprocess(spec: RunSpec) -> int:
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(MPL_CACHE)
    with spec.log_path.open("x", encoding="utf-8") as log:
        proc = subprocess.Popen(
            spec.command, cwd=HERE, env=env, shell=False,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                print(line, end="", flush=True)
                log.write(line)
                log.flush()
            return proc.wait()
        except BaseException:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
            raise


def execute_one(
    spec: RunSpec,
    manifest: dict[str, Any],
    index: int,
    total: int,
) -> None:
    print(f"\n===== [{index}/{total}] {spec.run_id} =====", flush=True)
    state = output_state(spec)
    if state == "complete":
        entry = manifest["runs"].setdefault(spec.run_id, {})
        validation = validate_spec(spec)
        assert_stored_artifact_hashes(spec, entry, validation)
        entry.update({
            "status": "validated_existing", "validated_at": now_iso(),
            "command": spec.command, "validation": validation,
        })
        save_manifest(manifest)
        print("完成成果物を再検証してskip", flush=True)
        return

    per_run_capacity = capacity_report([spec])
    enforce_capacity(per_run_capacity)
    source_before = verify_hashes(EXPECTED_SOURCE_SHA256, "source before run")
    reference_before = verify_hashes(
        EXPECTED_REFERENCE_SHA256, "reference before run"
    )
    started = time.time()
    entry: dict[str, Any] = {
        "status": "running", "started_at": now_iso(),
        "command": spec.command, "command_display": shlex.join(spec.command),
        "log": spec.log_path.name,
        "capacity_gate_before": per_run_capacity,
        "source_sha256_before": source_before,
        "reference_sha256_before": reference_before,
    }
    manifest["runs"][spec.run_id] = entry
    save_manifest(manifest)
    try:
        returncode = execute_subprocess(spec)
        entry["returncode"] = returncode
        entry["runtime_sec"] = time.time() - started
        if returncode != 0:
            raise RuntimeError(f"subprocess exit {returncode}")
        entry["source_sha256_after"] = verify_hashes(
            EXPECTED_SOURCE_SHA256, "source after run"
        )
        entry["reference_sha256_after"] = verify_hashes(
            EXPECTED_REFERENCE_SHA256, "reference after run"
        )
        entry["validation"] = validate_spec(spec)
    except BaseException as exc:
        entry.update({
            "status": "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
            "finished_at": now_iso(), "runtime_sec": time.time() - started,
            "error": f"{type(exc).__name__}: {exc}",
        })
        save_manifest(manifest)
        raise
    entry.update({"status": "completed", "finished_at": now_iso()})
    save_manifest(manifest)
    print(f"検証合格: {spec.run_id} ({entry['runtime_sec']:.1f}s)", flush=True)


def controls_ready(manifest: dict[str, Any]) -> bool:
    for spec in N4_CONTROLS:
        entry = manifest.get("runs", {}).get(spec.run_id, {})
        validation = entry.get("validation", {})
        if spec.mode in N4_CANONICAL:
            ok = validation.get("existing_T4000_exact_control", {}).get("byte_exact")
        else:
            ok = validation.get("new_T4000_control", {}).get("baseline_for_N4_T42000")
        if not ok:
            return False
    return True


def n4_ready(manifest: dict[str, Any]) -> bool:
    if not controls_ready(manifest):
        return False
    return all(
        manifest.get("runs", {}).get(spec.run_id, {})
        .get("validation", {}).get("first4000_vs_stage_control", {})
        .get("byte_exact", False)
        for spec in N4_LONG
    )


def n12_ready(manifest: dict[str, Any]) -> bool:
    for spec in N12_RUNS:
        validation = (
            manifest.get("runs", {}).get(spec.run_id, {}).get("validation", {})
        )
        if not validation.get("first42000_vs_existing", {}).get("byte_exact"):
            return False
        if spec.mode in N12_SHORT_REFERENCES and not (
            validation.get("first4000_vs_independent_existing", {}).get("byte_exact")
        ):
            return False
        if spec.mode == "fermion_family" and (
            validation.get("first4000_independent_evidence", {}).get("counted")
            is not False
        ):
            return False
    return True


def run_n4(manifest: dict[str, Any]) -> None:
    for index, spec in enumerate(N4_CONTROLS, 1):
        execute_one(spec, manifest, index, len(N4_CONTROLS))
    if not controls_ready(manifest):
        raise RuntimeError("N4 controls did not all pass; long runs are blocked")
    for index, spec in enumerate(N4_LONG, 1):
        execute_one(spec, manifest, index, len(N4_LONG))
    if not n4_ready(manifest):
        raise RuntimeError("N4 long prefix validation did not all pass")
    manifest["stage_N4_status"] = "complete"
    manifest["stage_N4_finished_at"] = now_iso()
    save_manifest(manifest)


def revalidate_n4_gate(manifest: dict[str, Any]) -> None:
    for spec in (*N4_CONTROLS, *N4_LONG):
        if output_state(spec) != "complete":
            raise RuntimeError(f"N12 blocked: {spec.run_id} is not complete")
        entry = manifest["runs"].setdefault(spec.run_id, {})
        validation = validate_spec(spec)
        assert_stored_artifact_hashes(spec, entry, validation)
        entry.update({
            "status": "validated_existing", "validated_at": now_iso(),
            "command": spec.command, "validation": validation,
        })
    if not n4_ready(manifest):
        raise RuntimeError("N12 blocked: N4 evidence gate failed")
    manifest["stage_N4_status"] = "complete_revalidated"
    save_manifest(manifest)


def run_n12(manifest: dict[str, Any]) -> None:
    revalidate_n4_gate(manifest)
    for index, spec in enumerate(N12_RUNS, 1):
        execute_one(spec, manifest, index, len(N12_RUNS))
    if not n12_ready(manifest):
        raise RuntimeError("N12 prefix validation did not all pass")
    manifest["stage_N12_status"] = "complete"
    manifest["stage_N12_finished_at"] = now_iso()
    save_manifest(manifest)


def smoke() -> dict[str, Any]:
    source = verify_hashes(EXPECTED_SOURCE_SHA256, "source smoke")
    references = verify_hashes(EXPECTED_REFERENCE_SHA256, "reference smoke")
    if len({spec.suffix for spec in ALL_SPECS}) != len(ALL_SPECS):
        raise AssertionError("output suffix collision")
    paths = [path for spec in ALL_SPECS for path in spec.all_outputs]
    if len(set(paths)) != len(paths):
        raise AssertionError("output path collision")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", spec.suffix)
           for spec in ALL_SPECS):
        raise AssertionError("invalid runner suffix")
    reference_json = validate_reference_jsons()
    reference_npz = validate_reference_npz_shapes()
    existing_prefix = validate_existing_reference_prefixes()
    capacity = capacity_report(list(ALL_SPECS))
    enforce_capacity(capacity)
    return {
        "smoke": "PASS", "experiments_executed": 0,
        "planned_runs": len(ALL_SPECS),
        "source_files_checked": len(source),
        "reference_files_checked": len(references),
        "reference_json": reference_json,
        "reference_npz": reference_npz,
        "existing_N12_T4000_vs_T42000": existing_prefix,
        "capacity": capacity,
    }


def selected_specs(args: argparse.Namespace) -> tuple[RunSpec, ...]:
    if args.run_n4:
        return (*N4_CONTROLS, *N4_LONG)
    if args.run_n12:
        return N12_RUNS
    if args.run_all:
        return ALL_SPECS
    return ALL_SPECS


def print_schedule() -> None:
    report = capacity_report(list(ALL_SPECS))
    estimates = {x["run_id"]: x for x in report["spec_estimates"]}
    print(f"planned_runs={len(ALL_SPECS)} sequential=true delta={DELTA}")
    for index, spec in enumerate(ALL_SPECS, 1):
        est = estimates[spec.run_id]
        print(
            f"{index:2d} {spec.run_id:34s} N={spec.n} T={spec.steps} "
            f"npz~{est['estimated_compressed_npz_bytes']/1024**2:.1f}MiB "
            f"RAM~{est['estimated_peak_ram_bytes']/1024**3:.2f}GiB"
        )
        print("   ", shlex.join(spec.command))
        print("   ->", spec.npz_path.name)
    available_text = (
        "unknown" if report["available_memory_bytes"] is None
        else f"{report['available_memory_bytes']/1024**3:.2f}GiB"
    )
    print(
        f"estimated_remaining_output={report['estimated_remaining_output_bytes']/1024**3:.2f}GiB "
        f"disk_required_gate={report['disk_required_bytes']/1024**3:.2f}GiB "
        f"disk_free={report['disk_free_bytes']/1024**3:.2f}GiB "
        f"peak_RAM={report['estimated_peak_ram_bytes']/1024**3:.2f}GiB "
        f"available_RAM={available_text} "
        f"runtime~{report['estimated_runtime_sec']/3600:.2f}h "
        f"gate_pass={report['gate_pass']}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--list", action="store_true", help="予定と見積りだけ表示")
    group.add_argument("--smoke", action="store_true", help="read-only短時間検査")
    group.add_argument("--run-n4", action="store_true", help="Stage N4を逐次実行/再開")
    group.add_argument("--run-n12", action="store_true", help="N4再検証後Stage N12を実行/再開")
    group.add_argument("--run-all", action="store_true", help="N4→N12を逐次実行/再開")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.list:
        print_schedule()
        return 0
    if args.smoke:
        print(json.dumps(smoke(), ensure_ascii=False, indent=2, default=float))
        return 0

    with campaign_lock():
        MPL_CACHE.mkdir(parents=True, exist_ok=True)
        source_before = verify_hashes(EXPECTED_SOURCE_SHA256, "source before campaign")
        reference_before = verify_hashes(
            EXPECTED_REFERENCE_SHA256, "reference before campaign"
        )
        validate_reference_jsons()
        validate_reference_npz_shapes()
        manifest = load_manifest()
        current_orchestrator_sha256 = sha256(Path(__file__).resolve())
        recorded_orchestrator_sha256 = manifest.get("orchestrator_sha256")
        if (
            recorded_orchestrator_sha256 is not None
            and recorded_orchestrator_sha256 != current_orchestrator_sha256
        ):
            raise RuntimeError(
                "orchestrator changed since manifest creation: "
                f"{recorded_orchestrator_sha256} != {current_orchestrator_sha256}"
            )
        report = capacity_report(list(selected_specs(args)))
        enforce_capacity(report)
        manifest.update({
            "orchestrator_sha256": current_orchestrator_sha256,
            "argv": sys.argv,
            "mplconfigdir": str(MPL_CACHE),
            "capacity_gate_before": report,
            "source_sha256_before": source_before,
            "reference_sha256_before": reference_before,
            "campaign_started_or_resumed_at": now_iso(),
        })
        save_manifest(manifest)

        if args.run_n4 or args.run_all:
            run_n4(manifest)
        if args.run_n12 or args.run_all:
            run_n12(manifest)

        manifest["source_sha256_after"] = verify_hashes(
            EXPECTED_SOURCE_SHA256, "source after campaign"
        )
        manifest["reference_sha256_after"] = verify_hashes(
            EXPECTED_REFERENCE_SHA256, "reference after campaign"
        )
        manifest["status"] = (
            "complete" if n4_ready(manifest) and n12_ready(manifest)
            else "N4_complete" if n4_ready(manifest)
            else "N12_complete" if n12_ready(manifest)
            else "incomplete"
        )
        manifest["finished_at"] = now_iso()
        save_manifest(manifest)
    print(f"Stage 4 complete: manifest={MANIFEST.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
