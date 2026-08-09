#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""不足しているシード型×δ長時間走行を、再現対照の合格後に逐次実行する。

科学条件は既存 mixed 長時間掃引と同じ N=12, T=42000, Nn=16,
Neta=8, window=[21000,42000]。既存成果物は上書きしない。

使い方:
  python3 run_missing_seed_sweeps_T42000_v1.py --list
  python3 run_missing_seed_sweeps_T42000_v1.py
  python3 run_missing_seed_sweeps_T42000_v1.py --controls-only

既に完成し検証に通る成果物は再開時に読み直して skip する。一部だけ存在する
成果物は上書きせず停止する。
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_nsweep_three_series_v2.py"
BASE = HERE / "run_tb_nsweep_1to20_v1.py"
UF = HERE.parent / "統一万能関数_v1"
MANIFEST = HERE / "manifest_不足シード型一括掃引_T42000_v1.json"
PREREGISTRATION = HERE / "事前登録_シード型別δ一括掃引_T42000_v1.md"
MPL_CACHE = Path(tempfile.gettempdir()) / "seed_type_delta_sweep_mplconfig_v1"

N = 12
T = 42000
CORE_DELTAS = (
    1e-15,
    1e-8,
    1e-4,
    1e-3,
    1e-2,
    0.03162277660168379,
    0.04357,
    0.1,
)

EXPECTED_SOURCE_SHA256 = {
    RUNNER: "c4b3b33a82657232734223f089cff79e103ce7a18fb5a064fc71a058d1d006fe",
    BASE: "1615b8769cf765a444a8e057718ff4a760f040376335f434f9903b8cf4f9dba4",
    UF / "unified_interaction_v1.py": (
        "b4db659b9bf958246d192c852a48fd97f5793f9bf2f891d5c75ce22ae5a5ed63"
    ),
    UF / "unified_dimension_v1.py": (
        "cfc9cf5d835e39ec6959eeeadff102e5dd23b93deed9dfeff52778e3ee48830a"
    ),
    UF / "unified_readout_v3.py": (
        "4dec07c9811c3c06b3c34941ade85593ac3feecbeb82e6902cecad1263958953"
    ),
    UF / "selection_v1.py": (
        "90d546e61e20ce58b4e9a962bdcd761ef96d46a9061a60decca2586ad324eed2"
    ),
}

REFERENCE_SHA256 = {
    HERE / "nsweep_mixed_T42000_N12_v2.npz": (
        "0de730a69d303fa7af0a88f09d768bea6f07159348a5c5b4b097cf60c1ad01ff"
    ),
    HERE / "nsweep_neutral_T42000_N12_v2.npz": (
        "6c78a39c37e50243ef78383d7ec8d94e7f741c207d5700d29bc55508f3be33c0"
    ),
    HERE / "nsweep_vacuum_T42000_N12_v2.npz": (
        "568ae5f6c576f7b78388a5498f1b5927a4f555fb2f841cefa783b6970d8b1c12"
    ),
}


@dataclass(frozen=True)
class RunSpec:
    mode: str
    delta: float
    suffix: str = ""
    control: bool = False

    @property
    def delta_tag(self) -> str:
        return f"{self.delta:g}"

    @property
    def tag(self) -> str:
        tag = f"_T{T}_d{self.delta_tag}"
        if self.suffix:
            tag += f"_rep-{self.suffix}"
        return tag

    @property
    def key(self) -> str:
        role = "control" if self.control else "science"
        return f"{role}:{self.mode}:d{self.delta_tag}:{self.suffix or '-'}"

    @property
    def command(self) -> list[str]:
        cmd = [sys.executable, "-u", str(RUNNER), self.mode,
               str(N), str(N), str(T), repr(self.delta)]
        if self.suffix:
            cmd.append(self.suffix)
        return cmd

    @property
    def json_path(self) -> Path:
        return HERE / f"result_nsweep_{self.mode}{self.tag}_v2.json"

    @property
    def npz_path(self) -> Path:
        return HERE / f"nsweep_{self.mode}{self.tag}_N{N}_v2.npz"

    @property
    def figure_paths(self) -> list[Path]:
        stem = f"fig_{self.mode}{self.tag}"
        return [
            HERE / f"{stem}_4panel_N{N}_v2.png",
            HERE / f"{stem}_mix_N{N}_v2.png",
            HERE / f"{stem}_ledger_N{N}_v2.png",
            HERE / f"{stem}_summary_v2.png",
            HERE / f"{stem}_birth_matrix_v2.png",
        ]

    @property
    def artifacts(self) -> list[Path]:
        return [self.json_path, self.npz_path, *self.figure_paths]

    @property
    def log_path(self) -> Path:
        return HERE / f"一括掃引ログ_{self.mode}{self.tag}_v1.txt"


CONTROLS = (
    RunSpec("mixed", 1e-2, "controlrep1", True),
    RunSpec("neutral", 1e-2, "controlrep1", True),
)


def science_specs() -> list[RunSpec]:
    specs: list[RunSpec] = []
    for mode in ("neutral", "electron", "fermion_family", "boson_family"):
        for delta in CORE_DELTAS:
            if mode == "neutral" and delta == 1e-2:
                continue
            specs.append(RunSpec(mode, delta))
    return specs


SCIENCE = tuple(science_specs())


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_hashes(expected: dict[Path, str], label: str) -> dict[str, str]:
    actual: dict[str, str] = {}
    bad: list[str] = []
    for path, wanted in expected.items():
        if not path.is_file():
            bad.append(f"missing: {path}")
            continue
        got = sha256(path)
        actual[str(path.relative_to(HERE) if path.is_relative_to(HERE) else path)] = got
        if got != wanted:
            bad.append(f"{path.name}: expected {wanted}, got {got}")
    if bad:
        raise RuntimeError(f"{label} SHA-256 mismatch: " + " | ".join(bad))
    return actual


def array_bytes(a: np.ndarray) -> bytes:
    return np.ascontiguousarray(a).tobytes(order="C")


def compare_npz_subset(reference: Path, candidate: Path, *,
                       prefix: bool = False,
                       keys: list[str] | None = None) -> dict:
    """reference の指定全配列を candidate の同名配列と byte exact 比較する。"""
    checked = 0
    with np.load(reference, allow_pickle=False) as ref, \
            np.load(candidate, allow_pickle=False) as cand:
        names = list(ref.files) if keys is None else list(keys)
        missing = [k for k in names if k not in cand.files]
        if missing:
            raise AssertionError(f"{candidate.name}: missing NPZ keys {missing}")
        for key in names:
            x = np.asarray(ref[key])
            y = np.asarray(cand[key])
            if prefix and x.shape != y.shape:
                if (x.ndim < 1 or y.ndim != x.ndim or
                        y.shape[1:] != x.shape[1:] or y.shape[0] < x.shape[0]):
                    raise AssertionError(
                        f"{key}: prefix-incompatible shape {x.shape} vs {y.shape}"
                    )
                y = y[:x.shape[0]]
            if x.dtype != y.dtype or x.shape != y.shape:
                raise AssertionError(
                    f"{key}: dtype/shape {x.dtype}/{x.shape} != {y.dtype}/{y.shape}"
                )
            if array_bytes(x) != array_bytes(y):
                raise AssertionError(f"{key}: array bytes differ")
            checked += 1
    return {"reference": reference.name, "candidate": candidate.name,
            "prefix": prefix, "arrays_checked": checked,
            "max_abs_difference": 0.0, "byte_exact": True}


def vacuum_keys(reference: Path) -> list[str]:
    # mode-dependent な primary巻き/相棒/target 計器は、同じ真空状態でも
    # 読むセルが違うため除外する。以下は状態そのものだけで決まる28配列。
    state_rec = {
        "r_mean", "r_med", "r_min", "r_max", "r_raw",
        "dist_alpha", "absdist_alpha",
        "r_nopump", "dist_alpha_nopump", "even_power_nopump",
        "odd_power", "odd_amp_max", "even_power", "total_power",
        "bands", "ledger", "ledger_t",
    }
    with np.load(reference, allow_pickle=False) as z:
        return [
            k for k in z.files
            if k.startswith("v_") or
            (k.startswith("rec_v_") and k.removeprefix("rec_v_") in state_rec)
        ]


def validate_json(spec: RunSpec) -> dict:
    data = json.loads(spec.json_path.read_text())
    env = data["env"]
    expected = {
        "mode": spec.mode,
        "Nn": 16,
        "Neta": 8,
        "T": T,
        "delta": spec.delta,
        "seed": 2,
        "cell": [2, 0],
        "order": 6,
        "window": [21000, 42000],
    }
    for key, wanted in expected.items():
        if env.get(key) != wanted:
            raise AssertionError(
                f"{spec.json_path.name}: env.{key}={env.get(key)!r}, expected {wanted!r}"
            )
    if env.get("output_suffix") != (spec.suffix or None):
        raise AssertionError(f"{spec.json_path.name}: output_suffix mismatch")
    rec = data.get("N", {}).get(str(N))
    if not rec or not rec.get("built"):
        raise AssertionError(f"{spec.json_path.name}: N=12 is not built")
    if data.get("failed_N") != []:
        raise AssertionError(f"{spec.json_path.name}: failed_N={data.get('failed_N')}")
    return {"env_fixed": True, "built_N12": True, "failed_N": []}


def validate_npz_shape(spec: RunSpec) -> dict:
    with np.load(spec.npz_path, allow_pickle=False) as z:
        for key in ("m_f2", "v_f2", "rec_m_r_mean", "rec_v_r_mean"):
            if key not in z.files or z[key].shape[0] != T:
                raise AssertionError(
                    f"{spec.npz_path.name}: {key} length is not {T}"
                )
        ledger_t = np.asarray(z["rec_m_ledger_t"])
        expected_t = np.array([1, *range(50, T + 1, 50)], dtype=float)
        if ledger_t.dtype != expected_t.dtype or array_bytes(ledger_t) != array_bytes(expected_t):
            raise AssertionError(f"{spec.npz_path.name}: ledger_t cadence mismatch")
        return {"npz_arrays": len(z.files), "time_steps": T,
                "ledger_snapshots": int(len(ledger_t))}


def validate_artifacts(spec: RunSpec) -> dict:
    missing = [p.name for p in spec.artifacts if not p.is_file()]
    if missing:
        raise AssertionError(f"{spec.key}: missing artifacts {missing}")
    out = {"json": validate_json(spec), "npz": validate_npz_shape(spec)}
    vac = HERE / "nsweep_vacuum_T42000_N12_v2.npz"
    out["internal_vacuum"] = compare_npz_subset(
        vac, spec.npz_path, keys=vacuum_keys(vac)
    )
    out["artifact_sha256"] = {
        p.name: sha256(p) for p in spec.artifacts
    }
    return out


def validate_control(spec: RunSpec) -> dict:
    if spec.mode == "mixed":
        long_ref = HERE / "nsweep_mixed_T42000_N12_v2.npz"
        short_ref = HERE / "nsweep_mixed_N12_v2.npz"
    elif spec.mode == "neutral":
        long_ref = HERE / "nsweep_neutral_T42000_N12_v2.npz"
        short_ref = HERE / "nsweep_neutral_N12_v2.npz"
    else:
        raise AssertionError(f"unknown control mode {spec.mode}")
    return {
        "all_existing_T42000_arrays": compare_npz_subset(long_ref, spec.npz_path),
        "all_existing_T4000_arrays_vs_prefix": compare_npz_subset(
            short_ref, spec.npz_path, prefix=True
        ),
    }


def load_manifest() -> dict:
    if MANIFEST.is_file():
        return json.loads(MANIFEST.read_text())
    return {
        "schema": "missing-seed-type-sweeps-T42000-v1",
        "created_at": now_iso(),
        "conditions": {
            "N": N, "T": T, "Nn": 16, "Neta": 8,
            "window": [21000, 42000], "parent_seed": 2,
            "cell": [2, 0], "order": 6,
            "core_deltas": list(CORE_DELTAS),
            "science_runs": len(SCIENCE), "control_runs": len(CONTROLS),
            "parallel": False,
        },
        "runs": {},
    }


def save_manifest(manifest: dict) -> None:
    manifest["updated_at"] = now_iso()
    MANIFEST.write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")


def existing_state(spec: RunSpec) -> str:
    count = sum(p.exists() for p in spec.artifacts)
    if count == 0:
        return "absent"
    if count == len(spec.artifacts):
        return "complete"
    present = [p.name for p in spec.artifacts if p.exists()]
    raise FileExistsError(
        f"{spec.key}: partial artifacts exist; refusing overwrite: {present}"
    )


def execute_one(spec: RunSpec, manifest: dict, index: int, total: int) -> None:
    print(f"\n===== [{index}/{total}] {spec.key} =====", flush=True)
    state = existing_state(spec)
    if state == "complete":
        validation = validate_artifacts(spec)
        if spec.control:
            validation["control"] = validate_control(spec)
        manifest["runs"][spec.key] = {
            "status": "validated_existing", "validated_at": now_iso(),
            "command": spec.command, "validation": validation,
        }
        save_manifest(manifest)
        print("既存の完成成果物を再検証し skip", flush=True)
        return
    if spec.log_path.exists():
        raise FileExistsError(
            f"{spec.key}: log exists without complete artifacts: {spec.log_path.name}"
        )

    verify_hashes(EXPECTED_SOURCE_SHA256, "source before run")
    started = time.time()
    entry = {
        "status": "running", "started_at": now_iso(),
        "command": spec.command, "log": spec.log_path.name,
    }
    manifest["runs"][spec.key] = entry
    save_manifest(manifest)

    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(MPL_CACHE)
    with spec.log_path.open("x", encoding="utf-8") as log:
        proc = subprocess.Popen(
            spec.command, cwd=HERE, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        rc = proc.wait()
    entry["returncode"] = rc
    entry["runtime_sec"] = time.time() - started
    if rc != 0:
        entry["status"] = "failed"
        entry["finished_at"] = now_iso()
        save_manifest(manifest)
        raise RuntimeError(f"{spec.key}: subprocess exit {rc}")

    validation = validate_artifacts(spec)
    if spec.control:
        validation["control"] = validate_control(spec)
    entry.update({
        "status": "completed", "finished_at": now_iso(),
        "validation": validation,
    })
    save_manifest(manifest)
    print(f"検証合格: {spec.key} ({entry['runtime_sec']:.1f}s)", flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--list", action="store_true", help="実行予定だけ表示する")
    p.add_argument("--controls-only", action="store_true",
                   help="完全再現対照2本の検証までで止める")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    specs = list(CONTROLS) + ([] if args.controls_only else list(SCIENCE))
    if args.list:
        print(f"controls={len(CONTROLS)} science={len(SCIENCE)} total={len(specs)}")
        for i, spec in enumerate(specs, 1):
            print(i, spec.key, "->", spec.npz_path.name)
        return

    MPL_CACHE.mkdir(parents=True, exist_ok=True)
    source_hashes = verify_hashes(EXPECTED_SOURCE_SHA256, "source")
    reference_hashes_before = verify_hashes(REFERENCE_SHA256, "reference before")
    manifest = load_manifest()
    manifest["orchestrator_sha256"] = sha256(Path(__file__).resolve())
    manifest["preregistration_sha256"] = sha256(PREREGISTRATION)
    manifest["source_sha256"] = source_hashes
    manifest["reference_sha256_before"] = reference_hashes_before
    manifest["mplconfigdir"] = str(MPL_CACHE)
    save_manifest(manifest)

    for index, spec in enumerate(specs, 1):
        execute_one(spec, manifest, index, len(specs))
        # 対照のどちらか一方でも不一致なら execute_one 内でここへ戻らない。

    manifest["source_sha256_after"] = verify_hashes(
        EXPECTED_SOURCE_SHA256, "source after batch"
    )
    manifest["reference_sha256_after"] = verify_hashes(
        REFERENCE_SHA256, "reference after batch"
    )
    manifest["status"] = "controls_complete" if args.controls_only else "complete"
    manifest["finished_at"] = now_iso()
    save_manifest(manifest)
    print(f"\n一括走行完了: {len(specs)} executions / manifest={MANIFEST.name}")


if __name__ == "__main__":
    main()
