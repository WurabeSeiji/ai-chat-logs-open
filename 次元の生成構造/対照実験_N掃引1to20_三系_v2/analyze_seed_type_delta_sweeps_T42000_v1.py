#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""シード型別δ一括掃引 (N=12, T=42000) の生 NPZ 後処理 v1。

このスクリプトは力学を実行しない。既存の JSON/NPZ を読み、

* core 8 δ × 5 mode の取得済み/部分/欠測を区別する。
* f2(1), tau_space, r_nopump/r_mean/r_raw, 通過・滞在、保存量、
  帯/帳簿の瞬間 CV と時間平均 CV を計算する。
* JSON + CSV + 同じ時間軸の比較図を、明示的に起動した時のみ生成する。

使い方:
  python3 analyze_seed_type_delta_sweeps_T42000_v1.py --discover-only
  python3 analyze_seed_type_delta_sweeps_T42000_v1.py --require-complete

注意:
  R=0.7, R_alpha, R_MZ, R_(620,117) は独立な比較中心である。
  数値的な近接や通過は、物理的同一性の証明ではない。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
N_TARGET = 12
T_TARGET = 42000
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
MODES = (
    "mixed",
    "neutral",
    "electron",
    "fermion_family",
    "boson_family",
)

# これらを一つの値として扱わない。由来と主張状態は REFERENCES に固定する。
R_DECIMAL_07 = 0.7
R_ALPHA = math.cos(23.0 * math.pi / 124.0) ** 2
R_MZ = 0.687822933884774
R_620_117 = math.cos(117.0 * math.pi / 620.0) ** 2
REFERENCES = {
    "R_decimal_0p7": {
        "value": R_DECIMAL_07,
        "definition": "decimal comparison level 0.7",
        "claim_status": "観測上の比較値。有限位数根や物理定数との同一性を定義しない。",
    },
    "R_alpha_124_23": {
        "value": R_ALPHA,
        "definition": "cos^2(23*pi/124)",
        "claim_status": (
            "有限位数根としては式から導出。alpha/charge との物理的対応は"
            "この時間発展実験が証明する定理ではない。"
        ),
    },
    "R_MZ_physical_correspondence": {
        "value": R_MZ,
        "definition": "supplied physical-correspondence comparator",
        "claim_status": (
            "外部から与えられた物理対応の比較値。R_(620,117) の式値ではなく、"
            "この後処理では独立な中心として扱う。"
        ),
    },
    "R_620_117": {
        "value": R_620_117,
        "definition": "cos^2(117*pi/620)",
        "claim_status": (
            "有限位数根としては式から導出。R_MZ との数値的近接は記録するが、"
            "物理的同一性として合并しない。"
        ),
    },
}
DEFAULT_BAND_HALFWIDTHS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)

REQUIRED_NPZ_KEYS = (
    "m_f2",
    "rec_m_r_nopump",
    "rec_m_r_mean",
    "rec_m_r_raw",
    "rec_m_total_power",
    "rec_m_bands",
    "rec_m_ledger",
    "rec_m_ledger_t",
)
EXPECTED_ENV = {
    "Nn": 16,
    "Neta": 8,
    "T": T_TARGET,
    "seed": 2,
    "cell": [2, 0],
    "order": 6,
    "window": [T_TARGET // 2, T_TARGET],
}

LEGACY_COUNTS = {
    "mixed": (5, 3),
    "neutral": (1, 0),
    "electron": (1, 0),
    "fermion_family": (5, 0),
    "boson_family": (0, 3),
}


def _finite_float(value: Any) -> float | None:
    """JSON 標準に入る有限 float だけを返す。"""
    if value is None:
        return None
    try:
        ans = float(value)
    except (TypeError, ValueError):
        return None
    return ans if math.isfinite(ans) else None


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return _finite_float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def _delta_label(delta: float) -> str:
    return f"{delta:.15g}"


def _match_core_delta(value: Any) -> float | None:
    observed = _finite_float(value)
    if observed is None or observed <= 0.0:
        return None
    for target in CORE_DELTAS:
        if math.isclose(observed, target, rel_tol=5e-12, abs_tol=1e-18):
            return target
    return None


def _n_record(payload: dict[str, Any]) -> dict[str, Any] | None:
    records = payload.get("N", {})
    rec = records.get(str(N_TARGET), records.get(N_TARGET))
    return rec if isinstance(rec, dict) else None


def _npz_path_from_json(path: Path) -> Path:
    stem = path.name
    if not (stem.startswith("result_nsweep_") and stem.endswith("_v2.json")):
        raise ValueError(f"非対応 JSON 名: {path.name}")
    core = stem[len("result_") : -len("_v2.json")]
    return path.with_name(f"{core}_N{N_TARGET}_v2.npz")


def _inspect_candidate(path: Path, expected_mode: str) -> dict[str, Any]:
    candidate: dict[str, Any] = {
        "json": path.name,
        "npz": None,
        "mode": expected_mode,
        "delta": None,
        "status": "partial",
        "issues": [],
        "output_suffix": None,
        "metadata_complete": False,
        "raw_complete": False,
        "conditions_valid": False,
    }
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # 同時書き込み中も partial として捕捉
        candidate["issues"].append(f"JSON read error: {type(exc).__name__}: {exc}")
        return candidate

    candidate["_json_payload"] = payload
    env = payload.get("env", {})
    candidate["output_suffix"] = env.get("output_suffix")
    candidate["delta"] = _match_core_delta(env.get("delta"))
    conditions_valid = True
    if env.get("mode") != expected_mode:
        conditions_valid = False
        candidate["issues"].append(
            f"env.mode={env.get('mode')!r} != filename mode={expected_mode!r}"
        )
    if env.get("T") != T_TARGET:
        conditions_valid = False
        candidate["issues"].append(f"env.T={env.get('T')!r} != {T_TARGET}")
    if candidate["delta"] is None:
        conditions_valid = False
        candidate["issues"].append(f"core8 外 delta={env.get('delta')!r}")
    for key, expected in EXPECTED_ENV.items():
        observed = env.get(key)
        if observed != expected:
            conditions_valid = False
            candidate["issues"].append(
                f"env.{key}={observed!r} != fixed {expected!r}"
            )
    functions = env.get("functions", [])
    required_functions = (
        "unified_interaction_v1",
        "unified_dimension_v1",
        "unified_readout_v3",
        "selection_v1",
    )
    if (
        not isinstance(functions, list)
        or len(functions) < len(required_functions)
        or any(required not in str(observed)
               for required, observed in zip(required_functions, functions))
    ):
        conditions_valid = False
        candidate["issues"].append(
            "env.functions does not declare fixed F/D/G/S in order"
        )

    rec = _n_record(payload)
    if rec is None:
        conditions_valid = False
        candidate["issues"].append(f"N={N_TARGET} record missing")
    elif rec.get("built") is not True:
        conditions_valid = False
        candidate["issues"].append(f"N={N_TARGET} built is not true")
    candidate["conditions_valid"] = conditions_valid

    try:
        npz_path = _npz_path_from_json(path)
    except ValueError as exc:
        candidate["issues"].append(str(exc))
        return candidate
    candidate["npz"] = npz_path.name
    candidate["_npz_path"] = npz_path
    if not npz_path.exists():
        candidate["issues"].append("NPZ missing")
        return candidate

    try:
        with np.load(npz_path, allow_pickle=False) as data:
            missing_keys = [key for key in REQUIRED_NPZ_KEYS if key not in data.files]
            if missing_keys:
                candidate["issues"].append("NPZ keys missing: " + ",".join(missing_keys))
            else:
                bad_lengths = {
                    key: int(data[key].shape[0])
                    for key in (
                        "m_f2",
                        "rec_m_r_nopump",
                        "rec_m_r_mean",
                        "rec_m_r_raw",
                        "rec_m_total_power",
                        "rec_m_bands",
                    )
                    if data[key].ndim == 0 or data[key].shape[0] != T_TARGET
                }
                if bad_lengths:
                    candidate["issues"].append(f"NPZ length mismatch: {bad_lengths}")
                else:
                    bands_shape = tuple(data["rec_m_bands"].shape)
                    ledger_shape = tuple(data["rec_m_ledger"].shape)
                    ledger_t = np.asarray(data["rec_m_ledger_t"], dtype=float)
                    shape_issues = []
                    if bands_shape != (T_TARGET, 16):
                        shape_issues.append(f"bands shape={bands_shape}")
                    if len(ledger_shape) != 3 or ledger_shape[1:] != (16, 8):
                        shape_issues.append(f"ledger shape={ledger_shape}")
                    if ledger_t.ndim != 1 or ledger_t.size != ledger_shape[0]:
                        shape_issues.append(
                            f"ledger_t shape={ledger_t.shape} vs ledger rows={ledger_shape[0]}"
                        )
                    elif (
                        ledger_t.size == 0
                        or ledger_t[0] != 1.0
                        or ledger_t[-1] != float(T_TARGET)
                    ):
                        shape_issues.append(
                            "ledger_t endpoints are not fixed 1..42000"
                        )
                    if shape_issues:
                        candidate["issues"].append(
                            "NPZ shape mismatch: " + "; ".join(shape_issues)
                        )
                    else:
                        candidate["raw_complete"] = True
    except Exception as exc:
        candidate["issues"].append(f"NPZ read error: {type(exc).__name__}: {exc}")

    candidate["metadata_complete"] = (
        "runtime_sec" in payload
        and "judgments" in payload
        and rec is not None
        and rec.get("built") is True
    )
    if (
        candidate["conditions_valid"]
        and candidate["raw_complete"]
        and candidate["metadata_complete"]
    ):
        candidate["status"] = "complete"
    elif candidate["conditions_valid"] and candidate["raw_complete"]:
        candidate["status"] = "raw_complete_metadata_partial"
    return candidate


def _public_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in candidate.items() if not k.startswith("_")}


def discover(data_dir: Path) -> dict[str, Any]:
    """core8 行列を探索し、物理条件を変えず欠測だけを列挙する。"""
    matrix: dict[str, dict[float, dict[str, Any]]] = {
        mode: {} for mode in MODES
    }
    ignored: list[dict[str, Any]] = []
    all_candidates: dict[tuple[str, float], list[dict[str, Any]]] = {}

    for mode in MODES:
        for path in sorted(data_dir.glob(f"result_nsweep_{mode}*_v2.json")):
            candidate = _inspect_candidate(path, mode)
            delta = candidate.get("delta")
            if delta is None or candidate.get("_json_payload", {}).get("env", {}).get("T") != T_TARGET:
                ignored.append(_public_candidate(candidate))
                continue
            all_candidates.setdefault((mode, delta), []).append(candidate)

    complete_count = raw_partial_count = partial_count = missing_count = 0
    selected_private: dict[tuple[str, float], dict[str, Any]] = {}
    for mode in MODES:
        for delta in CORE_DELTAS:
            candidates = all_candidates.get((mode, delta), [])
            # 解析可能な完了生データを優先し、同状態なら複製 suffix なしの
            # 正本を選ぶ。途中の正本が完了複製を隠さないようにする。
            ranked = sorted(
                candidates,
                key=lambda c: (
                    not c.get("raw_complete", False),
                    c.get("status") != "complete",
                    c.get("output_suffix") is not None,
                    c["json"],
                ),
            )
            selected = ranked[0] if ranked else None
            if selected is None:
                status = "missing"
                missing_count += 1
                entry = {
                    "mode": mode,
                    "delta": delta,
                    "status": status,
                    "selected": None,
                    "alternates": [],
                }
            else:
                status = selected["status"]
                if status == "complete":
                    complete_count += 1
                elif status == "raw_complete_metadata_partial":
                    raw_partial_count += 1
                else:
                    partial_count += 1
                if selected.get("raw_complete"):
                    selected_private[(mode, delta)] = selected
                entry = {
                    "mode": mode,
                    "delta": delta,
                    "status": status,
                    "selected": _public_candidate(selected),
                    "alternates": [_public_candidate(c) for c in ranked[1:]],
                }
            matrix[mode][delta] = entry

    rows = [
        matrix[mode][delta]
        for mode in MODES
        for delta in CORE_DELTAS
    ]
    return {
        "scope": {"N": N_TARGET, "T": T_TARGET, "modes": list(MODES),
                  "core_deltas": list(CORE_DELTAS)},
        "counts": {
            "expected": len(MODES) * len(CORE_DELTAS),
            "complete": complete_count,
            "raw_complete_metadata_partial": raw_partial_count,
            "partial": partial_count,
            "missing": missing_count,
        },
        "rows": rows,
        "ignored_json": ignored,
        "_selected_private": selected_private,
    }


def _series_summary(values: np.ndarray, late_start_index: int) -> dict[str, Any]:
    y = np.asarray(values, dtype=float).reshape(-1)
    finite = np.isfinite(y)
    finite_idx = np.flatnonzero(finite)
    result: dict[str, Any] = {
        "length": int(y.size),
        "finite_count": int(finite.sum()),
        "initial_t1": _finite_float(y[0]) if y.size else None,
        "final_tT": _finite_float(y[-1]) if y.size else None,
        "last_finite": None,
        "last_finite_step": None,
        "max": None,
        "argmax_step": None,
        "min": None,
        "argmin_step": None,
        "late_median": None,
        "late_mean": None,
        "late_slope_per_step": None,
    }
    if finite_idx.size == 0:
        return result
    last_idx = int(finite_idx[-1])
    result["last_finite"] = float(y[last_idx])
    result["last_finite_step"] = last_idx + 1
    valid_y = y[finite_idx]
    i_max = int(finite_idx[int(np.argmax(valid_y))])
    i_min = int(finite_idx[int(np.argmin(valid_y))])
    result.update({
        "max": float(y[i_max]),
        "argmax_step": i_max + 1,
        "min": float(y[i_min]),
        "argmin_step": i_min + 1,
    })

    late_idx = np.arange(y.size) >= late_start_index
    late_valid = late_idx & finite
    idx = np.flatnonzero(late_valid)
    if idx.size:
        vals = y[idx]
        result["late_median"] = float(np.median(vals))
        result["late_mean"] = float(np.mean(vals))
        if idx.size >= 2:
            x = idx.astype(float) + 1.0
            xc = x - x.mean()
            denom = float(np.dot(xc, xc))
            if denom > 0.0:
                result["late_slope_per_step"] = float(
                    np.dot(xc, vals - vals.mean()) / denom
                )
    return result


def _crossing_events(values: np.ndarray, center: float) -> dict[str, Any]:
    y = np.asarray(values, dtype=float).reshape(-1)
    if y.size < 2:
        return {"up": [], "down": [], "up_count": 0, "down_count": 0}
    a, b = y[:-1], y[1:]
    valid = np.isfinite(a) & np.isfinite(b)
    up_idx = np.flatnonzero(valid & (a < center) & (b >= center))
    down_idx = np.flatnonzero(valid & (a > center) & (b <= center))

    def events(indices: Iterable[int]) -> list[dict[str, Any]]:
        answer = []
        for raw_i in indices:
            i = int(raw_i)
            dy = float(b[i] - a[i])
            frac = float((center - a[i]) / dy) if dy != 0.0 else 0.0
            answer.append({
                "from_step": i + 1,
                "to_step": i + 2,
                "interpolated_step": float(i + 1 + frac),
                "local_slope_per_step": dy,
            })
        return answer

    ups, downs = events(up_idx), events(down_idx)
    return {
        "up": ups,
        "down": downs,
        "up_count": len(ups),
        "down_count": len(downs),
        "first_up_interpolated_step": ups[0]["interpolated_step"] if ups else None,
        "last_up_interpolated_step": ups[-1]["interpolated_step"] if ups else None,
        "first_down_interpolated_step": downs[0]["interpolated_step"] if downs else None,
        "last_down_interpolated_step": downs[-1]["interpolated_step"] if downs else None,
    }


def _true_intervals(mask: np.ndarray) -> list[tuple[int, int]]:
    flag = np.asarray(mask, dtype=bool).reshape(-1)
    if flag.size == 0:
        return []
    padded = np.concatenate(([False], flag, [False])).astype(np.int8)
    transitions = np.diff(padded)
    starts = np.flatnonzero(transitions == 1)
    ends_exclusive = np.flatnonzero(transitions == -1)
    return [(int(a), int(b - 1)) for a, b in zip(starts, ends_exclusive)]


def _dwell_summary(
    values: np.ndarray, center: float, halfwidth: float
) -> dict[str, Any]:
    y = np.asarray(values, dtype=float).reshape(-1)
    finite = np.isfinite(y)
    mask = finite & (np.abs(y - center) <= halfwidth)
    intervals = _true_intervals(mask)
    details = []
    for start, end in intervals:
        segment = y[start : end + 1]
        if segment.size >= 2:
            velocity = np.abs(np.diff(segment))
            mean_abs_velocity = float(np.mean(velocity))
            median_abs_velocity = float(np.median(velocity))
        else:
            mean_abs_velocity = median_abs_velocity = None
        entry_slope = (
            _finite_float(y[start] - y[start - 1])
            if start > 0 and np.isfinite(y[start - 1]) else None
        )
        exit_slope = (
            _finite_float(y[end + 1] - y[end])
            if end + 1 < y.size and np.isfinite(y[end + 1]) else None
        )
        details.append({
            "start_step": start + 1,
            "end_step": end + 1,
            "steps": end - start + 1,
            "entry_slope_per_step": entry_slope,
            "exit_slope_per_step": exit_slope,
            "inside_mean_abs_velocity": mean_abs_velocity,
            "inside_median_abs_velocity": median_abs_velocity,
        })
    total = int(mask.sum())
    finite_count = int(finite.sum())
    return {
        "halfwidth": halfwidth,
        "lower": center - halfwidth,
        "upper": center + halfwidth,
        "dwell_steps": total,
        "fraction_of_all_steps": total / y.size if y.size else None,
        "fraction_of_finite_steps": total / finite_count if finite_count else None,
        "interval_count": len(details),
        "longest_interval_steps": max((d["steps"] for d in details), default=0),
        "first_dwell_step": details[0]["start_step"] if details else None,
        "last_dwell_step": details[-1]["end_step"] if details else None,
        "intervals": details,
    }


def _reference_analysis(
    values: np.ndarray, halfwidths: tuple[float, ...]
) -> dict[str, Any]:
    y = np.asarray(values, dtype=float).reshape(-1)
    finite = np.isfinite(y)
    output: dict[str, Any] = {}
    for name, spec in REFERENCES.items():
        center = float(spec["value"])
        if finite.any():
            dist = np.full(y.shape, np.inf, dtype=float)
            dist[finite] = np.abs(y[finite] - center)
            nearest_index = int(np.argmin(dist))
            nearest = {
                "step": nearest_index + 1,
                "value": float(y[nearest_index]),
                "absolute_distance": float(dist[nearest_index]),
            }
        else:
            nearest = {"step": None, "value": None, "absolute_distance": None}
        output[name] = {
            "center": center,
            "crossings": _crossing_events(y, center),
            "nearest": nearest,
            "dwell_by_halfwidth": {
                f"{width:.1e}": _dwell_summary(y, center, width)
                for width in halfwidths
            },
        }
    return output


def _cv(vector: np.ndarray) -> float | None:
    values = np.asarray(vector, dtype=float).reshape(-1)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None
    mean = float(np.mean(values))
    if mean == 0.0:
        return 0.0 if np.all(values == 0.0) else None
    return float(np.std(values) / abs(mean))


def _matrix_cv_summary(
    matrix: np.ndarray,
    sample_steps: np.ndarray,
    late_first_step: int,
    feature_name: str,
) -> dict[str, Any]:
    values = np.asarray(matrix, dtype=float)
    if values.ndim < 2:
        return {"available": False, "reason": f"{feature_name} ndim < 2"}
    flat = values.reshape(values.shape[0], -1)
    steps = np.asarray(sample_steps, dtype=float).reshape(-1)
    if steps.size != flat.shape[0]:
        return {
            "available": False,
            "reason": f"sample steps {steps.size} != rows {flat.shape[0]}",
        }
    instant = np.array([
        np.nan if (item := _cv(row)) is None else item for row in flat
    ], dtype=float)
    late_mask = np.isfinite(steps) & (steps >= late_first_step)
    instant_late = instant[late_mask & np.isfinite(instant)]
    late_rows = flat[late_mask]
    if late_rows.size:
        with np.errstate(invalid="ignore"):
            time_average_features = np.nanmean(late_rows, axis=0)
        time_average_cv = _cv(time_average_features)
    else:
        time_average_cv = None
    return {
        "available": True,
        "feature_space": feature_name,
        "row_count": int(flat.shape[0]),
        "feature_count": int(flat.shape[1]),
        "late_first_step": late_first_step,
        "instantaneous_cv_definition": "population_std(features)/abs(mean(features)) per sample",
        "time_average_cv_definition": (
            "population_std(late-time mean of each feature)/"
            "abs(mean of late-time feature means)"
        ),
        "instantaneous_cv_initial": _finite_float(instant[0]) if instant.size else None,
        "instantaneous_cv_final": _finite_float(instant[-1]) if instant.size else None,
        "instantaneous_cv_late_median": (
            float(np.median(instant_late)) if instant_late.size else None
        ),
        "instantaneous_cv_late_mean": (
            float(np.mean(instant_late)) if instant_late.size else None
        ),
        "instantaneous_cv_late_min": (
            float(np.min(instant_late)) if instant_late.size else None
        ),
        "instantaneous_cv_late_max": (
            float(np.max(instant_late)) if instant_late.size else None
        ),
        "late_time_average_feature_cv": time_average_cv,
    }


def _power_drift(total_power: np.ndarray) -> dict[str, Any]:
    p = np.asarray(total_power, dtype=float).reshape(-1)
    finite = np.isfinite(p)
    if p.size == 0 or not finite.any() or not np.isfinite(p[0]):
        return {"available": False}
    initial = float(p[0])
    diff = p[finite] - initial
    final_diff = float(p[-1] - initial) if np.isfinite(p[-1]) else None
    scale = abs(initial)
    return {
        "available": True,
        "initial": initial,
        "final": _finite_float(p[-1]),
        "final_absolute_drift": final_diff,
        "max_absolute_drift": float(np.max(np.abs(diff))),
        "min_minus_initial": float(np.min(diff)),
        "max_minus_initial": float(np.max(diff)),
        "final_relative_drift": final_diff / scale if scale > 0 and final_diff is not None else None,
        "max_relative_absolute_drift": (
            float(np.max(np.abs(diff))) / scale if scale > 0 else None
        ),
    }


def _seed_strength(env: dict[str, Any], mode: str, delta: float) -> dict[str, Any]:
    strength = env.get("seed_strength")
    keys = ("nF", "nB", "PF", "PB", "Pseed", "Acoh")
    if isinstance(strength, dict) and all(key in strength for key in keys):
        return {
            "source": "JSON env.seed_strength",
            **{key: _finite_float(strength.get(key)) for key in keys},
            "n_cells": strength.get("n_cells"),
            "r_np_cell_ratio": _finite_float(strength.get("r_np_cell_ratio")),
            "r_np_power_ratio": _finite_float(strength.get("r_np_power_ratio")),
        }

    # 旧 mixed/neutral 正本は新メタデータ追加前。値は既存レシピから
    # 再構成し、JSON 実測値とは明示的に分ける。
    n_f, n_b = LEGACY_COUNTS[mode]
    p_f = n_f * delta * delta
    p_b = n_b * delta * delta
    n_total = n_f + n_b
    p_total = p_f + p_b
    return {
        "source": "legacy mode recipe inference (env.seed_strength absent)",
        "nF": n_f,
        "nB": n_b,
        "PF": p_f,
        "PB": p_b,
        "Pseed": p_total,
        "Acoh": n_total * delta,
        "n_cells": n_total,
        "r_np_cell_ratio": n_f / n_total if n_total else None,
        "r_np_power_ratio": p_f / p_total if p_total else None,
    }


def _optional_series(data: Any, key: str) -> np.ndarray | None:
    return np.asarray(data[key], dtype=float) if key in data.files else None


def _analyze_candidate(
    candidate: dict[str, Any],
    mode: str,
    delta: float,
    halfwidths: tuple[float, ...],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    payload = candidate["_json_payload"]
    env = payload.get("env", {})
    rec = _n_record(payload) or {}
    npz_path: Path = candidate["_npz_path"]
    late_start_index = T_TARGET // 2
    late_first_step = late_start_index + 1

    with np.load(npz_path, allow_pickle=False) as data:
        f2 = np.asarray(data["m_f2"], dtype=float)
        r_nopump = np.asarray(data["rec_m_r_nopump"], dtype=float)
        r_mean = np.asarray(data["rec_m_r_mean"], dtype=float)
        r_raw = np.asarray(data["rec_m_r_raw"], dtype=float)
        total_power = np.asarray(data["rec_m_total_power"], dtype=float)
        bands = np.asarray(data["rec_m_bands"], dtype=float)
        ledger = np.asarray(data["rec_m_ledger"], dtype=float)
        ledger_t = np.asarray(data["rec_m_ledger_t"], dtype=float)

        power_keys = (
            "rec_m_odd_power",
            "rec_m_even_power",
            "rec_m_even_power_nopump",
            "rec_m_total_power",
            "rec_m_target_power",
            "rec_m_nontarget_power",
            "rec_m_seed_power",
            "rec_m_primary_seed_power",
            "rec_m_pump_power",
        )
        power_series = {
            key.removeprefix("rec_m_"): _series_summary(data[key], late_start_index)
            for key in power_keys if key in data.files
        }

        balance = {
            "bands_16": _matrix_cv_summary(
                bands, np.arange(1, bands.shape[0] + 1), late_first_step, "16 bands"
            ),
            "ledger_128_cells": _matrix_cv_summary(
                ledger, ledger_t, late_first_step, "16x8=128 ledger cells"
            ),
        }
        targets = _optional_series(data, "rec_m_targets")
        if targets is not None and targets.ndim == 2 and targets.shape[1] > 0:
            balance["target_cells"] = _matrix_cv_summary(
                targets,
                np.arange(1, targets.shape[0] + 1),
                late_first_step,
                "declared target cells",
            )
        else:
            balance["target_cells"] = {
                "available": False,
                "reason": "no target-cell matrix",
            }

        run = {
            "mode": mode,
            "delta": delta,
            "status": candidate["status"],
            "source": {
                "json": candidate["json"],
                "npz": candidate["npz"],
                "output_suffix": candidate.get("output_suffix"),
                "npz_bytes": npz_path.stat().st_size,
            },
            "fixed_conditions_from_json": {
                key: env.get(key)
                for key in ("Nn", "Neta", "T", "seed", "cell", "order", "window")
            },
            "implementation_from_json": {
                "functions": env.get("functions"),
                "base_script_md5": env.get("base_script_md5"),
            },
            "seed_and_target_recipe_from_json": {
                "seed_cells": env.get("seed_cells"),
                "target_cells": env.get("target_cells"),
                "primary": env.get("primary"),
                "legacy_recipe": env.get("recipe"),
            },
            "seed_strength": _seed_strength(env, mode, delta),
            "f2_t1": _finite_float(f2[0]) if f2.size else None,
            "tau_space": rec.get("tau_space"),
            "readouts": {
                "r_nopump": _series_summary(r_nopump, late_start_index),
                "r_mean": _series_summary(r_mean, late_start_index),
                "r_raw": _series_summary(r_raw, late_start_index),
            },
            "r_nopump_reference_tests": _reference_analysis(r_nopump, halfwidths),
            "power_series": power_series,
            "total_power_drift": _power_drift(total_power),
            "distribution_balance": balance,
            "interpretation_limits": [
                "T=42000 内の未到達は、より長い時間での到達不可能を意味しない。",
                "参照値の通過、滞在、最近接は、参照値との物理的同一性を意味しない。",
                "同じ delta でもシードセル数が異なるため、総初期パワーは一致しない。",
            ],
        }
    plot = {"r_nopump": r_nopump, "r_mean": r_mean, "r_raw": r_raw}
    return run, plot


def _flatten_run_for_csv(run: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "mode": run["mode"],
        "delta": run["delta"],
        "status": run["status"],
    }
    source = run.get("source", {})
    row.update({"json": source.get("json"), "npz": source.get("npz")})
    if run["status"] not in {"complete", "raw_complete_metadata_partial"}:
        return row

    row["f2_t1"] = run.get("f2_t1")
    row["tau_space"] = run.get("tau_space")
    seed = run.get("seed_strength", {})
    row["seed_strength_source"] = seed.get("source")
    for key in ("nF", "nB", "PF", "PB", "Pseed", "Acoh"):
        row[key] = seed.get(key)

    for readout_name, stats in run.get("readouts", {}).items():
        for key in (
            "initial_t1", "max", "argmax_step", "final_tT", "last_finite",
            "last_finite_step", "late_median", "late_slope_per_step", "finite_count",
        ):
            row[f"{readout_name}_{key}"] = stats.get(key)

    drift = run.get("total_power_drift", {})
    for key in (
        "initial", "final", "final_absolute_drift", "max_absolute_drift",
        "final_relative_drift", "max_relative_absolute_drift",
    ):
        row[f"total_power_{key}"] = drift.get(key)

    for balance_name, summary in run.get("distribution_balance", {}).items():
        for key in (
            "instantaneous_cv_initial",
            "instantaneous_cv_final",
            "instantaneous_cv_late_median",
            "instantaneous_cv_late_mean",
            "late_time_average_feature_cv",
        ):
            row[f"{balance_name}_{key}"] = summary.get(key)

    for power_name, stats in run.get("power_series", {}).items():
        for key in ("initial_t1", "max", "argmax_step", "final_tT", "late_median",
                    "late_slope_per_step"):
            row[f"power_{power_name}_{key}"] = stats.get(key)

    for ref_name, result in run.get("r_nopump_reference_tests", {}).items():
        crossing = result["crossings"]
        nearest = result["nearest"]
        prefix = f"rnp_{ref_name}"
        for key in (
            "up_count", "down_count", "first_up_interpolated_step",
            "last_up_interpolated_step", "first_down_interpolated_step",
            "last_down_interpolated_step",
        ):
            row[f"{prefix}_{key}"] = crossing.get(key)
        row[f"{prefix}_nearest_step"] = nearest.get("step")
        row[f"{prefix}_nearest_value"] = nearest.get("value")
        row[f"{prefix}_nearest_absolute_distance"] = nearest.get("absolute_distance")
        for width_label, dwell in result["dwell_by_halfwidth"].items():
            wp = f"{prefix}_band_{width_label}"
            for key in (
                "dwell_steps", "fraction_of_all_steps", "interval_count",
                "longest_interval_steps", "first_dwell_step", "last_dwell_step",
            ):
                row[f"{wp}_{key}"] = dwell.get(key)
    return row


def _missing_csv_row(entry: dict[str, Any]) -> dict[str, Any]:
    selected = entry.get("selected") or {}
    return {
        "mode": entry["mode"],
        "delta": entry["delta"],
        "status": entry["status"],
        "json": selected.get("json"),
        "npz": selected.get("npz"),
        "issues": " | ".join(selected.get("issues", [])),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _figure_paths(data_dir: Path, prefix: str) -> list[Path]:
    return [
        data_dir / f"fig_{prefix}_r_nopump_by_delta.png",
        data_dir / f"fig_{prefix}_r_nopump_by_mode.png",
        data_dir / f"fig_{prefix}_r_mean_raw_by_mode.png",
    ]


def _prepare_output_paths(
    data_dir: Path, prefix: str, overwrite: bool
) -> tuple[Path, Path, list[Path]]:
    json_path = data_dir / f"{prefix}.json"
    csv_path = data_dir / f"{prefix}.csv"
    figures = _figure_paths(data_dir, prefix)
    occupied = [p for p in [json_path, csv_path, *figures] if p.exists()]
    if occupied and not overwrite:
        raise FileExistsError(
            "出力が既に存在する。--overwrite なしでは上書きしない: "
            + ", ".join(p.name for p in occupied)
        )
    return json_path, csv_path, figures


def _add_reference_lines(axis: Any, label: bool = False) -> None:
    styles = {
        "R_decimal_0p7": ("#111111", "--"),
        "R_alpha_124_23": ("#d62728", ":"),
        "R_MZ_physical_correspondence": ("#2ca02c", "-."),
        "R_620_117": ("#9467bd", (0, (1, 1))),
    }
    for name, spec in REFERENCES.items():
        color, linestyle = styles[name]
        axis.axhline(
            spec["value"], color=color, linestyle=linestyle, linewidth=0.9,
            alpha=0.8, label=(f"{name}={spec['value']:.9f}" if label else None),
        )


def _reference_handles(plt: Any) -> list[Any]:
    styles = {
        "R_decimal_0p7": ("#111111", "--"),
        "R_alpha_124_23": ("#d62728", ":"),
        "R_MZ_physical_correspondence": ("#2ca02c", "-."),
        "R_620_117": ("#9467bd", (0, (1, 1))),
    }
    return [
        plt.Line2D(
            [0], [0], color=styles[name][0], linestyle=styles[name][1],
            linewidth=1.0, label=f"{name}={spec['value']:.9f}",
        )
        for name, spec in REFERENCES.items()
    ]


def _save_figures(
    figure_paths: list[Path],
    plot_series: dict[tuple[str, float], dict[str, np.ndarray]],
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mode_colors = dict(zip(MODES, plt.get_cmap("tab10").colors[:len(MODES)]))
    delta_colors = dict(zip(CORE_DELTAS, plt.get_cmap("viridis")(
        np.linspace(0.05, 0.95, len(CORE_DELTAS))
    )))
    steps = np.arange(1, T_TARGET + 1)

    # 各 delta 内で mode を重ねる。8枚とも同じ横軸で欠測は描かない。
    fig, axes = plt.subplots(4, 2, figsize=(15, 13), sharex=True, sharey=True)
    for ax, delta in zip(axes.ravel(), CORE_DELTAS):
        for mode in MODES:
            item = plot_series.get((mode, delta))
            if item is None:
                continue
            y = item["r_nopump"]
            if np.isfinite(y).any():
                ax.plot(steps, y, color=mode_colors[mode], linewidth=0.8,
                        label=mode)
        _add_reference_lines(ax)
        ax.set_title(f"delta={_delta_label(delta)}")
        ax.set_xlim(1, T_TARGET)
        ax.grid(alpha=0.2)
    axes[-1, 0].set_xlabel("step (same x-axis: 1..42000)")
    axes[-1, 1].set_xlabel("step (same x-axis: 1..42000)")
    for ax in axes[:, 0]:
        ax.set_ylabel("r_nopump")
    handles = [plt.Line2D([0], [0], color=mode_colors[m], label=m) for m in MODES]
    handles.extend(_reference_handles(plt))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.975),
               ncol=5, frameon=False, fontsize=8)
    fig.suptitle("r_nopump: seed modes at each core delta (numerical comparison only)", y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(figure_paths[0], dpi=180)
    plt.close(fig)

    # 各 mode 内で delta を重ねる。
    fig, axes = plt.subplots(3, 2, figsize=(15, 11), sharex=True, sharey=True)
    flat_axes = axes.ravel()
    for ax, mode in zip(flat_axes, MODES):
        for delta in CORE_DELTAS:
            item = plot_series.get((mode, delta))
            if item is None:
                continue
            y = item["r_nopump"]
            if np.isfinite(y).any():
                ax.plot(steps, y, color=delta_colors[delta], linewidth=0.8,
                        label=_delta_label(delta))
        _add_reference_lines(ax)
        ax.set_title(mode)
        ax.set_xlim(1, T_TARGET)
        ax.grid(alpha=0.2)
    flat_axes[-1].axis("off")
    for ax in axes[-1, :]:
        if ax.axison:
            ax.set_xlabel("step (same x-axis: 1..42000)")
    for ax in axes[:, 0]:
        ax.set_ylabel("r_nopump")
    handles = [
        plt.Line2D([0], [0], color=delta_colors[d], label=_delta_label(d))
        for d in CORE_DELTAS
    ]
    handles.extend(_reference_handles(plt))
    flat_axes[-1].legend(handles=handles, title="delta / independent references",
                         loc="center", frameon=False, fontsize=8, ncol=2)
    fig.suptitle("r_nopump: core deltas within each seed mode (numerical comparison only)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(figure_paths[1], dpi=180)
    plt.close(fig)

    # r_mean と r_raw を合併せず、各 mode の左右に分ける。
    fig, axes = plt.subplots(len(MODES), 2, figsize=(16, 16), sharex=True, sharey=True)
    for row_index, mode in enumerate(MODES):
        for col_index, readout in enumerate(("r_mean", "r_raw")):
            ax = axes[row_index, col_index]
            for delta in CORE_DELTAS:
                item = plot_series.get((mode, delta))
                if item is None:
                    continue
                y = item[readout]
                if np.isfinite(y).any():
                    ax.plot(steps, y, color=delta_colors[delta], linewidth=0.7,
                            label=_delta_label(delta))
            _add_reference_lines(ax)
            ax.set_xlim(1, T_TARGET)
            ax.set_title(f"{mode}: {readout}")
            ax.grid(alpha=0.2)
            if col_index == 0:
                ax.set_ylabel("readout")
            if row_index == len(MODES) - 1:
                ax.set_xlabel("step (same x-axis: 1..42000)")
    handles = [
        plt.Line2D([0], [0], color=delta_colors[d], label=_delta_label(d))
        for d in CORE_DELTAS
    ]
    handles.extend(_reference_handles(plt))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.975),
               ncol=6, frameon=False,
               title="delta / independent references", fontsize=8)
    fig.suptitle("r_mean and r_raw kept as separate readouts", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(figure_paths[2], dpi=180)
    plt.close(fig)


def _reference_separations() -> dict[str, float]:
    return {
        "R_decimal_0p7_minus_R_alpha": R_DECIMAL_07 - R_ALPHA,
        "R_alpha_minus_R_MZ": R_ALPHA - R_MZ,
        "R_alpha_minus_R_620_117": R_ALPHA - R_620_117,
        "R_620_117_minus_R_MZ": R_620_117 - R_MZ,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_halfwidths(text: str) -> tuple[float, ...]:
    try:
        values = tuple(float(item.strip()) for item in text.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if not values or any((not math.isfinite(v) or v <= 0.0) for v in values):
        raise argparse.ArgumentTypeError("半幅は正の有限数をカンマ区切りで指定")
    return tuple(sorted(set(values), reverse=True))


def _print_discovery(report: dict[str, Any]) -> None:
    public = {k: v for k, v in report.items() if not k.startswith("_")}
    print(json.dumps(_json_safe(public), ensure_ascii=False, indent=2, allow_nan=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "N=12/T=42000 シード型別 core8 delta の生 NPZ を後処理する。"
            "力学の実行は行わない。"
        )
    )
    parser.add_argument("--data-dir", type=Path, default=HERE,
                        help="JSON/NPZ 格納ディレクトリ (既定: スクリプト所在地)")
    parser.add_argument("--discover-only", action="store_true",
                        help="取得/部分/欠測行列を stdout に出し、ファイルを生成しない")
    parser.add_argument("--require-complete", action="store_true",
                        help="40/40 の完了前は解析成果を生成せず終了する")
    parser.add_argument("--overwrite", action="store_true",
                        help="既存の後処理成果の上書きを許可する")
    parser.add_argument(
        "--output-prefix", default="seed_type_delta_sweeps_T42000_v1",
        help="JSON/CSV/図の共通出力名",
    )
    parser.add_argument(
        "--band-halfwidths",
        type=_parse_halfwidths,
        default=DEFAULT_BAND_HALFWIDTHS,
        help="滞在窓の半幅 (既定: 1e-2,1e-3,1e-4,1e-5,1e-6)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        print(f"data-dir が存在しない: {data_dir}", file=sys.stderr)
        return 2

    discovery = discover(data_dir)
    if args.discover_only:
        _print_discovery(discovery)
        return 0

    counts = discovery["counts"]
    if args.require_complete and counts["complete"] != counts["expected"]:
        print(
            f"未完了: complete={counts['complete']}/{counts['expected']}, "
            f"raw_metadata_partial={counts['raw_complete_metadata_partial']}, "
            f"partial={counts['partial']}, missing={counts['missing']}。"
            "出力は生成しない。",
            file=sys.stderr,
        )
        return 3

    selected = discovery["_selected_private"]
    if not selected:
        print("解析可能な生 NPZ が一本もない。出力は生成しない。", file=sys.stderr)
        return 4

    json_path, csv_path, figure_paths = _prepare_output_paths(
        data_dir, args.output_prefix, args.overwrite
    )
    runs: list[dict[str, Any]] = []
    plot_series: dict[tuple[str, float], dict[str, np.ndarray]] = {}
    run_by_key: dict[tuple[str, float], dict[str, Any]] = {}
    for mode in MODES:
        for delta in CORE_DELTAS:
            candidate = selected.get((mode, delta))
            if candidate is None:
                continue
            run, plot = _analyze_candidate(
                candidate, mode, delta, tuple(args.band_halfwidths)
            )
            runs.append(run)
            run_by_key[(mode, delta)] = run
            plot_series[(mode, delta)] = plot

    csv_rows: list[dict[str, Any]] = []
    for entry in discovery["rows"]:
        key = (entry["mode"], entry["delta"])
        run = run_by_key.get(key)
        csv_rows.append(_flatten_run_for_csv(run) if run else _missing_csv_row(entry))

    public_discovery = {k: v for k, v in discovery.items() if not k.startswith("_")}
    result = {
        "schema": "seed_type_delta_sweeps_T42000_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {
            "script": Path(__file__).name,
            "sha256": _sha256(Path(__file__)),
        },
        "scope": public_discovery["scope"],
        "late_window": {
            "array_slice_zero_based": [T_TARGET // 2, T_TARGET],
            "reported_steps_one_based": [T_TARGET // 2 + 1, T_TARGET],
        },
        "references": REFERENCES,
        "reference_separations": _reference_separations(),
        "dwell_band_halfwidths": list(args.band_halfwidths),
        "claim_boundaries": [
            "4参照値は由来の異なる独立値であり、近い値でも合并しない。",
            "通過時刻・滞在・最近接は数値的観測であり、物理的対応の証明ではない。",
            "シード数の異なる mode の同一 delta 比較には、総パワーとコヒーレント和の交絡が残る。",
            "T=42000 での未到達は右打切りであり、長時間極限を否定しない。",
            "瞬間 CV と時間平均後 CV は異なる統計量であり、一方で他方を代用しない。",
        ],
        "discovery": public_discovery,
        "runs": runs,
        "outputs": {
            "json": json_path.name,
            "csv": csv_path.name,
            "figures": [p.name for p in figure_paths],
        },
    }

    # 解析を先に完了させ、その後に明示的な成果物だけを書く。
    json_path.write_text(
        json.dumps(_json_safe(result), ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    _write_csv(csv_path, [_json_safe(row) for row in csv_rows])
    _save_figures(figure_paths, plot_series)
    print(
        f"後処理完了: analyzed={len(runs)}/{counts['expected']} / "
        f"{json_path.name}, {csv_path.name}, 図{len(figure_paths)}枚"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
