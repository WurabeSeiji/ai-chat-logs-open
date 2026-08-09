#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analyze passive pre-vertex relation-edge R_e probes.

The analysis is deliberately limited to root *visitation/visibility*:
distance, crossings, dwell, and instantaneous cyclotomic phase defect.  It
does not call these quantities a Jacobian, monodromy, or proof of
``U_rho^n = I`` in the N-body mother dynamics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any, Optional

import numpy as np


TOLERANCES = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)


def references() -> dict[str, dict[str, Any]]:
    rho_124 = math.cos(23.0 * math.pi / 124.0) ** 2
    rho_620 = math.cos(117.0 * math.pi / 620.0) ** 2
    rho_mz = 0.687822933884774
    rho_122 = math.cos(23.0 * math.pi / 122.0) ** 2
    return {
        "rho_124_23": {
            "value": rho_124, "kind": "finite_order_root", "n": 124, "m": 23,
        },
        "rho_620_117": {
            "value": rho_620, "kind": "finite_order_root", "n": 620, "m": 117,
        },
        "rho_MZ_physical": {
            "value": rho_mz,
            "kind": "physical_correspondence_not_finite_root",
            "n": None,
            "m": None,
        },
        "rho_122_23_neighbor_root": {
            "value": rho_122, "kind": "neighbor_finite_order_control", "n": 122, "m": 23,
        },
        "sham_124_minus_1e-3": {
            "value": rho_124 - 1e-3, "kind": "symmetric_numeric_sham", "n": None, "m": None,
        },
        "sham_124_plus_1e-3": {
            "value": rho_124 + 1e-3, "kind": "symmetric_numeric_sham", "n": None, "m": None,
        },
        "sham_620_minus_1e-3": {
            "value": rho_620 - 1e-3, "kind": "symmetric_numeric_sham", "n": None, "m": None,
        },
        "sham_620_plus_1e-3": {
            "value": rho_620 + 1e-3, "kind": "symmetric_numeric_sham", "n": None, "m": None,
        },
        "sham_midpoint_620_MZ": {
            "value": 0.5 * (rho_620 + rho_mz),
            "kind": "pair_resolution_midpoint_sham",
            "n": None,
            "m": None,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze Stage-3 pre-vertex edge-R probes without recurrence overclaim"
    )
    parser.add_argument("probe", nargs="+", type=Path, help="probe NPZ files")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stem", default="finite_order_root_probe_stage3_analysis_v1")
    parser.add_argument(
        "--ordinary-npz",
        action="append",
        type=Path,
        default=[],
        help="Optional mother NPZ in the same order as probe files; r_nopump stays separate",
    )
    return parser.parse_args()


def longest_true_run(mask: np.ndarray) -> int:
    if not np.any(mask):
        return 0
    padded = np.concatenate(([False], mask.astype(bool), [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return int(np.max(edges[1::2] - edges[::2]))


def crossing_summary(values: np.ndarray, center: float) -> dict[str, Any]:
    """Strict side-to-side crossings, independently for each relation edge.

    Exact-center plateaus are collapsed before comparing signs, so a tangency
    such as below -> center -> below is not mislabeled as a crossing.  A NaN
    breaks continuity and is never bridged.
    """

    if values.ndim != 2:
        raise ValueError("crossing_summary expects shape (time, series)")
    up: list[tuple[float, int]] = []
    down: list[tuple[float, int]] = []
    for edge in range(values.shape[1]):
        previous_index: Optional[int] = None
        previous_sign = 0
        for index, raw in enumerate(values[:, edge]):
            value = float(raw)
            if not math.isfinite(value):
                previous_index = None
                previous_sign = 0
                continue
            sign = -1 if value < center else (1 if value > center else 0)
            if sign == 0:
                continue
            if previous_index is not None and sign != previous_sign:
                y0 = float(values[previous_index, edge])
                if index == previous_index + 1:
                    fraction = (center - y0) / (value - y0)
                    crossing_step = float(previous_index + 1 + fraction)
                else:
                    # All intervening finite samples equal center; first plateau
                    # sample is an exact observed contact.
                    crossing_step = float(previous_index + 2)
                target = up if previous_sign < sign else down
                target.append((crossing_step, edge))
            previous_index = index
            previous_sign = sign
    return {
        "up_count": len(up),
        "down_count": len(down),
        "edges_with_up": len({x[1] for x in up}),
        "edges_with_down": len({x[1] for x in down}),
        "first_up_interpolated_step": min((x[0] for x in up), default=None),
        "last_up_interpolated_step": max((x[0] for x in up), default=None),
        "first_down_interpolated_step": min((x[0] for x in down), default=None),
        "last_down_interpolated_step": max((x[0] for x in down), default=None),
        "interpolation_note": (
            "strict side-to-side crossing; exact-center plateaus are collapsed, "
            "tangencies are excluded, and linear interpolation is not a substep measurement"
        ),
    }


def phase_diagnostics_at_target_nearest(
    observed: float,
    center: float,
    n: int,
    step: int,
    edge: int,
) -> dict[str, Any]:
    """Phase diagnostics at the sample already selected by distance to rho_(n,m).

    Minimizing the n-fold closure defect over all samples would be wrong here:
    every m with the same n is a zero of that defect.  Selecting by distance to
    the named rho first prevents a different n-th root from being reported as
    evidence for this target.
    """

    omega = math.pi + 2.0 * math.asin(math.sqrt(min(max(observed, 0.0), 1.0)))
    omega_target = math.pi + 2.0 * math.asin(math.sqrt(center))
    return {
        "step": step,
        "edge_index": edge,
        "R_e": observed,
        "absolute_R_distance_to_selected_root": abs(observed - center),
        "single_step_phase_distance_to_selected_root": float(
            abs(np.exp(1j * omega) - np.exp(1j * omega_target))
        ),
        "n_fold_closure_defect": float(abs(np.exp(1j * n * omega) - 1.0)),
        "status": (
            "evaluated only at the sample nearest the selected rho_(n,m); "
            "the n-fold scalar defect alone is m-degenerate and is not an "
            "N-step operator product"
        ),
    }


def analyze_target(
    values: np.ndarray,
    edge_ia: np.ndarray,
    edge_ib: np.ndarray,
    name: str,
    ref: dict[str, Any],
) -> dict[str, Any]:
    center = float(ref["value"])
    distance = np.abs(values - center)
    flat = int(np.argmin(distance))
    t_near, edge_near = np.unravel_index(flat, distance.shape)
    dwell: dict[str, Any] = {}
    for tolerance in TOLERANCES:
        mask = distance <= tolerance
        any_edge = np.any(mask, axis=1)
        per_edge_longest = [longest_true_run(mask[:, e]) for e in range(mask.shape[1])]
        dwell[f"{tolerance:.1e}"] = {
            "tolerance": tolerance,
            "edge_step_count": int(np.count_nonzero(mask)),
            "time_steps_with_any_edge": int(np.count_nonzero(any_edge)),
            "edge_count_ever_inside": int(np.count_nonzero(np.any(mask, axis=0))),
            "longest_any_edge_time_run": longest_true_run(any_edge),
            "longest_single_edge_run": int(max(per_edge_longest, default=0)),
        }
    result = {
        "name": name,
        **ref,
        "nearest": {
            "absolute_distance": float(distance[t_near, edge_near]),
            "step": int(t_near + 1),
            "edge_index": int(edge_near),
            "edge_endpoints": [int(edge_ia[edge_near]), int(edge_ib[edge_near])],
            "R_e": float(values[t_near, edge_near]),
        },
        "crossings": crossing_summary(values, center),
        "dwell": dwell,
    }
    if ref.get("n") is not None:
        result["phase_diagnostics_at_target_nearest_sample"] = (
            phase_diagnostics_at_target_nearest(
                float(values[t_near, edge_near]),
                center,
                int(ref["n"]),
                int(t_near + 1),
                int(edge_near),
            )
        )
    else:
        result["phase_diagnostics_at_target_nearest_sample"] = None
    return result


def series_summary(values: np.ndarray) -> dict[str, Any]:
    return {
        "shape": list(values.shape),
        "finite_count": int(np.count_nonzero(np.isfinite(values))),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "initial_edge_mean": float(np.mean(values[0])),
        "final_edge_mean": float(np.mean(values[-1])),
        "late_edge_time_median": float(np.median(values[len(values) // 2 :])),
    }


def analyze_scalar_target(
    values: np.ndarray,
    name: str,
    ref: dict[str, Any],
    timing_note: str,
) -> dict[str, Any]:
    """Analyze a one-dimensional diagnostic without edge semantics."""

    center = float(ref["value"])
    finite = np.isfinite(values)
    distance = np.where(finite, np.abs(values - center), np.inf)
    if np.any(finite):
        nearest_index = int(np.argmin(distance))
        nearest = {
            "absolute_distance": float(distance[nearest_index]),
            "step": nearest_index + 1,
            "value": float(values[nearest_index]),
        }
    else:
        nearest = {"absolute_distance": None, "step": None, "value": None}
    dwell: dict[str, Any] = {}
    for tolerance in TOLERANCES:
        mask = finite & (distance <= tolerance)
        dwell[f"{tolerance:.1e}"] = {
            "tolerance": tolerance,
            "sample_count": int(np.count_nonzero(mask)),
            "time_steps_inside": int(np.count_nonzero(mask)),
            "longest_run": longest_true_run(mask),
        }
    return {
        "name": name,
        **ref,
        "nearest": nearest,
        "crossings": crossing_summary(values[:, None], center),
        "dwell": dwell,
        "timing_note": timing_note,
    }


def finite_or_none(value: float) -> Optional[float]:
    return float(value) if math.isfinite(float(value)) else None


def summarize_scalar_r_nopump(
    arr: np.ndarray,
    refs: dict[str, dict[str, Any]],
    timing_note: str,
) -> dict[str, Any]:
    finite = np.isfinite(arr)
    finite_indices = np.flatnonzero(finite)
    if finite_indices.size:
        maximum_index = int(finite_indices[np.argmax(arr[finite])])
        late = arr[len(arr) // 2 :]
        late = late[np.isfinite(late)]
        maximum: Optional[float] = float(arr[maximum_index])
        argmax_step: Optional[int] = maximum_index + 1
        late_median: Optional[float] = (
            float(np.median(late)) if late.size else None
        )
    else:
        maximum = None
        argmax_step = None
        late_median = None
    return {
        "length": int(arr.size),
        "finite_count": int(finite_indices.size),
        "initial": finite_or_none(arr[0]),
        "maximum": maximum,
        "argmax_step": argmax_step,
        "final": finite_or_none(arr[-1]),
        "late_median": late_median,
        "timing_note": timing_note,
        "targets": [
            analyze_scalar_target(arr, name, ref, timing_note)
            for name, ref in refs.items()
        ],
    }


def load_separate_r_nopump(
    path: Optional[Path], refs: dict[str, dict[str, Any]], expected_length: int
) -> Optional[dict[str, Any]]:
    if path is None:
        return None
    with np.load(path, allow_pickle=False) as data:
        result: dict[str, Any] = {"source": str(path), "definition": (
            "global odd/(odd+even_nonpump) post-step passive diagnostic; "
            "not relation-edge R_e"
        )}
        for channel in ("m", "v"):
            key = f"rec_{channel}_r_nopump"
            if key not in data.files:
                result[channel] = None
                continue
            arr = np.asarray(data[key], dtype=float)
            if arr.shape != (expected_length,):
                raise AssertionError(
                    f"{path}: {key} shape {arr.shape} does not match probe time axis "
                    f"({expected_length},)"
                )
            result[channel] = summarize_scalar_r_nopump(
                arr,
                refs,
                "post-step diagnostic; not simultaneous with pre-vertex R_e",
            )
        return result


def rows_from_payload(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def flat_row(
        probe_path: str,
        channel: str,
        quantity: str,
        item: dict[str, Any],
        nearest_edge: Optional[int],
        nearest_edge_ia: Optional[int],
        nearest_edge_ib: Optional[int],
        nearest_value: Optional[float],
        dwell_kind: str,
    ) -> dict[str, Any]:
        row = {
            "probe": probe_path,
            "channel": channel,
            "quantity": quantity,
            "target": item["name"],
            "target_kind": item["kind"],
            "center": item["value"],
            "n": item["n"],
            "m": item["m"],
            "nearest_abs_distance": item["nearest"]["absolute_distance"],
            "nearest_step": item["nearest"]["step"],
            "nearest_edge": nearest_edge,
            "nearest_edge_ia": nearest_edge_ia,
            "nearest_edge_ib": nearest_edge_ib,
            "nearest_value": nearest_value,
            "up_count": item["crossings"]["up_count"],
            "down_count": item["crossings"]["down_count"],
            "dwell_observation_unit": dwell_kind,
        }
        for tolerance in TOLERANCES:
            d = item["dwell"][f"{tolerance:.1e}"]
            if dwell_kind == "edge-step":
                observations = d["edge_step_count"]
                time_steps = d["time_steps_with_any_edge"]
            else:
                observations = d["sample_count"]
                time_steps = d["time_steps_inside"]
            row[f"observations_within_{tolerance:.0e}"] = observations
            row[f"time_steps_within_{tolerance:.0e}"] = time_steps
        return row

    for probe in payload["probes"]:
        for channel in ("matter", "vacuum"):
            for item in probe["channels"][channel]["targets"]:
                rows.append(flat_row(
                    probe["probe"],
                    channel,
                    "prevertex_relation_edge_R_e",
                    item,
                    item["nearest"]["edge_index"],
                    item["nearest"]["edge_endpoints"][0],
                    item["nearest"]["edge_endpoints"][1],
                    item["nearest"]["R_e"],
                    "edge-step",
                ))
        for source_key, quantity in (
            ("prevertex_r_nopump", "prevertex_r_nopump"),
            ("ordinary_r_nopump", "poststep_r_nopump"),
        ):
            scalar_source = probe[source_key]
            if scalar_source is None:
                continue
            for short_channel, channel in (("m", "matter"), ("v", "vacuum")):
                scalar = scalar_source[short_channel]
                if scalar is None:
                    continue
                for item in scalar["targets"]:
                    rows.append(flat_row(
                        probe["probe"],
                        channel,
                        quantity,
                        item,
                        None,
                        None,
                        None,
                        item["nearest"]["value"],
                        "time-step",
                    ))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("no rows")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,127}", args.stem):
        raise ValueError("--stem must be 1..128 ASCII letters/digits/_/-")
    if args.ordinary_npz and len(args.ordinary_npz) != len(args.probe):
        raise ValueError("--ordinary-npz must be omitted or supplied once per probe")
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    refs = references()
    probes: list[dict[str, Any]] = []
    for index, path_arg in enumerate(args.probe):
        path = path_arg.resolve()
        ordinary = args.ordinary_npz[index].resolve() if args.ordinary_npz else None
        with np.load(path, allow_pickle=False) as data:
            required = {
                "step", "edge_index", "edge_ia", "edge_ib",
                "R_prevertex_matter", "R_prevertex_vacuum",
                "r_nopump_prevertex_matter", "r_nopump_prevertex_vacuum",
                "odd_power_prevertex_matter", "odd_power_prevertex_vacuum",
                "even_power_nopump_prevertex_matter",
                "even_power_nopump_prevertex_vacuum",
                "readout_calls_per_step_matter", "readout_calls_per_step_vacuum",
            }
            missing = sorted(required - set(data.files))
            if missing:
                raise KeyError(f"{path}: missing {missing}")
            step = np.asarray(data["step"], dtype=np.int64)
            edge_index = np.asarray(data["edge_index"], dtype=np.int64)
            edge_ia = np.asarray(data["edge_ia"], dtype=np.int64)
            edge_ib = np.asarray(data["edge_ib"], dtype=np.int64)
            if not np.array_equal(step, np.arange(1, len(step) + 1)):
                raise AssertionError(f"{path}: non-contiguous step axis")
            if not np.array_equal(edge_index, np.arange(len(edge_index))):
                raise AssertionError(f"{path}: non-contiguous edge index")
            if len(edge_ia) != len(edge_index) or len(edge_ib) != len(edge_index):
                raise AssertionError(f"{path}: inconsistent edge endpoint lengths")
            for channel in ("matter", "vacuum"):
                calls = np.asarray(data[f"readout_calls_per_step_{channel}"])
                if calls.shape != step.shape or not np.all(calls == 1):
                    raise AssertionError(f"{path}: bad {channel} pre-vertex call audit")
            channels: dict[str, Any] = {}
            synchronized_r_nopump: dict[str, Any] = {
                "source": str(path),
                "definition": (
                    "global odd/(odd+even_nonpump) diagnostic copied from the "
                    "same pre-vertex state as R_e; it remains a separate aggregate"
                ),
            }
            for channel in ("matter", "vacuum"):
                values = np.asarray(data[f"R_prevertex_{channel}"], dtype=float)
                if values.shape != (len(step), len(edge_ia)):
                    raise AssertionError(f"{path}: bad {channel} shape {values.shape}")
                if not np.all(np.isfinite(values)) or not np.all((0 <= values) & (values <= 1)):
                    raise AssertionError(f"{path}: invalid {channel} R_e")
                channels[channel] = {
                    "series": series_summary(values),
                    "targets": [
                        analyze_target(values, edge_ia, edge_ib, name, ref)
                        for name, ref in refs.items()
                    ],
                }
                rnp = np.asarray(data[f"r_nopump_prevertex_{channel}"], dtype=float)
                odd_power = np.asarray(data[f"odd_power_prevertex_{channel}"], dtype=float)
                even_np = np.asarray(
                    data[f"even_power_nopump_prevertex_{channel}"], dtype=float
                )
                for array_name, array in (
                    ("r_nopump", rnp),
                    ("odd_power", odd_power),
                    ("even_power_nopump", even_np),
                ):
                    if array.shape != step.shape:
                        raise AssertionError(
                            f"{path}: bad {channel} {array_name} shape {array.shape}"
                        )
                if not np.all(np.isfinite(odd_power)) or not np.all(odd_power >= 0.0):
                    raise AssertionError(f"{path}: invalid {channel} odd power")
                if not np.all(np.isfinite(even_np)) or not np.all(even_np >= 0.0):
                    raise AssertionError(f"{path}: invalid {channel} even nonpump power")
                finite_rnp = rnp[np.isfinite(rnp)]
                if not np.all((0.0 <= finite_rnp) & (finite_rnp <= 1.0)):
                    raise AssertionError(f"{path}: invalid {channel} pre-vertex r_nopump")
                short_channel = "m" if channel == "matter" else "v"
                synchronized_r_nopump[short_channel] = summarize_scalar_r_nopump(
                    rnp,
                    refs,
                    "same pre-vertex state as R_e; synchronized passive aggregate",
                )
            probes.append({
                "probe": str(path),
                "prevertex_r_nopump": synchronized_r_nopump,
                "ordinary_r_nopump": load_separate_r_nopump(
                    ordinary, refs, len(step)
                ),
                "edge_definition": "relation edge e=(edge_ia[e],edge_ib[e]); not graph node",
                "channels": channels,
            })
    rho_620 = refs["rho_620_117"]["value"]
    rho_mz = refs["rho_MZ_physical"]["value"]
    payload = {
        "schema": "finite_order_root_probe_stage3_analysis_v1",
        "claim_boundary": {
            "measured": (
                "pre-vertex relation-edge scalar R_e visitation plus a synchronized "
                "but definitionally separate global r_nopump aggregate"
            ),
            "not_measured": [
                "a two-channel invariant projection",
                "Jacobian or monodromy",
                "N-body realization of U_rho^n=I",
            ],
            "instantaneous_defect_note": (
                "phase diagnostics are evaluated at the sample selected nearest "
                "the named rho_(n,m); n-fold closure alone is m-degenerate"
            ),
        },
        "references": refs,
        "rho_620_minus_MZ": rho_620 - rho_mz,
        "half_separation_620_MZ": 0.5 * (rho_620 - rho_mz),
        "tolerances": list(TOLERANCES),
        "probes": probes,
    }
    json_path = outdir / f"{args.stem}.json"
    csv_path = outdir / f"{args.stem}.csv"
    if json_path.exists() or csv_path.exists():
        raise FileExistsError(f"analysis output already exists for stem {args.stem!r}")
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    write_csv(csv_path, rows_from_payload(payload))
    print(f"saved: {json_path}")
    print(f"saved: {csv_path}")


if __name__ == "__main__":
    main()
