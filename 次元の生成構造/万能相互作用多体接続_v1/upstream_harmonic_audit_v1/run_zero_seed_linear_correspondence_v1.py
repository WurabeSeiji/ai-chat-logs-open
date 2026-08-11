#!/usr/bin/env python3
"""追加事前登録に従う種なし線形インフレーション対応対照。"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT_PATH = HERE / "run_upstream_harmonic_robustness_v1.py"
spec = importlib.util.spec_from_file_location("harmonic_audit", AUDIT_PATH)
audit = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = audit
assert spec.loader is not None
spec.loader.exec_module(audit)


def posthoc_slow_full_reference(z: np.ndarray, wp0: np.ndarray) -> tuple[dict, np.ndarray]:
    """事前条件不成立後の診断：高速化しない親V2でE63を再走行する。"""
    dummy = np.zeros_like(z)
    c0, evens, _ = audit.initial_state("EO63-cl", 0.0, z, dummy, dummy)
    eng = audit.V2(audit.N, c0, wp0, vertex_on=True)
    fval = audit.fval_factory(z)
    fs = np.zeros(audit.T_LONG)
    crossing = None
    for t in range(audit.T_LONG):
        eng.step()
        z_even = np.sum(eng.C[:, evens], axis=1) / np.sqrt(len(evens))
        fs[t] = fval(z_even)
        if crossing is None and fs[t] > 0.05:
            crossing = t + 1
    return {"crossing": crossing, "final_closure_max": audit.pointwise_closure(eng.C)}, fs


def main() -> None:
    started = time.time()
    _, _, _, _, _, _, _, z, wp0 = audit.abl.build_init(audit.N, False)
    parent = audit.gen3.make_parent(audit.N, seed=2)
    csec = np.fft.fft(parent.relation_waves, axis=1) / audit.N
    old_seed = csec[:, 1] / np.linalg.norm(csec[:, 1])
    closed_seed = audit.make_closed_seed(z, old_seed)

    single, s_series = audit.run_one(
        "L-cl", 0.0, z, old_seed, closed_seed, wp0, audit.T_LONG)
    full, f_series = audit.run_one(
        "EO63-cl", 0.0, z, old_seed, closed_seed, wp0, audit.T_LONG)
    crossing = single["crossing_even_aggregate"]
    compare_hi = min(crossing or audit.T_LONG, audit.T_LONG)
    max_abs_f = float(np.max(np.abs(
        s_series["f_even"][:compare_hi] - f_series["f_even"][:compare_hi]
    )))
    criteria = {
        "crossing_both_1166": bool(
            single["crossing_even_aggregate"] == 1166
            and full["crossing_even_aggregate"] == 1166),
        "pre_crossing_f_max_abs_diff_le_1e-10": bool(max_abs_f <= 1e-10),
        "full_initial_closure_le_1e-12": bool(full["initial_closure_max"] <= 1e-12),
        "full_closure_change_le_1e-10": bool(full["closure_change_max"] <= 1e-10),
    }
    slow_diag, slow_full_f = posthoc_slow_full_reference(z, wp0)
    slow_diag["purpose"] = "posthoc diagnosis after preregistered trajectory criterion failed"
    slow_diag["fast_vs_slow_pre_crossing_max_abs_difference"] = float(np.max(np.abs(
        f_series["f_even"][:compare_hi] - slow_full_f[:compare_hi]
    )))
    slow_diag["single_vs_slow_pre_crossing_max_abs_difference"] = float(np.max(np.abs(
        s_series["f_even"][:compare_hi] - slow_full_f[:compare_hi]
    )))
    out = {
        "preregistration": "追加事前登録_種なし線形インフレーション対応_v1.md",
        "single": single,
        "full_even63": full,
        "pre_crossing_f_max_abs_difference": max_abs_f,
        "criteria": criteria,
        "all_pass": bool(all(criteria.values())),
        "posthoc_slow_engine_diagnosis": slow_diag,
        "runtime_sec": time.time() - started,
    }
    (HERE / "zero_seed_linear_correspondence_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    np.savez_compressed(
        HERE / "zero_seed_linear_correspondence_series_v1.npz",
        single_f=s_series["f_even"], full_even63_f=f_series["f_even"],
        full_even63_slow_reference_f=slow_full_f,
        single_closure=s_series["closure"], full_even63_closure=f_series["closure"],
    )
    print(json.dumps({"crossing_single": single["crossing_even_aggregate"],
                      "crossing_full": full["crossing_even_aggregate"],
                      "pre_crossing_f_max_abs_difference": max_abs_f,
                      "criteria": criteria, "all_pass": out["all_pass"],
                      "posthoc_slow_engine_diagnosis": slow_diag,
                      "runtime_sec": out["runtime_sec"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
