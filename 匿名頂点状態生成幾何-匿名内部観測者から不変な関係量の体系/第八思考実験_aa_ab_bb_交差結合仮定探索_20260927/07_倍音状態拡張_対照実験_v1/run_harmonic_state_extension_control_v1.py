#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 8 harmonic-state extension control experiment v1.

Purpose
-------
Extend the existing 123-component state Psi=(aa,ab,bb) by the minimum
harmonic-series state without enabling harmonics, then verify that the
original 123 physical components evolve identically.

Added persistent state (4 real components = 2 complex conserved values):
    B0 = B0_re + i B0_im   : series seed
    R  = R_re  + i R_im    : multiplicative series generator

Control condition (no harmonic structure enabled):
    B0 = 1 + 0 i
    R  = 1 + 0 i

The four new components are identity-mapped at every microstep and have no
cross-coupling to the original 123 components.  The existing physics and
state-dependent S(Psi) definitions are unchanged.
"""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
import numpy as np

import unified_engine as eng
import unified_init as init

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_harmonic_extension_control_v1"
PATTERN = {
    "pattern_id": "harmonic_extension_control_no_harmonics_v1",
    "charge_multiple_n": 1,
    "sign_a": 1,
    "sign_b": -1,
}
MACROSTEPS = 4000
HARMONIC_NAMES = ("B0_re", "B0_im", "R_re", "R_im")
HARMONIC_INITIAL = np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def build_extended(base_psi, base_kind, base_arg):
    n0 = base_psi.shape[0]
    n1 = n0 + 4
    psi = np.empty(n1, dtype=np.float64)
    psi[:n0] = base_psi
    psi[n0:] = HARMONIC_INITIAL

    kind = np.full((n1, n1), eng.K_ZERO, dtype=np.int64)
    arg = np.zeros((n1, n1, 4), dtype=np.int64)
    kind[:n0, :n0] = base_kind
    arg[:n0, :n0, :] = base_arg

    # Identity map for the two complex harmonic labels. No cross-coupling.
    for i in range(n0, n1):
        kind[i, i] = eng.K_ONE
    return psi, kind, arg


def run():
    OUT.mkdir(parents=True, exist_ok=True)

    base_psi, base_names, metas, normalization = init.initial_state(PATTERN)
    base_kind, base_arg = init.action_definition(base_psi.shape[0])
    ext_psi, ext_kind, ext_arg = build_extended(base_psi, base_kind, base_arg)

    n0 = base_psi.shape[0]
    base_next = np.empty_like(base_psi)
    ext_next = np.empty_like(ext_psi)
    base_S = np.empty((n0, n0), dtype=np.float64)
    ext_S = np.empty((n0 + 4, n0 + 4), dtype=np.float64)

    # Compile before the audited run.
    btmp = base_psi.copy(); etmp = ext_psi.copy()
    for _ in range(eng.MICROSTEPS_PER_MACRO):
        eng.microstep(btmp, base_next, base_S, base_kind, base_arg)
        eng.microstep(etmp, ext_next, ext_S, ext_kind, ext_arg)

    # Reinitialize exactly for audited run.
    base_psi, _, _, _ = init.initial_state(PATTERN)
    ext_psi, ext_kind, ext_arg = build_extended(base_psi, base_kind, base_arg)

    max_abs_diff = 0.0
    first_mismatch_macro = None
    first_mismatch_index = None
    harmonic_drift = np.zeros(4, dtype=np.float64)
    checkpoints = []

    def audit(step):
        nonlocal max_abs_diff, first_mismatch_macro, first_mismatch_index, harmonic_drift
        diff = np.abs(base_psi - ext_psi[:n0])
        dmax = float(np.max(diff))
        if dmax > max_abs_diff:
            max_abs_diff = dmax
        if first_mismatch_macro is None and not np.array_equal(base_psi, ext_psi[:n0]):
            neq = np.flatnonzero(base_psi != ext_psi[:n0])
            first_mismatch_macro = int(step)
            first_mismatch_index = int(neq[0]) if neq.size else -1
        harmonic_drift = np.maximum(harmonic_drift, np.abs(ext_psi[n0:] - HARMONIC_INITIAL))
        if step in (0, 1, 10, 100, 1000, MACROSTEPS):
            checkpoints.append({
                "macrostep": int(step),
                "bitwise_equal_original_123": bool(np.array_equal(base_psi, ext_psi[:n0])),
                "max_abs_diff_original_123": dmax,
                "harmonic_state": [float(x) for x in ext_psi[n0:]],
            })

    audit(0)
    for macro in range(1, MACROSTEPS + 1):
        for _ in range(eng.MICROSTEPS_PER_MACRO):
            eng.microstep(base_psi, base_next, base_S, base_kind, base_arg)
            eng.microstep(ext_psi, ext_next, ext_S, ext_kind, ext_arg)
        audit(macro)

    result = {
        "experiment": "Paper 8 harmonic-state extension control v1",
        "purpose": "Verify that adding conserved harmonic-series state with harmonics disabled does not change the existing 123-component dynamics.",
        "pattern": PATTERN,
        "normalization": normalization,
        "relation_meta": metas,
        "macrosteps": MACROSTEPS,
        "microsteps_per_macrostep": eng.MICROSTEPS_PER_MACRO,
        "original_state_dimension": int(n0),
        "extended_state_dimension": int(n0 + 4),
        "added_state_names": list(HARMONIC_NAMES),
        "added_state_initial": [float(x) for x in HARMONIC_INITIAL],
        "added_state_update": "identity map; zero cross-coupling",
        "control_interpretation": "B0=1+0i, R=1+0i; harmonic structure disabled",
        "bitwise_equal_all_macrosteps_original_123": first_mismatch_macro is None,
        "first_mismatch_macrostep": first_mismatch_macro,
        "first_mismatch_original_index": first_mismatch_index,
        "max_abs_diff_original_123": max_abs_diff,
        "max_abs_drift_added_state": [float(x) for x in harmonic_drift],
        "final_added_state": [float(x) for x in ext_psi[n0:]],
        "checkpoints": checkpoints,
    }

    (OUT / "control_result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    np.save(OUT / "final_original_123.npy", base_psi)
    np.save(OUT / "final_extended_127.npy", ext_psi)

    summary = [
        "# 第八思考実験：倍音状態拡張の対照実験 v1",
        "",
        "## 目的",
        "",
        "現行の `Psi=(aa,ab,bb)` 123成分に、倍音系列を将来生成するための最小保存状態として二つの複素量 `B0`, `R`（4実数）だけを追加し、倍音を無効化した対照条件で既存123成分が変化しないことを検証する。",
        "",
        "## 追加状態",
        "",
        "- `B0 = B0_re + i B0_im`: 倍音系列の初期要素",
        "- `R = R_re + i R_im`: 乗法的系列生成子",
        "- 対照条件: `B0=1+0i`, `R=1+0i`",
        "- 更新: 4成分とも恒等写像",
        "- 既存123成分との交差項: 厳密に0",
        "",
        "## 結果",
        "",
        f"- 実行: {MACROSTEPS} macrosteps = {MACROSTEPS * eng.MICROSTEPS_PER_MACRO} microsteps",
        f"- 既存123成分が全macrostepでbitwise一致: **{result['bitwise_equal_all_macrosteps_original_123']}**",
        f"- 既存123成分の最大絶対差: `{max_abs_diff:.17g}`",
        f"- 追加4状態の最大drift: `{result['max_abs_drift_added_state']}`",
        f"- 最終追加状態: `{result['final_added_state']}`",
        "",
        "## 判定",
        "",
        "この対照条件で完全一致すれば、倍音系列用の保存状態を追加すること自体は既存の力学を変更しない。次段階では `R != 1` として乗法反復だけで倍音系列を有効化し、その効果だけを比較できる。",
    ]
    (OUT / "README_ja.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    # provenance hashes
    prov = {
        "run_script_sha256": sha256_file(HERE / "run_harmonic_state_extension_control_v1.py"),
        "unified_engine_sha256": sha256_file(HERE / "unified_engine.py"),
        "unified_init_sha256": sha256_file(HERE / "unified_init.py"),
        "paper6_init_sha256": sha256_file(HERE / "paper6_supplement_program" / "paper6_init_G_q0_c1.py"),
    }
    (OUT / "SHA256SUMS.json").write_text(json.dumps(prov, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    run()
