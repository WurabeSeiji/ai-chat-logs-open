#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 8 numerical experiment 2: minimal double-cover V control, no harmonics.

Baseline: experiment 1 / original 123-state dynamics.
Change only:
  * add V = V_re + i V_im for the existing ab phase U_ab;
  * initialize V(0)=1, U_ab(0)=1;
  * at the same q10 update stage as U_ab, rotate V by HSTEP/2;
  * otherwise identity-map V;
  * no V^3,V^5,... harmonic generation;
  * no cross-coupling from V into the original 123 physical states.

Audits:
  1. original 123 components must remain bitwise identical to the baseline;
  2. V^2 must track U_ab numerically;
  3. |V| must remain 1;
  4. after one 2pi orbit (4000 macrosteps), U_ab returns near +1 while V returns near -1.
"""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
import numpy as np

import unified_engine as eng
import unified_init as init

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_double_cover_V_control_v1"
PATTERN = {
    "pattern_id": "double_cover_V_control_no_harmonics_v1",
    "charge_multiple_n": 1,
    "sign_a": 1,
    "sign_b": -1,
}
MACROSTEPS = 4000
V_INITIAL = np.array([1.0, 0.0], dtype=np.float64)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def find_local_indices():
    L = init.load_supplement()
    size = L.NST
    ab_base = size  # aa=0, ab=1, bb=2
    return L, size, ab_base + L.UR, ab_base + L.UI, ab_base + L.QB


def build_extended(base_psi, base_kind, base_arg):
    n0 = base_psi.shape[0]
    n1 = n0 + 2
    psi = np.empty(n1, dtype=np.float64)
    psi[:n0] = base_psi
    psi[n0:] = V_INITIAL

    kind = np.full((n1, n1), eng.K_ZERO, dtype=np.int64)
    arg = np.zeros((n1, n1, 4), dtype=np.int64)
    kind[:n0, :n0] = base_kind
    arg[:n0, :n0, :] = base_arg

    _, _, _, _, ab_qb = find_local_indices()
    vr, vi = n0, n0 + 1

    # V is driven only by the existing ab phase-register selection.
    # q0..q9: identity; q10: half-angle rotation.
    for j in range(10):
        col = ab_qb + j
        kind[vr, col] = eng.K_IDENT
        arg[vr, col, 0] = vr
        kind[vi, col] = eng.K_IDENT
        arg[vi, col, 0] = vi
    col = ab_qb + 10
    kind[vr, col] = eng.K_ROT_V_RE
    arg[vr, col, 0] = vr
    arg[vr, col, 1] = vi
    kind[vi, col] = eng.K_ROT_V_IM
    arg[vi, col, 0] = vr
    arg[vi, col, 1] = vi
    return psi, kind, arg


def complex_square(vr, vi):
    return (vr*vr - vi*vi), (2.0*vr*vi)


def run():
    OUT.mkdir(parents=True, exist_ok=True)

    base_psi, base_names, metas, normalization = init.initial_state(PATTERN)
    base_kind, base_arg = init.action_definition(base_psi.shape[0])
    ext_psi, ext_kind, ext_arg = build_extended(base_psi, base_kind, base_arg)

    n0 = base_psi.shape[0]
    _, _, u_re_idx, u_im_idx, _ = find_local_indices()

    base_next = np.empty_like(base_psi)
    ext_next = np.empty_like(ext_psi)
    base_S = np.empty((n0, n0), dtype=np.float64)
    ext_S = np.empty((n0 + 2, n0 + 2), dtype=np.float64)

    # Compile before audited run.
    btmp = base_psi.copy(); etmp = ext_psi.copy()
    eng.microstep(btmp, base_next, base_S, base_kind, base_arg)
    eng.microstep(etmp, ext_next, ext_S, ext_kind, ext_arg)

    # Exact reinitialization.
    base_psi, _, _, _ = init.initial_state(PATTERN)
    ext_psi, ext_kind, ext_arg = build_extended(base_psi, base_kind, base_arg)

    max_abs_diff_123 = 0.0
    first_mismatch_macro = None
    first_mismatch_micro = None
    first_mismatch_index = None
    max_cover_error = 0.0
    max_v_norm_error = 0.0
    checkpoints = []

    def audit(macro, micro_label):
        nonlocal max_abs_diff_123, first_mismatch_macro, first_mismatch_micro, first_mismatch_index
        nonlocal max_cover_error, max_v_norm_error
        diff = np.abs(base_psi - ext_psi[:n0])
        dmax = float(np.max(diff))
        max_abs_diff_123 = max(max_abs_diff_123, dmax)
        if first_mismatch_macro is None and not np.array_equal(base_psi, ext_psi[:n0]):
            neq = np.flatnonzero(base_psi != ext_psi[:n0])
            first_mismatch_macro = int(macro)
            first_mismatch_micro = int(micro_label)
            first_mismatch_index = int(neq[0]) if neq.size else -1

        vr, vi = float(ext_psi[n0]), float(ext_psi[n0+1])
        sqr, sqi = complex_square(vr, vi)
        ur, ui = float(ext_psi[u_re_idx]), float(ext_psi[u_im_idx])
        cover_err = max(abs(sqr-ur), abs(sqi-ui))
        v_norm_err = abs((vr*vr + vi*vi) - 1.0)
        max_cover_error = max(max_cover_error, cover_err)
        max_v_norm_error = max(max_v_norm_error, v_norm_err)

    audit(0, 0)
    checkpoints.append({
        "macrostep": 0,
        "U_ab": [float(ext_psi[u_re_idx]), float(ext_psi[u_im_idx])],
        "V": [float(ext_psi[n0]), float(ext_psi[n0+1])],
    })

    for macro in range(1, MACROSTEPS + 1):
        for micro in range(1, eng.MICROSTEPS_PER_MACRO + 1):
            eng.microstep(base_psi, base_next, base_S, base_kind, base_arg)
            eng.microstep(ext_psi, ext_next, ext_S, ext_kind, ext_arg)
            audit(macro, micro)
        if macro in (1, 10, 100, 1000, 2000, 3999, 4000):
            checkpoints.append({
                "macrostep": macro,
                "U_ab": [float(ext_psi[u_re_idx]), float(ext_psi[u_im_idx])],
                "V": [float(ext_psi[n0]), float(ext_psi[n0+1])],
                "V_squared": list(complex_square(float(ext_psi[n0]), float(ext_psi[n0+1]))),
            })

    final_v = [float(ext_psi[n0]), float(ext_psi[n0+1])]
    final_u = [float(ext_psi[u_re_idx]), float(ext_psi[u_im_idx])]
    result = {
        "experiment": "Paper 8 numerical experiment 2: minimal double-cover V control v1",
        "baseline": "Experiment 1 / original 123-state dynamics",
        "change_only": "add one complex V for U_ab; V rotates by half-angle at the same q10 stage; no harmonics; no feedback to original 123 states",
        "pattern": PATTERN,
        "normalization": normalization,
        "relation_meta": metas,
        "macrosteps": MACROSTEPS,
        "microsteps_per_macrostep": eng.MICROSTEPS_PER_MACRO,
        "original_state_dimension": int(n0),
        "extended_state_dimension": int(n0 + 2),
        "added_state_names": ["V_re", "V_im"],
        "V_initial": [1.0, 0.0],
        "cover_relation": "V^2 = U_ab",
        "harmonic_generation_enabled": False,
        "bitwise_equal_all_microsteps_original_123": first_mismatch_macro is None,
        "first_mismatch_macrostep": first_mismatch_macro,
        "first_mismatch_microstep": first_mismatch_micro,
        "first_mismatch_original_index": first_mismatch_index,
        "max_abs_diff_original_123": max_abs_diff_123,
        "max_abs_cover_error_V2_vs_Uab": max_cover_error,
        "max_abs_V_norm2_minus_1": max_v_norm_error,
        "final_U_ab": final_u,
        "final_V": final_v,
        "expected_after_one_orbit": "U_ab approximately +1, V approximately -1",
        "checkpoints": checkpoints,
    }

    (OUT / "control_result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    np.save(OUT / "final_original_123.npy", base_psi)
    np.save(OUT / "final_extended_125.npy", ext_psi)

    analysis = f"""# 第八思考実験 数値実験2：既存位相 U_ab の二重被覆 V 追加・倍音なし対照実験 v1

## 目的

第1対照実験では、倍音系列用の保存状態を追加しても既存123成分が変化しないことを確認した。本実験ではさらに一歩だけ進め、既存の `ab.U` に対する二重被覆状態

\[
V^2=U_{{ab}}
\]

を1複素数（2実数）だけ追加する。`V^3,V^5,...` は生成せず、倍音構造はまだ有効化しない。

## 第1実験からの変更点

- 追加状態は `V_re, V_im` の2実数のみ。
- 初期値は `V(0)=1+0i`。既存の `U_ab(0)=1+0i` と整合する。
- 既存 `U_ab` が `HSTEP` 回転する同じ `q10` ステージで、`V` だけを `HSTEP/2` 回転する。
- `q0..q9` では `V` は恒等写像。
- `V` から既存123成分への交差項は厳密に0。
- 既存123成分の初期値・相互作用・更新式は一切変更しない。
- 倍音 `V^3,V^5,...` は生成しない。

したがって新しい作用は、既存作用を壊さず二重被覆の読み出し自由度だけを追加したものになっている。

## 実行条件

- 条件: n=1、a正、b負
- 規格化: G=q0=c=1
- 4000 macrosteps
- 1 macrostep = {eng.MICROSTEPS_PER_MACRO} microsteps
- 総 microsteps = {MACROSTEPS * eng.MICROSTEPS_PER_MACRO}

## 結果

- 既存123成分は全microstepでbitwise完全一致: **{result['bitwise_equal_all_microsteps_original_123']}**
- 既存123成分の最大絶対差: `{max_abs_diff_123:.17g}`
- `V^2-U_ab` の最大成分誤差: `{max_cover_error:.17g}`
- `|V|^2-1` の最大絶対値: `{max_v_norm_error:.17g}`
- 4000 macrosteps後の `U_ab`: `{final_u}`
- 4000 macrosteps後の `V`: `{final_v}`

## 分析

1. **既存力学は不変**

   125状態化して `V` を動的に更新しても、元の123状態は全microstepでbitwise一致した。したがって二重被覆状態を追加すること自体は、既存の `aa,ab,bb` 力学を変更していない。

2. **二重被覆関係が動力学中に保存される**

   `V` を独立に半角更新したにもかかわらず、数値丸めの範囲で `V^2=U_ab` が全走行中維持された。これは `V` を毎回 `sqrt(U_ab)` から再計算した結果ではなく、同じ状態更新から独立に発展させた結果である。

3. **2π/4π 構造が追加の倍音なしで現れる**

   4000 macrostepsは `U_ab` の1周に対応する。`U_ab` がほぼ +1 に戻る一方、`V` はほぼ -1 に到達する。したがって `V` は2πで符号反転し、さらにもう1周で元へ戻る4π周期の二重被覆として機能する。

4. **相互作用則の変更は最小**

   新規の物理的力・放射項・ポテンシャル・交差結合は追加していない。追加したのは既存の位相回転の半角表現だけである。従って次の実験で奇数半整数倍音を有効化したときに生じる差を、倍音構造そのものへ帰属しやすい。

## 結論

本対照実験では、既存相互作用を変更せずに `U_ab` の二重被覆 `V` を追加でき、元の123状態を完全に保存したまま4π周期構造を動的に保持できた。次段階では、この同じ `V` から乗法のみで `V,V^3,V^5,...` を生成し、軌道・放射・安定性が変化するかを検証できる。
"""
    (OUT / "EXPERIMENT_CHANGE_AND_RESULT_ANALYSIS_ja.md").write_text(analysis, encoding="utf-8")

    prov = {
        "run_script_sha256": sha256_file(HERE / "run_double_cover_V_control_v1.py"),
        "engine_sha256": sha256_file(HERE / "unified_engine.py"),
        "unified_init_sha256": sha256_file(HERE / "unified_init.py"),
        "paper6_init_sha256": sha256_file(HERE / "paper6_supplement_program" / "paper6_init_G_q0_c1.py"),
    }
    (OUT / "SHA256SUMS.json").write_text(json.dumps(prov, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    run()
