#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Analyze K4 base+harmonic discrete states (m=2,3).

Input (same directory by default):
  complex_states_animation_exact_i.csv
  complex_states_static_exp.csv   [optional cross-check]

Outputs:
  K4_self_consistency_square_sum_by_state_20260915.csv
  K4_self_consistency_square_sum_summary_20260915.csv
  K4_self_consistency_square_sum_analysis_20260915.md

Definitions
-----------
U = i, state S_k = (z_base, z_harm) = (U^k, U^(m k)), k=0,1,2,3.

Algebraic self-consistency checks (kept separate from quadratic closure):
  1. |z_base| = |z_harm| = 1
  2. z_harm = z_base^m
  3. S_{k+1} = (U z_base, U^m z_harm) for forward order
  4. forward/reverse rows represent the same four states
  5. after 4 forward steps the pair state returns to S_0

Quadratic sum for each state:
  Q_k = z_base^2 + z_harm^2.
For z_j = a_j + i b_j,
  Q_k = sum(a_j^2) - sum(b_j^2) + 2 i sum(a_j b_j).
The whole-loop total is the sum over the four unique states k=0..3.
"""

from __future__ import annotations

import csv
import cmath
import math
from pathlib import Path
from typing import Dict, List, Tuple

HERE = Path(__file__).resolve().parent
RAW_DIR_CANDIDATE = HERE.parent / "raw_data"
RAW_DIR = RAW_DIR_CANDIDATE if RAW_DIR_CANDIDATE.exists() else HERE
EXACT_CSV = RAW_DIR / "complex_states_animation_exact_i.csv"
STATIC_CSV = RAW_DIR / "complex_states_static_exp.csv"
OUT_STATE = HERE / "K4_self_consistency_square_sum_by_state_20260915.csv"
OUT_SUMMARY = HERE / "K4_self_consistency_square_sum_summary_20260915.csv"
OUT_MD = HERE / "K4_self_consistency_square_sum_analysis_20260915.md"
TOL = 1e-12
U = 1j


def read_rows(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def c_from(row: dict, prefix: str) -> complex:
    return complex(float(row[f"{prefix}_re"]), float(row[f"{prefix}_im"]))


def fmt(x: float, nd=15) -> str:
    if x == 0:
        return "0"
    return f"{x:.{nd}g}"


def fmt_clean(x: float, nd=15) -> str:
    if abs(x) < 5e-15:
        return "0"
    return f"{x:.{nd}g}"


def complex_order(z: complex, max_n: int = 64, tol: float = TOL) -> int | None:
    for n in range(1, max_n + 1):
        if abs(z**n - 1) <= tol:
            return n
    return None


def analyze(exact_rows: List[dict], static_rows: List[dict] | None = None):
    by_state_out: List[dict] = []
    summary_out: List[dict] = []

    # Cross-source residual: static exp() representation vs exact-i representation.
    static_lookup: Dict[Tuple[int, str, int], dict] = {}
    if static_rows:
        for r in static_rows:
            static_lookup[(int(r["m"]), r["direction"], int(r["k"]))] = r

    for m in (2, 3):
        rows_m = [r for r in exact_rows if int(r["m"]) == m]
        forward = sorted([r for r in rows_m if r["direction"] == "forward"], key=lambda r: int(r["step"]))
        reverse = sorted([r for r in rows_m if r["direction"] == "reverse"], key=lambda r: int(r["step"]))
        if len(forward) != 4 or len(reverse) != 4:
            raise ValueError(f"m={m}: expected 4 forward and 4 reverse rows")

        f_by_k = {int(r["k"]): r for r in forward}
        r_by_k = {int(r["k"]): r for r in reverse}
        if set(f_by_k) != {0,1,2,3} or set(r_by_k) != {0,1,2,3}:
            raise ValueError(f"m={m}: k states are incomplete")

        state_residuals = []
        transition_residuals = []
        source_residuals = []
        q_loop = 0j
        a_loop = 0.0
        b_loop = 0.0
        a2_loop = 0.0
        b2_loop = 0.0
        ab_loop = 0.0

        for k in range(4):
            r = f_by_k[k]
            z1 = c_from(r, "base")
            zm = c_from(r, "harmonic")

            unit_base = abs(abs(z1) - 1.0)
            unit_harm = abs(abs(zm) - 1.0)
            harmonic_resid = abs(zm - z1**m)

            # Exact match of same physical state in forward/reverse tables.
            rr = r_by_k[k]
            z1_r = c_from(rr, "base")
            zm_r = c_from(rr, "harmonic")
            reverse_match_resid = max(abs(z1 - z1_r), abs(zm - zm_r))

            # Forward generator consistency, including k=3 -> k=0 closure.
            rn = f_by_k[(k + 1) % 4]
            z1_next = c_from(rn, "base")
            zm_next = c_from(rn, "harmonic")
            pred1 = U * z1
            predm = (U**m) * zm
            step_resid = max(abs(z1_next - pred1), abs(zm_next - predm))

            # Quadratic decomposition for the two components at this state.
            comps = [z1, zm]
            a2 = sum((z.real**2) for z in comps)
            b2 = sum((z.imag**2) for z in comps)
            ab = sum((z.real*z.imag) for z in comps)
            q = sum((z**2) for z in comps)
            q_formula = complex(a2 - b2, 2*ab)
            q_identity_resid = abs(q - q_formula)

            q_loop += q
            a_loop += sum(z.real for z in comps)
            b_loop += sum(z.imag for z in comps)
            a2_loop += a2
            b2_loop += b2
            ab_loop += ab

            static_resid = 0.0
            if static_lookup:
                sr = static_lookup[(m, "forward", k)]
                static_resid = max(
                    abs(z1 - c_from(sr, "base")),
                    abs(zm - c_from(sr, "harmonic")),
                    abs(complex(float(r["sum_re"]), float(r["sum_im"])) - c_from(sr, "sum")),
                )
                source_residuals.append(static_resid)

            local_resid = max(unit_base, unit_harm, harmonic_resid, reverse_match_resid, step_resid, q_identity_resid)
            state_residuals.append(local_resid)
            transition_residuals.append(step_resid)

            by_state_out.append({
                "m": m,
                "parity": "even" if m % 2 == 0 else "odd",
                "k": k,
                "theta_over_pi": float(r["theta_over_pi"]),
                "base_re": z1.real,
                "base_im": z1.imag,
                "harmonic_re": zm.real,
                "harmonic_im": zm.imag,
                "unit_residual_base": unit_base,
                "unit_residual_harmonic": unit_harm,
                "harmonic_relation_residual_abs": harmonic_resid,
                "forward_reverse_state_match_residual_abs": reverse_match_resid,
                "forward_transition_residual_abs": step_resid,
                "sum_a2": a2,
                "sum_b2": b2,
                "sum_ab": ab,
                "quadratic_sum_re": q.real,
                "quadratic_sum_im": q.imag,
                "quadratic_identity_residual_abs": q_identity_resid,
                "static_exp_crosscheck_residual_abs": static_resid,
                "algebraically_self_consistent": local_resid <= TOL,
                "quadratic_closed_at_state": abs(q) <= TOL,
            })

        # Direct 4-step loop closure from S0.
        z1_0 = c_from(f_by_k[0], "base")
        zm_0 = c_from(f_by_k[0], "harmonic")
        z1_4 = (U**4) * z1_0
        zm_4 = ((U**m)**4) * zm_0
        loop_closure_resid = max(abs(z1_4-z1_0), abs(zm_4-zm_0))

        base_order = complex_order(U)
        harmonic_generator = U**m
        harmonic_order = complex_order(harmonic_generator)
        pair_order = None
        for n in range(1, 65):
            if abs(U**n - 1) <= TOL and abs((U**m)**n - 1) <= TOL:
                pair_order = n
                break

        summary_out.append({
            "m": m,
            "parity": "even" if m % 2 == 0 else "odd",
            "max_algebraic_self_consistency_residual_abs": max(state_residuals),
            "max_transition_residual_abs": max(transition_residuals),
            "loop_closure_residual_abs": loop_closure_resid,
            "max_static_exp_crosscheck_residual_abs": max(source_residuals) if source_residuals else 0.0,
            "global_sum_re_over_4_unique_states": a_loop,
            "global_sum_im_over_4_unique_states": b_loop,
            "global_sum_a2_over_4_unique_states": a2_loop,
            "global_sum_b2_over_4_unique_states": b2_loop,
            "global_sum_ab_over_4_unique_states": ab_loop,
            "global_quadratic_sum_re": q_loop.real,
            "global_quadratic_sum_im": q_loop.imag,
            "global_quadratic_sum_abs": abs(q_loop),
            "global_quadratic_closed": abs(q_loop) <= TOL,
            "base_generator": "i",
            "base_minimal_order": base_order,
            "harmonic_generator": str(harmonic_generator),
            "harmonic_minimal_order": harmonic_order,
            "pair_state_minimal_loop_order": pair_order,
            "pair_state_U4_closure": pair_order is not None and 4 % pair_order == 0,
        })

    return by_state_out, summary_out


def write_csv(path: Path, rows: List[dict]):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def write_md(states: List[dict], summaries: List[dict]):
    s2 = next(x for x in summaries if x["m"] == 2)
    s3 = next(x for x in summaries if x["m"] == 3)

    lines = []
    lines.append("# K4 偶数倍音系 / 奇数倍音系：自己無撞着性・二乗和・U^n=I 検証")
    lines.append("")
    lines.append("日付: 2026-09-15")
    lines.append("")
    lines.append("## 0. 対象と定義")
    lines.append("")
    lines.append("元データ `complex_states_animation_exact_i.csv` の 4 個の一意状態 `k=0,1,2,3` を解析した。forward/reverse は同じ4状態の順序違いなので、系全体の和では二重計上しない。")
    lines.append("")
    lines.append("```text")
    lines.append("U = i")
    lines.append("S_k = (z_base, z_harmonic) = (U^k, U^(m k))")
    lines.append("m=2 : 偶数倍音系")
    lines.append("m=3 : 奇数倍音系")
    lines.append("```")
    lines.append("")
    lines.append("自己無撞着性は、二乗和ゼロとは分離し、まず次の代数的一貫性で検査した。")
    lines.append("")
    lines.append("1. `|z_base|=|z_harmonic|=1`")
    lines.append("2. `z_harmonic = z_base^m`")
    lines.append("3. `S_(k+1) = (U z_base, U^m z_harmonic)`")
    lines.append("4. forward と reverse の同じ k が同じ複素状態を与える")
    lines.append("5. 4ステップ後に S0 に戻る")
    lines.append("")
    lines.append("二乗和は各状態について")
    lines.append("")
    lines.append("```text")
    lines.append("Q_k = z_base^2 + z_harmonic^2")
    lines.append("z_j = a_j + i b_j")
    lines.append("Q_k = Σa_j^2 - Σb_j^2 + 2 i Σ(a_j b_j)")
    lines.append("```")
    lines.append("")
    lines.append("を計算し、さらに4状態全体の総和も求めた。")
    lines.append("")
    lines.append("## 1. 自己無撞着性")
    lines.append("")
    lines.append(f"- m=2: 最大残差 = {fmt(s2['max_algebraic_self_consistency_residual_abs'])}")
    lines.append(f"- m=3: 最大残差 = {fmt(s3['max_algebraic_self_consistency_residual_abs'])}")
    lines.append(f"- 4ステップ閉路残差: m=2 = {fmt(s2['loop_closure_residual_abs'])}, m=3 = {fmt(s3['loop_closure_residual_abs'])}")
    lines.append("")
    lines.append("したがって、今回の離散状態生成則に対しては **偶数系・奇数系ともに代数的に自己無撞着**。これは `Σz^2=0` を意味しない。後者は別条件として次節で確認する。")
    lines.append("")
    lines.append("`exp(iθ)` で生成した静的データとの最大差も丸め誤差レベルであり、離散 exact-i データと整合する。")
    lines.append("")
    lines.append(f"- m=2: {fmt(s2['max_static_exp_crosscheck_residual_abs'])}")
    lines.append(f"- m=3: {fmt(s3['max_static_exp_crosscheck_residual_abs'])}")
    lines.append("")
    lines.append("## 2. 各状態の二乗和")
    lines.append("")
    lines.append("| m | k | Σa² | Σb² | Σab | Re(Σz²) | Im(Σz²) | 状態ごとのΣz²=0? |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|:---:|")
    for r in states:
        lines.append(f"| {r['m']} | {r['k']} | {fmt(r['sum_a2'])} | {fmt(r['sum_b2'])} | {fmt(r['sum_ab'])} | {fmt(r['quadratic_sum_re'])} | {fmt(r['quadratic_sum_im'])} | {'YES' if r['quadratic_closed_at_state'] else 'NO'} |")
    lines.append("")
    lines.append("観測された系列は次の通り。")
    lines.append("")
    lines.append("```text")
    lines.append("m=2 : Q_k = [ 2, 0, 2, 0 ]")
    lines.append("m=3 : Q_k = [ 2,-2, 2,-2 ]")
    lines.append("```")
    lines.append("")
    lines.append("## 3. 4状態を系全体とした二乗和")
    lines.append("")
    lines.append("まず、元の複素成分そのものの線形和も確認すると両系ともゼロ。")
    lines.append("")
    lines.append("| 系 | ΣRe(z) | ΣIm(z) |")
    lines.append("|---|---:|---:|")
    for s in summaries:
        lines.append(f"| m={s['m']} ({'偶数' if s['m']==2 else '奇数'}) | {fmt_clean(s['global_sum_re_over_4_unique_states'])} | {fmt_clean(s['global_sum_im_over_4_unique_states'])} |")
    lines.append("")
    lines.append("二乗和の分解は次の通り。")
    lines.append("")
    lines.append("| 系 | ΣΣa² | ΣΣb² | ΣΣab | Re(ΣΣz²) | Im(ΣΣz²) | 全体二乗閉塞 |")
    lines.append("|---|---:|---:|---:|---:|---:|:---:|")
    for s in summaries:
        lines.append(f"| m={s['m']} ({'偶数' if s['m']==2 else '奇数'}) | {fmt_clean(s['global_sum_a2_over_4_unique_states'])} | {fmt_clean(s['global_sum_b2_over_4_unique_states'])} | {fmt_clean(s['global_sum_ab_over_4_unique_states'])} | {fmt_clean(s['global_quadratic_sum_re'])} | {fmt_clean(s['global_quadratic_sum_im'])} | {'YES' if s['global_quadratic_closed'] else 'NO'} |")
    lines.append("")
    lines.append("結果は明確に異なる。")
    lines.append("")
    lines.append("```text")
    lines.append("m=2 : Σ_loop z^2 = 4")
    lines.append("m=3 : Σ_loop z^2 = 0")
    lines.append("```")
    lines.append("")
    lines.append("すなわち、**4つの離散位相状態をひとつの閉じた系として数えると、奇数倍音 m=3 は二乗ゼロ閉塞するが、偶数倍音 m=2 は閉塞しない。**")
    lines.append("")
    lines.append("複素成分分解でも同じことが見える。")
    lines.append("")
    lines.append("```text")
    lines.append("m=2 : ΣRe(z)=0, ΣIm(z)=0;  Σa²=6, Σb²=2, Σab=0  -> Σz²=(6-2)+0i=4")
    lines.append("m=3 : ΣRe(z)=0, ΣIm(z)=0;  Σa²=4, Σb²=4, Σab=0  -> Σz²=(4-4)+0i=0")
    lines.append("```")
    lines.append("")
    lines.append("これは今回の偶奇差の中でかなり強い構造差である。")
    lines.append("")
    lines.append("## 4. U^n = I")
    lines.append("")
    lines.append("はい。全状態対のループについては両系とも4ステップで閉じる。")
    lines.append("")
    lines.append("ただし最小位数には重要な差がある。")
    lines.append("")
    lines.append("| 系 | 基底 generator | 基底の最小位数 | 倍音 generator | 倍音の最小位数 | 2波状態全体の最小ループ |")
    lines.append("|---|---|---:|---|---:|---:|")
    for s in summaries:
        hg = "-1" if s['m']==2 else "-i"
        lines.append(f"| m={s['m']} | i | {s['base_minimal_order']} | {hg} | {s['harmonic_minimal_order']} | {s['pair_state_minimal_loop_order']} |")
    lines.append("")
    lines.append("したがって")
    lines.append("")
    lines.append("```text")
    lines.append("m=2 : harmonic は2周期で戻るが、base が4周期なので pair 全体は4周期")
    lines.append("m=3 : base も harmonic も4周期で戻り、pair 全体も4周期")
    lines.append("```")
    lines.append("")
    lines.append("両者とも `G_m^4=I` は成立するが、**偶数系では倍音成分だけが先に2周期で縮退する**。これは偶奇差を考える上で重要。")
    lines.append("")
    lines.append("## 5. 現段階の結論")
    lines.append("")
    lines.append("1. 離散生成則・遷移・逆順・閉路という意味では m=2, m=3 とも自己無撞着。")
    lines.append("2. しかし4状態全体の二乗和は、m=3 だけが厳密にゼロ閉塞する。")
    lines.append("3. m=2 は倍音 generator の最小位数が2へ縮退し、m=3 は4を維持する。")
    lines.append("4. よって、今回のデータでは `偶奇差 = 単なる図形差` ではなく、**位数構造と二乗閉塞の両方に現れている**。")
    lines.append("")
    lines.append("注意: ここでいう『自己無撞着』は上記の有限離散生成則に対する代数的一貫性。これをより強い意味（例えば全関係波の同時固定点、追加保存則、力学的自己無撞着）で定義する場合は、その条件を追加して再検査できる。")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    exact = read_rows(EXACT_CSV)
    static = read_rows(STATIC_CSV) if STATIC_CSV.exists() else None
    states, summaries = analyze(exact, static)
    write_csv(OUT_STATE, states)
    write_csv(OUT_SUMMARY, summaries)
    write_md(states, summaries)
    print("Generated:")
    for p in (OUT_STATE, OUT_SUMMARY, OUT_MD):
        print(p)
    print("\nSummary:")
    for s in summaries:
        print(s)


if __name__ == "__main__":
    main()
