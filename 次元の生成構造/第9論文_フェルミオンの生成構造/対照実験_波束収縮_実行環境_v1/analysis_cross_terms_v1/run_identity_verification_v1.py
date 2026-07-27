#!/usr/bin/env python3
"""恒等式検証 v1: 報告書 cross_term_analysis_report_v1.md §1 の数値事実を
全 npz に対して系統的に検証し、恒久記録として保存する。

検証項目（各 (case, R) 系列、全記録衝突）:
  V1 実数性        max|Im C^A_mn(k)| / max|C^A_mn(k)|
  V2 η直交性       max_m |<a0_m, b0_m>_η|（A,B初期係数行の η 内積）
  V3 凸結合恒等式  C^A(k) = u(k)C^A(0) + v(k)C^B(0) の最大残差、max|u+v−1|、v の範囲
  V4 転写周期      v(k) の極値間隔から推定した周期（衝突数）
  V5 L アフィン恒等 R²(L~H), R²(L~H,X)（10桁精度）
  V6 Φ 時間不変性  Φ_A(k) の max−min
  V7 X 等分配天井  X_A_max と 1−1/N_modes の差（N_modes = 初期A∪B支持モード数）

出力: identity_verification_v1.csv（系列別全指標）
      identity_verification_summary_v1.json（全系列の最悪値集約と判定）

注意: 縮退対照点（R=0.0 / 1.0 など交換が起きない点）では v(k) が定数となり
      V4/V5 は未定義（NaN）になる。判定は非縮退系列のみで行う。
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROD = HERE.parent / "production_dump_v1"

RUNS = [
    "01_femtofocus_R137_B12",
    "02_B12_keyR",
    "03_oddN_B1_keyR",
    "03_oddN_B2_keyR_fullM",
    "03_oddN_B3_keyR",
    "03_oddN_B5_keyR",
    "03_oddN_B15_keyR",
    "03_oddN_B63_keyR",
]

SUPPORT_EPS = 1.0e-12


def gram(c: np.ndarray) -> np.ndarray:
    return c @ c.conj().T


def load_rows_lookup(output_dir: Path) -> dict:
    rows_csv = next(output_dir.glob("*_rows_v1.csv"))
    table = {}
    with open(rows_csv) as f:
        for row in csv.DictReader(f):
            key = (row["case_id"], float(row["R_input"]), int(row["collision"]), row["channel"])
            table[key] = float(row["L"])
    return table


def ols_r2(y: np.ndarray, cols: list[np.ndarray]) -> float:
    X = np.column_stack([np.ones_like(y)] + cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_tot = float((y - y.mean()) @ (y - y.mean()))
    return 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else float("nan")


def estimate_period(v: np.ndarray) -> float:
    """v(k) の内部極値数から周期を推定。極値がなければ NaN（縮退系列）。"""
    if np.ptp(v) < 1e-12:
        return float("nan")
    dv = np.diff(v)
    extrema = int(np.sum(np.diff(np.sign(dv)) != 0))
    if extrema == 0:
        return float("nan")
    return 2.0 * (len(v) - 1) / extrema


def verify_npz(npz_path: Path, rows_table: dict, run_name: str) -> dict:
    z = np.load(npz_path)
    meta = json.loads(str(z["meta"]))
    coeffs, colls = z["coeffs"], z["collisions"]
    case_id, r_input = meta["case_id"], meta["R_input"]
    n_rec = len(colls)

    a0, b0 = coeffs[0, 0], coeffs[0, 1]
    CA0, CB0 = gram(a0), gram(b0)

    # V2 η直交性
    v2_eta = float(np.abs(np.einsum("me,me->m", a0.conj(), b0)).max())

    # V7 用: 初期 A∪B 支持モード数
    p_union = np.sum(np.abs(a0) ** 2, axis=1) + np.sum(np.abs(b0) ** 2, axis=1)
    n_modes = int(np.sum(p_union > SUPPORT_EPS))

    basis = np.column_stack([CA0.real.ravel(), CB0.real.ravel()])
    v1_im, v1_abs = 0.0, 0.0
    v3_resid, v3_uv = 0.0, 0.0
    v_series = np.zeros(n_rec)
    x_series = np.zeros(n_rec)
    phi_series = np.zeros(n_rec)
    h_series = np.zeros(n_rec)
    L_series = np.zeros(n_rec)
    supp_a0 = np.sum(np.abs(a0) ** 2, axis=1) > SUPPORT_EPS

    for k in range(n_rec):
        C = gram(coeffs[k, 0])
        v1_im = max(v1_im, float(np.abs(C.imag).max()))
        v1_abs = max(v1_abs, float(np.abs(C).max()))
        sol, *_ = np.linalg.lstsq(basis, C.real.ravel(), rcond=None)
        v3_resid = max(v3_resid, float(np.abs(basis @ sol - C.real.ravel()).max()))
        v3_uv = max(v3_uv, abs(float(sol.sum()) - 1.0))
        v_series[k] = float(sol[1])
        off = C - np.diag(np.diag(C))
        x_series[k] = float(np.sum(np.abs(off) ** 2))
        denom = float(np.sum(np.abs(off)))
        phi_series[k] = float(abs(off.sum()) / denom) if denom > 0 else np.nan
        p = np.real(np.diag(C))
        h_series[k] = float(p[~supp_a0].sum())
        L_series[k] = rows_table[(case_id, r_input, int(colls[k]), "A_channel")]

    # L の分散が機械ノイズ水準の系列（例: B1 の R=0 全交換スワップ）では回帰は無意味
    l_degenerate = float(np.ptp(L_series)) < 1e-10 * max(abs(float(L_series.max())), 1e-300)
    if l_degenerate:
        r2_h = r2_hx = float("nan")
    else:
        r2_h = ols_r2(L_series, [h_series])
        r2_hx = ols_r2(L_series, [h_series, x_series])
    phi_span = float(np.nanmax(phi_series) - np.nanmin(phi_series)) if not np.all(np.isnan(phi_series)) else float("nan")
    x_ceiling_gap = float((1.0 - 1.0 / n_modes) - x_series.max())

    return {
        "run": run_name, "case_id": case_id, "R_input": r_input, "n_records": n_rec,
        "V1_max_im_ratio": v1_im / v1_abs if v1_abs > 0 else float("nan"),
        "V2_eta_overlap": v2_eta,
        "V3_convex_resid": v3_resid,
        "V3_uv_minus_1": v3_uv,
        "V3_v_min": float(v_series.min()), "V3_v_max": float(v_series.max()),
        "V4_period": estimate_period(v_series),
        "V5_R2_L_on_H": r2_h, "V5_R2_L_on_HX": r2_hx,
        "V5_affine_defect": 1.0 - r2_hx if np.isfinite(r2_hx) else float("nan"),
        "V6_phi_span": phi_span,
        "V7_n_modes": n_modes, "V7_x_ceiling_gap": x_ceiling_gap,
    }


def main() -> None:
    results = []
    for run_name in RUNS:
        out_dir = PROD / run_name / "output"
        rows_table = load_rows_lookup(out_dir)
        npzs = sorted((out_dir / "harmonic_dump_v1").glob("*.npz"))
        print(f"=== {run_name}: {len(npzs)} npz ===", flush=True)
        for p in npzs:
            results.append(verify_npz(p, rows_table, run_name))

    with open(HERE / "identity_verification_v1.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        for r in results:
            writer.writerow({k: (repr(float(v)) if isinstance(v, float) else v) for k, v in r.items()})

    def worst(key, agg=max, finite_only=True):
        vals = [r[key] for r in results]
        vals = [v for v in vals if np.isfinite(v)] if finite_only else vals
        return agg(vals) if vals else float("nan")

    nondegenerate = [r for r in results if np.isfinite(r["V4_period"])]
    affine_defects = [r["V5_affine_defect"] for r in nondegenerate if np.isfinite(r["V5_affine_defect"])]
    # Φ 時間不変性はラン系列ごとに判定（奇数N系列に固有の性質のため）
    phi_by_run: dict = {}
    for r in results:
        span = r["V6_phi_span"]
        if np.isfinite(span):
            phi_by_run.setdefault(r["run"], []).append(span)
    phi_invariant_runs = sorted(run for run, spans in phi_by_run.items() if max(spans) < 1e-6)
    phi_varying_runs = sorted(run for run, spans in phi_by_run.items() if max(spans) >= 1e-6)
    summary = {
        "n_series": len(results),
        "n_nondegenerate": len(nondegenerate),
        "V1_worst_im_ratio": worst("V1_max_im_ratio"),
        "V2_worst_eta_overlap": worst("V2_eta_overlap"),
        "V3_worst_convex_resid": worst("V3_convex_resid"),
        "V3_worst_uv_minus_1": worst("V3_uv_minus_1"),
        "V4_period_range_nondegenerate": [min(r["V4_period"] for r in nondegenerate), max(r["V4_period"] for r in nondegenerate)] if nondegenerate else None,
        "V5_worst_affine_defect": max(affine_defects) if affine_defects else None,
        "V5_n_L_degenerate_excluded": sum(1 for r in results if not np.isfinite(r["V5_affine_defect"])),
        "V6_phi_invariant_runs": phi_invariant_runs,
        "V6_phi_varying_runs": phi_varying_runs,
        "V7_note": "X 天井 1-1/n は初期網が union 支持上で等重みの場合（奇数N系列）のみ厳密。他は診断値。",
        "V7_worst_x_ceiling_gap_nondegenerate": max(abs(r["V7_x_ceiling_gap"]) for r in nondegenerate) if nondegenerate else None,
        "verdict": {
            "C_real_all_series": worst("V1_max_im_ratio") < 1e-12,
            "eta_orthogonal_all_series": worst("V2_eta_overlap") < 1e-12,
            "convex_identity_all_series": worst("V3_convex_resid") < 1e-12 and worst("V3_uv_minus_1") < 1e-6,
            "L_affine_in_H_X_nondegenerate": (max(affine_defects) < 1e-8) if affine_defects else None,
            "phi_time_invariant_odd_series_only": all(("oddN_B1" in r or "oddN_B3" in r or "oddN_B5" in r or "oddN_B15" in r or "oddN_B63" in r) for r in phi_invariant_runs) and len(phi_invariant_runs) > 0,
        },
    }
    (HERE / "identity_verification_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str), flush=True)


if __name__ == "__main__":
    main()
