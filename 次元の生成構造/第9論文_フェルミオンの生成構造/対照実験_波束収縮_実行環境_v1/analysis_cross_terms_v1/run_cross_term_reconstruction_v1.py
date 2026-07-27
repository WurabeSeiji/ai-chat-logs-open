#!/usr/bin/env python3
"""実行順序⑥⑦: 交差相関行列の再構成と局在性との時系列比較。

入力: production_dump_v1 の npz（倍音別複素係数）＋ 同ラン出力の rows CSV（L, N_eff）
出力: series/*.csv       (case,R)ごとの全衝突時系列
      cross_term_summary_v1.csv  (case,R)ごとの集約量
      analysis_meta_v1.json      解析条件の記録

定義（調査指示 §5 準拠。モード基底は符号付き chi 周波数、eta は Frobenius 縮約）:
  C_mn(t)   = Σ_η c_m(η,t) c_n(η,t)^*          交差相関行列（対角=モード強度）
  X(t)      = Σ_{m≠n} |C_mn|^2                  交差項総量（状態は単位ノルムなので既に規格化済）
  Φ(t)      = |Σ_{m≠n} C_mn| / Σ_{m≠n} |C_mn|  位相整列度 ∈ [0,1]
  G_d(t)    = Σ_{m−n=d} C_mn                    差周波数別集約
  H_acq(t)  = Σ_{n ∉ supp(自初期状態)} p_n(t)    獲得倍音量（相手側から移った強度）

時系列比較（§6）:
  時間遅れ相関 ρ(τ) = corr(x(t), L(t+τ))。x=X または Φ。原系列と一階差分の両方。
  交差項が収縮の原因なら最大相関は τ≥0 側（交差項が先行）に出るはず。
  説明力分離回帰: L ~ H / L ~ H+X+Φ / L ~ X+Φ の R² 比較（§6.3）。
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ENV = HERE.parent
PROD = ENV / "production_dump_v1"
SERIES_DIR = HERE / "series"
SERIES_DIR.mkdir(exist_ok=True)

MAX_LAG = 16
SUPPORT_EPS = 1.0e-12

RUNS = [
    "01_femtofocus_R137_B12",
    "02_B12_keyR",
    "03_oddN_B1_keyR",
    "03_oddN_B2_keyR_fullM",   # B2 は全スペクトル版を正とする（部分版 03_oddN_B2_keyR は解析対象外）
    "03_oddN_B3_keyR",
    "03_oddN_B5_keyR",
    "03_oddN_B15_keyR",
    "03_oddN_B63_keyR",
]


def load_rows_lookup(output_dir: Path) -> dict:
    """rows CSV → (case_id, R_input, collision, channel) → (L, N_eff)"""
    rows_csv = next(output_dir.glob("*_rows_v1.csv"))
    table = {}
    with open(rows_csv) as f:
        for row in csv.DictReader(f):
            key = (row["case_id"], float(row["R_input"]), int(row["collision"]), row["channel"])
            table[key] = (float(row["L"]), float(row["N_eff"]))
    return table


def gram(c: np.ndarray) -> np.ndarray:
    return c @ c.conj().T


def cross_metrics(C: np.ndarray) -> tuple[float, float]:
    off = C - np.diag(np.diag(C))
    x = float(np.sum(np.abs(off) ** 2))
    denom = float(np.sum(np.abs(off)))
    phi = float(abs(off.sum()) / denom) if denom > 0.0 else float("nan")
    return x, phi


def g_d_values(C: np.ndarray, harms: np.ndarray, d_list: list[int]) -> list[float]:
    out = []
    idx = {int(n): i for i, n in enumerate(harms)}
    for d in d_list:
        total = 0.0 + 0.0j
        for n in harms:
            m = int(n) + d
            if m in idx:
                total += C[idx[m], idx[int(n)]]
        out.append(abs(total))
    return out


def lag_correlation(x: np.ndarray, y: np.ndarray, max_lag: int) -> tuple[list[float], int, float]:
    """ρ(τ)=corr(x(t), y(t+τ)), τ=-max_lag..max_lag。返り値: (系列, 最良τ, 最良値)"""
    vals = []
    for tau in range(-max_lag, max_lag + 1):
        if tau >= 0:
            a, b = x[: len(x) - tau], y[tau:]
        else:
            a, b = x[-tau:], y[: len(y) + tau]
        if len(a) < 8 or np.std(a) < 1e-15 or np.std(b) < 1e-15:
            vals.append(float("nan"))
            continue
        vals.append(float(np.corrcoef(a, b)[0, 1]))
    arr = np.array(vals)
    if np.all(np.isnan(arr)):
        return vals, 0, float("nan")
    best = int(np.nanargmax(np.abs(arr)))
    return vals, best - max_lag, float(arr[best])


def ols_r2(y: np.ndarray, X_cols: list[np.ndarray]) -> float:
    X = np.column_stack([np.ones_like(y)] + X_cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) @ (y - y.mean())))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def analyze_npz(npz_path: Path, rows_table: dict, run_name: str) -> dict:
    z = np.load(npz_path)
    meta = json.loads(str(z["meta"]))
    coeffs, colls, harms = z["coeffs"], z["collisions"], z["harmonics"]
    case_id, r_input = meta["case_id"], meta["R_input"]
    n_rec = len(colls)

    p0 = {ch: np.sum(np.abs(coeffs[0, i]) ** 2, axis=1) for i, ch in enumerate("AB")}
    support = {ch: p0[ch] > SUPPORT_EPS for ch in "AB"}

    series = {k: np.zeros(n_rec) for k in
              ["X_A", "X_B", "Phi_A", "Phi_B", "H_acq_A", "H_acq_B",
               "L_A", "L_B", "N_eff_A", "N_eff_B",
               "Gd1_A", "Gd2_A", "Gd3_A", "Gd4_A"]}

    for k in range(n_rec):
        coll = int(colls[k])
        for i, ch in enumerate("AB"):
            c = coeffs[k, i]
            C = gram(c)
            x, phi = cross_metrics(C)
            series[f"X_{ch}"][k] = x
            series[f"Phi_{ch}"][k] = phi
            p = np.real(np.diag(C))
            series[f"H_acq_{ch}"][k] = float(p[~support[ch]].sum())
            L, n_eff = rows_table[(case_id, r_input, coll, f"{ch}_channel")]
            series[f"L_{ch}"][k] = L
            series[f"N_eff_{ch}"][k] = n_eff
            if ch == "A":
                for d, val in zip((1, 2, 3, 4), g_d_values(C, harms, [1, 2, 3, 4])):
                    series[f"Gd{d}_A"][k] = val

    # --- 時系列 CSV 保存 ---
    stem = npz_path.stem.replace("harmonic_coeffs_", "").replace("_v1", "")
    out_csv = SERIES_DIR / f"series_{run_name}_{stem}_v1.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["collision"] + list(series.keys()))
        for k in range(n_rec):
            writer.writerow([int(colls[k])] + [repr(float(series[key][k])) for key in series])

    # --- ⑦ 時間遅れ相関と回帰（A チャネル = 収縮側） ---
    LA, XA, PhiA, HA = series["L_A"], series["X_A"], series["Phi_A"], series["H_acq_A"]
    _, tau_XL_raw, rho_XL_raw = lag_correlation(XA, LA, MAX_LAG)
    _, tau_PL_raw, rho_PL_raw = lag_correlation(np.nan_to_num(PhiA), LA, MAX_LAG)
    dX, dL = np.diff(XA), np.diff(LA)
    dPhi = np.diff(np.nan_to_num(PhiA))
    _, tau_XL_diff, rho_XL_diff = lag_correlation(dX, dL, MAX_LAG)
    _, tau_PL_diff, rho_PL_diff = lag_correlation(dPhi, dL, MAX_LAG)

    phi_fit = np.nan_to_num(PhiA)
    r2_H = ols_r2(LA, [HA])
    r2_HXPhi = ols_r2(LA, [HA, XA, phi_fit])
    r2_XPhi = ols_r2(LA, [XA, phi_fit])

    return {
        "run": run_name, "case_id": case_id, "R_input": r_input,
        "series_csv": out_csv.name,
        "n_records": n_rec,
        "L_A_max": float(LA.max()), "L_A_final": float(LA[-1]),
        "X_A_max": float(XA.max()), "X_A_final": float(XA[-1]),
        "Phi_A_mean": float(np.nanmean(PhiA)),
        "H_acq_A_max": float(HA.max()),
        "corr_L_X_raw": float(np.corrcoef(XA, LA)[0, 1]) if np.std(XA) > 1e-15 and np.std(LA) > 1e-15 else float("nan"),
        "corr_L_H_raw": float(np.corrcoef(HA, LA)[0, 1]) if np.std(HA) > 1e-15 and np.std(LA) > 1e-15 else float("nan"),
        "tau_XL_raw": tau_XL_raw, "rho_XL_raw": rho_XL_raw,
        "tau_XL_diff": tau_XL_diff, "rho_XL_diff": rho_XL_diff,
        "tau_PhiL_raw": tau_PL_raw, "rho_PhiL_raw": rho_PL_raw,
        "tau_PhiL_diff": tau_PL_diff, "rho_PhiL_diff": rho_PL_diff,
        "R2_L_on_H": r2_H, "R2_L_on_HXPhi": r2_HXPhi, "R2_L_on_XPhi": r2_XPhi,
        "dR2_cross_terms": r2_HXPhi - r2_H,
    }


def main() -> None:
    results = []
    for run_name in RUNS:
        out_dir = PROD / run_name / "output"
        dump_dir = out_dir / "harmonic_dump_v1"
        if not dump_dir.exists():
            print(f"skip (no dump): {run_name}")
            continue
        rows_table = load_rows_lookup(out_dir)
        npzs = sorted(dump_dir.glob("*.npz"))
        print(f"=== {run_name}: {len(npzs)} npz ===", flush=True)
        for p in npzs:
            results.append(analyze_npz(p, rows_table, run_name))

    keys = list(results[0].keys())
    with open(HERE / "cross_term_summary_v1.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)

    (HERE / "analysis_meta_v1.json").write_text(json.dumps({
        "runs": RUNS, "max_lag": MAX_LAG, "support_eps": SUPPORT_EPS,
        "mode_basis": "signed chi frequency (carrier included), eta contracted via Frobenius inner product",
        "note_B2": "B2 は fullM(M=255) 版のみ解析。M=4 部分版は対象外。",
        "n_series": len(results),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"series: {len(results)}  summary: cross_term_summary_v1.csv", flush=True)


if __name__ == "__main__":
    main()
