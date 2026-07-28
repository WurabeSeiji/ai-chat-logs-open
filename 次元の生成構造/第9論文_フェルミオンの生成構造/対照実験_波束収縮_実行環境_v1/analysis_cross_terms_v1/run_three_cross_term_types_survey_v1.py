#!/usr/bin/env python3
"""三種交差項の全系列調査 v1

系統A（波束収縮系）のダンプ済み複素係数（production_dump_v1、184系列×全257衝突）から、
「交差項」と読める全種類を系統的に計測し、ケース間（ボゾン×ボゾン対照 B1 を含む）で比較する。

計測する交差項の種類:
  T1 チャネル内エルミート型   X_ch = Σ_{m≠n}|C_mn|²,  C_mn = Σ_η c_m(η) c_n(η)*
                               （⑥⑦で解析済みの型。参照用に再計算）
  T3 状態間エルミート型       Y_AB = ‖a b†‖_F,  (a b†)_mn = Σ_η a_m(η) b_n(η)*
                               （調査指針§3.2「状態間交差相関 C_AB=ab*」に対応）
  T4 毛間相関（チャネル内）    h12_ch = |<c̃(毛1), c̃(毛2)>_chi|
                               （η軸FFT後の毛モード1×毛モード2のchi縮約内積）
  T5 双線形・ηスライス別      C(η) = (A_η² − B2_η)/2,  A_η=Σ_m x_m(η), B2_η=Σ_m x_m(η)²
                               指標: ‖C(η)‖₂ と max_η|C(η)|（η和を取る前の生の双線形構造）
  T6 双線形・ペア別η和        max_{m≠n} |Σ_η x_m x_n| （構造的ゼロの検証、毛直交による）

構造検証:
  凸結合転送位相 v(k)（⑥の方法で取得）に対し、
    corr( Y_AB(k),  sqrt(v(1−v)) )   状態間交差項は転送の中間相か
    corr( h12_A(k), sqrt(v(1−v)) )   毛間相関は転送の中間相か
  を全非縮退系列で検証する。

出力:
  series_three_types/*.csv                        系列別時系列（全衝突）
  three_cross_term_types_summary_v1.csv           系列別集約
  three_cross_term_types_aggregate_v1.json        ケース別・検証項目の集約と判定
"""

from __future__ import annotations

import csv
import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROD = HERE.parent / "production_dump_v1"
OUT_SERIES = HERE / "series_three_types"
OUT_SERIES.mkdir(exist_ok=True)

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

HAIR_A, HAIR_B = 1, 2  # エンジン既定の毛モード（m_A=1, m_B=2、20260713 JSON params で確認済み）


def series_metrics(vals: np.ndarray) -> dict:
    m = float(vals.mean())
    return {
        "max": float(vals.max()), "min": float(vals.min()), "mean": m,
        "cv": float(vals.std() / m) if m > 1e-300 else float("nan"),
    }


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def analyze_npz(npz_path: Path, run_name: str) -> dict:
    z = np.load(npz_path)
    meta = json.loads(str(z["meta"]))
    coeffs, colls = z["coeffs"], z["collisions"]
    K = coeffs.shape[0]

    a0, b0 = coeffs[0, 0], coeffs[0, 1]
    CA0 = (a0 @ a0.conj().T).real
    CB0 = (b0 @ b0.conj().T).real
    basis = np.column_stack([CA0.ravel(), CB0.ravel()])

    cols = {k: np.zeros(K) for k in
            ["X_A", "X_B", "Y_AB", "h12_A", "h12_B",
             "Cslice_norm_A", "Cslice_max_A", "Cslice_norm_B", "Cslice_max_B",
             "pairwise_bilinear_max_A", "v"]}

    for k in range(K):
        a, b = coeffs[k, 0], coeffs[k, 1]
        # T1
        Ga = a @ a.conj().T
        Gb = b @ b.conj().T
        offa = Ga - np.diag(np.diag(Ga))
        offb = Gb - np.diag(np.diag(Gb))
        cols["X_A"][k] = float(np.sum(np.abs(offa) ** 2))
        cols["X_B"][k] = float(np.sum(np.abs(offb) ** 2))
        # T3
        Cab = a @ b.conj().T
        cols["Y_AB"][k] = float(np.sqrt(np.sum(np.abs(Cab) ** 2)))
        # T4
        ea = np.fft.fft(a, axis=1, norm="ortho")
        eb = np.fft.fft(b, axis=1, norm="ortho")
        cols["h12_A"][k] = float(np.abs(np.vdot(ea[:, HAIR_A], ea[:, HAIR_B])))
        cols["h12_B"][k] = float(np.abs(np.vdot(eb[:, HAIR_A], eb[:, HAIR_B])))
        # T5（両チャネル）
        for ch, mat, prefix in ((0, a, "A"), (1, b, "B")):
            As = mat.sum(axis=0)
            B2s = (mat * mat).sum(axis=0)
            Cs = (As * As - B2s) / 2.0
            cols[f"Cslice_norm_{prefix}"][k] = float(np.sqrt(np.sum(np.abs(Cs) ** 2)))
            cols[f"Cslice_max_{prefix}"][k] = float(np.abs(Cs).max())
        # T6
        P = a @ a.T
        cols["pairwise_bilinear_max_A"][k] = float(np.abs(P - np.diag(np.diag(P))).max())
        # 転送位相 v（凸結合フィット、⑥と同一手法）
        sol, *_ = np.linalg.lstsq(basis, Ga.real.ravel(), rcond=None)
        cols["v"][k] = float(sol[1])

    # 系列別時系列 CSV
    stem = npz_path.stem.replace("harmonic_coeffs_", "").replace("_v1", "")
    out_csv = OUT_SERIES / f"three_types_{run_name}_{stem}_v1.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["collision"] + list(cols.keys()))
        for k in range(K):
            w.writerow([int(colls[k])] + [repr(float(cols[key][k])) for key in cols])

    v = cols["v"]
    bridge = np.sqrt(np.clip(v * (1.0 - v), 0.0, None))
    row = {
        "run": run_name, "case_id": meta["case_id"], "R_input": meta["R_input"],
        "series_csv": out_csv.name,
    }
    for key in ["X_A", "Y_AB", "h12_A", "h12_B", "Cslice_norm_A", "Cslice_max_A",
                "Cslice_norm_B", "Cslice_max_B", "pairwise_bilinear_max_A"]:
        for stat, val in series_metrics(cols[key]).items():
            row[f"{key}_{stat}"] = val
    row["corr_YAB_bridge"] = safe_corr(cols["Y_AB"], bridge)
    row["corr_h12A_bridge"] = safe_corr(cols["h12_A"], bridge)
    row["corr_h12B_bridge"] = safe_corr(cols["h12_B"], bridge)
    row["corr_XA_v"] = safe_corr(cols["X_A"], v)
    row["v_min"] = float(v.min())
    row["v_max"] = float(v.max())
    return row


def main() -> None:
    rows = []
    for run_name in RUNS:
        npzs = sorted((PROD / run_name / "output" / "harmonic_dump_v1").glob("*.npz"))
        print(f"=== {run_name}: {len(npzs)} npz ===", flush=True)
        for p in npzs:
            rows.append(analyze_npz(p, run_name))

    keys = list(rows[0].keys())
    with open(HERE / "three_cross_term_types_summary_v1.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})

    # ケース別集約（非縮退系列のみで橋検証を判定）
    nond = [r for r in rows if np.isfinite(r["corr_YAB_bridge"])]
    agg = {
        "n_series": len(rows),
        "n_nondegenerate": len(nond),
        "T6_structural_zero_worst": max(r["pairwise_bilinear_max_A_max"] for r in rows),
        "bridge_check_YAB": {
            "min_corr": min(r["corr_YAB_bridge"] for r in nond),
            "median_corr": float(np.median([r["corr_YAB_bridge"] for r in nond])),
        },
        "bridge_check_h12A": {
            "min_corr": min(r["corr_h12A_bridge"] for r in nond),
            "median_corr": float(np.median([r["corr_h12A_bridge"] for r in nond])),
        },
        "by_case_at_R137": {},
    }
    for r in rows:
        if abs(r["R_input"] - 0.6971778791282474) < 1e-12:
            agg["by_case_at_R137"][f"{r['run']}"] = {
                "X_A_max": r["X_A_max"],
                "Y_AB_max": r["Y_AB_max"],
                "h12_A_max": r["h12_A_max"],
                "Cslice_norm_A_max": r["Cslice_norm_A_max"],
                "Cslice_norm_B_max": r["Cslice_norm_B_max"],
                "corr_YAB_bridge": r["corr_YAB_bridge"],
                "corr_h12A_bridge": r["corr_h12A_bridge"],
            }
    (HERE / "three_cross_term_types_aggregate_v1.json").write_text(
        json.dumps(agg, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps(agg, ensure_ascii=False, indent=2, default=str), flush=True)


if __name__ == "__main__":
    main()
