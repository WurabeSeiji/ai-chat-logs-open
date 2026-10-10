"""厳密根 R_{124,23} = cos^2(23π/124) における二波交換の完全回帰を、0713 の元プログラムの更新式で記録する。

目的
  第十三思考実験 §6.1 は R = 0.70 で「ω/2π が無理数なので完全回帰は起きない」と記した。
  R = 0.70 は、0713 の反射率掃引で二波の局在性の差が最小になった値である（[S3] §7）。その後の 0.7 前後の
  細かい掃引（0715 稿 [S5]、0718 稿 [S4] §3・§4.1）で、ピークの中心が有限位数根
  R_{124,23} = cos^2(23π/124) = 0.697177927556659… に収束した。[S4] §6.5 のとおり
      λ_a^31 = i,  λ_a^62 = −1,  λ_a^93 = −i,  λ_a^124 = 1
  となる。本スクリプトは R_VALUE をこの厳密根に置き、0713 の元プログラムの関数だけを使って
  衝突 0..N_COLL（二完全周期 = 248）の全ステップを記録する。

方法（検証_R070交換の正規化と振幅_20261009/measure_normalization_trace.py と同じ骨格）
  original_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
  （0713 の原本と一字一句同一、SHA256 = f815320f…）をモジュールとして読み込み、
  原本 recursive_snapshot_states（920〜983 行）と同じ更新
      a_next = normalize(r * a + t * b)
      b_next = normalize(t * a + r * b)
  を行い、記録用の行だけを足す。評価量（N_eff, L）は原本の recursive_state_metrics で読む。

記録する量（毎衝突）
  - 交換位相 kω（閉形式、deg）、閉形式の移行率 sin^2(kω/2)
  - N_eff_A, N_eff_B, L_A, L_B（原本の関数）、閉形式 N_A = 1 + 31 sin^2(kω/2)
  - ピークの絶対高さ（総和で割らない密度の最大値）
  - 正規化前のノルム^2（normalize が取り除く量）、A・B の重なり
  - 複素状態そのものの初期状態からの距離 ε（第十三思考実験 §12.2）:
        eps_A = ‖a_k − a_0‖,  eps_B = ‖b_k − b_0‖,  eps_state = sqrt((eps_A^2 + eps_B^2)/2)
        eps_state_gp: 共通の大域位相を除いた同じ量
        eps_swap: 入れ替わり状態 (b_0, a_0) からの距離（λ_a^62 = −1 で 0 になるべき量）
  - 座標 (alpha, beta, gamma, delta): a_k = alpha a_0 + beta b_0, b_k = gamma a_0 + delta b_0
        初期二状態は直交しているので射影で求まり、復元残差 recon_residual で完全性を確かめる。
        eta 方向の識別振動が異なる（m_A=1, m_B=2）ため chi 密度の干渉項は厳密に 0 で、
        rho_A,k(chi) = |alpha|^2 rho_a0(chi) + |beta|^2 rho_b0(chi) が成り立つ（density_recon_max_error で確認）。

出力（results/）
  trace_R12423.csv               全ステップの上記の量（原本の normalize 付き更新）
  trace_pure_unitary_R12423.csv  normalize を外した純粋 U_R の同じ記録（対照）
  coefficients_R12423.csv        全ステップの (alpha, beta, gamma, delta) と復元残差 — 状態の完全記録（小さい）
  basis_density_R12423.csv       chi, rho_a0, rho_b0（総和で割らない密度）。coefficients と合わせて全ステップの波形を復元できる
  rho_all_steps_R12423.npz       全ステップの chi 密度（A, B）そのもの
  states_R12423.npz              全ステップの複素状態そのもの（8192 成分 × 2 チャネル × (N_COLL+1)）
  control_R070_vs_S7.json        対照: 同じループを R=0.70 で回し、[S7] の summary.json の値と一致するか
  summary.json                   要点
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE / "original_copy" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

spec = importlib.util.spec_from_file_location("orig", SRC)
orig = importlib.util.module_from_spec(spec)
sys.modules["orig"] = orig
spec.loader.exec_module(orig)

N_ROOT, M_ROOT = 124, 23
R_VALUE = math.cos(math.pi * M_ROOT / N_ROOT) ** 2      # R_{124,23}
N_COLL = 2 * N_ROOT                                      # 二完全周期まで記録（完全回帰は k = 124, 248）
N_A_HARM, N_B_HARM = 1, 63                               # 原本の図 1 と同じ（A: 次数 1、B: 奇数倍音 1..63）
KEY_STEPS = [0, 31, 42, 62, 93, 103, 124, 155, 186, 217, 248]
CONTROL_R = 0.70
CONTROL_COLLISIONS = [0, 1, 2, 3, 5, 10, 20, 42]
REFERENCE_S7 = OUT / "reference_summary_R070_S7.json"   # 検証_R070交換の正規化と振幅_20261009/results/summary.json の写し


def raw_density(params, v):
    """総和で割らない chi 密度 sum_eta |psi|^2（状態が正規化されていれば総和 1）。"""
    psi = orig.reshape_state(params, v)
    return np.sum(np.abs(psi) ** 2, axis=1)


def omega_of(r_value: float) -> float:
    """交換位相の一回あたりの進み: U_R の固有値 1 と -e^{i delta} の位相差。"""
    return math.pi + 2.0 * math.asin(math.sqrt(r_value))


def run_trace(r_value: float, n_coll: int, renormalize: bool, keep_states: bool):
    params = orig.Params()
    delta_f = orig.delta_from_reflection_rate(r_value)
    t, r, T, R = orig.scattering_coefficients(delta_f)
    hair_enabled = True
    a = orig.make_state(params, N_A_HARM, params.q_A, params.m_A, hair_enabled, params.A_A)
    b = orig.make_state(params, N_B_HARM, params.q_B, params.m_B, hair_enabled, params.A_B)
    a0, b0 = a.copy(), b.copy()
    rho_a0, rho_b0 = raw_density(params, a0), raw_density(params, b0)
    omega = omega_of(r_value)
    label = f"R{int(round(r_value * 100)):03d}"
    rows, coeffs = [], []
    rho_a_all = np.empty((n_coll + 1, params.chi_grid_n))
    rho_b_all = np.empty((n_coll + 1, params.chi_grid_n))
    states_a = np.empty((n_coll + 1, a.size), dtype=complex) if keep_states else None
    states_b = np.empty((n_coll + 1, b.size), dtype=complex) if keep_states else None
    density_recon_err = 0.0
    for k in range(n_coll + 1):
        ma = orig.recursive_state_metrics(params, "trace", label, delta_f, T, R, N_A_HARM, N_B_HARM,
                                          hair_enabled, k, "A_channel", a, a0, b0)
        mb = orig.recursive_state_metrics(params, "trace", label, delta_f, T, R, N_A_HARM, N_B_HARM,
                                          hair_enabled, k, "B_channel", b, a0, b0)
        ov = orig.inner(a, b)
        na_pre = r * a + t * b
        nb_pre = t * a + r * b
        # 座標（初期二状態への射影）と復元残差
        alpha, beta = orig.inner(a0, a), orig.inner(b0, a)
        gamma, delta = orig.inner(a0, b), orig.inner(b0, b)
        recon = max(float(np.linalg.norm(a - (alpha * a0 + beta * b0))),
                    float(np.linalg.norm(b - (gamma * a0 + delta * b0))))
        # 状態距離
        da, db = float(np.linalg.norm(a - a0)), float(np.linalg.norm(b - b0))
        phi = math.atan2((orig.inner(a0, a) + orig.inner(b0, b)).imag, (orig.inner(a0, a) + orig.inner(b0, b)).real)
        gp = complex(math.cos(phi), math.sin(phi))
        da_gp, db_gp = float(np.linalg.norm(a - gp * a0)), float(np.linalg.norm(b - gp * b0))
        swap = math.sqrt((float(np.linalg.norm(a - b0)) ** 2 + float(np.linalg.norm(b - a0)) ** 2) / 2.0)
        rho_a, rho_b = raw_density(params, a), raw_density(params, b)
        rho_a_all[k], rho_b_all[k] = rho_a, rho_b
        density_recon_err = max(density_recon_err,
                                float(np.max(np.abs(rho_a - (abs(alpha) ** 2 * rho_a0 + abs(beta) ** 2 * rho_b0)))),
                                float(np.max(np.abs(rho_b - (abs(gamma) ** 2 * rho_a0 + abs(delta) ** 2 * rho_b0)))))
        if keep_states:
            states_a[k], states_b[k] = a, b
        phase_deg = math.degrees(k * omega) % 360.0
        frac = math.sin(0.5 * k * omega) ** 2
        rows.append({
            "collision": k,
            "exchange_phase_deg": phase_deg,
            "closed_form_transfer_fraction": frac,
            "closed_form_N_eff_A": 1.0 + 31.0 * frac,
            "N_eff_A": ma["N_eff"], "N_eff_B": mb["N_eff"], "N_eff_sum": ma["N_eff"] + mb["N_eff"],
            "L_A": ma["L"], "L_B": mb["L"],
            "peak_raw_A": float(np.max(rho_a)), "peak_raw_B": float(np.max(rho_b)),
            "norm2_A": orig.norm2(a), "norm2_B": orig.norm2(b),
            "next_prenorm_norm2_A": orig.norm2(na_pre), "next_prenorm_norm2_B": orig.norm2(nb_pre),
            "overlap_abs": abs(ov),
            "weight_b0_in_A": abs(beta) ** 2, "weight_a0_in_B": abs(gamma) ** 2,
            "eps_A": da, "eps_B": db, "eps_state": math.sqrt((da * da + db * db) / 2.0),
            "eps_state_gp": math.sqrt((da_gp * da_gp + db_gp * db_gp) / 2.0),
            "eps_swap": swap,
            "recon_residual": recon,
        })
        coeffs.append({
            "collision": k,
            "alpha_re": alpha.real, "alpha_im": alpha.imag, "beta_re": beta.real, "beta_im": beta.imag,
            "gamma_re": gamma.real, "gamma_im": gamma.imag, "delta_re": delta.real, "delta_im": delta.imag,
            "recon_residual": recon,
        })
        if k >= n_coll:
            break
        if renormalize:
            a_next = orig.normalize(r * a + t * b)   # 原本 recursive_snapshot_states と同一
            b_next = orig.normalize(t * a + r * b)
        else:
            a_next = r * a + t * b
            b_next = t * a + r * b
        a, b = a_next, b_next
    chi, _ = orig.make_grids(params)
    return {
        "rows": rows, "coeffs": coeffs, "omega": omega,
        "coefficients": {"t": [t.real, t.imag], "r": [r.real, r.imag], "T": T, "R": R, "delta_f": delta_f},
        "chi": chi, "rho_a0": rho_a0, "rho_b0": rho_b0, "rho_a_all": rho_a_all, "rho_b_all": rho_b_all,
        "states_a": states_a, "states_b": states_b, "density_recon_max_error": density_recon_err,
    }


def write_csv(path: Path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    t0 = time.time()
    # 対照: R = 0.70 で同じループを回し、[S7] の記録値と比較する
    control = {"R": CONTROL_R, "reference": str(REFERENCE_S7.name), "available": REFERENCE_S7.exists()}
    ctrl = run_trace(CONTROL_R, max(CONTROL_COLLISIONS), True, False)
    control["omega_over_2pi"] = ctrl["omega"] / (2 * math.pi)
    if REFERENCE_S7.exists():
        ref = json.loads(REFERENCE_S7.read_text(encoding="utf-8"))
        diffs = {}
        for k in CONTROL_COLLISIONS:
            row = ctrl["rows"][k]
            diffs[str(k)] = {
                "peak_raw_A_diff": row["peak_raw_A"] - ref["peak_raw_A_at"][str(k)],
                "peak_raw_B_diff": row["peak_raw_B"] - ref["peak_raw_B_at"][str(k)],
            }
        control["peak_diffs_vs_S7"] = diffs
        control["max_abs_peak_diff_vs_S7"] = max(abs(v[d]) for v in diffs.values() for d in v)
        control["coefficients_match_S7"] = (abs(ctrl["coefficients"]["R"] - ref["coefficients"]["R"]) == 0.0
                                            and abs(ctrl["coefficients"]["T"] - ref["coefficients"]["T"]) == 0.0)
        control["omega_over_2pi_diff_vs_S7"] = control["omega_over_2pi"] - ref["closed_form_check"]["omega_over_2pi"]
    control["N_eff_A_at_42"] = ctrl["rows"][42]["N_eff_A"]
    control["N_eff_B_at_42"] = ctrl["rows"][42]["N_eff_B"]
    (OUT / "control_R070_vs_S7.json").write_text(json.dumps(control, ensure_ascii=False, indent=2), encoding="utf-8")

    # 本番: 厳密根 R_{124,23}
    main_run = run_trace(R_VALUE, N_COLL, True, True)
    pure_run = run_trace(R_VALUE, N_COLL, False, False)
    rows, pure_rows = main_run["rows"], pure_run["rows"]

    write_csv(OUT / "trace_R12423.csv", rows)
    write_csv(OUT / "trace_pure_unitary_R12423.csv", pure_rows)
    write_csv(OUT / "coefficients_R12423.csv", main_run["coeffs"])
    with open(OUT / "basis_density_R12423.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["chi", "chi_over_pi", "rho_a0", "rho_b0"])
        for x, ra, rb in zip(main_run["chi"], main_run["rho_a0"], main_run["rho_b0"]):
            w.writerow([repr(float(x)), repr(float(x) / math.pi), repr(float(ra)), repr(float(rb))])
    np.savez_compressed(OUT / "rho_all_steps_R12423.npz", chi=main_run["chi"],
                        rho_A=main_run["rho_a_all"], rho_B=main_run["rho_b_all"])
    np.savez_compressed(OUT / "states_R12423.npz", a=main_run["states_a"], b=main_run["states_b"],
                        chi_grid_n=orig.Params().chi_grid_n, eta_grid_n=orig.Params().eta_grid_n)

    eps = [x["eps_state"] for x in rows]
    first_window = [k for k in range(1, N_ROOT + 20) if k <= N_COLL]
    k_min_first = min(first_window, key=lambda k: eps[k])
    neff_diff = max(max(abs(p["N_eff_A"] - q["N_eff_A"]), abs(p["N_eff_B"] - q["N_eff_B"])) for p, q in zip(rows, pure_rows))
    summary = {
        "root": {"n": N_ROOT, "m": M_ROOT, "R": R_VALUE, "R_repr": repr(R_VALUE),
                 "N_of_R_4pi_over_1mR2": 4 * math.pi / (1 - R_VALUE) ** 2},
        "collisions": N_COLL, "coefficients": main_run["coefficients"],
        "omega_over_2pi": main_run["omega"] / (2 * math.pi), "omega_over_2pi_exact": (N_ROOT - M_ROOT) / N_ROOT,
        "control_R070": {k: control[k] for k in control if k != "peak_diffs_vs_S7"},
        "max_abs_overlap": max(x["overlap_abs"] for x in rows),
        "max_prenorm_channel_norm2_deviation_from_1": max(max(abs(x["next_prenorm_norm2_A"] - 1), abs(x["next_prenorm_norm2_B"] - 1)) for x in rows),
        "max_N_eff_difference_normalized_vs_pure": neff_diff,
        "max_abs_N_eff_A_minus_closed_form": max(abs(x["N_eff_A"] - x["closed_form_N_eff_A"]) for x in rows),
        "max_abs_N_eff_sum_minus_33": max(abs(x["N_eff_sum"] - 33.0) for x in rows),
        "max_recon_residual": max(x["recon_residual"] for x in rows),
        "density_recon_max_error": main_run["density_recon_max_error"],
        "max_abs_weight_b0_in_A_minus_closed_form": max(abs(x["weight_b0_in_A"] - x["closed_form_transfer_fraction"]) for x in rows),
        "first_period": {"k_of_min_eps_state_in_1..143": k_min_first, "min_eps_state": eps[k_min_first]},
        "key_steps": {str(k): {kk: rows[k][kk] for kk in ("exchange_phase_deg", "N_eff_A", "N_eff_B", "peak_raw_A", "peak_raw_B",
                                                       "eps_state", "eps_state_gp", "eps_swap", "weight_b0_in_A")}
                      for k in KEY_STEPS if k <= N_COLL},
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"elapsed_sec {time.time() - t0:.1f}", file=sys.stderr)  # 所要時間は出力ファイルに残さない（再現性）


if __name__ == "__main__":
    main()
