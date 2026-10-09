"""R=0.70 の再帰交換（0713 の図 1 の元）で、各衝突の正規化が何をしているかを記録する。

元プログラム original_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
（0713 の原本と一字一句同一）をモジュールとして読み込み、その関数だけを使う。
ループは原本の recursive_snapshot_states（920〜983 行）と同じ更新

    a_next = normalize(r * a + t * b)
    b_next = normalize(t * a + r * b)

を行い、記録用の行だけを足す。renormalize=False の対照では normalize を外し、更新を純粋な U_R にする。

段階:
  1) 対照: renormalize=True の記録が、原本の recursive_snapshot_states と同じ値を出すか（衝突 0,1,2,3,5,10,20,42）
  2) 記録: 衝突 0..128 の正規化前ノルム、A・B の重なり、ピークの絶対高さ、N_eff、L
  3) 比較: renormalize=False（純粋な U_R）との差
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
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

R_VALUE = 0.70
N_COLL = 128
CONTROL_COLLISIONS = [0, 1, 2, 3, 5, 10, 20, 42]


def raw_peak(params, v):
    """ピークの絶対高さ: 総和で割らない chi 密度 sum_eta |psi|^2 の最大値。"""
    psi = orig.reshape_state(params, v)
    return float(np.max(np.sum(np.abs(psi) ** 2, axis=1)))


def trace(renormalize: bool):
    params = orig.Params()
    delta_f = orig.delta_from_reflection_rate(R_VALUE)
    t, r, T, R = orig.scattering_coefficients(delta_f)
    hair_enabled = True
    a = orig.make_state(params, 1, params.q_A, params.m_A, hair_enabled, params.A_A)
    b = orig.make_state(params, 63, params.q_B, params.m_B, hair_enabled, params.A_B)
    initial_a, initial_b = a.copy(), b.copy()
    rows = []
    for k in range(N_COLL + 1):
        ma = orig.recursive_state_metrics(params, "trace", f"R{int(round(R_VALUE*100)):03d}", delta_f, T, R, 1, 63,
                                          hair_enabled, k, "A_channel", a, initial_a, initial_b)
        mb = orig.recursive_state_metrics(params, "trace", f"R{int(round(R_VALUE*100)):03d}", delta_f, T, R, 1, 63,
                                          hair_enabled, k, "B_channel", b, initial_a, initial_b)
        ov = orig.inner(a, b)
        na_pre = r * a + t * b
        nb_pre = t * a + r * b
        rows.append({
            "collision": k,
            "norm2_A": orig.norm2(a), "norm2_B": orig.norm2(b),
            "overlap_abs": abs(ov), "overlap_re": ov.real, "overlap_im": ov.imag,
            "next_prenorm_norm2_A": orig.norm2(na_pre), "next_prenorm_norm2_B": orig.norm2(nb_pre),
            "next_prenorm_norm2_sum": orig.norm2(na_pre) + orig.norm2(nb_pre),
            "peak_raw_A": raw_peak(params, a), "peak_raw_B": raw_peak(params, b),
            "N_eff_A": ma["N_eff"], "N_eff_B": mb["N_eff"], "L_A": ma["L"], "L_B": mb["L"],
        })
        if renormalize:
            a_next = orig.normalize(r * a + t * b)
            b_next = orig.normalize(t * a + r * b)
        else:
            a_next = r * a + t * b
            b_next = t * a + r * b
        a, b = a_next, b_next
    return rows, {"t": [t.real, t.imag], "r": [r.real, r.imag], "T": T, "R": R}


def main():
    norm_rows, coeff = trace(True)
    pure_rows, _ = trace(False)

    # 1) 対照: 原本の関数と一致するか
    snaps = orig.recursive_snapshot_states(orig.Params(), R_VALUE, CONTROL_COLLISIONS)
    control = []
    for k in CONTROL_COLLISIONS:
        s = snaps[k]
        row = norm_rows[k]
        d = {
            "collision": k,
            "N_eff_A_orig": s["A_channel"]["N_eff"], "N_eff_A_trace": row["N_eff_A"],
            "N_eff_B_orig": s["B_channel"]["N_eff"], "N_eff_B_trace": row["N_eff_B"],
            "L_A_orig": s["A_channel"]["L"], "L_A_trace": row["L_A"],
            "L_B_orig": s["B_channel"]["L"], "L_B_trace": row["L_B"],
        }
        d["identical"] = (d["N_eff_A_orig"] == d["N_eff_A_trace"] and d["N_eff_B_orig"] == d["N_eff_B_trace"]
                          and d["L_A_orig"] == d["L_A_trace"] and d["L_B_orig"] == d["L_B_trace"])
        control.append(d)

    # 3) 比較
    dev_pre = max(max(abs(x["next_prenorm_norm2_A"] - 1.0), abs(x["next_prenorm_norm2_B"] - 1.0)) for x in norm_rows)
    sum_dev = max(abs(x["next_prenorm_norm2_sum"] - 2.0) for x in norm_rows)
    neff_diff = max(max(abs(p["N_eff_A"] - q["N_eff_A"]), abs(p["N_eff_B"] - q["N_eff_B"])) for p, q in zip(norm_rows, pure_rows))
    l_diff = max(max(abs(p["L_A"] - q["L_A"]), abs(p["L_B"] - q["L_B"])) for p, q in zip(norm_rows, pure_rows))
    pure_norm_range = (min(min(x["norm2_A"], x["norm2_B"]) for x in pure_rows), max(max(x["norm2_A"], x["norm2_B"]) for x in pure_rows))
    pure_sum_dev = max(abs(x["norm2_A"] + x["norm2_B"] - 2.0) for x in pure_rows)

    summary = {
        "R": R_VALUE, "collisions": N_COLL, "coefficients": coeff,
        "control_all_identical": all(c["identical"] for c in control),
        "max_abs_overlap": max(x["overlap_abs"] for x in norm_rows),
        "max_prenorm_channel_norm2_deviation_from_1": dev_pre,
        "max_prenorm_norm2_sum_deviation_from_2": sum_dev,
        "pure_unitary_channel_norm2_range": pure_norm_range,
        "pure_unitary_norm2_sum_max_deviation_from_2": pure_sum_dev,
        "max_N_eff_difference_normalized_vs_pure": neff_diff,
        "max_L_difference_normalized_vs_pure": l_diff,
        "peak_raw_A_at": {k: norm_rows[k]["peak_raw_A"] for k in CONTROL_COLLISIONS},
        "peak_raw_B_at": {k: norm_rows[k]["peak_raw_B"] for k in CONTROL_COLLISIONS},
    }

    def write(name, rows):
        with open(OUT / name, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write("trace_renormalized_R070.csv", norm_rows)
    write("trace_pure_unitary_R070.csv", pure_rows)
    write("control_vs_original_snapshots.csv", control)
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ks = [x["collision"] for x in norm_rows]
    fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    ax[0].plot(ks, [x["next_prenorm_norm2_A"] - 1 for x in norm_rows], label="A: |rA+tB|^2 - 1 (before normalize)")
    ax[0].plot(ks, [x["next_prenorm_norm2_B"] - 1 for x in norm_rows], label="B: |tA+rB|^2 - 1 (before normalize)")
    ax[0].set_ylabel("pre-normalization norm^2 - 1"); ax[0].legend(fontsize=8)
    ax[0].set_title("R=0.70: what the per-collision normalize removes")
    ax[1].plot(ks, [x["N_eff_A"] for x in norm_rows], label="A, normalized (original)")
    ax[1].plot(ks, [x["N_eff_B"] for x in norm_rows], label="B, normalized (original)")
    ax[1].plot(ks, [x["N_eff_A"] for x in pure_rows], "--", label="A, pure U_R")
    ax[1].plot(ks, [x["N_eff_B"] for x in pure_rows], "--", label="B, pure U_R")
    ax[1].set_ylabel("N_eff"); ax[1].legend(fontsize=8)
    ax[2].plot(ks, [x["peak_raw_A"] for x in norm_rows], label="A peak of sum_eta |psi|^2 (absolute)")
    ax[2].plot(ks, [x["peak_raw_B"] for x in norm_rows], label="B peak of sum_eta |psi|^2 (absolute)")
    ax[2].set_ylabel("absolute peak height"); ax[2].set_xlabel("collision"); ax[2].legend(fontsize=8)
    fig.savefig(OUT / "normalization_trace_R070.png", dpi=150)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
