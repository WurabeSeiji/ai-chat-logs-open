#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A 主張の数値をNPZから再算出し、既存報告の値と照合する v1

背景（この器を作った理由）:
  論文A の主張2（支持セル数・N_eff）、主張4（ピーク高）、主張5（コヒーレント振幅の
  対数則）は、既存の `結果報告_シード型別δ一括掃引_T42000_v1.md` に数値としては
  あるが、**それを算出する登録済みプログラムが存在しなかった**。
  シリーズ規約「各論文はシリーズ内プログラムで独立検算」を満たすため、
  正本NPZ から全主張の数値を再算出し、報告書の値と機械照合する器を用意する。

方針:
  - 新規走行はしない。既存の正本NPZ と result JSON のみを読む（read-only）
  - 報告書の値は REFERENCE に**先に**書き込んである（本ファイルの改変履歴で確認可能）。
    走行後に書き換えない
  - τ_space / τ_time は NPZ から母体と同じ規約（`first_true` = 0基点+1）で
    再計算し、**result JSON の値と一致するか**を独立に検査する
  - 一致しない項目は FAIL として記録し、黙って合わせない

出力: result_paperA_claims_v1.json

使い方: python3 aggregate_paperA_claims_v1.py
"""
from __future__ import annotations
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# ---- 母体と同一の宣言値 ----
T_LONG = 42000
WIN = (T_LONG // 2, T_LONG)          # 後半窓 [21000, 42000)
LEDGER_EVERY = 50
NN, NETA = 16, 8
K_PUMP = 2

DELTAS = [1e-15, 1e-08, 1e-04, 1e-03, 1e-02,
          0.03162277660168379, 0.04357, 0.1]
MODES = ["neutral", "electron", "fermion_family", "boson_family", "mixed"]
# 各シード型のセル数（母体 MODES 定義より）
N_CELLS = {"neutral": 1, "electron": 1, "fermion_family": 5,
           "boson_family": 3, "mixed": 8}
# 主張5 の弱域（報告書と同じ 3 点）
WEAK_DELTAS = [1e-08, 1e-04, 1e-03]

# ---- 照合基準（既存報告書からの転記・走行前に固定）----
REFERENCE = {
    "claim2_support": {          # δ=0.1・後半・正の支持セル数と N_eff
        "neutral": [16, 15.9994], "electron": [64, 63.9972],
        "fermion_family": [128, 127.9944], "mixed": [128, 127.9950],
        "boson_family": [4, 1.0606],
    },
    "claim4_peak_r_nopump_mixed": {   # 通常mixed の r_nopump 最大値
        "0.01": 0.633907, "0.0316228": 0.695137,
        "0.04357": 0.716139, "0.1": 0.677707,
    },
    "claim5_fit_Acoh": {"a": 9.892, "b": -48.611, "r2": 0.999868,
                        "rmse": 2.82, "max_resid": 7.60},
    "claim5_fit_Pseed": {"r2": 0.99220, "rmse": 21.66},
    "claim6_tau_space": {         # 通常mixed / 位相相殺 / F5
        "0.01": [74, 116, 116], "0.0316228": [29, 43, 43],
        "0.04357": [22, 33, 33], "0.1": [12, 17, 17],
    },
    "claim4_tau_space_vacuum": 1624,
}
TOL_REL = 1e-3       # 報告書の丸め桁に合わせた相対許容
TOL_ABS_INT = 0      # 整数量（更新回数・セル数）は完全一致を要求


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def first_true(mask) -> int | None:
    """母体 run_tb_nsweep_1to20_v1.first_true と同一規約（0基点 +1）。"""
    idx = np.flatnonzero(mask)
    return int(idx[0]) + 1 if len(idx) else None


def npz_path(mode: str, delta: float, suffix: str = "") -> Path:
    """正本の命名規約。δ=0.01 の mixed/neutral は当時の既定でタグ無し。"""
    tag = f"_rep-{suffix}" if suffix else ""
    if not suffix and abs(delta - 0.01) < 1e-18 and mode in ("mixed", "neutral"):
        return HERE / f"nsweep_{mode}_T{T_LONG}_N12_v2.npz"
    return HERE / f"nsweep_{mode}_T{T_LONG}_d{delta:g}{tag}_N12_v2.npz"


def result_json_path(mode: str, delta: float) -> Path:
    if abs(delta - 0.01) < 1e-18 and mode in ("mixed", "neutral"):
        return HERE / f"result_nsweep_{mode}_T{T_LONG}_v2.json"
    return HERE / f"result_nsweep_{mode}_T{T_LONG}_d{delta:g}_v2.json"


def measure(path: Path) -> dict:
    """1走行のNPZから、主張が使う量をすべて算出する。"""
    z = np.load(path)
    f2 = z["m_f2"]
    acq = z["m_acq"]
    tau_space = first_true(np.nan_to_num(f2, nan=-1) > 0.05)
    tau_time = first_true(acq)
    tau_space_vac = first_true(np.nan_to_num(z["v_f2"], nan=-1) > 0.05)

    # --- 主張2: 後半の 128 セル支持 ---
    led = z["rec_m_ledger"]                       # (snapshots, 16, 8)
    led_t = z["ledger_t"] if "ledger_t" in z.files else None
    if led_t is not None:
        sel = (led_t >= WIN[0]) & (led_t < WIN[1])
    else:
        t_axis = np.arange(led.shape[0]) * LEDGER_EVERY
        sel = (t_axis >= WIN[0]) & (t_axis < WIN[1])
    cell_mean = led[sel].mean(axis=0)              # (16, 8) セル別後半時間平均
    support = int(np.count_nonzero(cell_mean > 0.0))
    s1 = float(cell_mean.sum())
    s2 = float((cell_mean ** 2).sum())
    n_eff = (s1 * s1 / s2) if s2 > 0 else float("nan")

    # --- 主張3: 奇数帯パワーの全区間最大 ---
    odd_max = float(np.max(z["rec_m_odd_power"]))

    # --- 主張4: 混合率のピーク高 ---
    r_np = z["rec_m_r_nopump"]
    peak = float(np.nanmax(r_np))
    r_np0 = float(r_np[0])

    return {"file": path.name, "n_snapshots_in_window": int(sel.sum()),
            "tau_space": tau_space, "tau_time": tau_time,
            "tau_space_vacuum": tau_space_vac,
            "support_cells": support, "N_eff": n_eff,
            "odd_power_max": odd_max,
            "r_nopump_peak": peak, "r_nopump_initial": r_np0}


def linfit(x, y):
    """y = a + b·x の最小二乗（numpy.polyfit を使わず正規方程式で明示）。"""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    sx, sy = x.sum(), y.sum()
    sxx, sxy = (x * x).sum(), (x * y).sum()
    den = n * sxx - sx * sx
    b = (n * sxy - sx * sy) / den
    a = (sy - b * sx) / n
    pred = a + b * x
    resid = y - pred
    ss_res = float((resid ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {"a": float(a), "b": float(b),
            "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
            "rmse": float(math.sqrt(ss_res / n)),
            "max_abs_resid": float(np.abs(resid).max()), "n": int(n)}


def cmp_rel(got, ref, tol=TOL_REL):
    if got is None or ref is None:
        return {"got": got, "ref": ref, "ok": False}
    rel = abs(got - ref) / max(abs(ref), 1e-300)
    return {"got": got, "ref": ref, "rel": rel, "ok": bool(rel <= tol)}


def cmp_int(got, ref):
    return {"got": got, "ref": ref, "ok": bool(got == ref)}


def main() -> None:
    t0 = time.time()
    out = {"generator": {"script": Path(__file__).name,
                         "sha256": sha256(Path(__file__).resolve())},
           "declared": {"T": T_LONG, "window": list(WIN), "Nn": NN,
                        "Neta": NETA, "deltas": DELTAS, "modes": MODES,
                        "weak_deltas": WEAK_DELTAS, "tol_rel": TOL_REL},
           "reference_source": "結果報告_シード型別δ一括掃引_T42000_v1.md",
           "runs": {}, "checks": {}, "missing": []}

    # ---------------- 全 40 条件の測定 ----------------
    print("=== 40 条件の再算出 ===")
    print(f"{'mode':<16}{'δ':>12}{'τ_space':>9}{'τ_time':>8}"
          f"{'支持':>6}{'N_eff':>10}{'odd_max':>12}{'peak':>9}")
    for mode in MODES:
        for d in DELTAS:
            p = npz_path(mode, d)
            if not p.exists():
                out["missing"].append(p.name)
                continue
            m = measure(p)
            out["runs"][f"{mode}|{d:g}"] = m
            print(f"{mode:<16}{d:>12g}{str(m['tau_space']):>9}"
                  f"{str(m['tau_time']):>8}{m['support_cells']:>6}"
                  f"{m['N_eff']:>10.4f}{m['odd_power_max']:>12.3e}"
                  f"{m['r_nopump_peak']:>9.6f}")

    # ---------------- τ を result JSON と独立照合 ----------------
    tau_cross = []
    for mode in MODES:
        for d in DELTAS:
            key = f"{mode}|{d:g}"
            if key not in out["runs"]:
                continue
            rj = result_json_path(mode, d)
            if not rj.exists():
                continue
            rec = json.loads(rj.read_text())["N"]["12"]
            tau_cross.append({
                "key": key,
                "tau_space": cmp_int(out["runs"][key]["tau_space"],
                                     rec.get("tau_space")),
                "tau_time": cmp_int(out["runs"][key]["tau_time"],
                                    rec.get("tau_time")),
            })
    n_bad = sum(1 for c in tau_cross
                if not c["tau_space"]["ok"] or not c["tau_time"]["ok"])
    out["checks"]["tau_vs_result_json"] = {
        "n_compared": len(tau_cross), "n_mismatch": n_bad,
        "pass": n_bad == 0, "details": tau_cross}
    print(f"\n[照合] NPZ再計算 vs result JSON の τ: "
          f"{len(tau_cross)} 条件・不一致 {n_bad} → "
          f"{'PASS' if n_bad == 0 else 'FAIL'}")

    # ---------------- 主張2: 支持セル数と N_eff ----------------
    c2 = {}
    for mode, (ref_sup, ref_neff) in REFERENCE["claim2_support"].items():
        key = f"{mode}|0.1"
        if key not in out["runs"]:
            continue
        c2[mode] = {"support": cmp_int(out["runs"][key]["support_cells"], ref_sup),
                    "N_eff": cmp_rel(out["runs"][key]["N_eff"], ref_neff)}
    out["checks"]["claim2_support"] = {
        "pass": all(v["support"]["ok"] and v["N_eff"]["ok"] for v in c2.values()),
        "detail": c2}
    print(f"[主張2] 支持セル数・N_eff (δ=0.1): "
          f"{'PASS' if out['checks']['claim2_support']['pass'] else 'FAIL'}")

    # ---------------- 主張3: 純ボゾンの奇数帯は厳密に 0 ----------------
    b3 = {f"{d:g}": out["runs"][f"boson_family|{d:g}"]["odd_power_max"]
          for d in DELTAS if f"boson_family|{d:g}" in out["runs"]}
    out["checks"]["claim3_odd_exact_zero"] = {
        "max_over_all_deltas": max(b3.values()) if b3 else None,
        "pass": bool(b3) and all(v == 0.0 for v in b3.values()),
        "per_delta": b3}
    print(f"[主張3] 純ボゾンの奇数帯パワー最大 = "
          f"{max(b3.values()) if b3 else 'n/a'} → "
          f"{'PASS' if out['checks']['claim3_odd_exact_zero']['pass'] else 'FAIL'}")

    # ---------------- 主張4: ピーク高 ----------------
    c4 = {k: cmp_rel(out["runs"][f"mixed|{k}"]["r_nopump_peak"], v)
          for k, v in REFERENCE["claim4_peak_r_nopump_mixed"].items()
          if f"mixed|{k}" in out["runs"]}
    out["checks"]["claim4_peak"] = {
        "pass": all(v["ok"] for v in c4.values()), "detail": c4}
    print(f"[主張4] 混合率ピーク高: "
          f"{'PASS' if out['checks']['claim4_peak']['pass'] else 'FAIL'}")

    # ---------------- 主張5: コヒーレント振幅の対数則 ----------------
    xs_a, xs_p, ys = [], [], []
    pts = []
    for mode in MODES:
        n = N_CELLS[mode]
        for d in WEAK_DELTAS:
            key = f"{mode}|{d:g}"
            if key not in out["runs"]:
                continue
            ts = out["runs"][key]["tau_space"]
            if ts is None:
                continue
            a_coh = n * d
            p_seed = n * d * d
            xs_a.append(math.log(a_coh))
            xs_p.append(math.log(p_seed))
            ys.append(ts)
            pts.append({"mode": mode, "delta": d, "n_cells": n,
                        "A_coh": a_coh, "P_seed": p_seed, "tau_space": ts})
    fit_a = linfit(xs_a, ys)
    fit_p = linfit(xs_p, ys)
    ref_a = REFERENCE["claim5_fit_Acoh"]
    out["checks"]["claim5_fit"] = {
        "points": pts, "fit_A_coh": fit_a, "fit_P_seed": fit_p,
        "cmp": {"a": cmp_rel(fit_a["a"], ref_a["a"], 1e-3),
                "b": cmp_rel(fit_a["b"], ref_a["b"], 1e-3),
                "r2": cmp_rel(fit_a["r2"], ref_a["r2"], 1e-5),
                "rmse": cmp_rel(fit_a["rmse"], ref_a["rmse"], 5e-3),
                "P_seed_r2": cmp_rel(fit_p["r2"],
                                     REFERENCE["claim5_fit_Pseed"]["r2"], 1e-4),
                "P_seed_rmse": cmp_rel(fit_p["rmse"],
                                       REFERENCE["claim5_fit_Pseed"]["rmse"], 5e-3)},
    }
    out["checks"]["claim5_fit"]["pass"] = all(
        v["ok"] for v in out["checks"]["claim5_fit"]["cmp"].values())
    print(f"\n[主張5] τ_space = {fit_a['a']:.4f} {fit_a['b']:+.4f}·ln(A_coh)  "
          f"R²={fit_a['r2']:.6f} RMSE={fit_a['rmse']:.3f} n={fit_a['n']}")
    print(f"         P_seed で整理すると R²={fit_p['r2']:.6f} "
          f"RMSE={fit_p['rmse']:.3f} → 対数則は A_coh 側が優位")
    print(f"[主張5] 報告書との照合: "
          f"{'PASS' if out['checks']['claim5_fit']['pass'] else 'FAIL'}")

    # ---------------- 主張6: 位相相殺 ----------------
    c6 = {}
    for dk, (ref_ctl, ref_pb, ref_f5) in REFERENCE["claim6_tau_space"].items():
        d = float(dk)
        got = {}
        for label, suffix, mode in (("control", "cohmixv1-ctl1", "mixed"),
                                    ("phase_balanced", "pbmixv1-r1", "mixed"),
                                    ("F5", "", "fermion_family")):
            p = npz_path(mode, d, suffix)
            got[label] = measure(p)["tau_space"] if p.exists() else None
        c6[dk] = {"control": cmp_int(got["control"], ref_ctl),
                  "phase_balanced": cmp_int(got["phase_balanced"], ref_pb),
                  "F5": cmp_int(got["F5"], ref_f5),
                  "pb_equals_F5": bool(got["phase_balanced"] == got["F5"]),
                  "pb_differs_from_control": bool(
                      got["phase_balanced"] != got["control"])}
    out["checks"]["claim6_phase_balanced"] = {
        "pass": all(v["control"]["ok"] and v["phase_balanced"]["ok"]
                    and v["F5"]["ok"] and v["pb_equals_F5"]
                    and v["pb_differs_from_control"] for v in c6.values()),
        "detail": c6}
    print(f"[主張6] 位相相殺 τ_space（通常/相殺/F5）: "
          f"{'PASS' if out['checks']['claim6_phase_balanced']['pass'] else 'FAIL'}")

    # ---------------- 主張7: 分解能 N ----------------
    c7 = {}
    for mode in ("mixed", "neutral", "electron", "vacuum"):
        f = HERE / f"result_nsweep_{mode}_v2.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        rows = []
        for k in sorted(d["N"], key=int):
            r = d["N"][k]
            rows.append({"N": int(k), "M": r.get("M"),
                         "built": r.get("built", False),
                         "n_eff_med": r.get("n_eff_med"),
                         "align_med": r.get("align_med"),
                         "tau_space": r.get("tau_space"),
                         "matter_born": r.get("matter_born"),
                         "time_born": r.get("time_born")})
        c7[mode] = {"failed_N": d.get("failed_N"), "failed": d.get("failed"),
                    "rows": rows, "env": d.get("env", {})}
    out["checks"]["claim7_resolution"] = c7
    if "mixed" in c7:
        deg = [r["N"] for r in c7["mixed"]["rows"]
               if r["built"] and r["n_eff_med"] is not None and r["n_eff_med"] > 1.8]
        out["checks"]["claim7_degenerate_N"] = deg
        print(f"[主張7] 平面が縮退している N（n_eff>1.8・mixed）: {deg}")
        print(f"[主張7] 構築に失敗した N: {c7['mixed']['failed_N']}")

    # ---------------- 主張1: 三つの誕生の 4 行 ----------------
    def born(mode, d):
        k = f"{mode}|{d:g}"
        if k not in out["runs"]:
            return None
        m = out["runs"][k]
        return {"space": m["tau_space"] is not None,
                "matter": m["odd_power_max"] > 1e-30,
                "clock": m["tau_time"] is not None,
                "tau_space": m["tau_space"], "tau_time": m["tau_time"]}
    vac = HERE / f"nsweep_vacuum_T{T_LONG}_N12_v2.npz"
    row_vac = None
    if vac.exists():
        mv = measure(vac)
        row_vac = {"space": mv["tau_space"] is not None,
                   "matter": mv["odd_power_max"] > 1e-30,
                   "clock": mv["tau_time"] is not None,
                   "tau_space": mv["tau_space"]}
    out["checks"]["claim1_birth_table"] = {
        "シードなし": row_vac,
        "偶数kだけ(B3, δ=0.1)": born("boson_family", 0.1),
        "奇数kあり・強度1e-8(F5)": born("fermion_family", 1e-08),
        "奇数kあり・強度1e-4(F5)": born("fermion_family", 1e-04),
    }
    print("\n[主張1] 三つの誕生")
    print(f"  {'条件':<28}{'空間':>6}{'物質':>6}{'時計':>6}")
    for k, v in out["checks"]["claim1_birth_table"].items():
        if v is None:
            continue
        print(f"  {k:<28}{'○' if v['space'] else '×':>6}"
              f"{'○' if v['matter'] else '×':>6}{'○' if v['clock'] else '×':>6}")

    # ---------------- 総合 ----------------
    named = ["tau_vs_result_json", "claim2_support", "claim3_odd_exact_zero",
             "claim4_peak", "claim5_fit", "claim6_phase_balanced"]
    out["all_pass"] = all(out["checks"][k].get("pass") for k in named
                          if k in out["checks"])
    out["runtime_sec"] = time.time() - t0
    (HERE / "result_paperA_claims_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(f"\n=== 総合: {'全項目 PASS' if out['all_pass'] else 'FAIL あり'} ===")
    if out["missing"]:
        print(f"欠損NPZ {len(out['missing'])} 件: {out['missing'][:5]}")
    print(f"完了 {out['runtime_sec']:.0f}s → result_paperA_claims_v1.json")


if __name__ == "__main__":
    main()
