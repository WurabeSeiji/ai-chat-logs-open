#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第2予備実験：初期期間 N–分解能掃引（完全無seed）。解釈なし。

完全無seed（Z0=v）。制御変数は N と 有限分解能 Δ のみ。
1 step 順序（確定）: Cayley → [測定 before_quant] → 量子化 Q_Δ → [測定 after_quant] →
                     polar retraction 1回（第7論文 retract, 不変更） → [測定 after_reproj]。
Cayley 直後に retract は挿入しない（測定のみ）。retraction 補正量 ‖Zr−Zq‖_2 を生データ保存。
resolution_operator=OFF 基準は量子化・retract を適用しない（Cayley のみ、条件A と同一）。

早期停止(§8.2): f_outside≥1e-2 / NaN,Inf / closure_after_reproj>1e-8 / norm_after_reproj∉[1±1e-10] / max_step。
保存(§9): step 0..1000 毎step, 1001.. 5step毎, 停止step 必ず。全Z: {0,1,2,5,10,20,50,100,200,500,1000,停止}。

使い方（単一 run）: python3 run_initial_resolution_scan_v1.py N p Dref exec_idx    (Dref=OFF で基準run)
       全 run:    python3 run_initial_resolution_scan_v1.py ALL
"""
import csv
import json
import platform
import sys
import traceback
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
P2 = CODE.parent
PAPER8 = P2.parent
REPO = PAPER8.parent.parent
ENGINE = REPO / "時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1"
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
for pth in (ENGINE, V2 / "code"):
    sys.path.insert(0, str(pth))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_plane_flow_approx_v1 import parent_plane_split_approx
from run_n300_dimension_saturation_v2 import gram_reduce, dominant_plane
from run_transverse_stability_v1 import retract           # 第7論文 polar retraction（不変更）

N_REF = 40
M_REF = N_REF * (N_REF - 1) // 2
SIG_REL = 1e-6
Q_REL_TAU = 1e-8
MAX_STEP = {5: 2500, 40: 4500, 300: 10000}
P_VALUES = [0.0, 0.5, 1.0, 1.5, 2.0]
DREF_VALUES = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
N_VALUES = [5, 40, 300]
FULL_Z_STEPS = {0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000}
CLOSURE_LIMIT = 1e-8
NORM_LO, NORM_HI = 1 - 1e-10, 1 + 1e-10
F_LIMIT = 1e-2


def M_of(n):
    return n * (n - 1) // 2


def delta_actual(n, dref, p):
    return dref * (M_of(n) / M_REF) ** (-p / 2.0)


def Q_delta(Z, d):
    return d * np.round(Z.real / d) + 1j * (d * np.round(Z.imag / d))   # round = banker's (half to even)


def occ(B, Z):
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))


def qsv4(B0, Bd):
    ev = np.clip(np.linalg.eigvalsh(np.column_stack([B0, Bd]).T @ np.column_stack([B0, Bd]))[::-1], 0, None)
    return np.sqrt(ev)


def build_v(n):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    if n <= 40:
        _, B_p1, B_rot, _ = parent_plane_split_exact(sys_lr, v)
    else:
        _, B_p1, B_rot, _, _ = parent_plane_split_approx(sys_lr, v, SIG_REL)
    gr0 = gram_reduce(sys_lr, v); _, B0, _, _, _ = dominant_plane(sys_lr, gr0)
    wp = rng.normal(size=sys_lr.m)
    return sys_lr, v, B_p1, B0, wp


def observe(sys_lr, B_p1, B0, Z):
    totZ = float(np.real(np.conj(Z) @ Z))
    E_P1 = occ(B_p1, Z)
    f_out = max(0.0, 1 - E_P1 / totZ)
    a_out = float(np.sqrt(f_out))
    gr = gram_reduce(sys_lr, Z); _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
    E_dom = occ(Bdom, Z) / totZ
    qs = qsv4(B0, Bdom)
    rankQ = int(np.sum(qs > Q_REL_TAU * qs[0]))
    absZ = np.abs(Z)
    nz = absZ[absZ > 0]
    return {"f_outside": f_out, "a_outside": a_out, "q3": float(qs[2]), "q4": float(qs[3]),
            "rank_Q": rankQ, "E_dom": E_dom,
            "nonzero_real_count": int(np.sum(Z.real != 0.0)), "nonzero_imag_count": int(np.sum(Z.imag != 0.0)),
            "min_nonzero_abs_component": (float(nz.min()) if nz.size else 0.0),
            "max_abs_component": float(absZ.max())}


def f_outside_cheap(B_p1, Z):
    totZ = float(np.real(np.conj(Z) @ Z))
    return max(0.0, 1 - occ(B_p1, Z) / totZ)


def run_one(n, p, dref, exec_idx, res_off=False):
    if res_off:
        run_id = f"N{n:05d}_baseline_exec{exec_idx}"
        d = None
    else:
        run_id = f"N{n:05d}_p{p:.1f}_dref{dref:.0e}_exec{exec_idx}"
        d = delta_actual(n, dref, p)
    outdir = P2 / "raw" / run_id; outdir.mkdir(parents=True, exist_ok=True)
    fi = np.finfo(np.float64)
    import datetime  # noqa (timestamp via time not allowed; use env-free)
    cfg = {"run_id": run_id, "execution_index": exec_idx, "N": n, "M": M_of(n), "N_ref": N_REF, "M_ref": M_REF,
           "p": (None if res_off else p), "Delta_ref": (None if res_off else dref),
           "Delta_actual": d, "initial_seed": False, "metastable_seed": False, "random_kick": False,
           "external_noise": False, "resolution_operator": ("OFF" if res_off else "ON"),
           "rounding_mode": "half_to_even", "max_step": MAX_STEP[n],
           "save_schedule": "0..1000 every step; 1001.. every 5; stop step always",
           "stop_conditions": ["f_outside>=1e-2", "nan_inf", "closure_after_reproj>1e-8",
                               "norm_after_reproj outside [1-1e-10,1+1e-10]", "max_step"],
           "python_version": platform.python_version(), "numpy_version": np.__version__,
           "scipy_version": __import__("scipy").__version__, "platform": platform.platform(),
           "float_type": "float64", "mantissa_bits": int(fi.nmant), "machine_epsilon": float(fi.eps),
           "smallest_normal": float(fi.tiny), "smallest_subnormal": float(np.nextafter(0, 1)),
           "source_hashes_ref": "config/source_file_hashes.json"}
    with open(outdir / "run_config.json", "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2, ensure_ascii=False)

    sout = open(outdir / "stdout.log", "w"); serr = open(outdir / "stderr.log", "w")
    hdr = ["run_id", "N", "M", "p", "Delta_ref", "Delta_actual", "resolution_operator", "step", "time",
           "f_outside", "a_outside", "log_a_outside", "q3", "q4", "rank_Q", "E_dom",
           "closure_residual_before_quantization", "closure_residual_after_quantization",
           "closure_residual_after_reprojection", "norm_before_quantization", "norm_after_quantization",
           "norm_after_reprojection", "quantization_l2", "retraction_correction_l2",
           "nonzero_real_count", "nonzero_imag_count", "min_nonzero_abs_component", "max_abs_component", "stop_reason"]
    f_ts = open(outdir / "timeseries.csv", "w", newline=""); w = csv.writer(f_ts); w.writerow(hdr)
    fmt = "%.12e"
    save_rows = []
    zsteps = FULL_Z_STEPS
    z_saved = {}

    sys_lr, v, B_p1, B0, wp = build_v(n)
    Z = v.copy()
    stop_reason = None
    max_step = MAX_STEP[n]

    def is_save(t):
        return (t <= 1000) or (t % 5 == 0)

    def write_row(t, meas):
        ztz = complex(Z @ Z)
        a = meas["obs"]["a_outside"]
        loga = ("" if a <= 0 else fmt % np.log(a))
        w.writerow([run_id, n, M_of(n), (("" if res_off else fmt % p)), (("" if res_off else fmt % dref)),
                    (("" if d is None else fmt % d)), ("OFF" if res_off else "ON"), t, t,
                    fmt % meas["obs"]["f_outside"], fmt % a, loga, fmt % meas["obs"]["q3"], fmt % meas["obs"]["q4"],
                    meas["obs"]["rank_Q"], fmt % meas["obs"]["E_dom"],
                    fmt % meas["cbq"], fmt % meas["caq"], fmt % meas["car"],
                    fmt % meas["nbq"], fmt % meas["naq"], fmt % meas["nar"],
                    fmt % meas["ql2"], fmt % meas["rcorr"],
                    meas["obs"]["nonzero_real_count"], meas["obs"]["nonzero_imag_count"],
                    fmt % meas["obs"]["min_nonzero_abs_component"], fmt % meas["obs"]["max_abs_component"],
                    (stop_reason or "")])

    try:
        # t=0: Z0=v。step ではないので before/after/reproj は Z0 の値、ql2=rcorr=0
        obs0 = observe(sys_lr, B_p1, B0, Z)
        c0 = abs(complex(Z @ Z)); n0 = float(np.linalg.norm(Z))
        meas0 = {"obs": obs0, "cbq": c0, "caq": c0, "car": c0, "nbq": n0, "naq": n0, "nar": n0, "ql2": 0.0, "rcorr": 0.0}
        write_row(0, meas0); z_saved[0] = Z.copy()
        t = 0
        while t < max_step:
            # 1 step
            sys_lr.set_theta(np.angle(Z)); se, wp = sys_lr.sigma_max_power(wp)
            Zc = sys_lr.cayley_step(Z, se)                      # Cayley
            cbq = abs(complex(Zc @ Zc)); nbq = float(np.linalg.norm(Zc))   # before quant（測定のみ, retractなし）
            if res_off:
                Zr = Zc; caq = cbq; naq = nbq; ql2 = 0.0; car = cbq; nar = nbq; rcorr = 0.0
            else:
                Zq = Q_delta(Zc, d)                             # 量子化
                caq = abs(complex(Zq @ Zq)); naq = float(np.linalg.norm(Zq)); ql2 = float(np.linalg.norm(Zq - Zc))
                Zr = retract(Zq)                                # polar retraction 1回
                car = abs(complex(Zr @ Zr)); nar = float(np.linalg.norm(Zr)); rcorr = float(np.linalg.norm(Zr - Zq))
            Z = Zr; t += 1
            # 早期停止判定（数値例外）
            if not (np.all(np.isfinite(Z.real)) and np.all(np.isfinite(Z.imag))):
                stop_reason = "numerical_exception"
            f_now = f_outside_cheap(B_p1, Z)
            if stop_reason is None and f_now >= F_LIMIT:
                stop_reason = "f_outside_limit"
            if stop_reason is None and car > CLOSURE_LIMIT:
                stop_reason = "closure_residual_limit"
            if stop_reason is None and not (NORM_LO <= nar <= NORM_HI):
                stop_reason = "norm_limit"
            if stop_reason is None and t >= max_step:
                stop_reason = "max_step"
            save = is_save(t) or (stop_reason is not None)
            if save:
                obs = observe(sys_lr, B_p1, B0, Z)
                meas = {"obs": obs, "cbq": cbq, "caq": caq, "car": car, "nbq": nbq, "naq": naq, "nar": nar,
                        "ql2": ql2, "rcorr": rcorr}
                write_row(t, meas)
                if t in zsteps:
                    z_saved[t] = Z.copy()
            if stop_reason is not None:
                if t not in z_saved:
                    z_saved[t] = Z.copy()
                break
        if stop_reason is None:
            stop_reason = "max_step"
    except Exception as e:
        stop_reason = "numerical_exception"
        serr.write(traceback.format_exc())
    f_ts.close()

    # local_growth.csv（隣接保存点, 両端 a>0 のみ）
    rows = list(csv.DictReader(open(outdir / "timeseries.csv")))
    with open(outdir / "local_growth.csv", "w", newline="") as fh:
        wl = csv.writer(fh); wl.writerow(["step_from", "step_to", "a_from", "a_to", "gamma_local"])
        for r1, r2 in zip(rows[:-1], rows[1:]):
            a1 = float(r1["a_outside"]); a2 = float(r2["a_outside"]); t1 = int(r1["step"]); t2 = int(r2["step"])
            if a1 > 0 and a2 > 0 and t2 > t1:
                wl.writerow([t1, t2, "%.12e" % a1, "%.12e" % a2, "%.12e" % ((np.log(a2) - np.log(a1)) / (t2 - t1))])

    # 全Z保存（npz, 指定step+停止）
    if z_saved:
        np.savez_compressed(outdir / "state_vectors.npz",
                            steps=np.array(sorted(z_saved)),
                            **{f"Z_{s}": z_saved[s] for s in sorted(z_saved)})
    zh = None
    if stop_reason and (max(z_saved) if z_saved else None) is not None:
        import hashlib
        zf = z_saved[max(z_saved)]
        zh = hashlib.sha256(np.ascontiguousarray(zf).tobytes()).hexdigest()

    diag = {"run_id": run_id, "stop_reason": stop_reason, "stop_step": (max(z_saved) if z_saved else None),
            "n_saved_rows": len(rows), "final_Z_sha256": zh,
            "max_quantization_l2": max((float(r["quantization_l2"]) for r in rows), default=0.0),
            "max_retraction_correction_l2": max((float(r["retraction_correction_l2"]) for r in rows), default=0.0),
            "max_closure_after_reproj": max((float(r["closure_residual_after_reprojection"]) for r in rows), default=0.0),
            "max_norm_dev_after_reproj": max((abs(float(r["norm_after_reprojection"]) - 1) for r in rows), default=0.0),
            "final_f_outside": (float(rows[-1]["f_outside"]) if rows else None),
            "final_a_outside": (float(rows[-1]["a_outside"]) if rows else None)}
    with open(outdir / "run_diagnostics.json", "w", encoding="utf-8") as fh:
        json.dump(diag, fh, indent=2, ensure_ascii=False)
    sout.write(f"{run_id} stop={stop_reason} step={diag['stop_step']} rows={len(rows)}\n"); sout.close(); serr.close()
    return diag


def all_run_ids():
    ids = []
    for n in N_VALUES:
        for p in P_VALUES:
            for dref in DREF_VALUES:
                ids.append((n, p, dref, False))
    for n in N_VALUES:
        ids.append((n, None, None, True))
    return ids


def main():
    if sys.argv[1].upper() == "ALL":
        ids = all_run_ids()
        for (n, p, dref, off) in ids:
            for ex in (1, 2):
                dg = run_one(n, p, dref, ex, res_off=off)
                tag = (f"N{n} base" if off else f"N{n} p{p} d{dref:.0e}") + f" ex{ex}"
                print(f"[{tag}] stop={dg['stop_reason']} step={dg['stop_step']} "
                      f"final_a={dg['final_a_outside']:.2e} qL2max={dg['max_quantization_l2']:.1e} "
                      f"rcorrmax={dg['max_retraction_correction_l2']:.1e}")
        print("=== ALL DONE ===")
    else:
        n = int(sys.argv[1]); off = (sys.argv[3].upper() == "OFF")
        p = None if off else float(sys.argv[2]); dref = None if off else float(sys.argv[3])
        ex = int(sys.argv[4])
        dg = run_one(n, p, dref, ex, res_off=off)
        print(dg)


if __name__ == "__main__":
    main()
