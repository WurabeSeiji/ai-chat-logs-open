#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 第1予備実験：二段階seed除去対照（条件A/B/D）。指示書に忠実・解釈なし。

条件A: initial seed OFF / metastable seed OFF（Z0 = v。kernel seed 生成で乱数を消費しない）
条件B: initial seed ON  / metastable seed OFF（Z0 = (v+δg)/‖·‖, δ=1e-15。以後自然発展のみ）
条件D: initial seed ON  / metastable seed ON （B と t1 直前までビット一致。t1=crossing+3000 で
       単一横摂動 ε η_⊥ を一回注入 → 規格化 → 以後 B と同一の自然発展。Benettin/再注入/再正規化なし）
       ε=1e-8, transverse seed index=0（run_paper7_transverse.py seed=0 の η_⊥ 生成をそのまま流用、
       初回生成のみ。Benettin ループは使用しない）。

第7論文コードは read-only import（不変更）。seed の ON/OFF は本ラッパーで明示切替。
共通最終 step = 55000。crossing = f>0.05 の最初、t0/t1 = crossing+GUARD(3000)。

使い方: python3 run_preliminary_seed_ablation_v1.py 5 A
        python3 run_preliminary_seed_ablation_v1.py 40 B
        python3 run_preliminary_seed_ablation_v1.py 300 D
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
PAPER8 = CODE.parent
REPO = PAPER8.parent.parent
ENGINE = REPO / "時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1"
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
PL = V2 / "paper7_longtime"
for pth in (ENGINE, V2 / "code", PL / "code"):
    sys.path.insert(0, str(pth))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_plane_flow_approx_v1 import parent_plane_split_approx
from run_n300_dimension_saturation_v2 import gram_reduce, dominant_plane
from run_paper7_5color_timeseries import s4_new_dirs
from run_paper7_transverse import s4_basis, perp

DELTA = 1e-15
XMAX = 55000                 # COMMON_FINAL_STEP
GUARD = 3000                 # crossing → t1（第7論文と同一）
SIG_REL = 1e-6
Q_REL_TAU = 1e-8             # rank_Q 判定閾値（第7論文と同一）
SAMPLE = {5: 25, 40: 25, 300: 100}
D_EPS = 1e-8                 # 条件D 固定
D_SEED_INDEX = 0            # 条件D 固定

CONDITIONS = {
    "A": {"initial_seed": False, "metastable_seed": False, "file": "condition_A_no_seed"},
    "B": {"initial_seed": True, "metastable_seed": False, "file": "condition_B_initial_only"},
    "D": {"initial_seed": True, "metastable_seed": True, "file": "condition_D_existing_two_seed"},
}


def occ(B, Z):
    if B is None or B.shape[1] == 0:
        return 0.0
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))


def qsv4(B0, Bd):
    Q4 = np.column_stack([B0, Bd])
    ev = np.clip(np.linalg.eigvalsh(Q4.T @ Q4)[::-1], 0, None)
    return np.sqrt(ev)


def build_init(n, initial_seed):
    """第7論文 build と同一。ただし initial_seed=False では kernel seed g を生成せず乱数を消費しない。"""
    sys_lr = LowRankSystem(n); M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    if n <= 40:
        p1s, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    else:
        p1s, B_p1, B_rot, smax, thr = parent_plane_split_approx(sys_lr, v, SIG_REL)
    gr0 = gram_reduce(sys_lr, v)
    _, B0, _, _, _ = dominant_plane(sys_lr, gr0)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p; q = q / np.linalg.norm(q)
    if initial_seed:
        g = zero_closure_kernel_seed(sys_lr, rng)         # 乱数消費（B/D）
        Z0 = v + DELTA * g; Z0 = Z0 / np.linalg.norm(Z0)
    else:
        Z0 = v.copy()                                      # 無seed（乱数を消費しない）
    wp = rng.normal(size=M)
    return sys_lr, v, B_p1, B_rot, B0, p, q, Z0, wp


def evolve(sys_lr, Z, wp):
    sys_lr.set_theta(np.angle(Z)); se, wp = sys_lr.sigma_max_power(wp)
    return sys_lr.cayley_step(Z, se), wp


def run(n, cond):
    c = CONDITIONS[cond]
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = build_init(n, c["initial_seed"])
    M = sys_lr.m

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    # v の規格化・零二乗閉鎖 診断（条件A用にも）
    v_diag = {"v_norm_error": abs(float(np.real(np.conj(v) @ v)) - 1.0),
              "v_zero_square_abs": abs(complex(v @ v))}

    outdir = PAPER8 / "raw" / f"N{n:05d}"; outdir.mkdir(parents=True, exist_ok=True)
    f_ts = open(outdir / f"{c['file']}.csv", "w", newline=""); w = csv.writer(f_ts)
    w.writerow(["step", "time", "N", "condition", "initial_seed_enabled", "metastable_seed_enabled",
                "initial_seed_amplitude", "metastable_seed_amplitude", "parent_plane_occupation",
                "f_outside_parent", "q1", "q2", "q3", "q4", "rank_Q", "dominant_plane_occupation",
                "non_dominant_occupation", "kernel_occupation", "residual_occupation", "norm_Z",
                "dagger_norm_error", "zero_square_real", "zero_square_imag", "zero_square_abs",
                "projection_closure_error", "crossing_detected", "metastable_start_detected"])
    fmt = "%.10e"
    se_ev = SAMPLE[n]
    init_amp = DELTA if c["initial_seed"] else 0.0
    meta_amp = D_EPS if c["metastable_seed"] else 0.0

    # 診断集約
    dg = {"max_norm_error": 0.0, "max_zero_square_abs": 0.0, "max_closure_error": 0.0,
          "max_antisym_error": 0.0, "crossing_step": None, "metastable_start_step": None,
          "injected_at": None, "injection_eps": (D_EPS if cond == "D" else None),
          "injection_seed_index": (D_SEED_INDEX if cond == "D" else None)}

    crossing = None
    t = 0
    while True:
        f = fval(Z)
        if crossing is None and f > 0.05:
            crossing = t; dg["crossing_step"] = t
        t1 = (crossing + GUARD) if crossing is not None else None
        # 条件D：t1 で単一横摂動を一回だけ注入
        if cond == "D" and t1 is not None and t == t1 and dg["injected_at"] is None:
            rng_dir = np.random.default_rng(70000 + n)          # 第7論文 seed=0 と同一
            eta_r = rng_dir.normal(size=M); eta_i = rng_dir.normal(size=M)   # 初回生成（index 0）
            S4_t1 = s4_basis(sys_lr, B0, Z)
            eta_r = eta_r - S4_t1 @ (S4_t1.T @ eta_r); eta_i = eta_i - S4_t1 @ (S4_t1.T @ eta_i)
            eta = (eta_r + 1j * eta_i) / np.sqrt(eta_r @ eta_r + eta_i @ eta_i)
            Z = Z + D_EPS * eta
            Z = Z / np.linalg.norm(Z)                            # 注入直後に一度だけ規格化
            dg["injected_at"] = t; dg["metastable_start_step"] = t
            f = fval(Z)
        elif t1 is not None and dg["metastable_start_step"] is None:
            dg["metastable_start_step"] = t1

        if t % se_ev == 0 or t == XMAX:
            totZ = float(np.real(np.conj(Z) @ Z))
            E_P1 = occ(B_p1, Z); E_other = occ(B_rot, Z); E_ker = totZ - E_P1 - E_other
            gr = gram_reduce(sys_lr, Z); _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
            E_dom = occ(Bdom, Z)
            qs = qsv4(B0, Bdom)
            rankQ = int(np.sum(qs > Q_REL_TAU * qs[0]))
            ztz = complex(Z @ Z)
            closure = abs(totZ - E_P1 - E_other - E_ker)      # 恒等的に≈0
            K = None  # 反対称誤差は診断で別途（軽量化のため記録時のみ）
            residual = closure
            met_start = int(dg["metastable_start_step"] is not None and t >= dg["metastable_start_step"])
            w.writerow([t, t, n, cond, int(c["initial_seed"]), int(c["metastable_seed"]),
                        fmt % init_amp, fmt % meta_amp, fmt % (E_P1 / totZ), fmt % (1 - E_P1 / totZ),
                        fmt % qs[0], fmt % qs[1], fmt % qs[2], fmt % qs[3], rankQ,
                        fmt % (E_dom / totZ), fmt % (E_other / totZ), fmt % (E_ker / totZ),
                        fmt % (residual / totZ), fmt % np.sqrt(totZ),
                        fmt % abs(totZ - 1.0), fmt % ztz.real, fmt % ztz.imag, fmt % abs(ztz),
                        fmt % (residual / totZ), int(crossing is not None), met_start])
            dg["max_norm_error"] = max(dg["max_norm_error"], abs(totZ - 1.0))
            dg["max_zero_square_abs"] = max(dg["max_zero_square_abs"], abs(ztz))
            dg["max_closure_error"] = max(dg["max_closure_error"], residual / totZ)
        if t >= XMAX:
            break
        Z, wp = evolve(sys_lr, Z, wp); t += 1
    f_ts.close()

    # §13 反対称誤差（終端で1回, 軽量）
    Kfull = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)]) if M <= 800 else None
    antisym = (float(np.linalg.norm(Kfull + Kfull.T) / max(1.0, np.linalg.norm(Kfull, 2)))
               if Kfull is not None else None)  # N=300は密行列不可のため None

    diag = {"N": n, "condition": cond, "M": M,
            "initial_seed_enabled": c["initial_seed"], "metastable_seed_enabled": c["metastable_seed"],
            "v_normalized_zero_square": v_diag, "crossing_step": dg["crossing_step"],
            "metastable_start_step": dg["metastable_start_step"], "injection": {
                "injected_at": dg["injected_at"], "eps": dg["injection_eps"],
                "seed_index": dg["injection_seed_index"]},
            "max_norm_error": dg["max_norm_error"], "max_zero_square_abs": dg["max_zero_square_abs"],
            "max_projection_closure_error": dg["max_closure_error"],
            "antisymmetry_error_at_final": antisym, "xmax": XMAX, "sample_every": se_ev}
    (PAPER8 / "diagnostics").mkdir(exist_ok=True)
    with open(PAPER8 / "diagnostics" / f"N{n:05d}_condition_{cond}.json", "w", encoding="utf-8") as fh:
        json.dump(diag, fh, indent=2, ensure_ascii=False)
    print(f"[N={n} 条件{cond}] crossing={dg['crossing_step']} t1/準安定開始={dg['metastable_start_step']} "
          f"注入={dg['injected_at']} 規格化誤差max={dg['max_norm_error']:.1e} 零二乗max={dg['max_zero_square_abs']:.1e} "
          f"閉鎖max={dg['max_closure_error']:.1e}")
    return diag


if __name__ == "__main__":
    n = int(sys.argv[1]); cond = sys.argv[2].upper()
    assert cond in CONDITIONS, "condition must be A/B/D"
    run(n, cond)
