#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 図1/2/3 用 5色占有時系列（共通横軸 0..55000）。解釈なし。

論文6の固定親基底3分類（P1 / その他回転 other / 核 kernel）を維持し、確定指示により
「その他回転 other」を【新方向3・新方向4・残余その他回転】へ分解して計5色。
  新方向 = 時間依存 S4(t)=正規直交化[B0(=P1) | B_dom(t)] の B0 直交補2方向(e3,e4)を、
           固定 other 空間 B_rot へ射影・正規直交化した f3,f4。
  E_dir3=|Π_f3 Z|², E_dir4=|Π_f4 Z|², 残余other=E_other-E_dir3-E_dir4。P1・核は不変。
黒線 = 分裂量 f = 1 - E_P1（論文6と同一）。固有値順の色割当はしない。縮退平面は連続基底固定。

N=5,40: 固定親基底=parent_plane_split_exact, B_dom=gram(密行列一致検証済)。
N=300: 固定親基底=parent_plane_split_approx（σ_rel=1e-6, 論文6と同一）, B_dom=gram低ランク。

使い方: python3 run_paper7_5color_timeseries.py 5
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
P7 = CODE.parent                 # paper7_longtime/
V2 = P7.parent                   # exact_lowN_eigenspectrum_v2/
ENGINE = V2.parent               # 第5論文原本_..._v1/
sys.path.insert(0, str(ENGINE))
sys.path.insert(0, str(V2 / "code"))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_plane_flow_approx_v1 import parent_plane_split_approx
from run_n300_dimension_saturation_v2 import gram_reduce, dominant_plane

DELTA = 1e-15
XMAX = 55000                     # 共通横軸（絶対step）
SIG_REL = 1e-6
SAMPLE = {5: 25, 40: 25, 300: 100}


def occ(B, Z):
    if B is None or (hasattr(B, "shape") and B.shape[1] == 0):
        return 0.0
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))


def build(n):
    sys_lr = LowRankSystem(n); M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    if n <= 40:
        p1s, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    else:
        p1s, B_p1, B_rot, smax, thr = parent_plane_split_approx(sys_lr, v, SIG_REL)
    gr0 = gram_reduce(sys_lr, v)
    _, B0, _, _, _ = dominant_plane(sys_lr, gr0)   # = 親支配平面(gram)。B_p1 と同一平面
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g; Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p; q = q / np.linalg.norm(q)
    wp = rng.normal(size=M)
    return sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp


def s4_new_dirs(B0, Bdom):
    """S4=orthonormalize[B0|Bdom] の B0 直交補2方向 e3,e4 を返す。"""
    R = Bdom - B0 @ (B0.T @ Bdom)
    Qr, _ = np.linalg.qr(R)
    return Qr[:, :2]


def align_2d(f_prev, f_new):
    """f_new(M×2) を前時刻 f_prev へ 2×2 回転で整列（連続基底固定・色反転防止）。"""
    if f_prev is None:
        return f_new
    Ov = f_prev.T @ f_new                # 2×2
    U, _, Vt = np.linalg.svd(Ov)
    Rot = U @ Vt                          # 直交 2×2
    return f_new @ Rot.T


def run(n):
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = build(n)
    M = sys_lr.m
    method = "exact" if n <= 40 else "approx"

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    # crossing
    Zc = Z.copy(); wpc = wp.copy(); crossing = None; t = 0
    while True:
        if fval(Zc) > 0.05:
            crossing = t; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc); Zc = sys_lr.cayley_step(Zc, se); t += 1

    outdir = P7 / "raw" / f"N{n:05d}"; outdir.mkdir(parents=True, exist_ok=True)
    f_ts = open(outdir / "paper7_long_timeseries.csv", "w", newline=""); w = csv.writer(f_ts)
    w.writerow(["step", "time", "crossing_flag", "splitting_fraction",
                "direction_1_occupation", "direction_2_occupation",
                "direction_3_occupation", "direction_4_occupation",
                "other_rotating_occupation", "kernel_occupation", "occupation_sum",
                "plane_1_occupation", "plane_2_occupation",
                "norm_error", "conservation_error", "projection_closure_error"])
    fmt = "%.10e"
    se_ev = SAMPLE[n]
    f_prev = None
    max_close = 0.0
    Zr = Z.copy(); wpr = wp.copy(); t = 0
    while True:
        if t % se_ev == 0 or t == XMAX:
            totZ = float(np.real(np.conj(Zr) @ Zr))
            E_P1 = occ(B_p1, Zr)
            E_other = occ(B_rot, Zr)
            E_ker = totZ - E_P1 - E_other
            f = 1.0 - E_P1 / totZ
            # 支配平面(gram) → 新方向 e3,e4 → other空間へ射影 f3,f4
            gr = gram_reduce(sys_lr, Zr)
            _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
            e34 = s4_new_dirs(B0, Bdom)              # M×2, P1直交補
            proj = B_rot @ (B_rot.T @ e34)           # other空間へ射影
            fq, _ = np.linalg.qr(proj)
            f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
            f34 = align_2d(f_prev, f34); f_prev = f34
            E_d3 = occ(f34[:, [0]], Zr)
            E_d4 = occ(f34[:, [1]], Zr) if f34.shape[1] > 1 else 0.0
            E_rem_other = max(0.0, E_other - E_d3 - E_d4)
            # P1の2軸（bookkeeping）
            E_a1 = occ(B_p1[:, [0]], Zr); E_a2 = occ(B_p1[:, [1]], Zr)
            osum = (E_P1 + E_d3 + E_d4 + E_rem_other + E_ker) / totZ
            close = abs(osum - 1.0)
            max_close = max(max_close, close)
            plane1 = E_P1 / totZ                      # P1 平面占有
            plane2 = (E_d3 + E_d4) / totZ             # 新2次元平面占有
            w.writerow([t, t, int(t >= crossing), fmt % f,
                        fmt % (E_a1 / totZ), fmt % (E_a2 / totZ),
                        fmt % (E_d3 / totZ), fmt % (E_d4 / totZ),
                        fmt % (E_rem_other / totZ), fmt % (E_ker / totZ), fmt % osum,
                        fmt % plane1, fmt % plane2,
                        fmt % abs(totZ - 1.0), fmt % abs(totZ - 1.0), fmt % close])
        if t >= XMAX:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr); Zr = sys_lr.cayley_step(Zr, se); t += 1
    f_ts.close()

    summary = {"N": n, "M": M, "crossing": crossing, "xmax": XMAX, "sample_every": se_ev,
               "method_parent_basis": method, "dims_P1": int(B_p1.shape[1]),
               "dims_other": int(B_rot.shape[1]) if B_rot is not None else 0,
               "dims_kernel": int(M - B_p1.shape[1] - (B_rot.shape[1] if B_rot is not None else 0)),
               "max_projection_closure_error": max_close}
    (P7 / "summary").mkdir(exist_ok=True)
    with open(P7 / "summary" / f"N{n:05d}_5color_meta.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"[5色 N={n}] M={M} crossing={crossing} xmax={XMAX} sample_every={se_ev} "
          f"親基底dim(P1/other/ker)={summary['dims_P1']}/{summary['dims_other']}/{summary['dims_kernel']} "
          f"閉鎖誤差max={max_close:.1e}")
    return summary


if __name__ == "__main__":
    for a in sys.argv[1:]:
        run(int(a))
