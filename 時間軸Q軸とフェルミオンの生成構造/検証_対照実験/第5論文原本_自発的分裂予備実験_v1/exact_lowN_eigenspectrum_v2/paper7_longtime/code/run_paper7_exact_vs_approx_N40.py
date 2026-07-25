#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 §12 N=40 厳密法 vs 低ランク近似法 比較（5色占有・f・保存・射影閉鎖）。解釈なし。

固定親基底を厳密 eig(parent_plane_split_exact) と低ランク JG(parent_plane_split_approx) で
それぞれ構成し、同一軌道の 5色占有・分裂量 f・保存誤差・射影閉鎖を全サンプル時刻で比較。
B_dom(gram) は両者共通。横安定性は gram（低ランク）で両者同一。
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
P7 = CODE.parent; V2 = P7.parent; ENGINE = V2.parent
sys.path.insert(0, str(ENGINE)); sys.path.insert(0, str(V2 / "code")); sys.path.insert(0, str(CODE))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_plane_flow_approx_v1 import parent_plane_split_approx
from run_n300_dimension_saturation_v2 import gram_reduce, dominant_plane
from run_paper7_5color_timeseries import occ, s4_new_dirs

DELTA = 1e-15; XMAX = 55000; SAMPLE = 100; SIG_REL = 1e-6


def five_color(sys_lr, Zr, B_p1, B_rot, B0):
    totZ = float(np.real(np.conj(Zr) @ Zr))
    E_P1 = occ(B_p1, Zr); E_other = occ(B_rot, Zr); E_ker = totZ - E_P1 - E_other
    gr = gram_reduce(sys_lr, Zr); _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
    e34 = s4_new_dirs(B0, Bdom); proj = B_rot @ (B_rot.T @ e34)
    fq, _ = np.linalg.qr(proj); f34 = fq[:, :2]
    E_d3 = occ(f34[:, [0]], Zr); E_d4 = occ(f34[:, [1]], Zr)
    E_rem = max(0.0, E_other - E_d3 - E_d4)
    return np.array([E_P1, E_d3, E_d4, E_rem, E_ker]) / totZ, 1 - E_P1 / totZ, abs(totZ - 1.0)


def run():
    n = 40
    sys_lr = LowRankSystem(n); M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    _, Bp1_e, Brot_e, spec = parent_plane_split_exact(sys_lr, v)
    _, Bp1_a, Brot_a, smax, thr = parent_plane_split_approx(sys_lr, v, SIG_REL)
    gr0 = gram_reduce(sys_lr, v); _, B0, _, _, _ = dominant_plane(sys_lr, gr0)
    g = zero_closure_kernel_seed(sys_lr, rng); Z = v + DELTA * g; Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real); q = v.imag - (v.imag @ p) * p; q = q / np.linalg.norm(q)
    wp = rng.normal(size=M)

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))
    Zc = Z.copy(); wpc = wp.copy(); crossing = None; t = 0
    while True:
        if fval(Zc) > 0.05:
            crossing = t; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc); Zc = sys_lr.cayley_step(Zc, se); t += 1

    max_dev = {"bands": 0.0, "f": 0.0, "closure_exact": 0.0, "closure_approx": 0.0}
    Zr = Z.copy(); wpr = wp.copy(); t = 0
    while True:
        if t % SAMPLE == 0 or t == XMAX:
            be, fe, ce = five_color(sys_lr, Zr, Bp1_e, Brot_e, B0)
            ba, fa, ca = five_color(sys_lr, Zr, Bp1_a, Brot_a, B0)
            max_dev["bands"] = max(max_dev["bands"], float(np.max(np.abs(be - ba))))
            max_dev["f"] = max(max_dev["f"], abs(fe - fa))
            max_dev["closure_exact"] = max(max_dev["closure_exact"], ce)
            max_dev["closure_approx"] = max(max_dev["closure_approx"], ca)
        if t >= XMAX:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr); Zr = sys_lr.cayley_step(Zr, se); t += 1

    res = {"N": n, "crossing": crossing,
           "dims_exact": {"P1": int(Bp1_e.shape[1]), "other": int(Brot_e.shape[1])},
           "dims_approx": {"P1": int(Bp1_a.shape[1]), "other": int(Brot_a.shape[1])},
           "max_5color_band_deviation_exact_vs_approx": max_dev["bands"],
           "max_f_deviation": max_dev["f"],
           "max_conservation_error_exact": max_dev["closure_exact"],
           "max_conservation_error_approx": max_dev["closure_approx"],
           "note": "B_dom(gram) と横安定性は両法で同一(低ランク)。差は固定親基底の eig vs JG のみ。"}
    with open(P7 / "summary" / "N00040_exact_vs_approx.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print(f"[§12 N=40 厳密vs近似] crossing一致={crossing}")
    print(f"  5色占有 最大偏差={max_dev['bands']:.2e}  f最大偏差={max_dev['f']:.2e}")
    print(f"  保存誤差 厳密={max_dev['closure_exact']:.2e} 近似={max_dev['closure_approx']:.2e}")
    return res


if __name__ == "__main__":
    run()
