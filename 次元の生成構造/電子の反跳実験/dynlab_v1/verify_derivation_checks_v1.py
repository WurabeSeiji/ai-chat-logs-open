#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""導出ノート群の数値検証の正本 v1 — 使い捨てスクリプトで取った検証値の再現・保存

経緯: 2026-08-19 の導出セッションで、以下の検証値をインラインの使い捨てスクリプトで
取得し、導出ノートに引用していた（規約違反）。本プログラムはそれらを全て再現し、
result_verification_checks_v1.json に保存する。ノートの検証表は本プログラムを出所とする。

検証項目:
  V1  DL0  正単体定理の数値（λ=1/(2M), √λ, tr(B)=1/N, r_rms, 占有率）
  V2  Dyn-1 距離レート式 vs 頂点実装（乱数系・全辺）
  V3  Dyn-1a Σḋ²=0（エルミート保存の距離版）
  V4  Dyn-1b 真空 R=0 で頂点レート厳密零
  V5  Dyn-3 固有値・固有ベクトル摂動（一次精度）
  V6  H4   毛ゲージ不変性（(+3,+3) vs (+3,−3) の χ周辺化パワー、±8 エイリアシング域）
  V7  DL6  静的場方程式の自己無撞着固定点写像の固有値 (N−2)/2（真空自明性）

全て決定的（固定シード）。使い方: python3 verify_derivation_checks_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


# ----------------------------------------------------------------------
def v1_dl0_simplex():
    N = 16
    M = N * (N - 1) // 2
    lam = 1.0 / (2 * M)
    return {
        "N": N, "M": M,
        "d2": 1.0 / M,
        "lambda": lam,
        "sqrt_lambda": float(np.sqrt(lam)),
        "trB": (N - 1) * lam,
        "trB_eq_1_over_N": abs((N - 1) * lam - 1.0 / N),
        "r_rms": float(np.sqrt((N - 1) * lam / N)),
        "top3_occupancy": 3.0 / (N - 1),
    }


# ----------------------------------------------------------------------
def _vertex_rate_toy(x, R, ia, ib, n):
    a2 = np.abs(x) ** 2
    z2 = x ** 2
    A = np.zeros(n); B = np.zeros(n, complex)
    AR = np.zeros(n); BR = np.zeros(n, complex)
    np.add.at(A, ia, a2); np.add.at(A, ib, a2)
    np.add.at(B, ia, z2); np.add.at(B, ib, z2)
    np.add.at(AR, ia, R * a2); np.add.at(AR, ib, R * a2)
    np.add.at(BR, ia, R * z2); np.add.at(BR, ib, R * z2)
    cA = A[ia] + A[ib] - 2 * a2
    cB = B[ia] + B[ib] - 2 * z2
    cAR = AR[ia] + AR[ib] - 2 * R * a2
    cBR = BR[ia] + BR[ib] - 2 * R * z2
    rate = 0.5j * (R * (cA * x - cB * np.conj(x)) + (cAR * x - cBR * np.conj(x)))
    return rate, B, BR


def v2_v3_v4_dyn1():
    rng = np.random.default_rng(7)
    n = 6
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    ia = np.array([e[0] for e in edges]); ib = np.array([e[1] for e in edges])
    m = len(edges)
    x = rng.normal(size=m) + 1j * rng.normal(size=m)
    R = rng.random(m)

    rate, B, BR = _vertex_rate_toy(x, R, ia, ib, n)
    dd2_num = 2 * np.real(np.conj(x) * rate)
    xb2 = np.conj(x) ** 2
    dd2_theory = R * np.imag((B[ia] + B[ib]) * xb2) + np.imag((BR[ia] + BR[ib]) * xb2)

    rate0, _, _ = _vertex_rate_toy(x, np.zeros(m), ia, ib, n)
    return {
        "V2_max_abs_num_minus_theory": float(np.max(np.abs(dd2_num - dd2_theory))),
        "V3_sum_dd2": float(dd2_num.sum()),
        "V4_vacuum_max_rate": float(np.max(np.abs(rate0))),
    }


# ----------------------------------------------------------------------
def v5_eigen_perturbation():
    rng = np.random.default_rng(11)
    n = 8
    Q = rng.normal(size=(n, n)); B0 = (Q + Q.T) / 2
    P = rng.normal(size=(n, n)); dB = (P + P.T) / 2 * 1e-6
    lam0, V0 = np.linalg.eigh(B0)
    lam1, V1 = np.linalg.eigh(B0 + dB)
    dlam_th = np.einsum("ij,jk,ki->i", V0.T, dB, V0)
    i = n - 1
    w = np.array([(V0[:, j] @ dB @ V0[:, i]) / (lam0[i] - lam0[j]) if j != i else 0.0
                  for j in range(n)])
    v_pred = V0[:, i] + V0 @ w
    v_meas = V1[:, i] * np.sign(V1[:, i] @ V0[:, i])
    return {
        "V5_eigenvalue_err": float(np.max(np.abs((lam1 - lam0) - dlam_th))),
        "V5_eigenvector_err": float(np.max(np.abs(v_pred - v_meas))),
    }


# ----------------------------------------------------------------------
def v6_hair_gauge_invariance():
    _uni = _load("uni_vfy", UNI / "unified_interaction_v1.py")
    _cr0 = _load("cr0_vfy", EXP / "run_cr0_control_no_theta_v2.py")
    base = _uni.two_body_base
    step = _uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    sl, ic, _ = _cr0.calibrate_shift(sp, nc, ne)
    P = tuple(range(1, 18))
    eta = 2 * np.pi * np.arange(ne) / ne

    def mk(mA, mB):
        case = base.explicit_packet_case(
            mode=f"vfy{mA}_{mB}", packet_a=P, packet_b=P,
            packet_a_shift=_cr0.shift_for_deg(-30, sl, ic),
            packet_b_shift=_cr0.shift_for_deg(30, sl, ic))
        a = base.make_case_state(sp, case, "A", hair_enabled=True)
        b = base.make_case_state(sp, case, "B", hair_enabled=True)
        a = (a.reshape(nc, ne) * np.exp(1j * mA * eta)).reshape(-1)
        b = (b.reshape(nc, ne) * np.exp(1j * mB * eta)).reshape(-1)
        return a, b

    def chi_power(psi):
        return np.sum(np.abs(psi.reshape(nc, ne)) ** 2, axis=1)

    def run_pair(mA1, mB1, mA2, mB2, T=50):
        a1, b1 = mk(mA1, mB1); a2, b2 = mk(mA2, mB2)
        for _ in range(T):
            a1, b1, _ = step(a1, b1, sp)
            a2, b2, _ = step(a2, b2, sp)
        return float(np.max(np.abs(chi_power(a1) - chi_power(a2))))

    _, _, wA = _cr0.winding_spectrum(mk(3, 3)[0], nc, ne)
    return {
        "V6_diff_pp_vs_pm_m3": run_pair(3, 3, 3, -3),
        "V6_diff_pp_vs_pm_m8_alias": run_pair(8, 8, 8, -8),
        "V6_winding_A_after_imprint3": float(wA),
    }


# ----------------------------------------------------------------------
def v7_field_eq_fixed_point():
    out = {}
    for N in (3, 4, 5, 6, 8, 16):
        A = np.zeros((N, N))
        for v in range(N):
            for u in range(N):
                if u != v:
                    A[v, v] += 0.5
                    A[v, u] += 0.5
        J = np.eye(N) - np.ones((N, N)) / N
        lam = np.linalg.eigvalsh(J @ A @ J)
        nz = [float(x) for x in sorted(set(np.round(lam, 10))) if abs(x) > 1e-9]
        out[f"N{N}"] = {"nontrivial_eigs": nz, "theory": (N - 2) / 2.0}
    return out


# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    res = {
        "V1_dl0_simplex": v1_dl0_simplex(),
        "V2_V3_V4_dyn1": v2_v3_v4_dyn1(),
        "V5_eigen_perturbation": v5_eigen_perturbation(),
        "V6_hair_gauge_invariance": v6_hair_gauge_invariance(),
        "V7_field_eq_fixed_point": v7_field_eq_fixed_point(),
    }
    res["elapsed_sec"] = time.time() - t0
    (HERE / "result_verification_checks_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
