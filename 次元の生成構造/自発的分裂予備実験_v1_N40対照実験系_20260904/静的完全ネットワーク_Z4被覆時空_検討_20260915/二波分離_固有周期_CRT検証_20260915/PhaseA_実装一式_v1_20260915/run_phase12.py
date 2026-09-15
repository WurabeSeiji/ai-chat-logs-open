#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 / Phase 2 走行（仕様書 v1.2 §12、ChatGPT Phase 0 承認指示 2026-09-15）。

Phase 1: 単一モード初期値 Z_A(0)=c_a v_2, Z_B(0)=c_b v_3 の軌道観測。
Phase 2: 二モード Z_AB(0)=Z_A(0)+Z_B(0) を別走行し、非線形相互作用差
    Delta_NL(t) = ||F^t(Z_AB) - [F^t(Z_A)+F^t(Z_B)]|| / ||Z_AB(0)||
  を診断量として保存（重ね合わせ成立の期待量ではない）。

系: L in {6,12} x den in {L(主), 40(副)} x init in {A, B, AB}。T=500 step。
主軌道は疎 matrix-free（Lanczos tol=1e-14）、密正本コピーを影走行として
dense_sparse_drift(t) を別監査量で保存（閾値には算入しない）。

閾値（Phase 0 確定・変更禁止）:
  eps_floor = 7.60908171255871e-14
  eps_close = 7.60908171255871e-12
  eps_move  = 7.60908171255871e-8
  a_min/||Z|| = 7.60908171255871e-11
"""
import csv
import json
import math
import os

import numpy as np

from kernel_canonical_dense import adjacency, one_step
from kernel_sparse_matrixfree import SparseKernel
from generator_phaseA import generate, analysis_basis
from kernel_canonical_dense import edges

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'results', 'phase12')
os.makedirs(OUTDIR, exist_ok=True)

EPS_FLOOR = 7.60908171255871e-14
EPS_CLOSE = 7.60908171255871e-12
EPS_MOVE = 7.60908171255871e-8
A_MIN_OVER_NORM = 7.60908171255871e-11
LANCZOS_TOL = 1e-14
T_STEPS = 500
M_A, M_B = 2, 3


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def single_mode_states(L):
    ea, eb = edges(L)
    delta = (eb - ea).astype(np.float64)
    U = np.exp(2j * np.pi / L)
    c_a = np.exp(-1j * np.pi / (2 * L))
    c_b = np.exp(+1j * np.pi / (2 * L))
    zA = (c_a * U ** (M_A * delta)).astype(np.complex128)
    zB = (c_b * U ** (M_B * delta)).astype(np.complex128)
    return zA, zB


def run_one(L, den, init_name, z0, sk, A, Binfo):
    """1 走行。疎主軌道＋密影走行。時系列 CSV と要約 dict を返す。"""
    B, G_inv_Bh, ord_a, ord_b = Binfo
    M = len(z0)
    normZ = np.linalg.norm(z0)
    a_min = A_MIN_OVER_NORM * normZ
    H0 = float(np.vdot(z0, z0).real)
    Q20 = complex(np.sum(z0 ** 2))

    zs = z0.copy()
    zd = z0.copy()
    traj = np.empty((T_STEPS + 1, M), np.complex128)
    traj[0] = zs
    rows = []
    theta_ref = [None, None]
    cell_visits = set()
    n_valid_cov = 0
    n_excluded = 0
    stats = dict(drift_max=0.0, eta_max=0.0, gen_min=np.inf, gen_max=0.0,
                 dH_max=0.0, dQ2_max=0.0,
                 aA_min=np.inf, aB_min=np.inf,
                 rA_max=0.0, rB_max=0.0, rA_sum=0.0, rB_sum=0.0, r_cnt=0)

    def readout(t, z, drift):
        c = G_inv_Bh @ z
        aA, aB = abs(c[0]), abs(c[1])
        eta = float(np.linalg.norm(z - B @ c) / np.linalg.norm(z))
        gen = float(np.linalg.norm(sk.K_apply(np.exp(1j * np.angle(z)), z)) / np.linalg.norm(z))
        vA, vB = aA > a_min, aB > a_min
        thA = float(np.angle(c[0])) if vA else np.nan
        thB = float(np.angle(c[1])) if vB else np.nan
        nA = nB = -1
        rA = rB = np.nan
        if vA and theta_ref[0] is not None:
            tt = wrap(thA - theta_ref[0])
            nA = int(round(tt / (2 * np.pi / ord_a))) % ord_a
            rA = float(abs(wrap(tt - 2 * np.pi * nA / ord_a)))
        if vB and theta_ref[1] is not None:
            tt = wrap(thB - theta_ref[1])
            nB = int(round(tt / (2 * np.pi / ord_b))) % ord_b
            rB = float(abs(wrap(tt - 2 * np.pi * nB / ord_b)))
        return c, aA, aB, eta, gen, vA, vB, thA, thB, nA, nB, rA, rB

    # t=0: 位相基準の設定（有効な場合のみ）
    c0 = G_inv_Bh @ z0
    if abs(c0[0]) > a_min:
        theta_ref[0] = float(np.angle(c0[0]))
    if abs(c0[1]) > a_min:
        theta_ref[1] = float(np.angle(c0[1]))

    eps_abs = np.empty(T_STEPS + 1)
    eps_u1 = np.empty(T_STEPS + 1)
    phis = np.empty(T_STEPS + 1)
    eps_abs[0] = eps_u1[0] = 0.0
    phis[0] = 0.0

    for t in range(T_STEPS + 1):
        z = traj[t] if t == 0 else zs
        drift = float(np.linalg.norm(zs - zd) / normZ) if t > 0 else 0.0
        c, aA, aB, eta, gen, vA, vB, thA, thB, nA, nB, rA, rB = readout(t, z, drift)
        H = float(np.vdot(z, z).real)
        Q2 = complex(np.sum(z ** 2))
        if t > 0:
            d = z - z0
            eps_abs[t] = float(np.linalg.norm(d) / normZ)
            ip = np.vdot(z0, z)
            phis[t] = float(np.angle(ip))
            eps_u1[t] = float(np.sqrt(max(np.linalg.norm(z) ** 2 + normZ ** 2
                                          - 2 * abs(ip), 0.0)) / normZ)
        stats['drift_max'] = max(stats['drift_max'], drift)
        stats['eta_max'] = max(stats['eta_max'], eta)
        stats['gen_min'] = min(stats['gen_min'], gen)
        stats['gen_max'] = max(stats['gen_max'], gen)
        stats['dH_max'] = max(stats['dH_max'], abs(H - H0) / H0)
        stats['dQ2_max'] = max(stats['dQ2_max'], abs(Q2 - Q20) / max(H0, abs(Q20)))
        stats['aA_min'] = min(stats['aA_min'], aA)
        stats['aB_min'] = min(stats['aB_min'], aB)
        if vA and vB and nA >= 0 and nB >= 0:
            cell_visits.add((nA, nB))
            n_valid_cov += 1
        else:
            n_excluded += 1
        if not math.isnan(rA):
            stats['rA_max'] = max(stats['rA_max'], rA)
            stats['rA_sum'] += rA
        if not math.isnan(rB):
            stats['rB_max'] = max(stats['rB_max'], rB)
            stats['rB_sum'] += rB
        if not (math.isnan(rA) and math.isnan(rB)):
            stats['r_cnt'] += 1
        rows.append([t, H, Q2.real, Q2.imag, gen,
                     c[0].real, c[0].imag, c[1].real, c[1].imag,
                     aA, aB, eta, int(vA), int(vB), nA, nB, rA, rB,
                     eps_abs[t], eps_u1[t], phis[t], drift])
        if t < T_STEPS:
            zs = sk.one_step(zs, den, tol=LANCZOS_TOL)
            zd = one_step(zd, A, den)
            traj[t + 1] = zs

    # q_min（自明固定点排除つき、仕様 §6/§7）
    def qmin_of(eps):
        for q in range(2, T_STEPS + 1):
            if eps[q] < EPS_CLOSE:
                D = float(np.max(eps_abs[1:q]))
                if D > EPS_MOVE:
                    return q, D
        return None, None

    qa, Da = qmin_of(eps_abs)
    qu, Du = qmin_of(eps_u1)
    D_overall = float(np.max(eps_abs[1:]))

    run_id = f'L{L}_den{den}_{init_name}'
    with open(os.path.join(OUTDIR, f'{run_id}_timeseries.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t', 'H', 'Q2_re', 'Q2_im', 'gen_norm',
                    'cA_re', 'cA_im', 'cB_re', 'cB_im', 'a_A', 'a_B', 'eta_perp',
                    'valid_A', 'valid_B', 'n_A', 'n_B', 'r_A', 'r_B',
                    'eps_abs', 'eps_u1', 'phi_star', 'dense_sparse_drift'])
        w.writerows(rows)

    denom_cells = ord_a * ord_b
    summary = dict(
        run_id=run_id, L=L, den=den, init=init_name, T=T_STEPS,
        H0=H0, Q2_0=[Q20.real, Q20.imag],
        q_min_abs=qa, D_orbit_at_q_abs=Da,
        q_min_u1=qu, D_orbit_at_q_u1=Du,
        phi_star_at_q_u1=(float(phis[qu]) if qu else None),
        eps_abs_min_q=float(np.min(eps_abs[2:])), eps_abs_argmin=int(np.argmin(eps_abs[2:]) + 2),
        eps_u1_min_q=float(np.min(eps_u1[2:])), eps_u1_argmin=int(np.argmin(eps_u1[2:]) + 2),
        D_orbit_overall=D_overall,
        coverage=(len(cell_visits) / denom_cells),
        cells_visited=len(cell_visits), cells_total=denom_cells,
        n_valid_cov=n_valid_cov, n_excluded=n_excluded,
        rA_max=stats['rA_max'], rB_max=stats['rB_max'],
        rA_mean=(stats['rA_sum'] / max(stats['r_cnt'], 1)),
        rB_mean=(stats['rB_sum'] / max(stats['r_cnt'], 1)),
        aA_min=stats['aA_min'], aB_min=stats['aB_min'],
        eta_perp_max=stats['eta_max'],
        gen_norm_min=stats['gen_min'], gen_norm_max=stats['gen_max'],
        dH_rel_max=stats['dH_max'], dQ2_rel_max=stats['dQ2_max'],
        dense_sparse_drift_max=stats['drift_max'])
    return traj, summary


def main():
    summaries = []
    for L in (6, 12):
        A = adjacency(L)
        sk = SparseKernel(L)
        B = analysis_basis(L)
        G_inv_Bh = np.linalg.solve(B.conj().T @ B, B.conj().T)
        ord_a, ord_b = L // math.gcd(L, M_A), L // math.gcd(L, M_B)
        Binfo = (B, G_inv_Bh, ord_a, ord_b)
        zA0, zB0 = single_mode_states(L)
        zAB0 = generate(L)
        assert np.allclose(zAB0, zA0 + zB0, atol=1e-15)
        for den in (L, 40):
            trajs = {}
            for name, z0 in (('A', zA0), ('B', zB0), ('AB', zAB0)):
                traj, s = run_one(L, den, name, z0, sk, A, Binfo)
                trajs[name] = traj
                summaries.append(s)
                print(f"{s['run_id']}: eps_u1_min={s['eps_u1_min_q']:.3e}@q={s['eps_u1_argmin']} "
                      f"q_min_u1={s['q_min_u1']} cov={s['coverage']:.3f} "
                      f"eta_max={s['eta_perp_max']:.3e} drift={s['dense_sparse_drift_max']:.3e}", flush=True)
            # Phase 2: 非線形相互作用差 Delta_NL(t)
            nAB = np.linalg.norm(trajs['AB'][0])
            dNL = np.linalg.norm(trajs['AB'] - (trajs['A'] + trajs['B']), axis=1) / nAB
            with open(os.path.join(OUTDIR, f'L{L}_den{den}_delta_NL.csv'), 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['t', 'delta_NL'])
                w.writerows([[t, float(dNL[t])] for t in range(T_STEPS + 1)])
            for s in summaries:
                if s['L'] == L and s['den'] == den and s['init'] == 'AB':
                    s['delta_NL_max'] = float(np.max(dNL))
                    s['delta_NL_final'] = float(dNL[-1])
                    s['delta_NL_first_exceed_1e-3'] = (int(np.argmax(dNL > 1e-3))
                                                       if np.any(dNL > 1e-3) else None)
            print(f"L={L} den={den}: Delta_NL max={np.max(dNL):.3e} final={dNL[-1]:.3e}", flush=True)

    meta = dict(eps_floor=EPS_FLOOR, eps_close=EPS_CLOSE, eps_move=EPS_MOVE,
                a_min_over_normZ=A_MIN_OVER_NORM, lanczos_tol=LANCZOS_TOL,
                T=T_STEPS, m_a=M_A, m_b=M_B,
                note='閾値は Phase 0 確定値・変更禁止。dense_sparse_drift は別監査量（閾値に不算入）。')
    with open(os.path.join(HERE, 'results', 'phase12_summary.json'), 'w', encoding='utf-8') as f:
        json.dump({'meta': meta, 'runs': summaries}, f, indent=1)
    print('PHASE12 DONE')


if __name__ == '__main__':
    main()
