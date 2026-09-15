#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""状態 Z(t) の決定論的再生成と scaled 一致ゲート（仕様書 v1.1 §2）。

凍結済み Phase A コード（../PhaseA_実装一式_v1_20260915/、無変更で import）から
同一の初期条件・同一の stepping（疎 Lanczos 主軌道＋密正本影走行）を再実行し、
今回は全状態 Z(t) を results/states/*.npz に保存する。

観測量の計算式は run_phase12.py の readout の忠実複製である（旧データは読み取り専用。
旧 CSV には一切書き込まない）。

合格ゲート（凍結・変更禁止）:
  float:  |x_regen - x_old| <= 1e-12 * max(1, |x_old|)   （NaN は双方 NaN で一致）
  int 列（t, valid_A, valid_B, n_A, n_B）: 厳密一致
1 項目でも不合格なら SVD 解析へ進まず停止する（exit 1）。

対象走行（仕様書 §8 の解析対象）:
  L6_den6_AB, L6_den40_AB, L12_den12_AB, L12_den40_AB（第一対象＋対照）
  L6_den6_B（正値対照）
"""
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PHASEA = os.path.abspath(os.path.join(HERE, '..', 'PhaseA_実装一式_v1_20260915'))
sys.path.insert(0, PHASEA)

from kernel_canonical_dense import adjacency, one_step          # noqa: E402
from kernel_sparse_matrixfree import SparseKernel               # noqa: E402
from generator_phaseA import generate, analysis_basis           # noqa: E402
from run_phase12 import single_mode_states, wrap, M_A, M_B      # noqa: E402
from run_phase12 import (EPS_FLOOR, A_MIN_OVER_NORM,            # noqa: E402
                         LANCZOS_TOL, T_STEPS)

STATES_DIR = os.path.join(HERE, 'results', 'states')
os.makedirs(STATES_DIR, exist_ok=True)
OLD_CSV_DIR = os.path.join(PHASEA, 'results', 'phase12')

GATE = 1e-12
INT_COLS = {'t', 'valid_A', 'valid_B', 'n_A', 'n_B'}
TARGETS = [(6, 6, 'AB'), (6, 40, 'AB'), (12, 12, 'AB'), (12, 40, 'AB'), (6, 6, 'B')]


def regen_run(L, den, init_name):
    """run_phase12.run_one と同一手順の再走行（観測式は忠実複製）。状態と行列を返す。"""
    A = adjacency(L)
    sk = SparseKernel(L)
    B = analysis_basis(L)
    G_inv_Bh = np.linalg.solve(B.conj().T @ B, B.conj().T)
    ord_a, ord_b = L // math.gcd(L, M_A), L // math.gcd(L, M_B)
    zA0, zB0 = single_mode_states(L)
    z0 = {'A': zA0, 'B': zB0, 'AB': generate(L)}[init_name]

    M = len(z0)
    normZ = np.linalg.norm(z0)
    a_min = A_MIN_OVER_NORM * normZ
    zs = z0.copy()
    zd = z0.copy()
    traj = np.empty((T_STEPS + 1, M), np.complex128)
    traj[0] = zs
    rows = []
    theta_ref = [None, None]
    c0 = G_inv_Bh @ z0
    if abs(c0[0]) > a_min:
        theta_ref[0] = float(np.angle(c0[0]))
    if abs(c0[1]) > a_min:
        theta_ref[1] = float(np.angle(c0[1]))

    for t in range(T_STEPS + 1):
        z = traj[t] if t == 0 else zs
        drift = float(np.linalg.norm(zs - zd) / normZ) if t > 0 else 0.0
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
        H = float(np.vdot(z, z).real)
        Q2 = complex(np.sum(z ** 2))
        if t > 0:
            d = z - z0
            eps_abs = float(np.linalg.norm(d) / normZ)
            ip = np.vdot(z0, z)
            phi = float(np.angle(ip))
            eps_u1 = float(np.sqrt(max(np.linalg.norm(z) ** 2 + normZ ** 2
                                       - 2 * abs(ip), 0.0)) / normZ)
        else:
            eps_abs = eps_u1 = phi = 0.0
        rows.append([t, H, Q2.real, Q2.imag, gen,
                     c[0].real, c[0].imag, c[1].real, c[1].imag,
                     aA, aB, eta, int(vA), int(vB), nA, nB, rA, rB,
                     eps_abs, eps_u1, phi, drift])
        if t < T_STEPS:
            zs = sk.one_step(zs, den, tol=LANCZOS_TOL)
            zd = one_step(zd, A, den)
            traj[t + 1] = zs
    return traj, rows


def gate_compare(run_id, rows):
    """旧 CSV（読み取り専用）と scaled gate で全列比較。"""
    path = os.path.join(OLD_CSV_DIR, f'{run_id}_timeseries.csv')
    with open(path, newline='') as f:
        r = csv.reader(f)
        header = next(r)
        old_rows = [row for row in r]
    assert len(old_rows) == len(rows)
    n_checked = 0
    worst = 0.0
    failures = []
    for i, (old, new) in enumerate(zip(old_rows, rows)):
        for j, col in enumerate(header):
            xo = float(old[j])
            xn = float(new[j])
            if col in INT_COLS:
                if int(xo) != int(xn):
                    failures.append((i, col, xo, xn))
            else:
                if math.isnan(xo) or math.isnan(xn):
                    if not (math.isnan(xo) and math.isnan(xn)):
                        failures.append((i, col, xo, xn))
                else:
                    tol = GATE * max(1.0, abs(xo))
                    d = abs(xn - xo)
                    worst = max(worst, d / max(1.0, abs(xo)))
                    if d > tol:
                        failures.append((i, col, xo, xn))
            n_checked += 1
    return n_checked, worst, failures


def main():
    report = {'gate': 'abs(x_new-x_old) <= 1e-12*max(1,abs(x_old)); int exact; NaN==NaN',
              'runs': []}
    ok_all = True
    for L, den, init in TARGETS:
        run_id = f'L{L}_den{den}_{init}'
        traj, rows = regen_run(L, den, init)
        n, worst, fails = gate_compare(run_id, rows)
        ok = len(fails) == 0
        ok_all &= ok
        np.savez_compressed(os.path.join(STATES_DIR, f'{run_id}_states.npz'),
                            Z=traj, L=np.int64(L), den=np.int64(den),
                            init=init, T=np.int64(T_STEPS))
        report['runs'].append({'run_id': run_id, 'n_values_checked': n,
                               'worst_scaled_diff': worst, 'n_failures': len(fails),
                               'gate_pass': ok,
                               'failures_head': [[int(i), c, xo, xn] for i, c, xo, xn in fails[:5]]})
        print(f'{run_id}: checked={n} worst_scaled_diff={worst:.3e} '
              f'{"PASS" if ok else "FAIL(" + str(len(fails)) + ")"}', flush=True)
    with open(os.path.join(HERE, 'results', 'gate_report.json'), 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=1)
    if not ok_all:
        print('GATE FAILED — SVD 解析へ進まず停止する（仕様書 v1.1 §2）')
        sys.exit(1)
    print('ALL GATES PASSED')


if __name__ == '__main__':
    main()
