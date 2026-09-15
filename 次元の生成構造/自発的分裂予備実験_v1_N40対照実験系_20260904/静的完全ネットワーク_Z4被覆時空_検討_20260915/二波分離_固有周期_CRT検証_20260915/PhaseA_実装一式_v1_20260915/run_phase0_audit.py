#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 0: 力学カーネル監査（仕様書 v1.2 §12 Phase 0・§10・§14）。

検査項目:
  (1) A=0 ⇒ F=Id
  (2) 全位相同一 ⇒ F=Id
  (3) θ差 ∈ {0,π}（実状態）⇒ F=Id
  (4) S_P 頂点置換同変性
  (5) U(1) 同変性
  (6) 符号反転同変性
  (7) Q2, H の保存（500 step ドリフト）
  (8) 生成器の Q2(0), H(0), min|z| を §5/§4.2 事前登録値と照合
  (9) 密（正本コピー）vs 疎（matrix-free Lanczos）の一致
 (10) ε_floor の実測と閾値候補の確定

対象: L=6, 12（主系列サイズ）× den ∈ {L（主）, 40（副 Δτ0=2π/40）}。
出力: results/phase0_results.json（全数値）と標準出力サマリ。
"""
import json
import math
import os

import numpy as np

from kernel_canonical_dense import edges, adjacency, one_step
from kernel_sparse_matrixfree import SparseKernel
from generator_phaseA import generate

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'results')
os.makedirs(RESULTS, exist_ok=True)

LANCZOS_TOL = 1e-14
STEPS_DRIFT = 500
N_PERM = 20
N_PHASE = 20
RNG = np.random.default_rng(20260915)

# §5 事前登録値（仕様書 v1.2、独立検算済み）
PREREG = {
    6: dict(Q2=complex(3.526279, -2.892305), H=19.607695, min_z=0.517638),
    12: dict(Q2=complex(-22.694534, 6.561456), H=108.817780, min_z=0.261052),
}


def edge_perm_from_vertex_perm(P, pi):
    """頂点置換 pi が誘導する辺置換（canonical 再整列、値は共役なしで移動）。"""
    ea, eb = edges(P)
    idx = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(ea, eb))}
    sigma = np.empty(len(ea), dtype=np.int64)
    for i, (a, b) in enumerate(zip(ea, eb)):
        na, nb = int(pi[a]), int(pi[b])
        if na > nb:
            na, nb = nb, na
        sigma[i] = idx[(na, nb)]
    return sigma  # 新índice: (P_sigma z)[sigma[i]] = z[i]


def apply_edge_perm(sigma, z):
    out = np.empty_like(z)
    out[sigma] = z
    return out


def relres(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def audit(L, den):
    P = L
    A = adjacency(P)
    sk = SparseKernel(P)
    M = sk.M
    res = {'L': L, 'den': den, 'M': M}

    z0 = generate(L)
    zr = RNG.standard_normal(M) + 1j * RNG.standard_normal(M)
    zr /= np.linalg.norm(zr)
    tests = {'generator': z0 / np.linalg.norm(z0), 'random': zr}

    # (1) A=0 -> F=Id（密: A=0 行列 / 疎: A_apply を 0 に差し替え）
    zA = tests['random']
    r_dense = relres(one_step(zA, np.zeros_like(A), den), zA)
    sk0 = SparseKernel(P)
    sk0.A_apply = lambda y: np.zeros_like(y)
    r_sparse = relres(sk0.one_step(zA, den, tol=LANCZOS_TOL), zA)
    res['t1_A0_identity'] = {'dense': r_dense, 'sparse': r_sparse}

    # (2) 全位相同一 -> F=Id
    amp = np.abs(RNG.standard_normal(M)) + 0.1
    zsame = amp * np.exp(1j * 0.7)
    res['t2_same_phase'] = {'dense': relres(one_step(zsame, A, den), zsame),
                            'sparse': relres(sk.one_step(zsame, den, tol=LANCZOS_TOL), zsame)}

    # (3) θ差 ∈ {0,π}（実状態）-> F=Id
    zreal = (RNG.standard_normal(M) + np.sign(RNG.standard_normal(M)) * 0.2).astype(np.complex128)
    res['t3_real_state'] = {'dense': relres(one_step(zreal, A, den), zreal),
                            'sparse': relres(sk.one_step(zreal, den, tol=LANCZOS_TOL), zreal)}

    # (4)(5)(6) 同変性（密で実施。疎は (9) で密と突合するため間接検証）
    eq_perm, eq_u1, eq_sign = [], [], []
    for name, z in tests.items():
        Fz = one_step(z, A, den)
        for _ in range(N_PERM):
            pi = RNG.permutation(P)
            sg = edge_perm_from_vertex_perm(P, pi)
            lhs = one_step(apply_edge_perm(sg, z), A, den)
            eq_perm.append(relres(lhs, apply_edge_perm(sg, Fz)))
        for _ in range(N_PHASE):
            ph = np.exp(1j * RNG.uniform(0, 2 * np.pi))
            eq_u1.append(relres(one_step(ph * z, A, den), ph * Fz))
        eq_sign.append(relres(one_step(-z, A, den), -Fz))
    res['t4_perm_equiv_max'] = float(np.max(eq_perm))
    res['t5_u1_equiv_max'] = float(np.max(eq_u1))
    res['t6_sign_equiv_max'] = float(np.max(eq_sign))

    # (7) 保存量ドリフト（500 step、密・疎それぞれ）
    drift = {}
    for label, stepper in (('dense', lambda z: one_step(z, A, den)),
                           ('sparse', lambda z: sk.one_step(z, den, tol=LANCZOS_TOL))):
        z = (z0 / np.linalg.norm(z0)).copy()
        H0 = float(np.vdot(z, z).real)
        Q20 = complex(np.sum(z ** 2))
        dH = dQ = 0.0
        for _ in range(STEPS_DRIFT):
            z = stepper(z)
            dH = max(dH, abs(np.vdot(z, z).real - H0) / H0)
            dQ = max(dQ, abs(complex(np.sum(z ** 2)) - Q20) / max(H0, abs(Q20)))
        drift[label] = {'dH_rel_max': float(dH), 'dQ2_rel_max': float(dQ)}
    res['t7_conservation_500step'] = drift

    # (8) 生成器の事前登録値照合（正規化前の生成式そのもの）
    zg = generate(L)
    got = dict(Q2=complex(np.sum(zg ** 2)), H=float(np.vdot(zg, zg).real),
               min_z=float(np.min(np.abs(zg))))
    reg = PREREG[L]
    res['t8_preregistration'] = {
        'Q2_diff': abs(got['Q2'] - reg['Q2']),
        'H_diff': abs(got['H'] - reg['H']),
        'min_z_diff': abs(got['min_z'] - reg['min_z']),
        'Q2_measured': [got['Q2'].real, got['Q2'].imag],
        'H_measured': got['H'], 'min_z_measured': got['min_z']}

    # (9) 密 vs 疎の一致（1 step および 50 step）
    z = (z0 / np.linalg.norm(z0)).copy()
    zd, zs = z.copy(), z.copy()
    one_diff = None
    for t in range(50):
        zd = one_step(zd, A, den)
        zs = sk.one_step(zs, den, tol=LANCZOS_TOL)
        if t == 0:
            one_diff = relres(zs, zd)
    res['t9_dense_vs_sparse'] = {'step1': one_diff, 'step50': relres(zs, zd)}

    return res


def main():
    all_res = []
    for L in (6, 12):
        for den in (L, 40):
            r = audit(L, den)
            all_res.append(r)
            print(f"L={L} den={den}: perm={r['t4_perm_equiv_max']:.2e} "
                  f"u1={r['t5_u1_equiv_max']:.2e} sign={r['t6_sign_equiv_max']:.2e} "
                  f"dH={r['t7_conservation_500step']['dense']['dH_rel_max']:.2e} "
                  f"dQ2={r['t7_conservation_500step']['dense']['dQ2_rel_max']:.2e} "
                  f"dvs={r['t9_dense_vs_sparse']['step50']:.2e}", flush=True)

    # (10) ε_floor 実測（仕様書 §10 の定義。全系・密疎の最大値で保守的に取る）
    comp = {}
    Mmax = max(r['M'] for r in all_res)
    comp['sqrtM_epsmach'] = float(np.sqrt(Mmax) * np.finfo(np.float64).eps)
    comp['eps_expm'] = LANCZOS_TOL
    comp['eps_H'] = float(max(r['t7_conservation_500step'][k]['dH_rel_max']
                              for r in all_res for k in ('dense', 'sparse')))
    comp['eps_Q2'] = float(max(r['t7_conservation_500step'][k]['dQ2_rel_max']
                               for r in all_res for k in ('dense', 'sparse')))
    comp['eps_equiv'] = float(max(max(r['t4_perm_equiv_max'], r['t5_u1_equiv_max'],
                                      r['t6_sign_equiv_max']) for r in all_res))
    eps_floor = max(comp.values())
    thresholds = {
        'eps_floor_components': comp,
        'eps_floor': eps_floor,
        'C_close': 1e2, 'C_move': 1e6, 'C_a': 1e3,
        'eps_close': 1e2 * eps_floor,
        'eps_move': min(1e6 * eps_floor, 1e-3),
        'a_min_over_normZ': 1e3 * eps_floor,
        'dtau0_subseries': '2*pi/40 (den=40)',
        'note': '採用後は本実験の結果を見て変更禁止（仕様書 v1.2 §10）',
    }
    out = {'phase0_audits': all_res, 'thresholds': thresholds}
    with open(os.path.join(RESULTS, 'phase0_results.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    print('\neps_floor components:', json.dumps(comp, indent=1))
    print(f"eps_floor={eps_floor:.3e}  eps_close={thresholds['eps_close']:.3e}  "
          f"eps_move={thresholds['eps_move']:.3e}  a_min/||Z||={thresholds['a_min_over_normZ']:.3e}")
    print('PHASE0 DONE')


if __name__ == '__main__':
    main()
