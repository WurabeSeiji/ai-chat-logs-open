#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Control experiment: internalize the external symmetric mass-ratio nu as
an explicit sixth persistent multiplicative-scale state N.

Scope is deliberately narrow:
  * Legacy 5-state Method-B dynamics are the immutable control.
  * Only the location of nu information is changed.
  * Legacy path: micro_fast(s, nu)
  * New path:    micro_fast6(s), reading N from the OLD state.
  * N itself is propagated by the same one-hot product-sum selector with
    identity coefficient 1 in every phase: N' = sum_j q_j (1*N) = N.
  * No charged interaction is introduced in this experiment.

The broader pre-existing additive RK4 stage expressions are intentionally not
changed here. This experiment tests only whether nu can be internalized without
changing the previous experiment's physical/numerical results.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import importlib.util
import inspect
import json
import math
from pathlib import Path
from typing import Dict, Any

import numpy as np

HERE = Path(__file__).resolve().parent
REF_BASE = HERE / "reference_audit_two_paths_rk4_snapshot.py"
REF_B = HERE / "reference_audit_method_B_strict_snapshot.py"
SAVED_C1 = HERE / "experiment02_sn_2p5pn_rr_C1_N12.csv"

# Load immutable snapshots.
spec_base = importlib.util.spec_from_file_location("base_ref", REF_BASE)
base = importlib.util.module_from_spec(spec_base)
spec_base.loader.exec_module(base)

# The historical B snapshot imports /mnt/data/audit_two_paths_rk4.py. Instead of
# executing that module here, reproduce its exact micro_fast code locally below
# as legacy_micro_fast to avoid path-dependent behavior.
H, PH, PH2, TAU0 = base.H, base.PH, base.PH2, base.TAU0
rr, gN, ref_step = base.rr, base.gN, base.ref_step
init_B, coreB = base.init_B, base.coreB

IDENTITY_N_PHASE_COEFFICIENTS = np.ones(11, dtype=float)


def clone_state(s: Dict[str, Any]) -> Dict[str, Any]:
    return {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in s.items()}


def blend(q, cands):
    """Exact historical one-hot product-sum selector."""
    x = cands[0] * q[0]
    for j in range(1, len(cands)):
        x = x + cands[j] * q[j]
    return x


def legacy_micro_fast(s, nu):
    """Byte-for-byte algebraic transcription of historical audit_method_B_strict.py."""
    j = int(np.argmax(s['q']))
    o = clone_state(s)
    n = clone_state(s)
    U, P, E = o['U'], o['P'], o['E']
    if j == 0:
        n['k1'] = rr(U, P, E, nu)
    elif j == 1:
        n['Us'] = U * PH2; n['Ps'] = P + 0.5 * H * o['k1'][0]; n['Es'] = E + 0.5 * H * o['k1'][1]
    elif j == 2:
        n['k2'] = rr(o['Us'], o['Ps'], o['Es'], nu)
    elif j == 3:
        n['Us'] = U * PH2; n['Ps'] = P + 0.5 * H * o['k2'][0]; n['Es'] = E + 0.5 * H * o['k2'][1]
    elif j == 4:
        n['k3'] = rr(o['Us'], o['Ps'], o['Es'], nu)
    elif j == 5:
        n['Us'] = U * PH; n['Ps'] = P + H * o['k3'][0]; n['Es'] = E + H * o['k3'][1]
    elif j == 6:
        n['k4'] = rr(o['Us'], o['Ps'], o['Es'], nu)
    elif j == 7:
        n['d'] = (H / 6) * (o['k1'] + 2 * o['k2'] + 2 * o['k3'] + o['k4'])
    elif j == 8:
        n['Um'] = U * PH2; n['Pm'] = P + 0.5 * o['d'][0]; n['Em'] = E + 0.5 * o['d'][1]
    elif j == 9:
        n['dphi'] = H * gN(o['Um'], o['Pm'], o['Em'])
    elif j == 10:
        n['U'] = U * PH; n['P'] = P + o['d'][0]; n['E'] = E + o['d'][1]
        n['H'] = o['H'] * complex(math.cos(o['dphi']), math.sin(o['dphi']))
        n['Q'] = o['Q'] * math.exp(o['d'][2] / TAU0)
        n['k1'] = np.zeros(3); n['k2'] = np.zeros(3); n['k3'] = np.zeros(3); n['k4'] = np.zeros(3)
        n['Us'] = U * PH; n['Ps'] = P + o['d'][0]; n['Es'] = E + o['d'][1]; n['d'] = np.zeros(3)
        n['Um'] = U * PH; n['Pm'] = P + o['d'][0]; n['Em'] = E + o['d'][1]; n['dphi'] = 0.0
    n['q'] = np.roll(o['q'], 1)
    return n


def legacy_macro_fast(s, nu):
    for _ in range(11):
        s = legacy_micro_fast(s, nu)
    return s


def init_B6(U, P, E, Hphi, Q, nu0):
    """Initialize the six-state system.

    Mathematical parameterization: N0 = exp(log(nu0)) = nu0 for nu0>0.
    Machine representation stores the positive group element directly as the
    exact float nu0, deliberately avoiding an unnecessary transcendental
    exp(log(.)) round-trip that could alter low bits.
    """
    if not (nu0 > 0.0):
        raise ValueError("nu0 must be positive")
    s = init_B(U, P, E, Hphi, Q)
    s['N'] = float(nu0)
    return s


def micro_fast6(s):
    """New 5+1 microstep. No external nu argument exists.

    RR reads N from the OLD state. N is written by the same one-hot product-sum
    selector used to represent the phase-conditioned interaction:
        N_next = sum_j q_j * (c_j * N_old),  c_j = 1 for all 11 phases.
    Thus the sixth row is an identity row, but it is not an external copy.
    """
    j = int(np.argmax(s['q']))
    o = clone_state(s)
    n = clone_state(s)
    U, P, E, N = o['U'], o['P'], o['E'], o['N']

    if j == 0:
        n['k1'] = rr(U, P, E, N)
    elif j == 1:
        n['Us'] = U * PH2; n['Ps'] = P + 0.5 * H * o['k1'][0]; n['Es'] = E + 0.5 * H * o['k1'][1]
    elif j == 2:
        n['k2'] = rr(o['Us'], o['Ps'], o['Es'], N)
    elif j == 3:
        n['Us'] = U * PH2; n['Ps'] = P + 0.5 * H * o['k2'][0]; n['Es'] = E + 0.5 * H * o['k2'][1]
    elif j == 4:
        n['k3'] = rr(o['Us'], o['Ps'], o['Es'], N)
    elif j == 5:
        n['Us'] = U * PH; n['Ps'] = P + H * o['k3'][0]; n['Es'] = E + H * o['k3'][1]
    elif j == 6:
        n['k4'] = rr(o['Us'], o['Ps'], o['Es'], N)
    elif j == 7:
        n['d'] = (H / 6) * (o['k1'] + 2 * o['k2'] + 2 * o['k3'] + o['k4'])
    elif j == 8:
        n['Um'] = U * PH2; n['Pm'] = P + 0.5 * o['d'][0]; n['Em'] = E + 0.5 * o['d'][1]
    elif j == 9:
        n['dphi'] = H * gN(o['Um'], o['Pm'], o['Em'])
    elif j == 10:
        n['U'] = U * PH; n['P'] = P + o['d'][0]; n['E'] = E + o['d'][1]
        n['H'] = o['H'] * complex(math.cos(o['dphi']), math.sin(o['dphi']))
        n['Q'] = o['Q'] * math.exp(o['d'][2] / TAU0)
        n['k1'] = np.zeros(3); n['k2'] = np.zeros(3); n['k3'] = np.zeros(3); n['k4'] = np.zeros(3)
        n['Us'] = U * PH; n['Ps'] = P + o['d'][0]; n['Es'] = E + o['d'][1]; n['d'] = np.zeros(3)
        n['Um'] = U * PH; n['Pm'] = P + o['d'][0]; n['Em'] = E + o['d'][1]; n['dphi'] = 0.0

    # Sixth persistent state through the interaction selector itself.
    Nc = [IDENTITY_N_PHASE_COEFFICIENTS[k] * N for k in range(11)]
    n['N'] = float(blend(o['q'], Nc))
    n['q'] = np.roll(o['q'], 1)
    return n


def macro_fast6(s):
    for _ in range(11):
        s = micro_fast6(s)
    return s


def core6(s):
    return np.array([s['U'], s['P'], s['E'], s['H'], s['Q'], s['N']], dtype=np.complex128)


def pack_core5_bytes(z):
    return np.asarray(z, dtype=np.complex128).tobytes()


def float_bits_equal(a, b):
    return np.asarray([a], dtype=np.float64).tobytes() == np.asarray([b], dtype=np.float64).tobytes()


def readout_core(z):
    U, P, E, Hh, Q = z
    ce = base.trigU(U)[0]
    r = float(P.real) / (1.0 + float(E.real) * ce)
    x = r * float(Hh.real)
    y = r * float(Hh.imag)
    t = TAU0 * math.log(float(Q.real))
    return r, x, y, t


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def static_audit():
    src_micro = inspect.getsource(micro_fast6)
    src_macro = inspect.getsource(macro_fast6)
    tree_micro = ast.parse(src_micro)
    tree_macro = ast.parse(src_macro)
    micro_args = [a.arg for a in tree_micro.body[0].args.args]
    macro_args = [a.arg for a in tree_macro.body[0].args.args]
    names = {n.id for n in ast.walk(tree_micro) if isinstance(n, ast.Name)}
    # 'nu' must not be a runtime name in the new transition function.
    return {
        'micro_fast6_args': micro_args,
        'macro_fast6_args': macro_args,
        'runtime_name_nu_present_in_micro_fast6': 'nu' in names,
        'identity_N_phase_coefficients': IDENTITY_N_PHASE_COEFFICIENTS.tolist(),
        'all_identity_coefficients_exactly_one': bool(np.array_equal(IDENTITY_N_PHASE_COEFFICIENTS, np.ones(11))),
    }


def run_case(case_id, p0, e0, qmass, nsteps, raw_writer=None, sample_every=1):
    nu0 = qmass / (1.0 + qmass) ** 2
    x0 = np.array([1+0j, p0, e0, 1+0j, 1+0j], dtype=np.complex128)
    old = init_B(*x0)
    new = init_B6(*x0, nu0)

    max_core_abs = np.zeros(5, dtype=float)
    max_readout_abs = np.zeros(4, dtype=float)
    bitwise_core_mismatch = 0
    N_bitwise_mismatch = 0
    N_micro_hold_fail = 0
    phase_fail = 0

    def write_row(k):
        nonlocal bitwise_core_mismatch, N_bitwise_mismatch
        zo = coreB(old)
        zn = core6(new)[:5]
        ro = readout_core(zo)
        rn = readout_core(zn)
        core_eq = pack_core5_bytes(zo) == pack_core5_bytes(zn)
        n_eq = float_bits_equal(new['N'], nu0)
        if not core_eq:
            bitwise_core_mismatch += 1
        if not n_eq:
            N_bitwise_mismatch += 1
        dcore = np.abs(zn - zo)
        dread = np.abs(np.array(rn) - np.array(ro))
        max_core_abs[:] = np.maximum(max_core_abs, dcore)
        max_readout_abs[:] = np.maximum(max_readout_abs, dread)
        if raw_writer is not None and k % sample_every == 0:
            raw_writer.writerow([
                case_id, k, f"{p0:.17g}", f"{e0:.17g}", f"{qmass:.17g}", f"{nu0:.17g}", f"{new['N']:.17g}",
                int(n_eq), int(core_eq),
                f"{zo[0].real:.17g}", f"{zo[0].imag:.17g}", f"{zn[0].real:.17g}", f"{zn[0].imag:.17g}", f"{dcore[0]:.17g}",
                f"{zo[1].real:.17g}", f"{zn[1].real:.17g}", f"{dcore[1]:.17g}",
                f"{zo[2].real:.17g}", f"{zn[2].real:.17g}", f"{dcore[2]:.17g}",
                f"{zo[3].real:.17g}", f"{zo[3].imag:.17g}", f"{zn[3].real:.17g}", f"{zn[3].imag:.17g}", f"{dcore[3]:.17g}",
                f"{zo[4].real:.17g}", f"{zn[4].real:.17g}", f"{dcore[4]:.17g}",
                f"{ro[0]:.17g}", f"{rn[0]:.17g}", f"{dread[0]:.17g}",
                f"{ro[1]:.17g}", f"{rn[1]:.17g}", f"{dread[1]:.17g}",
                f"{ro[2]:.17g}", f"{rn[2]:.17g}", f"{dread[2]:.17g}",
                f"{ro[3]:.17g}", f"{rn[3]:.17g}", f"{dread[3]:.17g}",
            ])

    write_row(0)
    for k in range(1, nsteps + 1):
        old = legacy_macro_fast(old, nu0)
        # Inspect all 11 new microsteps so N identity is tested at the actual internal resolution.
        for _ in range(11):
            beforeN = new['N']
            q_expected = np.roll(new['q'], 1)
            new2 = micro_fast6(new)
            if not float_bits_equal(new2['N'], beforeN):
                N_micro_hold_fail += 1
            if not (np.array_equal(new2['q'], q_expected) and np.count_nonzero(new2['q'] == 1.0) == 1 and np.sum(new2['q']) == 1.0):
                phase_fail += 1
            new = new2
        write_row(k)

    return {
        'case_id': case_id,
        'p0': p0, 'e0': e0, 'qmass': qmass, 'nu0': nu0, 'nsteps': nsteps,
        'max_abs_U': float(max_core_abs[0]),
        'max_abs_P': float(max_core_abs[1]),
        'max_abs_E': float(max_core_abs[2]),
        'max_abs_H': float(max_core_abs[3]),
        'max_abs_Q': float(max_core_abs[4]),
        'max_abs_r': float(max_readout_abs[0]),
        'max_abs_x': float(max_readout_abs[1]),
        'max_abs_y': float(max_readout_abs[2]),
        'max_abs_t': float(max_readout_abs[3]),
        'bitwise_core_mismatch_count': int(bitwise_core_mismatch),
        'N_bitwise_mismatch_count': int(N_bitwise_mismatch),
        'N_micro_hold_fail_count': int(N_micro_hold_fail),
        'phase_fail_count': int(phase_fail),
        'final_N': float(new['N']),
        'final_N_minus_nu0': float(new['N'] - nu0),
    }


def compare_saved_c1_csv():
    data = np.genfromtxt(SAVED_C1, delimiter=',', names=True)
    p0 = float(data['p'][0]); e0 = float(data['e'][0]); nu0 = 0.25
    x0 = np.array([1+0j, p0, e0, 1+0j, 1+0j], dtype=np.complex128)
    old = init_B(*x0)
    new = init_B6(*x0, nu0)
    fields = ['p', 'e', 't', 'x', 'y']
    err_old = {k: 0.0 for k in fields}
    err_new = {k: 0.0 for k in fields}
    old_new_bitwise_mismatch = 0
    for i, row in enumerate(data):
        zo = coreB(old); zn = core6(new)[:5]
        if pack_core5_bytes(zo) != pack_core5_bytes(zn):
            old_new_bitwise_mismatch += 1
        for label, z, err in [('old', zo, err_old), ('new', zn, err_new)]:
            r, x, y, t = readout_core(z)
            err['p'] = max(err['p'], abs(float(z[1].real) - row['p']))
            err['e'] = max(err['e'], abs(float(z[2].real) - row['e']))
            err['t'] = max(err['t'], abs(t - row['t_osc']))
            err['x'] = max(err['x'], abs(x - row['x']))
            err['y'] = max(err['y'], abs(y - row['y']))
        if i + 1 < len(data):
            old = legacy_macro_fast(old, nu0)
            new = macro_fast6(new)
    return {
        'rows': int(len(data)), 'p0': p0, 'e0': e0, 'nu0': nu0,
        'legacy_max_abs_vs_saved_csv': err_old,
        'new6_max_abs_vs_saved_csv': err_new,
        'old_new_bitwise_mismatch_count': int(old_new_bitwise_mismatch),
        'error_dicts_exactly_equal': err_old == err_new,
    }


def main():
    outdir = HERE
    raw_path = outdir / 'nu_internal_state_control_timeseries.csv'
    cases_path = outdir / 'nu_internal_state_control_case_summary.csv'

    cases = []
    for p0 in (24.3, 60.0, 120.0):
        for e0 in (0.20, 0.45, 0.55):
            cases.append((f"p{p0:g}_e{e0:g}_q1", p0, e0, 1.0))
    for qmass in (2.0, 4.0, 10.0):
        cases.append((f"p60_e0.45_q{qmass:g}", 60.0, 0.45, qmass))

    header = [
        'case_id','step','p0','e0','qmass','nu0','N_state','N_bitwise_equal_nu0','core5_bitwise_equal',
        'U_old_re','U_old_im','U_new_re','U_new_im','abs_dU',
        'P_old','P_new','abs_dP','E_old','E_new','abs_dE',
        'H_old_re','H_old_im','H_new_re','H_new_im','abs_dH',
        'Q_old','Q_new','abs_dQ',
        'r_old','r_new','abs_dr','x_old','x_new','abs_dx','y_old','y_new','abs_dy','t_old','t_new','abs_dt'
    ]
    summaries = []
    with raw_path.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for case in cases:
            summaries.append(run_case(*case, nsteps=4000, raw_writer=w))

    # Representative long-run audit mirrors the historical 3 x 12000 check.
    long_cases = [
        ('C1_exact_3orbit', 24.3, 0.4358898943540673, 1.0),
        ('p60_e0.45_q4_3orbit', 60.0, 0.45, 4.0),
        ('p120_e0.55_q1_3orbit', 120.0, 0.55, 1.0),
    ]
    long_summaries = [run_case(*c, nsteps=12000, raw_writer=None) for c in long_cases]

    with cases_path.open('w', encoding='utf-8', newline='') as f:
        cols = list(summaries[0].keys())
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(summaries)

    static = static_audit()
    saved_c1 = compare_saved_c1_csv()
    all_summary = summaries + long_summaries
    global_max = {
        k: max(row[k] for row in all_summary)
        for k in ['max_abs_U','max_abs_P','max_abs_E','max_abs_H','max_abs_Q','max_abs_r','max_abs_x','max_abs_y','max_abs_t']
    }
    totals = {
        'bitwise_core_mismatch_count': sum(row['bitwise_core_mismatch_count'] for row in all_summary),
        'N_bitwise_mismatch_count': sum(row['N_bitwise_mismatch_count'] for row in all_summary),
        'N_micro_hold_fail_count': sum(row['N_micro_hold_fail_count'] for row in all_summary),
        'phase_fail_count': sum(row['phase_fail_count'] for row in all_summary),
    }

    result = {
        'experiment': 'nu internal sixth-state control v1',
        'scope': 'only external nu -> internal multiplicative-scale state N; charged interactions not introduced',
        'state_vector_persistent': ['U','P','E','H','Q','N'],
        'N_definition_math': 'N0 = exp(log(nu0)) = nu0, nu0>0',
        'N_machine_initialization': 'store exact float nu0 directly to avoid exp(log(.)) low-bit perturbation',
        'N_interaction_identity': 'N_next = sum_j q_j * (1 * N_old), j=0..10',
        'static_audit': static,
        'cases_12x4000': summaries,
        'representative_3x12000': long_summaries,
        'global_max_abs_old5_vs_new6_first5': global_max,
        'global_failure_counts': totals,
        'saved_C1_csv_check': saved_c1,
        'reference_sha256': {
            REF_BASE.name: sha256(REF_BASE),
            REF_B.name: sha256(REF_B),
            SAVED_C1.name: sha256(SAVED_C1),
        },
    }
    summary_path = outdir / 'nu_internal_state_control_summary.json'
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    coeff_path = outdir / 'nu_identity_interaction_coefficients.json'
    coeff_path.write_text(json.dumps({
        'phase_count': 11,
        'N_row_phase_coefficients': IDENTITY_N_PHASE_COEFFICIENTS.tolist(),
        'cross_couplings_into_N_row': {'U':0,'P':0,'E':0,'H':0,'Q':0},
        'formula': 'N_next = sum_j q_j * (1*N_old) = N_old'
    }, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    print(json.dumps({
        'global_max_abs': global_max,
        'failures': totals,
        'saved_C1': saved_c1,
        'static': static,
        'raw_rows': 1 + sum(4000+1 for _ in cases),
    }, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
