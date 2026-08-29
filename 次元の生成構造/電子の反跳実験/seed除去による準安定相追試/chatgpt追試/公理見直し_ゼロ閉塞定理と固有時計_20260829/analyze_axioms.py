#!/usr/bin/env python3
"""公理見直し：ゼロ閉塞定理・σ=N−1 のスケール不変形・τ と固有時計の乖離を、既存パッケージの保存データだけから再集計する。
新しい走行はしない（read-only）。入力は ../ 配下の 4 系統（original / fixed_baseline / fixed / fixed_equimodular）の step 0 と、
../N{5,8,10,16,20}_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828/ の時系列。"""
import csv, glob, gzip, os, sys, itertools
import numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); UP = os.path.dirname(HERE); OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
V1 = os.path.join(UP, '論文v1_全プログラム修正版_20260828'); ORIG = os.path.join(UP, '論文v1_全再現テスト_20260828', 'original')
SYSTEMS = [('original', ORIG), ('fixed_baseline', os.path.join(V1, 'fixed_baseline')), ('fixed', os.path.join(V1, 'fixed')), ('fixed_equimodular', os.path.join(V1, 'fixed_equimodular'))]
PK = {3:'N3_N4',4:'N3_N4',5:'N5',6:'N6_N7',7:'N6_N7',8:'N8_N9',9:'N8_N9',10:'N10_N11',11:'N10_N11',12:'N12_N13',13:'N12_N13',14:'N14_N15',15:'N14_N15',16:'N16'}
def edges(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def adjacency(N):
    E = edges(N); M = len(E); A = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            if i != j and set(E[i]) & set(E[j]): A[i, j] = 1.0
    return A
def step0(root, N):
    fs = sorted(glob.glob(os.path.join(root, f'{PK[N]}_complex_simplex_complete_analysis_20260826', f'N{N}_all_steps*')))
    if not fs: return None
    f = [x for x in fs if not x.endswith('.gz')] or fs; f = f[0]
    fh = gzip.open(f, 'rt') if f.endswith('.gz') else open(f)
    rd = csv.DictReader(fh)
    if 'edge_index' in rd.fieldnames:
        rows = []
        for x in rd:
            if x['step'] != '0': break
            rows.append(x)
        rows.sort(key=lambda x: int(x['edge_index']))
        return np.array([float(x['a']) + 1j * float(x['b']) for x in rows])
    x = next(rd); ac = [c for c in rd.fieldnames if c.endswith('_a')]; bc = [c[:-2] + '_b' for c in ac]
    return np.array([float(x[c]) for c in ac]) + 1j * np.array([float(x[c]) for c in bc])
# ---- 1. step 0 のゼロ閉塞（4 系統 × N）
w = csv.writer(open(os.path.join(OUT, 'closure_step0.csv'), 'w', newline=''))
w.writerow(['system', 'N', 'M', 'norm_v', 'r2_rel_spread', 'sum_a2', 'sum_b2', 'sum_ab', 'abs_sum_z2_over_H', 'a2_minus_b2_over_H', 'ab_over_H'])
print('## 1. step 0 ゼロ閉塞')
for name, root in SYSTEMS:
    for N in sorted(PK):
        z = step0(root, N)
        if z is None: continue
        a, b = z.real, z.imag; H = float((abs(z) ** 2).sum()); r2 = abs(z) ** 2
        A2, B2, AB = float((a * a).sum()), float((b * b).sum()), float((a * b).sum())
        row = [name, N, len(z), np.linalg.norm(z), (r2.max() - r2.min()) / r2.mean(), A2, B2, AB, abs((z * z).sum()) / H, (A2 - B2) / H, AB / H]
        w.writerow(row); print(f'{name:18s} N={N:2d} ‖v‖={row[3]:.6f} 幅={row[4]:.1e} Σa²={A2:.12f} Σb²={B2:.12f} Σab={AB:+.1e} |Σz²|/H={row[8]:.1e}')
# ---- 2. 走行中の保存（directHperp 時系列）＋ 位相進み
w2 = csv.writer(open(os.path.join(OUT, 'closure_conservation_and_phase_advance.csv'), 'w', newline=''))
w2.writerow(['N', 'step', 'Hperp_over_H', 'abs_ZT_Z', 'H_total', 'sigma1_of_step_generator', 'measured_phase_advance_rad'])
print('\n## 2. 走行中の |ZᵀZ| と τ 一歩あたりの位相進み')
summ = []
for N in [5, 8, 10, 16, 20]:
    f = os.path.join(UP, f'N{N}_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828', 'data', 'treatment_linear124_amplitude_aware_timeseries.csv')
    if not os.path.exists(f): continue
    r = list(csv.DictReader(open(f)))
    st = np.array([float(x['step']) for x in r]); zz = np.array([float(x['abs_ZT_Z']) for x in r]); H = np.array([float(x['H_total']) for x in r])
    fr = np.array([float(x['H_perp']) for x in r]) / H; s1 = np.array([float(x['sigma1_of_step_generator']) for x in r]); mp = np.array([float(x['measured_phase_advance']) for x in r])
    for t in [0, 1, 1000, 5000, 10000, 20000, 30000, 35000, 38000, 39999]:
        i = int(np.argmin(abs(st - t))); w2.writerow([N, int(st[i]), fr[i], zz[i], H[i], s1[i], mp[i]])
    ok = np.isfinite(mp); summ.append((N, zz.max(), int(st[zz.argmax()]), (H.max() - H.min()) / H.mean(), mp[ok].min(), mp[ok].max(), (mp[ok].max() - mp[ok].min()) / abs(mp[ok].mean()), fr.max()))
    print(f'N={N:2d}: |ZᵀZ| max={zz.max():.2e}(step {int(st[zz.argmax()])}) H_total 相対変動={summ[-1][3]:.1e}  位相進み min={mp[ok].min():.8f} max={mp[ok].max():.8f} 相対変動={summ[-1][6]:.1e}  H⊥/H max={fr.max():.2e}')
with open(os.path.join(OUT, 'phase_advance_summary.csv'), 'w', newline='') as fh:
    ww = csv.writer(fh); ww.writerow(['N', 'max_abs_ZT_Z', 'at_step', 'H_total_rel_var', 'phase_adv_min', 'phase_adv_max', 'phase_adv_rel_var', 'max_Hperp_over_H']); ww.writerows(summ)
# ---- 3. 等モジュラー親の σ スペクトル（スケール除去 σ/r²）
w3 = csv.writer(open(os.path.join(OUT, 'sigma_spectrum_equimodular_parent.csv'), 'w', newline=''))
w3.writerow(['N', 'M', 'r2', 'mu_over_r2', 'k', 'sigma_over_r2', 'ratio_to_sigma1', 'best_fraction_den_le_60', 'fraction_error'])
print('\n## 3. 等モジュラー親の K_amp スペクトル（σ/r²）と有理近似')
for N in range(3, 11):
    z = step0(os.path.join(V1, 'fixed_equimodular'), N)
    if z is None: continue
    A = adjacency(N); K = A * np.imag(np.conj(z)[:, None] * z[None, :]); r2 = float((abs(z) ** 2).mean())
    ev = np.linalg.eigvalsh(1j * K); s = np.sort(ev[ev > 1e-9])[::-1] / r2
    mu = float(np.vdot(z, 1j * (K @ z)).real / np.vdot(z, z).real) / r2
    line = []
    for k, x in enumerate(s):
        fr_ = Fraction(x / s[0]).limit_denominator(60); err = abs(x / s[0] - float(fr_))
        w3.writerow([N, len(z), r2, mu, k + 1, x, x / s[0], str(fr_), err]); line.append(f'{x:.5f}({fr_},{err:.0e})')
    print(f'N={N:2d} r²={r2:.6f} μ/r²={mu:+.6f} (N−1={N-1})  σ/r²: ' + ' '.join(line))
print('\nresults →', OUT)
