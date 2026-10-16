# -*- coding: utf-8 -*-
"""段 3c（README の事前固定のとおり）"""
import json, numpy as np
from s3_harmonics import *
N = 12; STEPS = 248
def resident_hw(q0):
    Vi = V(qinv(q0)); C = []
    for n in range(N+1):
        M = np.zeros((n+1, n+1), complex); D = Dn(Vi, n)
        M[:, 0] = 1j*D[0, :]/np.sqrt(n+1)        # f(q) = i Σ D_n(q0⁻¹ q)[0,0]
        C.append(M)
    return C
def Bpair(Cf, Cg):
    return sum(np.trace(Cf[k].T @ np.linalg.inv(J(k)) @ Cg[k] @ J(k)) for k in range(len(Cf)))
f0 = resident_hw(np.array([1., 0, 0, 0])); g0 = resident_hw(qexp(0.9, 'j').real)
out = dict(cases=[])
for case, (tp, to) in (('A_same', (2*np.pi*23/124, 2*np.pi*23/124)), ('B_different', (2*np.pi*23/124, 2*np.pi*29/156))):
    bp, bo = qexp(tp, 'k'), qexp(to, 'k')
    for u in 'ijk':
        for eta in (0.0, 0.1, 0.2, 0.3):
            Lam = qexp(1j*eta, u)
            f = translate_left(f0, Lam); g = g0; Bf0 = true_square(f)
            X = []; bmax = 0.0; f_late = f0
            for k in range(STEPS):
                X.append(Bpair(f, g))
                bmax = max(bmax, abs(true_square(f)-Bf0)/abs(Bf0))
                f = rotate_right(f, bp); g = rotate_right(g, bo); f_late = rotate_right(f_late, bp)
            # 可換の検査：ブーストを最後に掛けた場合の最終の交差項
            X_late = Bpair(translate_left(f_late, Lam), g)
            X = np.array(X)
            spec = np.fft.fft(X)/STEPS; top = np.argsort(-np.abs(spec))[:4]
            out['cases'].append(dict(case=case, u=u, eta=eta, B_self_max_rel=float(bmax),
                commute_rel=float(abs(X_late - Bpair(f, g))/max(abs(X_late), 1e-300)),
                X_std_over_mean=float(np.std(X)/max(abs(np.mean(X)), 1e-300)),
                top_freqs_cycles_per_step=[float(t/STEPS) for t in top], top_amps=[float(abs(spec[t])) for t in top]))
judge = dict(B_self_preserved=all(c['B_self_max_rel'] <= 1e-12 for c in out['cases']),
             boost_commutes_with_right=all(c['commute_rel'] <= 1e-10 for c in out['cases']))
out['judge'] = judge
json.dump(out, open('stage3c_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for c in out['cases']:
    if c['u'] in 'ij':
        print(c['case'], c['u'], c['eta'], 'X var/mean %.2e' % c['X_std_over_mean'], 'top freqs', np.round(c['top_freqs_cycles_per_step'], 5), 'amps', np.round(c['top_amps'], 4))
