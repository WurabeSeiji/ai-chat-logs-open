# -*- coding: utf-8 -*-
"""段 4a（README の事前固定のとおり）"""
import json, numpy as np
from s3_harmonics import *
N = 12; STEPS = 248; s0 = 0.01
def resident_hw():
    return [np.array([[1j/np.sqrt(n+1) if (a == 0 and b == 0) else 0 for b in range(n+1)] for a in range(n+1)], complex) for n in range(N+1)]
def coeffs(R):
    th = np.arcsin(np.sqrt(R)); ph = np.exp(1j*th)
    return ph*np.cos(th), -1j*ph*np.sin(th)          # t, r（論文 §2.1）
def nrm(C): return norm_h(C)
rows = []
for Rname, R in (('R_124_23', np.cos(23*np.pi/124)**2), ('R_0.7', 0.7)):
    t, r = coeffs(R)
    for eta in (0.0, 0.3):
        base = translate_left(resident_hw(), qexp(1j*eta, 'i')); w = np.sqrt(nrm(base))
        A = [np.sqrt(0.5+s0)*c/w for c in base]; Bm = [np.sqrt(0.5-s0)*c/w for c in base]
        A0, B0 = A, Bm; H0 = nrm(A)+nrm(Bm); shape0 = [c/np.sqrt(nrm(A)) for c in A]
        S = []; Z = []; shape_dev = 0.0; hmax = 0.0
        for k in range(STEPS):
            a2, b2 = nrm(A), nrm(Bm); S.append((a2-b2)/(a2+b2)); Z.append(true_square(A)+true_square(Bm))
            # A の形：A の係数を、初期の形の何倍か（複素の比）で割ったときのずれ
            ratio = sum(np.vdot(s, c) for s, c in zip(shape0, A))
            shape_dev = max(shape_dev, max(np.max(np.abs(c - ratio*s)) for s, c in zip(shape0, A))/max(np.sqrt(a2), 1e-300))
            hmax = max(hmax, abs(a2+b2-H0)/H0)
            A, Bm = [r*a + t*b for a, b in zip(A, Bm)], [t*a + r*b for a, b in zip(A, Bm)]
            if k == 123: A124, B124 = A, Bm
        ret = max(max(np.max(np.abs(x-y)) for x, y in zip(A124, A0)), max(np.max(np.abs(x-y)) for x, y in zip(B124, B0)))
        closed = [0.02*np.cos(k*(np.pi + 2*np.arcsin(np.sqrt(R)))) for k in range(STEPS)]
        rows.append(dict(R=Rname, eta=eta, norm_max_rel_change=float(hmax), return_err_124=float(ret),
                         S_vs_closed_form_max=float(max(abs(a-b) for a, b in zip(S, closed))),
                         sumz2_first=[float(np.real(Z[0])), float(np.imag(Z[0]))], sumz2_abs_min=float(min(abs(z) for z in Z)), sumz2_abs_max=float(max(abs(z) for z in Z)),
                         sumz2_rel_spread=float((max(abs(z) for z in Z)-min(abs(z) for z in Z))/max(abs(z) for z in Z)),
                         A_shape_max_dev=float(shape_dev)))
judge = dict(norm_preserved=all(x['norm_max_rel_change'] <= 1e-12 for x in rows),
             returns_after_124_for_R124_23=all(x['return_err_124'] <= 1e-12 for x in rows if x['R'] == 'R_124_23'))
json.dump(dict(rows=rows, judge=judge), open('stage4a_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for x in rows: print({k: (round(v, 15) if isinstance(v, float) else v) for k, v in x.items()})
