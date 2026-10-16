# -*- coding: utf-8 -*-
"""段 4b（README の事前固定のとおり。段 4a からの変更は交換を実数の角の回転にしたことだけ）"""
import json, numpy as np
from s3_harmonics import *
N = 12; STEPS = 248; s0 = 0.01
def resident_hw():
    return [np.array([[1j/np.sqrt(n+1) if (a == 0 and b == 0) else 0 for b in range(n+1)] for a in range(n+1)], complex) for n in range(N+1)]
rows = []
for name, th in (('theta_23pi_124', 23*np.pi/124), ('theta_R0.7', np.arccos(np.sqrt(0.7)))):
    c, s = np.cos(th), np.sin(th)
    for eta in (0.0, 0.3):
        base = translate_left(resident_hw(), qexp(1j*eta, 'i')); w = np.sqrt(norm_h(base))
        A = [np.sqrt(0.5+s0)*x/w for x in base]; Bm = [np.sqrt(0.5-s0)*x/w for x in base]
        A0, B0 = A, Bm; Z0 = true_square(A)+true_square(Bm); H0 = norm_h(A)+norm_h(Bm)
        shape0 = [x/np.sqrt(norm_h(A)) for x in A]
        S = []; zmax = hmax = sdev = 0.0
        for k in range(1, STEPS+1):
            A, Bm = [c*a - s*b for a, b in zip(A, Bm)], [s*a + c*b for a, b in zip(A, Bm)]
            a2, b2 = norm_h(A), norm_h(Bm); S.append((a2-b2)/(a2+b2))
            zmax = max(zmax, abs(true_square(A)+true_square(Bm)-Z0)/abs(Z0)); hmax = max(hmax, abs(a2+b2-H0)/H0)
            ratio = sum(np.vdot(sh, x) for sh, x in zip(shape0, A))
            sdev = max(sdev, max(np.max(np.abs(x - ratio*sh)) for sh, x in zip(shape0, A))/max(np.sqrt(a2), 1e-300))
            if k == 124: A124, B124 = A, Bm
        err124 = max(max(np.max(np.abs(x+y)) for x, y in zip(A124, A0)), max(np.max(np.abs(x+y)) for x, y in zip(B124, B0)))
        err248 = max(max(np.max(np.abs(x-y)) for x, y in zip(A, A0)), max(np.max(np.abs(x-y)) for x, y in zip(Bm, B0)))
        a0, b0 = np.sqrt(0.5+s0), np.sqrt(0.5-s0)
        closed = [((a0*np.cos(k*th)-b0*np.sin(k*th))**2 - (a0*np.sin(k*th)+b0*np.cos(k*th))**2) for k in range(1, STEPS+1)]
        rows.append(dict(case=name, eta=eta, sumz2_max_rel_change=float(zmax), norm_max_rel_change=float(hmax),
                         err_124_minus_identity=float(err124), err_248_identity=float(err248),
                         S_vs_closed_max=float(max(abs(a-b) for a, b in zip(S, closed))), S_min=float(min(S)), S_max=float(max(S)),
                         A_shape_max_dev=float(sdev)))
judge = dict(sumz2_preserved=all(r['sumz2_max_rel_change'] <= 1e-12 for r in rows), norm_preserved=all(r['norm_max_rel_change'] <= 1e-12 for r in rows),
             minus_identity_at_124=all(r['err_124_minus_identity'] <= 1e-12 for r in rows if r['case'] == 'theta_23pi_124'),
             identity_at_248=all(r['err_248_identity'] <= 1e-12 for r in rows if r['case'] == 'theta_23pi_124'))
json.dump(dict(rows=rows, judge=judge), open('stage4b_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for r in rows: print({k: (float('%.6g' % v) if isinstance(v, float) else v) for k, v in r.items()})
