# -*- coding: utf-8 -*-
"""段 3d（README の事前固定のとおり）"""
import json, numpy as np
def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2, a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2])
def qexp(phi, u):
    vec = {'i': (1,0,0), 'j': (0,1,0), 'k': (0,0,1)}[u]
    return np.array([np.cos(phi)] + [np.sin(phi)*c for c in vec], dtype=complex)
UIDX = {'i': 1, 'j': 2, 'k': 3}
TH = 2*np.pi*23/124; DELTA = 1.0; K = 20
rows = []
for u in 'ijk':
    for eta in (0.0, 0.1, 0.3, 0.6, 1.0):
        Lam = qexp(1j*eta, u)
        A = [np.array([1j*k*DELTA, 0, 0, 0], dtype=complex) for k in range(K)]
        Ap = [qmul(Lam, a) for a in A]
        inv_err = max(abs(ap @ ap - a @ a)/max(abs(a @ a), 1e-300) for a, ap in zip(A[1:], Ap[1:]))
        t = np.array([ap[0].imag for ap in Ap]); x = np.array([ap[UIDX[u]].real for ap in Ap])
        dt = np.diff(t); dx = np.diff(x)
        omega_ratio = (TH/dt.mean())/(TH/DELTA); v = (dx/dt).mean()
        # 右の時計の回転（内部の位相、右から掛ける単位四元数）と左のブーストの可換
        b = qexp(TH, 'k'); comm = max(np.max(np.abs(qmul(qmul(Lam, a), b) - qmul(Lam, qmul(a, b)))) for a in A)
        rows.append(dict(u=u, eta=eta, invariant_rel_err=float(inv_err), commute_err=float(comm),
                         omega_ratio=float(omega_ratio), inv_cosh=float(1/np.cosh(eta)),
                         v=float(v), tanh=float(np.tanh(eta)), sqrt_1_v2=float(np.sqrt(1 - v*v))))
judge = dict(invariant_preserved=all(r['invariant_rel_err'] <= 1e-12 for r in rows), boost_commutes_with_right_clock=all(r['commute_err'] <= 1e-12 for r in rows))
json.dump(dict(rows=rows, judge=judge), open('stage3d_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for r in rows:
    if r['u'] == 'i': print(f"η={r['eta']}: ω'/ω0={r['omega_ratio']:.6f}  1/coshη={r['inv_cosh']:.6f}  v'={r['v']:+.6f}  tanhη={r['tanh']:.6f}  √(1−v'²)={r['sqrt_1_v2']:.6f}")
print('max |ω\'/ω0 − 1/coshη| over all:', max(abs(r['omega_ratio']-r['inv_cosh']) for r in rows), '| max ||v\'|−tanhη|:', max(abs(abs(r['v'])-r['tanh']) for r in rows))
