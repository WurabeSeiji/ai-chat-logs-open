# -*- coding: utf-8 -*-
"""段 3：局在した住民と観測者（README の事前固定のとおり）"""
import json, numpy as np
from s3_harmonics import *
N = 12; A = 1.0; TH = 2*np.pi*23/124
def resident(q0):
    Vi = V(qinv(q0)); return [1j*A*Dn(Vi, n).T/np.sqrt(n+1) for n in range(N+1)]   # f(q)=iA Σ χ_n(q0⁻¹ q)
def rel(Ca, Cb):
    return max(np.max(np.abs(a-b)) for a, b in zip(Ca, Cb))/max(np.max(np.abs(b)) for b in Cb)
def Bpair(Cf, Cg, n=None):
    rng_ = range(len(Cf)) if n is None else [n]
    return sum(np.trace(Cf[k].T @ np.linalg.inv(J(k)) @ Cg[k] @ J(k)) for k in rng_)
out = {}
# 住民の確認：中心で山
q_p = np.array([1., 0, 0, 0]); q_o = qexp(0.9, 'j').real
f, g = resident(q_p), resident(q_o)
out['peak_check'] = dict(f_at_center=complex(field(f, q_p)).__repr__(), f_at_observer=complex(field(f, q_o)).__repr__(), g_at_observer=complex(field(g, q_o)).__repr__())
# (a) 並進：右の交換の回転を毎ステップ。予測は中心 q0 b^{-k}
b = qexp(TH, 'i'); C = f; worst = 0.0
for k in range(1, 249):
    C = rotate_right(C, b)
    pred_center = qmul(q_p, qinv(np.array(qexp(k*TH, 'i'))))
    worst = max(worst, rel(C, resident(pred_center.real)))
out['a_translation_max_rel_diff'] = float(worst)
# 時計：観測者の山に同じ右の回転、124 ステップで戻る
Cg = g
for k in range(124): Cg = rotate_right(Cg, b)
out['clock_observer_return_124'] = float(rel(Cg, g))
# (c) 観測者による読み出し：粒子に一度だけブースト
rows = []; B_f0 = true_square(f)
for u in 'ijk':
    X0 = Bpair(f, g); X0n = [Bpair(f, g, n) for n in range(N+1)]
    for eta in (0.0, 0.05, 0.1, 0.2, 0.3):
        Lf = translate_left(f, qexp(1j*eta, u))
        X = Bpair(Lf, g); Xn = [Bpair(Lf, g, n) for n in range(N+1)]
        rows.append(dict(u=u, eta=eta, B_self_rel_change=float(abs(true_square(Lf)-B_f0)/abs(B_f0)),
                         X_ratio=[float((X/X0).real), float((X/X0).imag)], cosh_eta=float(np.cosh(eta)),
                         Xn_ratio=[[float((a/b_).real), float((a/b_).imag)] if abs(b_) > 1e-12 else None for a, b_ in zip(Xn, X0n)]))
out['c_rows'] = rows
out['judge'] = dict(a_translation=out['a_translation_max_rel_diff'] <= 1e-12,
                    c_self_B_invariant=all(r['B_self_rel_change'] <= 1e-12 for r in rows),
                    clock_observer_return=out['clock_observer_return_124'] <= 1e-12)
json.dump(out, open('stage3_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(out['judge'], indent=1)); print('a', out['a_translation_max_rel_diff'], 'clock', out['clock_observer_return_124'])
print('peak', out['peak_check'])
for r in rows: print(r['u'], r['eta'], 'X/X0 =', np.round(r['X_ratio'], 6), 'cosh =', round(r['cosh_eta'], 6), 'B_self', f"{r['B_self_rel_change']:.1e}")
