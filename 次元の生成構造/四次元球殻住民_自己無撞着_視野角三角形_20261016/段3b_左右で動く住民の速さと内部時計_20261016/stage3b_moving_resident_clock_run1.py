# -*- coding: utf-8 -*-
"""段 3b（README の事前固定のとおり）"""
import json, numpy as np
from s3_harmonics import *
N = 12; TH = 2*np.pi*23/124; STEPS = 124
f0 = [np.zeros((n+1, n+1), complex) for n in range(N+1)]
for n in range(N+1): f0[n][0, 0] = 1j/np.sqrt(n+1)          # f(q) = i Σ D_n(q)[0,0]
S = np.linspace(-np.pi, np.pi, 4001)
basis = None
def along_circle(C):
    global basis
    if basis is None:
        basis = [[np.sqrt(n+1)*Dn(V(qexp(s, 'i')), n) for n in range(N+1)] for s in S]
    return np.array([sum(np.sum(C[n]*Bs[n]) for n in range(N+1)) for Bs in basis])
def peak(C):
    v = along_circle(C); a = np.abs(v); k = int(np.argmax(a)); k = min(max(k, 1), len(S)-2)
    y0, y1, y2 = a[k-1], a[k], a[k+1]; d = 0.5*(y0-y2)/(y0-2*y1+y2) if (y0-2*y1+y2) != 0 else 0.0
    s = S[k] + d*(S[1]-S[0])
    val = sum(np.sqrt(n+1)*np.sum(C[n]*Dn(V(qexp(s, 'i')), n)) for n in range(N+1))
    return s, val
rows = []; B0 = true_square(f0)
for frac in (1.0, 0.75, 0.5, 0.25, 0.0):
    thL = frac*TH; a = qexp(-thL, 'i'); b = qexp(TH, 'i')
    C = f0; s_prev, v_prev = peak(C); disp = []; phase = []; bmax = 0.0
    for k in range(1, STEPS+1):
        C = rotate_right(translate_left(C, a), b)
        if k <= 8:                                   # 速さと位相の進みは最初の 8 ステップで読む（山が大円を回り込む前）
            s, val = peak(C)
            ds = (s - s_prev + np.pi) % (2*np.pi) - np.pi
            dph = (np.angle(val) - np.angle(v_prev) + np.pi) % (2*np.pi) - np.pi
            disp.append(ds); phase.append(dph); s_prev, v_prev = s, val
        bmax = max(bmax, abs(true_square(C)-B0)/abs(B0))
    rows.append(dict(thetaL_over_thetaR=frac, step_disp_mean=float(np.mean(disp)), step_disp_std=float(np.std(disp)),
                     step_phase_mean=float(np.mean(phase)), step_phase_std=float(np.std(phase)), B_max_rel_change_124=float(bmax)))
c_disp = [r for r in rows if r['thetaL_over_thetaR'] == 0.0][0]['step_disp_mean']
w_rest = [r for r in rows if r['thetaL_over_thetaR'] == 1.0][0]['step_phase_mean']
for r in rows:
    r['v_read'] = r['step_disp_mean']/c_disp
    r['omega_over_rest'] = r['step_phase_mean']/w_rest if w_rest != 0 else None
    r['sqrt_1_minus_v2'] = float(np.sqrt(max(0.0, 1 - r['v_read']**2)))
judge = dict(B_preserved=all(r['B_max_rel_change_124'] <= 1e-12 for r in rows),
             rest_does_not_move=abs([r for r in rows if r['thetaL_over_thetaR'] == 1.0][0]['step_disp_mean']) <= 1e-9,
             massless_moves_thetaR=abs(c_disp/TH - 1) <= 1e-6 or abs(c_disp/TH + 1) <= 1e-6)
out = dict(rows=rows, judge=judge, thetaR=TH, massless_step_disp=c_disp, rest_step_phase=w_rest)
json.dump(out, open('stage3b_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for r in rows: print(f"θL/θR={r['thetaL_over_thetaR']:.2f}  Δs={r['step_disp_mean']:+.6f}  v_read={r['v_read']:.4f}  ω={r['step_phase_mean']:+.6f}  ω/ω_rest={r['omega_over_rest']:.4f}  √(1−v²)={r['sqrt_1_minus_v2']:.4f}  B={r['B_max_rel_change_124']:.1e}")
print('θR =', TH)
