# -*- coding: utf-8 -*-
"""段 5b（README の事前固定のとおり）"""
import json, numpy as np
src = open('stage5a_geometry_copy.py').read().split("TH = 23*np.pi/124")[0]
exec(src)                                   # P, edges, tets, d_geo, dihedrals(), curvature() を段 5a と同じ定義で得る
TH = 23*np.pi/124; STEPS = 248; s0 = 0.01; sig = 0.6
d0 = d_geo[0]; prof = np.exp(-d0**2/sig**2); far = int(np.argmax(d0))
Lap = np.zeros((NV, NV))
for i, j in edges: Lap[i, j] -= 1; Lap[j, i] -= 1; Lap[i, i] += 1; Lap[j, j] += 1
Lp = np.linalg.pinv(Lap)
K0 = 2*np.pi - 5*np.arccos(1/3)
shells = sorted(set(np.round(np.degrees(d0), 6)))
def solve(s):
    rhs = -(s - s.mean()); u = Lp @ rhs; u = u - u[far]
    return u, float(np.max(np.abs(Lap @ u - rhs)))
u_zero, _ = solve(np.zeros(NV))
rows = []
for kappa in (1e-4, 3e-4):
    A = 1j*np.sqrt(0.5+s0)*prof; B = 1j*np.sqrt(0.5-s0)*prof
    rho = A**2 + B**2; w = 1 + kappa*(np.abs(rho)-np.abs(rho).mean())/np.abs(rho).mean()
    delta, K = curvature(w); s = (K - K0)/K0
    u, res = solve(s); W = 1 - u + 0.5*u**2
    Z0 = np.sum(A**2+B**2); H0 = np.sum(np.abs(A)**2+np.abs(B)**2); zmax = hmax = 0.0
    for k in range(STEPS):
        c_, s_ = np.cos(W*TH), np.sin(W*TH); A, B = c_*A - s_*B, s_*A + c_*B
        zmax = max(zmax, abs(np.sum(A**2+B**2)-Z0)/abs(Z0)); hmax = max(hmax, abs(np.sum(np.abs(A)**2+np.abs(B)**2)-H0)/H0)
    prof_u = []; prof_W = []
    for sh in shells:
        m = np.isclose(np.degrees(d0), sh); prof_u.append(float(u[m].mean())); prof_W.append(float(W[m].mean()))
    dd = np.radians(np.array(shells)); g = np.where(dd > 1e-9, (np.pi-dd)/np.tan(np.where(dd > 1e-9, dd, 1)), np.nan); g = g - g[-1] if np.isfinite(g[-1]) else g
    mask = dd > 1e-9; uu = np.array(prof_u)
    amp = float(np.sum(uu[mask]*g[mask])/np.sum(g[mask]**2)); resid = float(np.max(np.abs(uu[mask]-amp*g[mask]))/max(np.max(np.abs(uu[mask])), 1e-300))
    rows.append(dict(kappa=kappa, source_sum=float(s.sum()), source_abs_sum=float(np.abs(s).sum()), poisson_residual=res,
                     shells_deg=[float(x) for x in shells], u_by_shell=prof_u, W_by_shell=prof_W,
                     clock_center_over_far=float(W[0]/W[far]), green_fit_amp=amp, green_fit_rel_resid_excl_center=resid,
                     sumz2_max_rel=float(zmax), norm_max_rel=float(hmax)))
judge = dict(zero_source_gives_zero_potential=float(np.max(np.abs(u_zero))) <= 1e-12,
             poisson_residual_small=all(r['poisson_residual'] <= 1e-10 for r in rows),
             sumz2_norm_preserved=all(r['sumz2_max_rel'] <= 1e-12 and r['norm_max_rel'] <= 1e-12 for r in rows))
json.dump(dict(rows=rows, judge=judge), open('stage5b_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for r in rows:
    print('κ', r['kappa'], '| source sum', '%.3e' % r['source_sum'], '(abs sum %.3e)' % r['source_abs_sum'], '| clock center/far %.6f' % r['clock_center_over_far'], '| Green fit resid %.3f' % r['green_fit_rel_resid_excl_center'])
    for sh, uu, ww in zip(r['shells_deg'], r['u_by_shell'], r['W_by_shell']): print('   d=%6.1f°  u=%+.6f  W=%.6f' % (sh, uu, ww))
