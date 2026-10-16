# -*- coding: utf-8 -*-
"""段 5c（README の事前固定のとおり）"""
import json, numpy as np
exec(open('stage5a_geometry_copy.py').read().split("TH = 23*np.pi/124")[0])     # P, edges, tets, d_geo, curvature()
from s3_harmonics import qmul
TH = 23*np.pi/124; s0 = 0.01; sig = 0.6; phi_ = (1+5**0.5)/2
d0 = d_geo[0]; far = int(np.argmax(d0)); K0 = 2*np.pi - 5*np.arccos(1/3)
Lap = np.zeros((NV, NV))
for i, j in edges: Lap[i, j] -= 1; Lap[j, i] -= 1; Lap[i, i] += 1; Lap[j, j] += 1
Lp = np.linalg.pinv(Lap)
g = P[[k for k in range(NV) if np.isclose(P[k][0], phi_/2)][0]]                  # 隣の頂点の群の元（実部 φ/2）
key = {tuple(np.round(p, 6)): i for i, p in enumerate(P)}
Pg = np.zeros((NV, NV))
for i, p in enumerate(P): Pg[key[tuple(np.round(qmul(p, g), 6))], i] = 1.0       # 右から g を掛ける置換
D = (Pg - Pg.T)/2
def build_K(W):
    n = NV; K = np.zeros((2*n, 2*n)); Wd = np.diag(W); sym = 0.5*(Wd @ D + D @ Wd)
    K[:n, :n] = TH*sym; K[n:, n:] = TH*sym                                      # 右の並進（A、B それぞれ）
    for p in range(n): K[p, n+p] = -TH*W[p]; K[n+p, p] = TH*W[p]               # 各点の交換（時計）
    return K
def geometry(v, kappa):
    A, B = v[:NV], v[NV:]; rho = A**2 + B**2
    w = 1 + kappa*(np.abs(rho)-np.abs(rho).mean())/np.abs(rho).mean()
    delta, Kc = curvature(w); s = (Kc - K0)/K0
    u = Lp @ (s - s.mean()); u = u - u[far]
    return s, u, 1 - u + 0.5*u**2, rho
rows = []
for kappa in (0.0, 1e-4, 3e-4):
    prof = np.exp(-d0**2/sig**2)
    v = np.concatenate([1j*np.sqrt(0.5+s0)*prof, 1j*np.sqrt(0.5-s0)*prof]); v = v/np.linalg.norm(v)
    hist = []; conv = None
    for it in range(1, 201):
        s, u, W, rho = geometry(v, kappa); K = build_K(W)
        mu, V = np.linalg.eig(K)
        ov = np.abs(V.conj().T @ v); k = int(np.argmax(ov)); vn = V[:, k]/np.linalg.norm(V[:, k])
        vn = vn*np.exp(-1j*np.angle(np.vdot(v, vn)))                                 # 位相をそろえる
        change = 1 - abs(np.vdot(v, vn)); hist.append(float(change)); v = vn
        if change < 1e-10: conv = it; break
    s, u, W, rho = geometry(v, kappa); K = build_K(W)
    shells = sorted(set(np.round(np.degrees(d0), 6)))
    rows.append(dict(kappa=kappa, converged_at=conv, last_change=hist[-1], change_history_head=hist[:6],
                     K_antisym=float(np.max(np.abs(K + K.T))), eig=[float(mu[k].real), float(mu[k].imag)],
                     sumz2_over_norm=float(abs(v @ v)/np.vdot(v, v).real),
                     center_rho_fraction=float(np.abs(rho[0])/np.abs(rho).sum()), rho_abs_max_vertex=int(np.argmax(np.abs(rho))),
                     source_sum=float(s.sum()), source_abs_sum=float(np.abs(s).sum()),
                     u_by_shell=[float(u[np.isclose(np.degrees(d0), sh)].mean()) for sh in shells],
                     W_by_shell=[float(W[np.isclose(np.degrees(d0), sh)].mean()) for sh in shells], shells_deg=[float(x) for x in shells]))
judge = dict(K_real_antisymmetric=all(r['K_antisym'] == 0.0 for r in rows), eigenmode_closure=all(r['sumz2_over_norm'] <= 1e-12 for r in rows),
             kappa0_W_is_one=all(abs(x-1) < 1e-12 for x in rows[0]['W_by_shell']))
json.dump(dict(rows=rows, judge=judge, g_used=[float(x) for x in g]), open('stage5c_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for r in rows:
    print('κ', r['kappa'], '| converged at', r['converged_at'], '| last change %.1e' % r['last_change'], '| eig', np.round(r['eig'], 6), '| |Σz²|/norm %.1e' % r['sumz2_over_norm'])
    print('   center |ρ| fraction %.4f  max at vertex %d | source sum %.3e (abs %.3e)' % (r['center_rho_fraction'], r['rho_abs_max_vertex'], r['source_sum'], r['source_abs_sum']))
    print('   W by shell', np.round(r['W_by_shell'], 6))
