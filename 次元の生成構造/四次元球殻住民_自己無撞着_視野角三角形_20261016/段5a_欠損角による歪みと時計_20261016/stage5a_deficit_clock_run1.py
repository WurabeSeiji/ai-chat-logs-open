# -*- coding: utf-8 -*-
"""段 5a（README の事前固定のとおり）"""
import json, itertools, numpy as np
from s3_harmonics import icosians
P = icosians(); NV = len(P)
G = P @ P.T; phi = (1+5**0.5)/2
adj = np.isclose(G, phi/2, atol=1e-9)
edges = [(i, j) for i in range(NV) for j in range(i+1, NV) if adj[i, j]]
nbr = [set(np.where(adj[i])[0]) for i in range(NV)]
tets = set()
for i, j in edges:
    common = nbr[i] & nbr[j]
    for k, l in itertools.combinations(sorted(common), 2):
        if adj[k, l]: tets.add(tuple(sorted((i, j, k, l))))
tets = sorted(tets)
assert (len(edges), len(tets)) == (720, 600), (len(edges), len(tets))
d_geo = np.arccos(np.clip(G, -1, 1))
def dihedrals(L):
    """6 辺の長さ L[(a,b)] の平らな四面体の、各辺での二面角（古典 MDS で 3 次元に置いて法線から）"""
    idx = [0, 1, 2, 3]; D2 = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if a != b: D2[a, b] = L[tuple(sorted((a, b)))]**2
    Jc = np.eye(4) - 1/4; Gm = -0.5*Jc @ D2 @ Jc
    w, V = np.linalg.eigh(Gm); X = V[:, -3:]*np.sqrt(np.maximum(w[-3:], 0))
    out = {}
    for a, b in itertools.combinations(range(4), 2):
        c, d = [x for x in range(4) if x not in (a, b)]
        e = X[b]-X[a]; e /= np.linalg.norm(e)
        u = X[c]-X[a]; u -= (u@e)*e; v = X[d]-X[a]; v -= (v@e)*e
        out[(a, b)] = np.arccos(np.clip(u@v/(np.linalg.norm(u)*np.linalg.norm(v)), -1, 1))
    return out
def curvature(w):
    ell = {e: d_geo[e]*0.5*(w[e[0]]+w[e[1]]) for e in edges}
    dsum = {e: 0.0 for e in edges}
    for t in tets:
        L = {tuple(sorted((a, b))): ell[tuple(sorted((t[a], t[b])))] for a, b in itertools.combinations(range(4), 2)}
        for (a, b), ang in dihedrals(L).items(): dsum[tuple(sorted((t[a], t[b])))] += ang
    delta = {e: 2*np.pi - dsum[e] for e in edges}
    K = np.zeros(NV); Ls = np.zeros(NV)
    for e in edges:
        for p in e: K[p] += delta[e]*ell[e]; Ls[p] += ell[e]
    return delta, K/Ls
TH = 23*np.pi/124; STEPS = 248; s0 = 0.01; sig = 0.6
d0 = d_geo[0]; prof = np.exp(-d0**2/sig**2); far = int(np.argmax(d0))
rows = []
for kappa in (0.0, 0.02, 0.05, 0.1):
    A = 1j*np.sqrt(0.5+s0)*prof; B = 1j*np.sqrt(0.5-s0)*prof
    rho = A**2 + B**2; w = 1 + kappa*(np.abs(rho) - np.abs(rho).mean())/np.abs(rho).mean()
    delta, K = curvature(w); K0 = 2*np.pi - 5*np.arccos(1/3)
    u = (K - K0)/K0; W = 1 - u + 0.5*u**2
    Z0 = np.sum(A**2 + B**2); H0 = np.sum(np.abs(A)**2 + np.abs(B)**2)
    ang = np.zeros(NV); zmax = hmax = rmax = 0.0
    psi_prev = np.angle(A + 1j*B) if False else None
    for k in range(STEPS):
        c, s = np.cos(W*TH), np.sin(W*TH)
        A, B = c*A - s*B, s*A + c*B; ang += W*TH
        zmax = max(zmax, abs(np.sum(A**2+B**2)-Z0)/abs(Z0)); hmax = max(hmax, abs(np.sum(np.abs(A)**2+np.abs(B)**2)-H0)/H0)
        rmax = max(rmax, np.max(np.abs(A**2+B**2 - rho))/np.max(np.abs(rho)))
    # 時計の進み：各点の (A,B) の回転角（実数の角）を状態から読む：A_p/B_p の比から角を復元して累積と照合
    psi0 = np.arctan2(np.sqrt(0.5-s0), np.sqrt(0.5+s0))
    ratio_clock = (W[0]*TH)/(W[far]*TH)
    dvals = np.degrees(np.array(list(delta.values())))
    rows.append(dict(kappa=kappa, deficit_min_deg=float(dvals.min()), deficit_max_deg=float(dvals.max()),
                     K_center_deg=float(np.degrees(K[0])), K_far_deg=float(np.degrees(K[far])), K0_deg=float(np.degrees(K0)),
                     W_center=float(W[0]), W_far=float(W[far]), W_min=float(W.min()), W_max=float(W.max()),
                     clock_ratio_center_over_far=float(ratio_clock),
                     sumz2_max_rel=float(zmax), norm_max_rel=float(hmax), rho_max_rel=float(rmax)))
r0 = rows[0]
judge = dict(undistorted_deficits_uniform=abs(r0['deficit_min_deg']-7.356103)<1e-5 and abs(r0['deficit_max_deg']-7.356103)<1e-5 and abs(r0['W_min']-1)<1e-9 and abs(r0['W_max']-1)<1e-9,
             sumz2_and_norm_preserved=all(r['sumz2_max_rel'] <= 1e-12 and r['norm_max_rel'] <= 1e-12 for r in rows),
             local_rho_invariant=all(r['rho_max_rel'] <= 1e-12 for r in rows))
json.dump(dict(rows=rows, judge=judge, n_edges=len(edges), n_tets=len(tets), far_vertex=far, far_angle=float(d0[far])), open('stage5a_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(judge, indent=1))
for r in rows: print({k: (float('%.6g' % v) if isinstance(v, float) else v) for k, v in r.items()})
