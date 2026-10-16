"""S3（v2：止まっている粒子はスピン未観測の混合・等方の絞り）背景の伸びと三角形の曲率：静止した重い波束を源に Φ、四面体の辺の欠損角、時計の比。"""
import json, itertools
import numpy as np
import m600_model as m

import sys
LW = float(sys.argv[1]) if len(sys.argv) > 1 else m.LAMBDA_WIDTH
KAPPA = 1.0
H_P = 2 * m.H_E
r = {'kappa': KAPPA, 'h_source': H_P, 'lambda_width': LW}
src = m.rest_mixed(H_P, 0, LW)
mu = m.mass_density_mixed(src, H_P)
Phi = m.solve_phi(mu, KAPPA)
N, s = m.stretch(Phi)
r['Phi_min'] = float(Phi.min()); r['Phi_max'] = float(Phi.max()); r['Phi_at_source'] = float(Phi[0])
# 連続の曲率 R s^3（左不変の正規直交枠 Y_i、ラプラシアンは H5 上）
lns = np.log(s)
lap_lns = m.LAP @ (m.PI @ lns)
grad2 = sum((Y @ (m.PI @ lns)) ** 2 for Y in m.YR)
R = s ** -2 * (6 - 4 * lap_lns - 2 * grad2)
C = R * s ** 3 - 6
# 四面体（4 頂点がすべて隣接）と辺
adj = np.abs(m.GRAM - m.PHI / 2) < 1e-9
edges = [(i, j) for i in range(m.NV) for j in range(i + 1, m.NV) if adj[i, j]]
eidx = {e: k for k, e in enumerate(edges)}
tets = []
nb = [set(np.nonzero(adj[i])[0]) for i in range(m.NV)]
for i in range(m.NV):
    for j in nb[i]:
        if j <= i: continue
        for k in nb[i] & nb[j]:
            if k <= j: continue
            for l in nb[i] & nb[j] & nb[k]:
                if l <= k: continue
                tets.append((i, j, k, l))
r['n_edges'] = len(edges); r['n_tets'] = len(tets)


def dihedrals(L):
    """辺の長さ L[(a,b)] の平坦な四面体の 6 本の辺の二面角。頂点 0..3。"""
    D = np.zeros((4, 4))
    for (a, b), v in L.items():
        D[a, b] = D[b, a] = v ** 2
    G = np.array([[(D[0, a] + D[0, b] - D[a, b]) / 2 for b in (1, 2, 3)] for a in (1, 2, 3)])
    w, U = np.linalg.eigh(G)
    X = np.vstack([np.zeros(3), (U * np.sqrt(np.clip(w, 0, None))).T @ np.eye(3) if False else U @ np.diag(np.sqrt(np.clip(w, 0, None)))])
    P = np.vstack([np.zeros(3), X[1:]])
    out = {}
    for a, b in itertools.combinations(range(4), 2):
        c, d = [x for x in range(4) if x not in (a, b)]
        e = P[b] - P[a]; e /= np.linalg.norm(e)
        u = P[c] - P[a]; u -= (u @ e) * e
        v = P[d] - P[a]; v -= (v @ e) * e
        out[(a, b)] = np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))
    return out


def deficits(svec):
    ell = {e: m.STEP * (svec[e[0]] + svec[e[1]]) / 2 for e in edges}
    tot = np.zeros(len(edges))
    for t in tets:
        L = {(a, b): ell[tuple(sorted((t[a], t[b])))] for a, b in itertools.combinations(range(4), 2)}
        dh = dihedrals(L)
        for (a, b), ang in dh.items():
            tot[eidx[tuple(sorted((t[a], t[b])))]] += ang
    dl = 2 * np.pi - tot
    return dl, np.array([ell[e] for e in edges])


d0, l0 = deficits(np.ones(m.NV))
d1, l1 = deficits(s)
r['flat_deficit_deg'] = float(np.degrees(d0.mean())); r['flat_deficit_spread'] = float(np.ptp(d0))
dDL = d1 * l1 - d0 * l0
Ce = np.array([(C[i] + C[j]) / 2 for i, j in edges])
mut = m.PI @ mu - (m.PI @ mu).mean()
mue = np.array([(mut[i] + mut[j]) / 2 for i, j in edges])
r['corr_regge_vs_continuum_RS3'] = float(np.corrcoef(dDL, Ce)[0, 1])
r['corr_regge_vs_mu_minus_mean_reference'] = float(np.corrcoef(dDL, mue)[0, 1])
r['sum_dDL'] = float(dDL.sum()); r['sum_abs_dDL'] = float(np.abs(dDL).sum())
r['sum_ratio'] = float(abs(dDL.sum()) / np.abs(dDL).sum())
# 時計の比：源からの距離ごとに、静止した電子役の波束を置く
ang = np.degrees(np.arccos(np.clip(m.V @ m.V[0], -1, 1)))
shells = sorted(set(np.round(ang, 6)))
K0 = m.K_of(m.H_E); K1 = m.K_of(m.H_E, Phi)
clock = []
for a in shells:
    i = int(np.nonzero(np.abs(ang - a) < 1e-6)[0][0])
    pr = m.rest_mixed(m.H_E, i, LW)
    w0 = sum(np.real(c.conj() @ (-1j * K0) @ c) for c in pr); w1 = sum(np.real(c.conj() @ (-1j * K1) @ c) for c in pr)
    clock.append(dict(distance_deg=float(a), clock_ratio=float(w1 / w0), Phi_at=float(Phi[i])))
r['clock'] = clock
r['Phi_max_spread_within_shell'] = float(max(np.ptp(Phi[np.abs(ang - a) < 1e-6]) for a in shells))
cr = [c['clock_ratio'] for c in clock]
r['clock_monotone_increasing_with_distance'] = bool(np.all(np.diff(cr) > 0))
r['PASS'] = bool(r['corr_regge_vs_continuum_RS3'] > 0.9 and r['sum_ratio'] <= 0.05 and r['clock_monotone_increasing_with_distance'])
json.dump(r, open('results/S3_stretch_triangles%s.json' % ('' if LW == m.LAMBDA_WIDTH else '_lw%g' % LW), 'w'), indent=1, ensure_ascii=False)
for k, v in r.items():
    if k != 'clock': print(k, v)
for c in clock: print(c)
