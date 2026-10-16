"""S4（v2）一つの波の自己無撞着：最も低い正の振動数の準位をクラマース対ごと占め（閉殻）、
μ → Φ → K → 準位 を収束まで回す。出発点の Φ は止まっている粒子（スピン未観測の混合）から作る。κ を掃引。"""
import json
import numpy as np
import m600_model as m

H = m.H_E
KAPPAS = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0]
MAXIT = 300
r = {'h': H, 'runs': []}


def pr(rho):
    return float(rho.sum() ** 2 / (m.NV * (rho ** 2).sum()))


def lowest_shell(K):
    w, Vv = np.linalg.eig(K)
    pos = w.imag > 1e-9
    w, Vv = w[pos].imag, Vv[:, pos]
    o = np.argsort(w); w, Vv = w[o], Vv[:, o]
    b = 1
    while b < len(w) and w[b] - w[0] < 1e-7:
        b += 1
    Qg, _ = np.linalg.qr(Vv[:, :b])
    return [Qg[:, j] / np.sqrt(b) for j in range(b)], float(w[:b].mean()), b


for kappa in KAPPAS:
    start = m.rest_mixed(H, 0)
    Phi = m.solve_phi(m.mass_density_mixed(start, H), kappa)
    hist = []; conv = False
    for it in range(MAXIT):
        comps, om, deg = lowest_shell(m.K_of(H, Phi))
        Phi_new = m.solve_phi(m.mass_density_mixed(comps, H), kappa)
        dP = float(np.abs(Phi_new - Phi).max())
        hist.append(dict(it=it, omega=om, degeneracy=deg, dPhi=dP, PR=pr(m.density_mixed(comps))))
        Phi = Phi_new
        if dP < 1e-8:
            conv = True
            break
    run = dict(kappa=kappa, converged=conv, iterations=len(hist), PR_start=pr(m.density_mixed(start)),
               PR_end=hist[-1]['PR'], omega_end=om, bound=bool(om < H - 1e-9), degeneracy_end=deg,
               Phi_range=[float(Phi.min()), float(Phi.max())],
               closure_end=float(max(abs((c ** 2).sum()) for c in comps)), self_consistency_residual=hist[-1]['dPhi'])
    r['runs'].append(run)
    print({k: (float('%.6g' % v) if isinstance(v, float) else v) for k, v in run.items()})
r['PASS'] = bool(all(x['converged'] and x['closure_end'] < 1e-12 for x in r['runs']))
print('PASS', r['PASS'])
json.dump(r, open('results/S4_self_consistency.json', 'w'), indent=1, ensure_ascii=False)
