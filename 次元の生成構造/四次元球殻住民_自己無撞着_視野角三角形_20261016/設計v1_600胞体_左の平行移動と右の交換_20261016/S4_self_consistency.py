"""S4 一つの波の自己無撞着：固有モード → μ → Φ → K → 固有モード を収束まで回す。κ を掃引。"""
import json
import numpy as np
import m600_model as m

H = m.H_E
KAPPAS = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0]
MAXIT = 300
r = {'h': H, 'runs': []}


def pr(psi):
    rho = m.density(psi)
    return float(rho.sum() ** 2 / (m.NV * (rho ** 2).sum()))


def select(K, prev):
    w, Vv = np.linalg.eig(K)
    pos = w.imag > 1e-9
    w, Vv = w[pos].imag, Vv[:, pos]
    order = np.argsort(w); w, Vv = w[order], Vv[:, order]
    groups = []; start = 0
    for i in range(1, len(w) + 1):
        if i == len(w) or w[i] - w[i - 1] > 1e-7:
            groups.append((start, i)); start = i
    best = None
    for a, b in groups:
        Qg, _ = np.linalg.qr(Vv[:, a:b])
        proj = Qg @ (Qg.conj().T @ prev)
        nrm = np.linalg.norm(proj)
        if best is None or nrm > best[0]:
            best = (nrm, proj / nrm, w[a:b].mean(), b - a)
    return best


for kappa in KAPPAS:
    psi, _ = m.packet(H, 0, 0)
    pr0 = pr(psi)
    Phi = m.solve_phi(m.mass_density(psi, H), kappa)
    hist = []
    conv = False
    for it in range(MAXIT):
        ov, new, om, deg = select(m.K_of(H, Phi), psi)
        Phi_new = m.solve_phi(m.mass_density(new, H), kappa)
        dP = float(np.abs(Phi_new - Phi).max())
        hist.append(dict(it=it, overlap=float(ov), omega=float(om), degeneracy=int(deg), dPhi=dP, PR=pr(new)))
        psi, Phi = new, Phi_new
        if dP < 1e-8 and it > 0:
            conv = True
            break
    run = dict(kappa=kappa, converged=conv, iterations=len(hist), PR_start=pr0, PR_end=pr(psi),
               omega_end=hist[-1]['omega'], degeneracy_end=hist[-1]['degeneracy'], overlap_first=hist[0]['overlap'],
               Phi_range=[float(Phi.min()), float(Phi.max())],
               closure_end=float(abs((psi ** 2).sum())), self_consistency_residual=hist[-1]['dPhi'])
    r['runs'].append(run)
    print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in run.items()})
r['PASS'] = bool(all(x['converged'] and x['closure_end'] < 1e-12 for x in r['runs']))
print('PASS', r['PASS'])
json.dump(r, open('results/S4_self_consistency.json', 'w'), indent=1, ensure_ascii=False)
