"""S5 二体（陽子役＋電子役）：同じファイバー上にゲージ 3 目盛（108°）離して置き、共有の Φ の中で発展させる。"""
import json
import numpy as np
from scipy.linalg import expm
import m600_model as m

H_E = m.H_E; H_P = 2 * m.H_E
DT = 0.05; T = 3.0; T_JUDGE = 1.25
KAPPAS = [0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0]
g3 = m.qmul(m.G10, m.qmul(m.G10, m.G10))     # e^{u·108°}
i_p = m.index_of(g3)
r = {'h_e': H_E, 'h_p': H_P, 'theta_e': m.THETA_E, 'theta_p': H_P * m.STEP, 'R_p=cos^2theta_p': float(np.cos(H_P * m.STEP) ** 2),
     'start_separation_deg': float(np.degrees(np.arccos(np.clip(m.V[0] @ m.V[i_p], -1, 1)))), 'runs': []}


def run(kappa):
    pe, _ = m.packet(H_E, 0, 0)
    pp, _ = m.packet(H_P, 0, i_p)
    s0 = [abs((pe ** 2).sum()), abs((pp ** 2).sum())]
    seps = []; clocks = []
    for st in range(int(round(T / DT)) + 1):
        mu = m.mass_density(pe, H_E) + m.mass_density(pp, H_P)
        Phi = m.solve_phi(mu, kappa)
        a, b = m.position(pe), m.position(pp)
        seps.append(float(np.arccos(np.clip(a @ b, -1, 1))))
        Ke = m.K_of(H_E, Phi); Kp = m.K_of(H_P, Phi)
        clocks.append(float(np.real(pe.conj() @ (-1j * Ke) @ pe)))
        pe = expm(DT * Ke) @ pe; pp = expm(DT * Kp) @ pp
    return dict(kappa=kappa, sep=seps, clock_e=clocks,
                closure_start=s0, closure_end=[float(abs((pe ** 2).sum())), float(abs((pp ** 2).sum()))],
                norm_end=[float(np.linalg.norm(pe)), float(np.linalg.norm(pp))])


for k in KAPPAS:
    x = run(k); r['runs'].append(x)
    print(k, np.round(x['sep'][::5], 4), 'closure', np.round(x['closure_end'], 16), 'norm', np.round(x['norm_end'], 13))
ts = np.arange(len(r['runs'][0]['sep'])) * DT
j = int(round(T_JUDGE / DT))
base = np.array(r['runs'][0]['sep'])
dsep = {x['kappa']: float(x['sep'][j] - base[j]) for x in r['runs']}
r['delta_sep_at_judge'] = {str(k): v for k, v in dsep.items()}
sign_ok = all(dsep[k] < 0 < dsep[-k] for k in (0.5, 1.0, 2.0))
lin = abs(dsep[1.0] / dsep[0.5] / 2 - 1)
sym = max(abs(dsep[k] + dsep[-k]) / abs(dsep[k]) for k in (0.5, 1.0))
r['sign_flip_ok'] = bool(sign_ok); r['linearity_dev_1_vs_0.5'] = float(lin); r['antisymmetry_dev'] = float(sym)
cons = max(max(x['closure_end']) for x in r['runs']); nrm = max(max(abs(np.array(x['norm_end']) - 1)) for x in r['runs'])
r['max_closure_end'] = float(cons); r['max_norm_dev'] = float(nrm)
r['PASS'] = bool(sign_ok and lin <= 0.3 and sym <= 0.3 and cons < 1e-12 and nrm < 1e-12)
for k in ['delta_sep_at_judge', 'sign_flip_ok', 'linearity_dev_1_vs_0.5', 'antisymmetry_dev', 'max_closure_end', 'max_norm_dev', 'PASS']:
    print(k, r[k])
json.dump(r, open('results/S5_two_body.json', 'w'), indent=1, ensure_ascii=False)
