#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""10,000歩終状態の最終検査（読み出しのみ）: 凍結性・パリティ・因数分解・無名性。
各 N の step 10000 について:
 (1) 剛体残差 r = min_φ‖Z[s+N]−e^{iφ}Z[s]‖/‖Z‖ と固定点残差 ‖Z[s+N]−Z[s]‖/‖Z‖
     （s=10000−N。ガラスか結晶かによらず「凍結」しているかの判定）
 (2) パリティ純度（6方向格子割り当て、E/ρ 信頼度つき）
 (3) Z[ω] 格子証明書つき因数分解（E<ρ/2 なら厳密判定、不成立なら 'no-cert'）
 (4) 無名性 1−S（振幅エントロピー）
出力: final_states_10000_table_v1.csv、check_final_states_10000_v1.json"""
import csv
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))

def parity_and_cert(z):
    s = z.astype(np.complex128) ** 2
    rho = float(np.median(np.abs(s)))
    u = np.angle(s) / (math.pi / 3.0)
    off = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    k = np.round(u - off).astype(int) % 6
    ideal = rho * np.exp(1j * (math.pi / 3.0) * (k + off))
    e_rel = float(np.sum(np.abs(s - ideal))) / rho
    a = [int(np.sum(k == j)) for j in range(6)]
    n_even = a[0] + a[2] + a[4]; n_odd = a[1] + a[3] + a[5]
    purity = abs(n_even - n_odd) / len(z)
    cert = bool(e_rel < 0.5)
    verdict = None
    if cert:
        Mn = len(z)
        pair = any(a[j] >= 1 and a[j + 3] >= 1 for j in range(3)) and Mn > 2
        tri = any(all(a[(j + 2 * i) % 6] >= 1 for i in range(3)) for j in (0, 1)) and Mn > 3
        verdict = 'prime' if not (pair or tri) else 'composite'
    return purity, e_rel, cert, verdict, a

rows = []
detail = {}
for N in range(3, 41):
    d = np.load(os.path.join(BASE, 'results', f'hm_N{N}_den_{N}_states_10000.npz'))
    Z = np.asarray(d['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    s = 10000 - N
    ip = np.vdot(Z[s], Z[s + N]); ph = float(np.angle(ip))
    nrm = float(np.linalg.norm(Z[s]))
    r_rigid = float(np.linalg.norm(Z[s + N] - np.exp(1j * ph) * Z[s])) / nrm
    r_fix = float(np.linalg.norm(Z[s + N] - Z[s])) / nrm
    purity, e_rel, cert, verdict, a = parity_and_cert(Z[10000])
    P = np.abs(Z[10000]) ** 2; P = P / P.sum()
    with np.errstate(divide='ignore', invalid='ignore'):
        S = float(-(P * np.where(P > 0, np.log(P), 0.0)).sum() / math.log(M))
    rows.append([N, M, M % 3, f'{r_rigid:.2e}', f'{r_fix:.2e}',
                 round(purity, 4), round(e_rel, 2), cert, verdict or '-', f'{1 - S:.2e}'])
    detail[N] = {'r_rigid': r_rigid, 'r_fixedpoint': r_fix, 'parity_purity': purity,
                 'E_over_rho': e_rel, 'certified': cert, 'closure_verdict': verdict,
                 'occ6': a, 'one_minus_S': 1 - S}
    print(f"N={N}: 剛体残差={r_rigid:.2e} 固定点残差={r_fix:.2e} 純度={purity:.3f} "
          f"E/ρ={e_rel:.2f} cert={cert} verdict={verdict} 1−S={1-S:.2e}", flush=True)

HEADER = ['N', 'M', 'M_mod_3', 'r_rigid', 'r_fixedpoint', 'parity_purity',
          'E_over_rho', 'certified', 'closure_verdict', 'one_minus_S']
with open(os.path.join(BASE, 'final_states_10000_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(HEADER); w.writerows(rows)
with open(os.path.join(BASE, 'check_final_states_10000_v1.json'), 'w') as f:
    json.dump(detail, f, indent=2, default=str)
print('ALL DONE')
