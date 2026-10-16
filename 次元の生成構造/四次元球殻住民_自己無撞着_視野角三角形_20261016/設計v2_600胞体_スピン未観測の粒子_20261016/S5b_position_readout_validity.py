"""S5 補足（自己査読用、合否に使わない）：位置の読み（密度の平均ベクトルを正規化）が有効な時間を、
平均ベクトルの大きさ |Σρq| で確かめる（κ=0、自由な広がり）。"""
import json
import numpy as np
from scipy.linalg import expm
import m600_model as m

DT = 0.05; T = 3.0
g3 = m.qmul(m.G10, m.qmul(m.G10, m.G10)); i_p = m.index_of(g3)
out = {}
for name, h, i0 in [('electron', m.H_E, 0), ('proton', 2 * m.H_E, i_p)]:
    c = m.rest_mixed(h, i0); U = expm(DT * m.K_of(h)); rows = []
    for st in range(int(round(T / DT)) + 1):
        rows.append(float(np.linalg.norm(m.density_mixed(c) @ m.V)))
        c = [U @ x for x in c]
    out[name] = rows
    print(name, 'mean-vector norm at beta=0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0:',
          [round(rows[int(round(b / DT))], 4) for b in (0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)])
json.dump(out, open('results/S5b_position_readout_validity.json', 'w'), indent=1)
