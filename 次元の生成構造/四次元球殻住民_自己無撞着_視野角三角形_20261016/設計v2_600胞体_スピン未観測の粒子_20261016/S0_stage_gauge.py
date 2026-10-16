"""S0 舞台とゲージ：群・11 次デザイン・P_g = exp(θ X_u)・左右の可換性・2T⊂2I・ファイバー・23/124。"""
import json, itertools
import numpy as np
from scipy.linalg import expm
import m600_model as m

r = {}
# 群として閉じる
mx = 0.0
for a in m.V:
    for b in m.V:
        q = m.qmul(a, b)
        mx = max(mx, np.min(np.linalg.norm(m.V - q, axis=1)))
r['group_closure_max_dist'] = mx
# 11 次デザイン：5 次までの調和関数が頂点上で正規直交（次数をまたいでも直交）
r['harmonics_dims'] = [B.shape[1] for B in m.BASES]
r['harmonics_orthonormal_err'] = float(np.abs(m.Q.T @ m.Q - np.eye(m.Q.shape[1])).max())
# P_g = exp(θ X_u) on H5、全要素について
err = 0.0
for g in m.V:
    c = np.clip(g[0], -1, 1); al = np.arccos(c)
    if np.sin(al) < 1e-9:
        continue
    u = g[1:] / np.sin(al)
    Xu = sum(ui * X for ui, X in zip(u, m.XL))
    err = max(err, np.abs(m.perm_left(g) @ m.PI - expm(al * Xu) @ m.PI).max())
r['perm_equals_exp_generator_err_all_elements'] = float(err)
r['commutator_XL_YR_max'] = float(max(np.abs(a @ b - b @ a).max() for a in m.XL for b in m.YR))
# su(2) 関係（左）：[X1,X2] = c X3 の c
c12 = m.XL[0] @ m.XL[1] - m.XL[1] @ m.XL[0]
coef = np.sum(c12 * m.XL[2]) / np.sum(m.XL[2] * m.XL[2])
r['su2_left_structure_constant'] = float(coef)
r['su2_left_residual'] = float(np.abs(c12 - coef * m.XL[2]).max())
# 2T ⊂ 2I
hur = [np.eye(4)[i] * s for i in range(4) for s in (1, -1)] + [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
r['hurwitz_units_in_600cell'] = int(sum(np.min(np.linalg.norm(m.V - h, axis=1)) < 1e-9 for h in hur))
# クリフォード・ファイバー（位数 10 の元の左の掛け算の軌道）
g10 = [v for v in m.V if abs(v[0] - m.PHI / 2) < 1e-9][0]
P = m.perm_left(g10)
seen = set(); orbits = []
for i in range(m.NV):
    if i in seen:
        continue
    orb = [i]; j = i
    while True:
        j = int(np.argmax(P[:, j]))
        if j == i:
            break
        orb.append(j)
    seen |= set(orb); orbits.append(len(orb))
r['clifford_fibres_orbit_sizes'] = sorted(set(orbits)); r['clifford_fibres_count'] = len(orbits)
# 2I の元の回転の半角と 23π/124
halfangles = sorted(set(np.round(np.degrees(np.arccos(np.clip(m.V[:, 0], -1, 1))), 9)))
r['element_half_angles_deg'] = halfangles
r['theta_23_124_deg'] = 23 * 180 / 124
r['theta_23_124_is_element_angle'] = bool(any(abs(a - 23 * 180 / 124) < 1e-6 for a in halfangles))
from fractions import Fraction
r['element_orders'] = sorted(set(Fraction(float(a) / 360).limit_denominator(100).denominator for a in halfangles))
crit = (r['group_closure_max_dist'] < 1e-12 and r['harmonics_orthonormal_err'] < 1e-12 and
        r['perm_equals_exp_generator_err_all_elements'] < 1e-12 and r['commutator_XL_YR_max'] < 1e-12 and
        r['hurwitz_units_in_600cell'] == 24 and r['clifford_fibres_count'] == 12 and r['clifford_fibres_orbit_sizes'] == [10]
        and not r['theta_23_124_is_element_angle'])
r['PASS'] = bool(crit)
json.dump(r, open('results/S0_stage_gauge.json', 'w'), indent=1, ensure_ascii=False)
for k, v in r.items():
    print(k, v)
