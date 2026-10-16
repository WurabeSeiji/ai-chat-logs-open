# -*- coding: utf-8 -*-
"""段 2：倍音の表現と並進の資格審査（設計書 v2 §5、検査は事前固定、すべて実装の正しさ）
N = 18（総 2,470 成分）、乱数種 20261016 の複素係数。"""
import json, numpy as np
from s3_harmonics import *
rng = np.random.default_rng(20261016); N = 18
C0 = [rng.normal(size=(n+1, n+1)) + 1j*rng.normal(size=(n+1, n+1)) for n in range(N+1)]
diff = lambda A, B: max(np.max(np.abs(a-b)) for a, b in zip(A, B))
res = {}
# 1 並進の合成
for u in 'ijk':
    A = translate_left(translate_left(C0, qexp(0.37, u)), qexp(1.11, u)); Bc = translate_left(C0, qexp(1.48, u))
    res[f'1_compose_{u}'] = float(diff(A, Bc))
# 2 奇偶：角 π（四元数 −1）で奇数次が −1 倍、偶数次は不変。角 2π で全次数不変
Cpi = translate_left(C0, qexp(np.pi, 'j'))
res['2_parity_pi'] = float(max(np.max(np.abs(Cpi[n] - (-1)**n*C0[n])) for n in range(N+1)))
res['2_full_turn_2pi'] = float(diff(translate_left(C0, qexp(2*np.pi, 'k')), C0))
# 3 左右の可換（並進・ブーストと交換の回転）
th = 2*np.pi*23/124
for name, a in (('translation', qexp(0.7, 'i')), ('boost', qexp(1j*0.8, 'i'))):
    for ub in 'ijk':
        b = qexp(th, ub)
        res[f'3_commute_{name}_right{ub}'] = float(diff(rotate_right(translate_left(C0, a), b), translate_left(rotate_right(C0, b), a)))
# 4 保存：並進・回転で B とノルムが保存、ブーストで B が保存しノルムは変わる
B0, H0 = true_square(C0), norm_h(C0)
Ct = translate_left(C0, qexp(0.9, 'k')); Cr = rotate_right(C0, qexp(th, 'j')); Cb = translate_left(C0, qexp(1j*0.8, 'j'))
res['4_B_translation'] = float(abs(true_square(Ct)-B0)/abs(B0)); res['4_H_translation'] = float(abs(norm_h(Ct)-H0)/H0)
res['4_B_rotation'] = float(abs(true_square(Cr)-B0)/abs(B0)); res['4_H_rotation'] = float(abs(norm_h(Cr)-H0)/H0)
res['4_B_boost'] = float(abs(true_square(Cb)-B0)/abs(B0)); res['4_H_boost_relative_change'] = float(abs(norm_h(Cb)-H0)/H0)
# 4b（記録）：B が S³ 上の f² の平均に等しいか（モンテカルロ 20,000 点、N=4 の部分で）
Cs = C0[:5]; qs = rng.normal(size=(20000, 4)); qs /= np.linalg.norm(qs, axis=1)[:, None]
mc = np.mean([field(Cs, q)**2 for q in qs]); res['4b_MC_mean_f2_over_B'] = [float((mc/true_square(Cs)).real), float((mc/true_square(Cs)).imag)]
# 5 標本点：600 胞体の頂点 g での並進と、頂点の値の並べ替えが一致するか（N=6 の部分で）
P = icosians(); Cs = C0[:7]; g = P[37]
Cg = translate_left(Cs, g)
vals = np.array([field(Cs, p) for p in P]); idx = {tuple(np.round(p, 6)): i for i, p in enumerate(P)}
gi = qinv(g).real; err = 0.0
for p in P[:40]:
    src = tuple(np.round(qmul(gi, p).real, 6)); err = max(err, abs(field(Cg, p) - vals[idx[src]]))
res['5_icosian_permutation'] = float(err); res['5_n_vertices'] = int(len(P))
# 6 有限位数：交換の回転 θ=2π·23/124 を 124 回で恒等、62 回で四元数 −1（対蹠点）
b = qexp(th, 'i'); C = C0
for k in range(124):
    C = rotate_right(C, b)
    if k == 61: C62 = C
res['6_right_124_identity'] = float(diff(C, C0))
res['6_right_62_parity'] = float(max(np.max(np.abs(C62[n] - (-1)**n*C0[n])) for n in range(N+1)))
TH = 1e-12
judge = {k: (v <= TH) for k, v in res.items() if isinstance(v, float) and not k.startswith('4_H_boost') and not k.startswith('4b')}
judge['4_H_boost_changes'] = res['4_H_boost_relative_change'] > 1e-6
out = dict(N=N, threshold=TH, results=res, judge=judge, all_pass=all(judge.values()))
json.dump(out, open('stage2_results.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
