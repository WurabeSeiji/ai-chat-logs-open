"""確認 2：フルヴィッツの格子 D4 の三重対称性（triality）。
D4 の自己同型をすべて数え、D4*/D4 の三つの剰余類（奇数の類 v、半整数で負号が偶数個の類 s、奇数個の類 c）
への作用を調べる。位数 3 で v→s→c と巡回させる元があるか、それが E8（同じ剰余類どうしの組）を保つか、
フルヴィッツの単位 24 個が v・s・c の三つの 8 個組に分かれるかを確かめる。"""
import json, itertools
import numpy as np

r = {}
roots = [np.array(v, float) for v in itertools.product((-1, 0, 1), repeat=4) if sum(abs(x) for x in v) == 2]
roots = np.array(roots); r['D4_roots'] = len(roots)
simple = np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 1, 1]], float)
G = simple @ simple.T
Sinv = np.linalg.inv(simple)
rootset = {tuple(np.round(x, 9)) for x in roots}


def coset(x):
    """D4*/D4 の類：0（D4）、v（整数・和が奇数）、s/c（半整数、負号の数の偶奇）"""
    x = np.asarray(x, float)
    if np.allclose(x, np.round(x)):
        return '0' if int(round(x.sum())) % 2 == 0 else 'v'
    y = x - 0.5
    if np.allclose(y, np.round(y)):
        # x = (1/2,1/2,1/2,1/2) + 整数ベクトル w。w の和の偶奇で s/c を分ける（D4 を法として）
        w = np.round(y)
        return 's' if int(w.sum()) % 2 == 0 else 'c'
    return '?'


autos = []
for imgs in itertools.product(range(len(roots)), repeat=4):
    R = roots[list(imgs)]
    if not np.allclose(R @ R.T, G):
        continue
    Mt = (Sinv @ R).T          # M simple_i = R_i  →  M = R^T S^{-T}
    M = R.T @ np.linalg.inv(simple.T)
    if not np.allclose(M.T @ M, np.eye(4)):
        continue
    if all(tuple(np.round(M @ x, 9)) in rootset for x in roots):
        autos.append(M)
r['automorphisms_of_D4'] = len(autos)
reps = {'v': np.array([1, 0, 0, 0.]), 's': np.array([.5, .5, .5, .5]), 'c': np.array([-.5, .5, .5, .5])}
perm_counts = {}
triality = None
for M in autos:
    p = tuple(coset(M @ reps[k]) for k in ('v', 's', 'c'))
    perm_counts[''.join(p)] = perm_counts.get(''.join(p), 0) + 1
    if p == ('s', 'c', 'v') and triality is None:
        if np.allclose(np.linalg.matrix_power(M, 3), np.eye(4)):
            triality = M
r['action_on_cosets_(images_of_v,s,c)'] = perm_counts
r['triality_found_order3_cycle_v_s_c'] = triality is not None
if triality is not None:
    r['triality_matrix'] = np.round(triality, 6).tolist()
    # 剰余類の最小ベクトル（各 8 個）を巡回させるか
    units_int = [np.eye(4)[i] * s for i in range(4) for s in (1, -1)]
    units_half = [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
    octets = {'v': units_int, 's': [u for u in units_half if coset(u) == 's'], 'c': [u for u in units_half if coset(u) == 'c']}
    r['hurwitz_units_split_v_s_c'] = {k: len(v) for k, v in octets.items()}
    def octet_image(src):
        imgs = [triality @ u for u in octets[src]]
        cls = {coset(x) for x in imgs}
        norms = {round(float(np.linalg.norm(x)), 9) for x in imgs}
        return dict(classes=sorted(cls), norms=sorted(norms))
    r['octet_images'] = {k: octet_image(k) for k in ('v', 's', 'c')}
    # Hurwitz の 24 単位の集合が保たれるか
    U24 = {tuple(np.round(u, 9)) for u in units_int + units_half}
    r['triality_preserves_24_units'] = all(tuple(np.round(triality @ u, 9)) in U24 for u in U24)
    # E8：(x, y) で x, y ∈ D4*、同じ類。ノルム² が 2 の根 240 個を作り、(T, T) で保たれるか
    cand = [np.array(v) / 2 for v in itertools.product(range(-2, 3), repeat=4)]
    cand = [x for x in cand if coset(x) != '?']
    E8 = []
    for x in cand:
        nx = x @ x
        if nx > 2: continue
        for y in cand:
            if abs(nx + y @ y - 2) < 1e-9 and coset(x) == coset(y):
                E8.append(np.concatenate([x, y]))
    E8set = {tuple(np.round(e, 9)) for e in E8}
    r['E8_roots'] = len(E8)
    TT = np.kron(np.eye(2), triality)
    r['triality_diag_preserves_E8'] = all(tuple(np.round(TT @ e, 9)) in E8set for e in E8)
    # 内積の行列が E8 の根系か（内積の値）
    Ea = np.array(E8); ip = set(np.round((Ea @ Ea.T).ravel(), 6))
    r['E8_inner_products'] = sorted(float(x) for x in ip)
    # 粒子の球の中心 (1/2,1/2,1/2,1/2) の類と、その像
    cpt = reps['s']
    r['ball_centre_class'] = coset(cpt)
    r['ball_centre_images_under_T,T2'] = [coset(triality @ cpt), coset(triality @ triality @ cpt)]
r['PASS'] = bool(r['automorphisms_of_D4'] == 1152 and r['triality_found_order3_cycle_v_s_c']
                 and r.get('triality_preserves_24_units') and r.get('E8_roots') == 240 and r.get('triality_diag_preserves_E8'))
json.dump(r, open('results/check2_triality_d4.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps(r, indent=1, ensure_ascii=False))
