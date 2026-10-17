"""確認 1：本当の二乗 Σz² とノルム Σ|z|² の対称性。
(a) 両方を保つ線形変換は実の直交変換に限られる（数値の例示。証明は README）。
(b) 項の置き場所の入替（座標に i を掛ける）は、符号の異なる群 O(p,q) を同じ O(n,C) の中へ写す。"""
import json
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(20261017)
n = 8
r = {}


def antisym(Z):
    return (Z - Z.T) / 2


# (a) 複素直交（本当の二乗を保つ）・ユニタリー（ノルムを保つ）・実直交（両方）
Kc = antisym(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))   # 複素反対称
Hh = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); Hh = (Hh + Hh.conj().T) / 2
Kr = antisym(rng.normal(size=(n, n)))                                   # 実反対称
T = {'complex_orthogonal': expm(Kc), 'unitary': expm(1j * Hh), 'real_orthogonal': expm(Kr)}
z = rng.normal(size=n) + 1j * rng.normal(size=n)
for name, M in T.items():
    w = M @ z
    r[name] = dict(true_square_change=float(abs((w @ w) - (z @ z))),
                   norm_change=float(abs(np.vdot(w, w).real - np.vdot(z, z).real)),
                   imag_part_max=float(np.abs(M.imag).max()))
# 生成子の交わり：複素反対称 ∩ 反エルミート = 実反対称（次元で確認）
basis = []
for i in range(n):
    for j in range(n):
        for part in (1, 1j):
            E = np.zeros((n, n), complex); E[i, j] = part; basis.append(E)
rows = []
for E in basis:
    rows.append(np.concatenate([(E + E.T).ravel().real, (E + E.T).ravel().imag,
                                (E + E.conj().T).ravel().real, (E + E.conj().T).ravel().imag]))
A = np.array(rows).T
null_dim = int(len(basis) - np.linalg.matrix_rank(A))
r['generator_intersection_dim'] = null_dim
r['real_antisymmetric_dim'] = n * (n - 1) // 2
# (b) 置き場所の入替：D = diag(1,1,1,i)（t に i を掛ける）。O(3,1) の元 Λ について DΛD⁻¹ ∈ O(4,C)
eta = np.diag([1., 1., 1., -1.])
rap = 0.83
boost = np.eye(4); boost[0, 0] = boost[3, 3] = np.cosh(rap); boost[0, 3] = boost[3, 0] = np.sinh(rap)
rot = np.eye(4); c, s = np.cos(0.4), np.sin(0.4); rot[0, 0] = rot[1, 1] = c; rot[0, 1] = -s; rot[1, 0] = s
Lam = rot @ boost
Dm = np.diag([1, 1, 1, 1j])
Mc = Dm @ Lam @ np.linalg.inv(Dm)
r['lorentz_preserves_eta'] = float(np.abs(Lam.T @ eta @ Lam - eta).max())
r['conjugated_is_complex_orthogonal'] = float(np.abs(Mc.T @ Mc - np.eye(4)).max())
r['conjugated_is_not_real'] = float(np.abs(Mc.imag).max())
# 同じ O(4,C) に符号 (4,0)・(3,1)・(2,2) の群がすべて入る
sig_ok = {}
for q in (0, 1, 2):
    et = np.diag([1.] * (4 - q) + [-1.] * q)
    Dq = np.diag([1] * (4 - q) + [1j] * q)
    # O(4-q,q) の元を生成子 X（X^T η + η X = 0）から作る
    X = rng.normal(size=(4, 4)); X = (X - et @ X.T @ et) / 2
    L = expm(0.5 * X)
    Mq = Dq @ L @ np.linalg.inv(Dq)
    sig_ok[f'({4 - q},{q})'] = dict(preserves_eta=float(np.abs(L.T @ et @ L - et).max()),
                                    in_O4C=float(np.abs(Mq.T @ Mq - np.eye(4)).max()))
r['signatures_inside_O4C'] = sig_ok
r['PASS'] = bool(r['real_orthogonal']['true_square_change'] < 1e-12 and r['real_orthogonal']['norm_change'] < 1e-12
                 and r['complex_orthogonal']['true_square_change'] < 1e-12 and r['complex_orthogonal']['norm_change'] > 1e-3
                 and r['unitary']['norm_change'] < 1e-12 and r['unitary']['true_square_change'] > 1e-3
                 and null_dim == r['real_antisymmetric_dim'] and r['conjugated_is_complex_orthogonal'] < 1e-12
                 and all(v['in_O4C'] < 1e-12 for v in sig_ok.values()))
json.dump(r, open('results/check1_true_square_and_norm.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps(r, indent=1, ensure_ascii=False))
