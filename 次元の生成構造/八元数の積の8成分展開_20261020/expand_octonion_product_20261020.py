#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""八元数の積を 8 成分すべて展開して書き出す（2026-10-20、木原の依頼「8 要素を完全に展開して描いて」）。
掛け算の定義は ../八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py の omul（ケーリー＝ディクソン構成、
基底 e0=1, e1..e3 = i,j,k, e4 = ℓ, e5..e7 = iℓ, jℓ, kℓ）を、記号（sympy）にそのまま適用して展開する。
出力：(1) 基底の掛け算表、(2) x·y の 8 成分、(3) ō_e·o_f の 8 成分、(4) 案2 の一項 Im(ō_e·o_f)·o_f の 8 成分、
（3）（4）は o_e = (a0..a7), o_f = (b0..b7)。結果は results/expansion.md に保存し、数値の omul と照合する。
"""
import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')
src = open(SRC).read().split("results = {")[0]
ns = {'__file__': __file__}
exec(src, ns)                                # 数値の omul（照合用）


def qmul(a, b):                              # 元ファイルの qmul と同じ式（記号用）
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return [w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2]


def qconj(a):
    return [a[0], -a[1], -a[2], -a[3]]


def omul_s(o1, o2):                          # 元ファイルの omul と同じ式：(a,b)(c,d) = (ac − d̄b, da + bc̄)
    a, b = o1[:4], o1[4:]; c, d = o2[:4], o2[4:]
    p = qmul(a, c); q = qmul(qconj(d), b); r = qmul(d, a); s = qmul(b, qconj(c))
    return [sp.expand(p[k] - q[k]) for k in range(4)] + [sp.expand(r[k] + s[k]) for k in range(4)]


def oconj_s(o):
    return [o[0]] + [-x for x in o[1:]]


a = sp.symbols('a0:8'); b = sp.symbols('b0:8')
names = ['1', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7']
E = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
lines = ['# 八元数の積の 8 成分展開（omul と同じ定義）', '']
lines += ['## (1) 基底の掛け算表（行 × 列）', '', '| × | ' + ' | '.join(names) + ' |', '|' + '---|' * 9]
for i in range(8):
    row = []
    for j in range(8):
        p = omul_s(E[i], E[j]); k = [n for n in range(8) if p[n] != 0][0]
        row.append(('−' if p[k] < 0 else '') + names[k])
    lines.append(f'| **{names[i]}** | ' + ' | '.join(row) + ' |')
xy = omul_s(list(a), list(b))
lines += ['', '## (2) x·y の 8 成分（x = (a0..a7), y = (b0..b7)）', '']
for k in range(8):
    lines.append(f'- 成分 {k}（{names[k]}）： {sp.sstr(xy[k])}')
w = omul_s(oconj_s(list(a)), list(b))
lines += ['', '## (3) 関係 w = ō_e·o_f の 8 成分（o_e = (a0..a7), o_f = (b0..b7)）', '']
for k in range(8):
    lines.append(f'- w{k}（{names[k]}）： {sp.sstr(w[k])}')
s = [0] + w[1:]
t = omul_s(s, list(b))
lines += ['', '## (4) 案2の一項 s·o_f、s = Im(w) = (0, w1..w7)', '',
          '各成分を w1..w7 と b0..b7 で書く（w を代入する前）：', '']
W = sp.symbols('w1:8')
t_w = omul_s([0] + list(W), list(b))
for k in range(8):
    lines.append(f'- 成分 {k}（{names[k]}）： {sp.sstr(t_w[k])}')
# 数値照合
rng = np.random.default_rng(0); err = 0.0
for _ in range(50):
    x, y = rng.normal(size=(2, 8)); sub = dict(zip(a, x)); sub.update(zip(b, y))
    err = max(err, float(np.max(np.abs(np.array([float(e.subs(sub)) for e in xy]) - ns['omul'](x, y)))))
lines += ['', f'数値照合：記号展開と omul の差の最大 = {err:.1e}（50 組）']
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
open(os.path.join(HERE, 'results', 'expansion.md'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
