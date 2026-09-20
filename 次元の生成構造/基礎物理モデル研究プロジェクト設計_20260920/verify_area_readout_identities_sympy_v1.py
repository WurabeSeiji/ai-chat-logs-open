#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_area_readout_identities_sympy_v1.py

第一思考実験「等価原理を無名な二自由度まで遡る」付録 A の恒等式を記号計算で確認する。
各行の右側が 0（または零行列）になれば、その恒等式は成り立っている。

依存: sympy
実行: python3 verify_area_readout_identities_sympy_v1.py
"""
import sympy as sp

p, q, r, a, b, c, d = sp.symbols("p q r a b c d", real=True)
Om = sp.Matrix([[0, 1], [-1, 0]])
P = sp.Matrix([[0, 1], [1, 0]])
X = sp.Matrix([a, b])
Y = sp.Matrix([c, d])


def w(U, V):
    """面積形式 ω(U,V)。"""
    return (U.T * Om * V)[0]


# 一般の S∈SL(2,R)：s22 = (1 + q r)/p で det S = 1 を厳密に満たす
S = sp.Matrix([[p, q], [r, (1 + q * r) / p]])
tau = S.trace()
kap = 2 - tau
N = (S - S.inv()) / 2
G = sp.simplify(Om * N)

A11, A12, A21, A22 = sp.symbols("A11 A12 A21 A22", real=True)
A = sp.Matrix([[A11, A12], [A21, A22]])

print("det S - 1                              :", sp.simplify(S.det() - 1))
print("T2  w(AX,AY) - det A * w(X,Y)          :", sp.simplify(w(A * X, A * Y) - A.det() * w(X, Y)))
print("T3  S + S^-1 - tau I                   :", sp.simplify(S + S.inv() - tau * sp.eye(2)))
print("T5  Q_S(SX) - Q_S(X)                   :", sp.simplify(w(S * X, S * S * X) - w(X, S * X)))
print("T5  G - G^T                            :", sp.simplify(G - G.T))
print("T5  X^T G X - w(X,SX)                  :", sp.simplify((X.T * G * X)[0] - w(X, S * X)))
print("T5  det G - (1 - tau^2/4)              :", sp.simplify(G.det() - (1 - tau**2 / 4)))
print("T5  N^2 + (1 - tau^2/4) I              :", sp.simplify(N * N + (1 - tau**2 / 4) * sp.eye(2)))
dX = (S - sp.eye(2)) * X
d2X = (S - 2 * sp.eye(2) + S.inv()) * X
print("T6  Q_S(dX) - kappa Q_S(X)             :", sp.simplify((dX.T * G * dX)[0] - kap * (X.T * G * X)[0]))
print("T6  d2X + kappa X                      :", sp.simplify(d2X + kap * X))

# 正規形での読み出し
H, th = sp.symbols("H theta", real=True)
R = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
B = sp.Matrix([[sp.cosh(H), sp.sinh(H)], [sp.sinh(H), sp.cosh(H)]])
Sq = sp.Matrix([[sp.exp(H), 0], [0, sp.exp(-H)]])
Sh = sp.Matrix([[1, 0], [p, 1]])
print("回転      w(X,RX) - sin(theta)(a^2+b^2):", sp.simplify(w(X, R * X) - sp.sin(th) * (a**2 + b**2)))
print("ブースト  w(X,BX) - sinh(H)(a^2-b^2)   :", sp.simplify(w(X, B * X) - sp.sinh(H) * (a**2 - b**2)))
print("圧搾      w(X,SqX) + 2 sinh(H) a b     :", sp.simplify(w(X, Sq * X) + 2 * sp.sinh(H) * a * b))
print("剪断      w(X,ShX) - p a^2             :", sp.simplify(w(X, Sh * X) - p * a**2))
print("回転      P R P - R^-1                 :", sp.simplify(P * R * P - R.inv()))
print("ブースト  P B P - B                    :", sp.simplify(P * B * P - B))

# T11：設計書 §6 の骨格を N=2・線形 f(x,y)=pp x+qq y で実装したときの S
pp, qq = sp.symbols("pp qq", real=True)
f = lambda x, y: pp * x + qq * y
S2 = sp.Matrix([[sp.diff(f(a, a) + f(a, b), a), sp.diff(f(a, a) + f(a, b), b)],
                [sp.diff(f(b, b) + f(b, a), a), sp.diff(f(b, b) + f(b, a), b)]])
print("T11 S - [[2pp+qq, qq],[qq, 2pp+qq]]    :", sp.simplify(S2 - sp.Matrix([[2 * pp + qq, qq], [qq, 2 * pp + qq]])))
print("T11 det S - 4 pp (pp+qq)               :", sp.simplify(S2.det() - 4 * pp * (pp + qq)))

# T3′：det S=1 を仮定しない一般の可逆な S（d = det S）。Y_k = X_{-k} は列を逆向きに読んだもの
p2, q2, r2, s2 = sp.symbols("p2 q2 r2 s2", real=True)
Sg = sp.Matrix([[p2, q2], [r2, s2]])
tg, dg = Sg.trace(), Sg.det()
Ym1, Y0, Y1 = Sg * X, X, Sg.inv() * X
print("T3' 逆向きの法則  Y1 - (tau/d) Y0 + (1/d) Y_{-1} :", sp.simplify(Y1 - (tg / dg) * Y0 + Ym1 / dg))
res = Y1 - tg * Y0 + dg * Ym1
print("T3' 残差 - (d-1)/d [-tau Y0 + (d+1) Y_{-1}]      :", sp.simplify(res - (dg - 1) / dg * (-tg * Y0 + (dg + 1) * Ym1)))
print("T3' tr S^2 - (tau^2 - 2d)                        :", sp.simplify((Sg * Sg).trace() - (tg**2 - 2 * dg)))
Gd = Om * (Sg - dg * Sg.inv()) / 2
print("T3' G_d - G_d^T                                  :", sp.simplify(Gd - Gd.T))
print("T3' X^T G_d X - w(X,SX)                          :", sp.simplify((X.T * Gd * X)[0] - w(X, Sg * X)))
print("T3' det G_d - (d - tau^2/4)                      :", sp.simplify(Gd.det() - (dg - tg**2 / 4)))
print("T3' w(SX,S^2X) - d w(X,SX)                       :", sp.simplify(w(Sg * X, Sg * Sg * X) - dg * w(X, Sg * X)))
print("T3' w(X,PX) - (a^2 - b^2)                        :", sp.simplify(w(X, P * X) - (a**2 - b**2)))
