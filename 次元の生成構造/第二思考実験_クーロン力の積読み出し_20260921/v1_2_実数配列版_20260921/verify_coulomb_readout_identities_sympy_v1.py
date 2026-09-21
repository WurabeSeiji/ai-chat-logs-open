#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_coulomb_readout_identities_sympy_v1.py

第二思考実験「クーロン型の逆二乗の力を、積の読み出しと面積の時計から導く」付録 A（A.13〜A.22）の
恒等式を記号計算で確認する。各行の右側が 0（または零行列）になれば、その恒等式は成り立っている。

第一思考実験の verify_area_readout_identities_sympy_v1.py と同じ体裁。

依存: sympy
実行: python3 verify_coulomb_readout_identities_sympy_v1.py
"""
import sympy as sp

a, b, c = sp.symbols("a b c", real=True)
p, q, r = sp.symbols("p q r", real=True)          # S の成分（s22 は det S=1 から決まる）
n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)    # 生成子 N（トレースゼロ）

Om = sp.Matrix([[0, 1], [-1, 0]])
P = sp.Matrix([[0, 1], [1, 0]])
D = sp.Matrix([[-1, 0], [0, 1]])
J0 = D * P
I2 = sp.eye(2)
X = sp.Matrix([a, b])


def omega(U, V):
    """面積形式 ω(U,V)。"""
    return (U.T * Om * V)[0]


def sym(M):
    return (M + M.T) / 2


def readout(U):
    """三つの面積読み出し (x, y, r)。"""
    return sp.Matrix([omega(U, P * U), omega(U, D * U), omega(U, J0 * U)])


def B(U):
    """B(U)=u1 I+u2 J0。実二次写像 Phi(U)=B(U)U を作る。"""
    return sp.Matrix([[U[0], -U[1]], [U[1], U[0]]])


def bilin(U, V):
    return B(U) * V


def Phi(U):
    return bilin(U, U)


print("=" * 78)
print("A.13（T13）読み出しの空間は三次元、核は cI に限る")
print("=" * 78)

m11, m12, m21, m22 = sp.symbols("m11 m12 m21 m22", real=True)
M = sp.Matrix([[m11, m12], [m21, m22]])
print("T13 omega(X,MX) - X^T sym(OmM) X                  :", sp.simplify(omega(X, M * X) - (X.T * sym(Om * M) * X)[0]))
print("T13 核: sym(OmM)=0 の解 M - m11 I             :",
      sp.simplify(M.subs(sp.solve([e for e in sym(Om * M)], [m12, m21, m22], dict=True)[0]) - m11 * I2))
print("T13 核: omega(X, cI X)                            :", sp.simplify(omega(X, c * I2 * X)))
print("T13 J0 - D P                                  :", sp.simplify(J0 - D * P))
print("T13 P^2 - I                                   :", sp.simplify(P * P - I2))
print("T13 D^2 - I                                   :", sp.simplify(D * D - I2))
print("T13 J0^2 + I                                  :", sp.simplify(J0 * J0 + I2))
print("T13 P D - (-J0)                               :", sp.simplify(P * D + J0))
print("T13 J0 + Om                                   :", sp.simplify(J0 + Om))
x_, y_, r_ = readout(X)
print("T13 omega(X,PX) - (a^2-b^2)                       :", sp.simplify(x_ - (a**2 - b**2)))
print("T13 omega(X,DX) - 2ab                             :", sp.simplify(y_ - 2 * a * b))
print("T13 omega(X,J0X) - (a^2+b^2)                      :", sp.simplify(r_ - (a**2 + b**2)))
AA, BB, CC = sp.symbols("A B C", real=True)
print("T13 像は二次形式の全体: (A-C)/2 x + B y + (A+C)/2 r - (A a^2+2B ab+C b^2) :",
      sp.expand((AA - CC) / 2 * x_ + BB * y_ + (AA + CC) / 2 * r_ - (AA * a**2 + 2 * BB * a * b + CC * b**2)))

print()
print("=" * 78)
print("A.14（T14）距離はゼロ閉塞として出る")
print("=" * 78)

print("T14 x^2 + y^2 - r^2                           :", sp.simplify(x_**2 + y_**2 - r_**2))
print("T14 B(X) - (a I + b J0)                       :", sp.simplify(B(X) - (a * I2 + b * J0)))
print("T14 Phi(X) - (x,y)                            :", sp.simplify(Phi(X) - sp.Matrix([x_, y_])))
print("T14 B(X)^T B(X) - r I                         :", sp.simplify(B(X).T * B(X) - r_ * I2))
print("T14 Phi(X).Phi(X) - r^2                       :", sp.simplify((Phi(X).T * Phi(X))[0] - r_**2))
print("T14 readout(PX) - (-x, y, r)                  :", sp.simplify(readout(P * X) - sp.Matrix([-x_, y_, r_])))
print("T14 readout(-X) - (x, y, r)                   :", sp.simplify(readout(-X) - sp.Matrix([x_, y_, r_])))
cc, dd = sp.symbols("cc dd", real=True)
Y = sp.Matrix([cc, dd])
print("T14 omega(PX,PY) + omega(X,Y)                         :", sp.simplify(omega(P * X, P * Y) + omega(X, Y)))
print("T14 P P P - P                                 :", sp.simplify(P * P * P - P))
print("T14 P D P + D                                 :", sp.simplify(P * D * P + D))
print("T14 P J0 P + J0                               :", sp.simplify(P * J0 * P + J0))

print()
print("=" * 78)
print("A.15（T15）掃く面積は積の読み出し")
print("=" * 78)

a0, b0, a1, b1 = sp.symbols("a0 b0 a1 b1", real=True)
X0, X1 = sp.Matrix([a0, b0]), sp.Matrix([a1, b1])
p_dot = (X0.T * X1)[0]
q_om = omega(X0, X1)
Y0, Y1 = Phi(X0), Phi(X1)
print("T15 omega(Phi(X0),Phi(X1))/2 - q p            :",
      sp.simplify(omega(Y0, Y1) / 2 - q_om * p_dot))
# 連続：X' = NX（トレースゼロ）。dA/ds = q r
N = sp.Matrix([[n1, n2], [n3, -n1]])
dX = N * X
Y = Phi(X)
dY = 2 * bilin(X, dX)
print("T15 dPhi/ds - 2 B(X)X'                         :",
      sp.simplify(Y.jacobian(X) * dX - dY))
print("T15 omega(Phi,dPhi/ds)/2 - q r                :",
      sp.simplify(omega(Y, dY) / 2 - omega(X, dX) * r_))

print()
print("=" * 78)
print("A.16（T16）逆二乗の力。dt = r ds。σ = -det N（楕円 -1 / 放物 0 / 双曲 +1 に規格化）")
print("=" * 78)

sig = -N.det()                                   # X'' = N^2 X = sigma X
print("T16 N^2 - sigma I                             :", sp.simplify(N * N - sig * I2))
print("T16 sigma - (n1^2 + n2 n3)                    :", sp.simplify(sig - (n1**2 + n2 * n3)))
mu = 2 * (dX.T * dX)[0] - 2 * sig * (X.T * X)[0]
E = 2 * sig
Y = Phi(X)
dY = 2 * bilin(X, dX)
ddY_s = 2 * Phi(dX) + 2 * sig * Y
rp = 2 * (X.T * dX)[0]
Yddot = sp.simplify((ddY_s * r_ - dY * rp) / r_**3)
print("T16 d^2Phi/dt^2 + mu Phi/r^3                  :", sp.simplify(Yddot + mu * Y / r_**3))
print("T16 mu - (2|X'|^2 - E|X|^2)                   :", sp.simplify(mu - (2 * (dX.T * dX)[0] - E * (X.T * X)[0])))
dmu = 4 * ((dX.T * (N * dX))[0]) - 4 * sig * (X.T * dX)[0]
print("T16 dmu/ds                                    :", sp.simplify(dmu))
Ydot = dY / r_
print("T16 Kepler energy |Ydot|^2/2 - mu/r - E       :",
      sp.simplify((Ydot.T * Ydot)[0] / 2 - mu / r_ - E))
print("T16 angular momentum omega(Phi,Ydot) - 2 q    :",
      sp.simplify(omega(Y, Ydot) - 2 * omega(X, dX)))
# sigma = n1^2 + n2 n3 を所定の値に固定する（n3 を解く）
sub_ell = {n3: (-1 - n1**2) / n2}
sub_hyp = {n3: (1 - n1**2) / n2}
sub_par = {n3: -n1**2 / n2}
print("T16 楕円型 N^2 = -I のとき E + 2               :", sp.simplify((E + 2).subs(sub_ell)))
print("T16 双曲型 N^2 = +I のとき E - 2               :", sp.simplify((E - 2).subs(sub_hyp)))
print("T16 放物型 det N = 0 のとき E                  :", sp.simplify(E.subs(sub_par)))

print()
print("=" * 78)
print("A.17（T17）円錐曲線は第一思考実験の保存量")
print("=" * 78)

S = sp.Matrix([[p, q], [r, (1 + q * r) / p]])     # det S = 1
tau = S.trace()
G = sp.simplify(Om * (S - S.inv()) / 2)
A_, B_, C_ = G[0, 0], G[0, 1], G[1, 1]
print("T17 det S - 1                                 :", sp.simplify(S.det() - 1))
print("T17 G - G^T                                   :", sp.simplify(G - G.T))
print("T17 X^T G X - omega(X,SX)                         :", sp.simplify((X.T * G * X)[0] - omega(X, S * X)))
print("T17 (A, B, C) - (r_S, (s-p)/2, -q_S)          :",
      sp.simplify(sp.Matrix([A_, B_, C_]) - sp.Matrix([r, (S[1, 1] - p) / 2, -q])))
print("T17 det G - (1 - tau^2/4)                     :", sp.simplify(G.det() - (1 - tau**2 / 4)))
qv = (X.T * G * X)[0]
print("T17 q - [(A-C)/2 x + B y + (A+C)/2 r]         :",
      sp.simplify(qv - ((A_ - C_) / 2 * x_ + B_ * y_ + (A_ + C_) / 2 * r_)))
pp_ = 2 * qv / (A_ + C_)
ee = sp.Matrix([(A_ - C_) / (A_ + C_), 2 * B_ / (A_ + C_)])
print("T17 r - [p - e.(x,y)]  （焦点が原点の円錐曲線）:",
      sp.simplify(r_ - (pp_ - (ee[0] * x_ + ee[1] * y_))))
print("T17 |e|^2 - (1 - 4 det G/(tr G)^2)            :",
      sp.simplify((ee.T * ee)[0] - (1 - 4 * G.det() / G.trace()**2)))
print("T17 ((A-C)^2+4B^2)/(A+C)^2 - (1-4(AC-B^2)/(A+C)^2) :",
      sp.simplify(((AA - CC)**2 + 4 * BB**2) / (AA + CC)**2 - (1 - 4 * (AA * CC - BB**2) / (AA + CC)**2)))

print()
print("=" * 78)
print("A.17（T17）Kepler の式との一致（生成子 N で書く。放物型 σ=0 も含む）")
print("=" * 78)

Gc = Om * N                                       # 連続版 G（離散の G = sin(Δs)・Gc、e は尺度不変）
qc = (X.T * Gc * X)[0]
print("T17 Gc - Gc^T                                 :", sp.simplify(Gc - Gc.T))
print("T17 X^T Gc X - omega(X,NX)                        :", sp.simplify(qc - omega(X, dX)))
print("T17 det Gc + sigma                            :", sp.simplify(Gc.det() + sig))
print("T17 mu - 2 q tr(Gc)   ★恒等式                 :", sp.simplify(mu - 2 * qc * Gc.trace()))
L = 2 * omega(X, dX)
print("T17 (1-4detGc/(trGc)^2) - (1+2 E L^2/mu^2)    :",
      sp.simplify((1 - 4 * Gc.det() / Gc.trace()**2) - (1 + 2 * E * L**2 / mu**2)))
print("T17 離散 G と連続 Gc の離心率の差（尺度不変）  :",
      sp.simplify((1 - 4 * (sp.Symbol('k', positive=True) * Gc).det() / (sp.Symbol('k', positive=True) * Gc).trace()**2)
                  - (1 - 4 * Gc.det() / Gc.trace()**2)))

print()
print("=" * 78)
print("A.18（T18）結合の符号は区画の開き角で決まる")
print("=" * 78)

f, g, al, be = sp.symbols("f g alpha beta", real=True)
vp = sp.Matrix([sp.cos(f), sp.sin(f)])            # 零方向（単位ベクトル）
vm = sp.Matrix([sp.cos(g), sp.sin(g)])
Ah = sp.simplify(sp.Matrix.hstack(vp, -vm) * sp.Matrix.hstack(vp, vm).inv())
print("T18 A vp - vp                                 :", sp.simplify(Ah * vp - vp))
print("T18 A vm + vm                                 :", sp.simplify(Ah * vm + vm))
print("T18 A^2 - I                                   :", sp.simplify(Ah * Ah - I2))
print("T18 tr A                                      :", sp.simplify(Ah.trace()))
Xh = al * vp + be * vm
dXh = Ah * Xh
print("T18 X' - (alpha vp - beta vm)                 :", sp.simplify(dXh - (al * vp - be * vm)))
muh = 2 * (dXh.T * dXh)[0] - 2 * (Xh.T * Xh)[0]   # 双曲型 E=+2
print("T18 mu + 8 alpha beta (vp.vm)                 :", sp.simplify(muh + 8 * al * be * (vp.T * vm)[0]))
print("T18 mu + 8 alpha beta cos(psi)                :", sp.simplify(muh + 8 * al * be * sp.cos(f - g)))
print("T18 vp.vm - cos(psi)                          :", sp.simplify((vp.T * vm)[0] - sp.cos(f - g)))
print("T18 名前の付け替え P で mu 不変（PX, PAP）    :",
      sp.simplify(2 * ((P * Ah * P * P * Xh).T * (P * Ah * P * P * Xh))[0]
                  - 2 * ((P * Xh).T * (P * Xh))[0] - muh))
print("T18 楕円型 mu = 2|X'|^2 + 2|X|^2 > 0 の確認（σ=-1 で mu - 2|X'|^2 - 2|X|^2）:",
      sp.simplify((mu - 2 * (dX.T * dX)[0] - 2 * (X.T * X)[0]).subs(sub_ell)))

print()
print("=" * 78)
print("A.19（T19）中性：四つの同値")
print("=" * 78)

print("T19 tr G - (r_S - q_S)                        :", sp.simplify(G.trace() - (r - q)))
print("T19 S 対称（q_S=r_S）のとき tr G              :", sp.simplify(G.trace().subs(r, q)))
Ao = sp.simplify(Ah.subs(g, f + sp.pi / 2))       # vp ⊥ vm
print("T19 (vp⊥vm) A - A^T                           :", sp.simplify(Ao - Ao.T))
print("T19 (vp⊥vm) tr(Om A)                          :", sp.simplify((Om * Ao).trace()))
print("T19 (vp⊥vm) mu（すべての配置で 0）             :", sp.simplify(muh.subs(g, f + sp.pi / 2)))
print("T19 (vp⊥vm) 読み出しの加速度 Phi¨（mu=0）      :",
      sp.simplify(Yddot.subs({n1: Ao[0, 0], n2: Ao[0, 1], n3: Ao[1, 0]})))
H = sp.symbols("H", real=True)
Ssym = sp.cosh(H) * I2 + sp.sinh(H) * Ao
print("T19 S=cosh(H)I+sinh(H)A は対称（S - S^T）      :", sp.simplify(Ssym - Ssym.T))
print("T19 det S - 1                                 :", sp.simplify(Ssym.det() - 1))
print("T19 tr G(S) （= 0）                            :", sp.simplify((Om * (Ssym - Ssym.inv()) / 2).trace()))
print("T19 mu = 0 ⟺ tr Gc = 0（mu - 2 q tr Gc の再掲）:", sp.simplify(mu - 2 * qc * Gc.trace()))
lam = sp.symbols("lambda", positive=True)
print("T19 対称な実行列は |tau|>=2：(lam+1/lam)^2-4-(lam-1/lam)^2 :",
      sp.simplify((lam + 1 / lam)**2 - 4 - (lam - 1 / lam)**2))

print()
print("=" * 78)
print("A.20（T20）衝突軌道は半直線の上")
print("=" * 78)

c0, c1, tau_ = sp.symbols("c0 c1 tau", real=True)
u = sp.Matrix([a, b])
Xk0, Xk1 = c0 * u, c1 * u
Xk2 = tau_ * Xk1 - Xk0
print("T20 omega(X0,X1)                                  :", sp.simplify(omega(Xk0, Xk1)))
print("T20 X2 が同じ直線上（omega(u, X2)）                :", sp.simplify(omega(u, Xk2)))
print("T20 omega(X1,X2)                                  :", sp.simplify(omega(Xk1, Xk2)))
ck = sp.symbols("c_k", real=True)
print("T20 Phi(c_k u) - c_k^2 Phi(u)                  :",
      sp.simplify(Phi(ck * u) - ck**2 * Phi(u)))
print("T20 omega(Phi(c_k u),Phi(u))（同じ半直線）     :",
      sp.simplify(omega(Phi(ck * u), Phi(u))))
print("T20 L = 2q = 0（角運動量ゼロ）                 :", sp.simplify(2 * omega(Xk0, Xk1)))
print("T20 S 法則で X0∥X1（固有値 lam）なら特性方程式 lam^2 - tau lam + 1 に tau=lam+1/lam :",
      sp.simplify((lam**2 - tau_ * lam + 1).subs(tau_, lam + 1 / lam)))
print("T20 その tau は |tau|>=2：(lam+1/lam)^2 - 4 - (lam-1/lam)^2 :",
      sp.simplify((lam + 1 / lam)**2 - 4 - (lam - 1 / lam)**2))
print("T20 E6 の設定（X1=(7/10)X0）が S 法則に要求する tau - 149/70 :",
      sp.nsimplify(sp.Rational(7, 10) + sp.Rational(10, 7) - sp.Rational(149, 70)))
print("T20 E6 が使う tau=17/10 との差 149/70 - 17/10 - 30/70（0 でない差の確認）:",
      sp.nsimplify(sp.Rational(149, 70) - sp.Rational(17, 10) - sp.Rational(30, 70)))

print()
print("=" * 78)
print("§3.1・§7 時計の規格化と刻みの式（T15 の帰結）")
print("=" * 78)

th, e_, u0 = sp.symbols("theta e u0", real=True)
ds = th                                            # 楕円型 Δs = arccos(tau/2), tau = 2 cos θ
tau_e = 2 * sp.cos(th)
print("clock Δs の規格化 (2Δs/tau) - θ/cosθ          :", sp.simplify(2 * ds / tau_e - th / sp.cos(th)))
T = sp.Matrix([[sp.sqrt(1 - e_), 0], [0, sp.sqrt(1 + e_)]])   # w-楕円（a_K = 1）
Rt = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
Se = sp.simplify(T * Rt * T.inv())
print("clock tr S - 2 cos θ                          :", sp.simplify(Se.trace() - tau_e))
print("clock det S - 1                               :", sp.simplify(Se.det() - 1))
Xu = lambda uu: T * sp.Matrix([sp.cos(uu), sp.sin(uu)])
Xk, Xk1_ = Xu(u0), Xu(u0 + th)
print("clock X_k.X_{k+1} - (cosθ - e cos(2u_k+θ))    :",
      sp.simplify((Xk.T * Xk1_)[0] - (sp.cos(th) - e_ * sp.cos(2 * u0 + th))))
print("clock r_k - (1 - e cos 2u_k)  （離心近点角=2u）:", sp.simplify((Xk.T * Xk)[0] - (1 - e_ * sp.cos(2 * u0))))
Ge = sp.simplify(Om * (Se - Se.inv()) / 2)
print("clock e_metric^2 - e^2  ★e は tau と独立      :",
      sp.simplify((1 - 4 * Ge.det() / Ge.trace()**2) - e_**2))
print("clock tr G - 2 sinθ/sqrt(1-e^2)  ★tau に依らぬ自由度 :",
      sp.simplify(Ge.trace() - 2 * sp.sin(th) / (sp.sqrt(1 - e_) * sp.sqrt(1 + e_))))
print("clock det G - (1 - tau^2/4)                   :", sp.simplify(Ge.det() - (1 - tau_e**2 / 4)))
Ae = sp.simplify(Se - sp.cos(th) * I2) / sp.sin(th)
mue = sp.simplify(2 * ((Ae * Xk).T * (Ae * Xk))[0] + 2 * (Xk.T * Xk)[0])
print("clock mu - 4  （a_K=1 の規格化）               :", sp.simplify(mue - 4))
for nn in (5, 6, 8):
    thn = sp.Rational(2, 1) * sp.pi / nn
    tot = sum((Xu(u0 + k * thn).T * Xu(u0 + (k + 1) * thn))[0] for k in range(nn))
    dtsum = sp.simplify((thn / sp.cos(thn)) * tot)
    print(f"clock n={nn}: Σ Δt_k - 2π mu/4 (= 2π)           :", sp.simplify(dtsum - 2 * sp.pi * 4 / 4))

print()
print("=" * 78)
print("A.21（T21）自由度の勘定と二体の結合")
print("=" * 78)

m1, m2, Gg, kc, s1, s2, s3 = sp.symbols("m1 m2 G_grav k_coul sigma1 sigma2 sigma3", real=True)
q1, q2 = s1 * m1, s2 * m2
F1 = -Gg * m1 * m2 + kc * q1 * q2                  # 相対座標の動径成分（引力を負に取る）/ r^2
a_rel = sp.simplify(F1 / m1 - (-F1) / m2)
print("T21 相対加速度の係数 + (m1+m2)(G - k s1 s2)   :", sp.simplify(a_rel + (m1 + m2) * (Gg - kc * s1 * s2)))
cik = m2 * (Gg - kc * s1 * s2)
cki = m1 * (Gg - kc * s1 * s2)
print("T21 c_ik m_i - c_ki m_k                       :", sp.simplify(cik * m1 - cki * m2))
print("T21 s12 s23 s31 - (s1 s2 s3)^2                :",
      sp.simplify((s1 * s2) * (s2 * s3) * (s3 * s1) - (s1 * s2 * s3)**2))
print("T21 形の自由度 2N-4: N=2 で 0                  :", sp.simplify(2 * 2 - 4 - 0))
print("T21 形の自由度 2N-4: N=3 で 2                  :", sp.simplify(2 * 3 - 4 - 2))
print("T21 モデルの量の数 (tau 1 + G の形 2 + X0 2) - 5 :", sp.simplify(1 + 2 + 2 - 5))
print("T21 連続軌道を決める数 (N 2 + X0 2) - 4        :", sp.simplify(2 + 2 - 4))

print()
print("=" * 78)
print("A.22（T22）D を一般の反転 D'_phi に取り替える")
print("=" * 78)
ph = sp.symbols("phi", real=True)
cph, sph = sp.cos(2 * ph), sp.sin(2 * ph)
Dg = sp.Matrix([[cph, sph], [sph, -cph]])
print("T22 D'^2 - I                                  :", sp.simplify(Dg * Dg - I2))
print("T22 phi=pi/2 で D' - D                        :", sp.simplify(Dg.subs(ph, sp.pi / 2) - D))
print("T22 omega(X,D'X) - (sin2phi x - cos2phi y)        :", sp.simplify(omega(X, Dg * X) - (sph * x_ - cph * y_)))
print("T22 omega(X,D'PX) + cos2phi r                     :", sp.simplify(omega(X, Dg * P * X) + cph * r_))
print("T22 基底の行列式 - cos^2(2phi)                 :", sp.simplify(sp.Matrix([[1, 0, 0], [sph, -cph, 0], [0, 0, -cph]]).det() - cph**2))
print("T22 r(D'X) - r                                :", sp.simplify(omega(Dg * X, J0 * Dg * X) - r_))
uu, vv, ww = sp.symbols("uu vv ww", real=True)
Qg = uu * x_ + vv * y_ + ww * r_
def co(e):
    e = sp.expand(e)
    return [e.coeff(a, 2).coeff(b, 0), e.coeff(a, 1).coeff(b, 1), e.coeff(a, 0).coeff(b, 2)]
QgP = Qg.subs({a: b, b: a}, simultaneous=True)
QgD = Qg.subs({a: cph * a + sph * b, b: sph * a - cph * b}, simultaneous=True)
solg = sp.solve([sp.simplify(t) for t in co(Qg - QgP) + co(Qg - QgD)], [uu, vv], dict=True)
print("T22 <P,D'> 不変な読み出しの (u, v)            :", sp.Matrix([solg[0][uu], solg[0][vv]]))
R4 = sp.Matrix([[sp.cos(4 * ph), sp.sin(4 * ph)],
                [sp.sin(4 * ph), -sp.cos(4 * ph)]])
print("T22 Phi(D'X) - R4 Phi(X)                      :", sp.simplify(Phi(Dg * X) - R4 * Phi(X)))
print("T22 R4^T R4 - I                               :", sp.simplify(R4.T * R4 - I2))
print("T22 det R4 + 1                                :", sp.simplify(R4.det() + 1))
