"""600 胞体（二項正二十面体群 2I）上の住民模型　設計書 v1 の共通部品。

左の平行移動（ゲージ）：群の左の掛け算 P_g と、その生成子 X_u。
運動の生成子：右の平行移動の生成子 Y_i（左の P_g と厳密に可換）。
右の交換（背景）：成分を回す M（実反対称）と、その一様な強さ h。
"""
import itertools, warnings
import numpy as np
warnings.filterwarnings('ignore')

PHI = (1 + 5 ** 0.5) / 2
NMAX = 5            # 分解能：5 次までの球面調和関数
STEP = np.pi / 5    # ゲージの一歩（最近接の角 36°）


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([a0*b0 - a1*b1 - a2*b2 - a3*b3,
                     a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1,
                     a0*b3 + a1*b2 - a2*b1 + a3*b0])


def vertices():
    V = []
    for i in range(4):
        for s in (1, -1):
            v = np.zeros(4); v[i] = s; V.append(v)
    for s in itertools.product((1, -1), repeat=4):
        V.append(np.array(s) / 2)
    even = [(0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2), (1, 0, 3, 2), (1, 2, 0, 3), (1, 3, 2, 0),
            (2, 0, 1, 3), (2, 1, 3, 0), (2, 3, 0, 1), (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0)]
    for p in even:
        for s in itertools.product((1, -1), repeat=3):
            base = np.array([PHI / 2 * s[0], 0.5 * s[1], 1 / (2 * PHI) * s[2], 0.0])
            v = np.zeros(4)
            for k in range(4):
                v[p[k]] = base[k]
            V.append(v)
    return np.array(V)


V = vertices()
NV = len(V)
GRAM = np.clip(V @ V.T, -1, 1)


def index_of(q):
    d = np.linalg.norm(V - q, axis=1)
    i = int(d.argmin())
    assert d[i] < 1e-9
    return i


def perm_left(g):
    """(P_g f)(q) = f(g^{-1} q)：左の掛け算による頂点の置換（クリフォード平行移動）。"""
    P = np.zeros((NV, NV))
    for i, q in enumerate(V):
        P[index_of(qmul(g, q)), i] = 1
    return P


def cheb_u(n, x):
    """第二種 Chebyshev 多項式 U_n とその微分（漸化式で厳密）。"""
    u0 = np.ones_like(x); d0 = np.zeros_like(x)
    if n == 0:
        return u0, d0
    u1 = 2 * x; d1 = 2 * np.ones_like(x)
    for _ in range(1, n):
        u0, u1 = u1, 2 * x * u1 - u0
        d0, d1 = d1, 2 * u0 + 2 * x * d1 - d0
    return u1, d1


def degree_bases():
    out = []
    for n in range(NMAX + 1):
        B = cheb_u(n, GRAM)[0]
        u_, s_, _ = np.linalg.svd(B)
        r = int((s_ > 1e-8 * s_[0]).sum())
        out.append(u_[:, :r])
    return out


BASES = degree_bases()
Q = np.hstack(BASES)                 # 91 列、正規直交
PI = Q @ Q.T                         # H5 への射影
LAP = sum(-n * (n + 2) * (B @ B.T) for n, B in enumerate(BASES))
LAP_INV = sum((-1.0 / (n * (n + 2))) * (B @ B.T) for n, B in enumerate(BASES) if n > 0)


def generator(u, side):
    """H5 上の平行移動の生成子。side='L'：f(e^{-u e} q) の微分、side='R'：f(q e^{u e}) の微分。"""
    uq = np.array([0.0, *u])
    if side == 'L':
        dq = np.array([-qmul(uq, q) for q in V])
    else:
        dq = np.array([qmul(q, uq) for q in V])
    X = np.zeros((NV, NV))
    tang = dq @ V.T
    for n, Bq in enumerate(BASES):
        B, dB = cheb_u(n, GRAM)
        X += (dB * tang) @ np.linalg.pinv(B) @ (Bq @ Bq.T)
    return X


E3 = np.eye(3)
XL = [generator(e, 'L') for e in E3]
YR = [generator(e, 'R') for e in E3]

_sx = np.array([[0., 1.], [1., 0.]]); _sz = np.array([[1., 0.], [0., -1.]])
_isy = np.array([[0., 1.], [-1., 0.]]); _I2 = np.eye(2)
A = [np.kron(_sx, _I2), np.kron(_sz, _I2), np.kron(_isy, _isy)]   # 実対称・反可換
M = np.kron(_isy, _sx)                                              # 実反対称・M^2=-1
PIF = np.kron(np.eye(4), PI)


def h_from_theta(theta):
    """右の交換の一歩あたりの角 θ（ゲージの一歩 36° の間）から強さ h。"""
    return theta / STEP


THETA_E = 23 * np.pi / 124
H_E = h_from_theta(THETA_E)


def stretch(Phi):
    N = 1 + Phi + Phi ** 2 / 2
    s = 1 - Phi + Phi ** 2 / 4
    return N, s


def K_of(h, Phi=None):
    if Phi is None:
        Phi = np.zeros(NV)
    N, s = stretch(Phi)
    W = np.diag(np.sqrt(N / s))
    K = sum(np.kron(a, PI @ W @ Y @ W @ PI) for a, Y in zip(A, YR))
    K = K + np.kron(M, PI @ np.diag(h * N) @ PI)
    return K


def mass_density(psi, h):
    P = psi.reshape(4, NV)
    return h * np.real(-1j * np.einsum('in,ij,jn->n', P.conj(), M, P))


def density(psi):
    return (np.abs(psi.reshape(4, NV)) ** 2).sum(0)


def solve_phi(mu, kappa):
    return kappa * (LAP_INV @ mu)


G10 = [v for v in V if abs(v[0] - PHI / 2) < 1e-9][0]     # 位数 10 の元（ゲージの一歩 36°）
U_AXIS = G10[1:] / np.linalg.norm(G10[1:])                  # 並進の方向（五回軸）
XU = sum(ui * X for ui, X in zip(U_AXIS, XL))


def joint_modes(h):
    """平坦な K と左の生成子 X_u（u は五回軸）の同時固有モード。戻り値：固有ベクトル、ω、k。"""
    K = K_of(h)
    Xu = np.kron(np.eye(4), XU)
    H = -1j * K
    P = -1j * Xu
    w, Vv = np.linalg.eigh(H + 1e-3 * np.pi * P)
    om = np.real(np.einsum('ij,ik,kj->j', Vv.conj(), H, Vv))
    kk = np.real(np.einsum('ij,ik,kj->j', Vv.conj(), P, Vv))
    return Vv, om, kk


def packet(h, k0, i0, dk=1.0, eta=None):
    """頂点 i0 に局在し、運動量 k0 のまわりに絞った、正の振動数だけの波束。"""
    if eta is None:
        eta = np.ones(4) / 2
    Vv, om, kk = joint_modes(h)
    sel = om > 1e-6
    loc = np.kron(eta, np.eye(NV)[i0])
    coef = (Vv[:, sel].conj().T @ loc) * np.exp(-(kk[sel] - k0) ** 2 / (2 * dk ** 2))
    psi = Vv[:, sel] @ coef
    nr = np.linalg.norm(psi)
    return psi / nr, dict(V=Vv[:, sel], om=om[sel], k=kk[sel], coef=coef / nr)


ZETA = V[:, 0] + 1j * (V[:, 1:] @ U_AXIS)   # u 方向の左の平行移動 e^{us} で e^{is} 倍になる関数


def position(psi):
    m = density(psi) @ V
    return m / np.linalg.norm(m)


def fibre_phase(psi):
    return np.angle(density(psi) @ ZETA)


# ---- v2：スピンを観測していない粒子（4 つのスピンの向きの等しい混合） ----

def packet_mixed(h, k0, i0, dk=1.0):
    """頂点 i0 のデルタ関数に η = e1..e4 を掛け、正の振動数・運動量 k0 のまわりに射影して等しく混ぜる。
    戻り値：成分のリスト（全体のノルムの二乗の和が 1）と、各成分のモード係数。"""
    Vv, om, kk = joint_modes(h)
    sel = om > 1e-6
    comps, coefs = [], []
    for e in np.eye(4):
        loc = np.kron(e, np.eye(NV)[i0])
        coef = (Vv[:, sel].conj().T @ loc) * np.exp(-(kk[sel] - k0) ** 2 / (2 * dk ** 2))
        comps.append(Vv[:, sel] @ coef); coefs.append(coef)
    tot = np.sqrt(sum(np.linalg.norm(c) ** 2 for c in comps))
    comps = [c / tot for c in comps]; coefs = [c / tot for c in coefs]
    return comps, dict(V=Vv[:, sel], om=om[sel], k=kk[sel], coefs=coefs)


def density_mixed(comps):
    return sum(density(c) for c in comps)


def mass_density_mixed(comps, h):
    return sum(mass_density(c, h) for c in comps)


def position_mixed(comps):
    mvec = density_mixed(comps) @ V
    return mvec / np.linalg.norm(mvec)


def fibre_phase_mixed(comps):
    return np.angle(density_mixed(comps) @ ZETA)


LAMBDA_WIDTH = 2.0   # 止まっている粒子の絞りの幅（v2 設計書 §2）


def rest_mixed(h, i0, lw=LAMBDA_WIDTH):
    """止まっている粒子：運動エネルギーの大きさ λ=√(ω²−h²) で等方に絞り、スピンの 4 向きを等しく混ぜる。"""
    Vv, om, kk = joint_modes(h)
    sel = om > 1e-6
    lam = np.sqrt(np.clip(om[sel] ** 2 - h ** 2, 0, None))
    filt = np.exp(-lam ** 2 / (2 * lw ** 2)) if lw else np.ones(int(sel.sum()))
    comps = [Vv[:, sel] @ ((Vv[:, sel].conj().T @ np.kron(e, np.eye(NV)[i0])) * filt) for e in np.eye(4)]
    tot = np.sqrt(sum(np.linalg.norm(c) ** 2 for c in comps))
    return [c / tot for c in comps]
