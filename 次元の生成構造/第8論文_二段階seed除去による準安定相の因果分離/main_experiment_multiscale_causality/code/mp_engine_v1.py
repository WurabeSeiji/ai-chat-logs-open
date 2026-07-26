#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 本実験 mpmath 高精度エンジン（Series 2）。解釈しない（力学の高精度移植のみ）。

第5論文 float64 エンジン run_n_scaling_lowrank_v1.LowRankSystem と同一アルゴリズムを
mpmath 任意精度で忠実移植する（更新則を変更しない）。σ は JG の固有分解で厳密算出（高精度版の
正規化。float64 版の warm-start 冪反復は σ_max の推定子であり、その厳密値を高精度で用いる）。
親は float64 make_parent で初期化後、高精度で自己無撞着精錬（同一親・高精度）。

観測（固有平面占有 q_j 等）も mpmath で算出し、最終スカラーのみ float64 で返す（sub-float64 の
振幅は値として float64 に格納可能）。
"""
import math

import mpmath as mp
import numpy as np

GAMMA_F64 = math.tan(math.pi / 144.0)


def edges(n):
    ea, eb = np.triu_indices(n, k=1)
    return list(zip(ea.tolist(), eb.tolist()))


class MPSystem:
    """位相差正弦生成子の低ランク表現（mpmath 版）。"""

    def __init__(self, n, prec_bits):
        self.n = n
        self.prec = prec_bits
        self.edges = edges(n)
        self.m = len(self.edges)
        with mp.workprec(prec_bits):
            self.gamma = mp.tan(mp.pi / 144)

    # ---- 位相設定: 辺角 theta(list mpf, len M) から c,s,G を構成 ----
    def set_theta(self, theta):
        n = self.n
        self.c = [mp.cos(t) for t in theta]
        self.s = [mp.sin(t) for t in theta]
        # n×n の辺角行列 T（対称, 対角0）
        T = [[mp.mpf(0)] * n for _ in range(n)]
        for e, (a, b) in enumerate(self.edges):
            T[a][b] = theta[e]; T[b][a] = theta[e]
        CT = [[mp.cos(T[i][j]) if i != j else mp.mpf(0) for j in range(n)] for i in range(n)]
        ST = [[mp.sin(T[i][j]) if i != j else mp.mpf(0) for j in range(n)] for i in range(n)]
        Gcc = [[CT[i][j] * CT[i][j] for j in range(n)] for i in range(n)]
        Gcs = [[CT[i][j] * ST[i][j] for j in range(n)] for i in range(n)]
        Gss = [[ST[i][j] * ST[i][j] for j in range(n)] for i in range(n)]
        for i in range(n):
            Gcc[i][i] = sum(Gcc[i]); Gcs[i][i] = sum(Gcs[i]); Gss[i][i] = sum(Gss[i])
        G = mp.zeros(2 * n)
        for i in range(n):
            for j in range(n):
                G[i, j] = Gcc[i][j]; G[i, n + j] = Gcs[i][j]
                G[n + i, j] = Gcs[i][j]; G[n + i, n + j] = Gss[i][j]
        self.G = G
        # J（2n×2n）
        J = mp.zeros(2 * n)
        for k in range(n):
            J[k, n + k] = mp.mpf(1); J[n + k, k] = mp.mpf(-1)
        self.J = J

    def _vsum(self, vals):
        """辺値(len M, 複素)を頂点(len n)へ散布加算。"""
        n = self.n
        out = [mp.mpc(0) for _ in range(n)]
        for e, (a, b) in enumerate(self.edges):
            out[a] += vals[e]; out[b] += vals[e]
        return out

    def wt(self, z):
        """W^T z（len 2n, 複素）"""
        cz = [self.c[e] * z[e] for e in range(self.m)]
        sz = [self.s[e] * z[e] for e in range(self.m)]
        vc = self._vsum(cz); vs = self._vsum(sz)
        return vc + vs

    def w(self, y):
        """W y（len M, 複素）。y = [yc(n)|ys(n)]"""
        n = self.n
        out = []
        for e, (a, b) in enumerate(self.edges):
            out.append(self.c[e] * (y[a] + y[b]) + self.s[e] * (y[n + a] + y[n + b]))
        return out

    def kmatvec(self, z):
        """K z = C(S^T z) - S(C^T z)（len M, 複素）"""
        sz = [self.s[e] * z[e] for e in range(self.m)]
        cz = [self.c[e] * z[e] for e in range(self.m)]
        vs = self._vsum(sz); vc = self._vsum(cz)
        out = []
        for e, (a, b) in enumerate(self.edges):
            out.append(self.c[e] * (vs[a] + vs[b]) - self.s[e] * (vc[a] + vc[b]))
        return out

    def jg_eig(self):
        """JG の固有値・固有ベクトル（mpmath eig）。戻り: (E list, EV matrix 2n×2n)。"""
        JG = self.J * self.G
        E, EV = mp.eig(JG)
        return E, EV

    def sigma_max(self):
        """σ_max = JG 固有値の最大正虚部。"""
        E, _ = self.jg_eig()
        ims = [mp.im(e) for e in E]
        return max(ims)

    def step_eig(self, kmax=8):
        """set_theta 済み前提で 1 回の eig から (σ_max, sigmas降順, Ws固有平面) を返す。

        1 step で σ_max（Cayley 用）と固有平面（観測用）を同一 eig から取得し計算量半減。
        """
        E, EV = self.jg_eig()
        items = [(mp.im(E[k]), k) for k in range(len(E))]
        items.sort(key=lambda x: x[0], reverse=True)
        n2 = 2 * self.n
        sig = []; Ws = []
        thr = mp.mpf(2) ** (-self.prec + 8)
        for im, k in items:
            if im <= thr:
                continue
            col = [EV[i, k] for i in range(n2)]
            w = self.w(col)
            nw = znorm(w)
            if nw == 0:
                continue
            if len(Ws) < kmax:
                w = [wi / nw for wi in w]
                Ws.append(w)
            sig.append(im)
        sigma_max = sig[0] if sig else mp.mpf(0)
        return sigma_max, sig, Ws

    def cayley_step(self, z, sigma):
        """z ← (I-γK̃)^{-1}(I+γK̃) z, K̃=K/σ。Woodbury（2n×2n 解, mp.lu_solve）。"""
        g = self.gamma
        gn = g / sigma
        Kz = self.kmatvec(z)
        r = [z[e] + gn * Kz[e] for e in range(self.m)]
        # A2 = (σ/γ)J + G
        A2 = mp.matrix(self.G.tolist())
        coef = sigma / g
        n2 = 2 * self.n
        for i in range(n2):
            for j in range(n2):
                A2[i, j] = A2[i, j] + coef * self.J[i, j]
        rhs = self.wt(r)
        rhs_v = mp.matrix(rhs)
        y = mp.lu_solve(A2, rhs_v)
        wy = self.w([y[i] for i in range(n2)])
        return [r[e] - wy[e] for e in range(self.m)]


# ---------------- 補助 ----------------

def zangle(z):
    """辺ごとの偏角（len M, mpf）。"""
    return [mp.atan2(mp.im(v), mp.re(v)) for v in z]


def znorm(z):
    return mp.sqrt(sum((mp.re(v) ** 2 + mp.im(v) ** 2) for v in z))


def zTz(z):
    """Z^T Z（複素, 零二乗閉鎖）。"""
    s = mp.mpc(0)
    for v in z:
        s += v * v
    return s


def cdot(a, b):
    """<a,b> = sum conj(a)·b（複素内積）。"""
    s = mp.mpc(0)
    for i in range(len(a)):
        s += mp.conj(a[i]) * b[i]
    return s


def f64(x):
    """mpf/mpc → float（実部）。桁溢れ回避のため mpf 経由。"""
    try:
        return float(mp.re(x))
    except (ValueError, OverflowError):
        return float("nan")


# ---------------- 親構成（高精度精錬） ----------------

def cast_parent(v_f64, prec_bits):
    """float64 親 v をそのまま高精度へキャスト（全系列で同一初期状態）。

    全系列 Series 1/2/3 が同一の Z0=v を共有する（Series 2 はビット等価なキャスト）ことで、
    軌道の分岐を純粋に演算精度の効果へ帰属させる。深い自己無撞着親の構成は不安定固定点＋σ縮退の
    ため Stage A では行わない（親自己無撞着残差 ~2e-13 = 全系列共通の初期 seed 床）。
    """
    with mp.workprec(prec_bits):
        return [mp.mpc(float(x.real), float(x.imag)) for x in v_f64]


def _eigenmode_at(sys_mp, theta, phase_ref_edge=0, phase_ref_val=None):
    """θ（len M）での σ_max 固有モード v（単位ノルム）。位相ゲージ: angle(v[ref])=phase_ref_val。"""
    sys_mp.set_theta(theta)
    E, EV = sys_mp.jg_eig()
    ims = [mp.im(E[k]) for k in range(len(E))]
    idx = min(range(len(ims)), key=lambda k: ims[k])   # λ=-iσ_max
    col = [EV[i, idx] for i in range(2 * sys_mp.n)]
    v = sys_mp.w(col)
    nv = znorm(v)
    v = [vi / nv for vi in v]
    if phase_ref_val is not None:
        # 大域位相を固定: v ← v·e^{-i(angle(v[ref])-phase_ref_val)}
        cur = mp.atan2(mp.im(v[phase_ref_edge]), mp.re(v[phase_ref_edge]))
        rot = mp.expjpi((phase_ref_val - cur) / mp.pi)
        v = [vi * rot for vi in v]
    return v


def refine_parent(sys_mp, theta0, iters=40, newton_iters=12):
    """高精度親精錬。減衰反復で近づけた後、Newton（位相ゲージ固定）で自己無撞着 θ*=angle(v(θ*))。

    自己無撞着写像は不安定固定点のため、Newton で F(θ_{1..M-1})=wrap(angle(v)_{1..M-1}-θ_{1..M-1})=0 を解く。
    ゲージ: θ[0] を固定し、v の大域位相を angle(v[0])=θ0[0] に合わせる。戻り: v, log10 residual。
    """
    n = sys_mp.n; M = sys_mp.m
    ref = 0
    th0_ref = mp.mpf(float(theta0[ref]))
    theta = [mp.mpf(float(t)) for t in theta0]

    def wrap(x):
        return mp.atan2(mp.sin(x), mp.cos(x))

    def residual_vec(theta_full):
        v = _eigenmode_at(sys_mp, theta_full, ref, th0_ref)
        ang = zangle(v)
        return [wrap(ang[e] - theta_full[e]) for e in range(M)], v

    # Levenberg–Marquardt on free components（θ0 近傍に留め別固定点への逸脱を防ぐ）
    free = list(range(M)); free.remove(ref)
    nf = len(free)
    h = mp.mpf(2) ** (-int(sys_mp.prec * 0.45))     # 有限差分ステップ
    STEP_CAP = mp.mpf("0.05")                        # 1反復の最大 |Δθ|（basin 逸脱防止）
    lam = mp.mpf("1e-6")
    F, v = residual_vec(theta)
    fnorm = max(abs(F[e]) for e in free)
    tol = mp.mpf(2) ** (-int(sys_mp.prec * 0.85))
    for _ in range(newton_iters):
        if fnorm < tol:
            break
        Ffree = mp.matrix([F[e] for e in free])
        J = mp.zeros(nf)
        for jj, ecol in enumerate(free):
            tp = list(theta); tp[ecol] = tp[ecol] + h
            Fp, _ = residual_vec(tp)
            for ii, erow in enumerate(free):
                J[ii, jj] = (Fp[erow] - F[erow]) / h
        JT = J.T
        JTJ = JT * J
        JTF = JT * Ffree
        accepted = False
        for _bt in range(30):                        # λ を増やしながら受理可能な step を探す
            A = mp.matrix(JTJ.tolist())
            for i in range(nf):
                A[i, i] = A[i, i] + lam * (JTJ[i, i] + mp.mpf(1))
            try:
                d = mp.lu_solve(A, JTF)
            except Exception:
                lam *= 10; continue
            dmax = max(abs(d[i]) for i in range(nf))
            scale = (STEP_CAP / dmax) if dmax > STEP_CAP else mp.mpf(1)
            tnew = list(theta)
            for jj, ecol in enumerate(free):
                tnew[ecol] = tnew[ecol] - scale * d[jj]
            Fn, vn = residual_vec(tnew)
            fnew = max(abs(Fn[e]) for e in free)
            if fnew < fnorm:                          # 残差減少なら受理, λ 減衰
                theta = tnew; F = Fn; v = vn; fnorm = fnew
                lam = max(lam / 3, mp.mpf(2) ** (-int(sys_mp.prec)))
                accepted = True
                break
            lam *= 10
        if not accepted:
            break
    v = _eigenmode_at(sys_mp, theta, ref, th0_ref)
    sys_mp.set_theta(zangle(v))
    r = eigenmode_residual(sys_mp, v)
    return v, float(mp.log(r + mp.mpf(2) ** (-4 * sys_mp.prec)) / mp.log(10))


def eigenmode_residual(sys_mp, v):
    """μ = Re<v, iKv>, ‖iKv - μv‖。"""
    kv = sys_mp.kmatvec(v)
    ikv = [mp.mpc(0, 1) * x for x in kv]
    mu = mp.re(cdot(v, ikv))
    r = [ikv[i] - mu * v[i] for i in range(len(v))]
    return znorm(r)


# ---------------- 観測（mpmath, 固有平面占有） ----------------

def eigenplanes_mp(sys_mp, z, kmax=8):
    """K(θ(z)) の固有平面 w_j（σ降順, 単位ノルム, len M 複素）上位 kmax。観測のみ。"""
    sys_mp.set_theta(zangle(z))
    E, EV = sys_mp.jg_eig()
    ims = [(mp.im(E[k]), k) for k in range(len(E))]
    ims.sort(key=lambda x: x[0], reverse=True)
    sig = []; Ws = []
    n2 = 2 * sys_mp.n
    for im, k in ims:
        if im <= mp.mpf(2) ** (-sys_mp.prec + 8):
            continue
        col = [EV[i, k] for i in range(n2)]
        w = sys_mp.w(col)
        nw = znorm(w)
        if nw == 0:
            continue
        w = [wi / nw for wi in w]
        sig.append(im); Ws.append(w)
        if len(Ws) >= kmax:
            break
    return sig, Ws


def plane_energies_mp(Ws, z):
    """E_j = |<w_j,z>|^2/|z|^2（float64 で返す）。"""
    tot = sum((mp.re(v) ** 2 + mp.im(v) ** 2) for v in z)
    out = []
    for w in Ws:
        pr = cdot(w, z)
        e = (mp.re(pr) ** 2 + mp.im(pr) ** 2) / tot
        out.append(f64(e))
    return np.array(out)
