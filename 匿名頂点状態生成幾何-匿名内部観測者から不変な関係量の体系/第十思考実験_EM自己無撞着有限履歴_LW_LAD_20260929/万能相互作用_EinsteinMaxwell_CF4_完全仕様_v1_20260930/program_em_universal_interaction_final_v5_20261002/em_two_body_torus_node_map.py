#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
em_two_body_torus_node_map.py — 仕様 最終版 v7（2026-10-03）
万能相互作用 Einstein–Maxwell（水素型二体）：2 次元トーラス（呼吸位相 θ_X、相対角 φ）上の一つの波 Ψ_{n,m}

仕様: ../coefficients_relational_full_list_ja_final_v7_20261003.md

表現（統一・匿名）
  状態 = Ψ_{n,m}（複素倍音係数）。n = 0…N は呼吸の倍音（動径作用 I_r = nħ、§10.3「呼吸は止まってよい」で n ≥ 0）、
  m = 1…M は回転の倍音（L = mħ、公理②「停止した電子は存在しない」で m ≥ 1）。Σ|Ψ|² = 1。ξ, π, B は状態にない。
  ほかに θ_a, θ_b（内部位相）、s_a, s_b（軌道軸に沿ったスピン、定数）、f（吸収率）。

初期化で E–M から作る表（相互作用中は固定）
  E(n,m)      : 作用 (I_r, L) = (2g n, 2g m)（μ で割った値）の軌道のエネルギー。ε(ξ,π,m) の 22 項を作用–角変数で一周求積して反転。
  Z_k(n̄,m̄)   : 軌道の Fourier 係数 ξ e^{i(φ − γθ_X)} = Σ_k Z_k e^{ikθ_X}。双極子 r e^{iφ} の行列要素
                ⟨n+k, m+1| r e^{iφ} |n, m⟩ = Z_k(n + k/2, m + ½)（Bohr の対応原理、Heisenberg 1925。中点の作用で評価）。
  W_k         : 同じく ξ² e^{2i(φ−γθ_X)}（四重極、Δm = 2）。  X_k : ξ（Δm = 0、呼吸の読み出し）。
  ⟨𝐏²⟩, ⟨1/ξ⟩, ⟨ξ⟩ : 軌道平均（時計因子に使う）。

節点写像（相互作用。割り算なし、積和と位相の乗算だけ）
  0. W = |Ψ|²、ω_X = Σ W Ω_X(n,m)、ρ̄ = Σ W ρ(n,m)、⟨v²⟩、⟨1/ξ⟩ → F_a、F_b
     Ω_X(n,m) = [E(n+1,m) − E(n,m)]/(2g)：呼吸の拍動（隣り合う n の Bohr 振動数）、ρ(n,m) = [E(n,m+1) − E(n,m)]/(2g)：回転の拍動。
     二つのチャネルで同じ定義。γ = ρ̄/ω_X。
  1. Δθ_k = Δθ_X · ω_X⁻¹（逆数は Newton 反復）
  2. D = ⟨Ψ|R|Ψ⟩（双極子）、Q = ⟨Ψ|Q̂|Ψ⟩（四重極）。直前 4・6 節点の後退差分で D‴、Q⁽⁵⁾。
  3. (VΨ) = −f (κ_C/3)[conj(D‴) RΨ + D‴ R†Ψ] + f (αμm_a/90)[Q⁽⁵⁾ Q̂†Ψ + conj(Q⁽⁵⁾) Q̂Ψ]
  4. Ψ ← [Ψ − i (1/2g) Δθ_k (VΨ)] · exp(−i E(n,m) Δθ_k /(2g))      （1/(2g) = μ/(2m_a) = a_a/ħ）
  5. θ_X += Δθ_X、θ_Y += ρ̄Δθ_k、θ_a += F_aΔθ_k、θ_b += (m_b/m_a)F_bΔθ_k
  6. 読み出し（発展に戻さない）
"""
import argparse
import cmath
import csv
import hashlib
import json
import math
import os
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# =============================================================================
# §0.4 結合定数、§1.2 単項式表（v6 と同じ）
# =============================================================================

class Const:
    pass


def build_constants(alpha=7.2973525643e-3, ma=7.0e-4, mb_over_ma=1836.15267343,
                    M=31, N=8, K=8, nodes_per_breath=62, nterms_exp=14, newton_only=False):
    C = Const()
    C.alpha, C.ma = alpha, ma
    C.mb = ma * mb_over_ma
    C.mb_over_ma = mb_over_ma
    C.Mtot = C.ma + C.mb
    C.mu = C.ma * C.mb / C.Mtot
    C.nu = C.mu / C.Mtot
    C.g = C.ma / C.mu
    C.Ja = 9.0 / (2.0 * alpha)
    C.aa = C.Ja / C.ma
    C.kC = 9.0 / (C.mu * C.aa)
    C.kG = C.Mtot / C.aa
    C.kappa = C.kC + C.kG
    C.C0m1 = 9.0 / (C.ma * C.mb)
    C.M, C.N, C.K = M, N, K
    C.dthX = 2.0 * math.pi / nodes_per_breath
    C.dthX_inv = 1.0 / C.dthX
    C.inv2g = 1.0 / (2.0 * C.g)                 # μ/(2m_a) = a_a/ħ
    C.g_inv2 = 1.0 / C.g ** 2
    C.r_over_g2 = ((C.g - 1.0) / C.g) ** 2
    C.Ua_coef = C.kG / C.g
    C.Ub_coef = C.kG * (C.g - 1.0) / C.g
    C.rr_em = C.kC / 3.0
    C.rr_gw = alpha * C.mu * C.ma / 90.0
    C.eloss_em = 2.0 * C.kC / 3.0
    C.eloss_gw = alpha * C.mu * C.ma / 45.0
    C.nterms_exp = nterms_exp
    C.inv_k = np.array([0.0] + [1.0 / k for k in range(1, nterms_exp + 1)])
    C.exp_nsq = 4
    C.exp_prescale = 1.0 / 2 ** C.exp_nsq
    C.newton_only = newton_only
    kC, kG, nu, g = C.kC, C.kG, C.nu, C.g
    rows = [
        (0, 2, 0, 0, 0, 0.5),
        (0, 4, 0, 0, 0, (3 * nu - 1) / 8),
        (0, 6, 0, 0, 0, (1 - 5 * nu + 5 * nu ** 2) / 16),
        (1, 0, 0, 0, 0, -(kC + kG)),
        (1, 2, 0, 0, 0, -nu * (kC + kG) - 1.5 * kG),
        (1, 4, 0, 0, 0, kC * nu * (1 - 2 * nu) / 2),
        (2, 0, 0, 0, 0, kG * (3 * kC + kG) / 2),
        (2, 0, 2, 0, 0, 2 * g ** 2),
        (2, 2, 0, 0, 0, kC ** 2 * nu / 4),
        (2, 2, 2, 0, 0, g ** 2 * (3 * nu - 1)),
        (2, 4, 2, 0, 0, 3 * g ** 2 * (1 - 5 * nu + 5 * nu ** 2) / 4),
        (3, 0, 0, 0, 0, -(kC + kG) * (g ** 2 - 2 * g + 2) / 2 + kC ** 3 * nu / 4),
        (3, 0, 0, 1, 1, -(g - 1) * (kC + kG)),
        (3, 0, 1, 0, 1, kC * (g ** 2 - 1) + kG * (3 * g ** 2 - 2 * g - 1)),
        (3, 0, 1, 1, 0, kC * (2 * g - 1) + kG * (4 * g - 1)),
        (3, 0, 2, 0, 0, -2 * g ** 2 * (nu * (kC + kG) + 3 * kG)),
        (3, 2, 2, 0, 0, g ** 2 * kC * nu * (3 - 4 * nu)),
        (4, 0, 2, 0, 0, g ** 2 * kC ** 2 * nu),
        (4, 0, 4, 0, 0, g ** 4 * (6 * nu - 2)),
        (4, 2, 4, 0, 0, g ** 4 * (3 - 15 * nu + 15 * nu ** 2)),
        (5, 0, 4, 0, 0, 2 * g ** 4 * kC * nu * (2 - 3 * nu)),
        (6, 0, 6, 0, 0, g ** 6 * (4 - 20 * nu + 20 * nu ** 2)),
    ]
    if newton_only:                      # 検証用：Newton＋クーロンだけ（Kepler の閉じた式と比較する）
        rows = [rows[0], rows[3], rows[7]]
        rows[1] = (1, 0, 0, 0, 0, -C.kappa)
    C.P = np.array([r[0] for r in rows], dtype=float)
    C.A = np.array([r[1] for r in rows], dtype=float)
    C.Bexp = np.array([r[2] for r in rows], dtype=float)
    C.C1 = np.array([r[3] for r in rows], dtype=float)
    C.C2 = np.array([r[4] for r in rows], dtype=float)
    C.coef = np.array([r[5] for r in rows], dtype=float)
    return C


def spin_factor(C, sa_z, sb_z):
    return np.where(C.C1 > 0, sa_z, 1.0) * np.where(C.C2 > 0, sb_z, 1.0)


def recip(x, y_seed, n=3):
    """1/x を Newton 反復 y ← y(2 − xy) で（積和のみ）。種は前節点の値。"""
    y = y_seed
    for _ in range(n):
        y = y * (2.0 - x * y)
    return y


def cexp_series(z, C):
    """exp(z)：z/2⁴ の級数（Horner、14 項）→ 4 回二乗 → 絶対値 1 に補正。積和のみ。"""
    w = z * C.exp_prescale
    acc = np.ones_like(w)
    for k in range(C.nterms_exp, 0, -1):
        acc = 1.0 + w * C.inv_k[k] * acc
    for _ in range(C.exp_nsq):
        acc = acc * acc
    return acc * (1.5 - 0.5 * (acc.real * acc.real + acc.imag * acc.imag))


# =============================================================================
# 初期化：E–M の表 ε(ξ, π, m) から作用–角変数の表を作る（ここでは割り算・求積を使ってよい）
# =============================================================================

class EpsM:
    """倍音番号 m（実数でよい）とスピンを固定した ε(ξ, π) とその微分。"""

    def __init__(self, C, m, sa, sb):
        sf = spin_factor(C, sa, sb)
        self.P, self.A = C.P, C.A
        self.cm = C.coef * sf * m ** C.Bexp
        self.cdm = C.coef * sf * C.Bexp * np.where(C.Bexp > 0, m ** np.maximum(C.Bexp - 1, 0), 0.0)
        self.Apos = C.A > 0
        self.Am1 = np.maximum(C.A - 1, 0)
        self.inv2g = C.inv2g
        self.m = m
        self.g = C.g

    def eps(self, xi, pi_):
        return float(np.sum(self.cm * pi_ ** self.A * xi ** (-self.P)))

    def d_xi(self, xi, pi_):
        return float(np.sum(-self.P * self.cm * pi_ ** self.A * xi ** (-self.P - 1)))

    def d_pi(self, xi, pi_):
        return float(np.sum(np.where(self.Apos, self.A * self.cm * pi_ ** self.Am1, 0.0) * xi ** (-self.P)))

    def d_m(self, xi, pi_):
        return float(np.sum(self.cdm * pi_ ** self.A * xi ** (-self.P)))

    def rhs(self, t, y):
        xi, pi_, phi, I = y
        xp = xi ** (-self.P)
        pa = pi_ ** self.A
        dxi = float(np.sum(np.where(self.Apos, self.A * self.cm * pi_ ** self.Am1, 0.0) * xp))
        dpi = -float(np.sum(-self.P * self.cm * pa * xp / xi))
        dphi = float(np.sum(self.cdm * pa * xp)) * self.inv2g
        return [dxi, dpi, dphi, pi_ * dxi]

    def circular(self):
        """∂ε/∂ξ = 0 の根 ξ_c、E_c、ω_c（小振幅の呼吸）、ρ_c（回転）。"""
        xi = (2.0 * self.g * self.m) ** 2 / (-self.cm[np.where((self.P == 1) & (self.A == 0))[0][0]])
        for _ in range(60):
            f = self.d_xi(xi, 0.0)
            a0 = (self.A == 0)
            df = float(np.sum(self.P[a0] * (self.P[a0] + 1) * self.cm[a0] * xi ** (-self.P[a0] - 2)))
            step = f / df
            xi -= step
            if abs(step) < 1e-15 * xi:
                break
        a0 = (self.A == 0)
        exx = float(np.sum(self.P[a0] * (self.P[a0] + 1) * self.cm[a0] * xi ** (-self.P[a0] - 2)))
        a2 = (self.A == 2)
        epp = float(np.sum(2.0 * self.cm[a2] * xi ** (-self.P[a2])))
        return xi, self.eps(xi, 0.0), math.sqrt(exx * epp), self.d_m(xi, 0.0) * self.inv2g


def orbit_at_energy(em, E, xi_c, K, Ns=2048, rtol=1e-12):
    """エネルギー E の軌道を近点→遠点まで求積し、一周（鏡映）の作用・振動数・Fourier 係数・平均を返す。"""
    f = lambda x: em.eps(x, 0.0) - E
    xi_p = brentq(f, xi_c * 0.1, xi_c, xtol=1e-14 * xi_c, rtol=1e-14)      # 近点は ξ_c/(1+e) ≥ 0.5 ξ_c

    def ev_apo(t, y):
        return y[1]
    ev_apo.direction = -1
    ev_apo.terminal = True
    kappa_est = -em.d_xi(xi_c, 0.0) * 0.0 + (-em.eps(xi_c, 0.0)) * 2.0 * xi_c      # Newton 次：E_c = −κ/(2ξ_c)
    T_est = 2.0 * math.pi * kappa_est / (-2.0 * E) ** 1.5                           # Kepler：ω = (−2E)^{3/2}/κ
    sol = solve_ivp(em.rhs, (0.0, 1.5 * T_est), [xi_p, 0.0, 0.0, 0.0], method="DOP853",
                    rtol=rtol, atol=1e-18, events=ev_apo, dense_output=True, max_step=T_est / 200)
    if not sol.t_events[0].size:
        raise RuntimeError("apocenter not found (E=%g, m=%g)" % (E, em.m))
    T_half = float(sol.t_events[0][0])
    y_half = sol.sol(T_half)
    omX = math.pi / T_half
    I_r = float(y_half[3]) / math.pi
    phi_half = float(y_half[2])
    gamma = phi_half / math.pi                       # 一周の回転 / 2π
    th = np.arange(Ns) * 2.0 * math.pi / Ns
    th_half = np.where(th <= math.pi, th, 2.0 * math.pi - th)
    yy = sol.sol(th_half / omX)
    xi_s = yy[0]
    pi_s = np.where(th <= math.pi, yy[1], -yy[1])
    phi_s = np.where(th <= math.pi, yy[2], 2.0 * phi_half - yy[2])
    rel = phi_s - gamma * th                          # 周期的な部分（近点で 0）
    fz = xi_s * np.exp(1j * rel)
    fw = xi_s ** 2 * np.exp(2j * rel)
    Zf = np.fft.fft(fz) / Ns
    Wf = np.fft.fft(fw) / Ns
    Xf = np.fft.fft(xi_s) / Ns
    ks = np.arange(-K, K + 1)
    Z = np.array([Zf[k % Ns] for k in ks])
    W = np.array([Wf[k % Ns] for k in ks])
    X = np.array([Xf[k % Ns] for k in ks])
    P2 = float(np.mean(pi_s ** 2 + (2.0 * em.g * em.m / xi_s) ** 2))
    return dict(E=E, I_r=I_r, omX=omX, gamma=gamma, xi_p=xi_p, xi_a=float(xi_s[Ns // 2]),
                Z=Z, W=W, X=X, P2=P2, invxi=float(np.mean(1.0 / xi_s)), xim=float(np.mean(xi_s)))


def orbit_circular(em, K):
    xi_c, E_c, om_c, rho_c = em.circular()
    Z = np.zeros(2 * K + 1, dtype=complex)
    W = np.zeros(2 * K + 1, dtype=complex)
    X = np.zeros(2 * K + 1, dtype=complex)
    Z[K] = xi_c
    W[K] = xi_c ** 2
    X[K] = xi_c
    return dict(E=E_c, I_r=0.0, omX=om_c, gamma=rho_c / om_c, xi_p=xi_c, xi_a=xi_c, Z=Z, W=W, X=X,
                P2=(2.0 * em.g * em.m / xi_c) ** 2, invxi=1.0 / xi_c, xim=xi_c)


def orbit_at_action(em, n_h, K, xi_c, E_c, rtol):
    """動径作用 I_r = 2g n_h の軌道。E を ∂E/∂I_r = ω_X の Newton 法で求める（Kepler の式を種に）。"""
    if n_h <= 0.0:
        return orbit_circular(em, K)
    target = 2.0 * em.g * n_h
    kappa_eff = -E_c * 2.0 * (2.0 * em.g * em.m)          # Newton: E_c = −κ²/(8g²m²) → κ = sqrt(...)；種には比で十分
    E = E_c * (em.m / (em.m + n_h)) ** 2
    orb = None
    for _ in range(12):
        E = min(E, -1e-12 * abs(E_c))
        E = max(E, E_c * (1.0 - 1e-12))
        orb = orbit_at_energy(em, E, xi_c, K, rtol=rtol)
        dI = target - orb["I_r"]
        if abs(dI) < 1e-10 * 2.0 * em.g:
            break
        E = E + orb["omX"] * dI
    return orb


def build_tables(C, sa_z, sb_z, rtol=1e-12, cache_dir="tables_cache", verbose=True):
    """整数格子 (n = 0..N, m = 1..M) の E、Ω_X、ρ、平均と、半整数格子の Z、W、X。キャッシュあり。"""
    key = hashlib.sha256(json.dumps(dict(alpha=C.alpha, ma=C.ma, r=C.mb_over_ma, M=C.M, N=C.N, K=C.K, sa=sa_z, sb=sb_z,
                                         rtol=rtol, newton=C.newton_only, coef=[float(c) for c in C.coef]), sort_keys=True).encode()).hexdigest()[:16]
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"tables_{key}.npz")
    if os.path.exists(path):
        d = np.load(path)
        return {k: d[k] for k in d.files}
    t0 = time.time()
    N, M, K = C.N, C.M, C.K
    nh = np.arange(0, 2 * N + 1) / 2.0                  # 0, 0.5, …, N
    mh = 1.0 + np.arange(0, 2 * M - 1) / 2.0            # 1, 1.5, …, M
    E_h = np.zeros((len(nh), len(mh)))
    om_h = np.zeros_like(E_h)
    gam_h = np.zeros_like(E_h)
    P2_h = np.zeros_like(E_h)
    invxi_h = np.zeros_like(E_h)
    xim_h = np.zeros_like(E_h)
    xip_h = np.zeros_like(E_h)
    xia_h = np.zeros_like(E_h)
    Ir_h = np.zeros_like(E_h)
    Z_h = np.zeros((2 * K + 1, len(nh), len(mh)), dtype=complex)
    W_h = np.zeros_like(Z_h)
    X_h = np.zeros_like(Z_h)
    for j, m in enumerate(mh):
        em = EpsM(C, m, sa_z, sb_z)
        xi_c, E_c, _, _ = em.circular()
        for i, n in enumerate(nh):
            o = orbit_at_action(em, n, K, xi_c, E_c, rtol)
            E_h[i, j], om_h[i, j], gam_h[i, j] = o["E"], o["omX"], o["gamma"]
            P2_h[i, j], invxi_h[i, j], xim_h[i, j] = o["P2"], o["invxi"], o["xim"]
            xip_h[i, j], xia_h[i, j], Ir_h[i, j] = o["xi_p"], o["xi_a"], o["I_r"]
            Z_h[:, i, j], W_h[:, i, j], X_h[:, i, j] = o["Z"], o["W"], o["X"]
        if verbose and (j % 10 == 0 or j == len(mh) - 1):
            print(f"  tables: m = {m:5.1f}  ({j + 1}/{len(mh)})  {time.time() - t0:6.1f} s", flush=True)
    T = dict(nh=nh, mh=mh, E_h=E_h, om_h=om_h, gam_h=gam_h, P2_h=P2_h, invxi_h=invxi_h, xim_h=xim_h,
             xip_h=xip_h, xia_h=xia_h, Ir_h=Ir_h, Z_h=Z_h, W_h=W_h, X_h=X_h, build_seconds=np.array(time.time() - t0))
    np.savez(path, **T)
    return T


class Tables:
    """整数格子の表と、演算子（R, R†, Q̂, Q̂†, X̂）の係数スライス。"""

    def __init__(self, C, T):
        N, M, K = C.N, C.M, C.K
        self.N, self.M, self.K = N, M, K
        ni = 2 * np.arange(N + 1)                     # 整数 n の半格子添字
        mi = 2 * (np.arange(1, M + 1) - 1)            # 整数 m の半格子添字
        self.E = T["E_h"][np.ix_(ni, mi)]             # [N+1, M]
        self.P2 = T["P2_h"][np.ix_(ni, mi)]
        self.invxi = T["invxi_h"][np.ix_(ni, mi)]
        self.xim = T["xim_h"][np.ix_(ni, mi)]
        self.gam_cl = T["gam_h"][np.ix_(ni, mi)]       # 古典の γ（読み出し比較用）
        self.om_cl = T["om_h"][np.ix_(ni, mi)]
        # 拍動（Bohr 振動数）：両チャネルで同じ定義。端は後退差分。
        E = self.E
        OmX = np.empty_like(E)
        OmX[:-1, :] = (E[1:, :] - E[:-1, :]) * C.inv2g
        OmX[-1, :] = OmX[-2, :]
        Rho = np.empty_like(E)
        Rho[:, :-1] = (E[:, 1:] - E[:, :-1]) * C.inv2g
        Rho[:, -1] = Rho[:, -2]
        self.OmX, self.Rho = OmX, Rho
        self.Ephase = E * C.inv2g                     # 位相の進み率
        Zh, Wh, Xh = T["Z_h"], T["W_h"], T["X_h"]
        # 演算子の係数スライス（k ごと）。(RΨ)[n,m] = Σ_k Z_k(n−k/2, m−½) Ψ[n−k, m−1]
        self.R_terms, self.Rd_terms, self.Q_terms, self.Qd_terms, self.X_terms = [], [], [], [], []
        for kk, k in enumerate(range(-K, K + 1)):
            # R（m を 1 上げる）：出力 (n, m)、入力 (n−k, m−1)。n ∈ [max(k,0), N+min(k,0)]、m ∈ [2, M]
            n_lo, n_hi = max(k, 0), N + min(k, 0)
            if n_hi >= n_lo:
                ns = np.arange(n_lo, n_hi + 1)
                ms = np.arange(2, M + 1)
                coef = Zh[kk][np.ix_(2 * ns - k, 2 * ms - 3)]       # nh = n − k/2 → 2n − k；mh = m − ½ → 2(m−½−1) = 2m − 3
                self.R_terms.append((k, n_lo, n_hi, coef))
                # R†（m を 1 下げる）：出力 (n, m)、入力 (n+k, m+1)、係数 Z_k(n + k/2, m + ½)。n ∈ [max(−k,0), N−max(k,0)]、m ∈ [1, M−1]
                n_lo2, n_hi2 = max(-k, 0), N - max(k, 0)
                ns2 = np.arange(n_lo2, n_hi2 + 1)
                ms2 = np.arange(1, M)
                coef2 = Zh[kk][np.ix_(2 * ns2 + k, 2 * ms2 - 1)]     # mh = m + ½ → 2(m+½−1) = 2m − 1
                self.Rd_terms.append((k, n_lo2, n_hi2, coef2))
                # Q̂（m を 2 上げる）：入力 (n−k, m−2)、係数 W_k(n − k/2, m − 1)。m ∈ [3, M]
                if M >= 3:
                    ms3 = np.arange(3, M + 1)
                    coef3 = Wh[kk][np.ix_(2 * ns - k, 2 * ms3 - 4)]  # mh = m − 1 → 2(m−2) = 2m − 4
                    self.Q_terms.append((k, n_lo, n_hi, coef3))
                    ms4 = np.arange(1, M - 1)
                    coef4 = Wh[kk][np.ix_(2 * ns2 + k, 2 * ms4)]     # mh = m + 1 → 2m
                    self.Qd_terms.append((k, n_lo2, n_hi2, coef4))
                # X̂（Δm = 0）：入力 (n−k, m)、係数 X_k(n − k/2, m)。
                ms5 = np.arange(1, M + 1)
                coef5 = Xh[kk][np.ix_(2 * ns - k, 2 * ms5 - 2)]
                self.X_terms.append((k, n_lo, n_hi, coef5))

    def apply_R(self, Psi):
        out = np.zeros_like(Psi)
        for k, n_lo, n_hi, coef in self.R_terms:
            out[n_lo:n_hi + 1, 1:] += coef * Psi[n_lo - k:n_hi - k + 1, :-1]
        return out

    def apply_Rd(self, Psi):
        out = np.zeros_like(Psi)
        for k, n_lo, n_hi, coef in self.Rd_terms:
            out[n_lo:n_hi + 1, :-1] += coef * Psi[n_lo + k:n_hi + k + 1, 1:]
        return out

    def apply_Q(self, Psi):
        out = np.zeros_like(Psi)
        for k, n_lo, n_hi, coef in self.Q_terms:
            out[n_lo:n_hi + 1, 2:] += coef * Psi[n_lo - k:n_hi - k + 1, :-2]
        return out

    def apply_Qd(self, Psi):
        out = np.zeros_like(Psi)
        for k, n_lo, n_hi, coef in self.Qd_terms:
            out[n_lo:n_hi + 1, :-2] += coef * Psi[n_lo + k:n_hi + k + 1, 2:]
        return out

    def apply_X(self, Psi):
        out = np.zeros_like(Psi)
        for k, n_lo, n_hi, coef in self.X_terms:
            out[n_lo:n_hi + 1, :] += coef * Psi[n_lo - k:n_hi - k + 1, :]
        return out


# =============================================================================
# 状態と初期化（§8）
# =============================================================================

class State:
    pass


def init_state(C, T, delta=0.065, f=1.0, sa_z=1.0, sb_z=-1.0):
    S = State()
    S.f = float(f)
    S.sa_z, S.sb_z = float(sa_z), float(sb_z)
    tb = Tables(C, T)
    S.tb = tb
    N, M = C.N, C.M
    # 初期条件：ξ₁ で接線速度 v₀ = α(1+δ) の対（STATUS 8.5：r、v/c ≈ α で置く）
    #   L/ħ = ⟨m⟩ = 1 + δ、エネルギー ε̄ = ε(ξ₁, 0, m = 1+δ)（22 項すべて）。
    #   最小の台 {(0,1), (0,2), (1,1)}：⟨m⟩ から w(0,2) = δ、ε̄ から w(1,1)。
    em1 = EpsM(C, 1.0, sa_z, sb_z)
    xi1, E01_c, _, _ = em1.circular()
    S.xi1 = xi1
    em_d = EpsM(C, 1.0 + delta, sa_z, sb_z)
    e_target = em_d.eps(xi1, 0.0)
    E01, E02, E11 = tb.E[0, 0], tb.E[0, 1], tb.E[1, 0]
    w02 = delta
    w11 = (e_target - E01 - w02 * (E02 - E01)) / (E11 - E01) if delta > 0 else 0.0
    w01 = 1.0 - w02 - w11
    if w11 < 0 or w01 < 0:
        raise ValueError(f"initial weights invalid: w01={w01} w02={w02} w11={w11}")
    Psi = np.zeros((N + 1, M), dtype=complex)
    Psi[0, 0] = math.sqrt(w01)
    Psi[0, 1] = math.sqrt(w02)
    Psi[1, 0] = math.sqrt(w11)
    S.Psi = Psi
    S.w0 = (w01, w02, w11)
    S.e_target = e_target
    S.delta = delta
    S.theta_X = S.theta_Y = S.theta_a = S.theta_b = 0.0
    S.cum_phase01 = 0.0
    S.k = 0
    S.histD, S.histQ = [], []
    W = np.abs(Psi) ** 2
    S.omX_inv = 1.0 / float(np.sum(W * tb.OmX))
    S.last = {}
    return S


# =============================================================================
# 節点写像（相互作用の唯一の関数。割り算なし）
# =============================================================================

def node_map(C, S):
    tb = S.tb
    Psi = S.Psi
    f = S.f
    W = Psi.real * Psi.real + Psi.imag * Psi.imag
    norm = float(np.sum(W))
    # 手順 0：拍動の平均と時計因子
    omX = float(np.sum(W * tb.OmX))
    rho = float(np.sum(W * tb.Rho))
    P2 = float(np.sum(W * tb.P2))
    invxi = float(np.sum(W * tb.invxi))
    va2 = P2 * C.g_inv2
    vb2 = P2 * C.r_over_g2
    Fa = 1.0 - 0.5 * va2 - C.Ua_coef * invxi - 0.125 * va2 * va2
    Fb = 1.0 - 0.5 * vb2 - C.Ub_coef * invxi
    # 手順 1
    omX_inv = recip(omX, S.omX_inv)
    dth = C.dthX * omX_inv
    dth_inv = omX * C.dthX_inv
    # 手順 2：二次の関係（双極子・四重極・半径）
    RPsi = tb.apply_R(Psi)
    D = complex(np.sum(np.conj(Psi) * RPsi))
    QPsi = tb.apply_Q(Psi)
    Q = complex(np.sum(np.conj(Psi) * QPsi))
    xi_mean = float(np.sum(np.conj(Psi) * tb.apply_X(Psi)).real)
    S.histD.append(D)
    S.histQ.append(Q)
    if len(S.histD) > 4:
        S.histD.pop(0)
    if len(S.histQ) > 6:
        S.histQ.pop(0)
    D3 = D2 = 0j
    Q5 = 0j
    if len(S.histD) == 4:
        h = S.histD
        D3 = (h[3] - 3.0 * h[2] + 3.0 * h[1] - h[0]) * dth_inv ** 3
        D2 = (h[3] - 2.0 * h[2] + h[1]) * dth_inv ** 2
    if len(S.histQ) == 6:
        q = S.histQ
        Q5 = (q[5] - 5.0 * q[4] + 10.0 * q[3] - 10.0 * q[2] + 5.0 * q[1] - q[0]) * dth_inv ** 5
    # 手順 3：放射反作用（両チャネルに同じ形で入る。R は Δm = +1 と Δn = k）
    VPsi = -f * C.rr_em * (np.conj(D3) * RPsi + D3 * tb.apply_Rd(Psi))
    if f != 0.0 and Q5 != 0j:
        VPsi = VPsi + f * C.rr_gw * (Q5 * tb.apply_Qd(Psi) + np.conj(Q5) * QPsi)
    # 手順 4：位相の乗算
    phase = cexp_series(-1j * tb.Ephase * dth, C)
    Psi_new = (Psi - 1j * C.inv2g * dth * VPsi) * phase
    # 手順 5：四つの位相
    dthY = rho * dth
    dtha = Fa * dth
    dthb = C.mb_over_ma * Fb * dth
    S.theta_X += C.dthX
    S.theta_Y += dthY
    S.theta_a += dtha
    S.theta_b += dthb
    S.cum_phase01 += -tb.Ephase[0, 0] * dth
    # 読み出しの記録
    ebar = float(np.sum(W * tb.E))
    S.last = dict(k=S.k, theta_X=S.theta_X - C.dthX, dth=dth, dthY=dthY, dtha=dtha, dthb=dthb,
                  omX=omX, rho=rho, gamma=rho * omX_inv, Fa=Fa, Fb=Fb,
                  v_over_c=float(np.sum(W * (2.0 * C.g * np.arange(1, C.M + 1)[None, :]) * tb.invxi)),
                  D_abs=abs(D), D_arg=cmath.phase(D), Q_abs=abs(Q), xi_mean=xi_mean,
                  ebar=ebar, norm=norm, n_mean=float(np.sum(W * np.arange(C.N + 1)[:, None])),
                  m_mean=float(np.sum(W * np.arange(1, C.M + 1)[None, :])),
                  w01=float(W[0, 0]), w02=float(W[0, 1]), w11=float(W[1, 0]),
                  w03=float(W[0, 2]) if C.M >= 3 else 0.0, w12=float(W[1, 1]), w21=float(W[2, 0]) if C.N >= 2 else 0.0,
                  w_rest=float(norm - W[0, 0] - W[0, 1] - W[1, 0]),
                  eloss_em=-f * C.eloss_em * abs(D2) ** 2, arg_psi01=cmath.phase(Psi[0, 0]))
    S.Psi = Psi_new
    S.omX_inv = omX_inv
    S.k += 1
    return S


READOUT_FIELDS = ["k", "theta_X", "theta_Y_over_2pi", "theta_a", "theta_b", "dth", "dthY", "dtha", "dthb",
                  "omX", "rho", "gamma", "Fa", "Fb", "v_over_c", "D_abs", "D_arg", "Q_abs", "xi_mean",
                  "ebar", "norm", "n_mean", "m_mean", "w01", "w02", "w11", "w03", "w12", "w21", "w_rest",
                  "eloss_em", "arg_psi01", "cum_phase01"]


def readout(C, S):
    row = dict(S.last)
    row["theta_Y_over_2pi"] = S.theta_Y / (2.0 * math.pi)
    row["theta_a"] = S.theta_a
    row["theta_b"] = S.theta_b
    row["cum_phase01"] = S.cum_phase01
    return {k: row[k] for k in READOUT_FIELDS}


# =============================================================================
# 自己検査：Newton のみの表を Kepler の閉じた式（Bessel）と比べる
# =============================================================================

def kepler_reference(C, n_h, m_h, K):
    from scipy.special import jv, jvp
    g, kappa = C.g, C.kappa
    L = 2.0 * g * m_h
    Ir = 2.0 * g * n_h
    Ltot = Ir + L
    E = -kappa ** 2 / (2.0 * Ltot ** 2)
    a = Ltot ** 2 / kappa
    e = math.sqrt(max(0.0, 1.0 - (L / Ltot) ** 2))
    omX = kappa ** 2 / Ltot ** 3
    # (x+iy)/a の Fourier 係数 c_j（e^{ijℓ}）
    def c(j):
        if j == 0:
            return -1.5 * e
        aj = abs(j)
        val = jvp(aj, aj * e) / aj
        val2 = math.sqrt(1 - e * e) * jv(aj, aj * e) / (aj * e) if e > 0 else (0.5 if aj == 1 else 0.0)
        return (val + val2) if j > 0 else (val - val2)
    Z = np.array([a * c(k + 1) for k in range(-K, K + 1)])        # Z_k = c_{k+1}
    return dict(E=E, omX=omX, e=e, a=a, Z=Z, gamma=1.0, P2=kappa / a, invxi=1.0 / a, xim=a * (1 + e * e / 2))


def selftest(C_full):
    C = build_constants(M=4, N=3, K=C_full.K, newton_only=True)
    print("selftest (Newton only, m=1..4, n=0..3): tables vs Kepler closed forms")
    T = build_tables(C, 1.0, -1.0, verbose=False, cache_dir="tables_cache_selftest")
    worst = dict(E=0.0, omX=0.0, Z=0.0, P2=0.0, invxi=0.0, gamma=0.0)
    for i, nh in enumerate(T["nh"]):
        for j, mh in enumerate(T["mh"]):
            ref = kepler_reference(C, nh, mh, C.K)
            worst["E"] = max(worst["E"], abs(T["E_h"][i, j] / ref["E"] - 1))
            worst["omX"] = max(worst["omX"], abs(T["om_h"][i, j] / ref["omX"] - 1))
            worst["gamma"] = max(worst["gamma"], abs(T["gam_h"][i, j] - 1.0))
            worst["P2"] = max(worst["P2"], abs(T["P2_h"][i, j] / ref["P2"] - 1))
            worst["invxi"] = max(worst["invxi"], abs(T["invxi_h"][i, j] / ref["invxi"] - 1))
            worst["Z"] = max(worst["Z"], float(np.max(np.abs(T["Z_h"][:, i, j] - ref["Z"])) / ref["a"]))
    for k, v in worst.items():
        print(f"  max |{k} deviation| = {v:.2e}")
    z = -1j * np.array([0.01, 0.1, 1.0, 10.0])
    print(f"  exp series max err (|z| ≤ 10) = {np.max(np.abs(cexp_series(z, C) - np.exp(z))):.2e}")


# =============================================================================
# 実行
# =============================================================================

def run(args):
    C = build_constants(M=args.M, N=args.N, K=args.K, nterms_exp=args.nterms_exp)
    if args.selftest:
        selftest(C)
        return
    spin = {"up": 1.0, "down": -1.0, "zero": 0.0}
    sa_z, sb_z = spin[args.sa], spin[args.sb]
    t0 = time.time()
    T = build_tables(C, sa_z, sb_z, rtol=args.rtol)
    t_tab = time.time() - t0
    S = init_state(C, T, delta=args.delta, f=args.f, sa_z=sa_z, sb_z=sb_z)
    tb = S.tb
    nodes = int(round(args.breaths * 62)) if args.nodes is None else args.nodes
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    init_info = dict(delta=args.delta, f=args.f, sa=args.sa, sb=args.sb, nodes=nodes, M=C.M, N=C.N, K=C.K, dthX=C.dthX,
                     xi1=S.xi1, e_target=S.e_target, w0=S.w0,
                     E01=float(tb.E[0, 0]), E02=float(tb.E[0, 1]), E11=float(tb.E[1, 0]), E03=float(tb.E[0, 2]), E12=float(tb.E[1, 1]),
                     E20=float(tb.E[2, 0]) if C.N >= 2 else None,
                     OmX01=float(tb.OmX[0, 0]), Rho01=float(tb.Rho[0, 0]), gamma01_beat=float(tb.Rho[0, 0] / tb.OmX[0, 0]),
                     gamma01_classical=float(tb.gam_cl[0, 0]), omX01_classical=float(tb.om_cl[0, 0]),
                     Z0_01=float(T["Z_h"][C.K, 0, 1].real), Zm1_halfhalf=float(T["Z_h"][C.K - 1, 1, 1].real),
                     table_build_seconds=float(T.get("build_seconds", t_tab)))
    t1 = time.time()
    with open(args.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=READOUT_FIELDS)
        wr.writeheader()
        first = None
        for k in range(nodes):
            node_map(C, S)
            if k % args.log_every == 0 or k == nodes - 1:
                row = readout(C, S)
                wr.writerow(row)
                if first is None:
                    first = row
    last = readout(C, S)
    elapsed = time.time() - t1
    W = np.abs(S.Psi) ** 2
    summary = dict(init=init_info, first=first, last=last,
                   weights_final=dict(by_n=W.sum(axis=1).tolist(), by_m=W.sum(axis=0).tolist(),
                                      top=[(int(i), int(j) + 1, float(W[i, j])) for i, j in
                                           zip(*np.unravel_index(np.argsort(W, axis=None)[::-1][:8], W.shape))]),
                   gamma_mean=S.theta_Y / (nodes * C.dthX),
                   a_cycles_per_rotation=(S.theta_a / max(S.theta_Y, 1e-300)),
                   theta_b_over_theta_a=(S.theta_b / S.theta_a if S.theta_a else None),
                   psi01_phase_per_breath_over_pi=(S.cum_phase01 / math.pi) / (nodes / 62.0),
                   elapsed_s=elapsed, nodes_per_s=nodes / max(elapsed, 1e-9))
    with open(os.path.splitext(args.out)[0] + "_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=float)
    print(f"done: {nodes} nodes in {elapsed:.2f} s ({nodes / max(elapsed, 1e-9):.0f} nodes/s); tables {init_info['table_build_seconds']:.1f} s -> {args.out}")
    print(f"  tables: E01={tb.E[0, 0]:.9e} E02={tb.E[0, 1]:.9e} E11={tb.E[1, 0]:.9e}  Ω_X(0,1)={tb.OmX[0, 0]:.6e} ρ(0,1)={tb.Rho[0, 0]:.6e} "
          f"γ_beat(0,1)−1={tb.Rho[0, 0] / tb.OmX[0, 0] - 1:.4e}  γ_cl(0,1)−1={tb.gam_cl[0, 0] - 1:.4e}")
    print(f"  init : xi1={S.xi1:.6f} w(0,1)={S.w0[0]:.5f} w(0,2)={S.w0[1]:.5f} w(1,1)={S.w0[2]:.5f} ebar={first['ebar']:.9e} |D|={first['D_abs']:.4f}")
    print(f"  last : w(0,1)={last['w01']:.6f} w(0,2)={last['w02']:.6f} w(1,1)={last['w11']:.6f} rest={last['w_rest']:.2e} "
          f"ebar={last['ebar']:.9e} |D|={last['D_abs']:.5f} norm-1={last['norm'] - 1:.1e} ⟨ξ⟩={last['xi_mean']:.4f}")
    print(f"  mean : γ_mean−1={summary['gamma_mean'] - 1:.4e}  θ_a/θ_Y={summary['a_cycles_per_rotation']:.4f}  "
          f"θ_b/θ_a={summary['theta_b_over_theta_a']:.7f}  Ψ01 phase/breath/π={summary['psi01_phase_per_breath_over_pi']:.6f}")


def main():
    ap = argparse.ArgumentParser(description="E–M 二体 トーラス波 Ψ_{n,m} の節点写像（最終版 v7）")
    ap.add_argument("--delta", type=float, default=0.065)
    ap.add_argument("--f", type=float, default=1.0)
    ap.add_argument("--sa", choices=["up", "down", "zero"], default="up")
    ap.add_argument("--sb", choices=["up", "down", "zero"], default="down")
    ap.add_argument("--breaths", type=float, default=10.0, help="走行長（節点数 = 62 × breaths。呼吸の拍動の周期数）")
    ap.add_argument("--nodes", type=int, default=None)
    ap.add_argument("--log-every", type=int, default=1)
    ap.add_argument("--M", type=int, default=31)
    ap.add_argument("--N", type=int, default=8)
    ap.add_argument("--K", type=int, default=8)
    ap.add_argument("--rtol", type=float, default=1e-12)
    ap.add_argument("--nterms-exp", type=int, default=14)
    ap.add_argument("--out", default="results_v7/run.csv")
    ap.add_argument("--selftest", action="store_true")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
