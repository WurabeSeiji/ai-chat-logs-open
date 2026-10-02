#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
em_two_body_node_map.py
万能相互作用 Einstein–Maxwell（水素型二体）— 節点写像の実装（仕様 最終版 v6）

仕様: ../coefficients_relational_full_list_ja_final_v6_20261002.md
      節番号は本ファイル内のコメントで「§」として参照する。

構成（統一関数は二つだけ）
  node_map(C, S)  : §2.1 手順 0〜8 を 1 節点ぶん進める。相互作用の唯一の関数。
  readout(C, S)   : §2.2・§7 の読み出し。発展へ帰還させない。
  init_state(...) : §8 初期化（ここだけ E–M の表を解く。割り算・平方根・反復・常微分方程式の求積を使ってよい）。

相互作用中（node_map）の規約（§6）
  * 割り算を使わない。逆数（1/ξ, 1/ω_X, 1/det J）は前節点の値を種にした Newton 反復 y ← y(2 − xy) で得る。
  * exp(iφ) は固定項数の級数（z/2^8 で級数 → 8 回二乗 → 絶対値を 1 に補正）。e^{iΔθ_X} は初期化で定数。
  * 保存部分は位相の乗算だけ。別の積分器を置かない。
  * 新しい係数は現れない：§1.2 の 22 単項式の指数をずらすだけで ∂ε/∂ξ、∂ε/∂π、∂ε/∂m、二階微分を作る。
  * 呼吸の軌道（§2.1 手順 2）は初期化で ε̄ の動径運動を一周だけ機械精度で解き、呼吸位相 θ_X の Fourier 級数
    ξ = ξ̄_c [Ĉ_0(s) + Σ_n Ĉ_n(s) Re((ζ_c B)^n)]、π = ω_X ξ̄_c Σ_n Ŝ_n(s) Im((ζ_c B)^n)、ω_X = ω_X0 Ω̂(s)、s = |B|² ζ_c²
    として保持する（Ĉ、Ŝ、Ω̂ は s の多項式。相互作用中は積和のみ）。打ち切りは読み出し基準 2.7×10⁻⁸ で決める。

実装上の注記（仕様と字面が異なる箇所。README にも記す）
  (a) 放射反作用の動径成分と、|ψ_m|² の変化による呼吸の中心・周波数の変化は、(ξ_k, π_k + F_r Δθ_k) を目標にして
      B を Newton 1〜2 歩で更新する（接触軌道要素の更新）。2×2 の Jacobian の逆は 1/det を逆数反復で。
  (b) ω_X0 = sqrt(ε̄_ξξ ε̄_ππ) に平方根を使う（割り算ではない）。
  (c) D''', Q^(5) は仕様どおり直前 4・6 節点の後退差分（1 次精度。位相の遅れ 1.5 節点・2.5 節点）。
"""
import argparse
import cmath
import csv
import json
import math
import os
import time

import numpy as np

# =============================================================================
# §0.4 結合定数、§1.2 単項式表（初期化でのみ計算）
# =============================================================================

class Const:
    pass


def build_constants(alpha=7.2973525643e-3, ma=7.0e-4, mb_over_ma=1836.15267343,
                    M=31, nodes_per_breath=62, nterms_exp=14, n_harm=12, poly_deg=5):
    C = Const()
    C.alpha = alpha
    C.ma = ma
    C.mb = ma * mb_over_ma
    C.mb_over_ma = mb_over_ma
    C.Mtot = C.ma + C.mb
    C.mu = C.ma * C.mb / C.Mtot
    C.nu = C.mu / C.Mtot
    C.g = C.ma / C.mu                       # = 1 + ma/mb
    C.Ja = 9.0 / (2.0 * alpha)              # ħ/2 = J_a = J_b
    C.aa = C.Ja / C.ma                      # リング半径（長さの単位）
    C.kC = 9.0 / (C.mu * C.aa)              # = 2 alpha g
    C.kG = C.Mtot / C.aa                    # = kC/(C0-1)
    C.kappa = C.kC + C.kG
    C.C0m1 = 9.0 / (C.ma * C.mb)
    C.M = M
    C.marr = np.arange(-M, M + 1, dtype=float)
    C.nm = 2 * M + 1
    C.dthX = 2.0 * math.pi / nodes_per_breath          # 節点 = 呼吸位相の固定増分（§6）
    C.dthX_inv = 1.0 / C.dthX
    C.eidthX = cmath.exp(1j * C.dthX)                   # 初期化で定数
    C.mu_over_2ma = C.mu / (2.0 * C.ma)                 # a_a/ħ 換算（位相の進み）
    C.kappa_inv = 1.0 / C.kappa
    C.g_inv = 1.0 / C.g
    C.g_inv2 = C.g_inv ** 2
    C.r_over_g2 = ((C.g - 1.0) / C.g) ** 2
    C.Ua_coef = C.kG * C.g_inv                          # U_a = kG/(g ξ)
    C.Ub_coef = C.kG * (C.g - 1.0) * C.g_inv            # U_b = kG (g-1)/(g ξ)
    C.rr_em = C.kC / 3.0                                # (V^EM)_m 係数
    C.rr_gw = alpha * C.mu * C.ma / 90.0                # (V^GW)_m 係数 = μ/(20 a_a)
    C.fr_em = 2.0 * C.kC / 3.0                          # F_r^EM
    C.fr_gw = 2.0 * alpha * C.mu * C.ma / 45.0          # F_r^GW
    C.eloss_em = 2.0 * C.kC / 3.0
    C.eloss_gw = alpha * C.mu * C.ma / 45.0
    # スピン歳差（§2.3）
    C.two_alpha = 2.0 * alpha
    C.two_alpha_over_C0m1 = 2.0 * alpha / C.C0m1
    C.two_alpha_ma2_over9 = 2.0 * alpha * C.ma ** 2 / 9.0
    # exp 級数の定数
    C.nterms_exp = nterms_exp
    C.inv_k = np.array([0.0] + [1.0 / k for k in range(1, nterms_exp + 1)])
    C.exp_nsq = 8                                       # 二乗の回数（|z| ≤ 256 まで）
    C.exp_prescale = 1.0 / 2 ** C.exp_nsq
    # 呼吸の級数
    C.n_harm = n_harm
    C.poly_deg = poly_deg
    C.readout_precision = 1.0 / 124.0 / 3.0e5           # εN < 1/124、N ≤ 3×10⁵ 周
    # 単項式表 §1.2：(p, a, b, c1, c2, coef)
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
    C.P = np.array([r[0] for r in rows], dtype=float)
    C.A = np.array([r[1] for r in rows], dtype=float)
    C.Bexp = np.array([r[2] for r in rows], dtype=float)
    C.C1 = np.array([r[3] for r in rows], dtype=float)
    C.C2 = np.array([r[4] for r in rows], dtype=float)
    C.coef = np.array([r[5] for r in rows], dtype=float)
    C.nmono = len(rows)
    C.is_newton = np.zeros(C.nmono, dtype=bool)
    C.is_newton[3] = True   # (1,0,0): -kappa/xi
    C.is_newton[7] = True   # (2,0,2): 2 g^2 m^2/xi^2
    C.mpow = np.array([C.marr ** b for b in range(7)])          # [7, nm]
    C.mono_m = np.array([C.mpow[int(b)] for b in C.Bexp])        # [nmono, nm]
    C.mono_dm = np.array([b * (C.mpow[int(b) - 1] if b >= 1 else 0.0 * C.marr) for b in C.Bexp])
    return C


# =============================================================================
# 割り算なしの道具（相互作用中に使う）
# =============================================================================

def recip(x, y_seed, n=2):
    """1/x を Newton 反復 y <- y(2 - x y) で得る。種は前節点の値。積和のみ。"""
    y = y_seed
    for _ in range(n):
        y = y * (2.0 - x * y)
    return y


def cexp_series(z, C):
    """exp(z) を固定項数の級数（Horner）＋固定回数の二乗で。z は複素 ndarray。
    倍音 m の 1 節点の位相 E_m Δθ_k は m² に比例し、m = 31 では約 50〜60 rad になるので、
    z を 2^{-n_sq} に縮めて級数を取り、n_sq 回二乗して戻す（積和のみ、|z| ≤ 2^{n_sq} で有効）。"""
    w = z * C.exp_prescale
    acc = np.ones_like(w)
    for k in range(C.nterms_exp, 0, -1):
        acc = 1.0 + w * C.inv_k[k] * acc
    for _ in range(C.exp_nsq):
        acc = acc * acc
    acc = acc * (1.5 - 0.5 * (acc.real * acc.real + acc.imag * acc.imag))   # |acc| → 1（Newton 1 歩）
    return acc


def spin_factor(C, sa_z, sb_z):
    fa = np.where(C.C1 > 0, sa_z, 1.0)
    fb = np.where(C.C2 > 0, sb_z, 1.0)
    return fa * fb


def eps_modes(C, zeta, pi_k, sa_z, sb_z):
    """ε(ξ_k, π_k, m, s) と ∂ε/∂m を倍音配列で返す。ζ = 1/ξ_k。"""
    zp = zeta ** C.P
    pa = pi_k ** C.A
    sf = spin_factor(C, sa_z, sb_z)
    w = C.coef * zp * pa * sf
    eps_m = np.dot(w, C.mono_m)
    deps_dm = np.dot(w, C.mono_dm)
    return eps_m, deps_dm


def moments_of(C, psi, sa_z, sb_z):
    """単項式ごとの ⟨m^b⟩ s 因子 Mom_j と、⟨m^b⟩ (b=0..6)。"""
    wgt = np.abs(psi) ** 2
    mom_b = np.dot(C.mpow, wgt)
    Mom = mom_b[C.Bexp.astype(int)] * spin_factor(C, sa_z, sb_z)
    return Mom, mom_b


def ebar_and_grad(C, Mom, xi, pi_):
    """ε̄(ξ, π) = Σ_j c_j π^{a_j} Mom_j ξ^{-p_j} と ∂ε̄/∂ξ、∂ε̄/∂π（初期化の軌道求積で使う。割り算可）。"""
    z = 1.0 / xi
    base = C.coef * Mom
    e = np.sum(base * pi_ ** C.A * z ** C.P)
    de_dxi = -np.sum(base * C.P * pi_ ** C.A * z ** (C.P + 1.0))
    a_pos = C.A >= 1
    de_dpi = np.sum(base[a_pos] * C.A[a_pos] * pi_ ** (C.A[a_pos] - 1.0) * z ** C.P[a_pos])
    return e, de_dxi, de_dpi


def node_constants(C, psi, sa_z, sb_z, zeta_seed):
    """§2.1 手順 0／5：|ψ_m|² のモーメントから ξ̄_c、ζ_c = 1/ξ̄_c、ω_X0 を作る（π = 0）。
    Newton 次 ζ_0 = κ/(4 g² ⟨m²⟩) に、非 Newton 単項式の 1 次補正を Newton 法 1 歩で加える。"""
    Mom, mom_b = moments_of(C, psi, sa_z, sb_z)
    a0 = (C.A == 0)
    xi0 = 4.0 * C.g ** 2 * mom_b[2] * C.kappa_inv
    zeta0 = recip(xi0, zeta_seed, n=3)
    sel = a0 & (~C.is_newton)
    P1 = -np.sum(C.P[sel] * C.coef[sel] * Mom[sel] * zeta0 ** (C.P[sel] + 1.0))
    zeta_c = zeta0 + P1 * C.kappa_inv * xi0
    xi_c = recip(zeta_c, xi0, n=3)
    eps_xx = np.sum(C.coef[a0] * C.P[a0] * (C.P[a0] + 1.0) * Mom[a0] * zeta_c ** (C.P[a0] + 2.0))
    a2 = (C.A == 2)
    eps_pp = np.sum(2.0 * C.coef[a2] * Mom[a2] * zeta_c ** C.P[a2])
    omX0 = math.sqrt(eps_xx * eps_pp)
    return xi_c, zeta_c, omX0, mom_b


def quad_relations(psi):
    S1 = np.sum(psi[:-1] * np.conj(psi[1:]))
    S2 = np.sum(psi[:-2] * np.conj(psi[2:]))
    return S1, S2


def shift(psi, d):
    out = np.zeros_like(psi)
    if d > 0:
        out[:-d] = psi[d:]
    elif d < 0:
        out[-d:] = psi[:d]
    else:
        out[:] = psi
    return out


def nn_tensor(S2):
    re, im = S2.real, S2.imag
    return np.array([[0.5 * (1 + re), 0.5 * im, 0.0],
                     [0.5 * im, 0.5 * (1 - re), 0.0],
                     [0.0, 0.0, 0.0]])


# =============================================================================
# 呼吸の級数（§2.1 手順 2）。評価は積和のみ。
# =============================================================================

def polyval_and_deriv(coefs, s):
    """coefs[n, d]: n ごとの s の多項式（d = 0..D）。値と s 微分を返す（Horner、積和）。"""
    D = coefs.shape[1] - 1
    val = coefs[:, D].copy()
    der = np.zeros_like(val)
    for d in range(D - 1, -1, -1):
        der = der * s + val
        val = val * s + coefs[:, d]
    return val, der


def breathing_eval(C, T, B, xi_c, zeta_c, omX0):
    """B、ξ̄_c、ζ_c、ω_X0 から ξ、π、ω_X と Jacobian ∂(ξ,π)/∂(Re B, Im B) を返す。"""
    x, y = B.real, B.imag
    s = (x * x + y * y) * zeta_c * zeta_c
    Ch, dCh = polyval_and_deriv(T["Chat"], s)      # n = 0..N
    Sh, dSh = polyval_and_deriv(T["Shat"], s)      # n = 0..N（n=0 は 0）
    Om, dOm = polyval_and_deriv(T["Ohat"], s)      # 1 行
    Om, dOm = Om[0], dOm[0]
    N = C.n_harm
    zB = zeta_c * B
    P = np.empty(N + 1, dtype=complex)
    P[0] = 1.0
    for n in range(1, N + 1):
        P[n] = P[n - 1] * zB
    Pr, Pi = P.real, P.imag
    nvec = np.arange(N + 1, dtype=float)
    sx, sy = 2.0 * x * zeta_c * zeta_c, 2.0 * y * zeta_c * zeta_c
    sum_c = np.dot(Ch, Pr)                                   # Ĉ_0 + Σ Ĉ_n Re P_n（P_0 = 1）
    sum_s = np.dot(Sh, Pi)
    xi = xi_c * sum_c
    omX = omX0 * Om
    pi_ = omX * xi_c * sum_s
    # 微分
    Prm1 = np.concatenate([[0.0], Pr[:-1]])                  # Re P_{n-1}
    Pim1 = np.concatenate([[0.0], Pi[:-1]])                  # Im P_{n-1}
    dsum_c_x = np.dot(dCh, Pr) * sx + np.dot(Ch * nvec * zeta_c, Prm1)
    dsum_c_y = np.dot(dCh, Pr) * sy - np.dot(Ch * nvec * zeta_c, Pim1)
    dsum_s_x = np.dot(dSh, Pi) * sx + np.dot(Sh * nvec * zeta_c, Pim1)
    dsum_s_y = np.dot(dSh, Pi) * sy + np.dot(Sh * nvec * zeta_c, Prm1)
    J11 = xi_c * dsum_c_x
    J12 = xi_c * dsum_c_y
    J21 = omX0 * xi_c * (dOm * sx * sum_s + Om * dsum_s_x)
    J22 = omX0 * xi_c * (dOm * sy * sum_s + Om * dsum_s_y)
    return xi, pi_, omX, (J11, J12, J21, J22), s


def build_breathing_table(C, Mom, xi_c0, zeta_c0, omX0, xi_peri0, e_max=None):
    """初期化：ε̄ の動径運動を一周だけ機械精度で解き、θ_X の Fourier 級数の係数を s の多項式として展開する。
    返り値 T: Chat[n,d], Shat[n,d], Ohat[1,d], 残差、初期の A_0、s_0。"""
    from scipy.integrate import solve_ivp

    N = C.n_harm
    e_K0 = xi_c0 / xi_peri0 - 1.0                       # Kepler 相当の離心率（pericenter と p から）
    circular0 = e_K0 < 1e-6                              # δ = 0：呼吸なし。級数は格子だけで作り、A_0 = 0
    if e_max is None:
        e_max = max(1.4 * e_K0, 0.05)
    e_grid = np.linspace(0.004, e_max, 18)
    if not circular0:
        e_grid = np.concatenate([e_grid, [e_K0]])
    Nsamp = 1024

    def rhs(t, yv):
        xi, pi_ = yv
        _, de_dxi, de_dpi = ebar_and_grad(C, Mom, xi, pi_)
        return [de_dpi, -de_dxi]

    def ev_apo(t, yv):
        return yv[1]
    ev_apo.direction = -1
    ev_apo.terminal = True

    rows = []
    for e in e_grid:
        xi_p = xi_c0 / (1.0 + e)
        T_est = 2.0 * math.pi / omX0
        sol = solve_ivp(rhs, (0.0, 1.5 * T_est), [xi_p, 0.0], method="DOP853", rtol=1e-13, atol=1e-18,
                        events=ev_apo, dense_output=True, max_step=T_est / 200)
        if not sol.t_events[0].size:
            raise RuntimeError("apocenter not found for e=%g" % e)
        T_half = float(sol.t_events[0][0])
        omX = math.pi / T_half
        # 平均近点角 θ_X ∈ [0, 2π) を等間隔に。前半 [0, π] は直接、後半は鏡映 ξ(2π−θ) = ξ(θ), π(2π−θ) = −π(θ)
        thX = np.arange(Nsamp) * 2.0 * math.pi / Nsamp
        th_half = np.where(thX <= math.pi, thX, 2.0 * math.pi - thX)
        yv = sol.sol(th_half / omX)
        xi_s = yv[0]
        pi_s = np.where(thX <= math.pi, yv[1], -yv[1])
        xi_apo = float(sol.sol(T_half)[0])
        A = 0.5 * (xi_apo - xi_p)
        q = A * zeta_c0
        s = q * q
        n = np.arange(N + 1)
        cosm = np.cos(np.outer(n, thX))
        sinm = np.sin(np.outer(n, thX))
        Cn = (2.0 / Nsamp) * np.dot(cosm, xi_s)
        Cn[0] *= 0.5
        Sn = (2.0 / Nsamp) * np.dot(sinm, pi_s)
        Chat = Cn / xi_c0 / np.where(n > 0, q ** n, 1.0)
        Shat = Sn / (omX * xi_c0) / np.where(n > 0, q ** n, 1.0)
        rows.append(dict(e=e, s=s, A=A, omX=omX, Cn=Cn, Sn=Sn, Chat=Chat, Shat=Shat, Ohat=omX / omX0, xi_p=xi_p, xi_apo=xi_apo))
    s_arr = np.array([r["s"] for r in rows])
    D = C.poly_deg

    def fit(key_hat, key_abs, scale, k):
        """Ĉ_n(s) を s の多項式で当てる。|C_n| が読み出し基準より小さい格子点（q^n が小さく Ĉ_n が丸めの雑音になる点）は除く。"""
        vals = np.array([r[key_hat][k] for r in rows])
        mags = np.array([abs(r[key_abs][k]) for r in rows]) / scale
        sel = mags >= 1e-8                         # 読み出し基準より小さい項は丸めの雑音：当てはめに使わない
        if sel.sum() == 0:
            return np.zeros(D + 1)                 # 全格子点で無視できる項は 0 にする（雑音を多項式にしない）
        deg = min(D, int(sel.sum()) - 1)
        c = np.polyfit(s_arr[sel], vals[sel], deg)[::-1]
        return np.concatenate([c, np.zeros(D + 1 - len(c))])
    Chat = np.array([fit("Chat", "Cn", xi_c0, k) for k in range(N + 1)])
    Shat = np.array([fit("Shat", "Sn", omX0 * xi_c0, k) for k in range(N + 1)])
    Ohat = np.array([np.polyfit(s_arr, np.array([r["Ohat"] for r in rows]), D)[::-1]])
    T = dict(Chat=Chat, Shat=Shat, Ohat=Ohat, s_max=float(s_arr.max()), e_grid=e_grid.tolist())
    # 残差：級数で再構成した ξ(θ_X) と求積の差、ε̄ の一周内の揺らぎ（検査する軌道 = 初期の軌道 e_K0。δ=0 なら格子の最大 e）
    r0 = rows[-1]
    T["e_K0"] = e_K0
    T["A0"] = 0.0 if circular0 else r0["A"]
    T["s0"] = 0.0 if circular0 else r0["s"]
    T["check_e"] = r0["e"]
    B_test = r0["A"] * np.exp(1j * np.arange(Nsamp) * 2.0 * math.pi / Nsamp)
    xi_rec = np.empty(Nsamp)
    pi_rec = np.empty(Nsamp)
    for i in range(Nsamp):
        xi_rec[i], pi_rec[i], _, _, _ = breathing_eval(C, T, B_test[i], xi_c0, zeta_c0, omX0)
    thX = np.arange(Nsamp) * 2.0 * math.pi / Nsamp
    th_half = np.where(thX <= math.pi, thX, 2.0 * math.pi - thX)
    # 再求積（最後の行 = e_K0）
    sol = None
    xi_p = r0["xi_p"]
    T_est = 2.0 * math.pi / omX0
    sol = solve_ivp(rhs, (0.0, 1.5 * T_est), [xi_p, 0.0], method="DOP853", rtol=1e-13, atol=1e-18,
                    events=ev_apo, dense_output=True, max_step=T_est / 200)
    T_half = float(sol.t_events[0][0])
    yv = sol.sol(th_half / (math.pi / T_half))
    xi_ex = yv[0]
    T["fit_resid_xi"] = float(np.max(np.abs(xi_rec - xi_ex)) / xi_c0)
    eb = np.array([ebar_and_grad(C, Mom, xi_rec[i], pi_rec[i])[0] for i in range(Nsamp)])
    T["energy_wobble"] = float((eb.max() - eb.min()) / abs(eb.mean()))
    # 必要次数の確認：q^{N+1} と読み出し基準
    T["q0_pow_Np1"] = float(math.sqrt(T["s0"]) ** (N + 1))
    return T


# =============================================================================
# 状態
# =============================================================================

class State:
    pass


def init_state(C, delta=0.065, f=1.0, sa=(0.0, 0.0, 1.0), sb=(0.0, 0.0, -1.0)):
    """§8 初期化。ここでは割り算・反復・求積を使ってよい（E–M の表を一度解いて展開する）。"""
    S = State()
    S.f = float(f)
    S.sa = np.array(sa, dtype=float)
    S.sb = np.array(sb, dtype=float)
    sa_z, sb_z = S.sa[2], S.sb[2]
    # §8 手順 2 ξ_1：単一倍音 m=1、π=0、走行のスピン値で ∂ε/∂ξ = 0（§3 の表はスピン 0 の値）
    psi1 = np.zeros(C.nm, dtype=complex)
    psi1[C.M + 1] = 1.0
    xi_c1, _, _, _ = node_constants(C, psi1, sa_z, sb_z, C.kappa / (4.0 * C.g ** 2))
    z = 1.0 / xi_c1
    a0 = (C.A == 0)
    Mom1, _ = moments_of(C, psi1, sa_z, sb_z)
    for _ in range(6):
        Pz = -np.sum(C.P[a0] * C.coef[a0] * Mom1[a0] * z ** (C.P[a0] + 1.0))
        dP = -np.sum(C.P[a0] * (C.P[a0] + 1.0) * C.coef[a0] * Mom1[a0] * z ** C.P[a0])
        z = z - Pz / dP
    S.xi1 = 1.0 / z
    S.xi1_residual = float(-np.sum(C.P[a0] * C.coef[a0] * Mom1[a0] * z ** (C.P[a0] + 1.0)))
    # §8 手順 4：倍音分布（台 {0,1,2}）、⟨m²⟩ = (1+δ)²
    S.delta = float(delta)
    m2 = (1.0 + delta) ** 2
    w02 = 0.5 * (m2 - 1.0)
    psi = np.zeros(C.nm, dtype=complex)
    psi[C.M + 0] = math.sqrt(w02)
    psi[C.M + 1] = math.sqrt(1.0 - 2.0 * w02)
    psi[C.M + 2] = math.sqrt(w02)
    S.psi = psi
    # 手順 0 の定数
    S.xi_c, S.zeta_c, S.omX0, S.mom = node_constants(C, psi, sa_z, sb_z, C.kappa / (4.0 * C.g ** 2 * m2))
    Mom, _ = moments_of(C, psi, sa_z, sb_z)
    # §8 手順 5：呼吸の軌道を初期化で解いて展開（pericenter = ξ_1）
    S.T = build_breathing_table(C, Mom, S.xi_c, S.zeta_c, S.omX0, S.xi1)
    S.B = complex(S.T["A0"], 0.0)                      # θ_X,0 = 0（pericenter）
    xi0, pi0, omX, J, s0 = breathing_eval(C, S.T, S.B, S.xi_c, S.zeta_c, S.omX0)
    S.omX_inv = 1.0 / omX
    S.zeta = 1.0 / xi0
    S.det_inv = 1.0 / (J[0] * J[3] - J[1] * J[2])
    S.theta_X = 0.0
    S.theta_a = 0.0
    S.theta_b = 0.0
    S.k = 0
    S.histD = []
    S.histQ = []
    S.cum_thY = 0.0
    S.cum_phase1 = 0.0
    S.last = {}
    S.xi0_check = xi0
    return S


# =============================================================================
# §2.1 節点写像（相互作用の唯一の関数。割り算なし）
# =============================================================================

def node_map(C, S):
    psi = S.psi
    f = S.f
    sa_z, sb_z = S.sa[2], S.sb[2]
    T = S.T
    xi_c, zeta_c, omX0 = S.xi_c, S.zeta_c, S.omX0
    # 手順 2：呼吸の級数から ξ_k、π_k、ω_X,k（積和）
    xi_k, pi_k, omX, J, s_k = breathing_eval(C, T, S.B, xi_c, zeta_c, omX0)
    # 手順 1：Δθ_k = Δθ_X ω_X^{-1}（中間量）
    omX_inv = recip(omX, S.omX_inv, n=3)
    dth = C.dthX * omX_inv
    dth_inv = omX * C.dthX_inv
    zeta = recip(xi_k, S.zeta, n=4)
    # ε_m、∂ε/∂m
    eps_m, deps_dm = eps_modes(C, zeta, pi_k, sa_z, sb_z)
    wgt = np.abs(psi) ** 2
    norm = float(np.sum(wgt))
    # 手順 3：四つの位相の進み
    vm = 2.0 * C.g * C.marr * zeta
    vbar2 = pi_k * pi_k + float(np.sum(wgt * vm * vm))
    va2 = vbar2 * C.g_inv2
    vb2 = vbar2 * C.r_over_g2
    Ua = C.Ua_coef * zeta
    Ub = C.Ub_coef * zeta
    Fa = 1.0 - 0.5 * va2 - Ua - 0.125 * va2 * va2
    Fb = 1.0 - 0.5 * vb2 - Ub
    rho = C.mu_over_2ma * float(np.sum(wgt * deps_dm))
    dthY = rho * dth
    dtha = Fa * dth
    dthb = C.mb_over_ma * Fb * dth
    # 手順 4：放射反作用（直前 4・6 節点の後退差分）と倍音の更新
    S1, S2 = quad_relations(psi)
    D = xi_k * S1
    Q = xi_k * xi_k * S2
    S.histD.append(D)
    S.histQ.append(Q)
    if len(S.histD) > 4:
        S.histD.pop(0)
    if len(S.histQ) > 6:
        S.histQ.pop(0)
    D3 = 0.0 + 0.0j
    Q5 = 0.0 + 0.0j
    D2 = 0.0 + 0.0j
    if len(S.histD) == 4:
        h = S.histD
        D3 = (h[3] - 3.0 * h[2] + 3.0 * h[1] - h[0]) * dth_inv ** 3
        D2 = (h[3] - 2.0 * h[2] + h[1]) * dth_inv ** 2
    if len(S.histQ) == 6:
        q = S.histQ
        Q5 = (q[5] - 5.0 * q[4] + 10.0 * q[3] - 10.0 * q[2] + 5.0 * q[1] - q[0]) * dth_inv ** 5
    Vpsi = (-f * C.rr_em * xi_k * (np.conj(D3) * shift(psi, -1) + D3 * shift(psi, +1))
            + f * C.rr_gw * xi_k * xi_k * (Q5 * shift(psi, +2) + np.conj(Q5) * shift(psi, -2)))
    phase = cexp_series(-1j * C.mu_over_2ma * eps_m * dth, C)
    psi_new = (psi - 1j * C.mu_over_2ma * dth * Vpsi) * phase
    # 動径力（手順 6）
    Fr = f * (C.fr_em * (np.conj(D3) * S1).real - C.fr_gw * xi_k * (Q5 * np.conj(S2)).real)
    # 手順 5：新しい定数
    xi_c_new, zeta_c_new, omX0_new, mom_new = node_constants(C, psi_new, sa_z, sb_z, zeta_c)
    # 手順 6：接触軌道要素の更新。目標 (ξ_k, π_k + F_r Δθ_k) を新しい定数で満たす B を Newton 2 歩（注記 (a)）
    xi_t = xi_k
    pi_t = pi_k + Fr * dth
    Bn = S.B
    det_inv = S.det_inv
    for _ in range(2):
        xi_e, pi_e, _, Je, _ = breathing_eval(C, T, Bn, xi_c_new, zeta_c_new, omX0_new)
        J11, J12, J21, J22 = Je
        det = J11 * J22 - J12 * J21
        det_inv = recip(det, det_inv, n=3)
        rx, ry = xi_e - xi_t, pi_e - pi_t
        dx = -(J22 * rx - J12 * ry) * det_inv
        dy = -(-J21 * rx + J11 * ry) * det_inv
        Bn = complex(Bn.real + dx, Bn.imag + dy)
    B_new = Bn * C.eidthX
    # 手順 7：内部位相・スピン
    S.theta_X += C.dthX
    S.theta_a += dtha
    S.theta_b += dthb
    m1 = float(np.dot(C.mpow[1], wgt))
    z3 = zeta ** 3
    nn = nn_tensor(S2)
    zhat = np.array([0.0, 0.0, 1.0])
    Om_a = (m1 * z3 * (C.two_alpha * (2 * C.g - 1) + C.two_alpha_over_C0m1 * (4 * C.g - 1))) * zhat \
        - z3 * (C.two_alpha * (C.g - 1) + C.two_alpha_ma2_over9) * (S.sb - 3.0 * np.dot(nn, S.sb))
    Om_b = (m1 * z3 * (C.two_alpha * ((C.g - 1) ** 2 + 2 * (C.g - 1))
                       + C.two_alpha_over_C0m1 * (4 * (C.g - 1) + 3 * (C.g - 1) ** 2))) * zhat \
        - z3 * (C.two_alpha * (C.g - 1) + C.two_alpha_ma2_over9) * (S.sa - 3.0 * np.dot(nn, S.sa))
    sa_new = S.sa + np.cross(Om_a, S.sa) * dth
    sb_new = S.sb + np.cross(Om_b, S.sb) * dth
    sa_new = sa_new * (1.5 - 0.5 * float(sa_new @ sa_new))
    sb_new = sb_new * (1.5 - 0.5 * float(sb_new @ sb_new))
    # 手順 8：何もしない
    # ---- 読み出し用の記録（発展には使わない）----
    ebar = float(np.sum(wgt * eps_m))
    eloss_em = -f * C.eloss_em * abs(D2) ** 2
    S.last = dict(k=S.k, theta_X=S.theta_X - C.dthX, dth=dth, dthY=dthY, dtha=dtha, dthb=dthb,
                  xi=xi_k, pi=pi_k, xi_c=xi_c, A=abs(S.B), q=abs(S.B) * zeta_c, s=s_k,
                  argB_minus_thX=cmath.phase(S.B) - (S.theta_X - C.dthX),
                  omX=omX, rho=rho, gamma=rho * omX_inv, Fa=Fa, Fb=Fb, v_over_c=xi_k * rho,
                  S1_abs=abs(S1), S1_arg=cmath.phase(S1), S2_abs=abs(S2), ebar=ebar, norm=norm,
                  m1=m1, m2=float(np.dot(C.mpow[2], wgt)), eloss_em=eloss_em, Fr=Fr,
                  arg_psi1=cmath.phase(psi[C.M + 1]), E1dth=-C.mu_over_2ma * eps_m[C.M + 1] * dth)
    S.cum_thY += dthY
    S.cum_phase1 += -C.mu_over_2ma * eps_m[C.M + 1] * dth
    # 状態の書き込み
    S.psi = psi_new
    S.B = B_new
    S.xi_c, S.zeta_c, S.omX0, S.mom = xi_c_new, zeta_c_new, omX0_new, mom_new
    S.omX_inv = omX_inv
    S.zeta = zeta
    S.det_inv = det_inv
    S.sa, S.sb = sa_new, sb_new
    S.k += 1
    return S


# =============================================================================
# §2.2・§7 読み出し
# =============================================================================

READOUT_FIELDS = ["k", "theta_X", "cum_thY_over_2pi", "theta_a", "theta_b", "dth", "dthY", "dtha", "dthb",
                  "xi", "pi", "xi_c", "A", "q", "argB_minus_thX", "omX", "rho", "gamma", "Fa", "Fb",
                  "v_over_c", "S1_abs", "S1_arg", "S2_abs", "ebar", "norm", "m1", "m2", "eloss_em", "Fr",
                  "arg_psi1", "cum_phase1", "sa_z", "sb_z"]


def readout(C, S):
    L = S.last
    row = dict(L)
    row["cum_thY_over_2pi"] = S.cum_thY / (2.0 * math.pi)
    row["theta_a"] = S.theta_a
    row["theta_b"] = S.theta_b
    row["cum_phase1"] = S.cum_phase1
    row["sa_z"] = S.sa[2]
    row["sb_z"] = S.sb[2]
    return {k: row[k] for k in READOUT_FIELDS}


def psi_weights(C, S):
    return (np.abs(S.psi) ** 2).tolist()


# =============================================================================
# 自己検査
# =============================================================================

def eps_direct(C, xi, pi_, m, sa, sb):
    kC, kG, nu, g = C.kC, C.kG, C.nu, C.g
    r = g - 1.0
    vm = 2.0 * m * g / xi
    P2 = pi_ ** 2 + vm ** 2
    e = P2 / 2 - kC / xi - kG / xi
    e += (3 * nu - 1) * P2 ** 2 / 8
    e += -nu * kC * (P2 + pi_ ** 2) / (2 * xi)
    e += -kG * ((3 + nu) * P2 + nu * pi_ ** 2) / (2 * xi) + kG ** 2 / (2 * xi ** 2)
    e += 1.5 * kC * kG / xi ** 2
    e += -kC / (8 * xi) * (3 * nu ** 2 * pi_ ** 4 + nu * (3 * nu - 2) * P2 ** 2 + 2 * nu * (nu - 1) * P2 * pi_ ** 2)
    e += (1 - 5 * nu + 5 * nu ** 2) * P2 ** 3 / 16 + nu * kC ** 2 * P2 / (4 * xi ** 2) + nu * kC ** 3 / (4 * xi ** 3)
    e += kC * (1 + 2 * r) * m * sa / xi ** 3 + kC * (r ** 2 + 2 * r) * m * sb / xi ** 3
    e += kG * (3 + 4 * r) * m * sa / xi ** 3 + kG * (4 * r + 3 * r ** 2) * m * sb / xi ** 3
    e += -kC * r * sa * sb / xi ** 3 - kG * r * sa * sb / xi ** 3
    e += -kC / (2 * xi ** 3) - kC * r ** 2 / (2 * xi ** 3) - kG / (2 * xi ** 3) - kG * r ** 2 / (2 * xi ** 3)
    return e


def selftest(C):
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(200):
        xi = rng.uniform(100.0, 2000.0)
        pi_ = rng.uniform(-0.01, 0.01)
        sa = rng.choice([-1.0, 0.0, 1.0])
        sb = rng.choice([-1.0, 0.0, 1.0])
        eps_m, _ = eps_modes(C, 1.0 / xi, pi_, sa, sb)
        for m in (-3, -1, 0, 1, 2, 5):
            d = eps_direct(C, xi, pi_, m, sa, sb)
            t = eps_m[C.M + m]
            worst = max(worst, abs(t - d) / max(abs(d), 1e-30))
    err_exp = 0.0
    for x in (0.05, 0.2, 1.0, 10.0, 60.0, 120.0):
        z = -1j * x
        err_exp = max(err_exp, abs(cexp_series(np.array([z]), C)[0] - cmath.exp(z)))
    return worst, err_exp


def mode_phase_per_node(C, S):
    eps_m, _ = eps_modes(C, 1.0 / S.xi1, 0.0, S.sa[2], S.sb[2])
    dth = C.dthX / S.omX0
    return np.abs(C.mu_over_2ma * eps_m * dth)


# =============================================================================
# 実行
# =============================================================================

def run(args):
    C = build_constants(M=args.M, nterms_exp=args.nterms_exp, n_harm=args.n_harm, poly_deg=args.poly_deg)
    if args.selftest:
        worst, err_exp = selftest(C)
        print(f"selftest: monomial table vs direct T1..T22  max rel diff = {worst:.3e}")
        print(f"selftest: exp series ({C.nterms_exp} terms, 2^{C.exp_nsq} scaling) max abs err for |z| ≤ 120 = {err_exp:.3e}")
        print(f"constants: kC={C.kC:.9e} kG={C.kG:.9e} nu={C.nu:.9e} g={C.g:.9f}")
        S = init_state(C, delta=0.065, f=0.0, sa=(0, 0, 1), sb=(0, 0, -1))
        T = S.T
        print(f"breathing table (δ=0.065): e_K0={T['e_K0']:.6f} A0={T['A0']:.5f} s0={T['s0']:.6f} "
              f"fit residual ξ/ξ̄_c={T['fit_resid_xi']:.2e}  ε̄ wobble along series={T['energy_wobble']:.2e}  "
              f"q0^(N+1)={T['q0_pow_Np1']:.2e} (criterion {C.readout_precision:.1e})  ξ(θ_X=0)−ξ1={S.xi0_check - S.xi1:.2e}")
        return
    sa = {"up": (0, 0, 1), "down": (0, 0, -1), "zero": (0, 0, 0)}[args.sa]
    sb = {"up": (0, 0, 1), "down": (0, 0, -1), "zero": (0, 0, 0)}[args.sb]
    S = init_state(C, delta=args.delta, f=args.f, sa=sa, sb=sb)
    nodes = int(round(args.orbits * 62)) if args.nodes is None else args.nodes
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    ph = mode_phase_per_node(C, S)
    m_alias = int(np.min(np.abs(C.marr[ph > math.pi]))) if np.any(ph > math.pi) else None
    T = S.T
    init_info = dict(delta=args.delta, f=args.f, sa=args.sa, sb=args.sb, nodes=nodes, M=C.M, n_harm=C.n_harm, poly_deg=C.poly_deg,
                     xi1=S.xi1, xi1_residual=S.xi1_residual, xi_c0=S.xi_c, A0=abs(S.B), q0=abs(S.B) * S.zeta_c,
                     e_K0=T["e_K0"], omX0=S.omX0, psi_w0=psi_weights(C, S), dthX=C.dthX,
                     breathing_fit_resid_xi=T["fit_resid_xi"], breathing_energy_wobble=T["energy_wobble"],
                     q0_pow_Np1=T["q0_pow_Np1"], readout_precision=C.readout_precision,
                     phase_per_node_m1=float(ph[C.M + 1]), phase_per_node_mM=float(ph[-1]),
                     smallest_m_with_phase_over_pi=m_alias)
    t0 = time.time()
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
    elapsed = time.time() - t0
    gamma_mean = S.cum_thY / (nodes * C.dthX)
    a_cycles_per_orbit = (S.theta_a / (2 * math.pi)) / max(S.cum_thY / (2 * math.pi), 1e-12)
    summary = dict(init=init_info, first=first, last=last, psi_w_final=psi_weights(C, S),
                   gamma_mean=gamma_mean, a_cycles_per_orbit=a_cycles_per_orbit,
                   theta_b_over_theta_a=(S.theta_b / S.theta_a if S.theta_a else None),
                   psi1_phase_per_orbit_over_pi=(S.cum_phase1 / math.pi) / max(S.cum_thY / (2 * math.pi), 1e-12),
                   elapsed_s=elapsed, nodes_per_s=nodes / max(elapsed, 1e-9))
    with open(os.path.splitext(args.out)[0] + "_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=float)
    print(f"done: {nodes} nodes in {elapsed:.2f} s ({nodes / max(elapsed, 1e-9):.0f} nodes/s) -> {args.out}")
    print(f"  init : xi1={S.xi1:.6f} (residual {S.xi1_residual:.1e}) xi_c0={init_info['xi_c0']:.6f} A0={init_info['A0']:.4f} "
          f"q0={init_info['q0']:.5f} e_K0={init_info['e_K0']:.5f} omX0={init_info['omX0']:.6e}")
    print(f"         breathing series: fit resid {T['fit_resid_xi']:.1e}, ε̄ wobble {T['energy_wobble']:.1e}, q0^(N+1)={T['q0_pow_Np1']:.1e}; "
          f"phase/node m=1 {init_info['phase_per_node_m1']:.4f} rad, m={C.M} {init_info['phase_per_node_mM']:.1f} rad, first |m| with phase>pi: {m_alias}")
    print(f"  last : orbits(thY/2pi)={last['cum_thY_over_2pi']:.4f} |S1|={last['S1_abs']:.5f} |S2|={last['S2_abs']:.5f} "
          f"ebar={last['ebar']:.9e} norm-1={last['norm'] - 1:.2e} m2={last['m2']:.7f} q={last['q']:.6f}")
    print(f"  mean : gamma_mean-1={gamma_mean - 1:.4e}  a cycles/orbit={a_cycles_per_orbit:.3f}  theta_b/theta_a={summary['theta_b_over_theta_a']:.7f}  "
          f"psi1 phase/orbit/pi={summary['psi1_phase_per_orbit_over_pi']:.7f}")


def main():
    ap = argparse.ArgumentParser(description="万能相互作用 E–M 二体 節点写像（最終版 v6 仕様）")
    ap.add_argument("--delta", type=float, default=0.065)
    ap.add_argument("--f", type=float, default=1.0, help="吸収率 f（B1=0, B2=設定値, B3=1）")
    ap.add_argument("--sa", choices=["up", "down", "zero"], default="up")
    ap.add_argument("--sb", choices=["up", "down", "zero"], default="down")
    ap.add_argument("--orbits", type=float, default=10.0, help="走行長（節点数 = 62 × orbits。orbits は呼吸の周期数）")
    ap.add_argument("--nodes", type=int, default=None)
    ap.add_argument("--log-every", type=int, default=1)
    ap.add_argument("--M", type=int, default=31)
    ap.add_argument("--nterms-exp", type=int, default=14)
    ap.add_argument("--n-harm", type=int, default=12, help="呼吸の Fourier 級数の倍音数")
    ap.add_argument("--poly-deg", type=int, default=5, help="係数の s 多項式の次数")
    ap.add_argument("--out", default="results/run.csv")
    ap.add_argument("--selftest", action="store_true")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
