#!/usr/bin/env python3
"""共通エンジン：系に依らない部分。準位表と辺（A 係数）を受け取り、生成行列・厳密解・時間発展・吸収体の三層を計算する。

系の部品（system_*.py）が与えるもの
  states  : list of dict(label, E_eV, J, M_kg[, 任意の付加情報])   E_eV は準位エネルギー、J は全角運動量（原子なら F、分子なら J）、M は反跳用の質量
  channels: list of (src, dst, kind, k, A)                          src → dst の自発放出率 A [s⁻¹]、kind は文字列、k は多重極の階数（E1・M1 = 1、E2・GW = 2、二光子 = 0）
エンジンが与えるもの
  generator(states, channels, T, beta, a2)      準位分解の生成行列 R(T)。下向き A(1+⟨n̄⟩)、上向き A ⟨n̄⟩ g_src/g_dst（g = 2J+1）。重力波と二光子に温度は掛けない
  generator_m(states, channels, T, beta, a2)    m 分解（Wigner–Eckart）の生成行列と m 基底。Δm = q の型で重みづけた ⟨n̄⟩_q、副準位は縮退なしの詳細釣り合い
  exact(R, p0)                                  滞在時間 τ = −R_TT⁻¹ p₀、終状態 p_∞（吸収状態がある場合）
  evolve(R, p0, t_max, n_times)                 Radau の時間発展
  fluxes(states, channels, tau)                 種別ごとの放出エネルギー Σ A τ ΔE
  simulate_atom(...)                            層 1：方向の印・重心速度・反跳の量子跳躍モンテカルロ（層 2 の温度を含む）
  角運動量代数 wigner3j / wigner6j / w3j000_sq、角分布 pattern(k, q, x)
出典：固定一式 ../exchange_rel_20261005/exchange_cascade.py と ../exchange_general_20261006/ の absorber_*.py から切り出し。
      regression_hydrogen.py で元の数値に戻ることを確認する。
"""
from fractions import Fraction
from math import factorial, sqrt
import numpy as np
from scipy.integrate import solve_ivp

C = 299792458.0
HBAR_SI = 1.054571817e-34
E_SI = 1.602176634e-19
ME_SI = 9.1093837139e-31
K_B_EV = 8.617333262e-5
HBAR_EVS = HBAR_SI / E_SI
GL_X, GL_W = np.polynomial.legendre.leggauss(96)
RANK_DEFAULT = {"e1": 1, "m1": 1, "e2": 2, "gw": 2, "2ph": 0}
THERMAL_KINDS = ("e1", "e2", "m1")      # 温度を掛ける種別（光子）。重力波・二光子には掛けない


# ---------- 角運動量代数（Racah、有理数で厳密） ----------
def _fact(x):
    return factorial(int(x))


def _tri(a, b, c):
    return (a + b + c).denominator == 1 and abs(a - b) <= c <= a + b


def _Delta(a, b, c):
    return Fraction(_fact(a + b - c) * _fact(a - b + c) * _fact(-a + b + c), _fact(a + b + c + 1))


def wigner6j(a, b, c, d, e, f):
    a, b, c, d, e, f = (Fraction(x) for x in (a, b, c, d, e, f))
    for tri in ((a, b, c), (a, e, f), (d, b, f), (d, e, c)):
        if not _tri(*tri):
            return 0.0
    pref = _Delta(a, b, c) * _Delta(a, e, f) * _Delta(d, b, f) * _Delta(d, e, c)
    kmin = int(max(a + b + c, a + e + f, d + b + f, d + e + c))
    kmax = int(min(a + b + d + e, a + c + d + f, b + c + e + f))
    s = Fraction(0)
    for k in range(kmin, kmax + 1):
        s += (-1) ** k * Fraction(_fact(k + 1), _fact(k - (a + b + c)) * _fact(k - (a + e + f)) * _fact(k - (d + b + f))
                                  * _fact(k - (d + e + c)) * _fact((a + b + d + e) - k) * _fact((a + c + d + f) - k) * _fact((b + c + e + f) - k))
    return float(np.sign(float(s)) * sqrt(float(pref * s * s)))


def wigner3j(j1, j2, j3, m1, m2, m3):
    j1, j2, j3, m1, m2, m3 = (Fraction(x) for x in (j1, j2, j3, m1, m2, m3))
    if m1 + m2 + m3 != 0 or not (abs(j1 - j2) <= j3 <= j1 + j2) or (j1 + j2 + j3).denominator != 1:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    for x in (j1 + m1, j1 - m1, j2 + m2, j2 - m2, j3 + m3, j3 - m3):
        if x.denominator != 1:
            return 0.0
    D = Fraction(_fact(j1 + j2 - j3) * _fact(j1 - j2 + j3) * _fact(-j1 + j2 + j3), _fact(j1 + j2 + j3 + 1))
    pref = D * _fact(j1 + m1) * _fact(j1 - m1) * _fact(j2 + m2) * _fact(j2 - m2) * _fact(j3 + m3) * _fact(j3 - m3)
    kmin = int(max(0, j2 - j3 - m1, j1 - j3 + m2))
    kmax = int(min(j1 + j2 - j3, j1 - m1, j2 + m2))
    s = Fraction(0)
    for k in range(kmin, kmax + 1):
        s += Fraction((-1) ** k, _fact(k) * _fact(j3 - j2 + k + m1) * _fact(j3 - j1 + k - m2) * _fact(j1 + j2 - j3 - k) * _fact(j1 - k - m1) * _fact(j2 - k + m2))
    return sqrt(float(pref)) * float(s) * (-1) ** int(j1 - j2 - m3)


def w3j000_sq(l1, l2, l3):
    J = l1 + l2 + l3
    if J % 2 or not (abs(l1 - l2) <= l3 <= l1 + l2):
        return 0.0
    g = J // 2
    D = Fraction(_fact(J - 2 * l1) * _fact(J - 2 * l2) * _fact(J - 2 * l3), _fact(J + 1))
    return float(D * Fraction(_fact(g), _fact(g - l1) * _fact(g - l2) * _fact(g - l3)) ** 2)


# ---------- 放射場 ----------
def planck_nbar(hw_eV, T):
    if T <= 0.0 or hw_eV <= 0.0:
        return 0.0
    x = hw_eV / (K_B_EV * T)
    return 0.0 if x > 700.0 else 1.0 / np.expm1(x)


def pattern(k, q, x):
    """Δm = q の放出の角分布 P_q(θ)（偏光で和、∫P dΩ = 1）。x = cos θ。"""
    if k == 1:
        return (3 / (8 * np.pi)) * (1 - x * x) if q == 0 else (3 / (16 * np.pi)) * (1 + x * x)
    if k == 2:
        if q == 0:
            return (15 / (8 * np.pi)) * (1 - x * x) * x * x
        if abs(q) == 1:
            return (5 / (16 * np.pi)) * (1 - 3 * x * x + 4 * x**4)
        return (5 / (16 * np.pi)) * (1 - x**4)
    return 1 / (4 * np.pi)


def T_of_x(T0, beta, a2, x):
    """原子の静止系で見た背景の温度。x = cos θ（印の軸 ẑ に対する源の方向）。双極子は Doppler、a₂ は背景自身の四重極。"""
    g = 1.0 / np.sqrt(1.0 - beta**2)
    return T0 / (g * (1.0 - beta * x)) * (1.0 + a2 * 0.5 * (3 * x * x - 1))


def nbar_avg(hw_eV, T0, beta=0.0, a2=0.0):
    """角度平均 ⟨n̄⟩_Ω。"""
    if T0 <= 0.0:
        return 0.0
    vals = np.array([planck_nbar(hw_eV, T_of_x(T0, beta, a2, x)) for x in GL_X])
    return 0.5 * float(np.dot(GL_W, vals))


def nbar_q(hw_eV, T0, beta, a2, k, q):
    """Δm = q の型で重みづけた占有数 ⟨n̄⟩_q = ∫ n̄ P_q dΩ。"""
    if T0 <= 0.0:
        return 0.0
    vals = np.array([planck_nbar(hw_eV, T_of_x(T0, beta, a2, x)) * pattern(k, q, x) for x in GL_X])
    return 2 * np.pi * float(np.dot(GL_W, vals))


# ---------- 準位分解の生成行列と厳密解 ----------
def generator(states, channels, T=0.0, beta=0.0, a2=0.0, thermal_kinds=THERMAL_KINDS):
    n = len(states)
    R = np.zeros((n, n))
    rows = []
    for src, dst, kind, k, A in channels:
        if A <= 0.0:
            continue
        hw = states[src]["E_eV"] - states[dst]["E_eV"]
        nb = nbar_avg(hw, T, beta, a2) if kind in thermal_kinds else 0.0
        down = A * (1.0 + nb)
        g_src, g_dst = 2 * states[src]["J"] + 1, 2 * states[dst]["J"] + 1
        up = A * nb * g_src / g_dst
        R[dst, src] += down
        R[src, src] -= down
        if up > 0.0:
            R[src, dst] += up
            R[dst, dst] -= up
        rows.append((src, dst, kind, k, A, nb, up, hw))
    return R, rows


def exact(R, p0):
    """τ = −R_TT⁻¹ p₀（過渡状態）、p_∞ = p₀_A + R_AT τ（吸収状態）。吸収状態が無い（T > 0）なら τ は定義されず None。"""
    n = R.shape[0]
    out_rate = -np.diag(R)
    trans = [i for i in range(n) if out_rate[i] > 0.0]
    absb = [i for i in range(n) if out_rate[i] == 0.0]
    if not absb:
        return None, None
    tau = np.zeros(n)
    tau[trans] = np.linalg.solve(-R[np.ix_(trans, trans)], p0[trans])
    p_inf = p0.copy()
    p_inf[trans] = 0.0
    p_inf[absb] += np.dot(R[np.ix_(absb, trans)], tau[trans])
    return tau, p_inf


def fluxes(states, channels, tau):
    out = {}
    for src, dst, kind, k, A in channels:
        out[kind] = out.get(kind, 0.0) + A * tau[src] * (states[src]["E_eV"] - states[dst]["E_eV"])
    return out


def evolve(R, p0, t_max=1.0e17, n_times=300, rtol=1e-11):
    times = np.concatenate([[0.0], np.logspace(-12, np.log10(t_max), n_times)])
    sol = solve_ivp(lambda t, p: np.dot(R, p), (0.0, t_max), p0, method="Radau", jac=lambda t, p: R,
                    t_eval=times, rtol=rtol, atol=1e-30)
    if not sol.success:
        raise RuntimeError("時間発展が解けない：" + sol.message)
    return times, sol.y.T


# ---------- m 分解（層 3） ----------
def m_basis(states):
    """(state index, m) の一覧。"""
    out = []
    for i, s in enumerate(states):
        J = s["J"]
        for mm in range(-int(J), int(J) + 1):
            out.append((i, mm))
    return out


def channels_m(states, channels):
    """m 分解した辺：(src_m, dst_m, kind, k, q, A_mm')。A_mm' = A (2J+1) (J' k J; m' q −m)²、q = m − m'。"""
    mb = m_basis(states)
    midx = {s: i for i, s in enumerate(mb)}
    out = []
    worst = 0.0
    for src, dst, kind, k, A in channels:
        if A <= 0.0:
            continue
        Js, Jd = states[src]["J"], states[dst]["J"]
        for m in range(-int(Js), int(Js) + 1):
            tot = 0.0
            for m2 in range(-int(Jd), int(Jd) + 1):
                q = m - m2
                if abs(q) > k:
                    continue
                if k == 0:
                    Amm = A if (m2 == m and Jd == Js) else 0.0
                else:
                    Amm = A * (2 * Js + 1) * wigner3j(Jd, k, Js, m2, q, -m) ** 2
                if Amm > 0.0:
                    out.append((midx[(src, m)], midx[(dst, m2)], kind, k, q, Amm))
                    tot += Amm
            worst = max(worst, abs(tot - A) / A)
    return mb, out, worst


def generator_m(states, channels, T=0.0, beta=0.0, a2=0.0, thermal_kinds=THERMAL_KINDS):
    mb, edm, worst = channels_m(states, channels)
    n = len(mb)
    R = np.zeros((n, n))
    for s, d, kind, k, q, A in edm:
        hw = states[mb[s][0]]["E_eV"] - states[mb[d][0]]["E_eV"]
        nq = nbar_q(hw, T, beta, a2, k, q) if kind in thermal_kinds else 0.0
        R[d, s] += A * (1.0 + nq)
        R[s, s] -= A * (1.0 + nq)
        if nq > 0.0:
            R[s, d] += A * nq
            R[d, d] -= A * nq
    return R, mb, edm, worst


def line_photons_m(mb, edm, tau):
    """線（準位対）ごとの N_q。"""
    acc = {}
    for s, d, kind, k, q, A in edm:
        key = (mb[s][0], mb[d][0], kind)
        acc.setdefault(key, {})
        acc[key][q] = acc[key].get(q, 0.0) + A * tau[s]
    return acc


def polarization_90(Nq):
    N0 = Nq.get(0, 0.0); Ns = Nq.get(1, 0.0) + Nq.get(-1, 0.0)
    return (2 * N0 - Ns) / (2 * N0 + Ns) if (2 * N0 + Ns) > 0 else float("nan")


# ---------- 層 1：方向・重心速度・反跳（量子跳躍モンテカルロ） ----------
def _gamma_minus_1(b):
    g = 1.0 / np.sqrt(1.0 - b * b)
    return b * b * g * g / (g + 1.0), g


def doppler_delta(n_lab, v):
    """δ = ω_lab/ω′ − 1 = 1/(γ(1 − β n_∥)) − 1 を桁落ちなしに。"""
    speed = np.linalg.norm(v)
    if speed == 0.0:
        return 0.0
    b = speed / C
    gm1, g = _gamma_minus_1(b)
    npar = np.dot(n_lab, v / speed)
    return (b * npar * g - gm1) / (g * (1.0 - b * npar))


def boost_dir(n_rest, v):
    speed = np.linalg.norm(v)
    if speed == 0.0:
        return n_rest
    vhat = v / speed; b = speed / C; _, g = _gamma_minus_1(b)
    npar = np.dot(n_rest, vhat)
    n = ((npar + b) / (1.0 + b * npar)) * vhat + (n_rest - npar * vhat) / (g * (1.0 + b * npar))
    return n / np.linalg.norm(n)


def isotropic(rng):
    u = rng.uniform(-1, 1); ph = rng.uniform(0, 2 * np.pi); s = np.sqrt(1 - u * u)
    return np.array([s * np.cos(ph), s * np.sin(ph), u])


def _T_dir(T0, v, n_src, a2=0.0):
    b = np.linalg.norm(v) / C
    if b == 0.0:
        x = n_src[2]
        return T0 * (1.0 + a2 * 0.5 * (3 * x * x - 1))
    _, g = _gamma_minus_1(b)
    x = np.dot(n_src, v / np.linalg.norm(v))
    return T0 / (g * (1.0 - b * x)) * (1.0 + a2 * 0.5 * (3 * x * x - 1))


def _sample_dir_weighted(rng, hw, T0, v, a2, src_weight=True):
    b = np.linalg.norm(v) / C
    _, g = _gamma_minus_1(b)
    nmax = planck_nbar(hw, T0 / (g * (1.0 - b)) * (1.0 + abs(a2)))
    while True:
        n = isotropic(rng)
        n_src = n if src_weight else -n
        if rng.uniform(0, nmax) <= planck_nbar(hw, _T_dir(T0, v, n_src, a2)):
            return n


def apply_kick(v, M_before, hw_lab_eV, n_lab, sign):
    b = np.linalg.norm(v) / C
    _, g = _gamma_minus_1(b)
    hw_J = hw_lab_eV * E_SI
    P = g * M_before * v
    Eat = g * M_before * C**2
    P2 = P + sign * (hw_J / C) * n_lab
    E2 = Eat + sign * hw_J
    return P2 * C**2 / E2


def build_channels_mc(states, channels, T0, thermal_kinds=THERMAL_KINDS):
    """モンテカルロ用：準位ごとのチャネル。放出は静止系の厳密な二体運動学 ħω′ = ΔE(1 − ΔE/2Mc²)、吸収は ΔE(1 + ΔE/2M_f c²)。"""
    ch = {i: [] for i in range(len(states))}
    for src, dst, kind, k, A in channels:
        if A <= 0.0:
            continue
        dE = states[src]["E_eV"] - states[dst]["E_eV"]
        Ms, Md = states[src]["M_kg"] * C**2 / E_SI, states[dst]["M_kg"] * C**2 / E_SI    # [eV]
        ch[src].append(dict(dst=dst, kind=kind, k=k, A=A, hw=dE * (1.0 - dE / (2.0 * Ms)), up=False))
        if kind in thermal_kinds and T0 > 0.0:
            g_s, g_d = 2 * states[src]["J"] + 1, 2 * states[dst]["J"] + 1
            ch[dst].append(dict(dst=src, kind=kind + "_abs", k=k, A=A * g_s / g_d, hw=dE * (1.0 + dE / (2.0 * Md)), up=True))
    return ch


def _rates_mc(ch_i, T0, v, a2, thermal_kinds):
    out = []
    for c in ch_i:
        base_kind = c["kind"].replace("_abs", "")
        if c["up"]:
            r = c["A"] * nbar_avg(c["hw"], T0, np.linalg.norm(v) / C, a2)
        elif base_kind in thermal_kinds:
            r = c["A"] * (1.0 + nbar_avg(c["hw"], T0, np.linalg.norm(v) / C, a2))
        else:
            r = c["A"]
        out.append(r)
    return np.array(out)


def simulate_atom(states, channels, T0, v0, rng, init_index, n_events=40, t_cap=1.0e17, a2=0.0,
                  thermal_kinds=THERMAL_KINDS, two_photon_kind="2ph"):
    """一個の原子の事象列。各事象：dict(t, kind, src, dst, hw_lab(eV), n_lab, delta, v_before, v_after)。
    二光子：和は厳密、方向は (1 + cos²θ₁₂) の相関、分配は y³(1−y)³（形は近似）。"""
    ch = build_channels_mc(states, channels, T0, thermal_kinds)
    def width(i, v):
        return float(np.sum(_rates_mc(ch[i], T0, v, a2, thermal_kinds)))
    i = init_index
    v = np.array(v0, dtype=float)
    t = 0.0
    events = []
    while len(events) < n_events and t < t_cap:
        r = _rates_mc(ch[i], T0, v, a2, thermal_kinds)
        Rtot = r.sum()
        if Rtot <= 0.0:
            break
        t += rng.exponential(1.0 / Rtot)
        if t > t_cap:        # 上限より後の事象は記録しない（⟨n̄⟩ ≈ 0 の吸収が t ~ 1e135 s に 1 件入るのを防ぐ）
            break
        k = rng.choice(len(r), p=r / Rtot)
        c = ch[i][k]
        f = c["dst"]
        rel_width = 0.5 * (width(i, v) + width(f, v)) * HBAR_EVS / c["hw"]
        eps = rng.standard_cauchy() * rel_width
        v_before = v.copy()
        M_i = states[i]["M_kg"]
        if c["kind"] == two_photon_kind:
            hw_tot = c["hw"] * (1.0 + eps)
            y = rng.beta(4, 4)
            n1 = isotropic(rng)
            while True:
                n2 = isotropic(rng)
                if rng.uniform(0, 2) <= 1 + np.dot(n1, n2) ** 2:
                    break
            n1l, n2l = boost_dir(n1, v), boost_dir(n2, v)
            hw1 = y * hw_tot * (1.0 + doppler_delta(n1l, v))
            hw2 = (1 - y) * hw_tot * (1.0 + doppler_delta(n2l, v))
            v = apply_kick(v, M_i, hw1, n1l, -1.0)
            v = apply_kick(v, M_i - hw1 * E_SI / C**2, hw2, n2l, -1.0)
            events.append(dict(t=t, kind=two_photon_kind, src=i, dst=f, hw_lab=hw1 + hw2, n_lab=n1l, delta=float("nan"),
                               v_before=v_before, v_after=v.copy()))
            i = f
            continue
        if c["up"]:
            n_lab = -_sample_dir_weighted(rng, c["hw"], T0, v, a2, src_weight=True)
            sign, kind = +1.0, c["kind"]
        else:
            stim = False
            if c["kind"] in thermal_kinds and T0 > 0.0:
                nav = nbar_avg(c["hw"], T0, np.linalg.norm(v) / C, a2)
                stim = rng.uniform(0, 1 + nav) > 1.0
            n_lab = _sample_dir_weighted(rng, c["hw"], T0, v, a2, src_weight=False) if stim else boost_dir(isotropic(rng), v)
            sign, kind = -1.0, c["kind"] + ("_stim" if stim else "")
        Dm1 = doppler_delta(n_lab, v)
        delta = Dm1 + eps * (1.0 + Dm1)
        hw_lab = c["hw"] * (1.0 + delta)
        v = apply_kick(v, M_i, hw_lab, n_lab, sign)
        events.append(dict(t=t, kind=kind, src=i, dst=f, hw_lab=hw_lab, n_lab=n_lab, delta=delta,
                           v_before=v_before, v_after=v.copy(), line_hw=c["hw"], rel_width=rel_width))
        i = f
    return events
