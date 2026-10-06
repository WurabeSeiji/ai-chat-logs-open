#!/usr/bin/env python3
"""二体（電荷 q_a、q_b、質量 m_a、m_b）の多段放出。電荷の符号を入力にし、準位を閉じた式でなく
径方向 Dirac 方程式の数値解で作る。引力なら束縛準位、斥力なら束縛準位なし。他の構造（反跳、超微細、重力 1PN、
Einstein の A 係数、マスター方程式）は固定した一式 ../exchange_rel_20261005/exchange_cascade.py と同じで、
共通の関数はそこから import する（シリーズ内再現性規約：コピー → 対照 → import）。

対照実験：SIGN = −1（水素）で閉じた式の一式と準位・A 係数・放出を比較する（control_general_vs_closed_form.py）。

単位：解法の内部は ħ = c = 1、質量 = μ（換算質量）。長さの単位は λ̄_μ = ħ/(μc)、エネルギーは μc²。
      Bohr 半径 a = 1/α、Coulomb ポテンシャル V = s·α/r（s = sign(q_a q_b)、水素は −1）。
      出力は固定一式と同じ m_e c² 単位（ħ = 9/α、e² = 9、m_e = 1）に換算する。
"""

import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
FIXED = HERE.parent / "exchange_rel_20261005"
sys.path.insert(0, str(FIXED))
from matplotlib import font_manager  # noqa: E402
font_manager.fontManager.addfont = lambda path: None   # 固定一式の import 時のフォント登録を無効化（本体は変更しない）
import exchange_cascade as X  # noqa: E402

ALPHA, MU_P, HBAR, ME_EV = X.ALPHA, X.MU_P, X.HBAR, X.ME_EV
MU_RED = MU_P / (1.0 + MU_P)
SIGN = -1.0            # q_a q_b の符号。−1 引力（水素）、+1 斥力（e⁻e⁻、p⁺p⁺）
N_MAX = 5
LAM = (9.0 / ALPHA) / MU_RED        # λ̄_μ を固定一式の長さ単位で表した値（ħ/μ、m_e = 1、ħ = 9/α）
A_BOHR = 1.0 / ALPHA                # Bohr 半径 [λ̄_μ]


# ---------- 径方向 Dirac 方程式（点電荷、質量 1、ħ = c = 1） ----------
# ψ = (G Ω_κ / r, i F Ω_{−κ} / r)。G' = −κG/r + (ε + 1 − V)F、F' = −(ε − 1 − V)G + κF/r、V = s α / r。
# 原点：G, F ∝ r^γ、γ = √(κ² − α²)、F/G = −s(γ + κ)/α。無限遠（束縛）：F/G → −√((1−ε)/(1+ε))、減衰 e^{−√(1−ε²) r}。
# 束縛エネルギー η = 1 − ε を変数にする（ε ≈ 1 の桁落ちを避け、η を相対精度 1e-15 で決める）。
def _rhs_dirac(x, y, kappa, eta, s):
    r = np.exp(x)
    G, F = y
    V = s * ALPHA / r
    dG = -kappa * G + r * (2.0 - eta - V) * F        # (ε + 1 − V) = 2 − η − V
    dF = -r * (-eta - V) * G + kappa * F             # (ε − 1 − V) = −η − V
    return [dG, dF]          # d/dx = r d/dr


def _shoot(kappa, eta, s, r0, rm, rmax, dense=False):
    gam = np.sqrt(kappa**2 - ALPHA**2)
    y0 = [r0**gam, -s * (gam + kappa) / ALPHA * r0**gam]
    out = solve_ivp(_rhs_dirac, (np.log(r0), np.log(rm)), y0, args=(kappa, eta, s), method="DOP853",
                    rtol=1e-13, atol=1e-300, dense_output=dense)
    yinf = [1.0, -np.sqrt(eta / (2.0 - eta))]        # F/G → −√((1−ε)/(1+ε))
    inn = solve_ivp(_rhs_dirac, (np.log(rmax), np.log(rm)), yinf, args=(kappa, eta, s), method="DOP853",
                    rtol=1e-13, atol=1e-300, dense_output=dense)
    return out, inn


def _mismatch(eta, kappa, s, r0, rm, rmax):
    out, inn = _shoot(kappa, eta, s, r0, rm, rmax)
    Go, Fo = out.y[:, -1]
    Gi, Fi = inn.y[:, -1]
    no, ni = np.hypot(Go, Fo), np.hypot(Gi, Fi)
    return (Go * Fi - Fo * Gi) / (no * ni)


def dirac_bound_state(n, kappa, s=SIGN, npts=6000):
    """主量子数 n、Dirac 量子数 κ の束縛状態。(ε, r, G, F)。斥力（s = +1）では束縛解が無いので None。
    G（大成分）の節点数は n − ℓ − 1（非相対論の R_nℓ と同じ）。"""
    if s > 0:
        return None
    ell, _ = lj_of_kappa(kappa)
    eta_nr = ALPHA**2 / (2.0 * n**2)
    half_gap = ALPHA**2 / (3.0 * n**3)                      # 隣の n との間隔の 1/3。κ による微細構造 α⁴ ≪ これ
    r0 = 1e-10 * A_BOHR                                     # 原点側の打ち切り。∫₀^{r0} の欠落は相対 r0/a = 1e-10
    rm = max(0.7 * n**2 * A_BOHR, 2.0 * A_BOHR)
    rmax = 40.0 / np.sqrt(2.0 * (eta_nr - half_gap)) + 2 * rm   # 減衰 √(1−ε²) ≈ √(2η)
    lo, hi = eta_nr - half_gap, eta_nr + half_gap
    f_lo, f_hi = _mismatch(lo, kappa, s, r0, rm, rmax), _mismatch(hi, kappa, s, r0, rm, rmax)
    if f_lo * f_hi > 0:
        grid = np.linspace(lo, hi, 41)
        vals = [_mismatch(e, kappa, s, r0, rm, rmax) for e in grid]
        k = next(i for i in range(40) if vals[i] * vals[i + 1] <= 0)
        lo, hi = grid[k], grid[k + 1]
    eta = brentq(_mismatch, lo, hi, args=(kappa, s, r0, rm, rmax), xtol=1e-300, rtol=1e-15, maxiter=200)
    eps = 1.0 - eta
    out, inn = _shoot(kappa, eta, s, r0, rm, rmax, dense=True)
    x_out = np.linspace(np.log(r0), np.log(rm), npts // 2)
    x_in = np.linspace(np.log(rm), np.log(rmax), npts // 2 + 1)[1:]
    yo = out.sol(x_out)
    yi = inn.sol(x_in)
    yi_m = inn.sol(np.log(rm))                                # 接続点 rm ちょうどで倍率を決める（格子点ではなく）
    scale = yo[0, -1] / yi_m[0] if abs(yi_m[0]) > abs(yi_m[1]) else yo[1, -1] / yi_m[1]
    x = np.concatenate([x_out, x_in])
    G = np.concatenate([yo[0], scale * yi[0]])
    F = np.concatenate([yo[1], scale * yi[1]])
    r = np.exp(x)
    norm = simpson((G**2 + F**2) * r, x=x)                  # ∫(G²+F²)dr = ∫(G²+F²) r dx
    G, F = G / np.sqrt(norm), F / np.sqrt(norm)
    nodes = int(np.sum(G[:-1] * G[1:] < 0))
    if nodes != n - ell - 1:
        raise RuntimeError("節点数が合わない n=%d κ=%d nodes=%d（期待 %d）" % (n, kappa, nodes, n - ell - 1))
    return eta, r, G, F


def integ(r, f):
    """∫ f dr を対数格子で。"""
    return simpson(f * r, x=np.log(r))


# ---------- 径方向 Schrödinger 方程式（重力 1PN の期待値用。質量 μ、ħ = 1、V = s α/r） ----------
def _rhs_schr(x, y, ell, E, s):
    r = np.exp(x)
    u, du = y                      # u = rR、du = u'（r 微分）
    V = s * ALPHA / r
    ddu = (ell * (ell + 1) / r**2 + 2.0 * (V - E)) * u
    return [r * du, r * ddu]


def schrodinger_bound_state(n, ell, s=SIGN, npts=6000):
    """(E, r, R, dR)。E は μc² 単位の非相対論束縛エネルギー（負）。"""
    if s > 0:
        return None
    E_nr = -ALPHA**2 / (2.0 * n**2)
    half_gap = ALPHA**2 / (3.0 * n**3)
    r0 = 1e-10 * A_BOHR
    rm = max(0.7 * n**2 * A_BOHR, 2.0 * A_BOHR)
    rmax = 40.0 / np.sqrt(-2.0 * (E_nr + half_gap)) + 2 * rm

    def shoot(E, dense=False):
        y0 = [r0 ** (ell + 1), (ell + 1) * r0**ell]
        out = solve_ivp(_rhs_schr, (np.log(r0), np.log(rm)), y0, args=(ell, E, s), method="DOP853",
                        rtol=1e-13, atol=1e-300, dense_output=dense)
        k = np.sqrt(-2.0 * E)
        inn = solve_ivp(_rhs_schr, (np.log(rmax), np.log(rm)), [1.0, -k], args=(ell, E, s), method="DOP853",
                        rtol=1e-13, atol=1e-300, dense_output=dense)
        return out, inn

    def mism(E):
        out, inn = shoot(E)
        uo, duo = out.y[:, -1]
        ui, dui = inn.y[:, -1]
        return (uo * dui - duo * ui) / (np.hypot(uo, duo) * np.hypot(ui, dui))

    lo, hi = E_nr - half_gap, E_nr + half_gap
    if mism(lo) * mism(hi) > 0:
        grid = np.linspace(lo, hi, 41)
        vals = [mism(e) for e in grid]
        k = next(i for i in range(40) if vals[i] * vals[i + 1] <= 0)
        lo, hi = grid[k], grid[k + 1]
    E = brentq(mism, lo, hi, xtol=1e-16, rtol=1e-15, maxiter=200)
    out, inn = shoot(E, dense=True)
    x_out = np.linspace(np.log(r0), np.log(rm), npts // 2)
    x_in = np.linspace(np.log(rm), np.log(rmax), npts // 2 + 1)[1:]
    yo, yi = out.sol(x_out), inn.sol(x_in)
    yi_m = inn.sol(np.log(rm))
    scale = yo[0, -1] / yi_m[0]                                # 接続点 rm ちょうどで倍率を決める
    x = np.concatenate([x_out, x_in])
    r = np.exp(x)
    u = np.concatenate([yo[0], scale * yi[0]])
    du = np.concatenate([yo[1], scale * yi[1]])
    norm = simpson(u**2 * r, x=x)
    u, du = u / np.sqrt(norm), du / np.sqrt(norm)
    R = u / r
    dR = du / r - u / r**2
    return E, r, R, dR


# ---------- 準位 ----------
def kappa_list(n):
    """n に対する κ の一覧（ℓ = 0..n−1、j = ℓ ± ½）。κ = −(ℓ+1) で j = ℓ+½、κ = ℓ で j = ℓ−½（ℓ ≥ 1）。"""
    out = []
    for ell in range(n):
        out.append(-(ell + 1))
        if ell > 0:
            out.append(ell)
    return out


def lj_of_kappa(kappa):
    j = abs(kappa) - 0.5
    ell = kappa if kappa > 0 else -kappa - 1
    return ell, j


def gravity_from_expectations(n, ell, j, F, r1, r2, r3, psi0sq, E_nr, I_dR2r):
    """固定一式 gravity_level と同じ項を、数値の期待値から作る。入力は固定一式の単位（m_e = 1、ħ = 9/α、e² = 9）。
    r1 = ⟨1/r⟩、r2 = ⟨1/r²⟩、r3 = ⟨1/r³⟩（ℓ ≥ 1）、psi0sq = |ψ(0)|²（ℓ = 0）、E_nr = Schrödinger の束縛エネルギー、
    I_dR2r = ∫R'² r dr（⟨p_i n_in_j/r p_j⟩/ħ²）。"""
    mu = MU_RED
    g = MU_P * X.MU**2 / 9.0
    G = 9.0 * g / MU_P
    hb = HBAR
    kin = 2.0 * mu * (E_nr * r1 - SIGN * 9.0 * r2)          # ⟨(1/r)p²⟩ = 2μ(E⟨1/r⟩ − V の係数 ⟨1/r²⟩)
    p_r_p = kin - 2.0 * np.pi * hb**2 * psi0sq
    p_nn_p = hb**2 * I_dR2r
    newton = -9.0 * g * r1
    self_pn = -(3.0 * G / 4.0) * (MU_P + 1.0 / MU_P) * 2.0 * kin
    darwin = (3.0 * np.pi * G * hb**2 / 2.0) * (MU_P + 1.0 / MU_P) * psi0sq
    cross = -(G / 2.0) * (7.0 * p_r_p + p_nn_p)
    g2 = 9.0 * g * G * (1.0 + MU_P) / 2.0 * r2
    mixed = 13.5 * G * (1.0 + MU_P) * r2 * (-SIGN)           # 固定一式は引力（q_eq_p = −9）：−q_aq_b(m_a+m_b) = +9M
    LS_e = 0.5 * hb**2 * (j * (j + 1) - ell * (ell + 1) - 0.75)
    LJ = 0.5 * hb**2 * (j * (j + 1) + ell * (ell + 1) - 0.75)
    JI = 0.5 * hb**2 * (F * (F + 1) - j * (j + 1) - 0.75)
    LI = LJ * JI / (hb**2 * j * (j + 1))
    so_e = G * (2.0 + 1.5 * MU_P) * r3 * LS_e
    so_p = G * (2.0 + 1.5 / MU_P) * r3 * LI
    if ell > 0:
        ss = G * r3 * (ell * (ell + 1) - j * (j + 1) + 0.75) / (2.0 * j * (j + 1)) * JI
    else:
        ss = G * (8.0 * np.pi / 3.0) * psi0sq * JI
    return {"newton": newton, "self_pn": self_pn, "darwin": darwin, "cross": cross, "g2": g2, "mixed": mixed,
            "so_e": so_e, "so_p": so_p, "ss": ss}


def hyperfine_dirac(kappa, j, F, r, G, Fr):
    """相対論的超微細（点磁気双極子、Breit 1930）。E = (g_p α/(m_p/m_e)) · κ/(j(j+1)) · ∫ G F /r² dr · [F(F+1) − I(I+1) − j(j+1)]/2。
    導出：H = e α·(μ_p × r)/r³ の Dirac 行列要素。角度因子 κ/(j(j+1))、動径積分 ∫GF/r² dr。定数は非相対論極限で Fermi の式
    E = (8π/3)(g_e μ_B)(g_p μ_N)|ψ(0)|²⟨S·I⟩（g_e = 2）に一致するように決める：1s で ∫GF/r² → −π|ψ(0)|²、κ/(j(j+1)) = −4/3 なので
    C·(4π/3)|ψ(0)|² = (8π/3)·2μ_B·g_pμ_N|ψ(0)|² → C = 4 g_p μ_N μ_B = g_p α/(m_p/m_e)（μ_B = √α/2、μ_N = √α/(2 m_p/m_e)、ħ = c = m_e = 1）。
    単位：∫GF/r² dr は次元 1/長さ² なので λ̄_μ 単位から m_e 単位へ μ_red² を掛ける。さらに質量 μ の Dirac 方程式は小成分 F ≈ (G'+κG/r)/(2μ)
    を通じて磁気モーメント e/(2μ) を与えるが、電子の実際のモーメントは e/(2m_e) なので μ_red を掛けて戻す。合計 μ_red³ で、
    非相対論の閉じた式の (m_r/m_e)³ と同じ質量依存。相対論補正は 1s で 1/(γ(2γ−1)) ≈ 1 + (3/2)α²。"""
    I_GF = integ(r, G * Fr / r**2)                           # [1/λ̄_μ²]（質量 μ の単位系）
    return (X.G_P * ALPHA / MU_P) * (kappa / (j * (j + 1))) * I_GF * MU_RED**3 * 0.5 * (F * (F + 1) - 0.75 - j * (j + 1))


def build_levels(s=SIGN, n_max=N_MAX, verbose=True):
    """束縛準位の表。各準位：(n, ℓ, s_e, s_p) のラベル（固定一式と同じ）、エネルギーの内訳、波動関数。"""
    levels = []
    cache_schr = {}
    for n in range(1, n_max + 1):
        for kappa in kappa_list(n):
            ell, j = lj_of_kappa(kappa)
            sol = dirac_bound_state(n, kappa, s)
            if sol is None:
                continue
            eta, r, G, Fd = sol
            eps = 1.0 - eta
            if (n, ell) not in cache_schr:
                cache_schr[(n, ell)] = schrodinger_bound_state(n, ell, s)
            E_nr, rs, R, dR = cache_schr[(n, ell)]
            # 期待値（固定一式の単位へ換算：長さ λ̄_μ = LAM）
            r1 = integ(rs, R**2 * rs) / LAM
            r2 = integ(rs, R**2) / LAM**2
            r3 = integ(rs, R**2 / rs) / LAM**3 if ell > 0 else 0.0
            psi0sq = (R[0] ** 2 / (4.0 * np.pi)) / LAM**3 if ell == 0 else 0.0
            I_dR2r = integ(rs, dR**2 * rs) / LAM**3            # 次元 1/長さ³（⟨1/r³⟩ と同じ）
            E_nr_prog = E_nr * MU_RED
            m_tot = 1.0 + MU_P
            breit = (ALPHA**4 * MU_RED**3 / (2.0 * n**3 * MU_P**2) * (1.0 / (j + 0.5) - 1.0 / (ell + 0.5)) if ell > 0 else 0.0)
            coul = -MU_RED * eta - MU_RED**2 * eta**2 / (2.0 * m_tot) + breit + sum(X.recoil_level(n, ell))   # η = 1 − ε（桁落ちなし）
            for s_p in (-1, 1):
                F = j + 0.5 * s_p
                s_e = 1 if ell == 0 else (1 if j > ell else -1)
                grav = gravity_from_expectations(n, ell, j, F, r1, r2, r3, psi0sq, E_nr_prog, I_dR2r)
                hfs = hyperfine_dirac(kappa, j, F, r, G, Fd)
                levels.append(dict(label=(n, ell, s_e, s_p), n=n, ell=ell, j=j, F=F, kappa=kappa, eps=eps, eta=eta,
                                   coul=coul, grav=sum(grav.values()), grav_terms=grav, hfs=hfs,
                                   r=r, G=G, Fd=Fd, E_nr=E_nr_prog,
                                   expect=dict(r1=r1, r2=r2, r3=r3, psi0sq=psi0sq, I_dR2r=I_dR2r)))
            if verbose:
                print("  n=%d κ=%+d ℓ=%d j=%.1f  1−ε = %.15e" % (n, kappa, ell, j, eta), flush=True)
    return levels


# ---------- 遷移率（数値の Dirac 行列要素） ----------
def radial_me(la, lb, power):
    """∫ (G_a G_b + F_a F_b) r^power dr を a（Bohr 半径）の単位で。両準位の格子が違うので共通の対数格子に補間。"""
    x = np.linspace(max(np.log(la["r"][0]), np.log(lb["r"][0])), min(np.log(la["r"][-1]), np.log(lb["r"][-1])), 8000)
    r = np.exp(x)
    Ga = np.interp(x, np.log(la["r"]), la["G"]); Fa = np.interp(x, np.log(la["r"]), la["Fd"])
    Gb = np.interp(x, np.log(lb["r"]), lb["G"]); Fb = np.interp(x, np.log(lb["r"]), lb["Fd"])
    return simpson((Ga * Gb + Fa * Fb) * r**power * r, x=x) / A_BOHR**power


def edges_general(levels, gw_weight):
    ETOT = np.array([L["coul"] + L["grav"] + L["hfs"] for L in levels])
    mu = MU_RED
    a_si = X.HBAR_SI / (X.ME_SI * mu * X.C_SI * ALPHA)
    out = []
    cache = {}
    for i, La in enumerate(levels):
        for f, Lb in enumerate(levels):
            dE = ETOT[i] - ETOT[f]
            if dE <= 0.0:
                continue
            w = dE * ME_EV * X.E_SI / X.HBAR_SI
            x = w * a_si / X.C_SI
            n, l, j, F = La["n"], La["ell"], La["j"], La["F"]
            n2, l2, j2, F2 = Lb["n"], Lb["ell"], Lb["j"], Lb["F"]
            if abs(l2 - l) == 1:
                key = (La["kappa"], n, Lb["kappa"], n2, 1)
                if key not in cache:
                    cache[key] = radial_me(La, Lb, 1)
                A = (4.0 / 3.0) * ALPHA * w * x**2 * cache[key] ** 2 * max(l, l2) / (2 * l + 1) * X.recouple(l, j, l2, j2, F, F2, 1)
                if A > 0.0:
                    out.append((i, f, "e1", A))
            elif abs(l2 - l) in (0, 2) and l + l2 >= 2:
                key = (La["kappa"], n, Lb["kappa"], n2, 2)
                if key not in cache:
                    cache[key] = radial_me(La, Lb, 2)
                A = (1.0 / 15.0) * ALPHA * w * x**4 * cache[key] ** 2 * (2 * l2 + 1) * X.w3j000_sq(l, 2, l2) * X.recouple(l, j, l2, j2, F, F2, 2)
                if A > 0.0:
                    out.append((i, f, "e2", A))
                    out.append((i, f, "gw", gw_weight * A))
            if (n, l) == (n2, l2) and abs(j2 - j) <= 1 and abs(F2 - F) <= 1 and F + F2 >= 1:
                A = X.m1_rate(l, j, F, j2, F2, w)
                if A > 0.0:
                    out.append((i, f, "m1", A))
            if (n, l) == (2, 0) and (n2, l2) == (1, 0):
                if F2 == F:
                    out.append((i, f, "2ph", X.A_2GAMMA))
                out.append((i, f, "m1", X.A_M1_2S * (2 * F2 + 1) * (2 * j + 1) * X.wigner6j(j, F, 0.5, F2, j2, 1) ** 2))
    return out, ETOT


def run_master(levels, ed, ETOT, init_label=(5, 1, 1, -1), t_max=1.0e17, n_times=300):
    """固定一式 run() と同じ：τ = −R_TT⁻¹p₀、p_∞、流量、Radau の時間履歴。"""
    DIM = len(levels)
    R = np.zeros((DIM, DIM))
    for src, dst, kind, A in ed:
        R[dst, src] += A
        R[src, src] -= A
    idx = {L["label"]: i for i, L in enumerate(levels)}
    p0 = np.zeros(DIM)
    p0[idx[init_label]] = 1.0
    out_rate = -np.diag(R)
    trans = [i for i in range(DIM) if out_rate[i] > 0.0]
    absb = [i for i in range(DIM) if out_rate[i] == 0.0]
    tau = np.zeros(DIM)
    tau[trans] = np.linalg.solve(-R[np.ix_(trans, trans)], p0[trans])
    p_inf = p0.copy()
    p_inf[trans] = 0.0
    p_inf[absb] += np.dot(R[np.ix_(absb, trans)], tau[trans])
    flux_E = {}
    for src, dst, kind, A in ed:
        flux_E[kind] = flux_E.get(kind, 0.0) + A * tau[src] * (ETOT[src] - ETOT[dst])
    times = np.concatenate([[0.0], np.logspace(-12, np.log10(t_max), n_times)])
    sol = solve_ivp(lambda t, p: np.dot(R, p), (0.0, t_max), p0, method="Radau", jac=lambda t, p: R,
                    t_eval=times, rtol=1e-11, atol=1e-30)
    return dict(R=R, tau=tau, p_inf=p_inf, p0=p0, flux_E=flux_E, times=times, hist=sol.y.T, idx=idx, absorbing=absb)


def main():
    import time
    global SIGN
    if len(sys.argv) > 1:
        SIGN = float(sys.argv[1])          # 例：python3 two_body_general.py +1（斥力）、−1（引力、既定）
    tag = "attract" if SIGN < 0 else "repel"
    t0 = time.time()
    print("SIGN = %+d（%s）" % (SIGN, "引力" if SIGN < 0 else "斥力"))
    levels = build_levels(SIGN)
    print("束縛準位数 %d（%.1f s）" % (len(levels), time.time() - t0))
    if not levels:
        (HERE / ("audit_general_%s.txt" % tag)).write_text(
            "SIGN=%+d\nbound_levels=0\n"
            "同符号の電荷：クーロンが斥力で、重力（比 %.3e）では束縛できず、束縛準位が一つも無い。\n"
            "準位と率の模型（この一式）には表現する対象が無い。E–M が予言するのは散乱（Rutherford／Mott の反跳、同種粒子では双極子が消え四重極の制動放射、\n"
            "重力波との比 4Gm²/q²）で、それには過程を時間に沿って追う計算が要る。\n" % (SIGN, MU_P * X.MU**2 / 9.0), encoding="utf-8")
        print((HERE / ("audit_general_%s.txt" % tag)).read_text())
        return
    ed, ETOT = edges_general(levels, X.W_GW)
    res = run_master(levels, ed, ETOT)
    emitted = sum(res["flux_E"].values()) * ME_EV
    drop = (res["p0"] @ ETOT - res["p_inf"] @ ETOT) * ME_EV
    kinds = {}
    for _, _, k, _ in ed:
        kinds[k] = kinds.get(k, 0) + 1
    lines = ["SIGN=%+d" % SIGN, "bound_levels=%d" % len(levels),
             "final_state=%s" % [(levels[i]["label"], float(res["p_inf"][i])) for i in res["absorbing"] if res["p_inf"][i] > 0],
             "emitted_eV=%.9f" % emitted, "energy_drop_eV=%.9f" % drop, "energy_balance=%.3e" % (emitted - drop),
             "edges=%s" % kinds] + ["emitted_%s_eV=%.6e" % (k, v * ME_EV) for k, v in sorted(res["flux_E"].items())] + [
             "grav_eV=%.6e" % ((res["p_inf"] - res["p0"]) @ np.array([L["grav"] for L in levels]) * ME_EV),
             "sum_p_at_t_max=%r" % float(res["hist"][-1].sum()),
             "max_abs_p_t_max_minus_p_inf=%.1e" % float(np.max(np.abs(res["hist"][-1] - res["p_inf"]))),
             "elapsed_s=%.1f" % (time.time() - t0)]
    (HERE / ("audit_general_%s.txt" % tag)).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((HERE / ("audit_general_%s.txt" % tag)).read_text())


if __name__ == "__main__":
    main()
