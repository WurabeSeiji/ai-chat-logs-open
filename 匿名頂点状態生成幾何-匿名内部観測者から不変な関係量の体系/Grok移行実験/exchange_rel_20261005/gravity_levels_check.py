#!/usr/bin/env python3
"""重力準位（gravity_level）と重力波比（W_GW）の検算。プログラム本体は触らず import するだけ。

1. 結合定数：g = G m_e m_p/e² を SI の CODATA から独立に計算し、プログラムの MU_P·MU²/9 と比較。
2. 期待値の恒等式：水素の動径関数 R_nℓ（Laguerre）で数値積分し、プログラムが使う閉じた式と比較。
     ⟨1/r⟩、⟨1/r²⟩、⟨1/r³⟩、|ψ(0)|²、
     ∫R'² r dr = (2μ/ħ²)[E_n⟨1/r⟩ + e²⟨1/r²⟩] − ℓ(ℓ+1)⟨1/r³⟩ − 2δ_ℓ0/(n³a³)      （⟨p_i n_i n_j/r p_j⟩/ħ²）
     ∫[R'² r + ℓ(ℓ+1)R²/r] dr = (2μ/ħ²)[E_n⟨1/r⟩ + e²⟨1/r²⟩] − 2δ_ℓ0/(n³a³)        （⟨p·(1/r)p⟩/ħ²）
3. Newton 項 = 2g × Bohr(n)。EIH の形 −(GMμ/2)[(3+ν)p²/μ² + ν(n·p)²/μ²] と自己項＋交差項の一致（ℓ ≥ 1）。
4. 重力波：W_GW を SI から独立に計算。3d→1s の E2 率 A = (ω⁵/90ħc⁵)(6/5)e²|⟨r²⟩₁₃|² を計算し、
   Γ_GW = W_GW·A と寿命を出す。
5. 準位ごとの重力項の内訳（eV）。
"""
import sys
from math import factorial, pi, sqrt
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import genlaguerre
from matplotlib import font_manager

font_manager.fontManager.addfont = lambda path: None   # 17 行目の Linux フォント登録を無効化（本体は変更しない）
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exchange_cascade as X  # noqa: E402

L = []
A = X.ALPHA
mu = X.MU_P / (1.0 + X.MU_P)
g_prog = X.MU_P * X.MU**2 / 9.0

# 1. 結合定数（SI、CODATA 2022）
G_SI, ME, MP, HB, C = 6.67430e-11, 9.1093837139e-31, 1.67262192595e-27, 1.054571817e-34, 299792458.0
E_SI, EPS0 = 1.602176634e-19, 8.8541878188e-12
K = 1.0 / (4.0 * pi * EPS0)
g_si = G_SI * ME * MP / (K * E_SI**2)
L.append("1. 結合定数")
L.append("   α_G = G m_e²/(ħc) = %.6e、MU = 3√(α_G/α) = %.6e" % (X.ALPHA_G, X.MU))
L.append("   g = G m_e m_p/e²：プログラム %.6e、SI から独立に %.6e、比 %.8f（元の入力は 4.0e-40）" % (g_prog, g_si, g_prog / g_si))

# 2. 動径関数と期待値の恒等式（a = ħ = μ = e = 1）
def R_and_dR(n, l):
    N = sqrt((2.0 / n) ** 3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    Lag = genlaguerre(n - l - 1, 2 * l + 1)
    dLag = Lag.deriv()
    def R(r):
        rho = 2.0 * r / n
        return N * np.exp(-rho / 2) * rho**l * Lag(rho)
    def dR(r):
        rho = 2.0 * r / n
        return N * np.exp(-rho / 2) * (2.0 / n) * (-0.5 * rho**l * Lag(rho) + (l * rho**(l - 1) * Lag(rho) if l > 0 else 0.0) + rho**l * dLag(rho))
    return R, dR

L.append("")
L.append("2. 期待値の恒等式（数値積分との差。a = ħ = μ = e = 1）")
worst = {"r1": 0, "r2": 0, "r3": 0, "psi0": 0, "nn": 0, "pp": 0, "norm": 0}
for n in range(1, 6):
    for l in range(n):
        R, dR = R_and_dR(n, l)
        I = lambda f: quad(f, 0, np.inf, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
        norm = I(lambda r: R(r) ** 2 * r**2)
        r1, r2 = 1.0 / n**2, 1.0 / (n**3 * (l + 0.5))
        r3 = 1.0 / (n**3 * l * (l + 0.5) * (l + 1)) if l > 0 else 0.0
        psi0 = (1.0 / (pi * n**3)) if l == 0 else 0.0
        E = -0.5 / n**2
        worst["norm"] = max(worst["norm"], abs(norm - 1))
        worst["r1"] = max(worst["r1"], abs(I(lambda r: R(r) ** 2 * r) - r1))
        worst["r2"] = max(worst["r2"], abs(I(lambda r: R(r) ** 2) - r2))
        if l > 0:
            worst["r3"] = max(worst["r3"], abs(I(lambda r: R(r) ** 2 / r) - r3))
        else:
            worst["psi0"] = max(worst["psi0"], abs(R(1e-12) ** 2 / (4 * pi) - psi0))
        nn_num = I(lambda r: dR(r) ** 2 * r)
        nn_cf = 2.0 * (E * r1 + r2) - l * (l + 1) * r3 - 2.0 * pi * psi0
        pp_num = nn_num + (l * (l + 1) * I(lambda r: R(r) ** 2 / r) if l > 0 else 0.0)
        pp_cf = 2.0 * (E * r1 + r2) - 2.0 * pi * psi0
        worst["nn"] = max(worst["nn"], abs(nn_num - nn_cf))
        worst["pp"] = max(worst["pp"], abs(pp_num - pp_cf))
L.append("   n = 1..5 全 (n,ℓ)：規格化 %.1e、⟨1/r⟩ %.1e、⟨1/r²⟩ %.1e、⟨1/r³⟩ %.1e、|ψ(0)|² %.1e、"
         "⟨p nn/r p⟩ %.1e、⟨p (1/r) p⟩ %.1e" % tuple(worst[k] for k in ("norm", "r1", "r2", "r3", "psi0", "nn", "pp")))

# 3. Newton = 2g·Bohr、EIH の形との一致
L.append("")
L.append("3. Newton 項と EIH の形")
hb = X.HBAR
G = 9.0 * g_prog / X.MU_P
nu = mu / (1.0 + X.MU_P)
M = 1.0 + X.MU_P
w_newton = 0.0
w_eih = 0.0
for n in range(1, 6):
    for l in range(n):
        for s_e in (-1, 1):
            j, F = X.jf_of(l, s_e, -1)
            t = X.gravity_level(n, l, j, F)
            bohr = -mu * A**2 / (2.0 * n**2)
            w_newton = max(w_newton, abs(t["newton"] - 2.0 * g_prog * bohr) / abs(t["newton"]))
            if l > 0:
                a = 9.0 / (A**2 * mu)
                En = bohr
                r1 = 1.0 / (n**2 * a); r2 = 1.0 / (n**3 * (l + 0.5) * a**2); r3 = 1.0 / (n**3 * l * (l + 0.5) * (l + 1) * a**3)
                kin = 2.0 * mu * (En * r1 + 9.0 * r2)
                pnn = kin - hb**2 * l * (l + 1) * r3
                eih = -(G * M * mu / 2.0) * ((3.0 + nu) * kin / mu**2 + nu * pnn / mu**2)
                w_eih = max(w_eih, abs((t["self_pn"] + t["cross"]) - eih) / abs(eih))
L.append("   Newton = 2g·E_Bohr(n)：最大相対差 %.1e" % w_newton)
L.append("   自己項＋交差項 = −(GMμ/2)[(3+ν)⟨p²/r⟩/μ² + ν⟨(n·p)²/r⟩/μ²]（ℓ ≥ 1）：最大相対差 %.1e" % w_eih)
L.append("   スピン軌道の係数：電子 2 + (3/2)m_p/m_e = %.3f（測地 3/2）、陽子 2 + (3/2)m_e/m_p = %.6f（Lense–Thirring 2）"
         % (2 + 1.5 * X.MU_P, 2 + 1.5 / X.MU_P))
# 重力–電磁交差項：v7 の T8 = (3/2)κ_Cκ_G/ξ²（κ_C = 9/(μa_a)、κ_G = GM/a_a、ξ = r/a_a、μ で割った形）→ μ を掛けて (27/2)GM/r²
t1s = X.gravity_level(1, 0, 0.5, 1.0)
a1 = 9.0 / (A**2 * mu)
kC, kG = 9.0 / (mu * a1), G * M / a1          # v7 §0.4 の κ_C、κ_G（a_a を a に置く）
r2_1s = 1.0 / (1.0 * 0.5 * a1**2)              # ⟨1/r²⟩_1s
L.append("   重力–電磁交差 (27/2)G(m_e+m_p)⟨1/r²⟩：1s で %.5e eV。v7 の T8 = (3/2)κ_Cκ_G⟨1/ξ²⟩ を μ 倍した値 %.5e eV（同一）"
         % (t1s["mixed"] * X.EV, 1.5 * kC * kG * mu * a1**2 * r2_1s * X.EV))

# 4. 重力波
L.append("")
L.append("4. 重力波（四重極）")
w_si = 4.0 * G_SI * ME**2 * MP**2 / (K * E_SI**2 * (MP - ME) ** 2)
L.append("   W_GW = 4Gμ²M²/(e²(m_p−m_e)²)：プログラム %.6e、SI から独立に %.6e（元の重みは Newton の比 4.0e-40 で、未使用だった）"
         % (X.W_GW, w_si))
mu_si = ME * MP / (ME + MP)
a_si = HB / (mu_si * C * A)
dE = (X.ORB_C[X.IDX[(1, 0, 1, -1)]] - X.ORB_C[X.IDX[(3, 2, 1, -1)]]) * X.EV * E_SI   # J（符号は大きさだけ使う）
omega = abs(dE) / HB
R1, _ = R_and_dR(1, 0); R3, _ = R_and_dR(3, 2)
r2_13 = quad(lambda r: R1(r) * R3(r) * r**4, 0, np.inf, limit=400)[0]        # a 単位
r2_13_cf = 8.0 * factorial(6) / (81.0 * sqrt(30.0)) * (3.0 / 4.0) ** 7
A_e2 = (6.0 / 5.0) * A * omega**5 * (r2_13 * a_si**2) ** 2 / (90.0 * C**4)
L.append("   3d→1s：⟨r²⟩₁₃ = %.5f a（閉じた式 %.5f）、ω = %.4e s⁻¹、A_E2 = (ω⁵/90ħc⁵)(6/5)e²|⟨r²⟩|² = %.1f s⁻¹"
         % (r2_13, r2_13_cf, omega, A_e2))
L.append("   Γ_GW(3d→1s) = W_GW·A_E2 = %.3e s⁻¹、寿命 %.2e 年" % (X.W_GW * A_e2, 1.0 / (X.W_GW * A_e2) / 3.15576e7))
ed = X.edges(X.W_GW)
kinds = {}
for a_, b_, k_, w_ in ed:
    kinds[k_] = kinds.get(k_, 0) + 1
L.append("   辺の数：%s" % ", ".join("%s %d" % kv for kv in sorted(kinds.items())))

# 5. 内訳
L.append("")
L.append("5. 準位ごとの重力項（eV、EV = m_e c²）。s_e → j、s_p → F")
L.append("   (n,ℓ,j,F)        Newton         1PN自己        Darwin         1PN交差        G²            重力–電磁交差  SO_e           SO_p           SS             合計")
for n in range(1, 4):
    for l in range(n):
        for s_e in (-1, 1):
            for s_p in (-1, 1):
                if l == 0 and s_e == -1:
                    continue
                j, F = X.jf_of(l, s_e, s_p)
                t = X.gravity_level(n, l, j, F)
                tot = X.ORB_G[X.IDX[(n, l, s_e, s_p)]]
                L.append("   (%d,%d,%.1f,%d)  " % (n, l, j, F) + " ".join("%+14.5e" % (t[k] * X.EV) for k in ("newton", "self_pn", "darwin", "cross", "g2", "mixed", "so_e", "so_p", "ss")) + " %+14.5e" % (tot * X.EV))
L.append("   1s の重力束縛 %.6e eV は Coulomb 13.5985 eV の %.3e 倍（= 2g）。1PN 項は Newton の α² 程度、スピン項は 1PN と同程度。" % (
    -X.gravity_level(1, 0, 0.5, 1.0)["newton"] * X.EV, X.gravity_level(1, 0, 0.5, 1.0)["newton"] / X.ORB_C[X.IDX[(1, 0, 1, -1)]]))

out = "\n".join(L) + "\n"
(HERE / "gravity_levels_check.txt").write_text(out, encoding="utf-8")
print(out)
