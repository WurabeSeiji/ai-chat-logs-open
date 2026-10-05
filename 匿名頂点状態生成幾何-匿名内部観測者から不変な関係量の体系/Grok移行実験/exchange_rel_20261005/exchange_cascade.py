#!/usr/bin/env python3
"""水素型二体（電子 q = −3、陽子 q = +3）の多段放出。準位は Einstein–Maxwell（Dirac–Coulomb 厳密解＋二体反跳
＋遅延の高次 (Zα)⁵・(Zα)⁶、二体 1PN 重力、超微細）、遷移は Einstein の A 係数、時間発展はマスター方程式 dp/dt = R p。

基底は (n, ℓ, s_e, s_p)、n = 1..5、50 準位。s_e は j − ℓ の符号（ℓ = 0 は j = ½ のみ）、s_p は F − j の符号。
遷移は全ての下向き対（ΔE > 0）：E1（Δℓ = ±1）、E2 と重力波（Δℓ = 0, ±2、同じ四重極演算子、比 W_GW）、
M1（同じ nℓ：超微細 ΔF = ±1 と微細構造間 Δj = ±1。μ = −μ_B(L+2S) + g_pμ_N I）、
2s→1s の二光子（8.2206 s⁻¹、ΔF = 0）と M1（2.496e-6 s⁻¹）。n → n−1 の拘束はない。
外からの吸収は 0 で、増えたら失敗にする。
電荷と静止質量は初期値と比較し、違えば失敗にする。
重力辺は残す。物理の重み W_GW と、重み 1 の実験を両方記録する。
放出は各辺の流量 A·τ（τ = 滞在時間の積分、−R⁻¹p₀ で厳密）に ΔE を掛けて数える。時間刻みの誤差はない。
"""

from fractions import Fraction
from math import factorial, pi, sqrt
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.integrate import quad, solve_ivp
from scipy.special import genlaguerre

# フォント：元の Linux 固定パスは他の OS で失敗するので、あれば登録し、無ければ OS にある日本語フォントへ落とす。
if Path("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc").exists():
    font_manager.fontManager.addfont("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc")
_FONT_PREF = ["WenQuanYi Zen Hei", "Hiragino Sans", "Noto Sans CJK JP", "IPAexGothic", "Yu Gothic"]
_FONTS = {f.name for f in font_manager.fontManager.ttflist}
plt.rcParams["font.family"] = next((f for f in _FONT_PREF if f in _FONTS), "DejaVu Sans")
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).resolve().parent
ALPHA = 1.0 / 137.035999177
HBAR = 9.0 / ALPHA
MU_P = 1836.152673426
# 重力結合：SI の CODATA 2022（G、m_e、ħ、c）から α_G = G m_e²/(ħc) を作り、G = c = q₀ = 1 の単位での
# 電子質量を MU = 3√(α_G/α) と導く。重力／クーロンの比 g = G m_e m_p/e² = MU_P·MU²/9 は入力でなく結果。
# （元は MU を「比 = 4.0e-40」から逆算していた。正しい比は 4.4078e-40。）
G_SI, ME_SI, HBAR_SI, C_SI = 6.67430e-11, 9.1093837139e-31, 1.054571817e-34, 299792458.0
E_SI = 1.602176634e-19
ALPHA_G = G_SI * ME_SI**2 / (HBAR_SI * C_SI)
MU = float(3.0 * np.sqrt(ALPHA_G / ALPHA))
# 相対論では単位は静止質量 m_e c²。1s を Rydberg（無限質量・非相対論）13.6057 eV に合わせる釘付けは
# 相対論の補正を絶対値から消すのでやめ、1s の値は入力でなく結果にする。CODATA 2022。
ME_EV = 510998.95069
NEUTRON_EV = 0.782e6
Q_E, Q_P = -3.0, 3.0
# 陽子の g 因子（CODATA 2022）。点 Dirac 粒子なら 2。陽子の内部構造による唯一の外部入力。電子は Dirac の g = 2（QED なし）。
G_P = 5.5856946893
# 2s→1s：二光子 8.2206 s⁻¹（Goldman 1989）、M1 2.496e-6 s⁻¹（Johnson 1972）。E1 も E2 も禁制なので定数で与える。
A_2GAMMA, A_M1_2S = 8.2206, 2.496e-6
# 重力波／電気四重極の放出比。同じ四重極演算子 r_i r_j なので行列要素は共通で、比は古典の
# 四重極公式の比に等しい（Landau–Lifshitz §110 と §71：(G/45)/(1/180) = 4G）。質量四重極は μ r r、
# 電荷四重極は e (m_p − m_e)/M · r r（重心系）。よって W_GW = 4 G μ² M² / (e² (m_p − m_e)²) = 4 g m_e m_p/(m_p − m_e)²。
W_GW = 4.0 * (MU_P * MU**2 / 9.0) * MU_P / (MU_P - 1.0) ** 2


def basis():
    """(n, ℓ, s_e, s_p)。ℓ = 0 は j = ½ だけなので s_e = +1 のみ（元の基底にあった s_e = −1 の重複は除いた）。50 準位。"""
    return [(n, ell, se, sp)
            for n in range(1, 6)
            for ell in range(n)
            for se in ((1,) if ell == 0 else (-1, 1))
            for sp in (-1, 1)]


ST = basis()
IDX = {s: i for i, s in enumerate(ST)}
DIM = len(ST)
VALID = list(range(DIM))


def jf_of(ell, s_e, s_p):
    """s_e は j − ℓ の符号（j = ℓ ± ½、ℓ = 0 は j = ½）、s_p は F − j の符号（F = j ± ½）。"""
    j = max(ell + 0.5 * s_e, 0.5)
    F = j + 0.5 * s_p
    return j, F


def gravity_level(n, ell, j, F):
    """二体 1PN の重力準位（m_e c² 単位）。Coulomb 固有状態での一次摂動。

    ハミルトニアン（Barker–O'Connell 1975 の二体 1PN。電子・陽子は Dirac 粒子で、量子化の順序は
    静的計量 ds² = V²dt² − W²dx² の Dirac ハミルトニアン H = βmV + ½{F, α·p}（Obukhov 2001）の
    Foldy–Wouthuysen 変換で固定：H_FW ∋ mΦ + (3/4m){p², Φ} + (3ħ²/8m)∇²Φ + (3/2m)(∇Φ × p)·S）：
      Newton        −G m_e m_p / r
      1PN 自己項    −(3G/4)(m_p/m_e + m_e/m_p) {p², 1/r}          （Schwarzschild 1PN −(3GM/2m r)p² の両粒子分）
      重力 Darwin   +(3πGħ²/2)(m_p/m_e + m_e/m_p) δ³(r)           （FW の (3ħ²/8m)∇²Φ）
      1PN 交差項    −(G/2) p_i [(7δ_ij + n_i n_j)/r] p_j           （EIH の −(7/2)p₁·p₂ − ½(p₁·n)(p₂·n)）
      G² 項         +G² m_e m_p (m_e + m_p)/(2r²)
      重力–電磁交差 +G[−q_e q_p (m_e+m_p) + (m_p q_e² + m_e q_p²)/2]/r² = (27/2) G (m_e+m_p)/r²
                    （Khalil–Buonanno–Steinhoff–Vines 2018、EFT 再導出 arXiv:2205.11591 式(4.8) の静的項 H = −L。
                      第十思考実験 v7 の T8 = (3/2)κ_Cκ_G/ξ²。調和座標と等方座標は 1PN で一致するので FW の項と足せる）
      スピン軌道    (G/r³)[(2 + (3/2)m_p/m_e) L·S_e + (2 + (3/2)m_e/m_p) L·S_p]  （測地 3/2、Lense–Thirring 2。g 因子なし）
      スピン・スピン G[(3(S_e·n)(S_p·n) − S_e·S_p)/r³ + (8π/3) S_e·S_p δ³(r)]      （B_g = ∇×A_g の接触項）
    自己項と交差項を合わせると EIH の −(GMμ/2r)[(3+ν)p²/μ² + ν(n·p)²/μ²] になる。
    KN 質量四重極（第十思考実験 v7 の T21、T22）はトレースレスの 2 階テンソルで、m_j を分解しない準位では
    一次の寄与が恒等的に 0（Σ_m ⟨P₂⟩ = 0）なので、この基底では項が無い。
    期待値：⟨1/r⟩ = 1/(n²a)、⟨1/r²⟩ = 1/(n³(ℓ+½)a²)、⟨1/r³⟩ = 1/(n³ℓ(ℓ+½)(ℓ+1)a³)、|ψ(0)|² = δ_ℓ0/(πn³a³)、
    ⟨(1/r)p²⟩ = 2μ(E_n⟨1/r⟩ + e²⟨1/r²⟩)（固有状態で厳密）、
    ⟨p·(1/r)p⟩ = ⟨(1/r)p²⟩ − 2πħ²|ψ(0)|²、
    ⟨p_i(n_in_j/r)p_j⟩ = ħ²∫R'² r dr = ⟨(1/r)p²⟩ − ħ²ℓ(ℓ+1)⟨1/r³⟩ − 2πħ²|ψ(0)|²（gravity_levels_check.py で数値積分と照合）。
    単位：m_e = 1、ħ = 9/α、e² = 9、a = 9/(α²μ)、E_n = −μα²/(2n²)、G m_e m_p = 9g。
    """
    mu = MU_P / (1.0 + MU_P)
    g = MU_P * MU**2 / 9.0
    G = 9.0 * g / MU_P
    hb = HBAR
    a = 9.0 / (ALPHA**2 * mu)
    En = -mu * ALPHA**2 / (2.0 * n**2)
    r1 = 1.0 / (n**2 * a)
    r2 = 1.0 / (n**3 * (ell + 0.5) * a**2)
    r3 = 1.0 / (n**3 * ell * (ell + 0.5) * (ell + 1) * a**3) if ell > 0 else 0.0
    psi0sq = 1.0 / (np.pi * n**3 * a**3) if ell == 0 else 0.0
    kin = 2.0 * mu * (En * r1 + 9.0 * r2)
    p_r_p = kin - 2.0 * np.pi * hb**2 * psi0sq
    p_nn_p = kin - hb**2 * ell * (ell + 1) * r3 - 2.0 * np.pi * hb**2 * psi0sq
    newton = -9.0 * g * r1
    self_pn = -(3.0 * G / 4.0) * (MU_P + 1.0 / MU_P) * 2.0 * kin
    darwin = (3.0 * np.pi * G * hb**2 / 2.0) * (MU_P + 1.0 / MU_P) * psi0sq
    cross = -(G / 2.0) * (7.0 * p_r_p + p_nn_p)
    g2 = 9.0 * g * G * (1.0 + MU_P) / 2.0 * r2
    mixed = 13.5 * G * (1.0 + MU_P) * r2
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


# Bethe 対数 ln k0(n,ℓ)。CODATA 2010 Table V（Drake–Swainson 1990）の 10 桁、無いものは Jentschura–Mohr 2005
# （quant-ph/0504002、9 桁）。非相対論の水素固有状態だけで決まる量。
LN_K0 = {(1, 0): 2.984128556, (2, 0): 2.811769893, (3, 0): 2.767663612, (4, 0): 2.749811840, (5, 0): 2.74082373,
         (2, 1): -0.030016709, (3, 1): -0.0381902294, (4, 1): -0.0419548946, (5, 1): -0.0440346956,
         (3, 2): -0.00523214814, (4, 2): -0.00674093888, (5, 2): -0.00760075126,
         (4, 3): -0.00173366148, (5, 3): -0.00220216838,
         (5, 4): -0.000772098902}


def recoil_level(n, ell):
    """二体の遅延（反跳）の高次。Breit（1/c²）の先、電子–陽子間の遅延した光子交換が束縛状態で出す項。m_e c² 単位。
    CODATA 2010（Mohr–Taylor–Newell, RMP 84, 1527）式 (28)–(33)、Eides–Grotch–Shelyuto 2001 式 (144)(146)(147)。
      E_S（Salpeter 1952、(Zα)⁵ m²/M、m/M の全次数）
        = (m_r³/(m_e² m_N))(α⁵/(πn³)) m_e c² {(1/3)δ_ℓ0 ln α⁻² − (8/3) ln k0(n,ℓ) − (1/9)δ_ℓ0 − (7/3)a_n
                                               − (2/(m_N²−m_e²)) δ_ℓ0 [m_N² ln(m_e/m_r) − m_e² ln(m_N/m_r)]}
        a_n = −2[ln(2/n) + Σ_{i=1}^n 1/i + 1 − 1/(2n)] δ_ℓ0 + (1−δ_ℓ0)/(ℓ(ℓ+1)(2ℓ+1))
      E_R（(Zα)⁶ m²/M と (Zα)⁷ の先頭項）
        = (m_e/m_N)(α⁶/n³) m_e c² [D60 + D72 α ln²α⁻²]
        D60 = 4 ln 2 − 7/2（nS、Pachucki–Grotch 1995、Eides–Grotch 1997）
        D60 = [3 − ℓ(ℓ+1)/n²]·2/((4ℓ²−1)(2ℓ+3))（ℓ ≥ 1、Golosov ほか 1995、Jentschura–Pachucki 1996）
        D72 = −11/(60π)（nS、Pachucki–Karshenboim 1999、Melnikov–Yelkhovsky 1999。(Zα)⁷ はこの対数二乗項だけが既知）
    検証：1S の E_S = 2409.5 kHz、2S 341.3 kHz、E_R(1S) = −7.38 kHz（Eides 表 VIII。coulomb_levels_check.py）。
    放射補正（Lamb シフト、異常磁気能率）と核サイズは含まない。
    """
    mu = MU_P / (1.0 + MU_P)
    d0 = 1.0 if ell == 0 else 0.0
    Hn = sum(1.0 / i for i in range(1, n + 1))
    a_n = -2.0 * (np.log(2.0 / n) + Hn + 1.0 - 1.0 / (2.0 * n)) * d0 + (0.0 if ell == 0 else 1.0 / (ell * (ell + 1) * (2 * ell + 1)))
    masslog = (2.0 / (MU_P**2 - 1.0)) * d0 * (MU_P**2 * np.log(1.0 / mu) - np.log(MU_P / mu))
    E_S = (mu**3 / MU_P) * (ALPHA**5 / (np.pi * n**3)) * (
        (1.0 / 3.0) * d0 * np.log(ALPHA**-2) - (8.0 / 3.0) * LN_K0[(n, ell)] - d0 / 9.0 - (7.0 / 3.0) * a_n - masslog)
    if ell == 0:
        D60, D72 = 4.0 * np.log(2.0) - 3.5, -11.0 / (60.0 * np.pi)
    else:
        D60, D72 = (3.0 - ell * (ell + 1) / n**2) * 2.0 / ((4 * ell**2 - 1) * (2 * ell + 3)), 0.0
    E_R = (1.0 / MU_P) * (ALPHA**6 / n**3) * (D60 + D72 * ALPHA * np.log(ALPHA**-2) ** 2)
    return E_S, E_R


def hyperfine_level(n, ell, j, F):
    """超微細（Maxwell の磁気双極子相互作用：Fermi 接触＋スピン双極子＋軌道。μ の Bohr 半径、g_e = 2）。m_e c² 単位。
    E = K [F(F+1) − j(j+1) − ¾] / [4 n³ (ℓ+½) j(j+1)]、K = g_p α⁴ μ³ / (m_p/m_e)。
    1s の分裂 (4/3)K = 1418.8 MHz（実測 1420.4。差は電子の異常磁気能率 g_e/2 = 1.00116 と QED）、2s 177.4、2p₁/₂ 59.1、2p₃/₂ 23.6 MHz。
    元の se = s_e α³/(9n³) は電子のスピン軌道で Dirac 厳密解に含まれるので削除（二重計上）。sp はこの式に置換。
    """
    mu = MU_P / (1.0 + MU_P)
    K = G_P * ALPHA**4 * mu**3 / MU_P
    return K * (F * (F + 1) - j * (j + 1) - 0.75) / (4.0 * n**3 * (ell + 0.5) * j * (j + 1))


def parts():
    mu_red = MU_P / (1.0 + MU_P)
    orb_c, orb_g, hfs = [], [], []
    m_tot = 1.0 + MU_P
    for n, ell, s_e, s_p in ST:
        # クーロン：Bohr（非相対論）→ 二体の相対論的クーロン準位（Barker–Glover 1955、CODATA の E_M と同じ式）。
        # 結合は q_e q_p/ħ = 9/ħ = α。
        #   第 1 項：電子の Dirac–Coulomb 厳密解 f(n,j)（Sommerfeld の微細構造を全次数で含む）
        #   第 2 項：陽子の反跳。二体の相対論的運動学 W² = m² + M² + 2mM·f の展開 −μ²(f−1)²/(2(m+M))
        #   第 3 項：Breit 相互作用の軌道–軌道・スピン–他軌道（磁気・遅延）α⁴μ³/(2n³M²)[1/(j+½) − 1/(ℓ+½)]、ℓ ≥ 1
        # 展開：−μα²/(2n²) − μα⁴/(2n³)[1/(j+½) − 3/(4n) + μ/(4n(m+M))] + α⁴μ³/(2n³M²)[1/(j+½) − 1/(ℓ+½)](1−δ_ℓ0)
        # これは CODATA の E_M（式 (26)、(1−δ_ℓ0)/(κ(2ℓ+1)) は上の [1/(j+½) − 1/(ℓ+½)] と同じ）。
        #   第 4 項：遅延の高次（recoil_level）。Salpeter の (Zα)⁵ と (Zα)⁶・(Zα)⁷ の既知項。
        j, F = jf_of(ell, s_e, s_p)
        kappa = j + 0.5
        delta = kappa - np.sqrt(kappa**2 - ALPHA**2)
        f_dirac = (1.0 + (ALPHA / (n - delta))**2) ** -0.5
        breit = (ALPHA**4 * mu_red**3 / (2.0 * n**3 * MU_P**2)
                 * (1.0 / (j + 0.5) - 1.0 / (ell + 0.5)) if ell > 0 else 0.0)
        coul = mu_red * (f_dirac - 1.0) - mu_red**2 * (f_dirac - 1.0)**2 / (2.0 * m_tot) + breit + sum(recoil_level(n, ell))
        orb_c.append(coul)
        # 重力：二体 1PN（gravity_level）。クーロンと桁が 1e-40 違うので ETOT の中では倍精度で消えるが、
        # 配列 ORB_G として別に保持し、重力の pool はこの配列で数える。
        orb_g.append(sum(gravity_level(n, ell, j, F).values()))
        hfs.append(hyperfine_level(n, ell, j, F))
    return np.array(orb_c), np.array(orb_g), np.array(hfs)


ORB_C, ORB_G, HFS = parts()
ORB = ORB_C + ORB_G
ETOT = ORB + HFS
EV = ME_EV


# ---------- 角運動量の代数（Racah の式。厳密な有理数） ----------
def _tri(a, b, c):
    return (a + b + c).denominator == 1 and abs(a - b) <= c <= a + b


def _Delta(a, b, c):
    return Fraction(factorial(int(a + b - c)) * factorial(int(a - b + c)) * factorial(int(-a + b + c)),
                    factorial(int(a + b + c + 1)))


def wigner6j(a, b, c, d, e, f):
    """{a b c; d e f}（符号つき）。"""
    a, b, c, d, e, f = (Fraction(x) for x in (a, b, c, d, e, f))
    for tri in ((a, b, c), (a, e, f), (d, b, f), (d, e, c)):
        if not _tri(*tri):
            return 0.0
    pref = _Delta(a, b, c) * _Delta(a, e, f) * _Delta(d, b, f) * _Delta(d, e, c)
    kmin = int(max(a + b + c, a + e + f, d + b + f, d + e + c))
    kmax = int(min(a + b + d + e, a + c + d + f, b + c + e + f))
    s = Fraction(0)
    for k in range(kmin, kmax + 1):
        s += (-1) ** k * Fraction(factorial(k + 1),
                                  factorial(k - int(a + b + c)) * factorial(k - int(a + e + f))
                                  * factorial(k - int(d + b + f)) * factorial(k - int(d + e + c))
                                  * factorial(int(a + b + d + e) - k) * factorial(int(a + c + d + f) - k)
                                  * factorial(int(b + c + e + f) - k))
    return float(np.sign(float(s)) * sqrt(float(pref * s * s)))


def w3j000_sq(l1, l2, l3):
    """(l1 l2 l3; 0 0 0)²。"""
    J = l1 + l2 + l3
    if J % 2 or not (abs(l1 - l2) <= l3 <= l1 + l2):
        return 0.0
    g = J // 2
    D = Fraction(factorial(J - 2 * l1) * factorial(J - 2 * l2) * factorial(J - 2 * l3), factorial(J + 1))
    return float(D * Fraction(factorial(g), factorial(g - l1) * factorial(g - l2) * factorial(g - l3)) ** 2)


def recouple(l, j, l2, j2, F, F2, k):
    """スピン ½（電子）と核スピン ½ を傍観者にした k 階テンソル遷移の角度因子。
    (2j'+1)(2ℓ+1){ℓ j ½; j' ℓ' k}² × (2F'+1)(2j+1){j F ½; F' j' k}²。j'、F' について和をとると 1。"""
    return ((2 * j2 + 1) * (2 * l + 1) * wigner6j(l, j, 0.5, j2, l2, k) ** 2
            * (2 * F2 + 1) * (2 * j + 1) * wigner6j(j, F, 0.5, F2, j2, k) ** 2)


# ---------- 水素の動径積分（μ の Bohr 半径 a = 1 の単位） ----------
def radial_fn(n, l):
    N = sqrt((2.0 / n) ** 3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    Lag = genlaguerre(n - l - 1, 2 * l + 1)
    return lambda r: N * np.exp(-r / n) * (2.0 * r / n) ** l * Lag(2.0 * r / n)


_RAD = {}


def radial(n, l, n2, l2, power):
    """∫ R_{n'ℓ'} r^power R_{nℓ} r² dr（a 単位）。"""
    key = (n, l, n2, l2, power)
    if key not in _RAD:
        R, R2 = radial_fn(n, l), radial_fn(n2, l2)
        _RAD[key] = quad(lambda r: R(r) * R2(r) * r ** (power + 2), 0, 60.0 * max(n, n2), limit=800, epsabs=0.0, epsrel=1e-12)[0]
    return _RAD[key]


# ---------- Einstein の A 係数 ----------
def lande_g(ell, j):
    return 1.0 + (j * (j + 1) + 0.75 - ell * (ell + 1)) / (2.0 * j * (j + 1))


def m1_reduced(ell, j, F, j2, F2):
    """同じ nℓ の二準位間の M1 の換算行列要素（μ_B 単位）。μ = −μ_B(L + 2S)/ħ + g_p μ_N I/ħ = −μ_B(J + S)/ħ + g_p μ_N I/ħ。
    結合 ((ℓ s) j I) F。Edmonds 7.1.7／7.1.8：
      ⟨(j' I)F'||J||(j I)F⟩ = δ_jj' (−1)^{j+I+F+1} √((2F+1)(2F'+1)) {j F' I; F j 1} √(j(j+1)(2j+1))
      ⟨((ℓs)j' I)F'||S||((ℓs)j I)F⟩ = (−1)^{j'+I+F+1} √((2F+1)(2F'+1)) {j' F' I; F j 1} ⟨(ℓs)j'||S||(ℓs)j⟩、
         ⟨(ℓs)j'||S||(ℓs)j⟩ = (−1)^{ℓ+s+j'+1} √((2j+1)(2j'+1)) {s j' ℓ; j s 1} √(s(s+1)(2s+1))
      ⟨(j I)F'||I||(j I)F⟩ = δ_jj' (−1)^{j+I+F'+1} √((2F+1)(2F'+1)) {I F' j; F I 1} √(I(I+1)(2I+1))
    j' = j では ⟨J⟩ + ⟨S⟩ = g_j⟨J⟩（射影定理。rates_check.py で確認）。j' ≠ j は微細構造間の M1（⟨S⟩ だけ）。"""
    I, s = 0.5, 0.5
    pref = sqrt((2 * F + 1) * (2 * F2 + 1))
    redJ = 0.0 if j2 != j else (-1) ** int(j + I + F + 1) * pref * wigner6j(j, F2, I, F, j, 1) * sqrt(j * (j + 1) * (2 * j + 1))
    redS_j = (-1) ** int(ell + s + j2 + 1) * sqrt((2 * j + 1) * (2 * j2 + 1)) * wigner6j(s, j2, ell, j, s, 1) * sqrt(s * (s + 1) * (2 * s + 1))
    redS = (-1) ** int(j2 + I + F + 1) * pref * wigner6j(j2, F2, I, F, j, 1) * redS_j
    redI = 0.0 if j2 != j else (-1) ** int(j + I + F2 + 1) * pref * wigner6j(I, F2, j, F, I, 1) * sqrt(I * (I + 1) * (2 * I + 1))
    return -(redJ + redS) + (G_P / MU_P) * redI


def m1_rate(ell, j, F, j2, F2, omega):
    """同じ nℓ の二準位間の M1。A = (4ω³/3ħc³)|⟨F'||μ||F⟩|²/(2F+1) = (α/3)ω³(ħ/m_ec²)²|M|²/(2F+1)、M は μ_B 単位。
    超微細（j' = j、ΔF = ±1）：1s F=1→0（21 cm）2.868e-15 s⁻¹（実測 2.884e-15、差は g_e/2 の 2 乗と ω³）。
    微細構造間（j' = j ± 1）：2p₃/₂→2p₁/₂ は 1e-12 s⁻¹ 台。"""
    return (ALPHA / 3.0) * omega**3 * (HBAR_SI / (ME_SI * C_SI**2)) ** 2 * m1_reduced(ell, j, F, j2, F2) ** 2 / (2 * F + 1)


def edges(gw_weight):
    """全ての下向き対の A 係数 [s⁻¹]。(src, dst, kind, A)。"""
    mu = MU_P / (1.0 + MU_P)
    a_si = HBAR_SI / (ME_SI * mu * C_SI * ALPHA)
    out = []
    for i in VALID:
        n, l, se, sp = ST[i]
        j, F = jf_of(l, se, sp)
        for f in VALID:
            n2, l2, se2, sp2 = ST[f]
            j2, F2 = jf_of(l2, se2, sp2)
            dE = ETOT[i] - ETOT[f]
            if dE <= 0.0:
                continue
            w = dE * ME_EV * E_SI / HBAR_SI
            x = w * a_si / C_SI
            if abs(l2 - l) == 1:
                A = (4.0 / 3.0) * ALPHA * w * x**2 * radial(n, l, n2, l2, 1) ** 2 * max(l, l2) / (2 * l + 1) * recouple(l, j, l2, j2, F, F2, 1)
                if A > 0.0:
                    out.append((i, f, "e1", A))
            elif abs(l2 - l) in (0, 2) and l + l2 >= 2:
                A = (1.0 / 15.0) * ALPHA * w * x**4 * radial(n, l, n2, l2, 2) ** 2 * (2 * l2 + 1) * w3j000_sq(l, 2, l2) * recouple(l, j, l2, j2, F, F2, 2)
                if A > 0.0:
                    out.append((i, f, "e2", A))
                    out.append((i, f, "gw", gw_weight * A))
            if (n, l) == (n2, l2) and abs(j2 - j) <= 1 and abs(F2 - F) <= 1 and F + F2 >= 1:
                A = m1_rate(l, j, F, j2, F2, w)
                if A > 0.0:
                    out.append((i, f, "m1", A))
            if (n, l) == (2, 0) and (n2, l2) == (1, 0):
                if F2 == F:
                    out.append((i, f, "2ph", A_2GAMMA))
                out.append((i, f, "m1", A_M1_2S * (2 * F2 + 1) * (2 * j + 1) * wigner6j(j, F, 0.5, F2, j2, 1) ** 2))
    return out


def generator(ed):
    R = np.zeros((DIM, DIM))
    for src, dst, kind, A in ed:
        R[dst, src] += A
        R[src, src] -= A
    return R


def run(gw_weight, t_max=1.0e17, n_times=300):
    """マスター方程式。滞在時間 τ = −R_TT⁻¹ p₀（厳密）から各辺の流量と終状態。時間発展（図と時刻ごとの pool）は
    硬い系なので陰的 Radau（Jacobian = R）で解く。率は 6e8 〜 1e-15 s⁻¹、時刻は 1e-12 〜 1e17 s にまたがる。"""
    ed = edges(gw_weight)
    R = generator(ed)
    p0 = np.zeros(DIM)
    p0[IDX[(5, 1, 1, -1)]] = 1.0
    charge = Q_E + Q_P
    mass_e, mass_p = MU, MU_P * MU
    e_in = 0.0
    if abs((Q_E + Q_P) - charge) > 0.0 or e_in != 0.0:
        raise RuntimeError("電荷または外からの吸収が動いた")
    if abs(mass_e - MU) > 0.0 or abs(mass_p - MU_P * MU) > 0.0:
        raise RuntimeError("静止質量が動いた")
    out_rate = -np.diag(R)
    trans = [i for i in VALID if out_rate[i] > 0.0]
    absb = [i for i in VALID if out_rate[i] == 0.0]
    tau = np.zeros(DIM)
    tau[trans] = np.linalg.solve(-R[np.ix_(trans, trans)], p0[trans])
    p_inf = p0.copy()
    p_inf[trans] = 0.0
    p_inf[absb] += np.dot(R[np.ix_(absb, trans)], tau[trans])
    flux_E = {}
    for src, dst, kind, A in ed:
        flux_E[kind] = flux_E.get(kind, 0.0) + A * tau[src] * (ETOT[src] - ETOT[dst])
    times = np.concatenate([[0.0], np.logspace(-12, np.log10(t_max), n_times)])
    # np.dot を使う：この環境（numpy 2.0.2 + Accelerate）では 50×50 以上の `@` が有限な乱数行列でも
    # "divide by zero in matmul" の偽の警告を出す（np.dot は出さず、値は一致。README 参照）。
    sol = solve_ivp(lambda t, p: np.dot(R, p), (0.0, t_max), p0, method="Radau", jac=lambda t, p: R,
                    t_eval=times, rtol=1e-11, atol=1e-30)
    if not sol.success:
        raise RuntimeError("時間発展が解けない：" + sol.message)
    hist = sol.y.T
    pools = []
    for k in range(1, len(times)):
        p, nxt = hist[k - 1], hist[k]
        d_hfs = nxt @ HFS - p @ HFS
        d_coul = nxt @ ORB_C - p @ ORB_C
        d_grav = nxt @ ORB_G - p @ ORB_G
        d_out = -(d_hfs + d_coul + d_grav) + e_in
        pools.append((d_hfs, d_coul, d_grav, d_out))
    return dict(times=times, hist=hist, pools=np.array(pools), ed=ed, R=R, tau=tau, p_inf=p_inf, p0=p0,
                flux_E=flux_E, absorbing=absb)


def main():
    res = run(W_GW)
    res_g = run(1.0)
    times, hist, pools = res["times"], res["hist"], res["pools"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0), dpi=140)
    ax = axes[0]
    tplot = np.where(times > 0, times, 1e-13)
    for n in range(1, 6):
        ax.plot(tplot, sum(hist[:, i] for i, s in enumerate(ST) if s[0] == n), lw=1.6, label=f"n={n}")
    ax.set_xscale("log")
    ax.set_xlabel("時刻 [s]")
    ax.set_ylabel("占有")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.set_title("マスター方程式の多段")
    ax = axes[1]
    arr = pools * EV
    tmid = tplot[1:]
    ax.plot(tmid, np.cumsum(arr[:, 3]), lw=1.6, label="管への放出（累積）")
    ax.plot(tmid, np.cumsum(arr[:, 1]), lw=1.6, label="クーロン（累積）")
    ax.plot(tmid, np.cumsum(arr[:, 2]) * 1e38, lw=1.6, label="重力（累積）×1e38")
    ax.set_xscale("log")
    ax.set_xlabel("時刻 [s]")
    ax.set_ylabel("eV")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("三つの変化と放出")
    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "cascade_result.png")
    p0, p_inf = res["p0"], res["p_inf"]
    emitted = float(sum(res["flux_E"].values()) * EV)
    e_drop = float((p0 @ ETOT - p_inf @ ETOT) * EV)
    grav = float((p_inf @ ORB_G - p0 @ ORB_G) * EV)
    final = [float(sum(p_inf[i] for i, s in enumerate(ST) if s[0] == n)) for n in range(1, 6)]
    kinds = {}
    for _, _, k, _ in res["ed"]:
        kinds[k] = kinds.get(k, 0) + 1
    lines = [
        f"final_n={final}",
        f"final_state={[(ST[i], float(p_inf[i])) for i in res['absorbing'] if p_inf[i] > 0]}",
        f"emitted_eV={emitted:.9f}",
        f"energy_drop_eV={e_drop:.9f}",
        f"energy_balance={emitted - e_drop:.3e}",
        f"neutron_open={emitted >= NEUTRON_EV}",
        f"charge={Q_E + Q_P}",
        f"mass_e={MU}", f"mass_p={MU_P * MU}",
        "e_in=0",
        f"t_max_s={times[-1]:.1e}",
        f"sum_p_at_t_max={float(hist[-1].sum())!r}",
        f"max_abs_p_t_max_minus_p_inf={float(np.max(np.abs(hist[-1] - p_inf))):.1e}",
        f"min_p_t={float(hist.min()):.1e}",
        f"edges={kinds}",
    ] + [f"emitted_{k}_eV={v * EV:.6e}" for k, v in sorted(res["flux_E"].items())] + [
        f"gw_out_phys_eV={res['flux_E'].get('gw', 0.0) * EV:.6e}",
        f"gw_out_weight1_eV={res_g['flux_E'].get('gw', 0.0) * EV:.6e}",
        f"grav_eV={grav:.6e}",
        f"pool_sum_vs_out={float(np.max(np.abs(pools.sum(axis=1))))}",
    ]
    (OUT / "audit.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((OUT / "audit.txt").read_text())


if __name__ == "__main__":
    main()
