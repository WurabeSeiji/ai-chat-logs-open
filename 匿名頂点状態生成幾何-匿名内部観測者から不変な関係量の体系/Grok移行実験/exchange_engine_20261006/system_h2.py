#!/usr/bin/env python3
"""系の部品：水素分子 H₂（基底電子状態 X¹Σg⁺ の振動回転準位 (v, J)）。
エンジン engine.py に渡すもの：states（準位：E_eV、J、M_kg）と channels（src→dst の自発放出率 A：E2・M1・GW）。

準位
  核の動径方程式  −(1/2μ_n) f″ + [V(R) + J(J+1)/(2μ_n R²)] f = E f、μ_n = m_p/2（断熱近似、Roueff+2019 と同じ）
  V(R) = E_el(R) + E_a(R)
    E_el : BO ポテンシャル。Pachucki, PRA 82, 032509 (2010) に基づく解析フィット（G. Łach）。H2SPECTRE 7.4 data/h2spectr_pot.f90 の FUNCTION V を移植。
           V(∞) = 0（2 個の H 原子を原点）。h2_data/h2_V_BO_pachucki2010.dat（原論文の表、82 点）との差を validation_h2.py で検査。
    E_a  : 断熱補正。Pachucki, Komasa, JCP 141, 224103 (2014) のフィット（H2SPECTRE 7.4 の Eaint = 2μ_n E_a）。E_a(∞) = 0（2H を原点）。
  離散化  sinc-DVR（Colbert–Miller；H2SPECTRE と同じ行列要素）、R_i = i·dr、i = 1..N。束縛準位 = E < 0 の固有値。
放出率（原子単位、ħ = m_e = e = 1、c = 1/α；時間単位 2.4188843265857e-17 s）
  E2 : A = (α⁵ ω⁵/15) (2J_f+1) (J_i 2 J_f; 0 0 0)² |⟨f_f|Θ|f_i⟩|²   ΔJ = 0, ±2（J_i+J_f 偶）
       Θ(R) = 分子四重極能率 = Q_WSD(R)/2（Wolniewicz, Simbotin, Dalgarno, ApJS 115, 293 (1998) 表 3；Roueff+2019 §3.1 の注）
       Q 枝で (2J+1)(J 2 J;000)² = J(J+1)/((2J−1)(2J+3))：Pachucki–Komasa 2011 式 (14) と一致
  M1 : A = (α⁵ ω³/3) (m_e/m_p)² J(J+1) |⟨f_f|g|f_i⟩|²               ΔJ = 0, v_i ≠ v_f
       g(R) = 回転 g 因子（Pachucki, Komasa, PRA 83, 032501 (2011) 表 I）。Roueff+2019 式：8.00e-18 σ³ J(J+1)|⟨g⟩|² と同値
  GW : A = 4 (α_Gp/α) |⟨f_f|q_m|f_i⟩|² / |⟨f_f|Θ|f_i⟩|² × A_E2        （LL §110 / §71 の比、水素一式の W_GW と同じ導出）
       q_m(R) = 質量四重極/m_p = R²/2 + (m_e/m_p)(R²/2 − Θ(R))        [m_p a₀²]、α_Gp = G m_p²/(ħc)
核スピン：J 偶 = パラ（I = 0）、J 奇 = オルト（I = 1）。E2・M1 は J の偶奇を保つので両系列は結ばれない。
         I の縮退は src・dst で共通なので詳細釣り合いの g_src/g_dst = (2J_s+1)/(2J_d+1) はそのまま。
質量：M = 2(m_p + m_e) + E/c²（反跳用）。
"""
import sys
from pathlib import Path
import numpy as np
from scipy.special import gammainc
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import engine as EN  # noqa: E402

# ---------- 定数（CODATA 2018） ----------
MPROT = 1836.15267343          # m_p/m_e
ALPHA = 7.2973525693e-3
HARTREE_EV = 27.211386245988
HARTREE_CM = 219474.6313632    # cm⁻¹ / hartree
AU_TIME_S = 2.4188843265857e-17
G_SI = 6.67430e-11
MP_SI = 1.67262192369e-27
ALPHA_GP = G_SI * MP_SI**2 / (EN.HBAR_SI * EN.C)   # G m_p²/(ħc) = 5.906e-39
MU_N = MPROT / 2.0             # 核の換算質量 [m_e]
E_B_H_EV = 13.598434599702     # 水素原子の束縛エネルギー（質量欠損）
M_H_ME = MPROT + 1.0 - E_B_H_EV / (EN.ME_SI * EN.C**2 / EN.E_SI)   # 水素原子の質量 [m_e] = m_p + m_e − E_b/c²
DATA = HERE / "h2_data"


# ---------- BO ポテンシャル E_el(R)（H2SPECTRE 7.4 FUNCTION V の移植） ----------
_EHE = -2.903724377034119
_BB = 1.84024290611525
_C = {6: 6.4990267054058393131, 7: 0.0, 8: 124.3990835836223436, 9: 0.0, 10: 3285.828414967421697,
      11: -3474.898037882, 12: 122727.6087007, 13: -326986.9240441, 14: 6361736.045092, 15: -28395580.63300,
      16: 441205192.2739, 17: -2739281653.140, 18: 39352477334.60, 19: -307082459389.3, 20: 4363762779418.0,
      21: -40360383996600.0, 22: 586030516404600.0, 23: -6200210182655000.0, 24: 93433616652420000.0,
      25: -1105388541580000000.0, 26: 17414122541340000000.0, 27: -226839565464200000000.0,
      28: 3747403197985000000000.0, 29: -53148855088090000000000.0, 30: 921579197325600000000000.0}
_B = [0.0, -4.57323700588624794726213559410563, 13.47638719366600967683673663479999, -28.51530241703605266394033745126115,
      45.78475015811747771910888974954924, -58.19615893167528525478814306565867, 58.48593413367164587893380279407919,
      -45.35039381018086925859899857147845, 26.47049312447772391118586407437590, -11.37752030083640364078327219266986,
      3.50402820299159351290056244127085, -0.73796135013579258240464270522277, 0.09631412647157079742433223552984,
      -0.00576418363156423769040137255530, -0.00014088867788010427947798439255, 0.00002845940482698706139384277587]
_A = [0.0] * 23
_A[1] = _B[1]
_A[2] = _B[2] + 2 + (_EHE + 1)
_A[3] = _B[3] + _A[2] * _B[1] - _B[1] * _B[2]
_A[4] = _B[4] + 2 * (_A[2] - _B[2] - 1) + _B[2] * (_A[2] - _B[2])
_A[5] = _B[5] + 2 * _B[1] * (_A[2] - _B[2] - 1) + _B[3] * (_A[2] - _B[2])
_A[6:] = [39.91391743339803899150896228191386, -5.75063757195175050513832273504327, -33.46588381771842371982158706320471,
          55.32089352568815354527823571430595, -51.42284107211464070094798429642712, 32.71776925369248112208919743003516,
          -15.66016122411670563080887658243106, 7.50377730360456765232166233248396, -5.40306648644156257500314013983694,
          4.47594305162345652142325458988993, -2.95983991352400754745996259257977, 1.42669802756450786769491630785461,
          -0.49306757735153614794116509830620, 0.11952677096537517464084697055123, -0.01932309354268259675592903647346,
          0.00186812017042699683819835094121, -0.00008146220049329661292603862566]


def tang_toennies(n, x):
    """F_n(x) = 1 − Γ(n+1, x)/Γ(n+1) = P(n+1, x)（正則化下側不完全ガンマ）。"""
    return gammainc(n + 1, x)


def v_el(R):
    """BO ポテンシャル E_el(R) [hartree]、E_el(∞) = 0。"""
    R = np.asarray(R, dtype=float)
    sr = np.sqrt(R)
    t1 = sum(_A[i] * sr**i for i in range(1, 23))
    t2 = sum(_B[i] * sr**i for i in range(1, 16))
    t3 = -sum(tang_toennies(i + 1, _BB * R) * _C[i] / R**i for i in range(6, 31))
    return (1 + t1) / (1 + t2) * np.exp(-2 * R) / R + t3


# ---------- 断熱補正 E_a(R)（H2SPECTRE 7.4 Eaint = 2μ_n E_a の移植、Pachucki–Komasa 2014） ----------
_AD_A, _AD_B, _AD_C = 1.1114055293534456784057320121522, 1.2095560917705853428726384695475, 1.01490515805955384474660168480738
_AD_E0 = 0.5313969260599797962744365
_AD_P = [2.57640985968790283123151781491906e4, -4.24766902765405008898541826259253e4, 3.23914101790994120219160120634889e4,
         -1.51487736617849601347139327355656e4, 4.85513784701655588907660942038169e3, -1.12817968586590251891697217394290e3,
         1.96022115958782784119305229987615e2, -2.58727150904267703054234804239174e1, 2.60684632154178626761090277540466e0,
         -1.99532537310733713612913792065941e-1, 1.14207904980296991756950227054134e-2, -4.74017092281636746493356165838112e-4,
         1.34857652007195080104183709321439e-5, -2.35400021232008399283943817140989e-7, 1.90202486167843259334945031635423e-9]
_AD_Q = [0.0, 3.99485699616875103429516544309691e4, -5.44880006664402627114573753359117e4, 5.26679590123087863666126645727291e4,
         -4.53948569025655989594911809433954e4, 3.33135367104566309045409601884433e4, -2.25922524785635131591735975775914e4,
         1.56343244224922385234291816278842e4, -1.60616168851563262719458549830300e4, 2.73228697168824991946522297552660e4,
         -5.16780392742578518472987095005810e4, 8.66782887179522293314013939749623e4, -1.23395802058510929643696213689536e5,
         1.49032029889160972833144538424413e5, -1.53587670024503730585958769007304e5, 1.35761994668913273821118458600299e5,
         -1.03311683808839637492445134792404e5, 6.78372389519334445654343720586628e4, -3.84767394910856453874838129483377e4,
         1.88499722967195307004803452140572e4, -7.96690860669712153936199181841002e3, 2.89804982597848873059044487914495e3,
         -9.03966166931006732358453963336031e2, 2.40513643411929918936145610873524e2, -5.41915909139083735974492686116770e1,
         1.02397672252399842815397491089614e1, -1.60129916283328421242788004210142e0, 2.03505145495990689451837401340804e-1,
         -2.04824339024431528193748270566111e-2, 1.57093891708355488534614263836794e-3, -8.62727960846911255035254246400885e-5,
         3.02264161462310781046446491735097e-6, -5.07959453408875288585515001119814e-8]
_AD_MP, _AD_A6, _AD_A8, _AD_A10 = 1836.15267247, 1.7699e-2, 0.144, 2.28


def eaint(R):
    """2μ_n E_a(R)（H2SPECTRE の Eaint）。"""
    R = np.asarray(R, dtype=float)
    w1 = np.polyval(_AD_P[::-1], R)
    w2 = (_AD_E0 - _AD_P[0]) + R * np.polyval(_AD_Q[:0:-1], R)
    near = np.exp(-_AD_A * R) * w1 + np.exp(-(_AD_B + _AD_C * R) * R) * w2
    rinv2 = 1.0 / (R * R)
    far = -_AD_MP * (_AD_A6 + (_AD_A8 + _AD_A10 * rinv2) * rinv2) * rinv2**3
    return np.where(R <= 11.5, near, far)


def e_ad(R):
    """断熱補正 E_a(R) [hartree]、E_a(∞) = 0（2H を原点）。"""
    return eaint(R) / MPROT


# ---------- Θ(R)、g(R)（文献の表を補間） ----------
def _load(name):
    d = np.loadtxt(DATA / name)
    return d[:, 0], d[:, 1]


_RQ, _QW = _load("h2_Q_WSD1998.dat")
_THETA_SPL = CubicSpline(_RQ, 0.5 * _QW)
_RG, _GG = _load("h2_g_PK2011.dat")
_G_SPL = CubicSpline(_RG, _GG)


def theta(R):
    """分子四重極能率 Θ(R) [e a₀²] = Q_WSD/2。表の外：R < 0.2 は ∝ R²、R > 20 は ∝ R⁻⁶（表の末尾の減衰 (19/20)⁶ に一致）。"""
    R = np.asarray(R, dtype=float)
    lo, hi = _RQ[0], _RQ[-1]
    inner = _THETA_SPL(np.clip(R, lo, hi))
    return np.where(R < lo, _THETA_SPL(lo) * (R / lo) ** 2, np.where(R > hi, _THETA_SPL(hi) * (hi / R) ** 6, inner))


def g_rot(R):
    """回転 g 因子 g(R)。表の外（R > 12）は末尾 2 点の指数減衰で外挿（値 5e-6 以下）。"""
    R = np.asarray(R, dtype=float)
    hi = _RG[-1]
    k = np.log(_GG[-2] / _GG[-1]) / (_RG[-1] - _RG[-2])
    return np.where(R > hi, _GG[-1] * np.exp(-k * (R - hi)), _G_SPL(np.clip(R, _RG[0], hi)))


def q_mass(R):
    """質量四重極/m_p [m_p a₀²]：陽子 R²/2 ＋ 電子 (m_e/m_p)(R²/2 − Θ)。"""
    R = np.asarray(R, dtype=float)
    return 0.5 * R * R + (0.5 * R * R - theta(R)) / MPROT


# ---------- 準位（sinc-DVR） ----------
def dvr_levels(J, dr=0.025, rmax=50.0, pot="ad", emax=0.0):
    """J を固定した準位。戻り値：E（hartree、2H 原点、昇順、E < emax）、F（各列が DVR 固有ベクトル、Σ F² = 1）、R 格子。
    emax = 0 で束縛準位だけ。emax > 0 なら箱（R ≤ rmax）で規格化した離散化連続状態も返す（放射会合に使う）。"""
    N = int(round(rmax / dr))
    i = np.arange(1, N + 1)
    R = i * dr
    V = v_el(R) + (e_ad(R) if pot == "ad" else 0.0) + J * (J + 1) / (2 * MU_N * R * R)
    I, K = np.meshgrid(i, i, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        T = (-1.0) ** ((I - K) % 2) / (MU_N * dr * dr) * (1.0 / (I - K) ** 2 - 1.0 / (I + K) ** 2)
    T[np.arange(N), np.arange(N)] = (np.pi**2 / 3 - 0.5 / i**2) / (2 * MU_N * dr * dr)
    H = T + np.diag(V)
    E, F = eigh(H, subset_by_value=(-np.inf, emax))
    return E, F, R


def all_levels(dr=0.025, rmax=50.0, pot="ad", jmax=31):
    """全 J の束縛準位。戻り値：list of dict(v, J, E_h, f) と R 格子。"""
    out = []
    R = None
    for J in range(jmax + 1):
        E, F, R = dvr_levels(J, dr, rmax, pot)
        for v in range(len(E)):
            out.append(dict(v=v, J=J, E_h=float(E[v]), f=F[:, v]))
    return out, R


# ---------- 原子間の重力（別勘定） ----------
def gravity_level(lev, R):
    """準位 (v,J) の原子間重力の一次摂動 ⟨f| −G m_H²/R |f⟩ [hartree]。
    点質量の原子（m_H = m_p + m_e − E_b/c²）どうしの Newton 項だけ。入れていないもの：電子–核・電子–電子の重力の交差項
    （Σ⟨1/r_ia⟩(R)、⟨1/r_12⟩(R) が要る。核–核項に対して ~1e-3、~1e-7）、1PN（(v/c)² ~ 1e-10）。
    大きさは 1e-31 cm⁻¹ で、準位エネルギーの倍精度の床（1e-12 cm⁻¹）の 20 桁下。E_eV に加えても値は変わらないので別勘定で持つ。"""
    inv_r = float(np.dot(lev["f"] * lev["f"], 1.0 / R))
    return -(ALPHA_GP / ALPHA) * (M_H_ME / MPROT) ** 2 * inv_r


# ---------- 放出率 ----------
def rates(lev_i, lev_f, R):
    """準位 i → f の (A_E2, A_M1, A_GW) [s⁻¹]。選択則を満たさない成分は 0。"""
    Ji, Jf = lev_i["J"], lev_f["J"]
    w = lev_i["E_h"] - lev_f["E_h"]
    if w <= 0.0:
        return 0.0, 0.0, 0.0
    A_e2 = A_m1 = A_gw = 0.0
    if abs(Ji - Jf) in (0, 2):
        ang = (2 * Jf + 1) * EN.w3j000_sq(Ji, 2, Jf)
        me_th = float(np.dot(lev_f["f"], lev_i["f"] * theta(R)))
        me_qm = float(np.dot(lev_f["f"], lev_i["f"] * q_mass(R)))
        A_e2 = ALPHA**5 * w**5 / 15.0 * ang * me_th**2 / AU_TIME_S
        A_gw = 4.0 * ALPHA_GP * ALPHA**4 * w**5 / 15.0 * ang * me_qm**2 / AU_TIME_S
    if Ji == Jf and lev_i["v"] != lev_f["v"]:
        me_g = float(np.dot(lev_f["f"], lev_i["f"] * g_rot(R)))
        A_m1 = ALPHA**5 * w**3 / 3.0 * (1.0 / MPROT) ** 2 * Ji * (Ji + 1) * me_g**2 / AU_TIME_S
    return A_e2, A_m1, A_gw


def build(init=(1, 3), dr=0.025, rmax=50.0, pot="ad", jmax=31, with_gw=True):
    """エンジン用の states, channels, init_index（と補助情報 aux）。"""
    levels, R = all_levels(dr, rmax, pot, jmax)
    levels.sort(key=lambda d: (d["J"], d["v"]))
    states = []
    for d in levels:
        I = 1 if d["J"] % 2 else 0
        Eg = gravity_level(d, R)
        states.append(dict(label=(d["v"], d["J"]), name="v=%d J=%d (%s)" % (d["v"], d["J"], "ortho" if I else "para"),
                           E_eV=d["E_h"] * HARTREE_EV, J=d["J"], I=I, E_grav_eV=Eg * HARTREE_EV, E_grav_cm=Eg * HARTREE_CM,
                           M_kg=(2 * (MPROT + 1) + d["E_h"] * ALPHA**2) * EN.ME_SI, v=d["v"]))
    channels = []
    for a, la in enumerate(levels):
        for b, lb in enumerate(levels):
            if la["E_h"] <= lb["E_h"]:
                continue
            e2, m1, gw = rates(la, lb, R)
            if e2 > 0.0:
                channels.append((a, b, "e2", 2, e2))
            if m1 > 0.0:
                channels.append((a, b, "m1", 1, m1))
            if with_gw and gw > 0.0:
                channels.append((a, b, "gw", 2, gw))
    idx = {s["label"]: k for k, s in enumerate(states)}
    if tuple(init) not in idx:
        raise SystemExit("初期準位 (v,J)=%s は束縛準位に無い" % (tuple(init),))
    aux = dict(R=R, levels=levels, idx=idx, dr=dr, rmax=rmax, pot=pot)
    return states, channels, idx[tuple(init)], aux


if __name__ == "__main__":
    states, channels, init, aux = build()
    print("levels %d, channels %d (e2 %d, m1 %d, gw %d), init %s" % (
        len(states), len(channels), sum(c[2] == "e2" for c in channels), sum(c[2] == "m1" for c in channels),
        sum(c[2] == "gw" for c in channels), states[init]["name"]))
    print("D0 [cm^-1] = %.4f" % (-states[aux["idx"][(0, 0)]]["E_eV"] / HARTREE_EV * HARTREE_CM))
