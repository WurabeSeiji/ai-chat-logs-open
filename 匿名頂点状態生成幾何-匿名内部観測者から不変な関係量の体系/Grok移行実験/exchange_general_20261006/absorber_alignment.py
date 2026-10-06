#!/usr/bin/env python3
"""層 3：m 副準位の整列と偏光（STATUS §12 の三層のうち、状態空間・生成子・読み出しの全部が変わる層）。層 2（温度）を含む。

基底：固定一式の 50 準位 (n, ℓ, j, F) を m_F まで分解 → Σ(2F+1) = 220 状態。量子化軸 ẑ = 吸収体の印（双極子軸）。
率：Wigner–Eckart。多重極の階数 k（E1、M1 は 1、E2・重力波は 2、二光子は 0）に対し
      A_{(F,m)→(F',m')} = A_{F→F'} (2F+1) (F' k F; m' q −m)²、q = m − m'（Σ_{m'} = A_{F→F'}）。
角分布（光子の偏光で和）：k = 1：P₀(θ) = (3/8π) sin²θ、P_{±1}(θ) = (3/16π)(1 + cos²θ)。
                        k = 2：P₀ = (15/8π) sin²θ cos²θ、P_{±1} = (5/16π)(1 − 3cos²θ + 4cos⁴θ)、P_{±2} = (5/16π)(1 − cos⁴θ)。
偏光（k = 1、θ = 90°）：π（q = 0）は ẑ に平行な直線偏光、σ（q = ±1）は垂直。直線偏光度 P = (2N₀ − N₊ − N₋)/(2N₀ + N₊ + N₋)。
背景放射：T(n̂_src) = T₀/(γ(1 − β n̂_src·ẑ)) · (1 + a₂ P₂(n̂_src·ẑ))。β = 370 km/s の Doppler 双極子（その 2 次が四重極を作る）、
          a₂ = 4e-6 は観測されている背景放射の四重極の目安（実際の四重極の軸は双極子と違うが、ここでは同じ軸に置く）。
          Δm = q の遷移が感じる占有数は型で重みづけた ⟨n̄⟩_q = ∫ n̄(ω; T(n̂)) P_q(n̂) dΩ。m 分解では副準位が縮退していないので
          詳細釣り合いは A_{m'm}(1 + ⟨n̄⟩_q)（下向き）と A_{m'm}⟨n̄⟩_q（上向き）、g 因子なし。等方極限で m 分解前の率に戻る。
読み出し（記録だけから）：
  (A) 初期 5p₃/₂ F=1 を m_F = 0（または +1）に置く。多段の光子の線ごとの (N₀, N₊, N₋) → 角分布 a + b cos²θ と偏光度。
      m が記録に出る。m 平均の初期状態なら等方・無偏光に戻る（検算）。
  (B) T = 2.7255 K の非等方背景での 1s F=1 の定常 m 分布 → 整列 A₂ = [p₊ + p₋ − 2p₀]/Σp と 21 cm 放出の偏光度。
      等方背景（β = 0、a₂ = 0）で 0 に戻る（検算）。
出力：absorber_alignment.md
"""
import sys
from fractions import Fraction
from math import factorial, sqrt
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
FIXED = HERE.parent / "exchange_rel_20261005"
sys.path.insert(0, str(FIXED))
from matplotlib import font_manager  # noqa: E402
font_manager.fontManager.addfont = lambda path: None
import exchange_cascade as X  # noqa: E402

C = X.C_SI
K_B_EV = 8.617333262e-5
T_CMB = 2.7255
V0 = 370.0e3
A2_CMB = 4.0e-6
EV = X.EV
GL_X, GL_W = np.polynomial.legendre.leggauss(96)
RANK = {"e1": 1, "m1": 1, "e2": 2, "gw": 2, "2ph": 0}


# ---------- 3j 記号（Racah、有理数で厳密） ----------
def _fact(x):
    return factorial(int(x))


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
    val = sqrt(float(pref)) * float(s) * (-1) ** int(j1 - j2 - m3)
    return val


# ---------- m 分解の基底 ----------
STATES = []           # (level index i, F, m)
for i, (n, l, se, sp) in enumerate(X.ST):
    j, F = X.jf_of(l, se, sp)
    for mm in range(-int(F), int(F) + 1):
        STATES.append((i, F, mm))
DIM = len(STATES)
SIDX = {s: k for k, s in enumerate(STATES)}
LEVEL_OF = np.array([s[0] for s in STATES])
ETOT_M = X.ETOT[LEVEL_OF]


def pattern(k, q, x):
    """角分布 P_q(θ)、x = cos θ。∫P dΩ = 1。"""
    if k == 1:
        return (3 / (8 * np.pi)) * (1 - x * x) if q == 0 else (3 / (16 * np.pi)) * (1 + x * x)
    if k == 2:
        if q == 0:
            return (15 / (8 * np.pi)) * (1 - x * x) * x * x
        if abs(q) == 1:
            return (5 / (16 * np.pi)) * (1 - 3 * x * x + 4 * x**4)
        return (5 / (16 * np.pi)) * (1 - x**4)
    return 1 / (4 * np.pi)


def nbar(hw_eV, T):
    if T <= 0.0:
        return 0.0
    x = hw_eV / (K_B_EV * T)
    return 0.0 if x > 700.0 else 1.0 / np.expm1(x)


def T_of_x(T0, beta, a2, x):
    g = 1.0 / np.sqrt(1.0 - beta**2)
    return T0 / (g * (1.0 - beta * x)) * (1.0 + a2 * 0.5 * (3 * x * x - 1))


def nbar_q(hw_eV, T0, beta, a2, k, q):
    """型で重みづけた占有数 ⟨n̄⟩_q = ∫ n̄(ω; T(n̂_src)) P_q(n̂) dΩ（P_q は cos θ の偶関数なので源の向きと進行方向の区別は要らない）。"""
    if T0 <= 0.0:
        return 0.0
    vals = np.array([nbar(hw_eV, T_of_x(T0, beta, a2, x)) * pattern(k, q, x) for x in GL_X])
    return 2 * np.pi * float(np.dot(GL_W, vals))


def build_edges_m():
    """m 分解した辺：(src_state, dst_state, kind, k, q, A_mm')。Σ_{m'} A_mm' = A_{F→F'} を確認。"""
    ed = X.edges(X.W_GW)
    out = []
    worst = 0.0
    for src, dst, kind, A in ed:
        if A <= 0.0:
            continue
        k = RANK[kind]
        Fs, Fd = X.jf_of(*X.ST[src][1:])[1], X.jf_of(*X.ST[dst][1:])[1]
        for m in range(-int(Fs), int(Fs) + 1):
            tot = 0.0
            for m2 in range(-int(Fd), int(Fd) + 1):
                q = m - m2
                if abs(q) > k:
                    continue
                if k == 0:
                    Amm = A if (m2 == m and Fd == Fs) else 0.0
                else:
                    Amm = A * (2 * Fs + 1) * wigner3j(Fd, k, Fs, m2, q, -m) ** 2
                if Amm > 0.0:
                    out.append((SIDX[(src, Fs, m)], SIDX[(dst, Fd, m2)], kind, k, q, Amm))
                    tot += Amm
            worst = max(worst, abs(tot - A) / A)
    return out, worst


def generator_m(edm, T0, beta, a2):
    R = np.zeros((DIM, DIM))
    for s, d, kind, k, q, A in edm:
        hw = (ETOT_M[s] - ETOT_M[d]) * EV
        nq = nbar_q(hw, T0, beta, a2, k, q) if kind in ("e1", "e2", "m1") else 0.0
        R[d, s] += A * (1.0 + nq)
        R[s, s] -= A * (1.0 + nq)
        if nq > 0.0:
            R[s, d] += A * nq
            R[d, d] -= A * nq
    return R


def occupation_time(R, p0):
    out_rate = -np.diag(R)
    trans = [i for i in range(DIM) if out_rate[i] > 0.0]
    absb = [i for i in range(DIM) if out_rate[i] == 0.0]
    tau = np.zeros(DIM)
    tau[trans] = np.linalg.solve(-R[np.ix_(trans, trans)], p0[trans])
    p_inf = p0.copy()
    p_inf[trans] = 0.0
    p_inf[absb] += np.dot(R[np.ix_(absb, trans)], tau[trans])
    return tau, p_inf


def line_photons(edm, tau):
    """線（準位対）ごとの (N₀, N₊, N₋)（k = 1）、k = 2 は (N₀, N_{±1}, N_{±2})。"""
    acc = {}
    for s, d, kind, k, q, A in edm:
        key = (LEVEL_OF[s], LEVEL_OF[d], kind)
        acc.setdefault(key, {})
        acc[key][q] = acc[key].get(q, 0.0) + A * tau[s]
    return acc


def lab(i):
    n, l, se, sp = X.ST[i]; j, F = X.jf_of(l, se, sp)
    return "%d%s_{%d/2}F%d" % (n, {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}[l], int(2 * j), int(F))


def angular_coeffs(k, Nq):
    """I(θ) = Σ_q N_q P_q(θ) を a + b x² (+ c x⁴) で。90° の直線偏光度（k = 1）。"""
    if k == 1:
        N0 = Nq.get(0, 0.0); Ns = Nq.get(1, 0.0) + Nq.get(-1, 0.0)
        a = (3 / (8 * np.pi)) * N0 + (3 / (16 * np.pi)) * Ns
        b = -(3 / (8 * np.pi)) * N0 + (3 / (16 * np.pi)) * Ns
        P90 = (2 * N0 - Ns) / (2 * N0 + Ns) if (2 * N0 + Ns) > 0 else float("nan")
        return a, b, P90
    return None


def main():
    L = []
    L.append("# 層 3：m 副準位の整列と偏光（層 2 を含む）")
    L.append("")
    edm, worst = build_edges_m()
    L.append("基底 220 状態（50 準位 × (2F+1)）、m 分解した辺 %d 本。和則 Σ_{m'} A_{mm'} = A_{F→F'} の最大相対差 %.1e。量子化軸 ẑ = 吸収体の印。" % (len(edm), worst))
    L.append("3j の検算：(1 1 0; 0 0 0)² = %.6f（1/3）、(1 1 2; 1 −1 0)² = %.6f（1/15）。" % (wigner3j(1, 1, 0, 0, 0, 0) ** 2, wigner3j(1, 1, 2, 1, -1, 0) ** 2))
    L.append("")
    # (A) 多段：T = 0、初期 m_F
    R0 = generator_m(edm, 0.0, 0.0, 0.0)
    L.append("## (A) 多段の光子の角分布と偏光（T = 0、初期 5p₃/₂ F=1 の m_F を指定）")
    L.append("")
    i5p = X.IDX[(5, 1, 1, -1)]
    results = {}
    for tag, weights in (("m_F = 0", {0: 1.0}), ("m_F = +1", {1: 1.0}), ("m 平均（分解前と同じ）", {-1: 1 / 3, 0: 1 / 3, 1: 1 / 3})):
        p0 = np.zeros(DIM)
        for mm, w in weights.items():
            p0[SIDX[(i5p, 1.0, mm)]] = w
        tau, p_inf = occupation_time(R0, p0)
        ph = line_photons(edm, tau)
        emitted = sum(A * tau[s] * (ETOT_M[s] - ETOT_M[d]) for s, d, kind, k, q, A in edm) * EV
        results[tag] = (tau, p_inf, ph, emitted)
        L.append("### 初期 %s：放出合計 %.9f eV、終状態 %s" % (tag, emitted, [(lab(LEVEL_OF[i]), STATES[i][2], round(float(p_inf[i]), 6)) for i in range(DIM) if p_inf[i] > 1e-12]))
        L.append("")
        L.append("| 線 | 種別 | 光子数の期待値 N | N₀ : N₊ : N₋ | 角分布 I(θ) ∝ 1 + (b/a) cos²θ | 90° の直線偏光度 |")
        L.append("|---|---|---|---|---|---|")
        rows = sorted(ph.items(), key=lambda kv: -sum(kv[1].values()))
        for (ls, ld, kind), Nq in rows[:8]:
            Ntot = sum(Nq.values())
            if Ntot < 1e-6 or RANK[kind] != 1:
                continue
            a, b, P90 = angular_coeffs(1, Nq)
            L.append("| %s → %s | %s | %.6f | %.4f : %.4f : %.4f | 1 %+.4f cos²θ | %+.4f |"
                     % (lab(ls), lab(ld), kind, Ntot, Nq.get(0, 0) / Ntot, Nq.get(1, 0) / Ntot, Nq.get(-1, 0) / Ntot, b / a, P90))
        L.append("")
    L.append("m_F = 0 の 5p₃/₂F1 → 1s F0 は q = 0 だけ（純 π）なので I ∝ sin²θ、偏光度 +1。m 平均では全ての線が等方（b/a = 0）・無偏光（P = 0）で、放出合計は m 分解前の 13.054539268 eV に戻る。"
             "初期の m が、放出光の角分布と偏光として記録に出る。")
    L.append("")
    # (B) 定常整列
    L.append("## (B) 非等方背景（T = 2.7255 K）での 1s F=1 の定常整列と 21 cm の偏光")
    L.append("")
    i1s1 = X.IDX[(1, 0, 1, 1)]
    hw21 = (X.ETOT[i1s1] - X.ETOT[X.IDX[(1, 0, 1, -1)]]) * EV
    L.append("| 背景 | ⟨n̄⟩₀（π） | ⟨n̄⟩₊₁（σ） | 1s F=1 の p(m=−1), p(0), p(+1) | 整列 A₂ = (p₊+p₋−2p₀)/Σ | 21 cm 放出の 90° 偏光度 |")
    L.append("|---|---|---|---|---|---|")
    p0 = np.zeros(DIM)
    p0[SIDX[(i5p, 1.0, 0)]] = 1.0
    times = np.array([0.0, 1.0e17])
    for name, beta, a2 in (("等方", 0.0, 0.0), ("双極子 β のみ（2 次で四重極）", V0 / C, 0.0), ("双極子 β ＋ 四重極 a₂ = 4e-6", V0 / C, A2_CMB)):
        R = generator_m(edm, T_CMB, beta, a2)
        sol = solve_ivp(lambda t, p: np.dot(R, p), (0.0, 1.0e17), p0, method="Radau", jac=lambda t, p: R, t_eval=times, rtol=1e-12, atol=1e-30)
        pT = sol.y[:, -1]
        pm = {mm: pT[SIDX[(i1s1, 1.0, mm)]] for mm in (-1, 0, 1)}
        S = sum(pm.values())
        A2 = (pm[1] + pm[-1] - 2 * pm[0]) / S
        # 21 cm 放出：F=1,m → F=0,0、q = m。偏光度 = (2N₀ − N₊ − N₋)/(2N₀ + N₊ + N₋)、N_q ∝ p_m A（A は m によらない）
        P90 = (2 * pm[0] - pm[1] - pm[-1]) / (2 * pm[0] + pm[1] + pm[-1])
        L.append("| %s | %.8f | %.8f | %.9f, %.9f, %.9f | %+.3e | %+.3e |"
                 % (name, nbar_q(hw21, T_CMB, beta, a2, 1, 0), nbar_q(hw21, T_CMB, beta, a2, 1, 1), pm[-1], pm[0], pm[1], A2, P90))
    L.append("")
    L.append("等方背景では整列 0（倍精度）。無偏光の強度双極子は 1 次では整列を作らず、2 次 β² ≈ 1.5e-6 の四重極成分と背景の四重極 a₂ が作る。"
             "定常の副準位の占有は p_m ∝ ⟨n̄⟩_q/(1 + ⟨n̄⟩_q) なので、占有数の相対差 Δn̄/n̄ ≈ 10⁻⁶ は 1/(n̄(1+n̄)) ≈ 1/1600 だけ抑えられ、整列は 10⁻⁸ 台になる"
             "（n̄ ≈ 40 で遷移が飽和に近いため。事前の見積もり 10⁻⁶ はこの抑制を見落としていた）。21 cm 放出の直線偏光度 2×10⁻⁸ を読むには 10¹⁶ 個の光子が要る。")
    L.append("")
    L.append("## 何が増えたか")
    L.append("")
    L.append("- 状態空間が m_F まで分解され（220）、生成子は Δm = q ごとの率と、背景の型ごとの占有数 ⟨n̄⟩_q になった。")
    L.append("- 記録の各光子に Δm の型（角分布・偏光）が付く。初期の m は多段光子の角分布と偏光にそのまま出る（(A)）。")
    L.append("- 非等方な背景は m 副準位を不均等に汲み上げ、定常の整列と 21 cm の偏光を作る（(B)）。大きさは背景の四重極成分 10⁻⁶ を飽和因子 1/(n̄(1+n̄)) で抑えた 10⁻⁸ 台。")
    L.append("- 入れていないもの：E2・重力波の型での汲み上げ（n̄ = 0 の振動数なので効かない）、二光子の偏光相関（対の相関は k = 0 の演算子で ΔF = 0 のまま）、円偏光（背景に円偏光は無い）。")
    out = "\n".join(L) + "\n"
    (HERE / "absorber_alignment.md").write_text(out, encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
