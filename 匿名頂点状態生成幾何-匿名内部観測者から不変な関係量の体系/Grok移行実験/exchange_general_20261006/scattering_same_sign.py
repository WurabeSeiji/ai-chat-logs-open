#!/usr/bin/env python3
"""同符号の対（ee、pp）の散乱の枝：in–out の関係量だけを計算する。軌道は出力しない。

条件は水素の一式と揃える：相対運動エネルギー E = |E(5p₃/₂ F=1)|（固定一式の ETOT から）、
スペクトルの衝突径数は ℓ = 1 の角運動量 L = ħ√2 から b = L/(μv)。

1. 断面積（厳密、クーロン）：Rutherford の振幅 f(θ) = −(η/2k sin²(θ/2)) exp[−iη ln sin²(θ/2) + 2iσ₀]、
   η = (q_aq_b/ħc)/(v/c)（Sommerfeld 径数）、σ₀ = arg Γ(1+iη)。同種スピン ½ 非偏極（Mott 1930）：
   dσ/dΩ = |f(θ)|² + |f(π−θ)|² − Re f(θ)f*(π−θ)
         = (η/2k)² [1/sin⁴(θ/2) + 1/cos⁴(θ/2) − cos(η ln tan²(θ/2))/(sin²(θ/2)cos²(θ/2))]。
   区別できる粒子（古典）の両粒子を数えた値 (η/2k)²[1/sin⁴ + 1/cos⁴] と比べる。θ = 90° で比 1/2（反対称化）。
2. スペクトル（古典、四重極）：同種粒子は重心系の電気双極子が恒等的に 0。放射は四重極
   P = (1/180c⁵) D⃛_{αβ}D⃛_{αβ}、D_{αβ} = Σ q(3x_αx_β − r²δ_αβ)（Landau–Lifshitz §71）。重心系で位置 ±r/2 なので D = (q/2)(3rr − r²δ)。
   斥力の双曲線 x = a(e + cosh ξ)、y = a√(e²−1) sinh ξ、t = (a/v)(e sinh ξ + ξ)、a = q_aq_b/(μv²)、e = √(1 + (b/a)²)。
   D⃛ は運動方程式 r̈ = (q_aq_b/μ) r/r³ から閉じた形で、dE/dω = (1/180πc⁵) Σ|D⃛_{αβ}(ω)|²、D⃛(ω) = ∫D⃛(t)e^{iωt}dt を ξ で数値積分。
   Parseval ∫P dt = ∫(dE/dω)dω を自己検査。重力波は質量四重極で同じ計算（比は 4Gm²/q² になるはず）。
   古典のスペクトルが妥当なのは η ≫ 1（pp）。ee（η ≈ 3.5）では目安（量子の Gaunt 因子は入れていない）。
3. 磁気（Darwin）項の大きさ v²/c²。

出力：scattering_same_sign.md
"""
import sys
from pathlib import Path
import numpy as np
from scipy.special import loggamma
from scipy.integrate import simpson

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import two_body_general as T  # noqa: E402
X = T.X
ALPHA, ME_EV = X.ALPHA, X.ME_EV
HBAR_EVS = X.HBAR_SI / X.E_SI               # ħ [eV s]
A0_M = X.HBAR_SI / (X.ME_SI * X.C_SI * ALPHA)   # Bohr 半径（電子質量）[m]
E_REL = abs(X.ETOT[X.IDX[(5, 1, 1, -1)]]) * ME_EV   # 相対運動エネルギー [eV]：水素の初期準位 5p₃/₂ F=1 の束縛エネルギー
ELL = 1                                          # 5p の ℓ
L_HBAR = np.sqrt(ELL * (ELL + 1))                # L/ħ


def mott(theta, eta, k):
    s2, c2 = np.sin(theta / 2) ** 2, np.cos(theta / 2) ** 2
    pref = (eta / (2 * k)) ** 2
    ruth_both = pref * (1 / s2**2 + 1 / c2**2)
    mott_ = pref * (1 / s2**2 + 1 / c2**2 - np.cos(eta * np.log(s2 / c2)) / (s2 * c2))
    return mott_, ruth_both


def quadrupole_spectrum(e, n_xi=20001, xi_max=14.0, omegas=None):
    """無次元（a = 1、v = 1、q = 1、位置 ±r/2）。戻り値：ω̃、Σ|D̃⃛(ω̃)|²、∫P̃dt、∫(1/π)Σ|D̃⃛|²dω̃。"""
    xi = np.linspace(-xi_max, xi_max, n_xi)
    c, s = np.cosh(xi), np.sinh(xi)
    g = np.sqrt(max(e**2 - 1.0, 0.0))
    x, y = e + c, g * s
    r = e * c + 1.0
    t = e * s + xi
    vx, vy = s / r, g * c / r                        # ẋ = (dx/dξ)/(dt/dξ)
    ax, ay = x / r**3, y / r**3                      # r̈ = r/r³（斥力、k/μ = 1）
    rdot = (x * vx + y * vy)
    jx = vx / r**3 - 3 * x * rdot / r**5             # r⃛
    jy = vy / r**3 - 3 * y * rdot / r**5
    # Q⃛_{αβ} = 3(r̈_α ṙ_β + ṙ_α r̈_β) + r⃛_α r_β + r_α r⃛_β、D⃛ = (1/2)(3Q⃛ − δ Tr Q⃛)
    Qxx = 3 * (2 * ax * vx) + 2 * jx * x
    Qyy = 3 * (2 * ay * vy) + 2 * jy * y
    Qxy = 3 * (ax * vy + vx * ay) + jx * y + x * jy
    tr = Qxx + Qyy
    Dxx, Dyy, Dxy = 0.5 * (3 * Qxx - tr), 0.5 * (3 * Qyy - tr), 1.5 * Qxy
    Dzz = -0.5 * tr                                  # z 成分：D_zz = (1/2)(0 − Tr)
    P_int = simpson((Dxx**2 + Dyy**2 + Dzz**2 + 2 * Dxy**2) * r, x=xi)   # ∫ Σ D⃛² dt、dt = r dξ
    if omegas is None:
        omegas = np.logspace(-4, 1.7, 400)
    S = np.empty_like(omegas)
    for i, w in enumerate(omegas):
        ph = np.exp(1j * w * t) * r
        F = [simpson(D * ph, x=xi) for D in (Dxx, Dyy, Dzz, Dxy)]
        S[i] = abs(F[0]) ** 2 + abs(F[1]) ** 2 + abs(F[2]) ** 2 + 2 * abs(F[3]) ** 2
    parseval = simpson(S, x=omegas) / np.pi
    return omegas, S, P_int, parseval


L = []
L.append("# 同符号の対の散乱：断面積と四重極制動放射のスペクトル（in–out の関係量のみ）")
L.append("")
L.append("条件：相対運動エネルギー E = |E(5p₃/₂ F=1)| = %.6f eV（水素の初期準位）。スペクトルの衝突径数は L = ħ√(ℓ(ℓ+1))、ℓ = %d。"
         % (E_REL, ELL))
L.append("")
for pair in ("ee", "pp"):
    p = T.PAIRS[pair]
    c = T.pair_couplings(p)
    mu = c["mu"]                                   # m_e 単位
    q2 = p["q_a"] * p["q_b"] / 9.0                 # e² 単位（正）
    v_c = np.sqrt(2.0 * E_REL / (mu * ME_EV))      # 非相対論
    eta = q2 * ALPHA / v_c
    k_a0 = mu * v_c / ALPHA                        # k [1/a₀]
    a_c = eta / k_a0                               # 距離の径数 a = q_aq_b/(μv²) [a₀]
    b = L_HBAR / k_a0                              # 衝突径数 [a₀]
    e = np.sqrt(1.0 + (b / a_c) ** 2)
    sigma0 = np.imag(loggamma(1 + 1j * eta))
    L.append("## %s（%s）" % (pair, p["name"]))
    L.append("")
    L.append("| 量 | 値 |")
    L.append("|---|---|")
    L.append("| 換算質量 μ | %.4f m_e |" % mu)
    L.append("| v/c | %.4e |" % v_c)
    L.append("| Sommerfeld 径数 η = (q_aq_b/ħc)/(v/c) | %.4f（%s） |" % (eta, "η ≫ 1：準古典" if eta > 20 else "η ≈ 1：量子領域（古典でも Born でもない）"))
    L.append("| 波数 k、de Broglie 波長 λ = 2π/k | %.4f /a₀、%.3f a₀ = %.3f nm |" % (k_a0, 2 * np.pi / k_a0, 2 * np.pi / k_a0 * A0_M * 1e9))
    L.append("| 距離の径数 a = q_aq_b/(μv²)、正面衝突の最接近距離 2a | %.3f a₀、%.3f nm |" % (a_c, 2 * a_c * A0_M * 1e9))
    L.append("| 衝突径数 b（L = ħ√2）、離心率 e、最接近距離 a(e+1) | %.3f a₀、%.4f、%.3f nm |" % (b, e, a_c * (e + 1) * A0_M * 1e9))
    L.append("| 磁気（Darwin）項の大きさ v²/c² | %.2e |" % v_c**2)
    L.append("| 重力/クーロン、重力波/電気四重極 4Gm²/q² | %.3e、%.3e |" % (c["g_ratio"], c["w_gw_identical"]))
    L.append("")
    L.append("### 断面積（重心系、同種スピン ½ 非偏極、Mott）。単位 nm²/sr")
    L.append("")
    L.append("| θ | Mott dσ/dΩ | 区別できる粒子（両方を数える、古典） | 比 |")
    L.append("|---|---|---|---|")
    for th_deg in (10, 20, 30, 45, 60, 75, 90, 105, 120, 135, 150):
        th = np.radians(th_deg)
        m_, r_ = mott(th, eta, k_a0)
        f = (A0_M * 1e9) ** 2
        L.append("| %3d° | %.4e | %.4e | %.4f |" % (th_deg, m_ * f, r_ * f, m_ / r_))
    L.append("")
    L.append("θ = 90° の比 1/2 はフェルミオンの反対称化。干渉項 cos(η ln tan²(θ/2)) は η が大きいほど細かく振れる。断面積は前方で発散（クーロンの長距離性）し、全断面積は定義されない。")
    L.append("")
    # スペクトル
    om, S, P_int, pars = quadrupole_spectrum(e)
    om0, S0, P0, pars0 = quadrupole_spectrum(1.0)
    pref_dEdw = (q2 * 27.211386245988 / a_c) * (a_c * A0_M / (v_c * X.C_SI)) * v_c**5 / (180 * np.pi)    # [eV·s]（q²/a × a/v × (v/c)⁵ /180π）
    unit_w = v_c * X.C_SI / (a_c * A0_M)           # ω の単位 v/a [1/s]
    dE_total = (q2 * 27.211386245988 / a_c) * v_c**5 / (180 * np.pi) * np.pi * pars          # ∫dE/dω dω = (q²/a)(v/c)⁵/(180π)·∫Σ|D̃|²dω̃
    dE_total_time = (q2 * 27.211386245988 / a_c) * v_c**5 / 180 * P_int
    dE_head = (q2 * 27.211386245988 / a_c) * v_c**5 / (180 * np.pi) * np.pi * pars0
    L.append("### 四重極制動放射のスペクトル（古典、b = L/(μv)。双極子は同種粒子で恒等的に 0）")
    L.append("")
    L.append("| ħω [eV] | ω/(v/a) | dE/dω [eV·s] | dN/dω = (dE/dω)/(ħω) [s] |")
    L.append("|---|---|---|---|")
    for wt in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
        i = int(np.argmin(abs(om - wt)))
        hw = HBAR_EVS * om[i] * unit_w
        dEdw = pref_dEdw * S[i]
        L.append("| %.3e | %.2f | %.3e | %.3e |" % (hw, om[i], dEdw, dEdw / hw))
    L.append("")
    L.append("放射エネルギーの合計（b = L/μv）：周波数積分 %.3e eV、時間積分（Parseval 検査）%.3e eV、比 %.6f。正面衝突（b = 0）：%.3e eV。"
             % (dE_total, dE_total_time, dE_total / dE_total_time, dE_head))
    L.append("運動エネルギー E = %.3f eV に対する比 %.1e。光子一個のエネルギーを ħω ~ ħv/a = %.2e eV とすると一回の衝突で出る光子数の期待値は %.1e 個。"
             "重力波はその %.1e 倍（4Gm²/q²）。つまり一回の反跳からは何も出ず、in と out の関係しか記録に残らない。"
             % (E_REL, dE_total / E_REL, HBAR_EVS * unit_w, dE_total / (HBAR_EVS * unit_w), c["w_gw_identical"]))
    ratio_q = HBAR_EVS * unit_w / E_REL
    L.append("スペクトルは ω → 0 で一定（放出前後の D̈ の差）、ω ≳ v/a で指数的に落ちる。ħ(v/a) = %.2e eV は E の %.1e 倍で、%s"
             % (HBAR_EVS * unit_w, ratio_q,
                "古典スペクトルの範囲 ħω ≪ E を満たす（η ≫ 1 と整合）。" if ratio_q < 0.1 else
                "古典スペクトルの前提 ħω ≪ E を満たさない。η ≈ 1 の量子領域で、この表は桁の目安（量子の Gaunt 因子は入れていない）。"))
    L.append("")

out = "\n".join(L) + "\n"
(HERE / "scattering_same_sign.md").write_text(out, encoding="utf-8")
print(out)
