#!/usr/bin/env python3
"""八状態の指数表現・固定交換係数 θ₀ = 23π/124 ― v2（割り当てを決定 B に直し、T を状態から読み、放射を乗法更新にする）

v1 からの変更（木原氏の指摘 10/4）
  1. 位相の枠 a1, b1 は体の時計の文字盤（状態、mod 2π）。a3, b3 はその一節点の進み（= F·ΔT。位相の進み率 = エネルギー）。
     軌道角は持たない（ゲージ。力学は一度も読まない）。
  2. 時間は状態から読む：この節点の ΔT（長さの対と ω⁻¹ の閉じた式）、周期内の位置 E（長さの対から）、時計の文字盤 a1, b1。
     経過時間 T の変数は無い（この系の量ではない）。積算も後処理も無い。
  3. エネルギーは時計の率の比から読む：b3/a3 = (m_b/m_a)·F_B/F_A → 𝐏² → a = κ/(2κ/ξ − 𝐏²)（反復）。
     角運動量は長さの増分 |g| = e^{a2} と a、ξ から：Re u = (|g|² − cos²θ₀ − sin²θ₀|u|²)/(2 sinθ₀cosθ₀)、|u|² = 2a/ξ − 1、
     Im u = +√(|u|² − Re u²)、L = 2ωξ·Im u。相対半角の増分 g = cosθ₀ + sinθ₀·u は導出量。
  4. 書き込みはすべて乗法更新（指数への加算）：a0 += a2、a1 += a3、a2 ← ln|g′|（g′ = 2cosθ₀ − 1/g）、
     a3 ← a3·[F_A′ΔT′/(F_AΔT)]。B も同じ規則で、A から上書きしない。
  5. 放射：エネルギーの損失は時計の率の枠へ乗法で（a3, b3 に因子）、角運動量の損失は長さの増分の枠へ乗法で（|g| に因子）。
     関係への加算と再構成（v1 の state_from_relations）は廃止。因子は掛けない（放出が消える条件は源の定義の問題）。
  6. 1PN〜（T4〜T22）：軌道角が無いので近点移動は現れない（ゲージ）。残るのは動径の振動数 ω の補正で、
     ε_ξξ·ε_ππ（辞書 §2）の Kepler 値との比を ω² に掛ける（係数。摂動）。
割り算・平方根は Newton 反復。exp/log は指数⇄振幅の境界だけ。cos/sin は θ₀ の定数。節点幅・窓・分母の入力はない。
"""
import argparse, csv, math, os

ALPHA = 7.2973525643e-3
M_A = 7.0e-4
RATIO = 1836.15267343
M_B = M_A * RATIO
MTOT = M_A + M_B
MU = M_A * M_B / MTOT
NU = MU / MTOT
GM = M_A / MU
KC = 2 * ALPHA * GM
C0M1 = 9.0 / (M_A * M_B)
KG = KC / C0M1
KAPPA = KC + KG
XI1 = 274.192027
CA2 = M_B / MTOT          # c_A²
CB2 = M_A / MTOT
THETA0 = 23 * math.pi / 124
C0, S0 = math.cos(THETA0), math.sin(THETA0)
C20, S20 = math.cos(2 * THETA0), math.sin(2 * THETA0)
TWO_PI = 2 * math.pi

MONO = [
    (0, 2, 0, 0, 0, 0.5), (0, 4, 0, 0, 0, (3 * NU - 1) / 8), (0, 6, 0, 0, 0, (1 - 5 * NU + 5 * NU**2) / 16),
    (1, 0, 0, 0, 0, -KC - KG), (1, 2, 0, 0, 0, -NU * (KC + KG) - 1.5 * KG), (1, 4, 0, 0, 0, KC * NU * (1 - 2 * NU) / 2),
    (2, 0, 0, 0, 0, KG * (3 * KC + KG) / 2), (2, 0, 2, 0, 0, 2 * GM**2), (2, 2, 0, 0, 0, KC**2 * NU / 4),
    (2, 2, 2, 0, 0, GM**2 * (3 * NU - 1)), (2, 4, 2, 0, 0, 3 * GM**2 * (1 - 5 * NU + 5 * NU**2) / 4),
    (3, 0, 0, 0, 0, -(KC + KG) * (GM**2 - 2 * GM + 2) / 2 + KC**3 * NU / 4), (3, 0, 0, 1, 1, -(GM - 1) * (KC + KG)),
    (3, 0, 1, 0, 1, KC * (GM**2 - 1) + KG * (3 * GM**2 - 2 * GM - 1)), (3, 0, 1, 1, 0, KC * (2 * GM - 1) + KG * (4 * GM - 1)),
    (3, 0, 2, 0, 0, -2 * GM**2 * (NU * (KC + KG) + 3 * KG)), (3, 2, 2, 0, 0, GM**2 * KC * NU * (3 - 4 * NU)),
    (4, 0, 2, 0, 0, GM**2 * KC**2 * NU), (4, 0, 4, 0, 0, GM**4 * (6 * NU - 2)), (4, 2, 4, 0, 0, GM**4 * (3 - 15 * NU + 15 * NU**2)),
    (5, 0, 4, 0, 0, 2 * GM**4 * KC * NU * (2 - 3 * NU)), (6, 0, 6, 0, 0, GM**6 * (4 - 20 * NU + 20 * NU**2)),
]
KEPLER_KEYS = {(0, 2, 0), (1, 0, 0), (2, 0, 2)}


# ---- 逆数・平方根（Newton 反復、種はビット操作。割り算なし）
def recip(x, iters=7):
    m, e = math.frexp(x)
    y = math.ldexp(1.0, -e) * (2.914 - 2.0 * abs(m)) * (1 if m > 0 else -1)
    for _ in range(iters):
        y = y * (2.0 - x * y)
    return y


def sqrt_newton(x, iters=8):
    if x <= 0.0:
        return 0.0
    m, e = math.frexp(x)
    if e % 2:
        m *= 2.0; e -= 1
    y = math.ldexp(0.5 + 0.5 * m, e // 2)
    for _ in range(iters):
        y = 0.5 * (y + x * recip(y))
    return y


def _pw(x, n):
    return 1.0 if n == 0 else x**n


def eps_terms(xi, pi, m, sa, sb, exclude_kepler=False):
    """ε とその二階偏微分 ε_ξξ、ε_ππ（単項式の積和。ξ の負冪は e^{−p ln ξ} で作る）"""
    lnxi = math.log(xi)
    e = exx = epp = 0.0
    for p, a, b, ca, cb, c in MONO:
        if exclude_kepler and (p, a, b) in KEPLER_KEYS:
            continue
        s = c * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb)
        xp = math.exp(-p * lnxi)
        e += s * _pw(pi, a) * xp
        if p > 0:
            exx += s * _pw(pi, a) * p * (p + 1) * math.exp(-(p + 2) * lnxi)
        if a >= 2:
            epp += s * a * (a - 1) * _pw(pi, a - 2) * xp
    return e, exx, epp


def clock_factors(xi, P2):
    rxi = recip(xi)
    va2 = P2 * recip(GM**2)
    Ua = (KG / GM) * rxi
    Fa = 1 - 0.5 * va2 - Ua - va2 * va2 / 8
    vb2 = ((GM - 1) / GM)**2 * P2
    Ub = KG * (GM - 1) / GM * rxi
    Fb = 1 - 0.5 * vb2 - Ub
    return Fa, Fb


def P2_from_clock_ratio(xi, rho, iters=8):
    """時計の率の比 ρ = b3/a3 = RATIO·F_B/F_A から 𝐏² を読む（Newton 反復。種は円軌道 κ/ξ）"""
    target = rho * recip(RATIO)
    P2 = KAPPA * recip(xi)
    for _ in range(iters):
        Fa, Fb = clock_factors(xi, P2)
        rFa = recip(Fa)
        h = Fb * rFa - target
        dFa = -0.5 * recip(GM**2) - P2 * recip(4 * GM**4)
        dFb = -0.5 * ((GM - 1) / GM)**2
        dh = (dFb * Fa - Fb * dFa) * rFa * rFa
        P2 = P2 - h * recip(dh)
    return P2


def wrap_pi(x):
    return (x + math.pi) % TWO_PI - math.pi


def readout(a, b, sa, sb, em):
    """八状態 → 関係（ξ、a_orb、ω、u、g、L、e、E、ΔT、F、x3）。すべて現在の状態から。"""
    rA = math.exp(2 * a[0]); rB = math.exp(2 * b[0])
    xi = rA + rB                                     # 重心の両側：ξ = r_A + r_B
    absg = math.exp(0.5 * (a[2] + b[2]))             # 長さの増分（両体で同じはず。幾何平均で読む）
    rho = b[3] * recip(a[3])                         # 時計の率の比
    P2 = P2_from_clock_ratio(xi, rho)                # 𝐏²（反復）
    inv_a = (2 * KAPPA * recip(xi) - P2) * recip(KAPPA)   # 1/a = (2κ/ξ − 𝐏²)/κ
    a_orb = recip(inv_a)
    om2 = 0.25 * KAPPA * inv_a                       # ω² = κ/(4a)
    # 1PN〜：動径の振動数の補正 corr = ω_r²(全項)/ω_r²(Kepler)（ε_ξξ ε_ππ、円軌道近傍の係数）。
    # 交換の運動学（ω、u、L、e の関係）は Kepler のまま保ち、補正は一節点の座標時間 ΔT に 1/√corr として入れる
    # （動径の周期が 1PN で変わる分だけ、一節点にかかる座標時間が変わる）。時計 x3 = F·ΔT がそれを受ける。
    corr = 1.0
    if em == "full":
        u2_tmp = 2 * a_orb * recip(xi) - 1
        L_tmp = 2 * sqrt_newton(om2) * xi * sqrt_newton(max(u2_tmp, 0.0))
        m_tmp = L_tmp * recip(2 * GM)
        _, exx_f, epp_f = eps_terms(a_orb, 0.0, m_tmp, sa, sb)
        exx_k = 12 * GM**2 * m_tmp**2 * math.exp(-4 * math.log(a_orb)) - 2 * KAPPA * math.exp(-3 * math.log(a_orb))
        corr = (exx_f * epp_f) * recip(exx_k)
    om = sqrt_newton(om2)
    u2 = 2 * a_orb * recip(xi) - 1                   # |u|² = 2a/ξ − 1
    Reu = (absg * absg - C0 * C0 - S0 * S0 * u2) * recip(2 * S0 * C0)
    Imu2 = u2 - Reu * Reu
    Imu = sqrt_newton(Imu2) if Imu2 > 0 else 0.0     # 順行
    L = 2 * om * xi * Imu
    e2 = 1 - L * L * recip(KAPPA * a_orb)            # e² = 1 − L²/(κa)
    ecc = sqrt_newton(e2) if e2 > 0 else 0.0
    cosE = (1 - xi * inv_a) * recip(ecc) if ecc > 1e-12 else 1.0
    cosE = max(-1.0, min(1.0, cosE))
    E = math.acos(cosE) if Reu >= 0 else TWO_PI - math.acos(cosE)   # 周期内の位置（長さの対から）
    n_mean = 2 * om * inv_a                          # 2ω/a = √(κ/a³)
    t_cyc = (E - ecc * math.sin(E)) * recip(n_mean)  # 周期内の時刻（状態から）
    dT = (xi * (THETA0 / 2 + S20 / 4) + u2 * xi * (THETA0 / 2 - S20 / 4) + xi * Reu * (1 - C20) / 2) * recip(om) * recip(sqrt_newton(corr))
    Fa, Fb = clock_factors(xi, P2)
    gr, gi = C0 + S0 * Reu, S0 * Imu                 # g = cosθ₀ + sinθ₀ u
    return dict(xi=xi, rA=rA, rB=rB, a_orb=a_orb, om=om, Reu=Reu, Imu=Imu, u2=u2, L=L, ecc=ecc, E=E, t_cyc=t_cyc,
                dT=dT, Fa=Fa, Fb=Fb, P2=P2, gr=gr, gi=gi, absg=absg, corr=corr, pi=2 * om * Reu * xi * recip(xi),
                epsK=-0.5 * KAPPA * inv_a)


def run(delta, f, nodes, sa, sb, em, out_csv):
    # ---- 初期化（Bohr 型：ξ₁、近点、m = 1 + δ）。両体は重心の両側。時計の文字盤は 0 から
    m0 = 1.0 + delta
    vt = 2 * GM * m0 / XI1
    L0 = XI1 * vt
    epsK0 = 0.5 * vt * vt - KAPPA / XI1
    a_orb0 = -KAPPA / (2 * epsK0)
    om0 = math.sqrt(abs(epsK0) / 2)
    u0 = L0 / (2 * om0 * XI1)                        # 近点：Re u = 0、Im u = L/(2ωξ)
    g0 = complex(C0, S0 * u0)
    dT0 = (XI1 * (THETA0 / 2 + S20 / 4) + u0 * u0 * XI1 * (THETA0 / 2 - S20 / 4)) / om0
    Fa0, Fb0 = clock_factors(XI1, vt * vt)
    a = [0.5 * math.log(CA2 * XI1), 0.0, math.log(abs(g0)), Fa0 * dT0]
    b = [0.5 * math.log(CB2 * XI1), 0.0, math.log(abs(g0)), RATIO * Fb0 * dT0]

    fields = ["node", "xi", "a_orb", "ecc", "L", "eps_K", "omega", "dtau", "dT_node", "E_cycle", "t_cycle",
              "a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3", "absg", "argg", "Reu", "Imu", "corr_1PN",
              "F_a", "F_b", "x3a", "x3b", "d_argz", "Delta_a", "Delta_b", "Delta_ab", "s_dir",
              "dE_node", "dL_node", "sin2_theta_rad"]
    with open(out_csv, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(fields)
        for k in range(nodes + 1):
            # ---- 読出し（状態だけから）
            r = readout(a, b, sa, sb, em)
            argg = math.atan2(r["gi"], r["gr"])
            d_argz = 2 * argg
            Da = wrap_pi(d_argz - a[3]); Db = wrap_pi(d_argz - b[3]); Dab = wrap_pi(b[3] - a[3])
            s_dir = (math.sin(Da / 2)**2 + math.sin(Db / 2)**2 + math.sin(Dab / 2)**2) / 3.0
            rxi = recip(r["xi"])
            dE = -f * (2.0 / 3.0) * KC * KAPPA**2 * rxi**4 * r["dT"]
            dL = -f * (2.0 / 3.0) * KC * KAPPA * r["L"] * rxi**3 * r["dT"]
            sin2_rad = -dE * recip(abs(r["epsK"]))
            wr.writerow([k, f"{r['xi']:.15e}", f"{r['a_orb']:.15e}", f"{r['ecc']:.12e}", f"{r['L']:.15e}", f"{r['epsK']:.15e}",
                         f"{r['om']:.15e}", f"{THETA0 * recip(r['om']):.12e}", f"{r['dT']:.12e}", f"{r['E']:.12e}", f"{r['t_cyc']:.12e}",
                         f"{a[0]:.15e}", f"{a[1]:.12e}", f"{a[2]:.15e}", f"{a[3]:.15e}", f"{b[0]:.15e}", f"{b[1]:.12e}", f"{b[2]:.15e}", f"{b[3]:.15e}",
                         f"{r['absg']:.15e}", f"{argg:.12e}", f"{r['Reu']:.12e}", f"{r['Imu']:.12e}", f"{r['corr']:.12e}",
                         f"{r['Fa']:.15e}", f"{r['Fb']:.15e}", f"{a[3]:.12e}", f"{b[3]:.12e}", f"{d_argz:.9e}",
                         f"{Da:.9e}", f"{Db:.9e}", f"{Dab:.9e}", f"{s_dir:.9e}", f"{dE:.6e}", f"{dL:.6e}", f"{sin2_rad:.6e}"])
            if k == nodes:
                break
            # ---- 書き込み（乗法更新 = 指数への加算）
            # 1. 現在値：長さと時計の文字盤が自身の増分で進む
            a[0] += a[2]; b[0] += b[2]
            a[1] = (a[1] + a[3]) % TWO_PI; b[1] = (b[1] + b[3]) % TWO_PI
            # 2. 長さの増分：g′ = 2cosθ₀ − 1/g（複素の積和。1/g は |g|⁻¹ e^{−i arg g}）
            gr, gi = r["gr"], r["gi"]
            inv_g2 = recip(gr * gr + gi * gi)
            gpr = 2 * C0 - gr * inv_g2; gpi = gi * inv_g2
            absg_new2 = gpr * gpr + gpi * gpi
            # 3. 放射（片側の交換）：エネルギーは時計の率の枠へ、角運動量は長さの増分の枠へ、どちらも因子で
            a_new = r["a_orb"]; L_new = r["L"]
            if f > 0:
                inv_a_new = recip(r["a_orb"]) - 2 * dE * recip(KAPPA)          # κ/a′ = κ/a − 2ΔE
                a_new = recip(inv_a_new)
                L_new = r["L"] * (1 + dL * recip(r["L"]))                        # L′ = L·(1 + ΔL/L)
            # 次の節点の配位（長さは進んだ後の値）から、次の増分と時計の率を作る
            xi_new = math.exp(2 * a[0]) + math.exp(2 * b[0])
            om_new = sqrt_newton(0.25 * KAPPA * recip(a_new))
            u2_new = 2 * a_new * recip(xi_new) - 1
            Reu_n = (gpr - C0) * recip(S0)                                      # 交換後の Re u（呼吸の向きの符号もここから）
            if f > 0:
                # 放射で位置 ξ は変わらず速度だけ変わる：|u′|² = 2a′/ξ′ − 1、Im u′ = L′/(2ω′ξ′)、Re u′ = ±√(|u′|² − Im u′²)（符号は交換後の向き）
                Imu_new = L_new * recip(2 * om_new * xi_new)
                Reu_n = math.copysign(sqrt_newton(max(u2_new - Imu_new * Imu_new, 0.0)), Reu_n)
                absg_new2 = (C0 + S0 * Reu_n)**2 + (S0 * Imu_new)**2           # 長さの増分 |g| への因子（a2 ← ln|g′|）
            lng_new = 0.5 * math.log(absg_new2)
            a[2] = lng_new; b[2] = lng_new
            # 4. 時計の率：次の節点の ΔT′ と F′ から、現在の値への因子で
            P2_new = 2 * KAPPA * recip(xi_new) - KAPPA * recip(a_new)
            corr_new = r["corr"]                                                 # 1PN の係数（次節点の読出しで更新される。ここでは現在値）
            dT_new = (xi_new * (THETA0 / 2 + S20 / 4) + u2_new * xi_new * (THETA0 / 2 - S20 / 4) + xi_new * Reu_n * (1 - C20) / 2) * recip(om_new) * recip(sqrt_newton(corr_new))
            Fa_n, Fb_n = clock_factors(xi_new, P2_new)
            a[3] = a[3] * (Fa_n * dT_new) * recip(r["Fa"] * r["dT"])
            b[3] = b[3] * (Fb_n * dT_new) * recip(r["Fb"] * r["dT"])
    return dict(theta0=THETA0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--f", type=float, default=0.0)
    ap.add_argument("--nodes", type=int, default=1240)
    ap.add_argument("--sa", type=int, default=1)
    ap.add_argument("--sb", type=int, default=-1)
    ap.add_argument("--em", choices=["kepler", "full"], default="full")
    ap.add_argument("--out", required=True)
    x = ap.parse_args()
    os.makedirs(os.path.dirname(x.out) or ".", exist_ok=True)
    run(x.delta, x.f, x.nodes, x.sa, x.sb, x.em, x.out)
    print(f"delta={x.delta} f={x.f} nodes={x.nodes} em={x.em} theta0=23pi/124 R={S0*S0:.6f} done")


if __name__ == "__main__":
    main()
