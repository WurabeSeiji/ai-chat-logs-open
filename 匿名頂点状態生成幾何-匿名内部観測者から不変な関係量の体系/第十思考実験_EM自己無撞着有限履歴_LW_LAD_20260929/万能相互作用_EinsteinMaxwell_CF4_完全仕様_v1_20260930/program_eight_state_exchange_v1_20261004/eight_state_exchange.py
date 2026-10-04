#!/usr/bin/env python3
"""八状態の指数表現・固定交換係数 θ₀ = 23π/124（固定点 5 の実装）

状態（隠しなし）：
  体 A (a0,a1,a2,a3)、体 B (b0,b1,b2,b3)：指数。W = e^{a0+i a1}（半角の対象）、G = e^{a2+i a3}（生成子）
  dtau：節点の目盛（補助状態 = θ₀·ω⁻¹、ω は状態から読む）
  T   ：座標時間（状態。一節点ごとに ΔT を乗算値として積む）
読出し：関係 z = W²、ξ、L、ε_K、a、e、F_A、F_B、時計の進み x3a・x3b と位相 x1a・x1b、三位相の差、124 節点の戻り。
一節点（読出し・演算・書き込み）：
  1. z1 ← z1 + z2（両体。W ← W·G、乗法更新）
  2. e^{z2'} = 2cosθ₀ − e^{−z2}（増分の書き戻し。係数固定の積和。回転 U_R(θ₀) と同じ写像）
  3. （--em full）T4〜T22：ε の偏微分 × ΔT を関係 (ξ, π, φ) への積和として入れ、(W, G) に写す（一次の摂動）
  4. （--f）放射：ΔE、ΔL（双極子 W²、Kepler の z̈）を (ε_K, L) へ入れ、(W, G) に写す。因子は掛けない（読み (a)(b) は出力）
  5. ω を不変量から読み、dtau = θ₀·ω⁻¹、ΔT = ∫|W|² を閉じた式で、T に積む。時計を読出しとして積算
割り算と平方根は Newton 反復（種はビット操作、割り算なし）。exp/log は指数⇄振幅の読み書きの境界だけ。cos/sin は θ₀ の定数。
E–M は係数だけ。走行の合否は書かない。節点幅・窓・分母の入力はない。
"""
import argparse, cmath, csv, math, os

# ---- 定数（v7）
ALPHA = 7.2973525643e-3
M_A = 7.0e-4
RATIO = 1836.15267343
M_B = M_A * RATIO
MTOT = M_A + M_B
MU = M_A * M_B / MTOT
NU = MU / MTOT
GM = M_A / MU                       # 質量因子 g（v7）。増分 G と区別するため GM
KC = 2 * ALPHA * GM
C0M1 = 9.0 / (M_A * M_B)
KG = KC / C0M1
KAPPA = KC + KG
XI1 = 274.192027
CA = math.sqrt(M_B / MTOT)          # 体 A の振幅係数（重心系）
CB = math.sqrt(M_A / MTOT)
THETA0 = 23 * math.pi / 124          # 固定の交換係数（量子化の拘束）
C0, S0 = math.cos(THETA0), math.sin(THETA0)
C20, S20 = math.cos(2 * THETA0), math.sin(2 * THETA0)

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
KEPLER_KEYS = {(0, 2, 0), (1, 0, 0), (2, 0, 2)}     # T1〜T3


def _pw(x, n):
    return 1.0 if n == 0 else x**n


def eps_full(xi, pi, m, sa, sb):
    return sum(c * _pw(pi, a) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p) for p, a, b, ca, cb, c in MONO)


def d_eps_full(xi, pi, m, sa, sb):
    """(∂ε/∂ξ, ∂ε/∂π, ∂ε/∂m) を単項式から（T1〜T3 を除いた摂動部分だけ）"""
    dxi = dpi = dm = 0.0
    for p, a, b, ca, cb, c in MONO:
        if (p, a, b) in KEPLER_KEYS:
            continue
        s = _pw(sa, ca) * _pw(sb, cb)
        if p > 0:
            dxi += -p * c * _pw(pi, a) * _pw(m, b) * s * xi**(-p - 1)
        if a > 0:
            dpi += c * a * _pw(pi, a - 1) * _pw(m, b) * s * xi**(-p)
        if b > 0:
            dm += c * _pw(pi, a) * b * _pw(m, b - 1) * s * xi**(-p)
    return dxi, dpi, dm


# ---- 割り算・平方根を使わない逆数と平方根（Newton 反復、種はビット操作）
def recip(x, iters=7):
    m, e = math.frexp(x)                 # x = m·2^e, 0.5 ≤ |m| < 1
    y = math.ldexp(1.0, -e) * (2.914 - 2.0 * abs(m)) * (1 if m > 0 else -1)   # 線形の種
    for _ in range(iters):
        y = y * (2.0 - x * y)
    return y


def sqrt_newton(x, iters=8):
    if x <= 0.0:
        return 0.0
    m, e = math.frexp(x)
    if e % 2:
        m *= 2.0; e -= 1
    y = math.ldexp(0.5 + 0.5 * m, e // 2)   # 種
    for _ in range(iters):
        y = 0.5 * (y + x * recip(y))
    return y


def clock_factors(xi, P2):
    rxi = recip(xi)
    va2 = P2 * recip(GM**2)
    Ua = (KG / GM) * rxi
    Fa = 1 - 0.5 * va2 - Ua - va2**2 / 8
    vb2 = ((GM - 1) / GM)**2 * P2
    Ub = KG * (GM - 1) / GM * rxi
    Fb = 1 - 0.5 * vb2 - Ub
    return Fa, Fb


def wrap_pi(x):
    return (x + math.pi) % (2 * math.pi) - math.pi


def state_from_relations(xi, phi, pi, L):
    """関係 (ξ, φ, π, L) から半角の対象 W、ω、増分 G を作る（初期化と、関係へのキックの書き戻しに使う）"""
    W = sqrt_newton(xi) * cmath.exp(1j * phi / 2)
    zdot = (pi + 1j * L * recip(xi)) * cmath.exp(1j * phi)
    epsK = 0.5 * abs(zdot)**2 - KAPPA * recip(xi)
    om = sqrt_newton(0.5 * abs(epsK))
    Wdot = zdot * xi * recip(2 * W.real**2 + 2 * W.imag**2) * W.conjugate()     # ż ξ /(2W) = ż ξ W̄ /(2|W|²)
    V = Wdot * recip(om)
    G = C0 + S0 * V * W.conjugate() * recip(abs(W)**2)   # cosθ₀ + sinθ₀ V/W
    return W, om, G, epsK


def run(delta, f, nodes, sa, sb, em, out_csv, tol_rec=1e-3):
    # ---- 初期化（Bohr 型：ξ₁、近点 π = 0、m = 1 + δ）。両体は重心系で W_A = c_A W、W_B = i c_B W
    m0 = 1.0 + delta
    vt = 2 * GM * m0 / XI1
    W, om, G, epsK = state_from_relations(XI1, 0.0, 0.0, XI1 * vt)
    a = [math.log(CA) + math.log(abs(W)), cmath.phase(W), math.log(abs(G)), cmath.phase(G)]
    b = [math.log(CB) + math.log(abs(W)), cmath.phase(W) + math.pi / 2, math.log(abs(G)), cmath.phase(G)]
    dtau = THETA0 * recip(om)
    T = 0.0
    x1a = x1b = 0.0
    z_init = None; clock_init = None
    rec124 = []

    fields = ["node", "T", "xi", "phi", "pi", "m", "L", "eps_K", "eps_full", "omega", "dtau", "a_orb", "ecc",
              "a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3", "absG", "argG", "R_fixed",
              "dT_node", "F_a", "F_b", "x3a", "x3b", "x1a", "x1b", "d_argz", "Delta_a", "Delta_b", "Delta_ab", "s_dir",
              "dE_node", "dL_node", "sin2_theta_rad", "varpi", "rec124_dz", "rec124_dclock_a", "rec124_dclock_b"]
    with open(out_csv, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(fields)
        for k in range(nodes + 1):
            # ---- 読出し（状態 → 振幅 → 関係）
            WA = cmath.exp(complex(a[0], a[1])); WB = cmath.exp(complex(b[0], b[1]))
            G = cmath.exp(complex(a[2], a[3]))
            z = WA * WA - WB * WB                       # 相対関係 = 二体の波の二乗の差（重心系で W²）
            W = WA * recip(CA)                           # 半角の対象（枝は A から）
            xi = abs(W)**2
            # W′/ω：増分 G は次へ進む増分（W_{n+1} = W_n·G = cosθ₀W_n + sinθ₀V_n）なので V = (G − cosθ₀)·W/sinθ₀
            V = (G - C0) * W * recip(S0)
            a_orb = 0.5 * (abs(W)**2 + abs(V)**2)
            om = sqrt_newton(KAPPA * 0.25 * recip(a_orb))
            L = 2 * om * (W.conjugate() * V).imag
            epsK = -0.5 * KAPPA * recip(a_orb)
            e2 = 1 - L * L * recip(KAPPA * a_orb)
            ecc = sqrt_newton(e2) if e2 > 0 else 0.0
            zdot = 2 * om * W * V * recip(xi)
            P2 = abs(zdot)**2
            phi = cmath.phase(z)
            pi_r = (z.conjugate() * zdot).real * recip(xi)
            m = L * recip(2 * GM)
            ef = eps_full(xi, pi_r, m, sa, sb)
            Fa, Fb = clock_factors(xi, P2)
            dtau = THETA0 * recip(om)
            dT = (abs(W)**2 * (THETA0 / 2 + S20 / 4) + abs(V)**2 * (THETA0 / 2 - S20 / 4) + (W.conjugate() * V).real * (1 - C20) / 2) * recip(om)
            x3a = Fa * dT; x3b = RATIO * Fb * dT
            # 近点方向（離心率ベクトル）
            evec = ((P2 - KAPPA * recip(xi)) * z - (z.conjugate() * zdot).real * zdot) * recip(KAPPA)
            varpi = cmath.phase(evec) if abs(evec) > 1e-300 else 0.0
            # 放射の量（双極子 W²、Kepler の z̈ = −κ z/ξ³）
            rxi = recip(xi)
            dE = -f * (2.0 / 3.0) * KC * KAPPA**2 * rxi**4 * dT
            dL = -f * (2.0 / 3.0) * KC * KAPPA * L * rxi**3 * dT
            sin2_rad = -dE * recip(abs(epsK)) if epsK != 0 else 0.0
            # 読み (a)：双極子の一節点の位相の進み = 2·arg G（この節点の増分）
            d_argz = 2 * a[3]
            Da = wrap_pi(d_argz - x3a); Db = wrap_pi(d_argz - x3b); Dab = wrap_pi(x3b - x3a)
            s_dir = (math.sin(Da / 2)**2 + math.sin(Db / 2)**2 + math.sin(Dab / 2)**2) / 3.0
            # 読み (b)：124 節点ごとの戻り（関係 z と二つの時計）
            r_dz = r_ca = r_cb = float("nan")
            if k % 124 == 0:
                if z_init is None:
                    z_init = z; clock_init = (x1a, x1b)
                else:
                    r_dz = abs(z - z_init) * recip(abs(z_init))
                    r_ca = wrap_pi(x1a - clock_init[0]); r_cb = wrap_pi(x1b - clock_init[1])
                    rec124.append((k, r_dz, r_ca, r_cb))
            wr.writerow([k, f"{T:.15e}", f"{xi:.15e}", f"{phi:.12e}", f"{pi_r:.9e}", f"{m:.12e}", f"{L:.15e}",
                         f"{epsK:.15e}", f"{ef:.12e}", f"{om:.15e}", f"{dtau:.12e}", f"{a_orb:.15e}", f"{ecc:.12e}",
                         f"{a[0]:.12e}", f"{a[1]:.12e}", f"{a[2]:.12e}", f"{a[3]:.12e}", f"{b[0]:.12e}", f"{b[1]:.12e}", f"{b[2]:.12e}", f"{b[3]:.12e}",
                         f"{abs(G):.12e}", f"{cmath.phase(G):.12e}", f"{S0*S0:.9f}",
                         f"{dT:.12e}", f"{Fa:.15e}", f"{Fb:.15e}", f"{x3a:.12e}", f"{x3b:.12e}", f"{x1a:.9e}", f"{x1b:.9e}",
                         f"{d_argz:.9e}", f"{Da:.9e}", f"{Db:.9e}", f"{Dab:.9e}", f"{s_dir:.9e}",
                         f"{dE:.6e}", f"{dL:.6e}", f"{sin2_rad:.6e}", f"{varpi:.9e}", f"{r_dz:.3e}", f"{r_ca:.6e}", f"{r_cb:.6e}"])
            if k == nodes:
                break
            # ---- 書き込み 1・2：現在値の乗法更新と増分の積和（両体）
            for st in (a, b):
                st[0] += st[2]; st[1] = (st[1] + st[3]) % (2 * math.pi)
                g_new = 2 * C0 - cmath.exp(complex(-st[2], -st[3]))
                st[2] = math.log(abs(g_new)); st[3] = cmath.phase(g_new)
            # ---- 書き込み 3・4：関係へのキック（1PN、放射）を (W, G) へ写す
            if em == "full" or f > 0:
                WA = cmath.exp(complex(a[0], a[1])); WB = cmath.exp(complex(b[0], b[1]))
                z = WA * WA - WB * WB; W = WA * recip(CA); xi = abs(W)**2
                G = cmath.exp(complex(a[2], a[3])); V = (G - C0) * W * recip(S0)
                a_orb = 0.5 * (abs(W)**2 + abs(V)**2); om = sqrt_newton(KAPPA * 0.25 * recip(a_orb))
                L = 2 * om * (W.conjugate() * V).imag; zdot = 2 * om * W * V * recip(xi)
                phi = cmath.phase(z); pi_r = (z.conjugate() * zdot).real * recip(xi); m = L * recip(2 * GM)
                epsK = 0.5 * abs(zdot)**2 - KAPPA * recip(xi)
                if em == "full":
                    dxi, dpi, dm = d_eps_full(xi, pi_r, m, sa, sb)
                    xi = xi + dpi * dT                      # Δξ = ∂δε/∂π·ΔT
                    pi_r = pi_r - dxi * dT                  # Δπ = −∂δε/∂ξ·ΔT
                    phi = phi + dm * recip(2 * GM) * dT     # Δφ = ∂δε/∂L·ΔT
                if f > 0:
                    epsK = epsK + dE; L = L + dL
                    v2 = 2 * (epsK + KAPPA * recip(xi))
                    vt2 = L * L * recip(xi * xi)
                    pr2 = v2 - vt2
                    pi_r = math.copysign(sqrt_newton(pr2), pi_r) if pr2 > 0 else 0.0
                W_new, om, G_new, epsK = state_from_relations(xi, phi, pi_r, L)
                # 半角の枝：A の現在の枝に近い方
                if abs(W_new * CA - cmath.exp(complex(a[0], a[1]))) > abs(-W_new * CA - cmath.exp(complex(a[0], a[1]))):
                    W_new = -W_new
                a[0] = math.log(CA) + math.log(abs(W_new)); a[1] = cmath.phase(W_new) % (2 * math.pi)
                b[0] = math.log(CB) + math.log(abs(W_new)); b[1] = (cmath.phase(W_new) + math.pi / 2) % (2 * math.pi)
                a[2] = b[2] = math.log(abs(G_new)); a[3] = b[3] = cmath.phase(G_new)
            # ---- 書き込み 5：時間と時計（読出しの積算）
            T += dT
            x1a = (x1a + x3a) % (2 * math.pi); x1b = (x1b + x3b) % (2 * math.pi)
    return dict(rec124=rec124, dtau0=THETA0 * recip(om))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--f", type=float, default=0.0)
    ap.add_argument("--nodes", type=int, default=1240)
    ap.add_argument("--sa", type=int, default=1)
    ap.add_argument("--sb", type=int, default=-1)
    ap.add_argument("--em", choices=["kepler", "full"], default="full")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    r = run(a.delta, a.f, a.nodes, a.sa, a.sb, a.em, a.out)
    rec = r["rec124"]
    print(f"delta={a.delta} f={a.f} nodes={a.nodes} em={a.em} theta0=23pi/124 R={S0*S0:.6f} "
          f"rec124: count={len(rec)} first={[(k, f'{dz:.1e}', f'{ca:.3f}', f'{cb:.3f}') for k, dz, ca, cb in rec[:4]]}")


if __name__ == "__main__":
    main()
