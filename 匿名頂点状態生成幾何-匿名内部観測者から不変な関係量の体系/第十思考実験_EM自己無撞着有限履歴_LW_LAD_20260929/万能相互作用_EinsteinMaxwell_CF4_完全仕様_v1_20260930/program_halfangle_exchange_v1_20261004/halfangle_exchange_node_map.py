#!/usr/bin/env python3
"""半角位相子の交換としての二体 E–M（固定点 4 の実装）

状態：p, q ∈ ℂ（半角の対象 w = √2 p とその率 w′ = √2 iω q）、二つの時計の位相 x1a, x1b、座標時間 t、節点時間 s。
節点：
  0. 読み出し：ξ = 2|p|²、ż = 2 w w′/ξ、L = 2 Im(w̄ w′)、ε_K = ½|ż|² − κ/ξ、ω = √(|ε_K|/2)、a、e、π、m、F_a、F_b、ε_full
  1. 一節点の座標時間 Δt = ∫|w|² ds（厳密）
  2. 交換 U_R(θ)、θ = ωΔs：p′ = cosθ p + i sinθ q、q′ = i sinθ p + cosθ q   ← S1。R = sin²θ がノルムを配る
  3. （--em full）近点移動：共通位相 Δα = ½·(∂⟨δε⟩/∂L)·Δt（δε = ε_full − ε_K の軌道平均。一次の摂動）
  4. 時計：x1a += F_a Δt、x1b += (m_b/m_a) F_b Δt
  5. 放射（f）：背景への片側の交換。ΔE = −(2/3)κ_C κ²/ξ⁴ Δt、ΔL = −(2/3)κ_C κ L/ξ³ Δt（双極子 z = w²、Kepler の z̈ = −κz/ξ³）。
     N_A + N_B（長半径）と N_A − N_B（L/2ω）を減らす。位相は保つ。--lock a なら源の因子 s_dir を掛ける（既定は掛けない）
  6. 読み出し：三位相の進み（x3a, x3b, 2Δarg w）の差と s_dir（読み a）、全配列の初期への戻り（読み b）
E–M は係数だけ。走行の合否は書かない。窓・分母の入力はない。節点幅は v7 の Δt = 10153（ξ₁ での値）を Δs = Δt/ξ₁ として継承。
"""
import argparse, cmath, csv, math, os

# ---- v7 の定数
ALPHA = 7.2973525643e-3
M_A = 7.0e-4
RATIO = 1836.15267343
M_B = M_A * RATIO
MTOT = M_A + M_B
MU = M_A * M_B / MTOT
NU = MU / MTOT
G = M_A / MU
KC = 2 * ALPHA * G
C0M1 = 9.0 / (M_A * M_B)
KG = KC / C0M1
KAPPA = KC + KG
XI1 = 274.192027
DTH0 = 10153.0

# ---- ε の 22 単項式（v7 §1.2）。T1〜T3 は Kepler 部（ε_K）として別に扱う
MONO = [
    (0, 2, 0, 0, 0, 0.5),
    (0, 4, 0, 0, 0, (3 * NU - 1) / 8),
    (0, 6, 0, 0, 0, (1 - 5 * NU + 5 * NU**2) / 16),
    (1, 0, 0, 0, 0, -KC - KG),
    (1, 2, 0, 0, 0, -NU * (KC + KG) - 1.5 * KG),
    (1, 4, 0, 0, 0, KC * NU * (1 - 2 * NU) / 2),
    (2, 0, 0, 0, 0, KG * (3 * KC + KG) / 2),
    (2, 0, 2, 0, 0, 2 * G**2),
    (2, 2, 0, 0, 0, KC**2 * NU / 4),
    (2, 2, 2, 0, 0, G**2 * (3 * NU - 1)),
    (2, 4, 2, 0, 0, 3 * G**2 * (1 - 5 * NU + 5 * NU**2) / 4),
    (3, 0, 0, 0, 0, -(KC + KG) * (G**2 - 2 * G + 2) / 2 + KC**3 * NU / 4),
    (3, 0, 0, 1, 1, -(G - 1) * (KC + KG)),
    (3, 0, 1, 0, 1, KC * (G**2 - 1) + KG * (3 * G**2 - 2 * G - 1)),
    (3, 0, 1, 1, 0, KC * (2 * G - 1) + KG * (4 * G - 1)),
    (3, 0, 2, 0, 0, -2 * G**2 * (NU * (KC + KG) + 3 * KG)),
    (3, 2, 2, 0, 0, G**2 * KC * NU * (3 - 4 * NU)),
    (4, 0, 2, 0, 0, G**2 * KC**2 * NU),
    (4, 0, 4, 0, 0, G**4 * (6 * NU - 2)),
    (4, 2, 4, 0, 0, G**4 * (3 - 15 * NU + 15 * NU**2)),
    (5, 0, 4, 0, 0, 2 * G**4 * KC * NU * (2 - 3 * NU)),
    (6, 0, 6, 0, 0, G**6 * (4 - 20 * NU + 20 * NU**2)),
]


def _pw(x, n):
    return 1.0 if n == 0 else x**n


def eps_full(xi, pi, m, sa, sb):
    return sum(c * _pw(pi, a) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p) for p, a, b, ca, cb, c in MONO)


def eps_kepler(xi, pi, m):
    return 0.5 * (pi**2 + (2 * G * m / xi)**2) - KAPPA / xi


def clock_factors(xi, P2):
    va2 = P2 / G**2
    Ua = (KG / G) / xi
    Fa = 1 - 0.5 * va2 - Ua - va2**2 / 8
    vb2 = ((G - 1) / G)**2 * P2
    Ub = KG * (G - 1) / (G * xi)
    Fb = 1 - 0.5 * vb2 - Ub
    return Fa, Fb


def wrap_pi(x):
    return (x + math.pi) % (2 * math.pi) - math.pi


def orbit_avg_delta_eps(a_orb, e, sa, sb, nq=64):
    """δε = ε_full − ε_K の Kepler 軌道平均（離心近点角 E で求積、重み dt/dE ∝ 1 − e cos E）。
    固定 a で e を動かすので、L = √(κa(1−e²))、m = L/(2g) も e に従って動かす（軌道上で整合）。"""
    e = abs(e)
    L = math.sqrt(KAPPA * a_orb * (1 - e * e))
    m = L / (2 * G)
    tot = 0.0
    for j in range(nq):
        E = 2 * math.pi * (j + 0.5) / nq
        r = a_orb * (1 - e * math.cos(E))
        rdot = math.sqrt(KAPPA / a_orb) * e * math.sin(E) / (1 - e * math.cos(E))
        tot += (eps_full(r, rdot, m, sa, sb) - eps_kepler(r, rdot, m)) * (1 - e * math.cos(E))
    return tot / nq


def precession_rate(a_orb, e, L, sa, sb):
    """近点移動の率 ∂⟨δε⟩/∂L（座標時間あたり、Delaunay の固定 a）。L ↔ e は e² = 1 − L²/(κa)。
    ⟨δε⟩ は e の偶関数なので (1/e)∂/∂e を差分で取る（e→0 で有限）。"""
    h = 1e-3
    f1 = orbit_avg_delta_eps(a_orb, e + h, sa, sb)
    f2 = orbit_avg_delta_eps(a_orb, abs(e - h), sa, sb)
    if e > h:
        one_over_e_dde = (f1 - f2) / (2 * h * e)
    else:
        # e ≪ h：⟨δε⟩(e+h) − ⟨δε⟩(h−e) ≈ 2e·⟨δε⟩′(h) ⇒ (1/e)∂/∂e ≈ ⟨δε⟩′(h)/h（偶関数の極限）
        fp = orbit_avg_delta_eps(a_orb, h + h * 0.5, sa, sb)
        fm = orbit_avg_delta_eps(a_orb, h - h * 0.5, sa, sb)
        one_over_e_dde = ((fp - fm) / h) / h
    de_dL = -L / (KAPPA * a_orb)      # 符号付き：de/dL = −L/(κ a e) の e を上で割っている
    return one_over_e_dde * de_dL


def run(delta, f, nodes, sa, sb, em, lock, out_csv, dth0=DTH0, tol_rec=1e-3):
    # ---- 初期化（v7 §8：近点 ξ₁、π = 0、m = 1 + δ）。半角の対象 w = √ξ₁、w′ = i v_t √ξ₁ / 2
    m0 = 1.0 + delta
    vt = 2 * G * m0 / XI1
    w = complex(math.sqrt(XI1), 0.0)
    wp = 1j * vt * math.sqrt(XI1) / 2
    epsK = eps_kepler(XI1, 0.0, m0)
    om = math.sqrt(abs(epsK) / 2)
    p = w / math.sqrt(2)
    q = wp / (math.sqrt(2) * 1j * om)
    ds = dth0 / XI1                      # 継承した節点幅（節点時間）
    t = 0.0; s = 0.0; x1a = 0.0; x1b = 0.0
    init_cfg = None
    recur = []

    fields = ["node", "t", "s", "xi", "pi", "m", "L", "eps_K", "eps_full", "omega", "a_orb", "ecc",
              "N_p", "N_q", "Re_pq", "Im_pq", "N_A", "N_B", "theta", "R_sin2theta", "dt_node",
              "F_a", "F_b", "x3a", "x3b", "x3b_over_x3a", "d_argz", "Delta_a", "Delta_b", "Delta_ab", "s_dir",
              "prec_rate", "dalpha_node", "dE_node", "dL_node", "sin2_theta_rad", "arg_z", "arg_pq", "x1a", "x1b"]
    with open(out_csv, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(fields)
        for k in range(nodes + 1):
            # ---- 0. 読み出し（関係から）
            w = math.sqrt(2) * p
            Np, Nq = abs(p)**2, abs(q)**2
            pq = p.conjugate() * q
            xi = 2 * Np
            # ω は前節点の値で w′ を復元し、ε_K を読み直して ω を更新（Kepler では不変）
            wp = math.sqrt(2) * 1j * om * q
            zdot = 2 * w * wp / xi
            z = w * w
            L = 2 * (w.conjugate() * wp).imag
            P2 = abs(zdot)**2
            pi_r = (z.conjugate() * zdot).real / xi
            m = L / (2 * G)
            epsK = 0.5 * P2 - KAPPA / xi
            om = math.sqrt(abs(epsK) / 2)
            q = wp / (math.sqrt(2) * 1j * om)          # ω を更新したので q を整合（Kepler では無変化）
            Nq = abs(q)**2; pq = p.conjugate() * q
            alpha = (p + q) / math.sqrt(2); beta = (p - q) / math.sqrt(2)
            NA, NB = abs(alpha)**2, abs(beta)**2
            a_orb = NA + NB
            ecc = 2 * math.sqrt(NA * NB) / a_orb if a_orb > 0 else 0.0
            Fa, Fb = clock_factors(xi, P2)
            ef = eps_full(xi, pi_r, m, sa, sb)
            theta = om * ds
            # ---- 1. 一節点の座標時間（厳密）
            dt = a_orb * ds + 2 * (alpha * beta.conjugate() * (cmath.exp(2j * om * ds) - 1) / (2j * om)).real
            x3a = Fa * dt
            x3b = RATIO * Fb * dt
            # ---- 3. 近点移動（em full）：共通位相。率は軌道平均の一次摂動
            prec = precession_rate(a_orb, ecc, L, sa, sb) if em == "full" else 0.0
            dalpha = 0.5 * prec * dt
            # ---- 5. 放射の量（双極子 z = w²、Kepler の z̈ = −κz/ξ³）
            dE = -f * (2.0 / 3.0) * KC * KAPPA**2 / xi**4 * dt
            dL = -f * (2.0 / 3.0) * KC * KAPPA * L / xi**3 * dt
            # ---- 6. 読み出し（a）：三位相の進みの差。双極子の進み 2Δarg w はこの節点の交換で決まる（先に計算）
            p_n = math.cos(theta) * p + 1j * math.sin(theta) * q
            q_n = 1j * math.sin(theta) * p + math.cos(theta) * q
            w_n = math.sqrt(2) * p_n * cmath.exp(1j * dalpha)
            d_argz = 2 * cmath.phase(w_n / w)          # 双極子 z = w² の一節点の位相の進み（(−2π, 2π]）
            Da = wrap_pi(d_argz - x3a); Db = wrap_pi(d_argz - x3b); Dab = wrap_pi(x3b - x3a)
            s_dir = (math.sin(Da / 2)**2 + math.sin(Db / 2)**2 + math.sin(Dab / 2)**2) / 3.0
            s_eff = s_dir if lock == "a" else 1.0
            dE *= s_eff; dL *= s_eff
            sin2_rad = -dE / abs(epsK) if epsK != 0 else 0.0
            # ---- 読み出し（b）：全配列の初期への戻り（双極子の位相、二つの時計）
            arg_z = cmath.phase(z); arg_pq = cmath.phase(pq) if abs(pq) > 0 else 0.0
            cfg = (arg_z, x1a, x1b)
            if init_cfg is None:
                init_cfg = cfg
            elif k > 0 and all(abs(wrap_pi(c - c0)) < tol_rec * 2 * math.pi for c, c0 in zip(cfg, init_cfg)):
                recur.append((k, t))
            wr.writerow([k, f"{t:.15e}", f"{s:.15e}", f"{xi:.15e}", f"{pi_r:.9e}", f"{m:.12e}", f"{L:.15e}",
                         f"{epsK:.15e}", f"{ef:.12e}", f"{om:.15e}", f"{a_orb:.15e}", f"{ecc:.15e}",
                         f"{Np:.12e}", f"{Nq:.12e}", f"{pq.real:.12e}", f"{pq.imag:.12e}", f"{NA:.12e}", f"{NB:.12e}",
                         f"{theta:.12e}", f"{math.sin(theta)**2:.9e}", f"{dt:.9e}",
                         f"{Fa:.15e}", f"{Fb:.15e}", f"{x3a:.9e}", f"{x3b:.9e}", f"{x3b / x3a:.12e}", f"{d_argz:.9e}",
                         f"{Da:.9e}", f"{Db:.9e}", f"{Dab:.9e}", f"{s_dir:.9e}",
                         f"{prec:.9e}", f"{dalpha:.9e}", f"{dE:.6e}", f"{dL:.6e}", f"{sin2_rad:.6e}",
                         f"{arg_z:.9e}", f"{arg_pq:.9e}", f"{x1a:.9e}", f"{x1b:.9e}"])
            if k == nodes:
                break
            # ---- 2. 交換 U_R(θ)（S1）＋ 3. 共通位相
            p, q = p_n * cmath.exp(1j * dalpha), q_n * cmath.exp(1j * dalpha)
            # ---- 4. 時計・時間
            x1a = (x1a + x3a) % (2 * math.pi)
            x1b = (x1b + x3b) % (2 * math.pi)
            t += dt; s += ds
            # ---- 5. 放射：ノルムだけを背景へ出す（位相は保つ）
            if f > 0 and (dE != 0 or dL != 0):
                epsK_n = epsK + dE
                L_n = L + dL
                a_n = -KAPPA / (2 * epsK_n)
                om_n = math.sqrt(abs(epsK_n) / 2)
                dN = L_n / (2 * om_n)                    # N_A − N_B
                NA_n = (a_n + dN) / 2; NB_n = (a_n - dN) / 2
                alpha = (p + q) / math.sqrt(2); beta = (p - q) / math.sqrt(2)
                alpha = alpha * math.sqrt(max(NA_n, 0.0) / abs(alpha)**2) if abs(alpha) > 0 else alpha
                beta = beta * math.sqrt(max(NB_n, 0.0) / abs(beta)**2) if abs(beta) > 0 else beta
                p = (alpha + beta) / math.sqrt(2); q = (alpha - beta) / math.sqrt(2)
                om = om_n
    return dict(recur=recur[:20], n_recur=len(recur), ds=ds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--f", type=float, default=0.0)
    ap.add_argument("--nodes", type=int, default=2000)
    ap.add_argument("--sa", type=int, default=1)
    ap.add_argument("--sb", type=int, default=-1)
    ap.add_argument("--em", choices=["kepler", "full"], default="full", help="kepler: T1〜T3 のみ / full: 近点移動（一次摂動）も")
    ap.add_argument("--lock", choices=["none", "a"], default="none", help="none: 放出に因子を掛けない（既定） / a: 読み (a) の s_dir を掛ける")
    ap.add_argument("--dth0", type=float, default=DTH0, help="継承した節点の座標時間幅（ξ₁ での値）")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    r = run(a.delta, a.f, a.nodes, a.sa, a.sb, a.em, a.lock, a.out, dth0=a.dth0)
    print(f"delta={a.delta} f={a.f} nodes={a.nodes} em={a.em} lock={a.lock} ds={r['ds']:.6f} "
          f"recur_count={r['n_recur']} first_recur(node,t)={r['recur'][:6]}")


if __name__ == "__main__":
    main()
