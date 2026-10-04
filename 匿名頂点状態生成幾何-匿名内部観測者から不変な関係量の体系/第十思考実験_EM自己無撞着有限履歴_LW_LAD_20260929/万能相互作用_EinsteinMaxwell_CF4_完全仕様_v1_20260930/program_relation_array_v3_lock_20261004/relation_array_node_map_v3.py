#!/usr/bin/env python3
"""関係配列写像 v3 ＝ v1 の忠実コピー ＋ 決定済みの引き込み（--lock off|on）

決定済み（チャット 10/1、STATUS §7.1・7.5）を v1 の上にそのまま実装する。窓・分母・刻みの入力はない。
  1. 三つの位相の進み：a の時計 x3a、b の時計 x3b、関係の位相 2θ。読めるのはその差だけ。
     Δa = 2θ − x3a、Δb = 2θ − x3b、Δab = x3b − x3a（mod 2π）。三つが整数関係（≡ 0）なら何も変化が読めない。
  2. 放出の源は位相配位の関数：s = (1/3)Σ sin²(Δ_i/2)。整数関係で s = 0 ⇒ 方向を失い放出が止まる。
     E–M の放射反作用（ALD の P_rad、Δm）はそのまま。背景は B3（等方な吸収体、戻らない）。
  3. 節点の時計 δτ は動的：a の位相が関係の位相と閉じる値 x3a = 2πK + 2θ（Δa ≡ 0）で各節点ごとに決まる。
     K は v1 の節点幅（v7 の 10153 ＝ 1615.86 周）から継承した整数の周数 1616。新しい入力ではない。
  --lock off は v1 と同一（対照走行で照合）。

---- 以下 v1 の説明 ----
関係配列写像 v1（辞書 10/4・固定点 1〜3 の最小実装）

状態は a, b の成分ではなく「関係の配列」とその率（固定点 2 §4, §6）:
  λ   : 対数距離（長さの対）          ξ = e^λ  [a_a 単位]
  Δλ  : λ の一節点の増分（長さの率）   π = ξ Δλ / Δθ_k
  ψ   : 関係の位相（背景に対する向き。mod 2π。配列の回転角の積算）
  m   : 反対称の組の不変量（角運動量 L/ħ）。θ はここから ε を介して読む
  x3a : a の位相の一節点の増分 = 節点の時計（固定点 3）。Δθ_k = x3a / F_a
  x1a, x1b : 二つの時計の位相（mod 2π。帳簿）
節点 = 読み出し → 角 θ, θ_X, η_a, η_b → 面②（呼吸: E–M の正準更新）→ 面①（ψ += 2θ）
       → 時計 → 放射（f: 一節点に出る割合 sin²θ_rad = ΔE/|ε|、Δm = −ΔE/(2g φ̇)）
E–M は係数（ε の 22 単項式, v7 §1.2）だけ。走行の合否は書かない。
"""
import argparse, csv, math, os, sys

# ---- v7 の定数（coefficients_relational_full_list_ja_final_v7_20261003.md §0, §8）
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
XI1 = 274.192027
DTH0 = 10153.0          # v7 の節点幅（状態 (0,1)）。初期の x3a = F_a·DTH0 を与えるだけ

# ---- ε の 22 単項式 (p, a, b, s_a, s_b, c)  v7 §1.2
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


def eps(xi, pi, m, sa, sb):
    return sum(c * _pw(pi, a) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p)
               for p, a, b, ca, cb, c in MONO)


def d_m(xi, pi, m, sa, sb):
    return sum(c * _pw(pi, a) * b * _pw(m, b - 1) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p)
               for p, a, b, ca, cb, c in MONO if b > 0)


def d_xi(xi, pi, m, sa, sb):
    return sum(-p * c * _pw(pi, a) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p - 1)
               for p, a, b, ca, cb, c in MONO if p > 0)


def d_pi(xi, pi, m, sa, sb):
    return sum(c * a * _pw(pi, a - 1) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p)
               for p, a, b, ca, cb, c in MONO if a > 0)


def d_xixi(xi, pi, m, sa, sb):
    return sum(p * (p + 1) * c * _pw(pi, a) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p - 2)
               for p, a, b, ca, cb, c in MONO if p > 0)


def d_pipi(xi, pi, m, sa, sb):
    return sum(c * a * (a - 1) * _pw(pi, a - 2) * _pw(m, b) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p)
               for p, a, b, ca, cb, c in MONO if a >= 2)


def d_mxi(xi, pi, m, sa, sb):
    return sum(-p * c * _pw(pi, a) * b * _pw(m, b - 1) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p - 1)
               for p, a, b, ca, cb, c in MONO if b > 0 and p > 0)


def d_mpi(xi, pi, m, sa, sb):
    return sum(c * a * _pw(pi, a - 1) * b * _pw(m, b - 1) * _pw(sa, ca) * _pw(sb, cb) * xi**(-p)
               for p, a, b, ca, cb, c in MONO if b > 0 and a > 0)


def wrap_pi(x):
    """角を (−π, π] に写す（v3：位相の進みの差を読むため）"""
    return (x + math.pi) % (2 * math.pi) - math.pi


def clock_factors(xi, P2):
    """v7 §3 手順 0（表の平均ではなく節点の値）"""
    va2 = P2 / G**2
    Ua = (KG / G) / xi
    Fa = 1 - 0.5 * va2 - Ua - va2**2 / 8
    vb2 = ((G - 1) / G)**2 * P2
    Ub = KG * (G - 1) / (G * xi)
    Fb = 1 - 0.5 * vb2 - Ub
    return Fa, Fb


def run(delta, f, nodes, sa, sb, out_csv, tol_turn=1e-3, lock="off"):
    # ---- 初期化（v7 §8：半径 ξ1、近点 π = 0、m = 1 + δ、節点幅 DTH0 を x3a に写す）
    m = 1.0 + delta
    xi = XI1
    pi = 0.0
    P2 = pi**2 + (2 * G * m / xi)**2
    Fa, Fb = clock_factors(xi, P2)
    x3a = Fa * DTH0            # a の位相の一節点の増分 = 節点の時計（以後、v1 では不変）
    K = int(round(x3a / (2 * math.pi)))   # v3: 継承した節点の周数（1616）。lock on では x3a = 2πK + 2θ で閉じる
    lam = math.log(xi)
    dlam = 0.0
    psi = 0.0                  # 関係の位相（配列の回転角の積算）
    x1a = 0.0
    x1b = 0.0
    eps0 = eps(xi, pi, m, sa, sb)
    cum_turns = 0.0
    recur = []                 # ψ が初期に戻る節点（一周の 1/1000 以内）

    fields = ["node", "xi", "pi", "m", "theta", "R_sin2theta", "two_theta_over_2pi", "theta_X",
              "gamma", "F_a", "F_b", "x3b_over_x3a", "eps", "psi", "cum_turns",
              "theta_rad", "sin2_theta_rad", "dE_node", "eta_a", "eta_b",
              "Delta_a", "Delta_b", "Delta_ab", "s_dir", "x3a", "K"]   # v3: 末尾 6 列（先頭 20 列は v1 と同一）
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for k in range(nodes + 1):
            # ---- 手順 0: 読み出し（長さの側 ξ, π；時計の側 Δθ_k, F；反対称の不変量 m）
            xi = math.exp(lam)
            Fa, Fb = clock_factors(xi, pi**2 + (2 * G * m / xi)**2)
            for _ in range(3):                      # Δθ_k と F の自己無撞着（一節点内の反復）
                if lock == "on":
                    # v3: 節点の時計は a の位相が関係の位相と閉じる値。x3a = 2πK + 2θ、θ = (x3a/(4gFa))∂ε/∂m
                    #     ⇒ x3a = 2πK / (1 − ∂ε/∂m /(2 g Fa))。Δa ≡ 0 が恒等的に成り立つ
                    x3a = 2 * math.pi * K / (1.0 - d_m(xi, pi, m, sa, sb) / (2 * G * Fa))
                dth = x3a / Fa
                pi = xi * dlam / dth if k > 0 else pi
                P2 = pi**2 + (2 * G * m / xi)**2
                Fa, Fb = clock_factors(xi, P2)
            if lock == "on":
                x3a = 2 * math.pi * K / (1.0 - d_m(xi, pi, m, sa, sb) / (2 * G * Fa))
            dth = x3a / Fa
            x3b = RATIO * Fb * dth
            e = eps(xi, pi, m, sa, sb)
            # ---- 手順 1: 角（面①②③④）
            dm = d_m(xi, pi, m, sa, sb)
            theta = dth / (4 * G) * dm
            exx = d_xixi(xi, pi, m, sa, sb)
            epp = d_pipi(xi, pi, m, sa, sb)
            theta_X = 0.5 * dth * math.sqrt(max(exx * epp, 0.0))
            eta_a = -math.log(Fa)
            eta_b = -math.log(Fb)
            # ---- 放射の読み（形の群：|D''|²。φ̇ = 2θ/Δθ_k, ξ̇ = ∂ε/∂π, ξ̈ = −∂ε/∂ξ, φ̈ は連鎖律）
            phidot = 2 * theta / dth
            xidot = d_pi(xi, pi, m, sa, sb)
            xiddot = -d_xi(xi, pi, m, sa, sb)
            phiddot = (d_mxi(xi, pi, m, sa, sb) * xidot + d_mpi(xi, pi, m, sa, sb) * xiddot) / (2 * G)
            D2 = (xiddot - xi * phidot**2)**2 + (2 * xidot * phidot + xi * phiddot)**2
            P_D = (2.0 / 3.0) * KC * D2
            # ---- v3: 三つの位相の進み（x3a, x3b, 2θ）の差（STATUS §7.1）。整数関係なら s_dir = 0 ⇒ 放出停止（§7.5）
            Da = wrap_pi(2 * theta - x3a)
            Db = wrap_pi(2 * theta - x3b)
            Dab = wrap_pi(x3b - x3a)
            s_dir = (math.sin(Da / 2)**2 + math.sin(Db / 2)**2 + math.sin(Dab / 2)**2) / 3.0
            s_eff = s_dir if lock == "on" else 1.0
            dE = f * P_D * dth * s_eff
            sin2_rad = dE / abs(e) if e != 0 else 0.0
            theta_rad = math.asin(math.sqrt(min(sin2_rad, 1.0)))
            w.writerow([k, f"{xi:.12e}", f"{pi:.12e}", f"{m:.12e}", f"{theta:.12e}",
                        f"{math.sin(theta)**2:.12e}", f"{theta / math.pi:.12e}", f"{theta_X:.12e}",
                        f"{theta / theta_X if theta_X else float('nan'):.12e}", f"{Fa:.15e}", f"{Fb:.15e}",
                        f"{x3b / x3a:.12e}", f"{e:.12e}", f"{psi:.12e}", f"{cum_turns:.9f}",
                        f"{theta_rad:.6e}", f"{sin2_rad:.6e}", f"{dE:.6e}", f"{eta_a:.6e}", f"{eta_b:.6e}",
                        f"{Da:.9e}", f"{Db:.9e}", f"{Dab:.9e}", f"{s_dir:.9e}", f"{x3a:.9f}", K])
            if k == nodes:
                break
            # ---- 手順 2: 面②（呼吸。E–M の正準更新。シンプレクティック Euler）
            pi = pi - d_xi(xi, pi, m, sa, sb) * dth
            xi_new = xi + d_pi(xi, pi, m, sa, sb) * dth
            dlam = math.log(xi_new) - lam
            lam = math.log(xi_new)
            # ---- 手順 3: 面①（配列の回転：関係の位相は倍角 2θ 進む）
            psi = (psi + 2 * theta) % (2 * math.pi)
            cum_turns += 2 * theta / (2 * math.pi)
            if k > 0 and min(psi, 2 * math.pi - psi) < tol_turn * 2 * math.pi:
                recur.append((k + 1, cum_turns))
            # ---- 時計（面③④の双曲角は F に入っている）
            x1a = (x1a + x3a) % (2 * math.pi)
            x1b = (x1b + x3b) % (2 * math.pi)
            # ---- 手順 4: 放射 = 背景との片側の交換（f）。Δm = −ΔE/(2g φ̇)（∂ε/∂m = 2g φ̇ なので ε は ΔE だけ減る）
            if dE > 0 and phidot > 0:
                m = m - dE / (2 * G * phidot)
    return dict(eps0=eps0, recur=recur[:12], n_recur=len(recur), x3a=x3a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--f", type=float, default=0.0)
    ap.add_argument("--nodes", type=int, default=2000)
    ap.add_argument("--sa", type=int, default=1)
    ap.add_argument("--sb", type=int, default=-1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lock", choices=["off", "on"], default="off",
                    help="off: v1 と同一 / on: 決定済みの引き込み（三位相の整数関係、δτ 動的、源は位相配位の関数）")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    r = run(a.delta, a.f, a.nodes, a.sa, a.sb, a.out, lock=a.lock)
    print(f"delta={a.delta} f={a.f} nodes={a.nodes} lock={a.lock} eps0={r['eps0']:.9e} x3a_final={r['x3a']:.6f} "
          f"recur_count={r['n_recur']} first_recur(node,turns)={r['recur']}")


if __name__ == "__main__":
    main()
