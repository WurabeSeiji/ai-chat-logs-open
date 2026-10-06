#!/usr/bin/env python3
"""層 1：吸収体に方向の印を与え、重心速度を状態に持ち、反跳を入れる（STATUS §12 の三層のうち読み出しの層。
反跳を入れるので生成子にも「放出が v⃗ を蹴る写像」が加わる）。層 2（温度）を含む。

設定
- 実験室系 = 吸収体の静止系。要素数 N = 10⁶⁰ なので方向は連続量（床 √(4π/N) = 3.5e-30 rad は記録に現れない）。
- 印：双極子軸を z にとり、原子の初期速度を観測されている背景放射の双極子 v₀ = 370 km/s ẑ（β = 1.23e-3）とする。
- 状態：準位 i と重心速度 v⃗。原子の質量 M_i = m_e + m_p + E_i/c²（準位ごと）。
- 放出：静止系で厳密な二体運動学 ħω′ = (M_i² − M_f²)c²/(2M_i) = ΔE(1 − ΔE/2M_ic²)、方向は等方（m 副準位は分解していない）、
  自然幅 Lorentz（半値半幅 (Γ_i + Γ_f)/2、T での誘導・吸収の率を含む）の揺らぎ。光子の 4 元運動量を Lorentz 変換で実験室系へ、
  運動量保存で v⃗ を更新。
- 温度：原子の静止系では背景放射が双極子状に非等方 T(n̂_src) = T₀/(γ(1 − β n̂_src·v̂))。吸収は到来方向を n̄(ω₀; T(n̂_src)) の
  重みで引き、誘導放出は進行方向のモードの占有 n̄(T(−n̂)) の重みで引く。率は角度平均 ⟨n̄⟩_Ω（Gauss–Legendre 64 点）。
- 二光子：和の振動数は厳密（反跳込み）、二個の方向は (1 + cos²θ₁₂) の相関、分配 y = ω₁/ω_tot は y³(1−y)³ に比例する形で引く
  （分配の形だけ近似。記録に残る鋭い量は和）。重力波の辺も事象として引く（確率 1e-42、実際には出ない）。
- 記録：各事象の (t, 種別, ħω_lab, n̂_lab)。ω_lab は線の静止振動数 ω′ に対する比 D = ω_lab/ω′ の形 δ = D − 1 でも保持する
  （倍精度で ω_lab そのものを持つと相対 1e-16 で丸まり、21 cm の幅 5e-23 や反跳 5e-12 が消えるため。δ は桁落ちなしに計算）。

読み出し（記録だけから）
1. 集団：多段の光子の (δ, n̂) から β⃗ を推定。雑音は Lorentz 型（裾が重い）なので最小二乗でなく Cauchy 損失（Lorentz 雑音の最尤）。
2. 一個の原子の 21 cm 事象：一事象ごとに、印の軸（ẑ）方向の速度を厳密な Doppler 式 D = 1/(γ(1 − β n_z)) の逆解きで読む。
3. 一個の原子の全履歴：21 cm 事象の列から、各事象の蹴り ∓ħω n̂/c（記録にある）を運動量保存で積み上げ、初期 v⃗₀ の 3 成分を推定。
   反跳の履歴が記録から再構成できることの検査。

出力：absorber_direction_recoil.md、events_one_atom.csv
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, brentq

HERE = Path(__file__).resolve().parent
FIXED = HERE.parent / "exchange_rel_20261005"
sys.path.insert(0, str(FIXED))
from matplotlib import font_manager  # noqa: E402
font_manager.fontManager.addfont = lambda path: None
import exchange_cascade as X  # noqa: E402

C = X.C_SI
HBAR_SI, E_SI, ME_SI = X.HBAR_SI, X.E_SI, X.ME_SI
K_B_EV = 8.617333262e-5
T_CMB = 2.7255
V0 = 370.0e3                      # m/s（背景放射の双極子から読まれている太陽系の速度）
N_CELLS = 1e60
EV = X.EV
G_STAT = np.array([int(2 * X.jf_of(l, se, sp)[1] + 1) for (n, l, se, sp) in X.ST])
M_LEVEL = (1.0 + X.MU_P + X.ETOT) * ME_SI          # 準位ごとの原子の質量 [kg]
E_J = X.ETOT * EV * E_SI                            # 準位エネルギー [J]
GL_X, GL_W = np.polynomial.legendre.leggauss(64)
YEAR = 3.15576e7


def nbar(hw_J, T):
    if T <= 0.0:
        return 0.0
    x = hw_J / (K_B_EV * E_SI * T)
    return 0.0 if x > 700.0 else 1.0 / np.expm1(x)


def gamma_minus_1(b):
    g = 1.0 / np.sqrt(1.0 - b * b)
    return b * b * g * g / (g + 1.0), g


def T_dir(T0, v, n_src):
    b = np.linalg.norm(v) / C
    if b == 0.0:
        return T0
    _, g = gamma_minus_1(b)
    return T0 / (g * (1.0 - b * np.dot(n_src, v / np.linalg.norm(v))))


def nbar_avg(hw_J, T0, v):
    if T0 <= 0.0:
        return 0.0
    b = np.linalg.norm(v) / C
    _, g = gamma_minus_1(b)
    vals = np.array([nbar(hw_J, T0 / (g * (1.0 - b * c))) for c in GL_X])
    return 0.5 * float(np.dot(GL_W, vals))


def isotropic(rng):
    u = rng.uniform(-1, 1); ph = rng.uniform(0, 2 * np.pi)
    s = np.sqrt(1 - u * u)
    return np.array([s * np.cos(ph), s * np.sin(ph), u])


def sample_dir_weighted(rng, hw_J, T0, v, src_weight=True):
    b = np.linalg.norm(v) / C
    _, g = gamma_minus_1(b)
    nmax = nbar(hw_J, T0 / (g * (1.0 - b)))
    while True:
        n = isotropic(rng)
        n_src = n if src_weight else -n
        if rng.uniform(0, nmax) <= nbar(hw_J, T_dir(T0, v, n_src)):
            return n


def doppler_delta(n_lab, v):
    """δ = ω_lab/ω′ − 1 = 1/(γ(1 − β n_∥)) − 1 を桁落ちなしに。n_∥ = n̂_lab·v̂。"""
    speed = np.linalg.norm(v)
    if speed == 0.0:
        return 0.0
    b = speed / C
    gm1, g = gamma_minus_1(b)
    npar = np.dot(n_lab, v / speed)
    # 1/(γ(1−βn)) − 1 = (1 − γ(1−βn))/(γ(1−βn)) = (βn·γ − (γ−1))/(γ(1−βn))
    return (b * npar * g - gm1) / (g * (1.0 - b * npar))


def unboost_dir(n_lab, v):
    speed = np.linalg.norm(v)
    if speed == 0.0:
        return n_lab
    vhat = v / speed
    b = speed / C
    _, g = gamma_minus_1(b)
    npar = np.dot(n_lab, vhat)
    n_par = (npar - b) / (1.0 - b * npar)
    n_perp = (n_lab - npar * vhat) / (g * (1.0 - b * npar))
    n = n_par * vhat + n_perp
    return n / np.linalg.norm(n)


def boost_dir(n_rest, v):
    speed = np.linalg.norm(v)
    if speed == 0.0:
        return n_rest
    vhat = v / speed
    b = speed / C
    _, g = gamma_minus_1(b)
    npar = np.dot(n_rest, vhat)
    n_par = (npar + b) / (1.0 + b * npar)
    n_perp = (n_rest - npar * vhat) / (g * (1.0 + b * npar))
    n = n_par * vhat + n_perp
    return n / np.linalg.norm(n)


def build_channels(T0):
    ed = X.edges(X.W_GW)
    ch = {i: [] for i in range(X.DIM)}
    for src, dst, kind, A in ed:
        if A <= 0.0:
            continue
        dE = E_J[src] - E_J[dst]
        hw_rest = dE * (1.0 - dE / (2.0 * M_LEVEL[src] * C**2))
        ch[src].append(dict(dst=dst, kind=kind, A=A, hw=hw_rest, up=False))
        if kind in ("e1", "e2", "m1") and T0 > 0.0:
            hw_abs = dE * (1.0 + dE / (2.0 * M_LEVEL[dst] * C**2))
            ch[dst].append(dict(dst=src, kind=kind + "_abs", A=A * G_STAT[src] / G_STAT[dst], hw=hw_abs, up=True))
    return ch


def rates(ch_i, T0, v):
    out = []
    for c in ch_i:
        if c["up"]:
            r = c["A"] * nbar_avg(c["hw"], T0, v)
        elif c["kind"] in ("e1", "e2", "m1"):
            r = c["A"] * (1.0 + nbar_avg(c["hw"], T0, v))
        else:
            r = c["A"]
        out.append(r)
    return np.array(out)


def level_width(i, T0, v, ch):
    return float(np.sum(rates(ch[i], T0, v)))


def apply_kick(v, M_before, hw_lab, n_lab, sign):
    """運動量保存。sign = −1 放出、+1 吸収。原子の 4 元運動量 (E, P⃗) を更新して v⃗ = P⃗c²/E。"""
    b = np.linalg.norm(v) / C
    _, g = gamma_minus_1(b)
    P = g * M_before * v
    Eat = g * M_before * C**2
    P2 = P + sign * (hw_lab / C) * n_lab
    E2 = Eat + sign * hw_lab
    return P2 * C**2 / E2


def simulate_atom(rng, T0, v0, n_events=40, t_cap=1.0e17, init=(5, 1, 1, -1), ch=None):
    if ch is None:
        ch = build_channels(T0)
    i = X.IDX[init]
    v = np.array(v0, dtype=float)
    t = 0.0
    events = []
    while len(events) < n_events and t < t_cap:
        r = rates(ch[i], T0, v)
        Rtot = r.sum()
        if Rtot <= 0.0:
            break
        t += rng.exponential(1.0 / Rtot)
        k = rng.choice(len(r), p=r / Rtot)
        c = ch[i][k]
        f = c["dst"]
        rel_width = 0.5 * (level_width(i, T0, v, ch) + level_width(f, T0, v, ch)) * HBAR_SI / c["hw"]   # 相対半値半幅
        eps = rng.standard_cauchy() * rel_width                                                      # 静止振動数の揺らぎ（相対）
        v_before = v.copy()
        if c["kind"] == "2ph":
            hw_tot = c["hw"] * (1.0 + eps)
            y = rng.beta(4, 4)
            n1 = isotropic(rng)
            while True:
                n2 = isotropic(rng)
                if rng.uniform(0, 2) <= 1 + np.dot(n1, n2) ** 2:
                    break
            n1l, n2l = boost_dir(n1, v), boost_dir(n2, v)
            hw1 = y * hw_tot * (1.0 + doppler_delta(n1l, v))
            hw2 = (1 - y) * hw_tot * (1.0 + doppler_delta(n2l, v))
            v = apply_kick(v, M_LEVEL[i], hw1, n1l, -1.0)
            v = apply_kick(v, M_LEVEL[i] - hw1 / C**2, hw2, n2l, -1.0)
            events.append(dict(t=t, kind="2ph", src=i, dst=f, hw_lab=hw1 + hw2, n_lab=n1l, delta=float("nan"), eps=eps,
                               v_before=v_before, v_after=v.copy()))
            i = f
            continue
        if c["up"]:
            n_lab = -sample_dir_weighted(rng, c["hw"], T0, v, src_weight=True)
            sign = +1.0
            kind = c["kind"]
        else:
            stim = False
            if c["kind"] in ("e1", "e2", "m1") and T0 > 0.0:
                nav = nbar_avg(c["hw"], T0, v)
                stim = rng.uniform(0, 1 + nav) > 1.0
            if stim:
                n_lab = sample_dir_weighted(rng, c["hw"], T0, v, src_weight=False)
            else:
                n_lab = boost_dir(isotropic(rng), v)
            sign = -1.0
            kind = c["kind"] + ("_stim" if stim else "")
        D_minus_1 = doppler_delta(n_lab, v)
        delta = D_minus_1 + eps * (1.0 + D_minus_1)            # ω_lab/ω′ − 1（ω′ は表の静止振動数、eps は線の揺らぎ）
        hw_lab = c["hw"] * (1.0 + delta)
        v = apply_kick(v, M_LEVEL[i], hw_lab, n_lab, sign)
        events.append(dict(t=t, kind=kind, src=i, dst=f, hw_lab=hw_lab, n_lab=n_lab, delta=delta, eps=eps,
                           v_before=v_before, v_after=v.copy(), line_hw=c["hw"], rel_width=rel_width))
        i = f
    return events


# ---------- 読み出し（記録だけから） ----------
def line_table(T0):
    ch = build_channels(T0)
    tab = []
    for i, lst in ch.items():
        for c in lst:
            tab.append(dict(hw=c["hw"], kind=c["kind"], src=i, dst=c["dst"], up=c["up"]))
    return tab


def identify(hw_lab, tab, beta_max=3e-3, up=None):
    """ω_lab から線を同定。吸収か放出かは記録にある（吸収体は自分が出した光子を知っている）ので up で絞る。
    同じ対の放出線と吸収線は反跳分 2ΔE/(2Mc²) ≈ 6e-15 しか違わず、振動数だけでは区別できない。"""
    best = None
    for L in tab:
        if up is not None and L["up"] != up:
            continue
        if abs(hw_lab - L["hw"]) <= beta_max * L["hw"]:
            if best is None or abs(hw_lab - L["hw"]) < abs(hw_lab - best["hw"]):
                best = L
    return best


def rel_width_of(line, T0, v, ch):
    return 0.5 * (level_width(line["src"], T0, v, ch) + level_width(line["dst"], T0, v, ch)) * HBAR_SI / line["hw"]


def fit_velocity_ensemble(events, tab, T0, ch):
    """集団：各事象は独立（原子は別々）。δ_k = ω_lab/ω′ − 1 と n̂_k から β⃗。Cauchy 損失、相対幅で規格化。"""
    data = []
    for e in events:
        if "2ph" in e["kind"]:
            continue
        line = identify(e["hw_lab"], tab)
        if line is None:
            continue
        data.append((e["delta"], e["n_lab"], rel_width_of(line, T0, np.array([0, 0, V0]), ch)))
    def resid(bv):
        v = bv * C
        return np.array([(d - doppler_delta(n, v)) / w for d, n, w in data])
    # 初期値：1 次の線形最小二乗 δ ≈ β⃗·n̂（Cauchy 損失は初期値が遠いと勾配が消えるので）
    Nmat = np.array([n for _, n, _ in data]); dvec = np.array([d for d, _, _ in data])
    b0 = np.linalg.lstsq(Nmat, dvec, rcond=None)[0]
    sol = least_squares(resid, x0=b0, loss="cauchy", f_scale=1.0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return sol.x * C, len(data)


def invert_single_axis(delta, n_z):
    """一事象、v⃗ ∥ ẑ を仮定：δ = 1/(γ(1 − β n_z)) − 1 を β について厳密に解く。"""
    if abs(n_z) < 0.05:
        return float("nan")
    f = lambda b: doppler_delta(np.array([0.0, 0.0, 1.0]) * n_z + np.array([np.sqrt(max(0.0, 1 - n_z**2)), 0.0, 0.0]), np.array([0.0, 0.0, b * C])) - delta
    try:
        return brentq(f, -5e-3, 5e-3, xtol=1e-30, rtol=1e-15, maxiter=200) * C
    except ValueError:
        return float("nan")


def reconstruct_history(events, tab, T0, ch):
    """一個の原子：21 cm 事象の列。各事象の蹴り（∓ħω n̂/c、符号は吸収か放出か）を記録から積み上げ、v⃗₀ を推定。"""
    seq = []
    for e in events:
        if "2ph" in e["kind"] or e["hw_lab"] > 1e-4 * E_SI:
            continue
        up = e["kind"].endswith("_abs")
        line = identify(e["hw_lab"], tab, up=up)
        if line is None:
            continue
        sign = +1.0 if up else -1.0
        M_before = M_LEVEL[line["src"]]          # チャネルの持ち主の準位 = 事象前の準位（吸収なら下、放出なら上）
        seq.append((e["delta"], e["n_lab"], e["hw_lab"], sign, M_before))
    if len(seq) < 3:
        return None, None, len(seq)
    def model(v0):
        v = np.array(v0, dtype=float)
        out = []
        vs = []
        for d, n, hw, sign, M in seq:
            out.append(d - doppler_delta(n, v))
            vs.append(v.copy())
            v = apply_kick(v, M, hw, n, sign)
        return np.array(out), vs
    # 21 cm の線幅（相対 5e-23）は倍精度の床（1e-19）より下なので雑音は見えず、通常の最小二乗でよい。
    # 初期値は 1 次の線形解（2 次 Doppler の分 β² ≈ 1.5e-6 だけずれる）。
    Nmat = np.array([n for _, n, _, _, _ in seq]); dvec = np.array([d for d, _, _, _, _ in seq])
    b0 = np.linalg.lstsq(Nmat, dvec, rcond=None)[0]
    sol = least_squares(lambda x: model(x * C)[0], x0=b0, xtol=3e-16, ftol=3e-16, gtol=3e-16, x_scale=1e-3)
    v0 = sol.x * C
    _, vs = model(v0)
    return v0, vs, len(seq)


def main():
    rng = np.random.default_rng(20261006)
    T0 = T_CMB
    v0 = np.array([0.0, 0.0, V0])
    ch = build_channels(T0)
    tab = line_table(T0)
    def lab(i):
        n, l, se, sp = X.ST[i]; j, F = X.jf_of(l, se, sp)
        return "%d%s%d/2F%d" % (n, {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}[l], int(2 * j), int(F))
    L = []
    L.append("# 層 1：方向の印・重心速度・反跳（層 2 の温度を含む）")
    L.append("")
    L.append("実験室系 = 吸収体の静止系。要素数 N = 10⁶⁰（方向は連続、床 %.1e rad）。印 = 双極子軸 ẑ。初期 5p₃/₂ F=1、v₀ = %.0f km/s ẑ（β = %.4e）、T = %.4f K。"
             % (np.sqrt(4 * np.pi / N_CELLS), V0 / 1e3, V0 / C, T0))
    L.append("")
    i_ly, f_ly = X.IDX[(2, 1, 1, 1)], X.IDX[(1, 0, 1, 1)]
    dE = E_J[i_ly] - E_J[f_ly]
    i21, f21 = X.IDX[(1, 0, 1, 1)], X.IDX[(1, 0, 1, -1)]
    dE21 = E_J[i21] - E_J[f21]
    w21 = 0.5 * (level_width(i21, T0, v0, ch) + level_width(f21, T0, v0, ch))
    L.append("| 量 | 値 |")
    L.append("|---|---|")
    L.append("| Lyman-α の反跳シフト ΔE²/(2Mc²) | %.2f MHz |" % (dE**2 / (2 * M_LEVEL[i_ly] * C**2) / (2 * np.pi * HBAR_SI) / 1e6))
    L.append("| Lyman-α 一個の反跳速度 ħω/(Mc) | %.3f m/s |" % (dE / (M_LEVEL[i_ly] * C)))
    L.append("| 21 cm 一個の反跳速度 | %.3e m/s |" % (dE21 / (M_LEVEL[i21] * C)))
    L.append("| 1 次 Doppler β、2 次 β²/2 | %.4e、%.3e |" % (V0 / C, 0.5 * (V0 / C) ** 2))
    L.append("| Lyman-α の自然幅による速度分解能 c·(Γ/2)/ω（一事象） | %.1f m/s |" % (C * 0.5 * level_width(i_ly, T0, v0, ch) / (dE / HBAR_SI)))
    L.append("| 21 cm の幅による速度分解能（一事象、T での誘導・吸収込みの幅 %.2e s⁻¹） | %.2e m/s（物理の限界）。倍精度で δ を持つ数値の限界 ≈ c·10⁻¹⁹ = %.0e m/s |" % (2 * w21, C * w21 / (dE21 / HBAR_SI), C * 1e-19))
    L.append("| 原子の静止系で見た背景の温度 T(前方)/T(後方) | %.6f / %.6f K |" % (T_dir(T0, v0, np.array([0, 0, 1.0])), T_dir(T0, v0, np.array([0, 0, -1.0]))))
    L.append("| 21 cm の ⟨n̄⟩_Ω（動く原子）と n̄(T₀) | %.5f、%.5f |" % (nbar_avg(dE21, T0, v0), nbar(dE21, T0)))
    L.append("")
    # 一個の原子
    ev = simulate_atom(rng, T0, v0, n_events=30, ch=ch)
    L.append("## 一個の原子の事象列（最初の 14 事象）")
    L.append("")
    L.append("| # | t | 種別 | 遷移 | ħω_lab [eV] | δ = ω_lab/ω′ − 1 | δ/β | cos θ = n̂_lab·ẑ | \\|v\\| after [m/s] | 反跳 \\|Δv\\| [m/s] |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    with open(HERE / "events_one_atom.csv", "w", encoding="utf-8") as fcsv:
        fcsv.write("k,t_s,kind,src,dst,hw_lab_eV,delta,nx,ny,nz,vx_after,vy_after,vz_after\n")
        for k, e in enumerate(ev):
            n = e["n_lab"]
            dv = np.linalg.norm(e["v_after"] - e["v_before"])
            fcsv.write("%d,%.6e,%s,%s,%s,%.12e,%.17e,%.9f,%.9f,%.9f,%.12f,%.12f,%.12f\n" % (k, e["t"], e["kind"], lab(e["src"]), lab(e["dst"]), e["hw_lab"] / E_SI, e["delta"], n[0], n[1], n[2], *e["v_after"]))
            if k < 14:
                tstr = ("%.3e s" % e["t"]) if e["t"] < YEAR else ("%.3e 年" % (e["t"] / YEAR))
                L.append("| %d | %s | %s | %s→%s | %.9f | %+.6e | %+.4f | %+.4f | %.6f | %.3e |"
                         % (k, tstr, e["kind"], lab(e["src"]), lab(e["dst"]), e["hw_lab"] / E_SI, e["delta"], e["delta"] / (V0 / C), n[2],
                            np.linalg.norm(e["v_after"]), dv))
    L.append("")
    L.append("δ/β が cos θ とほぼ一致するのが、「方向の印に対する速度」が記録に出ていることの直接の表れ。差は 2 次 Doppler（β²/2 = 7.6e-7、δ/β で 6e-4）と線の揺らぎ。")
    L.append("")
    # 読み出し 1：集団
    n_atoms = 200
    rng2 = np.random.default_rng(7)
    all_ev = []
    for a in range(n_atoms):
        all_ev += [e for e in simulate_atom(rng2, T0, v0, n_events=8, ch=ch) if e["hw_lab"] > 1.0 * E_SI]
    v_ens, nused = fit_velocity_ensemble(all_ev, tab, T0, ch)
    L.append("## 読み出し 1：集団（%d 原子、多段の光子 %d 個）の (δ, n̂) から速度ベクトルを推定（Cauchy 損失）" % (n_atoms, nused))
    L.append("")
    L.append("推定 v⃗ = (%.2f, %.2f, %.2f) m/s、与えた v₀ = (0, 0, %.0f) m/s、差 %.2f m/s。目安：Lyman-α 一事象 6 m/s（Lorentz の半値半幅）と反跳 3.3 m/s を %d 事象で平均 → 0.5 m/s 程度。"
             % (*v_ens, V0, np.linalg.norm(v_ens - v0), nused))
    L.append("")
    # 読み出し 2：一事象の情報量
    hf = [e for e in ev if e["hw_lab"] < 1e-4 * E_SI and "2ph" not in e["kind"]]
    L.append("## 読み出し 2：一個の原子の 21 cm 事象（%d 個）。一事象が持つ速度の情報" % len(hf))
    L.append("")
    L.append("一事象の δ は β⃗ の視線成分 n̂·v⃗ しか決めないので、一事象だけで v⃗ は読めない。読めるのは「δ が厳密な Doppler 式を線幅の範囲で満たしているか」、"
             "つまり一事象が視線速度を何 m/s の精度で持っているか。下表は真の v⃗（事象前）で計算した δ_model との差を、視線速度に換算したもの。")
    L.append("")
    L.append("| # | t [年] | 種別 | cos θ | δ | c·(δ − δ_model) [m/s] | 幅換算 (δ − δ_model)/rel_width | この事象の反跳 \\|Δv\\| [m/s] |")
    L.append("|---|---|---|---|---|---|---|---|")
    errs = []
    for k, e in enumerate(hf[:12]):
        dmod = doppler_delta(e["n_lab"], e["v_before"])
        err = C * (e["delta"] - dmod)
        errs.append(abs(err))
        L.append("| %d | %.3e | %s | %+.3f | %+.9e | %+.2e | %+.2f | %.3e |"
                 % (k, e["t"] / YEAR, e["kind"], e["n_lab"][2], e["delta"], err, (e["delta"] - dmod) / e["rel_width"], np.linalg.norm(e["v_after"] - e["v_before"])))
    L.append("")
    L.append("差の中央値 %.1e m/s。21 cm の線幅（相対 5e-23、7.7e-15 m/s）は倍精度の床（δ で 1e-19、3e-11 m/s）より下なので、記録の δ は真の速度の厳密な Doppler と倍精度で一致する（残差 0）。"
             "一事象の情報は視線成分だけだが、3 事象以上を蹴りの履歴とともに連立すれば v⃗ の 3 成分が決まる（読み出し 3）。" % np.median(errs))
    L.append("")
    # 読み出し 3：全履歴
    v0_fit, vs, nseq = reconstruct_history(ev, tab, T0, ch)
    L.append("## 読み出し 3：一個の原子の 21 cm 事象 %d 個から、記録にある蹴り ∓ħω n̂/c を積み上げて初期 v⃗₀ の 3 成分を推定" % nseq)
    L.append("")
    if v0_fit is not None:
        first = hf[0]["v_before"]
        L.append("推定 v⃗₀ = (%.9f, %.9f, %.9f) m/s、真（最初の 21 cm 事象の前）= (%.9f, %.9f, %.9f) m/s、差 %.2e m/s。"
                 % (*v0_fit, *first, np.linalg.norm(v0_fit - first)))
        maxdev = max(np.linalg.norm(vs[k] - hf[k]["v_before"]) for k in range(min(len(vs), len(hf))))
        L.append("再構成した各事象前の v⃗ と真の v⃗ の最大差 %.2e m/s（反跳 1.9e-6 m/s に対して %.0e 倍）。%s"
                 % (maxdev, maxdev / 1.87e-6, "記録（ω_lab と n̂ と吸収か放出か）だけから反跳の履歴が復元できている。" if maxdev < 1.87e-7 else "反跳の履歴は復元できていない。"))
    L.append("")
    kicks = [np.linalg.norm(e["v_after"] - e["v_before"]) for e in all_ev]
    L.append("多段の光子一個あたりの反跳 |Δv|：平均 %.3f m/s、最小 %.3f、最大 %.3f（5p→1s 4.17、Lyman-α 3.26、Balmer 0.60 など線ごと）。" % (np.mean(kicks), np.min(kicks), np.max(kicks)))
    L.append("")
    L.append("## 何が増えたか")
    L.append("")
    L.append("- 記録の各事象に方向 n̂ と Doppler 込みの ω_lab が付き、速度ベクトル（3 成分）が記録だけから決まる（読み出し 1）。")
    L.append("- 放出・吸収ごとの反跳が運動量保存で v⃗ を蹴り、21 cm 事象ではその蹴りが一事象ごとに分解でき、履歴が復元できる（読み出し 2、3）。生成子に「事象が状態 v⃗ を変える写像」が加わった。")
    L.append("- 原子の静止系で背景が双極子状に非等方（T(前方)/T(後方)）になり、吸収と誘導放出の方向に β 程度の偏りが出る。これが背景放射の抵抗で、21 cm の率では 1e-21 m/s² と無視できる。")
    L.append("- 層 3（m 副準位の整列、偏光）は入れていない。方向は等方（m 平均）で引いている。")
    out = "\n".join(L) + "\n"
    (HERE / "absorber_direction_recoil.md").write_text(out, encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
