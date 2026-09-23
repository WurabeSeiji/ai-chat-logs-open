#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数値実験02: 強い場の2体の観測データ生成（観測1 軌道観測・観測2 重力波観測）
（匿名頂点状態生成幾何 ― 匿名内部観測者から不変な関係量の体系）

設定
    単位は G = c = M = 1（M は全質量 m_A + m_B）。長さと時間は M を単位とする。
    天体 A, B は点とし、自転・電荷は無視する。A を重いほう（m_A >= m_B）とする。
    軌道の形のパラメータ a, b は、初期時刻の軌道の長半径・短半径（M 単位、調和座標）。
      近点距離 a(1-e) と遠点距離 a(1+e)（e = sqrt(1 - b^2/a^2)）が、放射反作用を除いた
      2PN の運動の実際の折り返し点になるよう、近点での接線速度を数値的に決める。
    観測者 c は十分に軽く十分に遠く、軌道面の法線方向にいる。ドップラー効果は無視する。

軌道の計算
    相対運動（x = x_A - x_B）を調和座標のポストニュートン近似で積分する。
      1PN・2PN の保存項と 2.5PN の放射反作用（Kidder, Phys. Rev. D 52, 821 (1995) 式 (2.2)）
    個々の天体の位置は 1PN の重心の関係式で求める。
    A, B の固有時間は、相手の場だけを考えた 1PN の計量で積分する
      dτ_A/dt = 1 - m_B/r - v_A^2/2,  dτ_B/dt = 1 - m_A/r - v_B^2/2

観測データ（input/。解析側に渡すのはこれだけ）
    観測1 obs1_Dk.csv : t_s, body1_x_nrad, body1_y_nrad, body2_x_nrad, body2_y_nrad
        時刻は観測者 c の時計 [s]、位置は天球上の角度 [nrad]。
        分解能は、軌道の大きさ a の 1/1000 に相当する正規分布の誤差。
    観測2 obs2_Dk.csv : t_s, h_plus, h_cross
        重力波の歪み（四重極公式、観測者は法線方向）。分解能は振幅の最大値の 1/1000。
    伏せているもの: 距離 D、M の秒換算 T_M、時刻の原点、天球上の回転・原点・鏡映、
                    どちらの列が重い天体か。

正解（truth/）
    Dk_meta.json       生成パラメータ・伏せたパラメータ・自己検証の結果
    Dk_timeseries.csv  正規化した単位（M 単位）での誤差なしの時系列
    manifest.json      データセットID と条件の対応表、全体設定

使い方
    python generate_binary_pn_observations.py
    python generate_binary_pn_observations.py --out 出力先フォルダ
"""

import argparse
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

GENERATOR_VERSION = "1.0.0"
MASTER_SEED = 2026092202

# (条件ID, ラベル, a, b, 質量比 q = m_A/m_B, 計算する周回数)
CASES = [
    ("C1", "strong",        30.0,  27.0, 1.0, 20),
    ("C2", "intermediate",  60.0,  50.0, 2.0, 30),
    ("C3", "moderate",     120.0, 112.0, 4.0, 30),
]

STEPS_PER_ORBIT = 4000          # 積分の刻み（局所的な力学時間 r^1.5 の 2π/4000）
R_STOP = 8.0                    # これより近づいたら打ち切る（ポストニュートン近似の限界の手前）
OBS1_SAMPLES_PER_ORBIT = 48     # 観測1の時間間隔 = 初期周期 / 48
OBS2_OVERSAMPLING = 4           # 観測2は観測1の4倍の細かさ（初期周期 / 192）
OBS1_RESOLUTION = 1.0e-3        # 観測1の分解能（a に対する比）
OBS2_RESOLUTION = 1.0e-3        # 観測2の分解能（振幅の最大値に対する比）
T_SUN_S = 4.925490947e-6        # GM_sun / c^3 [s]
M_TOTAL_SUN_RANGE = (20.0, 80.0)    # 全質量の秒換算に使う太陽質量の範囲（伏せる）
DISTANCE_RANGE = (3.0e5, 3.0e6)     # 観測者までの距離 D [M]（対数一様、伏せる）
SKY_OFFSET_RANGE = (-2.0e5, 2.0e5)  # 天球上の原点のずれ [nrad]（伏せる）
TIME_OFFSET_RANGE = (0.0, 1000.0)   # 時刻の原点のずれ [s]（伏せる）
VALIDATION_ORBITS = 3               # 自己検証で追加計算する周回数


# ---------------------------------------------------------------- 運動方程式

def pn_acceleration(x, y, vx, vy, nu, rr):
    """相対加速度（調和座標、G = c = M = 1）。Kidder (1995) 式 (2.2)。"""
    r = math.hypot(x, y)
    nx, ny = x / r, y / r
    v2 = vx * vx + vy * vy
    rd = nx * vx + ny * vy
    mr = 1.0 / r
    a1 = (1.0 + 3.0 * nu) * v2 - 2.0 * (2.0 + nu) * mr - 1.5 * nu * rd * rd
    b1 = -2.0 * (2.0 - nu) * rd
    a2 = (0.75 * (12.0 + 29.0 * nu) * mr * mr
          + nu * (3.0 - 4.0 * nu) * v2 * v2
          + 1.875 * nu * (1.0 - 3.0 * nu) * rd ** 4
          - 1.5 * nu * (3.0 - 4.0 * nu) * v2 * rd * rd
          - 0.5 * nu * (13.0 - 4.0 * nu) * mr * v2
          - (2.0 + 25.0 * nu + 2.0 * nu * nu) * mr * rd * rd)
    b2 = -0.5 * rd * (nu * (15.0 + 4.0 * nu) * v2
                      - (4.0 + 41.0 * nu + 8.0 * nu * nu) * mr
                      - 3.0 * nu * (3.0 + 2.0 * nu) * rd * rd)
    f = -mr * mr
    ax = f * ((1.0 + a1 + a2) * nx + (b1 + b2) * vx)
    ay = f * ((1.0 + a1 + a2) * ny + (b1 + b2) * vy)
    if rr:
        c = 1.6 * nu * mr ** 3
        an = rd * (18.0 * v2 + (2.0 / 3.0) * mr - 25.0 * rd * rd)
        bv = -(6.0 * v2 - 2.0 * mr - 15.0 * rd * rd)
        ax += c * (an * nx + bv * vx)
        ay += c * (an * ny + bv * vy)
    return ax, ay


def pn_energy(x, y, vx, vy, nu, order=2):
    """換算質量あたりの保存エネルギー（調和座標、2PN まで）。Kidder (1995) 式 (2.7)。"""
    r = np.hypot(x, y)
    mr = 1.0 / r
    v2 = vx * vx + vy * vy
    rd = (x * vx + y * vy) / r
    e = 0.5 * v2 - mr
    if order >= 1:
        e = e + (0.375 * (1.0 - 3.0 * nu) * v2 * v2 + 0.5 * (3.0 + nu) * v2 * mr
                 + 0.5 * nu * mr * rd * rd + 0.5 * mr * mr)
    if order >= 2:
        e = e + (0.3125 * (1.0 - 7.0 * nu + 13.0 * nu * nu) * v2 ** 3
                 + 0.125 * (21.0 - 23.0 * nu - 27.0 * nu * nu) * mr * v2 * v2
                 + 0.25 * nu * (1.0 - 15.0 * nu) * mr * v2 * rd * rd
                 - 0.375 * nu * (1.0 - 3.0 * nu) * mr * rd ** 4
                 - 0.25 * (2.0 + 15.0 * nu) * mr ** 3
                 + 0.125 * (14.0 - 55.0 * nu + 4.0 * nu * nu) * mr * mr * v2
                 + 0.125 * (4.0 + 69.0 * nu + 12.0 * nu * nu) * mr * mr * rd * rd)
    return e


def rhs(state, nu, m_a, m_b, rr):
    """状態 (x, y, vx, vy, τ_A, τ_B) の時間微分。"""
    x, y, vx, vy = state[0], state[1], state[2], state[3]
    ax, ay = pn_acceleration(x, y, vx, vy, nu, rr)
    r = math.hypot(x, y)
    v2 = vx * vx + vy * vy
    return (vx, vy, ax, ay,
            1.0 - m_b / r - 0.5 * m_b * m_b * v2,
            1.0 - m_a / r - 0.5 * m_a * m_a * v2)


def rk4_step(state, k1, h, nu, m_a, m_b, rr):
    k2 = rhs(tuple(s + 0.5 * h * k for s, k in zip(state, k1)), nu, m_a, m_b, rr)
    k3 = rhs(tuple(s + 0.5 * h * k for s, k in zip(state, k2)), nu, m_a, m_b, rr)
    k4 = rhs(tuple(s + h * k for s, k in zip(state, k3)), nu, m_a, m_b, rr)
    return tuple(s + h / 6.0 * (p + 2.0 * q + 2.0 * u + w)
                 for s, p, q, u, w in zip(state, k1, k2, k3, k4))


def next_apoapsis(rp, vp, nu, m_a, m_b, steps_per_orbit=STEPS_PER_ORBIT):
    """放射反作用を除いた運動で、近点 (rp, 0) から接線速度 vp で出発したときの次の遠点距離。"""
    state = (rp, 0.0, 0.0, vp, 0.0, 0.0)
    eta = 2.0 * math.pi / steps_per_orbit
    t, hist = 0.0, [(0.0, rp)]
    k1 = rhs(state, nu, m_a, m_b, False)
    for _ in range(50 * steps_per_orbit):
        r = math.hypot(state[0], state[1])
        state = rk4_step(state, k1, eta * r ** 1.5, nu, m_a, m_b, False)
        t += eta * r ** 1.5
        k1 = rhs(state, nu, m_a, m_b, False)
        r_new = math.hypot(state[0], state[1])
        hist.append((t, r_new))
        rd = (state[0] * state[2] + state[1] * state[3]) / r_new
        if len(hist) >= 3 and rd <= 0.0:
            if len(hist) == 3 and hist[1][1] <= hist[0][1]:
                return rp          # 出発点がすでに遠点（速度が足りない）
            (t0, r0), (t1, r1), (t2, r2) = hist[-3:]
            c2, c1, c0 = np.polyfit([t0 - t1, 0.0, t2 - t1], [r0, r1, r2], 2)
            return float(c0 - c1 * c1 / (4.0 * c2)) if c2 < 0 else float(max(r0, r1, r2))
    return float("inf")            # 遠点に達しない（束縛されていない）


def initial_state(a, b, nu, m_a, m_b):
    """近点から出発する初期状態。遠点距離が a(1+e) になる接線速度を二分法で決める。"""
    e = math.sqrt(1.0 - (b / a) ** 2)
    rp, ra = a * (1.0 - e), a * (1.0 + e)
    v_newton = math.sqrt((1.0 + e) / (a * (1.0 - e)))
    lo, hi = 0.9 * v_newton, v_newton
    while next_apoapsis(rp, lo, nu, m_a, m_b) > ra:
        lo *= 0.95
    while next_apoapsis(rp, hi, nu, m_a, m_b) < ra:
        hi *= 1.02
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if next_apoapsis(rp, mid, nu, m_a, m_b) < ra:
            lo = mid
        else:
            hi = mid
        if (hi - lo) < 1.0e-15 * hi:
            break
    vp = 0.5 * (lo + hi)
    info = {"r_periapsis": rp, "r_apoapsis_target": ra,
            "r_apoapsis_achieved": next_apoapsis(rp, vp, nu, m_a, m_b),
            "v_periapsis": vp, "v_periapsis_newtonian": v_newton}
    return (rp, 0.0, 0.0, vp, 0.0, 0.0), e, info


def integrate(a, b, nu, m_a, m_b, t_end, rr=True, steps_per_orbit=STEPS_PER_ORBIT, init=None):
    """4次のルンゲ・クッタ法（刻みは局所的な力学時間に比例）。全ステップを返す。"""
    state = init if init is not None else initial_state(a, b, nu, m_a, m_b)[0]
    eta = 2.0 * math.pi / steps_per_orbit
    t = 0.0
    ts, states, derivs = [t], [state], [rhs(state, nu, m_a, m_b, rr)]
    stop = "t_end"
    while t < t_end:
        r = math.hypot(state[0], state[1])
        if r < R_STOP:
            stop = "r_stop"
            break
        h = eta * r ** 1.5
        state = rk4_step(state, derivs[-1], h, nu, m_a, m_b, rr)
        t += h
        ts.append(t)
        states.append(state)
        derivs.append(rhs(state, nu, m_a, m_b, rr))
    return np.array(ts), np.array(states), np.array(derivs), stop


def interpolate(ts, states, derivs, tq):
    """各成分を、その時間微分を使った3次エルミート補間で求める。"""
    idx = np.clip(np.searchsorted(ts, tq, side="right") - 1, 0, len(ts) - 2)
    t0 = ts[idx]
    h = ts[idx + 1] - t0
    s = (tq - t0) / h
    h00 = 2 * s ** 3 - 3 * s ** 2 + 1
    h10 = s ** 3 - 2 * s ** 2 + s
    h01 = -2 * s ** 3 + 3 * s ** 2
    h11 = s ** 3 - s ** 2
    return (h00[:, None] * states[idx] + (h10 * h)[:, None] * derivs[idx]
            + h01[:, None] * states[idx + 1] + (h11 * h)[:, None] * derivs[idx + 1])


# ---------------------------------------------------------------- 自己検証

def periapsis_passages(t, x, y):
    """相対距離の極小（近点通過）の時刻と、そのときの位置角（連続化）を返す。"""
    r = np.hypot(x, y)
    phi = np.unwrap(np.arctan2(y, x))
    i = np.where((r[1:-1] < r[:-2]) & (r[1:-1] <= r[2:]))[0] + 1
    den = r[i - 1] - 2.0 * r[i] + r[i + 1]
    d = np.where(den != 0.0, 0.5 * (r[i - 1] - r[i + 1]) / den, 0.0)
    tp = t[i] + d * (t[i + 1] - t[i - 1]) * 0.5
    pp = phi[i] + d * (phi[i + 1] - phi[i - 1]) * 0.5
    ok = i > 1   # 出発点（近点）そのものは除く
    return np.concatenate(([t[0]], tp[ok])), np.concatenate(([phi[0]], pp[ok]))


def validate(a, b, nu, m_a, m_b, t0_period, ts, st, dv, samples, init):
    res = {}
    e0 = math.sqrt(1.0 - (b / a) ** 2)
    p0 = a * (1.0 - e0 * e0)
    t_val = VALIDATION_ORBITS * t0_period

    # (1) 刻みを半分にした計算との一致（収束）
    ts2, st2, dv2, _ = integrate(a, b, nu, m_a, m_b, t_val * 1.001, rr=True,
                                 steps_per_orbit=2 * STEPS_PER_ORBIT, init=init)
    q = np.array([t_val])
    y1 = interpolate(ts, st, dv, q)[0]
    y2 = interpolate(ts2, st2, dv2, q)[0]
    res["convergence_position_diff_over_a"] = float(math.hypot(y1[0] - y2[0], y1[1] - y2[1]) / a)

    # (2) 放射反作用なしでの保存エネルギーの変動（1PN までと 2PN までを比較）
    tc, sc, _, _ = integrate(a, b, nu, m_a, m_b, t_val, rr=False, init=init)
    e1 = pn_energy(sc[:, 0], sc[:, 1], sc[:, 2], sc[:, 3], nu, order=1)
    e2 = pn_energy(sc[:, 0], sc[:, 1], sc[:, 2], sc[:, 3], nu, order=2)
    res["conservative_energy_variation_1pn"] = float(np.max(np.abs(e1 - e1[0])) / abs(e1[0]))
    res["conservative_energy_variation_2pn"] = float(np.max(np.abs(e2 - e2[0])) / abs(e2[0]))

    # (3) 近点移動（1周あたり）と最低次の式 6π/p の比
    tp, pp = periapsis_passages(samples["t"], samples["x"], samples["y"])
    k = min(5, len(tp) - 1)
    dphi = (pp[k] - pp[0]) / k - 2.0 * math.pi if k >= 1 else float("nan")
    res["periapsis_passages_found"] = int(len(tp))
    res["precession_per_orbit_rad"] = float(dphi)
    res["precession_leading_order_rad"] = float(6.0 * math.pi / p0)
    res["precession_ratio_to_leading_order"] = float(dphi / (6.0 * math.pi / p0))

    # (4) エネルギーの減り方と Peters の式（最低次、周回平均）の比
    ee = pn_energy(samples["x"], samples["y"], samples["vx"], samples["vy"], nu, order=2)
    kk = min(3, len(tp) - 1)
    if kk >= 1:
        e_at = np.interp(tp[[0, kk]], samples["t"], ee)
        slope = (e_at[1] - e_at[0]) / (tp[kk] - tp[0])
        f_e = (1.0 + 73.0 / 24.0 * e0 ** 2 + 37.0 / 96.0 * e0 ** 4) / (1.0 - e0 ** 2) ** 3.5
        peters = -6.4 * nu * a ** -5 * f_e
        res["energy_loss_rate_per_mu"] = float(slope)
        res["energy_loss_rate_peters_per_mu"] = float(peters)
        res["energy_loss_ratio_to_peters"] = float(slope / peters)
    else:
        res["energy_loss_ratio_to_peters"] = float("nan")

    checks = {
        "convergence": res["convergence_position_diff_over_a"] < 1.0e-7,
        "conservative_energy_2pn_improves_on_1pn":
            res["conservative_energy_variation_2pn"] < 0.5 * res["conservative_energy_variation_1pn"],
        "precession_sane": 0.5 < res["precession_ratio_to_leading_order"] < 2.0,
        "energy_loss_sane": 0.5 < res["energy_loss_ratio_to_peters"] < 1.5,
    }
    res["checks"] = checks
    res["all_passed"] = bool(all(checks.values()))
    return res


def pn_scaling_test(nu=0.25, e=0.436, a_values=(60.0, 120.0, 240.0, 480.0)):
    """実装全体の検査: 放射反作用なしで1周積分し、保存エネルギーの変動が
    1PN まで打ち切ると x^2、2PN まで打ち切ると x^3 以上で小さくなることを確かめる。"""
    m_a = m_b = 0.5 if nu == 0.25 else None
    if m_a is None:
        disc = math.sqrt(1.0 - 4.0 * nu)
        m_a, m_b = 0.5 * (1.0 + disc), 0.5 * (1.0 - disc)
    rows = []
    for a in a_values:
        b = a * math.sqrt(1.0 - e * e)
        init, _, _ = initial_state(a, b, nu, m_a, m_b)
        _, sc, _, _ = integrate(a, b, nu, m_a, m_b, 2.0 * math.pi * a ** 1.5, rr=False, init=init)
        var = []
        for order in (1, 2):
            en = pn_energy(sc[:, 0], sc[:, 1], sc[:, 2], sc[:, 3], nu, order=order)
            var.append(float(np.max(np.abs(en - en[0])) / abs(en[0])))
        rows.append({"a": a, "x_periapsis": 1.0 / (a * (1.0 - e)),
                     "energy_variation_1pn": var[0], "energy_variation_2pn": var[1]})
    lx = np.log([r_["x_periapsis"] for r_ in rows])
    s1 = float(np.polyfit(lx, np.log([r_["energy_variation_1pn"] for r_ in rows]), 1)[0])
    s2 = float(np.polyfit(lx, np.log([r_["energy_variation_2pn"] for r_ in rows]), 1)[0])
    return {"nu": nu, "e": e, "rows": rows, "slope_1pn": s1, "slope_2pn": s2,
            "passed": bool(s1 > 1.5 and s2 > 2.7)}


# ---------------------------------------------------------------- 書き出し

def to_jsonable(obj):
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [to_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


def write_json(path, obj):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(to_jsonable(obj), f, ensure_ascii=False, indent=2)
        f.write("\n")


def write_csv(path, columns, fmt):
    names = list(columns.keys())
    fmts = ["{:d}" if np.issubdtype(np.asarray(columns[k]).dtype, np.integer) else fmt for k in names]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(names) + "\n")
        for row in zip(*(columns[k] for k in names)):
            f.write(",".join(g.format(v) for g, v in zip(fmts, row)) + "\n")


# ---------------------------------------------------------------- データセット1つ

def make_dataset(index, case, seed):
    case_id, label, a, b, q, n_orbits = case
    rng = np.random.default_rng(seed)
    m_a = q / (1.0 + q)
    m_b = 1.0 / (1.0 + q)
    nu = m_a * m_b
    e0 = math.sqrt(1.0 - (b / a) ** 2)
    t0_period = 2.0 * math.pi * a ** 1.5

    # 初期条件（近点・遠点が a(1-e), a(1+e) になる接線速度）と軌道の積分（放射反作用あり）
    init, _, init_info = initial_state(a, b, nu, m_a, m_b)
    ts, st, dv, stop = integrate(a, b, nu, m_a, m_b, n_orbits * t0_period, rr=True, init=init)
    dt2 = t0_period / (OBS1_SAMPLES_PER_ORBIT * OBS2_OVERSAMPLING)
    tq = np.arange(0.0, ts[-1], dt2)
    y = interpolate(ts, st, dv, tq)
    x, yy, vx, vy, tau_a, tau_b = (y[:, i] for i in range(6))
    r = np.hypot(x, yy)
    v2 = vx * vx + vy * vy
    delta = m_a - m_b
    com = 0.5 * nu * delta * (v2 - 1.0 / r)
    xa, ya = (m_b + com) * x, (m_b + com) * yy
    xb, yb = (-m_a + com) * x, (-m_a + com) * yy
    samples = {"t": tq, "x": x, "y": yy, "vx": vx, "vy": vy}

    # 伏せる観測条件
    m_total_sun = rng.uniform(*M_TOTAL_SUN_RANGE)
    t_m = m_total_sun * T_SUN_S
    dist = math.exp(rng.uniform(math.log(DISTANCE_RANGE[0]), math.log(DISTANCE_RANGE[1])))
    phi = 2.0 * math.pi * rng.random()
    flip = bool(rng.random() < 0.5)
    off_x = rng.uniform(*SKY_OFFSET_RANGE)
    off_y = rng.uniform(*SKY_OFFSET_RANGE)
    t_off = rng.uniform(*TIME_OFFSET_RANGE)
    swap = bool(rng.random() < 0.5)

    c, s = math.cos(phi), math.sin(phi)
    fy = -1.0 if flip else 1.0

    def to_sky(px, py):
        py = fy * py
        return c * px - s * py, s * px + c * py

    scale = 1.0e9 / dist   # [nrad / M]
    ax_s, ay_s = to_sky(xa, ya)
    bx_s, by_s = to_sky(xb, yb)
    ax_s, ay_s = off_x + scale * ax_s, off_y + scale * ay_s
    bx_s, by_s = off_x + scale * bx_s, off_y + scale * by_s

    # 重力波（四重極公式、法線方向の観測者、天球の座標軸で）
    rx, ry = to_sky(x, yy)
    ux, uy = to_sky(vx, vy)
    nxs, nys = rx / r, ry / r
    amp = 4.0 * nu / dist
    h_xx = amp * (ux * ux - nxs * nxs / r)
    h_yy = amp * (uy * uy - nys * nys / r)
    h_xy = amp * (ux * uy - nxs * nys / r)
    h_plus = 0.5 * (h_xx - h_yy)
    h_cross = h_xy

    t_obs = t_off + t_m * tq

    # 観測誤差
    sig1 = OBS1_RESOLUTION * a * scale
    sig2 = OBS2_RESOLUTION * float(max(np.max(np.abs(h_plus)), np.max(np.abs(h_cross))))
    k1 = np.arange(0, tq.size, OBS2_OVERSAMPLING)
    n1 = k1.size
    first = (bx_s, by_s) if swap else (ax_s, ay_s)
    second = (ax_s, ay_s) if swap else (bx_s, by_s)
    obs1 = {
        "t_s": t_obs[k1],
        "body1_x_nrad": first[0][k1] + sig1 * rng.standard_normal(n1),
        "body1_y_nrad": first[1][k1] + sig1 * rng.standard_normal(n1),
        "body2_x_nrad": second[0][k1] + sig1 * rng.standard_normal(n1),
        "body2_y_nrad": second[1][k1] + sig1 * rng.standard_normal(n1),
    }
    obs2 = {
        "t_s": t_obs,
        "h_plus": h_plus + sig2 * rng.standard_normal(tq.size),
        "h_cross": h_cross + sig2 * rng.standard_normal(tq.size),
    }

    truth_ts = {
        "k": np.arange(tq.size),
        "in_obs1": (np.arange(tq.size) % OBS2_OVERSAMPLING == 0).astype(int),
        "t_M": tq,
        "t_s": t_obs,
        "x_rel": x, "y_rel": yy, "vx_rel": vx, "vy_rel": vy,
        "r": r,
        "phi_unwrapped": np.unwrap(np.arctan2(yy, x)),
        "xA": xa, "yA": ya, "xB": xb, "yB": yb,
        "tau_A": tau_a, "tau_B": tau_b,
        "h_plus_true": h_plus, "h_cross_true": h_cross,
        "energy_2pn_per_mu": pn_energy(x, yy, vx, vy, nu, order=2),
    }

    validation = validate(a, b, nu, m_a, m_b, t0_period, ts, st, dv, samples, init)

    meta = {
        "case_id": case_id,
        "case_label": label,
        "generator": "generate_binary_pn_observations.py",
        "generator_version": GENERATOR_VERSION,
        "seed": seed,
        "units": "G = c = M = 1 (M = m_A + m_B); lengths and times in units of M",
        "orbit_shape_parameters": {
            "a": a, "b": b,
            "definition": ("semi-major / semi-minor axes of the orbit at t = 0 in harmonic coordinates: "
                           "the conservative 2PN motion has its turning points at a(1-e) and a(1+e), "
                           "e = sqrt(1 - b^2/a^2); start at periapsis with the tangential speed found by bisection"),
            "initial_condition": init_info,
            "e": e0, "p": a * (1.0 - e0 * e0), "r_periapsis": a * (1.0 - e0),
            "newtonian_period_T0": t0_period,
        },
        "masses": {"q_mA_over_mB": q, "m_A": m_a, "m_B": m_b, "nu": nu,
                   "note": "A is the heavier body"},
        "dynamics": ("relative motion in harmonic coordinates: Newtonian + 1PN + 2PN + 2.5PN radiation "
                     "reaction (Kidder 1995, Eq. 2.2); individual positions from the 1PN centre-of-mass "
                     "relation; proper times from the 1PN metric of the companion only"),
        "integration": {"method": "RK4 with step = (2*pi/steps_per_orbit) * r^1.5",
                        "steps_per_orbit": STEPS_PER_ORBIT, "r_stop": R_STOP,
                        "orbits_requested": n_orbits, "t_end_M": float(ts[-1]),
                        "orbits_covered_in_T0": float(ts[-1] / t0_period), "stop_reason": stop},
        "observation_1": {
            "file_columns": list(obs1.keys()),
            "cadence_M": dt2 * OBS2_OVERSAMPLING,
            "n_samples": int(n1),
            "resolution_fraction_of_a": OBS1_RESOLUTION,
            "noise_sigma_nrad": sig1,
            "geometry": "face-on, observer far away; Doppler and light-travel-time differences neglected",
        },
        "observation_2": {
            "file_columns": list(obs2.keys()),
            "cadence_M": dt2,
            "n_samples": int(tq.size),
            "resolution_fraction_of_max_amplitude": OBS2_RESOLUTION,
            "noise_sigma_strain": sig2,
            "waveform": "leading-order quadrupole h_ij = (4 nu / D)(v_i v_j - n_i n_j / r), TT in the sky plane",
        },
        "hidden_observer_parameters": {
            "distance_D_in_M": dist,
            "angular_scale_nrad_per_M": scale,
            "total_mass_solar": m_total_sun,
            "T_M_seconds_per_M": t_m,
            "time_offset_s": t_off,
            "sky_rotation_rad": phi,
            "mirror_flip": flip,
            "sky_offset_nrad": [off_x, off_y],
            "body1_is": "B" if swap else "A",
            "body2_is": "A" if swap else "B",
            "sky_transform": "X = offset + scale * R(rotation) * diag(1, -1 if mirror else 1) * x",
        },
        "validation": validation,
    }
    return meta, obs1, obs2, truth_ts


def main():
    parser = argparse.ArgumentParser(description="強い場の2体の観測データを生成する（数値実験02）")
    parser.add_argument("--out", default=None, help="出力先（既定: このスクリプトと同じフォルダの data）")
    args = parser.parse_args()
    out = Path(args.out) if args.out else Path(__file__).resolve().parent / "data"
    (out / "input").mkdir(parents=True, exist_ok=True)
    (out / "truth").mkdir(parents=True, exist_ok=True)

    id_rng = np.random.default_rng(MASTER_SEED)
    order = np.argsort(id_rng.random(len(CASES)))
    ids = {int(j): f"D{rank + 1}" for rank, j in enumerate(order)}

    summaries = []
    for j, case in enumerate(CASES):
        did = ids[j]
        seed = MASTER_SEED * 10 + j
        meta, obs1, obs2, truth_ts = make_dataset(j, case, seed)
        meta["dataset_id"] = did
        write_csv(out / "input" / f"obs1_{did}.csv", obs1, "{:.12g}")
        write_csv(out / "input" / f"obs2_{did}.csv", obs2, "{:.12g}")
        write_csv(out / "truth" / f"{did}_timeseries.csv", truth_ts, "{:.15g}")
        write_json(out / "truth" / f"{did}_meta.json", meta)
        v = meta["validation"]
        summaries.append({
            "dataset_id": did, "case_id": case[0], "case_label": case[1],
            "a": case[2], "b": case[3], "q_mA_over_mB": case[4],
            "orbits_covered_in_T0": meta["integration"]["orbits_covered_in_T0"],
            "validation_passed": v["all_passed"],
            "obs1_file": f"input/obs1_{did}.csv", "obs2_file": f"input/obs2_{did}.csv",
            "meta_file": f"truth/{did}_meta.json", "timeseries_file": f"truth/{did}_timeseries.csv",
        })
        print(f"{did} {case[0]} a={case[2]:g} b={case[3]:g} q={case[4]:g}  "
              f"orbits={meta['integration']['orbits_covered_in_T0']:.2f}  "
              f"precession/6pi_p={v['precession_ratio_to_leading_order']:.3f}  "
              f"dE/Peters={v['energy_loss_ratio_to_peters']:.3f}  "
              f"conv={v['convergence_position_diff_over_a']:.2e}  passed={v['all_passed']}")

    scaling = pn_scaling_test()
    print(f"PN consistency (energy variation vs x): slope 1PN = {scaling['slope_1pn']:.2f}, "
          f"slope 2PN = {scaling['slope_2pn']:.2f}, passed = {scaling['passed']}")

    summaries.sort(key=lambda d: d["dataset_id"])
    manifest = {
        "experiment": "数値実験02 強い場の2体の観測データ（観測1 軌道観測・観測2 重力波観測）",
        "project": "匿名頂点状態生成幾何 ― 匿名内部観測者から不変な関係量の体系",
        "generator": "generate_binary_pn_observations.py",
        "generator_version": GENERATOR_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "master_seed": MASTER_SEED,
        "units": "G = c = M = 1 (M = m_A + m_B)",
        "cases": [{"case_id": c[0], "label": c[1], "a": c[2], "b": c[3], "q_mA_over_mB": c[4],
                   "orbits": c[5]} for c in CASES],
        "solver_input": "data/input/ only (truth/ must not be given to the analysis)",
        "pn_consistency_scaling_test": scaling,
        "all_datasets_passed_validation": all(d["validation_passed"] for d in summaries),
        "datasets": summaries,
    }
    write_json(out / "truth" / "manifest.json", manifest)
    print(f"all datasets passed validation: {manifest['all_datasets_passed_validation']}")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
