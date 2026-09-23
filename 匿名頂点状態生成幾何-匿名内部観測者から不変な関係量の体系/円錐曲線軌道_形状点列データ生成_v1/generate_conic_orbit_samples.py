#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数値実験01: 円錐曲線軌道の「形状のみ」点列データ生成
（匿名頂点状態生成幾何 ― 匿名内部観測者から不変な関係量の体系）

目的
    時刻情報を含まない軌道点列（形と並び順だけ）を、円・楕円・放物線・双曲線の
    各パターンで生成する。焦点の位置・長軸方向・進行方向は乱数で与え、解く側に
    形の種類・位置・向きの先入観を与えない。生成に使った厳密なパラメータと
    生成式を保存し、後から完全に再現・採点できるようにする。

生成式（焦点を原点、近点方向を +x とする標準座標 -> データ座標）
    r(θ) = p / (1 + e cos θ)
    q(θ) = r(θ) (cos θ, sin θ)
    P    = O + R(φ) q(θ),   R(φ) = [[cos φ, -sin φ], [sin φ, cos φ]]
    形の2自由度は (p, e)。e = 0 円、0 < e < 1 楕円、e = 1 放物線、e > 1 双曲線。
    a, b は有限な場合だけ派生値として保存する。

時刻（正解データ側のみ。解く側には渡さない）
    τ = sqrt(GM/p^3) (t - t_p) = ∫_0^θ dθ' / (1 + e cos θ')^2
    GM は形からは決まらないので、p を長さの単位とする無次元時刻で保存する。
    t_p は近点通過時刻。

出力（既定: このスクリプトと同じフォルダの data/）
    data/input/orbit_Dxx.csv           解く側に渡す点列（n, x, y）。渡すのはこれだけ
    data/truth/orbit_Dxx_meta.json     生成パラメータ・派生幾何量・生成式・並び順の扱い
    data/truth/orbit_Dxx_truth.csv     各点の正解（θ、近点角、時刻 τ など。行は入力と同じ順）
    data/truth/manifest.json           データセットID とケースの対応表、全体設定
    data/truth/validation_report.json  生成データの自己検証結果

使い方
    python generate_conic_orbit_samples.py
    python generate_conic_orbit_samples.py --out 出力先フォルダ
"""

import argparse
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

GENERATOR_VERSION = "1.0.0"
MASTER_SEED = 20260922

P_RANGE = (0.5, 3.0)            # 半直弦 p の範囲（任意の長さ単位）
OFFSET_RANGE = (-20.0, 20.0)    # データ座標での焦点 O の各成分の範囲
R_MAX_OVER_P = 10.0             # 開いた軌道・部分弧の観測窓: r <= R_MAX_OVER_P * p
H_MAX_OVER_P = 0.2              # 隣接点の弧長間隔の上限（円錐曲線の最小曲率半径 p に対する比）
MIN_POINTS = 250                # 短い曲線でも点数を確保するため、間隔の上限を L / MIN_POINTS 以下にもする
ARC_TABLE_SIZE = 400001         # 弧長 -> θ の逆写像に使う表の分割数
GL_PANELS = 400                 # 検算用数値積分（合成ガウス・ルジャンドル 20 点則）の区間数

# (ケースID, ラベル, 離心率 e, 範囲)   範囲: "full" = 1周、"window" = r <= R_MAX_OVER_P * p の弧
CASES = [
    ("C01", "circle",                     0.00, "full"),
    ("C02", "ellipse_near_circular",      0.02, "full"),
    ("C03", "ellipse",                    0.60, "full"),
    ("C04", "ellipse_high_e",             0.97, "full"),
    ("C05", "ellipse_high_e_partial_arc", 0.97, "window"),
    ("C06", "parabola",                   1.00, "window"),
    ("C07", "hyperbola_near_parabolic",   1.03, "window"),
    ("C08", "hyperbola",                  1.80, "window"),
]

# 刻み方。どちらも時間とは無関係な、弧長を基準にした不規則な刻み。
# 揺らぎの幅は、連続する間隔が g[k+1] < g[k] + g[k-1] を満たすように選んである
# （最近傍を順にたどる並べ替えが、閉曲線ではどの点から始めても順序を復元するための条件）。
SAMPLING_MODES = {
    "A": {"name": "arc_jitter",
          "rule": "gap = h_max * U(0.55, 1.0)"},
    "B": {"name": "arc_cluster",
          "rule": "gap = h_max * exp(beta * (sin(2*pi*k*s/L + psi) - 1)) * U(0.7, 1.0)",
          "beta": 0.9,
          "k_choices": [1, 2, 3]},
}

GENERATION_EQUATION = (
    "P = O + R(phi) q(theta);  q(theta) = r(theta) (cos theta, sin theta);  "
    "r(theta) = p / (1 + e cos theta);  R(phi) = [[cos phi, -sin phi], [sin phi, cos phi]]"
)
TIME_DEFINITION = (
    "tau = sqrt(GM/p^3) (t - t_p) = integral_0^theta dtheta' / (1 + e cos theta')^2;  "
    "t_p = periapsis passage.  GM is not determined by the shape, so time is dimensionless "
    "in units of sqrt(p^3/GM)."
)
ANOMALY_TYPE = {
    "circle": "E: eccentric anomaly (equals theta for a circle)",
    "ellipse": "E: eccentric anomaly",
    "parabola": "D = tan(theta/2)",
    "hyperbola": "H: hyperbolic anomaly",
}

_GL_X, _GL_W = np.polynomial.legendre.leggauss(20)
_U = ((np.arange(GL_PANELS)[:, None] + 0.5 * (_GL_X[None, :] + 1.0)) / GL_PANELS).ravel()
_WU = np.tile(0.5 * _GL_W / GL_PANELS, GL_PANELS)


# ---------------------------------------------------------------- 円錐曲線の基本量

def pattern_of(e):
    if e == 0.0:
        return "circle"
    if e < 1.0:
        return "ellipse"
    if e == 1.0:
        return "parabola"
    return "hyperbola"


def conic_r(theta, p, e):
    return p / (1.0 + e * np.cos(theta))


def arc_density(theta, p, e):
    """ds/dθ = p sqrt(1 + 2e cos θ + e^2) / (1 + e cos θ)^2"""
    c = np.cos(theta)
    return p * np.sqrt(1.0 + 2.0 * e * c + e * e) / (1.0 + e * c) ** 2


def time_density(theta, e):
    """dτ/dθ = 1 / (1 + e cos θ)^2（面積速度一定を無次元化したもの）"""
    return 1.0 / (1.0 + e * np.cos(theta)) ** 2


def integrate_from_zero(func, theta):
    """各 θ について ∫_0^θ func(θ') dθ' を合成ガウス・ルジャンドル則で計算する（θ < 0 も可）。"""
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    out = np.empty_like(theta)
    for i, th in enumerate(theta):
        out[i] = th * np.dot(_WU, func(th * _U))
    return out


def kepler_time_closed_form(theta, e):
    """近点通過からの無次元時刻 K(θ) = ∫_0^θ dθ'/(1 + e cos θ')^2 を解析式で返す。

    θ は閉曲線なら (-π, π]、開いた曲線なら |θ| < θ∞ で与える。
    戻り値 (anomaly, K)。anomaly は 円・楕円: 離心近点角 E、放物線: D = tan(θ/2)、
    双曲線: 双曲線近点角 H。
    """
    theta = np.asarray(theta, dtype=float)
    kind = pattern_of(e)
    if kind == "circle":
        return theta.copy(), theta.copy()
    if kind == "ellipse":
        beta = e / (1.0 + math.sqrt(1.0 - e * e))
        E = theta - 2.0 * np.arctan(beta * np.sin(theta) / (1.0 + beta * np.cos(theta)))
        M = E - e * np.sin(E)
        return E, M / (1.0 - e * e) ** 1.5
    if kind == "parabola":
        D = np.tan(0.5 * theta)
        return D, 0.5 * (D + D ** 3 / 3.0)
    k = math.sqrt((e - 1.0) / (e + 1.0))
    H = 2.0 * np.arctanh(k * np.tan(0.5 * theta))
    return H, (e * np.sinh(H) - H) / (e * e - 1.0) ** 1.5


def stumpff_S(z):
    """Stumpff 関数 S(z) = Σ (-z)^k / (2k+3)!"""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.empty_like(z)
    small = np.abs(z) < 0.1
    zs = z[small]
    term = np.full_like(zs, 1.0 / 6.0)
    acc = term.copy()
    for k in range(1, 10):
        term = term * (-zs) / ((2 * k + 2) * (2 * k + 3))
        acc = acc + term
    out[small] = acc
    pos = z >= 0.1
    sq = np.sqrt(z[pos])
    out[pos] = (sq - np.sin(sq)) / sq ** 3
    neg = z <= -0.1
    sq = np.sqrt(-z[neg])
    out[neg] = (np.sinh(sq) - sq) / sq ** 3
    return out


def universal_time(anomaly, e):
    """普遍変数による無次元時刻（p = 1 単位）: τ = r_p χ + e χ^3 S(α χ^2)、α = 1 - e^2"""
    kind = pattern_of(e)
    alpha = 1.0 - e * e
    rp = 1.0 / (1.0 + e)
    if kind in ("circle", "ellipse"):
        chi = anomaly / math.sqrt(alpha)
    elif kind == "parabola":
        chi = anomaly
    else:
        chi = anomaly / math.sqrt(-alpha)
    return rp * chi + e * chi ** 3 * stumpff_S(alpha * chi ** 2)


def theta_range(e, span):
    if span == "full":
        return 0.0, 2.0 * math.pi
    th = math.acos((1.0 / R_MAX_OVER_P - 1.0) / e)
    return -th, th


def arc_table(p, e, th0, th1):
    grid = np.linspace(th0, th1, ARC_TABLE_SIZE)
    f = arc_density(grid, p, e)
    s = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(grid))))
    return grid, s


def rotate(phi, xc, yc):
    c, s = math.cos(phi), math.sin(phi)
    return c * xc - s * yc, s * xc + c * yc


# ---------------------------------------------------------------- サンプリング

def draw_gaps(L, h_max, mode, rng):
    """弧長の間隔列を作り、合計がちょうど L になるよう縮める（縮尺 <= 1 なので上限 h_max は保たれる）。"""
    params = {}
    if mode == "B":
        spec = SAMPLING_MODES["B"]
        params = {"beta": spec["beta"],
                  "k": spec["k_choices"][int(3 * rng.random())],
                  "psi": 2.0 * math.pi * rng.random()}
    gaps = []
    s = 0.0
    while s < L:
        u = rng.random()
        if mode == "A":
            g = h_max * (0.55 + 0.45 * u)
        else:
            mod = math.exp(params["beta"] * (math.sin(2.0 * math.pi * params["k"] * s / L + params["psi"]) - 1.0))
            g = h_max * mod * (0.7 + 0.3 * u)
        gaps.append(g)
        s += g
    gaps = np.asarray(gaps)
    scale = L / gaps.sum()
    return gaps * scale, params, scale


# ---------------------------------------------------------------- 派生幾何量

def derived_geometry(p, e, O, phi, th0, th1):
    kind = pattern_of(e)

    def to_data(xc, yc):
        dx, dy = rotate(phi, xc, yc)
        return [O[0] + dx, O[1] + dy]

    rp = p / (1.0 + e)
    g = {
        "pattern": kind,
        "p_semi_latus_rectum": p,
        "e_eccentricity": e,
        "a": None,
        "b": None,
        "r_periapsis": rp,
        "attracting_focus_xy": [O[0], O[1]],
        "periapsis_direction_rad": phi,
        "periapsis_point_xy": to_data(rp, 0.0),
    }
    if kind in ("circle", "ellipse"):
        a = p / (1.0 - e * e)
        b = p / math.sqrt(1.0 - e * e)
        ra = p / (1.0 - e)
        g.update({
            "a": a,
            "b": b,
            "r_apoapsis": ra,
            "center_xy": to_data(-a * e, 0.0),
            "empty_focus_xy": to_data(-2.0 * a * e, 0.0),
            "apoapsis_point_xy": to_data(-ra, 0.0),
            "period_tau": 2.0 * math.pi / (1.0 - e * e) ** 1.5,
        })
        if kind == "circle":
            g["note"] = "circle: focus = center; the periapsis direction phi is not identifiable from the shape"
    elif kind == "parabola":
        g.update({"focal_length": 0.5 * p, "vertex_xy": to_data(rp, 0.0)})
    else:
        a = p / (e * e - 1.0)
        b = p / math.sqrt(e * e - 1.0)
        th_inf = math.acos(-1.0 / e)
        g.update({
            "a": a,
            "b": b,
            "center_xy": to_data(a * e, 0.0),
            "empty_focus_xy": to_data(2.0 * a * e, 0.0),
            "theta_infinity_rad": th_inf,
            "asymptote_direction_angles_rad": [(phi + th_inf) % (2.0 * math.pi),
                                               (phi - th_inf) % (2.0 * math.pi)],
            "deflection_angle_rad": 2.0 * math.asin(1.0 / e),
        })
    # データ座標での一般形 A x^2 + B xy + C y^2 + D x + E y + F = 0
    c, s = math.cos(phi), math.sin(phi)
    A = 1.0 - e * e * c * c
    B = -2.0 * e * e * c * s
    C = 1.0 - e * e * s * s
    D0, E0, F0 = 2.0 * p * e * c, 2.0 * p * e * s, -p * p
    Ox, Oy = O
    D = D0 - 2.0 * A * Ox - B * Oy
    E = E0 - 2.0 * C * Oy - B * Ox
    F = A * Ox * Ox + B * Ox * Oy + C * Oy * Oy - D0 * Ox - E0 * Oy + F0
    g["general_conic_equation"] = "A x^2 + B x y + C y^2 + D x + E y + F = 0"
    g["general_conic_coefficients_ABCDEF"] = [A, B, C, D, E, F]
    g["discriminant_B2_minus_4AC"] = B * B - 4.0 * A * C
    g["theta_range_rad"] = [th0, th1]
    return g


# ---------------------------------------------------------------- 自己検証

def greedy_chain(D2, start):
    n = D2.shape[0]
    visited = np.zeros(n, dtype=bool)
    order = np.empty(n, dtype=int)
    cur = start
    visited[cur] = True
    order[0] = cur
    for k in range(1, n):
        d = np.where(visited, np.inf, D2[cur])
        cur = int(np.argmin(d))
        visited[cur] = True
        order[k] = cur
    return order


def validate_dataset(v):
    x, y, O, phi, p, e = v["x"], v["y"], v["O"], v["phi"], v["p"], v["e"]
    closed = v["closed"]
    res = {}

    # (1) 生成式の再現: 正解のパラメータで標準座標に戻し、円錐曲線の式を満たすか
    dx, dy = x - O[0], y - O[1]
    c, s = math.cos(phi), math.sin(phi)
    xc = c * dx + s * dy
    yc = -s * dx + c * dy
    r = np.hypot(xc, yc)
    res["conic_residual_max_over_p"] = float(np.max(np.abs(r - (p - e * xc))) / p)
    dth = np.angle(np.exp(1j * (np.arctan2(yc, xc) - v["theta"])))
    res["theta_recovery_max_abs_rad"] = float(np.max(np.abs(dth)))

    # (2) 最近傍による並び順の復元可能性（入力の順で評価）
    io = v["input_order"]
    X, Y, S = x[io], y[io], v["s_gen"][io]
    n = X.size
    D2 = (X[:, None] - X[None, :]) ** 2 + (Y[:, None] - Y[None, :]) ** 2
    np.fill_diagonal(D2, np.inf)
    rows = np.arange(n)
    #   実際の並べ替え: 閉曲線は8か所の始点から、開いた曲線は両端から最近傍をたどる
    if closed:
        starts = sorted({int(round(k * n / 8)) % n for k in range(8)})
        chain_ok = True
        for st in starts:
            ch = greedy_chain(D2, st)
            chain_ok &= (np.array_equal(ch, (st + rows) % n) or np.array_equal(ch, (st - rows) % n))
    else:
        starts = [0, n - 1]
        chain_ok = (np.array_equal(greedy_chain(D2, 0), rows) and
                    np.array_equal(greedy_chain(D2, n - 1), rows[::-1]))
    res["greedy_chain_starts_tested"] = starts
    res["greedy_nearest_neighbour_chain_recovers_order"] = bool(chain_ok)
    #   閉曲線で「どの点から始めても」復元できるための局所条件:
    #   d(i, i+1) < d(i, i-2) かつ d(i, i-1) < d(i, i+2)（始点のすぐ後ろへ戻る誤りが起きない）
    if closed:
        d_next = D2[rows, (rows + 1) % n]
        d_prev = D2[rows, (rows - 1) % n]
        local_viol = int(np.sum(d_next >= D2[rows, (rows - 2) % n]) +
                         np.sum(d_prev >= D2[rows, (rows + 2) % n]))
    else:
        local_viol = None
    res["local_any_start_condition_violations"] = local_viol

    # (3) 隣接点の間隔
    seg = np.hypot(np.diff(X), np.diff(Y))
    if closed:
        seg = np.append(seg, math.hypot(X[0] - X[-1], Y[0] - Y[-1]))
    res["chord_gap_max_over_p"] = float(seg.max() / p)
    res["chord_gap_min_over_p"] = float(seg.min() / p)
    res["chord_gap_mean_over_p"] = float(seg.mean() / p)
    res["arc_gap_upper_bound_over_p"] = float(v["h_max"] / p)
    #   曲線の離れた部分との余裕: 曲線に沿って p 以上離れた点までの最短距離 / 隣接点の最大間隔
    #   （1 を大きく超えていれば、最近傍をたどるときに曲線の別の部分へ飛び移ることはない）
    ds = np.abs(S[:, None] - S[None, :])
    if closed:
        ds = np.minimum(ds, v["L"] - ds)
    res["across_curve_clearance_ratio"] = float(math.sqrt(np.min(D2[ds >= p])) / seg.max())

    # (4) 時刻の検算: 解析式 vs 数値積分、解析式 vs 普遍変数
    K = v["K"]
    scale = max(1.0, float(np.max(np.abs(K))))
    res["tau_closed_form_vs_quadrature_max_rel"] = float(np.max(np.abs(K - v["K_quad"])) / scale)
    res["tau_closed_form_vs_universal_variable_max_rel"] = float(np.max(np.abs(K - v["K_uv"])) / scale)
    tau_t = v["tau_elapsed_t"]
    if closed:
        total = float(tau_t[-1] + v["closing"])
        res["period_tau_expected"] = v["T_tau"]
        res["period_tau_from_samples"] = total
        res["period_rel_error"] = float(abs(total - v["T_tau"]) / v["T_tau"])

    # (5) 時刻の刻みの不規則さ（点の密度に時間情報が焼き込まれていないことの目安）
    dt = np.diff(tau_t)
    if closed:
        dt = np.append(dt, v["closing"])
    res["dtau_coefficient_of_variation"] = float(dt.std() / dt.mean())
    res["dtau_max_over_min"] = float(dt.max() / dt.min())

    # (6) 進行方向: 引力中心まわりの角運動量の符号
    to = v["time_order"]
    xt, yt = x[to] - O[0], y[to] - O[1]
    lz = float(np.sum(xt[:-1] * np.diff(yt) - yt[:-1] * np.diff(xt)))
    res["motion_direction_check"] = bool(np.sign(lz) == v["direction"])

    # (7) 観測窓（開いた曲線・部分弧）
    if not closed:
        res["r_max_over_p"] = float(r.max() / p)

    checks = {
        "conic_residual": res["conic_residual_max_over_p"] < 1e-11,
        "theta_recovery": res["theta_recovery_max_abs_rad"] < 1e-10,
        "nearest_neighbour_order": bool(chain_ok) and (local_viol in (None, 0)),
        "across_curve_clearance": res["across_curve_clearance_ratio"] > 1.0,
        "gap_bound": res["chord_gap_max_over_p"] <= res["arc_gap_upper_bound_over_p"] * (1.0 + 1e-6),
        "time_quadrature": res["tau_closed_form_vs_quadrature_max_rel"] < 1e-10,
        "time_universal_variable": res["tau_closed_form_vs_universal_variable_max_rel"] < 1e-10,
        "motion_direction": res["motion_direction_check"],
    }
    if closed:
        checks["period"] = res["period_rel_error"] < 1e-11
    else:
        checks["observation_window"] = res["r_max_over_p"] <= R_MAX_OVER_P * (1.0 + 1e-12)
    res["checks"] = checks
    res["all_passed"] = bool(all(checks.values()))
    return res


# ---------------------------------------------------------------- データセット1つの生成

def make_dataset(job_index, case, mode, seed, direction, input_reversed):
    """direction（進行方向の符号）と input_reversed（入力の並びを進行方向に対して反転するか）は
    データセット全体で半数ずつになるよう main で割り当てる。"""
    case_id, label, e, span = case
    rng = np.random.default_rng(seed)

    # 形の大きさ・焦点の位置・近点方向
    p = P_RANGE[0] + (P_RANGE[1] - P_RANGE[0]) * rng.random()
    O = (OFFSET_RANGE[0] + (OFFSET_RANGE[1] - OFFSET_RANGE[0]) * rng.random(),
         OFFSET_RANGE[0] + (OFFSET_RANGE[1] - OFFSET_RANGE[0]) * rng.random())
    phi = 2.0 * math.pi * rng.random()

    # 弧長を基準にした不規則な刻みで θ を決める
    closed = (span == "full")
    th0, th1 = theta_range(e, span)
    grid, s_tab = arc_table(p, e, th0, th1)
    L = float(s_tab[-1])
    h_max = min(H_MAX_OVER_P * p, L / MIN_POINTS)
    gaps, mode_params, gap_scale = draw_gaps(L, h_max, mode, rng)
    if closed:
        s_start = L * rng.random()
        s_gen = np.mod(s_start + np.concatenate(([0.0], np.cumsum(gaps[:-1]))), L)
    else:
        s_start = 0.0
        s_gen = np.concatenate(([0.0], np.cumsum(gaps)))
        s_gen[-1] = L
    theta = np.interp(s_gen, s_tab, grid)          # 生成順 = θ の増える向き
    n_pts = theta.size

    # 並び順: 生成順 -> 時間順（進行方向） -> 入力順（向きを乱数で反転、閉曲線は始点もずらす）
    gen_index = np.arange(n_pts)
    time_order = gen_index.copy() if direction == 1 else gen_index[::-1].copy()
    shift = int(n_pts * rng.random()) if closed else 0
    input_order = np.roll(time_order, -shift) if closed else time_order.copy()
    if input_reversed:
        input_order = input_order[::-1].copy()

    # 点の座標（生成式そのもの）
    r = conic_r(theta, p, e)
    xc = r * np.cos(theta)
    yc = r * np.sin(theta)
    dx, dy = rotate(phi, xc, yc)
    x = O[0] + dx
    y = O[1] + dy

    # 正解の時刻
    theta_pm = np.where(theta > math.pi, theta - 2.0 * math.pi, theta) if closed else theta
    anomaly, K = kepler_time_closed_form(theta_pm, e)
    tau_sp = direction * K
    tau_t = tau_sp[time_order]
    if closed:
        T_tau = 2.0 * math.pi / (1.0 - e * e) ** 1.5
        d = np.diff(tau_t)
        d = np.where(d < 0.0, d + T_tau, d)
        closing = float((tau_t[0] - tau_t[-1]) % T_tau)
    else:
        T_tau = None
        d = np.diff(tau_t)
        closing = None
    tau_elapsed_t = np.concatenate(([0.0], np.cumsum(d)))
    tau_elapsed = np.empty(n_pts)
    tau_elapsed[time_order] = tau_elapsed_t
    time_rank = np.empty(n_pts, dtype=int)
    time_rank[time_order] = np.arange(n_pts)

    arc_peri = integrate_from_zero(lambda t: arc_density(t, p, e), theta_pm)
    K_quad = integrate_from_zero(lambda t: time_density(t, e), theta_pm)
    K_uv = universal_time(anomaly, e)

    validation = validate_dataset({
        "x": x, "y": y, "O": O, "phi": phi, "p": p, "e": e, "theta": theta, "closed": closed,
        "input_order": input_order, "time_order": time_order, "direction": direction,
        "K": K, "K_quad": K_quad, "K_uv": K_uv, "tau_elapsed_t": tau_elapsed_t,
        "closing": closing, "T_tau": T_tau, "h_max": h_max, "s_gen": s_gen, "L": L,
    })

    meta = {
        "dataset_id": None,
        "case_id": case_id,
        "case_label": label,
        "pattern": pattern_of(e),
        "generator": "generate_conic_orbit_samples.py",
        "generator_version": GENERATOR_VERSION,
        "seed": seed,
        "job_index": job_index,
        "shape_parameters": {"p": p, "e": e},
        "pose": {"focus_offset_O_xy": [O[0], O[1]],
                 "rotation_phi_rad": phi,
                 "rotation_phi_deg": math.degrees(phi)},
        "motion_direction": "counterclockwise" if direction == 1 else "clockwise",
        "motion_direction_sign": direction,
        "closed_curve": closed,
        "span": "full revolution" if closed else f"window r <= {R_MAX_OVER_P:g} p",
        "sampling": {
            "mode": mode,
            **SAMPLING_MODES[mode],
            "mode_parameters": mode_params,
            "h_max": h_max,
            "h_max_over_p": h_max / p,
            "gap_scale_factor": gap_scale,
            "arc_length_total": L,
            "arc_start_offset_from_periapsis": s_start if closed else None,
            "arc_measured_from": "periapsis (theta = 0)" if closed else "theta = -theta_max end",
            "n_points": n_pts,
        },
        "input_order": {
            "description": ("input rows carry only what nearest-neighbour sorting can recover: "
                            "orientation reversed at random relative to the motion; "
                            "closed curves also start at a random point"),
            "reversed_relative_to_motion": input_reversed,
            "cyclic_shift": shift,
        },
        "generation_equation": GENERATION_EQUATION,
        "time_definition": TIME_DEFINITION,
        "anomaly_type": ANOMALY_TYPE[pattern_of(e)],
        "derived_geometry": derived_geometry(p, e, O, phi, th0, th1),
        "truth_columns": {
            "n": "row index of the input file",
            "time_rank": "true order in time (0 = first)",
            "theta_rad": "true anomaly theta (closed curves: (-pi, pi])",
            "anomaly": "see anomaly_type",
            "tau_since_periapsis": "dimensionless time since periapsis passage",
            "tau_elapsed": "dimensionless time since the first point in time order",
            "r": "distance from the attracting focus",
            "x_canonical": "x in the canonical frame (focus at origin, periapsis on +x)",
            "y_canonical": "y in the canonical frame",
            "arc_from_periapsis": "signed arc length from periapsis (positive for theta > 0)",
        },
    }

    io = input_order
    truth_cols = {
        "n": np.arange(n_pts),
        "time_rank": time_rank[io],
        "theta_rad": theta_pm[io],
        "anomaly": anomaly[io],
        "tau_since_periapsis": tau_sp[io],
        "tau_elapsed": tau_elapsed[io],
        "r": r[io],
        "x_canonical": xc[io],
        "y_canonical": yc[io],
        "arc_from_periapsis": arc_peri[io],
    }
    return {"meta": meta, "x": x[io], "y": y[io], "truth": truth_cols, "validation": validation}


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


def write_csv(path, columns):
    names = list(columns.keys())
    fmts = ["{:d}" if np.issubdtype(np.asarray(columns[k]).dtype, np.integer) else "{:.17g}" for k in names]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(names) + "\n")
        for row in zip(*(columns[k] for k in names)):
            f.write(",".join(fmt.format(v) for fmt, v in zip(fmts, row)) + "\n")


def main():
    parser = argparse.ArgumentParser(description="円錐曲線軌道の形状点列データを生成する（数値実験01）")
    parser.add_argument("--out", default=None,
                        help="出力先フォルダ（既定: このスクリプトと同じフォルダの data）")
    args = parser.parse_args()
    out = Path(args.out) if args.out else Path(__file__).resolve().parent / "data"
    in_dir, tr_dir = out / "input", out / "truth"
    in_dir.mkdir(parents=True, exist_ok=True)
    tr_dir.mkdir(parents=True, exist_ok=True)

    jobs = [(case, mode) for case in CASES for mode in ("A", "B")]

    # データセットID はケースと無関係な乱数順に振る（ファイル名から種類が分からないように）
    id_rng = np.random.default_rng(MASTER_SEED)
    keys = id_rng.random(len(jobs))
    rank_of_job = np.empty(len(jobs), dtype=int)
    rank_of_job[np.argsort(keys)] = np.arange(len(jobs))
    # 進行方向と入力の並びの反転は、それぞれちょうど半数ずつを乱数で割り当てる
    half = len(jobs) // 2
    direction_of_job = np.ones(len(jobs), dtype=int)
    direction_of_job[np.argsort(id_rng.random(len(jobs)))[:half]] = -1
    reversed_of_job = np.zeros(len(jobs), dtype=bool)
    reversed_of_job[np.argsort(id_rng.random(len(jobs)))[:half]] = True

    summaries, report = [], {}
    for j, (case, mode) in enumerate(jobs):
        did = f"D{rank_of_job[j] + 1:02d}"
        seed = MASTER_SEED * 1000 + j
        ds = make_dataset(j, case, mode, seed, int(direction_of_job[j]), bool(reversed_of_job[j]))
        ds["meta"]["dataset_id"] = did

        f_in = f"input/orbit_{did}.csv"
        f_meta = f"truth/orbit_{did}_meta.json"
        f_truth = f"truth/orbit_{did}_truth.csv"
        write_csv(out / f_in, {"n": np.arange(ds["x"].size), "x": ds["x"], "y": ds["y"]})
        write_csv(out / f_truth, ds["truth"])
        write_json(out / f_meta, ds["meta"])

        report[did] = ds["validation"]
        summaries.append({
            "dataset_id": did,
            "case_id": case[0],
            "case_label": case[1],
            "pattern": ds["meta"]["pattern"],
            "e": case[2],
            "span": ds["meta"]["span"],
            "sampling_mode": mode,
            "n_points": int(ds["x"].size),
            "validation_passed": ds["validation"]["all_passed"],
            "input_file": f_in,
            "meta_file": f_meta,
            "truth_file": f_truth,
        })
        print(f"{did}  {case[0]} {case[1]:<28s} e={case[2]:<5.2f} mode={mode}  "
              f"N={ds['x'].size:4d}  passed={ds['validation']['all_passed']}")

    summaries.sort(key=lambda d: d["dataset_id"])
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    all_passed = all(v["all_passed"] for v in report.values())
    manifest = {
        "experiment": "数値実験01 円錐曲線軌道の形状点列データ（形と並び順のみ）",
        "project": "匿名頂点状態生成幾何 ― 匿名内部観測者から不変な関係量の体系",
        "generator": "generate_conic_orbit_samples.py",
        "generator_version": GENERATOR_VERSION,
        "generated_at_utc": now,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "master_seed": MASTER_SEED,
        "per_dataset_seed_rule": "seed = master_seed * 1000 + job_index",
        "settings": {
            "p_range": P_RANGE,
            "offset_range_each_axis": OFFSET_RANGE,
            "rotation_phi_range_rad": [0.0, 2.0 * math.pi],
            "motion_direction": "counterclockwise / clockwise, each assigned to exactly half of the datasets at random",
            "input_order_reversal": "input rows reversed relative to the motion in exactly half of the datasets (assigned at random)",
            "window_r_max_over_p": R_MAX_OVER_P,
            "arc_gap_upper_bound_over_p": H_MAX_OVER_P,
            "min_points_rule": f"arc gap upper bound <= L / {MIN_POINTS}",
        },
        "generation_equation": GENERATION_EQUATION,
        "time_definition": TIME_DEFINITION,
        "cases": [{"case_id": c[0], "case_label": c[1], "e": c[2], "span": c[3]} for c in CASES],
        "sampling_modes": SAMPLING_MODES,
        "solver_input": "data/input/ only (the truth/ folder must not be given to the solver)",
        "all_datasets_passed_validation": all_passed,
        "datasets": summaries,
    }
    write_json(tr_dir / "manifest.json", manifest)
    write_json(tr_dir / "validation_report.json", {
        "generated_at_utc": now,
        "all_datasets_passed": all_passed,
        "datasets": {k: report[k] for k in sorted(report)},
    })
    print(f"all datasets passed validation: {all_passed}")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
