#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_coulomb_readout_experiments_v1.py

第二思考実験（二つの値と線形の法則のまま、位置を積の読み出し、時計を面積保存から決めると
クーロン型の逆二乗の力が出る）の数値検証 E1〜E7。

約束（設計書 §12、第一思考実験と同じ）:
  * 法則は第一思考実験のまま X_{k+1} = S X_k, det S = 1。状態更新に丸め・クリッピング・正規化を入れない。
  * 落ちる条件を FAIL_CONDITIONS に事前に宣言し、実行結果と照合して PASS / FAIL を出す。
  * 恒等式は有理数の厳密計算でも確かめる。乱数の種は固定。パスは相対。
  * 出力: coulomb_readout_experiments_results_v1.json（このスクリプトと同じフォルダ）

読み出し（P: 入れ替え、D: 片方の符号反転、J0 = D P: 四分の一回転）:
  x = ω(X,PX) = a²−b²,  y = ω(X,DX) = 2ab,  r = ω(X,J0 X) = a²+b²
時計（読み出し平面で一歩に掃く面積が時計あたり一定、から決まる）:
  Δt_k = (2Δs/τ) · X_k·X_{k+1},   Δs = arccos(τ/2)（楕円型）, arccosh(|τ|/2)（双曲型）

依存: numpy。実行: python3 run_coulomb_readout_experiments_v1.py
"""
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS_JSON = HERE / "coulomb_readout_experiments_results_v1.json"
SEED = 20260921

FAIL_CONDITIONS = {
    "E1a": "有理数厳密計算で x²+y²−r² が一度でも 0 にならない。",
    "E1b": "有理数厳密計算で、第一論文の保存量 q=XᵀGX が読み出しの一次式 ((A−C)/2)x+By+((A+C)/2)r と一度でも食い違う。",
    "E2a": "計量から出した離心率 e²=1−4detG/(trG)² が、Kepler の式 e²=1+2EL²/μ² と 1e-9 を超えてずれる。",
    "E2b": "型と離心率の対応（楕円型⇔e<1、双曲型⇔e>1）が一件でも破れる。",
    "E3": "補間した軌道で、時計 dt=r ds のもとの加速度が −μz/r³ から相対 1e-9 を超えてずれる（楕円型・双曲型、引力・斥力）。",
    "E4a": "離散の時計で、一歩に掃く三角形の面積÷Δt_k が軌道に沿って相対 1e-12 を超えて変わる。",
    "E4b": "閉軌道の離散の時計の一周の和が Kepler の周期 πa の 2 倍から相対 1e-12 を超えてずれる。",
    "E4c": "「時計が単調（すべての Δt_k>0）⇔ cosθ>e」が一件でも破れる。",
    "E4c_prime": "【E4c が落ちた後に追加した診断。事前登録ではない】(i) cosθ>e なのに負の刻みが出る、または (ii) cosθ<e なのに、初期位相を細かく走査しても負の刻みが一度も出ない。",
    "E4d": "離散の時計での差分加速度と −μz/r³ の誤差の、歩幅 θ に対する次数が [1.8, 2.2] を外れる。",
    "E5a": "楕円型で μ≤0 が一件でも出る。",
    "E5b": "双曲型で μ が −8αβ(v₊·v₋) から相対 1e-9 を超えてずれる、または引力・斥力のどちらかが一件も出ない。",
    "E5d": "双曲型で、配置が入っている区画（零方向 v₊, v₋ の向きを α>0, β>0 となるよう取ったときの二方向のあいだ）の開き角 ψ について、「引力 ⇔ ψ>90°」が一件でも破れる。",
    "E5c": "厳密に無名な法則（ブースト）と圧搾の正規形で |μ| または読み出し点の加速度が 1e-12 を超える。",
    "E6": "q=0 の軌道（X_0∥X_1）で、q_k が厳密に 0 でない、読み出し点が原点からの一本の半直線を外れる、または値が有限でなくなる。",
    "E6b": "【E6 の後に追加した検証。事前登録ではない】τ=149/70=λ+1/λ（λ=7/10）で、(i) D2 の範囲の q=0 の軌道（固有軌道）で c_k が λ^k にならない、q_k が厳密に 0 でない、読み出しが半直線を外れる、または原点通過が 0 回でない、(ii) 二階則の一般解（係数が異符号。S の軌道ではない）で原点通過がちょうど 1 回でない。",
    "E7a": "記述の回転 X→RX, S→RSR⁻¹ で、e・歩ごとの角度差・距離の比が 1e-10 を超えて変わる。",
    "E7b": "入れ替え P で (x,y,r)→(−x,y,r) にならない、または X→−X で読み出しが変わる。",
}

P = np.array([[0.0, 1.0], [1.0, 0.0]])
D = np.array([[-1.0, 0.0], [0.0, 1.0]])
J0 = D @ P
OM = np.array([[0.0, 1.0], [-1.0, 0.0]])


def om(u, v):
    return float(u[0] * v[1] - u[1] * v[0])


def readout(X):
    """三つの面積読み出し (x, y, r)。"""
    return np.array([om(X, P @ X), om(X, D @ X), om(X, J0 @ X)])


def random_sl2(rng):
    while True:
        A = rng.standard_normal((2, 2))
        d = np.linalg.det(A)
        if abs(d) < 1e-2:
            continue
        if d < 0:
            A = A[::-1].copy()
            d = -d
        return A / math.sqrt(d)


def law(S):
    """(型, Δs, A)。A=(S−S⁻¹)/(2 sinΔs または 2 sinhΔs)。X'=AX, X''=∓X（s の単位で角振動数 1）。"""
    tau = float(np.trace(S))
    N = (S - np.linalg.inv(S)) / 2.0
    if abs(tau) < 2:
        ds = math.acos(tau / 2)
        return "elliptic", ds, N / math.sin(ds)
    ds = math.acosh(abs(tau) / 2)
    return "hyperbolic", ds, N / math.sinh(ds)


def kepler(X, dX, kind):
    """(E, μ, L)。E=∓2, μ=2|X'|²−E|X|², L=2ω(X,X')。"""
    E = -2.0 if kind == "elliptic" else 2.0
    return E, 2 * float(dX @ dX) - E * float(X @ X), 2 * om(X, dX)


# --------------------------------------------------------------------------- E1
def e1_exact(rng, n=300, steps=10):
    bad_null = bad_plane = 0
    F = lambda: Fraction(int(rng.integers(-6, 7)), int(rng.integers(1, 6)))
    done = 0
    while done < n:
        p, q_, r_ = F(), F(), F()
        if p == 0:
            continue
        S = [[p, q_], [r_, (1 + q_ * r_) / p]]
        Si = [[S[1][1], -S[0][1]], [-S[1][0], S[0][0]]]
        Nm = [[(S[i][j] - Si[i][j]) / 2 for j in range(2)] for i in range(2)]
        G = [[Nm[1][0], Nm[1][1]], [-Nm[0][0], -Nm[0][1]]]          # G = Ω N
        A_, B_, C_ = G[0][0], G[0][1], G[1][1]
        X = [F(), F()]
        if X == [0, 0]:
            continue
        for _ in range(steps):
            a, b = X
            x, y, r = a * a - b * b, 2 * a * b, a * a + b * b
            if x * x + y * y - r * r != 0:
                bad_null += 1
            qv = A_ * a * a + 2 * B_ * a * b + C_ * b * b
            if qv != (A_ - C_) / 2 * x + B_ * y + (A_ + C_) / 2 * r:
                bad_plane += 1
            X = [S[0][0] * a + S[0][1] * b, S[1][0] * a + S[1][1] * b]
        done += 1
    return {"laws": n, "steps": steps, "E1a_violations": bad_null, "E1b_violations": bad_plane}


# ----------------------------------------------------------------------- E2, E3, E5
def e235_float(rng, n=2000):
    out = {"laws": 0, "elliptic": 0, "hyperbolic": 0, "attractive_hyp": 0, "repulsive_hyp": 0,
           "E2a_max": 0.0, "E2b_violations": 0, "E3_max": 0.0, "E5a_violations": 0, "E5b_max": 0.0, "E5d_violations": 0}
    while out["laws"] < n:
        S = random_sl2(rng)
        tau = float(np.trace(S))
        if tau < -2 or abs(abs(tau) - 2) < 0.1 or tau > 8:
            continue
        kind, ds, A = law(S)
        X0 = rng.standard_normal(2)
        E, mu, L = kepler(X0, A @ X0, kind)
        if abs(mu) < 1e-3:
            continue
        G = OM @ (S - np.linalg.inv(S)) / 2
        e_metric = math.sqrt(max(0.0, 1 - 4 * np.linalg.det(G) / np.trace(G) ** 2))
        e_kepler = math.sqrt(max(0.0, 1 + 2 * E * L * L / mu ** 2))
        out["E2a_max"] = max(out["E2a_max"], abs(e_metric - e_kepler) / max(1.0, e_kepler))
        if (kind == "elliptic") != (e_metric < 1):
            out["E2b_violations"] += 1
        sgn = -1.0 if kind == "elliptic" else 1.0
        for s in (0.37, 1.23):                       # 補間した軌道 X(s)
            c, sh = (math.cos(s), math.sin(s)) if kind == "elliptic" else (math.cosh(s), math.sinh(s))
            X = c * X0 + sh * (A @ X0)
            dX = (sgn * sh) * X0 + c * (A @ X0)
            w, dw = complex(*X), complex(*dX)
            z, dz, ddz = w * w, 2 * w * dw, 2 * dw * dw + 2 * sgn * w * w
            r, dr = abs(w) ** 2, 2 * (w.conjugate() * dw).real
            acc = (ddz * r - dz * dr) / r ** 3       # d²z/dt²,  dt = r ds
            out["E3_max"] = max(out["E3_max"], abs(acc + mu * z / r ** 3) / abs(mu * z / r ** 3))
        if kind == "elliptic":
            out["elliptic"] += 1
            out["E5a_violations"] += mu <= 0
        else:
            out["hyperbolic"] += 1
            ev, V = np.linalg.eig(S)
            ip = int(np.argmax(np.abs(ev)))
            vp, vm = V[:, ip].real, V[:, 1 - ip].real
            al, be = np.linalg.solve(np.c_[vp, vm], X0)
            out["E5b_max"] = max(out["E5b_max"], abs(mu + 8 * al * be * float(vp @ vm)) / abs(mu))
            out["attractive_hyp" if mu > 0 else "repulsive_hyp"] += 1
            vp_, vm_ = (vp if al > 0 else -vp), (vm if be > 0 else -vm)          # 配置が入っている区画の二辺
            psi = math.degrees(math.acos(float(vp_ @ vm_) / (np.linalg.norm(vp_) * np.linalg.norm(vm_))))
            out["E5d_violations"] += (mu > 0) != (psi > 90.0)
        out["laws"] += 1
    return out


def e5c_neutral(rng):
    worst = 0.0
    for H in (0.2, 0.7, 1.5):
        for S in (np.array([[math.cosh(H), math.sinh(H)], [math.sinh(H), math.cosh(H)]]),
                  np.diag([math.exp(H), math.exp(-H)])):
            _, _, A = law(S)
            for _ in range(20):
                X = rng.standard_normal(2)
                dX = A @ X
                _, mu, _ = kepler(X, dX, "hyperbolic")
                w, dw = complex(*X), complex(*dX)
                z, dz, ddz = w * w, 2 * w * dw, 2 * dw * dw + 2 * w * w
                r, dr = abs(w) ** 2, 2 * (w.conjugate() * dw).real
                worst = max(worst, abs(mu) / float(X @ X), abs((ddz * r - dz * dr) / r ** 3) * r)
    return {"E5c_max": worst}


# --------------------------------------------------------------------------- E4
def closed_orbit(rng, n, e_target=None):
    """一周 n 歩の楕円型の法則と初期配置。読み出し側の離心率を e_target にできる。"""
    th = 2 * math.pi / n
    e = rng.uniform(0.0, 0.95) if e_target is None else e_target
    al, be = math.sqrt(1 - e), math.sqrt(1 + e)           # w-楕円の半軸（a_K = 1）
    T = np.diag([al, be])
    psi = rng.uniform(0, 2 * math.pi)
    Rm = np.array([[math.cos(psi), -math.sin(psi)], [math.sin(psi), math.cos(psi)]])
    T = Rm @ T
    Rt = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    S = T @ Rt @ np.linalg.inv(T)
    u0 = rng.uniform(0, 2 * math.pi)
    X0 = T @ np.array([math.cos(u0), math.sin(u0)])
    return S, X0, th, e


def e4_clock(rng):
    res = {"E4a_max": 0.0, "E4b_max": 0.0, "E4c_violations": 0, "E4c_cases": 0,
           "E4c_monotone_although_cos_le_e": 0, "E4c_nonmonotone_although_cos_gt_e": 0}
    for _ in range(400):
        n = int(rng.integers(5, 61))
        S, X0, th, e = closed_orbit(rng, n)
        Xs = [X0]
        for _k in range(n):
            Xs.append(S @ Xs[-1])
        dts = np.array([(th / math.cos(th)) * float(Xs[k] @ Xs[k + 1]) for k in range(n)])
        zs = [complex(*X) ** 2 for X in Xs]
        tri = np.array([0.5 * (zs[k].conjugate() * zs[k + 1]).imag for k in range(n)])
        monotone = bool(np.all(dts > 0))
        res["E4c_cases"] += 1
        if monotone != (math.cos(th) > e):
            res["E4c_violations"] += 1
            res["E4c_monotone_although_cos_le_e" if monotone else "E4c_nonmonotone_although_cos_gt_e"] += 1
        ratio = tri / dts
        res["E4a_max"] = max(res["E4a_max"], float(np.max(np.abs(ratio / ratio[0] - 1))))
        _, _, A = law(S)
        _, mu, _ = kepler(X0, A @ X0, "elliptic")
        res["E4b_max"] = max(res["E4b_max"], abs(float(dts.sum()) / (2 * math.pi * mu / 4) - 1))
    # E4c': 診断。X_k·X_{k+1} = a[cosθ − e cos(u_k+θ)] なので、cosθ>e は「どの初期位相でも単調」の条件
    bad_i = bad_ii = 0
    for _ in range(200):
        n = int(rng.integers(5, 61))
        th = 2 * math.pi / n
        e = float(rng.uniform(0.0, 0.95))
        if abs(math.cos(th) - e) < 1e-3:
            continue
        found_negative = False
        for u0 in np.linspace(0.0, 2 * th, 400, endpoint=False):
            ticks = [math.cos(th) - e * math.cos(u0 + 2 * k * th + th) for k in range(n)]
            if min(ticks) < 0:
                found_negative = True
                break
        if math.cos(th) > e and found_negative:
            bad_i += 1
        if math.cos(th) < e and not found_negative:
            bad_ii += 1
    res.update({"E4c_prime_bad_i": bad_i, "E4c_prime_bad_ii": bad_ii})
    # 刻みの式そのものの確認: X_k·X_{k+1} = a[cosθ − e cos(u_k+θ)]（a=1）
    S, X0, th, e = closed_orbit(rng, 24, e_target=0.6)
    Xs = [X0]
    for _k in range(24):
        Xs.append(S @ Xs[-1])
    dots = np.array([float(Xs[k] @ Xs[k + 1]) for k in range(24)])
    res["tick_formula_check"] = {"mean_over_cos_theta": float(dots.mean() / math.cos(th)),
                                 "amplitude_over_e": float((dots.max() - dots.min()) / 2 / e) }
    # E4d: 離散の時計での差分加速度の収束次数
    ns, errs = [32, 64, 128, 256, 512, 1024], []
    for n in ns:
        S, X0, th, e = closed_orbit(np.random.default_rng(SEED + 1), n, e_target=0.4)
        _, _, A = law(S)
        _, mu, _ = kepler(X0, A @ X0, "elliptic")
        Xs = [X0]
        for _k in range(n + 1):
            Xs.append(S @ Xs[-1])
        zs = [complex(*X) ** 2 for X in Xs]
        dt = [(th / math.cos(th)) * float(Xs[k] @ Xs[k + 1]) for k in range(n + 1)]
        worst = 0.0
        for k in range(1, n):
            acc = 2 * ((zs[k + 1] - zs[k]) / dt[k] - (zs[k] - zs[k - 1]) / dt[k - 1]) / (dt[k] + dt[k - 1])
            r = abs(zs[k])
            worst = max(worst, abs(acc + mu * zs[k] / r ** 3) / (mu / r ** 2))
        errs.append(worst)
    slope = float(np.polyfit(np.log([2 * math.pi / n for n in ns]), np.log(errs), 1)[0])
    res.update({"E4d_n": ns, "E4d_errors": errs, "E4d_order": slope})
    return res


# --------------------------------------------------------------------------- E6
def e6_collision():
    bad = 0
    tau = Fraction(17, 10)
    X = [[Fraction(3, 2), Fraction(-1, 3)]]
    X.append([Fraction(7, 10) * X[0][0], Fraction(7, 10) * X[0][1]])     # X_1 ∥ X_0
    for k in range(1, 40):
        X.append([tau * X[k][0] - X[k - 1][0], tau * X[k][1] - X[k - 1][1]])
    for k in range(40):
        if X[k][0] * X[k + 1][1] - X[k][1] * X[k + 1][0] != 0:
            bad += 1
    zs = [complex(float(a), float(b)) ** 2 for a, b in X]
    ray = zs[0] / abs(zs[0])
    off = max(abs(z / ray - abs(z)) for z in zs)             # 同じ半直線上なら 0
    crossings = sum(1 for k in range(40) if X[k][0] * X[k + 1][0] + X[k][1] * X[k + 1][1] < 0)
    return {"E6_q_nonzero": bad, "E6_max_off_ray": off, "E6_finite": bool(np.isfinite(zs).all()),
            "E6_origin_crossings_in_40_steps": crossings,
            "E6_r_min_over_r_max": min(abs(z) for z in zs) / max(abs(z) for z in zs)}


def e6b_collision_in_d2():
    """D2 の範囲にある衝突軌道。tau = lam + 1/lam なので |tau| >= 2。すべて有理数の厳密計算。"""
    lam = Fraction(7, 10)
    tau = lam + 1 / lam                                   # 149/70
    u = [Fraction(3, 2), Fraction(-1, 3)]
    u2 = [u[0] * u[0] - u[1] * u[1], 2 * u[0] * u[1]]     # u^2（読み出しの向き）

    def run(c0, c1, n=40):
        c = [c0, c1]
        for k in range(1, n):
            c.append(tau * c[k] - c[k - 1])
        return c

    def checks(c):
        X = [[ck * u[0], ck * u[1]] for ck in c]
        q_bad = sum(1 for k in range(len(c) - 1)
                    if X[k][0] * X[k + 1][1] - X[k][1] * X[k + 1][0] != 0)
        off = sum(1 for ck in c
                  if (ck * ck * u2[0]) * (-u2[1]) + (ck * ck * u2[1]) * u2[0] != 0)
        cross = sum(1 for k in range(len(c) - 1) if c[k] * c[k + 1] < 0)
        return q_bad, off, cross

    pure, mixed = run(Fraction(1), lam), run(Fraction(1), Fraction(-1))
    qa, offa, xa = checks(pure)
    qb, offb, xb = checks(mixed)
    return {"tau": str(tau), "tau_ge_2": tau >= 2,
            "pure_is_geometric": all(pure[k] == lam ** k for k in range(40)),
            "pure_q_nonzero": qa, "pure_off_ray": offa, "pure_origin_crossings": xa,
            "mixed_q_nonzero": qb, "mixed_off_ray": offb, "mixed_origin_crossings": xb}


# --------------------------------------------------------------------------- E7
def e7_readable(rng):
    worst = 0.0
    bad = 0
    for _ in range(200):
        S, X0, th, e = closed_orbit(rng, int(rng.integers(7, 40)))
        ps = rng.uniform(0, 2 * math.pi)
        Rm = np.array([[math.cos(ps), -math.sin(ps)], [math.sin(ps), math.cos(ps)]])
        def summary(S_, X_):
            G = OM @ (S_ - np.linalg.inv(S_)) / 2
            ecc = math.sqrt(max(0.0, 1 - 4 * np.linalg.det(G) / np.trace(G) ** 2))
            Xs = [X_]
            for _k in range(6):
                Xs.append(S_ @ Xs[-1])
            zs = [complex(*X) ** 2 for X in Xs]
            dphi = [np.angle(zs[k + 1] / zs[k]) for k in range(6)]
            rr = [abs(zs[k + 1]) / abs(zs[k]) for k in range(6)]
            return np.array([ecc] + dphi + rr)
        worst = max(worst, float(np.max(np.abs(summary(S, X0) - summary(Rm @ S @ Rm.T, Rm @ X0)))))
        x, y, r = readout(X0)
        xp, yp, rp = readout(P @ X0)
        if max(abs(xp + x), abs(yp - y), abs(rp - r)) > 1e-12 or np.max(np.abs(readout(-X0) - readout(X0))) > 0:
            bad += 1
    return {"E7a_max": worst, "E7b_violations": bad}


# ---------------------------------------------------------------------------
def main():
    rng = np.random.default_rng(SEED)
    r1, r235, r5c = e1_exact(rng), e235_float(rng), e5c_neutral(rng)
    r4, r6, r7 = e4_clock(rng), e6_collision(), e7_readable(rng)
    r6b = e6b_collision_in_d2()
    verdict = {
        "E1a": r1["E1a_violations"] == 0, "E1b": r1["E1b_violations"] == 0,
        "E2a": r235["E2a_max"] <= 1e-9, "E2b": r235["E2b_violations"] == 0,
        "E3": r235["E3_max"] <= 1e-9,
        "E4a": r4["E4a_max"] <= 1e-12, "E4b": r4["E4b_max"] <= 1e-12,
        "E4c": r4["E4c_violations"] == 0,
        "E4c_prime": r4["E4c_prime_bad_i"] == 0 and r4["E4c_prime_bad_ii"] == 0,
        "E4d": 1.8 <= r4["E4d_order"] <= 2.2,
        "E5a": r235["E5a_violations"] == 0,
        "E5b": r235["E5b_max"] <= 1e-9 and r235["attractive_hyp"] > 0 and r235["repulsive_hyp"] > 0,
        "E5c": r5c["E5c_max"] <= 1e-12, "E5d": r235["E5d_violations"] == 0,
        "E6": r6["E6_q_nonzero"] == 0 and r6["E6_max_off_ray"] <= 1e-12 and r6["E6_finite"],
        "E6b": (r6b["tau_ge_2"] and r6b["pure_is_geometric"]
                and r6b["pure_q_nonzero"] == 0 and r6b["mixed_q_nonzero"] == 0
                and r6b["pure_off_ray"] == 0 and r6b["mixed_off_ray"] == 0
                and r6b["pure_origin_crossings"] == 0 and r6b["mixed_origin_crossings"] == 1),
        "E7a": r7["E7a_max"] <= 1e-10, "E7b": r7["E7b_violations"] == 0,
    }
    results = {"seed": SEED, "fail_conditions": FAIL_CONDITIONS, "E1": r1, "E2_E3_E5": r235, "E5c": r5c,
               "E4": r4, "E6": r6, "E6b": r6b, "E7": r7, "verdict": {k: "PASS" if v else "FAIL" for k, v in verdict.items()}}
    with open(RESULTS_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=1, default=float)
    for k in verdict:
        print(f"{k}: {'PASS' if verdict[k] else 'FAIL'}  -- 落ちる条件: {FAIL_CONDITIONS[k]}")
    print(json.dumps({k: results[k] for k in ("E1", "E2_E3_E5", "E5c", "E4", "E6", "E7")}, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
