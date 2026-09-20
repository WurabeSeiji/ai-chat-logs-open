#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_area_readout_experiments_v1.py

第一思考実験「等価原理を無名な二自由度まで遡る」の数値検証 E1〜E6 と図の生成。

設計上の約束（設計書 §12 への適合）:
  * 状態更新の中で丸め・クリッピング・正規化を一切行わない。
  * 乱数の種は固定し、結果は JSON に保存する。
  * 各検証は「落ちる条件」をコード冒頭の FAIL_CONDITIONS に事前に宣言し、
    実行結果と照合して PASS / FAIL を出す。
  * パスはこのスクリプトの置き場所からの相対パス。特定環境の絶対パスを使わない。
  * 出力ファイル名は論文本文に書いた名前と一致させる。

出力（すべてこのスクリプトと同じフォルダ）:
  * area_readout_experiments_results_v1.json
  * fig01〜fig03 の SVG（依存なしの簡易描画。常に生成される）
  * fig01〜fig03 の PNG（matplotlib がある場合だけ生成される）

依存: numpy, mpmath。PNG も作る場合は matplotlib。
実行: python3 run_area_readout_experiments_v1.py
"""
from __future__ import annotations

import json
import math
import os
from fractions import Fraction
from pathlib import Path

import numpy as np
import mpmath as mp

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    HAVE_MPL = True
except ImportError:          # PNG は作らず、SVG と JSON だけを出力する
    HAVE_MPL = False

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE
RESULTS_JSON = HERE / "area_readout_experiments_results_v1.json"

FIG1 = "fig01_area_readout_three_classes_ja_v1"
FIG2 = "fig02_area_readout_closure_staircase_ja_v1"
FIG3 = "fig03_area_readout_metric_from_orbit_ja_v1"

SEED = 20260920

# ---------------------------------------------------------------------------
# 事前宣言：各検証が「落ちる」条件
# ---------------------------------------------------------------------------
FAIL_CONDITIONS = {
    "E1a": "有理数厳密計算で、面積読み出し q_k が一度でも q_0 と異なる。",
    "E1b": "浮動小数で、尺度で正規化した q のずれ max|q_k-q_0|/(|X_k||X_{k+1}|) が 1e-9 を超える。",
    "E1c": "誘導計量 G の符号型（固有値の符号）が τ による分類と一件でも食い違う。",
    "E1d": "（対照）楕円型の一般の記述で、素朴な a^2+b^2 が保存してしまう（中央値のゆらぎ < 1%）。",
    "E2a": "二階の無向則 X_{k+1}+X_{k-1}=τX_k の相対残差が 1e-9 を超える。",
    "E2b": "回転 (n>=3) で |PSP - S^-1| が 1e-12 を超える、または |PSP - S| が 1e-3 未満になる。",
    "E2c": "置換同変な実線形写像 γI+β11^T (N=2..8) に虚部 > 1e-12 の固有値が出る。",
    "E3a": "有限分解能での初回帰還数 n が、α の連分数近似分母でない値を一度でも取る。",
    "E3b": "初回帰還数が上界 n < 1/δ を破る。",
    "E4a": "有理数厳密計算で |ΔX|_G^2 = κ X^T G X または Δ^2 X = -κ X が破れる。",
    "E4b": "浮動小数で（離散曲率）×|C| が 1 から 1e-7 を超えてずれる。",
    "E4c": "同一の S の二軌道で、誘導固有長の比が |C1/C2| から 1e-7 を超えてずれる。",
    "E4b_prime": "【事後に追加した診断検証】初期点の前後一歩 (X_{-1}, X_0, X_1) だけで評価したとき、（離散曲率）×|C| が 1 から 1e-9 を超えてずれる。",
    "E6": "【E4b の失敗を受けて立てた予測】双曲型で r_k = |q|/(|X_k||X_{k+1}|) の対数の傾きが -2H (H=arccosh(|τ|/2)) から 1% を超えてずれる。",
    "E5a": "n>=3 の閉軌道で、二次モーメントからの計量読み出しの相対誤差が 1e-10 を超える。",
    "E5b": "n>=3 で、n が p を割らないのに冪和 |Σ z_k^p| / n が 1e-10 を超える。",
    "E5c": "n=2 で面積の自己読み出しが恒等的にゼロにならない（|G| > 1e-12）。",
}

OMEGA = np.array([[0.0, 1.0], [-1.0, 0.0]])
P_SWAP = np.array([[0.0, 1.0], [1.0, 0.0]])


# ---------------------------------------------------------------------------
# 基本演算
# ---------------------------------------------------------------------------
def omega(x: np.ndarray, y: np.ndarray) -> float:
    """面積形式 ω(X,Y) = a_X b_Y - b_X a_Y。"""
    return float(x[0] * y[1] - x[1] * y[0])


def random_sl2(rng: np.random.Generator) -> np.ndarray:
    """SL(2,R) の乱数元。ガウス行列を det で規格化（det<0 なら行を入れ替える）。"""
    while True:
        a = rng.standard_normal((2, 2))
        d = np.linalg.det(a)
        if abs(d) < 1e-3:
            continue
        if d < 0:
            a = a[::-1, :].copy()
            d = -d
        return a / math.sqrt(d)


def induced_metric(s: np.ndarray) -> np.ndarray:
    """G = Ω (S - S^-1)/2。X^T G X = ω(X, SX)。"""
    return OMEGA @ (s - np.linalg.inv(s)) / 2.0


def classify(tau: float) -> str:
    if abs(tau) < 2.0:
        return "elliptic"
    if abs(tau) > 2.0:
        return "hyperbolic"
    return "parabolic"


# ---------------------------------------------------------------------------
# E1: 面積読み出しの保存・符号型・素朴な二乗和との対照
# ---------------------------------------------------------------------------
def frac_rand(rng: np.random.Generator, lo=-6, hi=6, dmax=5) -> Fraction:
    n = int(rng.integers(lo, hi + 1))
    d = int(rng.integers(1, dmax + 1))
    return Fraction(n, d)


def exact_tests(rng: np.random.Generator, n_samples=300, steps=12):
    """有理数の厳密計算。E1a と E4a。"""
    e1a_fail = 0
    e4a_fail = 0
    law_fail = 0
    done = 0
    while done < n_samples:
        p = frac_rand(rng)
        q = frac_rand(rng)
        r = frac_rand(rng)
        if p == 0:
            continue
        s22 = (1 + q * r) / p          # det = 1 を厳密に満たす
        S = [[p, q], [r, s22]]
        Sinv = [[s22, -q], [-r, p]]
        tau = p + s22
        kap = 2 - tau
        x = [frac_rand(rng), frac_rand(rng)]
        if x[0] == 0 and x[1] == 0:
            continue

        def mul(m, v):
            return [m[0][0] * v[0] + m[0][1] * v[1], m[1][0] * v[0] + m[1][1] * v[1]]

        def om(u, v):
            return u[0] * v[1] - u[1] * v[0]

        # G = Ω (S - S^-1)/2
        n_ = [[(S[i][j] - Sinv[i][j]) / 2 for j in range(2)] for i in range(2)]
        G = [[n_[1][0], n_[1][1]], [-n_[0][0], -n_[0][1]]]

        def qf(v):
            gv = mul(G, v)
            return v[0] * gv[0] + v[1] * gv[1]

        q0 = om(x, mul(S, x))
        xs = [x]
        for _ in range(steps + 1):
            xs.append(mul(S, xs[-1]))
        for k in range(steps):
            qk = om(xs[k], xs[k + 1])
            if qk != q0:
                e1a_fail += 1
            if qf(xs[k]) != qk:
                e1a_fail += 1
        for k in range(1, steps):
            # 無向則
            for c in range(2):
                if xs[k + 1][c] + xs[k - 1][c] != tau * xs[k][c]:
                    law_fail += 1
            dx = [xs[k + 1][c] - xs[k][c] for c in range(2)]
            d2 = [xs[k + 1][c] - 2 * xs[k][c] + xs[k - 1][c] for c in range(2)]
            if qf(dx) != kap * qf(xs[k]):
                e4a_fail += 1
            for c in range(2):
                if d2[c] != -kap * xs[k][c]:
                    e4a_fail += 1
        done += 1
    return {"samples": n_samples, "steps": steps,
            "E1a_violations": e1a_fail, "E4a_violations": e4a_fail,
            "second_order_law_violations_exact": law_fail}


def float_tests(rng: np.random.Generator, n_samples=3000, steps=40):
    drift_max = 0.0
    law_res_max = 0.0
    sig_mismatch = 0
    counts = {"elliptic": 0, "hyperbolic": 0, "parabolic": 0}
    naive_fluct = []
    curv_err_max = 0.0
    curv_err_by_class = {"elliptic": 0.0, "hyperbolic": 0.0}
    curv_err_prime_max = 0.0
    amp_at_worst = 0.0
    worst_tau = 0.0
    ratio_err_max = 0.0
    n_curv = 0
    for _ in range(n_samples):
        S = random_sl2(rng)
        tau = float(np.trace(S))
        cls = classify(tau)
        counts[cls] += 1
        G = induced_metric(S)
        ev = np.linalg.eigvalsh((G + G.T) / 2)
        if cls == "elliptic":
            ok = (ev[0] > 0 and ev[1] > 0) or (ev[0] < 0 and ev[1] < 0)
        elif cls == "hyperbolic":
            ok = ev[0] < 0 < ev[1]
        else:
            ok = True
        if not ok:
            sig_mismatch += 1

        x0 = rng.standard_normal(2)
        xs = [x0]
        for _k in range(steps + 1):
            xs.append(S @ xs[-1])          # 状態更新に丸め・クリップなし
        q0 = omega(xs[0], xs[1])
        for k in range(steps):
            qk = omega(xs[k], xs[k + 1])
            scale = np.linalg.norm(xs[k]) * np.linalg.norm(xs[k + 1])
            drift_max = max(drift_max, abs(qk - q0) / scale)
        for k in range(1, steps):
            res = xs[k + 1] + xs[k - 1] - tau * xs[k]
            sc = np.linalg.norm(xs[k + 1]) + np.linalg.norm(xs[k - 1]) + abs(tau) * np.linalg.norm(xs[k])
            law_res_max = max(law_res_max, float(np.linalg.norm(res) / sc))
        if cls == "elliptic":
            naive = np.array([v @ v for v in xs[:steps]])
            naive_fluct.append(float((naive.max() - naive.min()) / naive.mean()))

        # E4b/E4c: 曲率恒等式（放物型近傍と零錐近傍は桁落ちするので除く）
        detG = float(np.linalg.det(G))
        if abs(detG) > 1e-2:
            Gn = G / math.sqrt(abs(detG))
            scale0 = np.linalg.norm(xs[0]) * np.linalg.norm(xs[1])
            if abs(q0) > 1e-2 * scale0:
                C = math.sqrt(abs(xs[0] @ Gn @ xs[0]))
                for k in range(1, 4):
                    dx = xs[k + 1] - xs[k]
                    d2 = xs[k + 1] - 2 * xs[k] + xs[k - 1]
                    ds2 = abs(dx @ Gn @ dx)
                    a2 = abs(d2 @ Gn @ d2)
                    curv = math.sqrt(a2) / ds2
                    err = abs(curv * C - 1.0)
                    curv_err_by_class[cls] = max(curv_err_by_class[cls], err)
                    if err > curv_err_max:
                        curv_err_max = err
                        # 桁落ちの増幅率の目安: |X_k|^2 |Gn| / C^2
                        amp_at_worst = float((xs[k] @ xs[k]) * np.linalg.norm(Gn) / (C * C))
                        worst_tau = tau
                # E4b': 成長を挟まない評価（X_{-1}, X_0, X_1）
                xm1 = np.linalg.inv(S) @ xs[0]
                dxp = xs[1] - xs[0]
                d2p = xs[1] - 2 * xs[0] + xm1
                curvp = math.sqrt(abs(d2p @ Gn @ d2p)) / abs(dxp @ Gn @ dxp)
                curv_err_prime_max = max(curv_err_prime_max, abs(curvp * C - 1.0))
                # 同じ S の別軌道との固有長の比
                y0 = rng.standard_normal(2)
                qy = omega(y0, S @ y0)
                if abs(qy) > 1e-2 * np.linalg.norm(y0) * np.linalg.norm(S @ y0):
                    C2 = math.sqrt(abs(y0 @ Gn @ y0))
                    dx1 = xs[1] - xs[0]
                    dy1 = S @ y0 - y0
                    l1 = math.sqrt(abs(dx1 @ Gn @ dx1))
                    l2 = math.sqrt(abs(dy1 @ Gn @ dy1))
                    ratio_err_max = max(ratio_err_max, abs((l1 / l2) / (C / C2) - 1.0))
                n_curv += 1
    naive_fluct = np.array(naive_fluct)
    return {
        "samples": n_samples, "steps": steps, "class_counts": counts,
        "E1b_max_scaled_drift": drift_max,
        "E1c_signature_mismatches": sig_mismatch,
        "E1d_naive_sumsq_fluct_median": float(np.median(naive_fluct)),
        "E1d_naive_sumsq_fluct_frac_gt_1pct": float(np.mean(naive_fluct > 0.01)),
        "E2a_max_rel_residual": law_res_max,
        "E4_samples_used": n_curv,
        "E4b_max_curvature_identity_error": curv_err_max,
        "E4b_max_error_by_class": curv_err_by_class,
        "E4b_amplification_at_worst_case": amp_at_worst,
        "E4b_tau_at_worst_case": worst_tau,
        "E4b_prime_max_error": curv_err_prime_max,
        "E4c_max_proper_length_ratio_error": ratio_err_max,
    }


# ---------------------------------------------------------------------------
# E2: 名前の付け替えと向き、置換同変な線形写像のスペクトル
# ---------------------------------------------------------------------------
def relabel_tests(rng: np.random.Generator):
    rot = {}
    for n in range(3, 13):
        th = 2 * math.pi / n
        S = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        psp = P_SWAP @ S @ P_SWAP
        rot[n] = {"|PSP-S^-1|": float(np.linalg.norm(psp - np.linalg.inv(S))),
                  "|PSP-S|": float(np.linalg.norm(psp - S))}
    max_im = 0.0
    for N in range(2, 9):
        for _ in range(300):
            g, b = rng.standard_normal(2) * 3
            A = g * np.eye(N) + b * np.ones((N, N))
            ev = np.linalg.eigvals(A)
            max_im = max(max_im, float(np.max(np.abs(ev.imag))))
    # 厳密に無名（P と可換）でユニモジュラ → |τ| >= 2
    strict_elliptic = 0
    for _ in range(2000):
        H = rng.uniform(-3, 3)
        sgn = rng.choice([-1.0, 1.0])
        S = sgn * (math.cosh(H) * np.eye(2) + math.sinh(H) * P_SWAP)
        if abs(np.trace(S)) < 2.0 - 1e-12:
            strict_elliptic += 1
    return {"rotation_relabel": rot,
            "E2c_max_imag_eig_perm_equivariant": max_im,
            "strict_nameless_elliptic_count": strict_elliptic}


# ---------------------------------------------------------------------------
# E3: 有限分解能のもとでの閉路数（初回帰還）と連分数
# ---------------------------------------------------------------------------
def convergent_denominators(alpha_mp, n_terms=40):
    mp.mp.dps = 80
    x = mp.mpf(alpha_mp)
    a = []
    for _ in range(n_terms):
        ai = int(mp.floor(x))
        a.append(ai)
        fr = x - ai
        if fr < mp.mpf(10) ** (-60):
            break
        x = 1 / fr
    qs = []
    q_m2, q_m1 = 1, 0                      # q_{-2}=1, q_{-1}=0
    for ai in a:
        qn = ai * q_m1 + q_m2
        qs.append(qn)
        q_m2, q_m1 = q_m1, qn
    return a, qs


def closure_tests():
    mp.mp.dps = 80
    alphas = {
        "golden (sqrt5-1)/2": (mp.sqrt(5) - 1) / 2,
        "sqrt2-1": mp.sqrt(2) - 1,
        "e-2": mp.e - 2,
        "pi-3": mp.pi - 3,
    }
    rhos = np.logspace(0, 5, 241)           # ρ = C/ε
    n_max = int(2 * math.pi * rhos[-1]) + 10
    out = {}
    e3a_viol = 0
    e3b_viol = 0
    curves = {}
    for name, al in alphas.items():
        a_cf, qs = convergent_denominators(al)
        qset = set(qs)
        al_f = float(al)
        n = np.arange(1, n_max + 1, dtype=np.float64)
        fracp = np.mod(n * al_f, 1.0)
        dist = np.minimum(fracp, 1.0 - fracp)
        runmin = np.minimum.accumulate(dist)
        ns = []
        for rho in rhos:
            delta = math.asin(min(1.0, 1.0 / (2.0 * rho))) / math.pi
            idx = int(np.argmax(runmin < delta))
            if not (runmin[idx] < delta):
                ns.append(None)
                continue
            n_first = idx + 1
            # 高精度で再確認
            d_hp = abs(n_first * al - mp.nint(n_first * al))
            assert d_hp < delta * (1 + 1e-6)
            if n_first not in qset:
                e3a_viol += 1
            if not (n_first < 1.0 / delta):
                e3b_viol += 1
            ns.append(n_first)
        distinct = sorted(set(v for v in ns if v is not None))
        ratio = [v / (2 * math.pi * r) for v, r in zip(ns, rhos) if v is not None and r >= 10]
        out[name] = {
            "partial_quotients_head": a_cf[:12],
            "convergent_denominators_head": qs[:16],
            "distinct_closure_numbers": distinct,
            "n_over_2pi_rho_min(rho>=10)": float(min(ratio)),
            "n_over_2pi_rho_max(rho>=10)": float(max(ratio)),
        }
        curves[name] = ns
    return {"rhos_spec": "numpy.logspace(0, 5, 241)", "per_alpha": out,
            "E3a_violations": e3a_viol, "E3b_violations": e3b_viol}, curves, rhos


# ---------------------------------------------------------------------------
# E5: 閉軌道の二次モーメントから計量を読む／冪和のゼロ閉塞
# ---------------------------------------------------------------------------
def metric_readout_tests(rng: np.random.Generator, trials=200):
    per_n = {}
    worst_pow = 0.0
    for n in range(3, 13):
        errs = []
        for _ in range(trials):
            ms = [m for m in range(1, n) if math.gcd(m, n) == 1]
            m = int(rng.choice(ms))
            th = 2 * math.pi * m / n
            R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
            T = random_sl2(rng)
            S = T @ R @ np.linalg.inv(T)
            G = induced_metric(S)
            detG = float(np.linalg.det(G))
            Gn = G / math.sqrt(abs(detG))
            if np.trace(Gn) < 0:
                Gn = -Gn
            x = rng.standard_normal(2)
            C2 = float(x @ Gn @ x)
            M2 = np.zeros((2, 2))
            xk = x.copy()
            for _k in range(n):
                M2 += np.outer(xk, xk)
                xk = S @ xk
            pred = (n * C2 / 2.0) * np.linalg.inv(Gn)
            errs.append(float(np.linalg.norm(M2 - pred) / np.linalg.norm(pred)))
        per_n[n] = {"max_rel_err": float(max(errs)), "median_rel_err": float(np.median(errs))}
        # 正規記述での冪和
        z = np.exp(2j * np.pi * np.arange(n) / n) * np.exp(1j * rng.uniform(0, 2 * np.pi))
        for p in range(1, 2 * n + 1):
            s = abs(np.sum(z ** p)) / n
            if p % n != 0:
                worst_pow = max(worst_pow, float(s))
    # n=2: S=-I → 面積の自己読み出しは恒等的にゼロ
    S2 = -np.eye(2)
    G2 = induced_metric(S2)
    x = rng.standard_normal(2)
    M2 = np.outer(x, x) + np.outer(-x, -x)
    return {"per_n": per_n,
            "E5b_worst_nonmultiple_power_sum": worst_pow,
            "n2_|G|": float(np.linalg.norm(G2)),
            "n2_rank_M2": int(np.linalg.matrix_rank(M2)),
            "n2_sum_z2_over_n": 1.0}


# ---------------------------------------------------------------------------
# E6: 双曲型での射影的収束（保存量の相対的な読めなさ）
# ---------------------------------------------------------------------------
def projective_attractor_tests(rng: np.random.Generator, n_samples=300):
    worst = 0.0
    used = 0
    while used < n_samples:
        S = random_sl2(rng)
        tau = float(np.trace(S))
        if not (2.2 < abs(tau) < 20.0):
            continue
        H = math.acosh(abs(tau) / 2.0)
        x = rng.standard_normal(2)
        q0 = omega(x, S @ x)
        if abs(q0) < 1e-2 * np.linalg.norm(x) * np.linalg.norm(S @ x):
            continue
        ks, lr = [], []
        xk = x.copy()
        for k in range(0, 200):
            xn = S @ xk
            r = abs(q0) / (np.linalg.norm(xk) * np.linalg.norm(xn))
            if k >= 6 and r > 1e-9:
                ks.append(k)
                lr.append(math.log(r))
            if r < 1e-9:
                break
            xk = xn
        if len(ks) < 4:
            continue
        slope = float(np.polyfit(ks, lr, 1)[0])
        worst = max(worst, abs(slope / (-2.0 * H) - 1.0))
        used += 1
    return {"samples": n_samples, "E6_max_rel_slope_error": worst}


# ---------------------------------------------------------------------------
# 図のデータ（PNG と SVG で共通）
# ---------------------------------------------------------------------------
def figure1_data():
    T = np.array([[1.6, 0.7], [0.2, 0.9]])
    T = T / math.sqrt(np.linalg.det(T))
    Ti = np.linalg.inv(T)
    th = 2 * math.pi * (math.sqrt(5) - 1) / 2 / 3.0
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    Se = T @ R @ Ti
    x = np.array([1.0, 0.2])
    pts = [x]
    for _ in range(60):
        pts.append(Se @ pts[-1])
    pts = np.array(pts)
    r0 = float(np.linalg.norm(Ti @ x))
    ph = np.linspace(0, 2 * math.pi, 73)
    ell = (T @ np.vstack([r0 * np.cos(ph), r0 * np.sin(ph)])).T
    H = 0.35
    hyp = {}
    for C in (0.5, 0.8, 1.1):
        hyp[C] = np.array([[C * math.cosh(k * H), C * math.sinh(k * H)] for k in range(-3, 4)])
    Sp = np.array([[1.0, 0.0], [0.4, 1.0]])
    par = {}
    for a0 in (0.5, 1.0, 1.5):
        q = [np.array([a0, 0.0])]
        for _ in range(8):
            q.append(Sp @ q[-1])
        par[a0] = np.array(q)
    naive = np.array([p @ p for p in pts[:60]])
    area = np.array([omega(pts[i], pts[i + 1]) for i in range(60)])
    return {"pts": pts, "ell": ell, "hyp": hyp, "par": par,
            "naive": naive / naive[0], "area": area / area[0]}


# ---------------------------------------------------------------------------
# 依存なしの簡易 SVG 描画（文字はテキストのまま埋め込む）
# ---------------------------------------------------------------------------
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def _n(v):
    s = f"{v:.1f}"
    return s[:-2] if s.endswith(".0") else s


def _esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Svg:
    def __init__(self, w, h):
        self.w, self.h, self.items = w, h, []

    def add(self, s):
        self.items.append(s)

    def text(self, x, y, t, size=12, anchor="start", rot=None):
        tr = f' transform="rotate({rot} {_n(x)} {_n(y)})"' if rot is not None else ""
        body = _esc(t)
        if "^" in t:                       # "10^-5" → 10 の右肩に指数（dy で持ち上げる）
            base, ex = t.split("^", 1)
            body = f'{_esc(base)}<tspan dy="-5" font-size="{size - 3}">{_esc(ex)}</tspan>'
        self.add(f'<text x="{_n(x)}" y="{_n(y)}" font-size="{size}" text-anchor="{anchor}"{tr}>{body}</text>')

    def save(self, path):
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" '
                f'font-family="Hiragino Sans, Yu Gothic, Noto Sans CJK JP, sans-serif">\n'
                f'<rect width="{self.w}" height="{self.h}" fill="#fff"/>\n')
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(head + "\n".join(self.items) + "\n</svg>\n")


class Axes:
    def __init__(self, svg, x0, y0, w, h, xlim, ylim, xlog=False, ylog=False):
        self.s, self.x0, self.y0, self.w, self.h = svg, x0, y0, w, h
        self.xlog, self.ylog = xlog, ylog
        self.xl = [math.log10(v) for v in xlim] if xlog else list(xlim)
        self.yl = [math.log10(v) for v in ylim] if ylog else list(ylim)

    def X(self, v):
        v = math.log10(v) if self.xlog else v
        return self.x0 + (v - self.xl[0]) / (self.xl[1] - self.xl[0]) * self.w

    def Y(self, v):
        v = math.log10(v) if self.ylog else v
        return self.y0 + self.h - (v - self.yl[0]) / (self.yl[1] - self.yl[0]) * self.h

    def frame(self, xticks, yticks, xfmt=None, yfmt=None, title="", xlabel="", ylabel="", ylab_dx=40):
        s = self.s
        g = []
        for v in xticks:
            g.append(f"M{_n(self.X(v))} {_n(self.y0)}v{_n(self.h)}")
        for v in yticks:
            g.append(f"M{_n(self.x0)} {_n(self.Y(v))}h{_n(self.w)}")
        s.add(f'<path d="{"".join(g)}" stroke="#ddd" fill="none"/>')
        s.add(f'<rect x="{_n(self.x0)}" y="{_n(self.y0)}" width="{_n(self.w)}" height="{_n(self.h)}" fill="none" stroke="#000"/>')
        for v in xticks:
            s.text(self.X(v), self.y0 + self.h + 15, (xfmt or _n)(v), 11, "middle")
        for v in yticks:
            s.text(self.x0 - 6, self.Y(v) + 4, (yfmt or _n)(v), 11, "end")
        if title:
            s.text(self.x0 + self.w / 2, self.y0 - 8, title, 13, "middle")
        if xlabel:
            s.text(self.x0 + self.w / 2, self.y0 + self.h + 32, xlabel, 12, "middle")
        if ylabel:
            s.text(self.x0 - ylab_dx, self.y0 + self.h / 2, ylabel, 12, "middle", rot=-90)

    def line(self, xs, ys, color, width=1.5, dash=None):
        d = " ".join(f"{_n(self.X(x))},{_n(self.Y(y))}" for x, y in zip(xs, ys))
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.s.add(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="{width}"{da}/>')

    def dots(self, xs, ys, color, r=2.5):
        d = "".join(f"M{_n(self.X(x))} {_n(self.Y(y))}h0" for x, y in zip(xs, ys))
        self.s.add(f'<path d="{d}" stroke="{color}" stroke-width="{2 * r}" stroke-linecap="round" fill="none"/>')

    def legend(self, entries, x, y, box_w=None):
        if box_w:
            self.s.add(f'<rect x="{_n(x - 6)}" y="{_n(y - 10)}" width="{_n(box_w)}" height="{16 * len(entries) + 6}" fill="#fff" fill-opacity="0.85" stroke="#ccc"/>')
        for i, (label, color, dash) in enumerate(entries):
            yy = y + 16 * i
            da = f' stroke-dasharray="{dash}"' if dash else ""
            self.s.add(f'<path d="M{_n(x)} {_n(yy)}h22" stroke="{color}" stroke-width="2"{da}/>')
            self.s.text(x + 28, yy + 4, label, 11)


def _pow10(v):
    e = int(round(math.log10(v)))
    if e == 0:
        return "1"
    if e == 1:
        return "10"
    return "10^" + str(e).replace("-", "−")


def svg_fig1(d):
    s = Svg(900, 800)
    # (a) 楕円型。縦横の尺度をそろえる
    a = Axes(s, 60, 40, 360, 300, (-1.2, 1.2), (-1.0, 1.0))
    a.frame([-1, -0.5, 0, 0.5, 1], [-1, -0.5, 0, 0.5, 1], title="(a) 楕円型 0<κ<4：一般の記述では楕円", xlabel="a", ylabel="b")
    a.line(d["ell"][:, 0], d["ell"][:, 1], COLORS[1], 1)
    a.dots(d["pts"][:, 0], d["pts"][:, 1], COLORS[0])
    # (b) 双曲型（ブースト）。縦横の尺度をそろえる
    b = Axes(s, 500, 40, 360, 300, (-0.6, 3.0), (-1.5, 1.5))
    b.frame([0, 1, 2, 3], [-1, 0, 1], title="(b) 双曲型 κ<0：a²−b²=C²（破線は零錐 a=±b）", xlabel="a", ylabel="b")
    b.line([0, 1.5], [0, 1.5], "#000", 1, "4 3")
    b.line([0, 1.5], [0, -1.5], "#000", 1, "4 3")
    ent = []
    for i, (C, arr) in enumerate(d["hyp"].items()):
        b.line(arr[:, 0], arr[:, 1], COLORS[i], 1)
        b.dots(arr[:, 0], arr[:, 1], COLORS[i])
        ent.append((f"C={C:g}", COLORS[i], None))
    b.legend(ent, 512, 56, box_w=80)
    # (c) 放物型
    c = Axes(s, 60, 440, 360, 300, (0, 2), (-0.2, 5.2))
    c.frame([0, 0.5, 1, 1.5, 2], [0, 1, 2, 3, 4, 5], title="(c) 放物型 κ=0：a が一定、b が k に比例して増える", xlabel="a", ylabel="b")
    for i, (a0, arr) in enumerate(d["par"].items()):
        c.line(arr[:, 0], arr[:, 1], COLORS[i], 1)
        c.dots(arr[:, 0], arr[:, 1], COLORS[i])
    # (d) 二つの読み出し
    e = Axes(s, 500, 440, 360, 300, (0, 60), (0, 1.5))
    e.frame([0, 10, 20, 30, 40, 50, 60], [0, 0.5, 1, 1.5], title="(d) (a) の軌道に沿った二つの読み出し（初期値で規格化）", xlabel="k")
    k = list(range(60))
    e.line(k, d["naive"], COLORS[0], 1.2)
    e.line(k, d["area"], COLORS[1], 2.2)
    e.legend([("素朴な a²+b²（保存しない）", COLORS[0], None), ("面積読み出し ω(X_k, X_{k+1})（保存）", COLORS[1], None)], 512, 456, box_w=250)
    s.save(FIG_DIR / (FIG1 + ".svg"))


def svg_fig2(rhos, curves):
    s = Svg(900, 600)
    a = Axes(s, 70, 40, 790, 490, (1, 1e5), (0.7, 1e6), xlog=True, ylog=True)
    a.frame([10 ** e for e in range(0, 6)], [10 ** e for e in range(0, 7)], _pow10, _pow10,
            title="有限分解能のもとでの閉路数：n は θ/2π の連分数近似分母だけを取る",
            xlabel="ρ = C/ε（曲率半径 ÷ 分解能）", ylabel="初回帰還数 n（＝有効な閉路数）")
    names = {"golden (sqrt5-1)/2": "α=(√5−1)/2", "sqrt2-1": "α=√2−1", "e-2": "α=e−2", "pi-3": "α=π−3"}
    ent = []
    for i, (name, ns) in enumerate(curves.items()):
        xs, ys = [rhos[0]], [ns[0]]
        for j in range(1, len(rhos)):
            if ns[j] != ns[j - 1]:            # 値が変わる点だけを残す（階段）
                xs += [rhos[j], rhos[j]]
                ys += [ns[j - 1], ns[j]]
        xs.append(rhos[-1]); ys.append(ns[-1])
        a.line(xs, ys, COLORS[i], 1.6)
        ent.append((names.get(name, name), COLORS[i], None))
    rr = np.logspace(0, 5, 26)
    a.line(rr, [math.pi / math.asin(min(1.0, 1.0 / (2.0 * r))) for r in rr], "#000", 1, "5 4")
    ent.append(("上界 1/δ ≈ 2πC/ε", "#000", "5 4"))
    a.legend(ent, 86, 60, box_w=170)
    s.save(FIG_DIR / (FIG2 + ".svg"))


def svg_fig3(res5):
    s = Svg(800, 520)
    a = Axes(s, 90, 40, 670, 410, (1.5, 12.5), (1e-17, 1e-8), ylog=True)
    a.frame(list(range(2, 13)), [10.0 ** e for e in range(-17, -7, 2)], None, _pow10,
            title="閉軌道の二次モーメントから計量を読む", xlabel="閉路数 n", ylabel="計量読み出しの相対誤差", ylab_dx=62)
    ns = sorted(int(k) for k in res5["per_n"].keys())
    mx = [res5["per_n"][n]["max_rel_err"] for n in ns]
    a.line([1.5, 12.5], [1e-10, 1e-10], "#000", 1, "5 4")
    a.line(ns, mx, COLORS[0], 1.6)
    a.dots(ns, mx, COLORS[0], 3.5)
    a.legend([("n≥3：相対誤差の最大値（200 試行）", COLORS[0], None), ("事前に宣言した許容値 1e-10", "#000", "5 4")], 500, 60, box_w=250)
    s.text(a.X(1.7), a.Y(3e-16), "n=2：S=−I で面積の自己読み出しが", 11)
    s.text(a.X(1.7), a.Y(3e-16) + 15, "恒等的にゼロ（計量が定義できない）", 11)
    s.save(FIG_DIR / (FIG3 + ".svg"))


# ---------------------------------------------------------------------------
# PNG（matplotlib がある場合だけ）
# ---------------------------------------------------------------------------
def pick_japanese_font():
    candidates = [
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
        "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
        "/System/Library/Fonts/Hiragino Sans GB.ttc",
        "/Library/Fonts/Arial Unicode.ttf",
        "C:/Windows/Fonts/YuGothR.ttc",
        "C:/Windows/Fonts/meiryo.ttc",
        "C:/Windows/Fonts/msgothic.ttc",
    ]
    for c in candidates:
        if os.path.exists(c):
            return fm.FontProperties(fname=c)
    return None


FP = pick_japanese_font() if HAVE_MPL else None
JA = FP is not None


def L(ja: str, en: str) -> str:
    return ja if JA else en


def tkw(size=11):
    return {"fontproperties": FP, "fontsize": size} if JA else {"fontsize": size}


def make_fig1(d):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.5), dpi=170)
    ax = axes[0, 0]
    ax.plot(d["ell"][:, 0], d["ell"][:, 1], "-", lw=0.9, alpha=0.7, color="C1")
    ax.plot(d["pts"][:, 0], d["pts"][:, 1], "o", ms=3.5, color="C0")
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax.set_title(L("(a) 楕円型 0<κ<4：一般の記述では楕円", "(a) elliptic 0<kappa<4: an ellipse in a generic description"), **tkw(11))
    ax.set_xlabel("a"); ax.set_ylabel("b")

    ax = axes[0, 1]
    for C, arr in d["hyp"].items():
        ax.plot(arr[:, 0], arr[:, 1], "o-", ms=3.5, lw=0.8, label=f"C={C:g}")
    ax.plot([0, 1.5], [0, 1.5], "k--", lw=0.7)
    ax.plot([0, 1.5], [0, -1.5], "k--", lw=0.7)
    ax.set_xlim(-0.6, 3.0); ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title(L("(b) 双曲型 κ<0：a²−b²=C²（破線は零錐 a=±b）", "(b) hyperbolic kappa<0: a^2-b^2=C^2 (dashed: null cone)"), **tkw(11))
    ax.set_xlabel("a"); ax.set_ylabel("b")

    ax = axes[1, 0]
    for a0, arr in d["par"].items():
        ax.plot(arr[:, 0], arr[:, 1], "o-", ms=3.5, lw=0.8)
    ax.set_xlim(0, 2); ax.grid(alpha=0.3)
    ax.set_title(L("(c) 放物型 κ=0：a が一定、b が k に比例して増える", "(c) parabolic kappa=0: a constant, b grows linearly in k"), **tkw(11))
    ax.set_xlabel("a"); ax.set_ylabel("b")

    ax = axes[1, 1]
    k = np.arange(60)
    ax.plot(k, d["naive"], "-", lw=1.2, label=L("素朴な a²+b²（保存しない）", "naive a^2+b^2 (not conserved)"))
    ax.plot(k, d["area"], "-", lw=2.0, label=L("面積読み出し ω(X_k, X_{k+1})（保存）", "area readout w(X_k, X_{k+1}) (conserved)"))
    ax.set_xlabel("k"); ax.grid(alpha=0.3)
    ax.legend(loc="upper right", prop=FP if JA else None, fontsize=9)
    ax.set_title(L("(d) (a) の軌道に沿った二つの読み出し（初期値で規格化）", "(d) two readouts along orbit (a), normalised"), **tkw(11))
    fig.tight_layout()
    fig.savefig(FIG_DIR / (FIG1 + ".png"))
    plt.close(fig)


def make_fig2(rhos, curves):
    fig, ax = plt.subplots(figsize=(9.5, 6.4), dpi=170)
    for name, ns in curves.items():
        y = np.array([np.nan if v is None else v for v in ns], dtype=float)
        ax.step(rhos, y, where="post", lw=1.4, label=name)
    delta = np.arcsin(np.minimum(1.0, 1.0 / (2.0 * rhos))) / math.pi
    ax.plot(rhos, 1.0 / delta, "k--", lw=1.0, label=L("上界 1/δ ≈ 2πC/ε", "upper bound 1/delta ~ 2 pi C/eps"))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(L("ρ = C/ε（曲率半径 ÷ 分解能）", "rho = C/eps"), **tkw(12))
    ax.set_ylabel(L("初回帰還数 n（＝有効な閉路数）", "first-return number n"), **tkw(12))
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop=FP if JA else None, fontsize=9, loc="upper left")
    ax.set_title(L("有限分解能のもとでの閉路数：n は θ/2π の連分数近似分母だけを取る",
                   "closure number under finite resolution: n takes only convergent denominators"), **tkw(12))
    fig.tight_layout()
    fig.savefig(FIG_DIR / (FIG2 + ".png"))
    plt.close(fig)


def make_fig3(res5):
    ns = sorted(int(k) for k in res5["per_n"].keys())
    mx = [res5["per_n"][n]["max_rel_err"] for n in ns]
    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=170)
    ax.semilogy(ns, mx, "o-", label=L("n≥3：相対誤差の最大値（200 試行）", "n>=3: max relative error (200 trials)"))
    ax.axhline(1e-10, color="k", ls="--", lw=0.9, label=L("事前に宣言した許容値 1e-10", "pre-declared tolerance 1e-10"))
    ax.set_xticks([2] + ns)
    ax.annotate(L("n=2：S=−I で面積の自己読み出しが\n恒等的にゼロ（計量が定義できない）",
                  "n=2: S=-I, area self-readout vanishes\nidentically (no metric)"),
                xy=(2, 3e-17), xytext=(1.65, 6e-17), **tkw(10))
    ax.set_xlim(1.5, 12.5)
    ax.set_ylim(1e-17, 1e-8)
    ax.set_xlabel(L("閉路数 n", "closure number n"), **tkw(12))
    ax.set_ylabel(L("計量読み出しの相対誤差", "relative error of metric readout"), **tkw(11))
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop=FP if JA else None, fontsize=9, loc="upper right")
    ax.set_title(L("閉軌道の二次モーメントから計量を読む", "reading the metric from the orbit's second moment"), **tkw(12))
    fig.tight_layout()
    fig.savefig(FIG_DIR / (FIG3 + ".png"))
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    rng = np.random.default_rng(SEED)
    results = {"seed": SEED, "fail_conditions": FAIL_CONDITIONS}

    ex = exact_tests(rng)
    fl = float_tests(rng)
    rl = relabel_tests(rng)
    cl, curves, rhos = closure_tests()
    mr = metric_readout_tests(rng)
    pa = projective_attractor_tests(rng)
    results.update({"exact": ex, "float": fl, "relabel": rl, "closure": cl,
                    "metric_readout": mr, "projective_attractor": pa})

    rot = rl["rotation_relabel"]
    verdict = {
        "E1a": ex["E1a_violations"] == 0,
        "E1b": fl["E1b_max_scaled_drift"] <= 1e-9,
        "E1c": fl["E1c_signature_mismatches"] == 0,
        "E1d": fl["E1d_naive_sumsq_fluct_median"] >= 0.01,
        "E2a": fl["E2a_max_rel_residual"] <= 1e-9 and ex["second_order_law_violations_exact"] == 0,
        "E2b": all(v["|PSP-S^-1|"] <= 1e-12 and v["|PSP-S|"] >= 1e-3 for v in rot.values()),
        "E2c": rl["E2c_max_imag_eig_perm_equivariant"] <= 1e-12,
        "E3a": cl["E3a_violations"] == 0,
        "E3b": cl["E3b_violations"] == 0,
        "E4a": ex["E4a_violations"] == 0,
        "E4b": fl["E4b_max_curvature_identity_error"] <= 1e-7,
        "E4c": fl["E4c_max_proper_length_ratio_error"] <= 1e-7,
        "E4b_prime": fl["E4b_prime_max_error"] <= 1e-9,
        "E5a": all(v["max_rel_err"] <= 1e-10 for v in mr["per_n"].values()),
        "E5b": mr["E5b_worst_nonmultiple_power_sum"] <= 1e-10,
        "E5c": mr["n2_|G|"] <= 1e-12,
        "E6": pa["E6_max_rel_slope_error"] <= 0.01,
    }
    results["verdict"] = {k: ("PASS" if v else "FAIL") for k, v in verdict.items()}

    d1 = figure1_data()
    svg_fig1(d1)
    svg_fig2(rhos, curves)
    svg_fig3(mr)
    if HAVE_MPL:
        make_fig1(d1)
        make_fig2(rhos, curves)
        make_fig3(mr)
    results["figures_svg"] = [f + ".svg" for f in (FIG1, FIG2, FIG3)]
    results["figures_png"] = [f + ".png" for f in (FIG1, FIG2, FIG3)] if HAVE_MPL else []
    results["png_japanese_font_found"] = JA

    with open(RESULTS_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    for k in sorted(verdict):
        print(f"{k}: {'PASS' if verdict[k] else 'FAIL'}  -- falls if: {FAIL_CONDITIONS[k]}")
    print(json.dumps({k: results[k] for k in ("exact", "float")}, ensure_ascii=False, indent=2))
    print(json.dumps(results["closure"]["per_alpha"], ensure_ascii=False, indent=2))
    print(json.dumps(results["metric_readout"], ensure_ascii=False, indent=2))
    print(json.dumps(pa, indent=2))
    print(json.dumps({"E2c": rl["E2c_max_imag_eig_perm_equivariant"],
                      "strict_elliptic": rl["strict_nameless_elliptic_count"],
                      "rot_n3": rot[3], "rot_n12": rot[12]}, indent=2))


if __name__ == "__main__":
    main()
