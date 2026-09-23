#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
readout_binary_pn.py  (v1.0.0)
数値実験02の読出し: 観測1（2天体の天球上の位置）と観測2（重力波の2偏光）から、
G = c = M = 1 の2体の状態を読み出し、最後に正解と照合する。

読み出す量
  観測者 c から見た量: 質量比 q とどちらが重いか、ν、軌道の形 a, b（M 単位）、
                       M の時間換算 T_M [s]、距離 D [M]、相対距離 r(t) [M]、経過時間 t [M]
  A, B それぞれの量  : 固有時間 τ_A(t), τ_B(t) [M]、
                       自分の時計で測るレーダー距離 d = r dτ/dt [M]（シャピロ遅延を除く1PN近似）

読出しの原理（物理法則は既知、条件の値は使わない）
  1. 観測1の重心: X1 = C + (m2 + 1PN補正) (X1 - X2) → 質量比と、どちらが重いか
  2. 近点通過ごとの特徴量の列を作る
       観測1: 近点時刻、近点の方向、近点距離、遠点距離（相対軌道への局所的な円錐曲線の当てはめ）
       観測2: 振幅 |h+ - i hx| の極大の時刻、その時の位相、極大値、極小値
  3. 生成と同じ運動方程式（2PN + 2.5PN）のモデル軌道から同じ特徴量の列を計算し、観測の列に当てはめる
       非線形パラメータ: a, e（, ν）  線形パラメータ: 時刻の原点, T_M, 方向の原点, 尺度（角度/M または 1/D）
       近点移動の大きさ → a、近点と遠点の比 → e、周期の減少 → ν、周期の秒数 → T_M、角度の大きさ・振幅 → D
  4. 読み出したパラメータでモデル軌道を再構成し、r(t), t, τ_A, τ_B, d_A, d_B を計算する
  読出しは data/input だけを使う。正解 data/truth は最後の照合（evaluate）でだけ読む。

使い方
  python readout_binary_pn.py [--gen 生成フォルダ] [--data データフォルダ] [--out 出力フォルダ]
  既定: 生成フォルダ = ../強い場の二体軌道_観測1観測2データ生成_v1、データ = その data/、出力 = ./results
"""
import argparse
import csv
import json
import math
import pickle
import sys
import time
from pathlib import Path

import numpy as np

VERSION = "1.0.0"
HERE = Path(__file__).resolve().parent
GEN_DIR_DEFAULT = HERE.parent / "強い場の二体軌道_観測1観測2データ生成_v1"
MODEL_STEPS_PER_ORBIT = 500      # モデル計算の刻み（生成の 4000 より粗いが、積分誤差は 1e-7 以下で観測誤差より十分小さい）
W_PERI, W_APO = 1.4, 1.2         # 近点・遠点付近の当てはめに使う位置角の半幅 [rad]
FD_STEP = (1.0e-4, 1.0e-4, 1.0e-4)   # 数値微分の刻み（ln a, e, ν）


# ------------------------------------------------------------------ 入出力

def load_generator(gen_dir):
    sys.path.insert(0, str(Path(gen_dir).resolve()))
    import generate_binary_pn_observations as gen
    return gen


def read_csv(path):
    arr = np.genfromtxt(path, delimiter=",", names=True)
    return {name: np.asarray(arr[name], dtype=float) for name in arr.dtype.names}


def write_csv(path, header, rows, digits=9):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow(["" if (isinstance(v, float) and not math.isfinite(v)) else
                        (f"{v:.{digits}g}" if isinstance(v, float) else v) for v in row])


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return float(f"{v:.7g}") if math.isfinite(v) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


# ------------------------------------------------------------------ 特徴量の抽出（観測データ）

def period_in_samples(y):
    """自己相関の最初の山から、動径周期のおおよそのサンプル数を求める。"""
    y = np.asarray(y) - np.mean(y)
    n = y.size
    f = np.fft.rfft(y, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac = ac / ac[0]
    zc = int(np.argmax(ac < 0.0))
    return zc + int(np.argmax(ac[zc: n // 2]))


def extrema_indices(y, n_per, kind):
    """周期の 1/3 の窓の中で最小（最大）になるサンプル。端で切れたものは使わない。"""
    ys = np.convolve(y, np.ones(3) / 3.0, mode="same")
    w = max(2, n_per // 3)
    out = []
    for i in range(w, y.size - w):
        seg = ys[i - w: i + w + 1]
        if (kind == "min" and ys[i] <= seg.min()) or (kind == "max" and ys[i] >= seg.max()):
            if not out or i - out[-1] > w:
                out.append(i)
    return np.array(out, dtype=int)


def taper(u):
    """端をなめらかに 0 に落とす窓（biweight）。"""
    return np.clip(1.0 - u * u, 0.0, None) ** 2


def local_conic(phi, r, center, half_width, kappa, apo=False):
    """1/r = A + B cos(κφ) + C sin(κφ) を、中心 center の前後 half_width の点に重みつき最小二乗で当てはめる。
    κ = 1/(1 + k) は近点移動を含む円錐曲線（k: 1動径周期あたりの近点移動の割合）。
    重みは 1/r の誤差（σ/r^2）と端をなめらかに落とす窓。求めた極値の方向を中心にもう一度当てはめる
    （標本の位置で結果が飛ばないようにするため）。"""
    for _ in range(2):
        u = (phi - center) / half_width
        sel = np.abs(u) < 1.0
        ph, rr = phi[sel], r[sel]
        w = rr ** 2 * np.sqrt(taper(u[sel]))
        mat = np.column_stack([np.ones_like(ph), np.cos(kappa * ph), np.sin(kappa * ph)]) * w[:, None]
        a_, b_, c_ = np.linalg.lstsq(mat, w / rr, rcond=None)[0]
        amp = math.hypot(b_, c_)
        psi = math.atan2(c_, b_) + (math.pi if apo else 0.0)     # 1/r が最大（近点）または最小（遠点）の位相
        m = round((kappa * center - psi) / (2.0 * math.pi))
        center = (psi + 2.0 * math.pi * m) / kappa
    return {"r_min": 1.0 / (a_ + amp), "r_max": 1.0 / (a_ - amp), "phi_ext": center, "n_points": int(sel.sum())}


def value_at(x, y, x0, n_pts=6, deg=3):
    """x0 に近い n_pts 点の多項式で y(x0) を求める（x は単調）。"""
    idx = np.sort(np.argsort(np.abs(x - x0))[:n_pts])
    scale = np.ptp(x[idx]) or 1.0
    return float(np.polyval(np.polyfit((x[idx] - x0) / scale, y[idx], deg), 0.0))


def poly_extremum(t, y, i0, half, deg, kind):
    """i0 の付近に多項式を当てはめて極値の時刻と値を求める。窓の端はなめらかに落とし、
    求めた極値の時刻を中心にもう一度当てはめる。"""
    lo, hi = max(1, i0 - 2), min(t.size - 1, i0 + 3)
    width = (half + 0.5) * float(np.median(np.diff(t[lo - 1:hi])))
    center, s0, c = float(t[i0]), 0.0, None
    for _ in range(2):
        u = (t - center) / width
        sel = np.abs(u) < 1.0
        c = np.polyfit(u[sel], y[sel], deg, w=np.sqrt(taper(u[sel])))
        roots = np.roots(np.polyder(c))
        roots = roots[np.isreal(roots)].real
        roots = roots[np.abs(roots) < 0.8]
        curv = np.polyval(np.polyder(c, 2), roots)
        roots = roots[(curv < 0) if kind == "max" else (curv > 0)]
        s0 = float(roots[np.argmin(np.abs(roots))]) if roots.size else 0.0
        center += s0 * width
    return center, float(np.polyval(c, s0))


def features_obs1(d):
    """観測1から、重心の回帰に使う量と近点通過ごとの特徴量を作る。"""
    t = d["t_s"]
    x1, y1, x2, y2 = d["body1_x_nrad"], d["body1_y_nrad"], d["body2_x_nrad"], d["body2_y_nrad"]
    rx, ry = x1 - x2, y1 - y2
    r = np.hypot(rx, ry)
    phi = np.unwrap(np.arctan2(ry, rx))
    sense = 1.0 if phi[-1] > phi[0] else -1.0     # 回る向き（鏡映）は読めないので、正の向きにそろえる
    phi = sense * phi
    n_per = period_in_samples(r)
    ip = extrema_indices(r, n_per, "min")
    ia = extrema_indices(r, n_per, "max")
    kappa = 1.0
    for _ in range(2):                            # 1回目で近点移動を見積もり、2回目は κ を入れて当てはめ直す
        rows = []
        for j, i0 in enumerate(ip):
            f = local_conic(phi, r, phi[i0], W_PERI, kappa)
            tp = value_at(phi, t, f["phi_ext"])
            upper = ip[j + 1] if j + 1 < ip.size else r.size
            cand = ia[(ia > i0) & (ia < upper)]
            ra = local_conic(phi, r, phi[cand[0]], W_APO, kappa, apo=True)["r_max"] if cand.size else math.nan
            rows.append((tp, f["phi_ext"], f["r_min"], ra))
        rows = np.array(rows)
        kprec = np.median(np.diff(rows[:, 1])) / (2.0 * math.pi) - 1.0
        kappa = 1.0 / (1.0 + kprec)
    p_est = float(np.median(np.diff(rows[:, 0])))
    n0 = int(round((rows[0, 0] - t[0]) / p_est))
    return {"channel": "obs1", "n": n0 + np.arange(len(rows)), "t": rows[:, 0], "ang": rows[:, 1],
            "peri": rows[:, 2], "apo": rows[:, 3], "kappa": kappa, "sense": sense,
            "samples_per_period": n_per,
            "rho": (rx, sense * ry), "x1": (x1, sense * y1), "t_s": t}


def features_obs2(d):
    """観測2から、近点通過ごとの特徴量（振幅の極大の時刻・位相・極大値、その後の極小値）を作る。"""
    t = d["t_s"]
    h = d["h_plus"] - 1j * d["h_cross"]
    amp = np.abs(h)
    ph = -np.unwrap(np.angle(h))
    if ph[-1] < ph[0]:
        ph = -ph
    n_per = period_in_samples(amp)
    ip = extrema_indices(amp, n_per, "max")
    ia = extrema_indices(amp, n_per, "min")
    rows = []
    for j, i0 in enumerate(ip):
        tp, ap = poly_extremum(t, amp, i0, 6, 4, "max")
        php = value_at(t, ph, tp, 6, 3)
        upper = ip[j + 1] if j + 1 < ip.size else amp.size
        cand = ia[(ia > i0) & (ia < upper)]
        aa = poly_extremum(t, amp, cand[0], 10, 4, "min")[1] if cand.size else math.nan
        rows.append((tp, php, ap, aa))
    rows = np.array(rows)
    p_est = float(np.median(np.diff(rows[:, 0])))
    n0 = int(round((rows[0, 0] - t[0]) / p_est))
    return {"channel": "obs2", "n": n0 + np.arange(len(rows)), "t": rows[:, 0], "ang": rows[:, 1],
            "peri": rows[:, 2], "apo": rows[:, 3], "samples_per_period": n_per, "t_s": t}


# ------------------------------------------------------------------ モデル（生成と同じ運動方程式）

def masses_from_nu(nu):
    d = math.sqrt(max(0.0, 1.0 - 4.0 * min(nu, 0.25)))
    return 0.5 * (1.0 + d), 0.5 * (1.0 - d)


def initial_state_fast(gen, a, e, nu):
    """生成と同じ定義（保存系 2PN 運動の折り返し点が a(1±e)）の初期速度を、Illinois 法で求める。"""
    rp, ra = a * (1.0 - e), a * (1.0 + e)
    m_a, m_b = masses_from_nu(nu)

    def f(v):
        return gen.next_apoapsis(rp, v, nu, m_a, m_b, MODEL_STEPS_PER_ORBIT) - ra

    vn = math.sqrt((1.0 + e) / (a * (1.0 - e)))
    lo, hi = 0.9 * vn, vn
    flo, fhi = f(lo), f(hi)
    while flo > 0.0:
        lo *= 0.95
        flo = f(lo)
    while not fhi > 0.0:
        hi *= 1.02
        fhi = f(hi)
    v, side = hi, 0
    for _ in range(100):
        v = (lo * fhi - hi * flo) / (fhi - flo) if math.isfinite(fhi) else 0.5 * (lo + hi)
        fv = f(v)
        if not math.isfinite(fv) or fv > 0.0:
            hi, fhi = v, fv
            if side == 1:
                flo *= 0.5
            side = 1
        else:
            lo, flo = v, fv
            if side == -1 and math.isfinite(fhi):
                fhi *= 0.5
            side = -1
        if abs(fv) < 1.0e-11 * ra or hi - lo < 1.0e-14 * hi:
            break
    return (rp, 0.0, 0.0, v, 0.0, 0.0)


def run_model(gen, a, e, nu, t_end):
    m_a, m_b = masses_from_nu(nu)
    init = initial_state_fast(gen, a, e, nu)
    ts, st, dv, stop = gen.integrate(a, a * math.sqrt(1.0 - e * e), nu, m_a, m_b, t_end, rr=True,
                                     steps_per_orbit=MODEL_STEPS_PER_ORBIT, init=init)
    return {"ts": ts, "st": st, "dv": dv, "stop": stop, "m_a": m_a, "m_b": m_b, "nu": nu}


def model_features(gen, run):
    """モデル軌道の近点・遠点通過（動径速度 = 0）の時刻、方向、距離、接線速度。n = 0 は出発点。"""
    ts, st, dv = run["ts"], run["st"], run["dv"]
    x, y, vx, vy = st[:, 0], st[:, 1], st[:, 2], st[:, 3]
    rdot = x * vx + y * vy
    i = np.where(np.sign(rdot[1:]) != np.sign(rdot[:-1]))[0]
    i = i[i > 0]
    lo, hi, glo = ts[i].copy(), ts[i + 1].copy(), rdot[i].copy()
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        s = gen.interpolate(ts, st, dv, mid)
        g = s[:, 0] * s[:, 2] + s[:, 1] * s[:, 3]
        same = np.sign(g) == np.sign(glo)
        lo, glo, hi = np.where(same, mid, lo), np.where(same, g, glo), np.where(same, hi, mid)
    tv = 0.5 * (lo + hi)
    s = gen.interpolate(ts, st, dv, tv)
    rv = np.hypot(s[:, 0], s[:, 1])
    vt = np.abs(s[:, 0] * s[:, 3] - s[:, 1] * s[:, 2]) / rv
    phi_steps = np.unwrap(np.arctan2(y, x))
    ang = phi_steps[i] + np.angle(np.exp(1j * (np.arctan2(s[:, 1], s[:, 0]) - phi_steps[i])))
    peri = rdot[i] < 0.0
    tp = np.concatenate(([0.0], tv[peri]))
    out = {"t": tp, "ang": np.concatenate(([0.0], ang[peri])),
           "rp": np.concatenate(([st[0, 0]], rv[peri])), "vtp": np.concatenate(([st[0, 3]], vt[peri])),
           "ra": np.full(tp.size, np.nan), "vta": np.full(tp.size, np.nan)}
    ta, ra_, va_ = tv[~peri], rv[~peri], vt[~peri]
    for n in range(tp.size):
        upper = tp[n + 1] if n + 1 < tp.size else np.inf
        k = np.where((ta > tp[n]) & (ta < upper))[0]
        if k.size:
            out["ra"][n], out["vta"][n] = ra_[k[0]], va_[k[0]]
    nu = run["nu"]
    out["Ap"] = 2.0 * nu * (out["vtp"] ** 2 + 1.0 / out["rp"])   # |h+ - i hx| × D（近点）
    out["Aa"] = 2.0 * nu * (out["vta"] ** 2 + 1.0 / out["ra"])   # 同（遠点）
    return out


# ------------------------------------------------------------------ 特徴量の列の当てはめ

def pseudo_features(gen, run, chan, t_model, t_extra):
    """モデル軌道を観測と同じ時刻（モデルの時間）で標本化し、観測と同じ手順で特徴量を取り出す。
    こうすると、特徴量の取り出し方に由来する偏りは観測とモデルで打ち消し合う。"""
    dt = np.median(np.diff(t_model))
    tq = np.concatenate([t_model, t_model[-1] + dt * np.arange(1, int(t_extra / dt) + 1)])
    tq = tq[(tq >= 0.0) & (tq <= run["ts"][-1])]
    y = gen.interpolate(run["ts"], run["st"], run["dv"], tq)
    x, yy, vx, vy = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
    if chan == "obs1":
        zero = np.zeros_like(x)
        return features_obs1({"t_s": tq, "body1_x_nrad": x, "body1_y_nrad": yy,
                              "body2_x_nrad": zero, "body2_y_nrad": zero})
    r = np.hypot(x, yy)
    h = 2.0 * run["nu"] * ((vx - 1j * vy) ** 2 - (x - 1j * yy) ** 2 / r ** 3)   # h+ - i hx（D = 1）
    return features_obs2({"t_s": tq, "h_plus": h.real, "h_cross": -h.imag})


def match_sequences(feat, pf):
    """観測の近点番号 n に対応するモデルの特徴量の列。そろわなければ None。"""
    pos = {int(k): j for j, k in enumerate(pf["n"])}
    if any(int(k) not in pos for k in feat["n"]):
        return None
    j = np.array([pos[int(k)] for k in feat["n"]])
    return pf["t"][j], pf["ang"][j], pf["peri"][j], pf["apo"][j]


def linear_solution(feat, seq, sig):
    """線形パラメータ（時刻の原点 t0, T_M, 方向の原点, 尺度）の重みつき最小二乗解。"""
    tm, am, pm, qm = seq
    mat = np.column_stack([np.ones_like(tm), tm])
    t0, tscale = np.linalg.lstsq(mat, feat["t"], rcond=None)[0]
    a0 = float(np.mean(feat["ang"] - am))
    ok = np.isfinite(feat["apo"]) & np.isfinite(qm)
    num = np.sum(feat["peri"] * pm) + np.sum(feat["apo"][ok] * qm[ok])
    den = np.sum(pm * pm) + np.sum(qm[ok] ** 2)
    return np.array([t0, tscale, a0, num / den])


def residual_blocks(feat, seq, beta):
    tm, am, pm, qm = seq
    t0, tscale, a0, scale = beta
    ok = np.isfinite(feat["apo"])
    qa = np.where(np.isfinite(qm[ok]), qm[ok], 1.0e3)
    return [feat["t"] - t0 - tscale * tm, feat["ang"] - a0 - am,
            feat["peri"] - scale * pm, feat["apo"][ok] - scale * qa]


class SequenceFit:
    """モデルの特徴量の列を観測の列に当てはめる（Levenberg-Marquardt、線形パラメータは消去）。"""

    def __init__(self, gen, feat, nu_fixed=None):
        self.gen, self.feat, self.nu_fixed = gen, feat, nu_fixed
        self.n = feat["n"].astype(int)
        self.nmax = int(self.n.max())
        self.sig = np.array([1.0e-4 * np.median(np.diff(feat["t"])), 2.0e-3,
                             1.0e-3 * np.nanmedian(feat["peri"]), 1.0e-3 * np.nanmedian(feat["peri"])])
        self.calls = 0
        self.size = 1
        self.tmap = None

    def params(self, theta):
        nu = self.nu_fixed if self.nu_fixed is not None else theta[2]
        return math.exp(theta[0]), theta[1], nu

    def model(self, theta):
        a, e, nu = self.params(theta)
        t0, tscale = self.tmap
        t_model = (self.feat["t_s"] - t0) / tscale
        period = 2.0 * math.pi * a ** 1.5
        t_end = max((self.nmax + 1.6) * period * 1.1, t_model[-1] + 2.0 * period)
        run = run_model(self.gen, a, e, nu, t_end)
        self.calls += 1
        if run["ts"][-1] < t_model[-1]:
            return None, run
        pf = pseudo_features(self.gen, run, self.feat["channel"], t_model, 1.5 * period)
        return match_sequences(self.feat, pf), run

    def init_time_map(self, theta):
        """初期値のモデル軌道の近点時刻と観測の近点時刻の直線関係から、時刻の原点と T_M の初期値を決める。"""
        a, e, nu = self.params(theta)
        run = run_model(self.gen, a, e, nu, (self.nmax + 1.6) * 2.0 * math.pi * a ** 1.5 * 1.1)
        mf = model_features(self.gen, run)
        n = self.n[self.n < mf["t"].size]
        tm = mf["t"][n]
        c = np.polyfit(tm, self.feat["t"][:n.size], 1)
        self.tmap = (float(c[1]), float(c[0]))

    def residuals(self, theta, beta=None):
        seq, run = self.model(theta)
        if seq is None:                             # モデル軌道が観測の周回数に届かない
            return np.full(self.size, 1.0e6), None, None
        if beta is None:
            beta = linear_solution(self.feat, seq, self.sig)
        blocks = residual_blocks(self.feat, seq, beta)
        res = np.concatenate([b / s for b, s in zip(blocks, self.sig)])
        self.size = res.size
        return res, beta, {"seq": seq, "run": run, "blocks": blocks}

    def clamp(self, theta):
        th = theta.copy()
        th[1] = min(max(th[1], 0.02), 0.85)
        if th.size > 2:
            th[2] = min(max(th[2], 0.01), 0.27)    # ν は運動方程式の係数として 0.25 を少し超えることを許す
        return th

    def lm(self, theta, max_iter=15):
        res, beta, info = self.residuals(theta)
        chi2, lam = float(res @ res), 1.0e-3
        for _ in range(max_iter):
            jac = np.empty((res.size, theta.size))
            for j in range(theta.size):
                th = theta.copy()
                th[j] += FD_STEP[j]
                jac[:, j] = (self.residuals(th)[0] - res) / FD_STEP[j]
            amat, grad = jac.T @ jac, jac.T @ res
            improved = False
            for _trial in range(6):
                step = -np.linalg.solve(amat + lam * np.diag(np.diag(amat)), grad)
                th_new = self.clamp(theta + step)
                r_new, b_new, i_new = self.residuals(th_new)
                c_new = float(r_new @ r_new)
                if c_new < chi2:
                    rel = (chi2 - c_new) / chi2
                    theta, res, beta, info, chi2 = th_new, r_new, b_new, i_new, c_new
                    lam = max(lam / 10.0, 1.0e-9)
                    improved = True
                    break
                lam *= 10.0
            if not improved or rel < 1.0e-10:
                break
        return theta, res, beta, info, chi2

    def solve(self, theta0):
        theta = self.clamp(np.array(theta0, dtype=float))
        self.init_time_map(theta)
        theta, res, beta, info, chi2 = self.lm(theta)
        for _ in range(2):      # 時刻の対応と各特徴量の誤差（残差の大きさ）を更新して当てはめ直す
            self.tmap = (beta[0], beta[1])
            self.sig = np.array([max(np.sqrt(np.mean(b ** 2)), 1.0e-30) for b in info["blocks"]])
            theta, res, beta, info, chi2 = self.lm(theta)
        # 共分散: 非線形パラメータ（線形パラメータを固定して数値微分）と線形パラメータの全体のヤコビ行列
        cols = []
        for j in range(theta.size):
            th = theta.copy()
            th[j] += FD_STEP[j]
            cols.append((self.residuals(th, beta)[0] - res) / FD_STEP[j])
        tm, am, pm, qm = info["seq"]
        ok = np.isfinite(self.feat["apo"])
        nt, na, npp, nq = tm.size, am.size, pm.size, int(ok.sum())
        z = np.zeros
        lin = [np.concatenate([-np.ones(nt) / self.sig[0], z(na + npp + nq)]),
               np.concatenate([-tm / self.sig[0], z(na + npp + nq)]),
               np.concatenate([z(nt), -np.ones(na) / self.sig[1], z(npp + nq)]),
               np.concatenate([z(nt + na), -pm / self.sig[2], -qm[ok] / self.sig[3]])]
        jac = np.column_stack(cols + lin)
        norm = np.sqrt(np.sum(jac ** 2, axis=0))
        jn = jac / norm
        dof = max(1, res.size - jac.shape[1])
        cov = np.linalg.inv(jn.T @ jn) / np.outer(norm, norm) * (chi2 / dof)
        a, e, nu = self.params(theta)
        return {"theta": theta, "a": a, "e": e, "nu": nu, "beta": beta, "cov": cov,
                "chi2_per_dof": chi2 / dof, "sigmas": self.sig, "info": info, "model_calls": self.calls}


# ------------------------------------------------------------------ 初期値

def precession_to_a(kprec, e, nu):
    """2PN の近点移動の式（Damour & Schafer 1988）から a のおおよその値を求める（初期値用）。"""
    c1 = 3.0 / (1.0 - e * e)
    c2 = ((78.0 - 28.0 * nu) + (51.0 - 26.0 * nu) * e * e) / (4.0 * (1.0 - e * e) ** 2)
    x = (-c1 + math.sqrt(c1 * c1 + 4.0 * c2 * kprec)) / (2.0 * c2)
    return 1.0 / x


def initial_guess(feat, nu_guess):
    n = feat["n"][:6]
    slope = np.polyfit(n, feat["ang"][:6], 1)[0]
    if feat["channel"] == "obs1":
        kprec = slope / (2.0 * math.pi) - 1.0
        e = (feat["apo"][0] - feat["peri"][0]) / (feat["apo"][0] + feat["peri"][0])
    else:
        kprec = slope / (4.0 * math.pi) - 1.0
        ratio = feat["peri"][0] / feat["apo"][0]    # (1+e)(2+e)/((1-e)(2-e))（ニュートン近似）
        e = (3.0 * (ratio + 1.0) - math.sqrt(9.0 * (ratio + 1.0) ** 2 - 8.0 * (ratio - 1.0) ** 2)) / (2.0 * (ratio - 1.0))
    a = precession_to_a(kprec, e, 0.2 if nu_guess is None else nu_guess)
    if nu_guess is None:
        c = np.polyfit(feat["n"], feat["t"], 2)
        pdot = 2.0 * c[0] / (c[1] + 2.0 * c[0] * feat["n"][0])
        pm = 2.0 * math.pi * a ** 1.5
        fe = (1.0 + 73.0 / 24.0 * e * e + 37.0 / 96.0 * e ** 4) / (1.0 - e * e) ** 3.5
        nu_guess = min(0.25, max(0.02, -pdot / (192.0 * math.pi / 5.0 * (2.0 * math.pi / pm) ** (5.0 / 3.0) * fe)))
        a = precession_to_a(kprec, e, nu_guess)
    return [math.log(a), e, nu_guess], kprec


# ------------------------------------------------------------------ 重心と質量比（観測1）

def com_mass_fraction(f1, com_factor=None, iterations=4):
    """X1 = C + (m2 + ν (m1 - m2) f(t)) ρ を解いて m2（観測1の2番目の天体の質量の割合）を求める。
    f(t) = (v^2 - 1/r)/2 は1PNの重心補正（com_factor が None ならニュートン近似）。"""
    rx, ry = f1["rho"]
    x1, y1 = f1["x1"]
    ones, zeros = np.ones_like(rx), np.zeros_like(rx)
    mat = np.vstack([np.column_stack([ones, zeros, rx]), np.column_stack([zeros, ones, ry])])
    mu2 = 0.5
    for _ in range(iterations if com_factor is not None else 1):
        corr = 0.0 if com_factor is None else mu2 * (1.0 - mu2) * (1.0 - 2.0 * mu2) * com_factor
        rhs = np.concatenate([x1 - corr * rx, y1 - corr * ry])
        sol, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
        mu2 = float(sol[2])
    resid = rhs - mat @ sol
    cov = np.linalg.inv(mat.T @ mat) * (resid @ resid) / (rhs.size - 3)
    return mu2, math.sqrt(cov[2, 2])


def model_at(gen, run, t_model):
    tq = np.clip(t_model, 0.0, run["ts"][-1])
    y = gen.interpolate(run["ts"], run["st"], run["dv"], tq)
    return {k: y[:, i] for i, k in enumerate(("x", "y", "vx", "vy", "tau_heavy", "tau_light"))}


# ------------------------------------------------------------------ 読出し（正解は使わない）

def readout_dataset(gen, data_dir, name, cache_dir):
    """読出し（正解は使わない）。当てはめごとに cache/ に保存し、途中で止まっても続きから再開する。"""
    f1 = features_obs1(read_csv(data_dir / "input" / f"obs1_{name}.csv"))
    f2 = features_obs2(read_csv(data_dir / "input" / f"obs2_{name}.csv"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{name}.pkl"
    state = pickle.loads(path.read_bytes()) if path.exists() else {}

    def step(key, func):
        if key not in state:
            t0 = time.time()
            value = func()
            if isinstance(value, dict) and "info" in value:
                value["info"].pop("run", None)
            state[key], state["time_" + key] = value, time.time() - t0
            path.write_bytes(pickle.dumps(state))
            print(f"[{name}] {key}: {state['time_' + key]:.0f} s", flush=True)
        return state[key]

    # 観測1: 重心（ニュートン近似）→ ν を固定して a, e, T_M, 尺度 → 1PN の重心補正で質量比を読み直して再当てはめ
    mu2_n = step("com_newtonian", lambda: com_mass_fraction(f1)[0])
    nu0 = mu2_n * (1.0 - mu2_n)
    th1, kp1 = initial_guess(f1, nu0)
    fit1a = step("fit_obs1_first", lambda: SequenceFit(gen, f1, nu_fixed=nu0).solve(th1[:2]))

    def com_1pn():
        t0_, ts_ = fit1a["beta"][0], fit1a["beta"][1]
        run = run_model(gen, fit1a["a"], fit1a["e"], fit1a["nu"], (f1["t_s"][-1] - t0_) / ts_ + 50.0)
        m = model_at(gen, run, (f1["t_s"] - t0_) / ts_)
        return com_mass_fraction(f1, 0.5 * (m["vx"] ** 2 + m["vy"] ** 2 - 1.0 / np.hypot(m["x"], m["y"])))

    mu2, mu2_sig = step("com_1pn", com_1pn)
    nu_com = mu2 * (1.0 - mu2)
    fit1 = step("fit_obs1", lambda: SequenceFit(gen, f1, nu_fixed=nu_com).solve(fit1a["theta"]))
    # 観測1: ν も自由にして、軌道（周期の減少）だけから ν を読む
    fit1n = step("fit_obs1_free_nu", lambda: SequenceFit(gen, f1).solve(list(fit1["theta"]) + [nu_com]))
    # 観測2: 重力波だけから a, e, ν を読む
    th2, kp2 = initial_guess(f2, None)
    fit2 = step("fit_obs2", lambda: SequenceFit(gen, f2).solve(th2))
    return {"dataset": name, "features": {"obs1": f1, "obs2": f2}, "mu2_newtonian": mu2_n, "mu2": mu2,
            "mu2_sigma": mu2_sig, "nu_com": nu_com, "fit_obs1": fit1, "fit_obs1_free_nu": fit1n,
            "fit_obs2": fit2, "initial_precession": {"obs1": kp1, "obs2": kp2},
            "elapsed_s": sum(v for k, v in state.items() if k.startswith("time_"))}


def summarize_fit(fit, chan):
    sd = np.sqrt(np.diag(fit["cov"]))
    ntheta = fit["theta"].size
    a, e = fit["a"], fit["e"]
    s = {"a_M": a, "a_sigma": a * sd[0], "e": e, "e_sigma": sd[1],
         "b_M": a * math.sqrt(1.0 - e * e),
         "nu": fit["nu"], "nu_sigma": sd[2] if ntheta > 2 else None,
         "t_ref_s": fit["beta"][0], "t_ref_sigma": sd[ntheta],
         "T_M_s": fit["beta"][1], "T_M_sigma": sd[ntheta + 1],
         "chi2_per_dof": fit["chi2_per_dof"], "feature_sigmas": fit["sigmas"],
         "model_calls": fit["model_calls"]}
    cab = fit["cov"][0, 1] * a
    s["b_sigma"] = math.sqrt(max(0.0, (math.sqrt(1 - e * e) * s["a_sigma"]) ** 2
                                 + (a * e / math.sqrt(1 - e * e) * s["e_sigma"]) ** 2
                                 - 2.0 * math.sqrt(1 - e * e) * a * e / math.sqrt(1 - e * e) * cab))
    scale, ssd = fit["beta"][3], sd[ntheta + 3]
    if chan == "obs1":
        s.update({"scale_nrad_per_M": scale, "scale_sigma": ssd, "D_M": 1.0e9 / scale, "D_sigma": 1.0e9 / scale * ssd / scale})
    else:
        s.update({"strain_scale_1_over_D": scale, "D_M": 1.0 / scale, "D_sigma": ssd / scale ** 2})
    return s


def q_from_nu(nu):
    nu = min(nu, 0.25)
    return (1.0 - 2.0 * nu + math.sqrt(max(0.0, 1.0 - 4.0 * nu))) / (2.0 * nu)


def reconstruct(gen, fit, t_s):
    """読み出したパラメータでモデル軌道を再構成し、観測時刻での r, t, τ, レーダー距離を返す。"""
    t0, tscale = fit["beta"][0], fit["beta"][1]
    tm = (t_s - t0) / tscale
    run = run_model(gen, fit["a"], fit["e"], fit["nu"], tm.max() + 50.0)
    m = model_at(gen, run, tm)
    r = np.hypot(m["x"], m["y"])
    v2 = m["vx"] ** 2 + m["vy"] ** 2
    mh, ml = run["m_a"], run["m_b"]
    return {"t_M": tm, "r": r, "tau_heavy": m["tau_heavy"], "tau_light": m["tau_light"],
            "d_heavy": r * (1.0 - ml / r - 0.5 * ml * ml * v2), "d_light": r * (1.0 - mh / r - 0.5 * mh * mh * v2),
            "stop": run["stop"]}


def blind_summary(ro):
    f1 = ro["features"]["obs1"]
    s1, s1n, s2 = (summarize_fit(ro["fit_obs1"], "obs1"), summarize_fit(ro["fit_obs1_free_nu"], "obs1"),
                   summarize_fit(ro["fit_obs2"], "obs2"))
    mu2 = ro["mu2"]
    heavy = "body1" if mu2 < 0.5 else "body2"
    lam = mu2 / (1.0 - mu2)
    q = 1.0 / lam if lam < 1 else lam
    q_sig = ro["mu2_sigma"] / (min(mu2, 1 - mu2) ** 2)
    return {
        "dataset": ro["dataset"],
        "obs1_centre_of_mass": {"m2_fraction": mu2, "m2_fraction_sigma": ro["mu2_sigma"],
                                "m2_fraction_newtonian": ro["mu2_newtonian"], "heavier": heavy,
                                "q_heavy_over_light": q, "q_sigma": q_sig, "nu": ro["nu_com"]},
        "obs1_orbit_fit_nu_fixed_to_com": s1,
        "obs1_orbit_fit_nu_free": s1n,
        "obs2_gw_fit": dict(s2, q_from_nu=q_from_nu(s2["nu"])),
        "features": {"obs1_passages": int(f1["n"].size), "obs2_passages": int(ro["features"]["obs2"]["n"].size),
                     "obs1_kappa": f1["kappa"], "initial_precession_estimate": ro["initial_precession"]},
        "elapsed_s": ro["elapsed_s"],
    }


# ------------------------------------------------------------------ 照合（ここで初めて正解を読む）

def rel(read, true):
    return (read - true) / true


def evaluate_dataset(gen, data_dir, ro, blind):
    name = ro["dataset"]
    meta = json.loads((data_dir / "truth" / f"{name}_meta.json").read_text(encoding="utf-8"))
    tr = read_csv(data_dir / "truth" / f"{name}_timeseries.csv")
    sel = tr["in_obs1"] == 1
    tr = {k: v[sel] for k, v in tr.items()}
    hid, orb, mas = meta["hidden_observer_parameters"], meta["orbit_shape_parameters"], meta["masses"]
    f1 = ro["features"]["obs1"]
    t_s = f1["t_s"]
    v2 = tr["vx_rel"] ** 2 + tr["vy_rel"] ** 2
    ma, mb = mas["m_A"], mas["m_B"]
    d_a_true = tr["r"] * (1.0 - mb / tr["r"] - 0.5 * mb * mb * v2)
    d_b_true = tr["r"] * (1.0 - ma / tr["r"] - 0.5 * ma * ma * v2)
    rec1, rec2 = reconstruct(gen, ro["fit_obs1"], t_s), reconstruct(gen, ro["fit_obs2"], t_s)
    r_data = np.hypot(*f1["rho"]) / ro["fit_obs1"]["beta"][3]

    rows = []

    def add(quantity, method, read, sigma, true, ref=None):
        # 正解が 0 の量（等質量のときの差）は、基準の量 ref（τ_A や d_A）で割った誤差にする
        err = rel(read, true) if true != 0 else (read - true) / ref
        rows.append({"quantity": quantity, "method": method, "readout": read, "sigma": sigma, "truth": true,
                     "relative_error": err, "normalized_by": "truth" if true != 0 else "reference",
                     "pull": (read - true) / sigma if sigma else None})

    c, s1, s1n, s2 = (blind["obs1_centre_of_mass"], blind["obs1_orbit_fit_nu_fixed_to_com"],
                      blind["obs1_orbit_fit_nu_free"], blind["obs2_gw_fit"])
    add("q", "obs1 centre of mass", c["q_heavy_over_light"], c["q_sigma"], mas["q_mA_over_mB"])
    add("q", "obs2 GW (from nu)", s2["q_from_nu"], None, mas["q_mA_over_mB"])
    add("nu", "obs1 centre of mass", c["nu"], None, mas["nu"])
    add("nu", "obs1 orbit decay", s1n["nu"], s1n["nu_sigma"], mas["nu"])
    add("nu", "obs2 GW", s2["nu"], s2["nu_sigma"], mas["nu"])
    for key, label, true in (("a_M", "a", orb["a"]), ("b_M", "b", orb["b"]), ("e", "e", orb["e"]),
                             ("T_M_s", "T_M", hid["T_M_seconds_per_M"]), ("D_M", "D", hid["distance_D_in_M"])):
        sk = {"a_M": "a_sigma", "b_M": "b_sigma", "e": "e_sigma", "T_M_s": "T_M_sigma", "D_M": "D_sigma"}[key]
        add(label, "obs1", s1[key], s1[sk], true)
        add(label, "obs2 GW", s2[key], s2[sk], true)
    heavy_true = "body1" if hid["body1_is"] == "A" else "body2"

    def rms(x):
        return float(np.sqrt(np.mean(x ** 2)))

    derived = {
        "heavier_body_identified_correctly": None if mas["q_mA_over_mB"] == 1 else c["heavier"] == heavy_true,
        "time_origin_t_ref_minus_t_offset_s": {"obs1": s1["t_ref_s"] - hid["time_offset_s"],
                                               "obs2": s2["t_ref_s"] - hid["time_offset_s"]},
        "r_rms_relative_error": {"obs1_data_scaled": rms(rel(r_data, tr["r"])),
                                 "obs1_model": rms(rel(rec1["r"], tr["r"])), "obs2_model": rms(rel(rec2["r"], tr["r"]))},
        "t_M_end": {"truth": tr["t_M"][-1], "obs1": rec1["t_M"][-1], "obs2": rec2["t_M"][-1]},
    }
    end = -1
    for lab, rec in (("obs1", rec1), ("obs2 GW", rec2)):
        add("tau_A end", lab, rec["tau_heavy"][end], None, tr["tau_A"][end])
        add("tau_B end", lab, rec["tau_light"][end], None, tr["tau_B"][end])
        add("tau_A - tau_B end", lab, rec["tau_heavy"][end] - rec["tau_light"][end], None,
            tr["tau_A"][end] - tr["tau_B"][end], ref=tr["tau_A"][end])
        add("radar d_A mean", lab, float(np.mean(rec["d_heavy"])), None, float(np.mean(d_a_true)))
        add("radar d_B mean", lab, float(np.mean(rec["d_light"])), None, float(np.mean(d_b_true)))
        add("radar d_A - d_B mean", lab, float(np.mean(rec["d_heavy"] - rec["d_light"])), None,
            float(np.mean(d_a_true - d_b_true)), ref=float(np.mean(d_a_true)))
    recon = {"t_s": t_s, "t_M_true": tr["t_M"], "t_M_obs1": rec1["t_M"], "t_M_obs2": rec2["t_M"],
             "r_true": tr["r"], "r_obs1_data": r_data, "r_obs1_model": rec1["r"], "r_obs2_model": rec2["r"],
             "dtau_true": tr["tau_A"] - tr["tau_B"], "dtau_obs1": rec1["tau_heavy"] - rec1["tau_light"],
             "dtau_obs2": rec2["tau_heavy"] - rec2["tau_light"]}
    return {"dataset": name, "case_label": meta["case_label"], "rows": rows, "derived": derived}, recon


# ------------------------------------------------------------------ 出力

def feature_rows(ro):
    rows = []
    for chan, key in (("obs1", "fit_obs1"), ("obs2", "fit_obs2")):
        f, fit = ro["features"][chan], ro[key]
        tm, am, pm, qm = fit["info"]["seq"]
        t0, tscale, a0, scale = fit["beta"]
        for j in range(f["n"].size):
            rows.append([chan, int(f["n"][j]), f["t"][j], f["ang"][j], f["peri"][j], f["apo"][j],
                         t0 + tscale * tm[j], a0 + am[j], scale * pm[j], scale * qm[j]])
    return rows


REC_KEYS = ["t_M_true", "r_true", "r_obs1_data", "r_obs1_model", "r_obs2_model", "dtau_true", "dtau_obs1", "dtau_obs2"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gen", default=str(GEN_DIR_DEFAULT), help="生成プログラムのあるフォルダ")
    ap.add_argument("--data", default=None, help="data フォルダ（既定: 生成フォルダ/data）")
    ap.add_argument("--out", default=str(HERE / "results"), help="出力フォルダ")
    ap.add_argument("--datasets", default="D1,D2,D3")
    args = ap.parse_args()
    gen = load_generator(args.gen)
    data_dir = Path(args.data) if args.data else Path(args.gen) / "data"
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    blind_all, evals, feat_rows, rec_rows = [], [], [], []
    for name in args.datasets.split(","):
        print(f"[{name}] 読出し ...", flush=True)
        ro = readout_dataset(gen, data_dir, name, out_dir / "cache")
        blind = blind_summary(ro)
        blind_all.append(blind)
        feat_rows += [[name] + row for row in feature_rows(ro)]
        ev, recon = evaluate_dataset(gen, data_dir, ro, blind)          # ここで初めて正解を読む
        evals.append(ev)
        idx = np.unique(np.linspace(0, recon["t_s"].size - 1, 101).astype(int))
        rec_rows += [[name] + [float(recon[k][i]) for k in REC_KEYS] for i in idx]
        for row in ev["rows"]:
            print(f"   {row['quantity']:>20s} {row['method']:<22s} read {row['readout']:.7g}  truth {row['truth']:.7g}"
                  f"  rel.err {row['relative_error']:+.2e}" + (f"  pull {row['pull']:+.2f}" if row["pull"] is not None else ""))
    common = {"program": "readout_binary_pn.py", "version": VERSION, "units": "G = c = M = 1",
              "model_steps_per_orbit": MODEL_STEPS_PER_ORBIT}
    (out_dir / "readout_results.json").write_text(
        json.dumps(to_jsonable(dict(common, readout_blind=blind_all, evaluation=evals)), ensure_ascii=False, indent=1),
        encoding="utf-8")
    write_csv(out_dir / "summary.csv", ["dataset", "quantity", "method", "readout", "sigma", "truth",
                                        "relative_error", "normalized_by", "pull"],
              [[ev["dataset"], r["quantity"], r["method"], float(r["readout"]),
                float(r["sigma"]) if r["sigma"] is not None else math.nan, float(r["truth"]),
                float(r["relative_error"]), r["normalized_by"], float(r["pull"]) if r["pull"] is not None else math.nan]
               for ev in evals for r in ev["rows"]], digits=7)
    write_csv(out_dir / "features.csv", ["dataset", "channel", "n", "t_s", "angle_rad", "peri", "apo",
                                         "t_fit", "angle_fit", "peri_fit", "apo_fit"], feat_rows)
    write_csv(out_dir / "reconstruction.csv", ["dataset"] + REC_KEYS, rec_rows, digits=8)
    print("出力:", out_dir)


if __name__ == "__main__":
    main()
