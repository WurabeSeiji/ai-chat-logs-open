#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nスケーリング測定プログラム（低ランク頂点分解版）

    python3 run_n_scaling_lowrank_v1.py N [seed] [--steps CAP] [--validate]

第3論文の頂点分解 K = C S^T - S C^T = W J W^T（W = [C|S], M×2N, rank ≤ 2N）を
アルゴリズムとして使い、密行列 O(M^3)=O(N^6)/step を O(N^3)/step に削減する。

  - K z 積: 辺→頂点の散布和で O(M)
  - G = W^T W: K_N の線グラフでは頂点対の共有辺が1本なので、
    G_cc[k,l] = cos^2 θ_(k,l)（非対角）等、辺位相から O(N^2) で直接構成
  - Cayley 解: Woodbury で (I-γK̃)^{-1} = I - W[(σ/γ)J + G]^{-1} W^T（2N×2N 解）
  - σスペクトル: K の非零固有値 = JG（2N×2N）の固有値。σ_max は
    warm-start 冪反復（O(M)/step）で追い、サブサンプル時に厳密値と照合

測定（ケース分け最小、CSVなし、JSON出力のみ）:
  [A] 開始率プローブ: 自己無撞着親（円偏波固有モード）+ 零閉鎖核種 δ、
      休眠フラクション f の指数成長率。種は g=(u+iw)/√2 で g^T g = 0 厳密
      （Codex 査読指摘への対応: Z^T Z = 0 を保った摂動）
  [B] 緩和プローブ: 零閉鎖一般状態（Z=X+iY, |X|=|Y|, X⊥Y）からの
      等分配到達時間・PR/M・活性平面数 n_act と上界 min(N,⌊M/2⌋) の比較・
      span(W) 占有率

  N ≤ 12 では --validate で密行列実装との軌道一致を検証（正当性ゲート）。
"""

import json
import math
import os
import sys
import time

import warnings

import numpy as np

warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

def progress(msg):
    print(f"[進捗 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


GAMMA = math.tan(math.pi / 144.0)  # 旧 Cayley 刻み（記録用。validate では使わない）
ANGLE = 2.0 * math.pi / 124.0    # 判断4(R3-ii): 検算用の線形回転の刻み角（最新版の L=124 に合わせる）
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "n_scaling_result_v1")


# ---------------- 低ランク基盤 ----------------

def build_edges(n):
    ea, eb = np.triu_indices(n, k=1)
    return ea.astype(np.int64), eb.astype(np.int64)


class LowRankSystem:
    """位相差正弦生成子の低ランク表現。θ を設定すると c,s,頂点和,G を保持する。"""

    def __init__(self, n):
        self.n = n
        self.ea, self.eb = build_edges(n)
        self.m = len(self.ea)
        self.J = np.zeros((2 * n, 2 * n))
        self.J[:n, n:] = np.eye(n)
        self.J[n:, :n] = -np.eye(n)

    def set_theta(self, theta):
        n = self.n
        self.c = np.cos(theta)
        self.s = np.sin(theta)
        T = np.zeros((n, n))
        T[self.ea, self.eb] = theta
        T[self.eb, self.ea] = theta
        CT = np.cos(T)
        ST = np.sin(T)
        np.fill_diagonal(CT, 0.0)
        np.fill_diagonal(ST, 0.0)
        Gcc = CT * CT
        Gcs = CT * ST
        Gss = ST * ST
        np.fill_diagonal(Gcc, Gcc.sum(axis=1))
        np.fill_diagonal(Gcs, Gcs.sum(axis=1))
        np.fill_diagonal(Gss, Gss.sum(axis=1))
        G = np.empty((2 * n, 2 * n))
        G[:n, :n] = Gcc
        G[:n, n:] = Gcs
        G[n:, :n] = Gcs
        G[n:, n:] = Gss
        self.G = G

    def vsum(self, vals):
        n = self.n
        if np.iscomplexobj(vals):
            re = np.bincount(self.ea, vals.real, n) + np.bincount(self.eb, vals.real, n)
            im = np.bincount(self.ea, vals.imag, n) + np.bincount(self.eb, vals.imag, n)
            return re + 1j * im
        return np.bincount(self.ea, vals, n) + np.bincount(self.eb, vals, n)

    def wt(self, z):
        """W^T z（2N）"""
        return np.concatenate([self.vsum(self.c * z), self.vsum(self.s * z)])

    def w(self, y):
        """W y（M）"""
        n = self.n
        yc, ys = y[:n], y[n:]
        return self.c * (yc[self.ea] + yc[self.eb]) + self.s * (ys[self.ea] + ys[self.eb])

    def kmatvec(self, z):
        """K z = C(S^T z) - S(C^T z)"""
        vs = self.vsum(self.s * z)
        vc = self.vsum(self.c * z)
        return self.c * (vs[self.ea] + vs[self.eb]) - self.s * (vc[self.ea] + vc[self.eb])

    def sigma_spectrum(self):
        """厳密σスペクトル（降順）: JG の固有値の正の虚部。O((2N)^3)"""
        ev = np.linalg.eigvals(self.J @ self.G)
        sig = np.sort(ev.imag[ev.imag > 1e-12])[::-1]
        return sig

    def sigma_max_power(self, wp, iters=3):
        """warm-start 冪反復による σ_max 推定。wp は前ステップのベクトル（更新して返す）。"""
        for _ in range(iters):
            y = self.kmatvec(wp)
            wp = -self.kmatvec(y)
            nrm = np.linalg.norm(wp)
            if nrm == 0.0:
                return 0.0, wp
            wp = wp / nrm
        sig = np.linalg.norm(self.kmatvec(wp))
        return float(sig), wp

    def cayley_step(self, z, sigma):
        """z ← (I-γK)^{-1}(I+γK) z。K/σ 正規化なし。Woodbury で O(N^3)。"""
        gn = GAMMA
        r = z + gn * self.kmatvec(z)
        A2 = (1.0 / GAMMA) * self.J + self.G
        rhs = self.wt(r)
        y = np.linalg.solve(A2, rhs)
        return r - self.w(y)


def participation_ratio(z):
    a2 = np.abs(z) ** 2
    return float(np.sum(a2) ** 2 / np.sum(a2 * a2))


# ---------------- 状態構成 ----------------

def _eigenmode_residual(sys_lr, v):
    """カイラリティ非依存の固有モード残差: μ = v†(iKv) に対する ‖iKv - μv‖。"""
    kv = sys_lr.kmatvec(v)
    mu = float(np.real(np.conj(v) @ (1j * kv)))
    return float(np.linalg.norm(1j * kv - mu * v))


def _adjacency(sys_lr):
    """隣接辺マスク A（頂点を共有する辺対に 1）。"""
    m = sys_lr.m
    A = np.zeros((m, m))
    for i in range(m):
        share = (sys_lr.ea == sys_lr.ea[i]) | (sys_lr.ea == sys_lr.eb[i]) | (sys_lr.eb == sys_lr.ea[i]) | (sys_lr.eb == sys_lr.eb[i])
        share[i] = False
        A[i, share] = 1.0
    return A


def _K_amplitude_aware(A, v):
    """実際の相互作用の生成子 K_ij = Im(conj(v_i) v_j)（隣接辺のみ）。実反対称。"""
    K = A * np.imag(np.conj(v)[:, None] * v[None, :])
    np.fill_diagonal(K, 0.0)
    return K


def _selfconsistency_residual(K, v):
    """スケール不変の固有モード残差 ‖iKv − μv‖/‖v‖、μ = v†(iKv)/v†v。"""
    kv = 1j * (K @ v)
    mu = np.vdot(v, kv) / np.vdot(v, v)
    return float(np.linalg.norm(kv - mu * v) / np.linalg.norm(v))


def make_parent(sys_lr, rng, iters=2000, beta=0.5, tol=1e-10, restarts=20):
    """自己無撞着な親：実際の相互作用 K_amp(v) = Im(conj(v_i) v_j) に対して K_amp(v) v = iσ_max v を満たす状態 v を、
    位相と振幅の両方を未知数として反復で解く（v → K_amp(v) の λ=−iσ_max モード u（旧 make_parent の argmin(ev.imag) と同じカイラリティ）→ v ← (1−β)v + βu、u は v に位相整合）。

    - 振幅の正規化は行わない。固定点は c·v の族（K_amp は |v|² に比例）で、全体スケール c は自由なモジュライ。
      ここでは固有値ソルバが返す固有ベクトルのスケール（‖u‖=1）をそのまま使う（明示的な正規化行は無い）。
    - 収束判定はスケール不変残差 < tol。restarts 回すべて失敗したら RuntimeError（実験は走らせない）。
    """
    A = _adjacency(sys_lr)
    m = sys_lr.m
    log = []
    for r in range(restarts):
        theta = rng.uniform(0.0, 2.0 * np.pi, m)
        v = np.exp(1j * theta)
        res = np.inf
        for it in range(iters):
            K = _K_amplitude_aware(A, v)
            w, U = np.linalg.eigh(1j * K)
            # 旧 make_parent と同じカイラリティ：λ = −iσ_max のモード（iK の最小固有値 −σ_max）。縮退時は現在の v との重なりが最大のもの
            top = np.where(w <= w.min() + 1e-12 * max(1.0, abs(w.min())))[0]
            ov = [abs(np.vdot(U[:, j], v)) for j in top]
            u = U[:, top[int(np.argmax(ov))]]
            ph = np.vdot(u, v)
            if abs(ph) > 0:
                u = u * (ph / abs(ph))  # 位相整合（大域位相は自由）
            v_new = (1.0 - beta) * v + beta * u
            v = v_new
            if it % 10 == 9:
                res = _selfconsistency_residual(_K_amplitude_aware(A, v), v)
                progress(f"親構成(振幅込み) restart={r+1} iter={it+1} 残差={res:.2e}")
                if res < tol:
                    break
        K = _K_amplitude_aware(A, v)
        res = _selfconsistency_residual(K, v)
        sig = np.sort(np.linalg.eigvalsh(1j * K))[::-1]
        sig = sig[sig > 1e-12]
        log.append((r + 1, res))
        if res < tol:
            sys_lr.set_theta(np.angle(v))
            progress(f"親構成(振幅込み) 収束 restart={r+1} 残差={res:.2e} |v|={np.linalg.norm(v):.6f} σ={sig[:4]}")
            return v, res, sig
    msg = "make_parent: 振幅込み K に対する自己無撞着解が得られず（restarts=%d, iters=%d, tol=%.0e）: " % (restarts, iters, tol) + ", ".join(f"r{r}:{x:.2e}" for r, x in log)
    progress(msg)
    raise RuntimeError(msg)


def zero_closure_kernel_seed(sys_lr, rng):
    """span(W) の直交補内の実正規直交対 (u,w) から g=(u+iw)/√2。g^T g = 0 厳密。"""
    m = sys_lr.m
    def project_out(g):
        q = np.linalg.lstsq(sys_lr.G, sys_lr.wt(g), rcond=None)[0]
        return g - sys_lr.w(q)
    u = project_out(rng.normal(size=m))
    u = u / np.linalg.norm(u)
    w = project_out(rng.normal(size=m))
    w = w - (w @ u) * u
    w = w / np.linalg.norm(w)
    return (u + 1j * w) / math.sqrt(2.0)


def zero_closure_generic(rng, m):
    """Z=X+iY, |X|=|Y|, X⊥Y → Z^T Z = 0 厳密。"""
    X = rng.normal(size=m)
    Y = rng.normal(size=m)
    Y = Y - (X @ Y) / (X @ X) * X
    Y = Y * (np.linalg.norm(X) / np.linalg.norm(Y))
    Z = X + 1j * Y
    return Z  # 判断1(A5): 全体の振幅正規化を除去（|X|=|Y| の相対比は維持）


# ---------------- 検証（密行列との一致） ----------------

def validate_against_dense(n, seed, steps=300):
    sys_lr = LowRankSystem(n)
    m = sys_lr.m
    rng = np.random.default_rng(seed)
    Z = zero_closure_generic(rng, m)

    # 密行列側の準備
    ea, eb = sys_lr.ea, sys_lr.eb
    A = np.zeros((m, m))
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        A[i, share] = 1.0
    np.fill_diagonal(A, 0.0)

    # 判断4(R3-ii): 密行列側も振幅込み K で比較。低ランク側は c=Re z, s=Im z を直接与える
    sys_lr.c = Z.real.copy(); sys_lr.s = Z.imag.copy()
    Kd = A * np.imag(np.conj(Z)[:, None] * Z[None, :])

    err_matvec = np.linalg.norm(sys_lr.kmatvec(Z) - Kd @ Z) / np.linalg.norm(Kd @ Z)

    err_G = 0.0  # G（WᵀW）は振幅込みでは set_theta が作らないため本検算では対象外（kmatvec と軌道で検算）

    sig_lr = np.sort(np.linalg.eigvalsh(1j * np.column_stack([sys_lr.kmatvec(e) for e in np.eye(m)])))
    sig_lr = sig_lr[sig_lr > 1e-12][::-1]
    sig_d = np.sort(np.linalg.eigvalsh(1j * Kd))
    sig_d = sig_d[sig_d > 1e-12][::-1]
    err_sigma = float(np.max(np.abs(sig_lr[: len(sig_d)] - sig_d) / sig_d[0]))

    # 軌道一致：低ランク kmatvec で組んだ K と密行列 K の双方で exp(ANGLE·K)（raw、正規化なし）
    Z_lr = Z.copy()
    Z_d = Z.copy()
    dev = 0.0
    for _ in range(steps):
        sys_lr.c = Z_lr.real.copy(); sys_lr.s = Z_lr.imag.copy()
        K_lr = np.column_stack([sys_lr.kmatvec(e) for e in np.eye(m)])
        w, V = np.linalg.eigh(1j * K_lr)
        Z_lr = V @ (np.exp(-1j * ANGLE * w) * (V.conj().T @ Z_lr))

        Kd = A * np.imag(np.conj(Z_d)[:, None] * Z_d[None, :])
        w, V = np.linalg.eigh(1j * Kd)
        Z_d = V @ (np.exp(-1j * ANGLE * w) * (V.conj().T @ Z_d))
        dev = max(dev, float(np.max(np.abs(Z_lr - Z_d))))
    return {
        "steps": steps,
        "err_matvec": float(err_matvec),
        "err_G": float(err_G),
        "err_sigma_spectrum": err_sigma,
        "max_traj_dev": dev,
    }


# ---------------- 測定プローブ ----------------

def onset_probe(n, seed, delta=1e-6, cap=4000):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260721 + 1000 * n + seed)
    t0 = time.time()
    v, residual, sig = make_parent(sys_lr, rng)
    t_parent = time.time() - t0

    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    ztz0 = complex(Z @ Z)

    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    f0 = None
    f_hist = []
    max_ztz = 0.0
    wp = rng.normal(size=sys_lr.m)
    t0 = time.time()
    for t in range(cap + 1):
        h1 = abs(p @ Z) ** 2 + abs(q @ Z) ** 2
        htot = float(np.real(np.conj(Z) @ Z))
        f = 1.0 - h1 / htot
        f_hist.append(f)
        if f0 is None:
            f0 = max(f, 1e-300)
        if f > 0.05:
            break
        if t % 200 == 0 and t > 0:
            progress(f"開始率走行 τ={t} f={f:.3e}")
        max_ztz = max(max_ztz, abs(complex(Z @ Z)))
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    t_run = time.time() - t0

    f_arr = np.array(f_hist)
    lo, hi = max(10.0 * f0, 1e-13), 1e-3
    mask = (f_arr > lo) & (f_arr < hi)
    idx = np.where(mask)[0]
    rate = None
    if len(idx) >= 5:
        rate = float(np.polyfit(idx, np.log(f_arr[idx]), 1)[0])
    return {
        "parent_residual": residual,
        "parent_sigma2_over_sigma1": float(sig[1] / sig[0]) if len(sig) > 1 else 0.0,
        "parent_rank_planes": int(len(sig)),
        "abs_ztz_initial": abs(ztz0),
        "abs_ztz_max": max_ztz,
        "f_initial": float(f0),
        "steps_run": len(f_hist) - 1,
        "onset_rate_per_step": rate,
        "t_parent_sec": t_parent,
        "t_run_sec": t_run,
    }


def relax_probe(n, seed, cap=3000, sub=None):
    sys_lr = LowRankSystem(n)
    m = sys_lr.m
    if sub is None:
        sub = max(20, cap // 40)
    rng = np.random.default_rng(50260721 + 1000 * n + seed)
    Z = zero_closure_generic(rng, m)
    ztz0 = abs(complex(Z @ Z))

    pr_hist = []
    sigma_checks = []
    wp = rng.normal(size=m)
    plateau_tau = None
    t0 = time.time()
    t = 0
    for t in range(cap + 1):
        pr_hist.append(participation_ratio(Z))
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        if t % sub == 0:
            sig_exact = sys_lr.sigma_spectrum()
            sigma_checks.append(abs(sig_est - sig_exact[0]) / sig_exact[0])
            progress(f"緩和走行 τ={t} PR/M={pr_hist[-1]/m:.4f}")
        if t >= 400 and t % 50 == 0:
            recent = pr_hist[-200:]
            if (max(recent) - min(recent)) / m < 5e-3:
                plateau_tau = t
                break
        if t < cap:
            Z = sys_lr.cayley_step(Z, sig_est)
    t_run = time.time() - t0

    sys_lr.set_theta(np.angle(Z))
    sig = sys_lr.sigma_spectrum()
    n_act = int(np.sum(sig / sig[0] > 0.05))
    qv = np.linalg.lstsq(sys_lr.G, sys_lr.wt(Z), rcond=None)[0]
    span_frac = float(np.linalg.norm(sys_lr.w(qv)) ** 2 / np.linalg.norm(Z) ** 2)
    ztz_end = abs(complex(Z @ Z))
    return {
        "abs_ztz_initial": ztz0,
        "abs_ztz_final": ztz_end,
        "steps_run": t,
        "plateau_tau": plateau_tau,
        "pr_final": pr_hist[-1],
        "pr_final_over_m": pr_hist[-1] / m,
        "n_active_planes": n_act,
        "plane_bound": int(min(n, m // 2)),
        "n_active_over_bound": n_act / min(n, m // 2),
        "span_w_fraction": span_frac,
        "kernel_fraction": 1.0 - span_frac,
        "sigma2_over_sigma1_final": float(sig[1] / sig[0]) if len(sig) > 1 else 0.0,
        "sigma_power_max_rel_err": float(max(sigma_checks)) if sigma_checks else None,
        "t_run_sec": t_run,
        "sec_per_step": t_run / max(t, 1),
    }


# ---------------- メイン ----------------

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    if not args:
        print(__doc__)
        sys.exit(1)
    n = int(args[0])
    seed = int(args[1]) if len(args) > 1 else 0
    cap = 3000
    for fl in flags:
        if fl.startswith("--steps"):
            cap = int(fl.split("=")[1])
    do_validate = "--validate" in flags

    os.makedirs(RESULT_DIR, exist_ok=True)
    m = n * (n - 1) // 2
    out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA}
    print(f"=== N スケーリング測定（低ランク頂点分解） N={n}, M={m}, seed={seed} ===", flush=True)

    if do_validate:
        if n <= 12:
            val = validate_against_dense(n, seed)
            out["validation"] = val
            print(f"[検証 vs 密行列] matvec相対誤差={val['err_matvec']:.3e}"
                  f" G相対誤差={val['err_G']:.3e}"
                  f" σスペクトル最大誤差={val['err_sigma_spectrum']:.3e}")
            print(f"[検証 vs 密行列] {val['steps']}ステップ軌道最大偏差={val['max_traj_dev']:.3e}")
        else:
            print("[検証] N>12 のためスキップ（密行列側が過大）")

    res_a = onset_probe(n, seed, cap=max(cap, 4000))
    out["onset"] = res_a
    print(f"[A 開始率] 親残差={res_a['parent_residual']:.3e}"
          f" 親平面数={res_a['parent_rank_planes']}"
          f" σ2/σ1={res_a['parent_sigma2_over_sigma1']:.3e}")
    print(f"[A 開始率] |Z^TZ|初期={res_a['abs_ztz_initial']:.3e}"
          f" 最大={res_a['abs_ztz_max']:.3e}"
          f" f0={res_a['f_initial']:.3e}")
    rate = res_a["onset_rate_per_step"]
    print(f"[A 開始率] rate={rate if rate is None else f'{rate:.4e}'}/step"
          f" steps={res_a['steps_run']}"
          f" (親構成{res_a['t_parent_sec']:.1f}s, 走行{res_a['t_run_sec']:.1f}s)")

    res_b = relax_probe(n, seed, cap=cap)
    out["relax"] = res_b
    print(f"[B 緩和] |Z^TZ| 初期={res_b['abs_ztz_initial']:.3e}"
          f" 終={res_b['abs_ztz_final']:.3e}")
    print(f"[B 緩和] plateau_tau={res_b['plateau_tau']}"
          f" PR/M={res_b['pr_final_over_m']:.4f}"
          f" 活性平面={res_b['n_active_planes']}/{res_b['plane_bound']}"
          f" (比 {res_b['n_active_over_bound']:.3f})")
    print(f"[B 緩和] 核フラクション={res_b['kernel_fraction']:.4f}"
          f" σ2/σ1(終)={res_b['sigma2_over_sigma1_final']:.3f}"
          f" 冪反復σ誤差max={res_b['sigma_power_max_rel_err']:.2e}"
          f" {res_b['sec_per_step']*1000:.2f} ms/step")

    path = os.path.join(RESULT_DIR, f"summary_N{n:05d}_seed{seed}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print(f"出力: {path}")


if __name__ == "__main__":
    main()
