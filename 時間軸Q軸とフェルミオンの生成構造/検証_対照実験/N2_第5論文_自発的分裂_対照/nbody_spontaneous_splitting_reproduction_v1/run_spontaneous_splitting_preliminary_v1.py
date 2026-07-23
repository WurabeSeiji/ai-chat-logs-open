#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
自発的分裂 予備実験 v1（実験N統制 + 実験O1）

目的:
  固定生成子（実験N統制）と逐次再構成正弦生成子（実験O1）のもとで、
  参照平面レジスタの移動（レジスタ再配分としての分裂開始）が
  自発的に生じるかを予備測定する。

設計要点:
  - 状態は複素ベクトル Z ∈ C^M（M = C(N,2) 関係波）。
  - 更新は実直交 Cayley 行列 U = (I - γK)^{-1}(I + γK)。
    実直交更新は複素双線形二乗和 Z^T Z（閉鎖量、§13の実部・虚部二条件）と
    エルミートノルム |Z|^2 を同時に厳密保存する。
  - 統制（実験N）: K を τ=0 の位相 arg(Z0) で評価して凍結。
    参照平面レジスタは機械精度で個別保存されるはず（帰結16.5/16.6）。
  - 実験O1: K を毎ステップ現在の位相 arg(Z(τ)) から再構成（凍結解除）。
    生成子は位相差のみの関数なので零次同次（公理0.5 スケール無名性に適合）。
  - 親状態: K(arg Z0) の最大σ平面に集中、他平面は相対種振幅 δ。
    自己無撞着不動点はゲージ固定（平面内基底を参照軸へピン留め）付きの
    減衰反復で追い込み、残差と達成集中度を記録する。
  - 自発性判定の予備データ: δ 系列で休眠フラクション f(τ) の成長を測定し、
    f(0) 比の成長時刻の δ 依存性を調べる。
  - スケール不変性検証（公理0.5）: λ=1024（2冪、まるめ誤差ゼロ）で
    比率観測量が bitwise 一致することを確認。非2冪 λ・共通位相回転・
    微小双子摂動の三者の軌道発散率を比較し、発散がスケール対称性の破れでなく
    カオス的鋭敏性によることを切り分ける。

出力: spontaneous_splitting_result_v1/ に JSON・CSV・図。
"""

import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "spontaneous_splitting_result_v1")

GAMMA = math.tan(math.pi / 144.0)
N_BODIES = 4
STEPS = 4320
DELTAS = [1e-2, 1e-3, 1e-4, 1e-5]
SEEDS = list(range(6))
SC_ITERS = 400
SC_BETA = 0.5
ACTIVE_EPS2 = 0.0025  # 相対レジスタ h_j/h_total の活性判定閾値 (= 0.05^2)
F_ABS_LEVELS = [1e-2, 1e-1]  # 休眠フラクションの絶対閾値
F_REL_LEVELS = [10.0, 100.0]  # f(0) 比の相対閾値
SUBSAMPLE = 8


def complete_graph_edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def line_graph_adjacency(n):
    edges = complete_graph_edges(n)
    m = len(edges)
    A = np.zeros((m, m))
    for a in range(m):
        for b in range(m):
            if a != b and len(set(edges[a]) & set(edges[b])) > 0:
                A[a, b] = 1.0
    return A


def sine_generator(theta, A):
    """K_ef = A_ef sin(θ_f - θ_e)、スペクトルノルム正規化。位相差のみの関数（零次同次）。"""
    K = A * np.sin(theta[None, :] - theta[:, None])
    s = np.linalg.norm(K, 2)
    if s < 1e-300:
        return K
    return K / s


def cayley(K, gamma=GAMMA):
    I = np.eye(K.shape[0])
    return np.linalg.solve(I - gamma * K, I + gamma * K)


def plane_decomposition(K, tol=1e-10):
    """実反対称 K の回転平面 (σ_j, p_j, q_j) と核の実正規直交基底を抽出。

    向き規約 K p = +σ q。平面内 SO(2) ゲージは、平面への射影が最大となる
    標準基底ベクトルへ p をピン留めして固定する（決定的にするため）。
    """
    mu, V = np.linalg.eigh(1j * K)
    order = np.argsort(-mu)
    planes = []
    for idx in order:
        m_val = float(mu[idx])
        if m_val <= tol:
            continue
        v = V[:, idx]
        p = np.sqrt(2.0) * v.real
        q = np.sqrt(2.0) * v.imag
        p = p / np.linalg.norm(p)
        q = q - (q @ p) * p
        q = q / np.linalg.norm(q)
        # ゲージ固定: 平面射影が最大の標準基底軸に p を揃える
        proj = p[None, :] * p[:, None] + q[None, :] * q[:, None]  # P = p p^T + q q^T
        norms = np.linalg.norm(proj, axis=1)
        ref = int(np.argmax(norms))
        p_new = proj[ref] / norms[ref]
        q_new = (K @ p_new) / m_val
        q_new = q_new / np.linalg.norm(q_new)
        planes.append({"sigma": m_val, "p": p_new, "q": q_new})
    # 核: SVD の零特異値方向（実正規直交）
    _, s, Vt = np.linalg.svd(K)
    kernel = [Vt[i] for i in range(K.shape[0]) if s[i] <= tol]
    return planes, kernel


def plane_count(K, rel_tol=0.05):
    """瞬時生成子の活性回転平面数。成分数増加の直接読出し。

    正規化後は σ_max = 1 なので、σ > rel_tol は比率閾値である
    （公理0.5: 絶対閾値は置けない。数値ダスト σ = O(δ) を平面と数えない）。"""
    mu = np.linalg.eigvalsh(1j * K)
    return int(np.sum(mu > rel_tol))


def registers(Z, planes):
    """参照平面ごとのエルミートレジスタ h_j と双線形（閉鎖）レジスタ c_j。"""
    h = np.empty(len(planes))
    c = np.empty(len(planes), dtype=complex)
    for j, pl in enumerate(planes):
        a = pl["p"] @ Z
        b = pl["q"] @ Z
        h[j] = abs(a) ** 2 + abs(b) ** 2
        c[j] = a * a + b * b
    return h, c


def phase_err(t1, t2):
    """共通位相をmod outした位相差の最大絶対値（位相差収束の真の残差）。"""
    d = np.angle(np.exp(1j * (t1 - t2)))
    common = np.angle(np.sum(np.exp(1j * d)))
    return float(np.max(np.abs(np.angle(np.exp(1j * (d - common))))))


def prepare_initial_state(rng, A, m, delta, iters=SC_ITERS, beta=SC_BETA):
    """親状態の自己無撞着反復。

    親状態は K(arg Z) の最大σ平面の固有モード v = (p + iq)/√2（円偏波）。
    実係数の平面内状態は全成分位相が共線（ψ, ψ+π）となり sin 位相差が消えて
    生成子が退化するため、自己無撞着な親は固有モードでなければならない。
    固有モードなら更新は共通位相回転（U v = e^{-iθ} v）なので位相差が不変、
    したがって K(arg Z) が不変となり、δ=0 で厳密な周期軌道になる。
    不動点反復: θ ← arg(v(K(θ)))。平面内 SO(2) ゲージは共通位相にしか
    効かないため、この反復はゲージ選択に依存しない。

    数値的知見: 自己無撞着な K(arg v) は rank 2（回転平面が一つだけ）で、
    残り m-2 次元は核（静的部分空間）となる。したがって休眠種 δ は
    核方向へ複素係数で置く。
    """
    theta = rng.uniform(0.0, 2.0 * np.pi, m)
    g1 = rng.normal(size=m)  # 休眠種の生データ（実、決定的）
    g2 = rng.normal(size=m)
    residual = None
    v = None
    for _ in range(iters):
        K = sine_generator(theta, A)
        planes, _ = plane_decomposition(K)
        pl0 = planes[0]
        v = (pl0["p"] + 1j * pl0["q"]) / math.sqrt(2.0)
        theta_new = np.angle(v)
        residual = phase_err(theta_new, theta)
        mix = (1.0 - beta) * np.exp(1j * theta) + beta * np.exp(1j * theta_new)
        theta = np.angle(mix)
    K_fin = sine_generator(np.angle(v), A)
    # 真の自己無撞着残差: v が K(arg v) の固有モードであることの残差
    mu_v = float(np.real(np.conj(v) @ (1j * K_fin @ v)))
    residual = float(np.linalg.norm(1j * K_fin @ v - mu_v * v))
    _, kernel = plane_decomposition(K_fin)
    # 零閉鎖種: 核内の実正規直交対 (u,w) から g=(u+iw)/√2（g^T g = 0 厳密）
    def kproj(x):
        out = np.zeros(m)
        for k_vec in kernel:
            out += (k_vec @ x) * k_vec
        return out
    u = kproj(g1)
    w = kproj(g2)
    if delta > 0.0 and np.linalg.norm(u) > 0.0:
        u = u / np.linalg.norm(u)
        w = w - (w @ u) * u
        w = w / np.linalg.norm(w)
        g = (u + 1j * w) / math.sqrt(2.0)
        Z = v + delta * g
    else:
        Z = v.copy()
    Z = Z / np.linalg.norm(Z)
    return Z, residual


def run_dynamics(Z0, A, planes_ref, steps, sequential):
    """力学を走らせ、参照レジスタ等の時系列を返す。"""
    Z = Z0.copy()
    K_ref = sine_generator(np.angle(Z0), A)
    U_fixed = cayley(K_ref)
    n_pl_ref = plane_count(K_ref)
    n_pl = len(planes_ref)
    h_hist = np.empty((steps + 1, n_pl))
    c_total_hist = np.empty(steps + 1, dtype=complex)
    h_total_hist = np.empty(steps + 1)
    k_dist_hist = np.empty(steps + 1)
    n_planes_hist = np.empty(steps + 1, dtype=int)
    sigma2_hist = np.empty(steps + 1)
    mu_ref = np.linalg.eigvalsh(1j * K_ref)
    sigma2_ref = float(np.sort(mu_ref)[-2]) if len(mu_ref) >= 2 else 0.0
    for t in range(steps + 1):
        h, _ = registers(Z, planes_ref)
        h_hist[t] = h
        c_total_hist[t] = Z @ Z
        h_total_hist[t] = float(np.real(np.conj(Z) @ Z))
        if sequential:
            K_now = sine_generator(np.angle(Z), A)
            k_dist_hist[t] = float(np.linalg.norm(K_now - K_ref, "fro"))
            mu_now = np.linalg.eigvalsh(1j * K_now)
            n_planes_hist[t] = int(np.sum(mu_now > 0.05))
            sigma2_hist[t] = float(np.sort(mu_now)[-2])
            U = cayley(K_now)
        else:
            k_dist_hist[t] = 0.0
            n_planes_hist[t] = n_pl_ref
            sigma2_hist[t] = sigma2_ref
            U = U_fixed
        if t < steps:
            Z = U @ Z
    return {
        "h": h_hist,
        "c_total": c_total_hist,
        "h_total": h_total_hist,
        "k_dist": k_dist_hist,
        "n_planes": n_planes_hist,
        "sigma2": sigma2_hist,
    }


def dormant_fraction(rec):
    h_total = rec["h_total"]
    return 1.0 - rec["h"][:, 0] / h_total


def conservation_errors(rec):
    c = rec["c_total"]
    h = rec["h_total"]
    return {
        "max_abs_dev_re_c": float(np.max(np.abs(c.real - c.real[0]))),
        "max_abs_dev_im_c": float(np.max(np.abs(c.imag - c.imag[0]))),
        "max_rel_dev_h_total": float(np.max(np.abs(h - h[0])) / h[0]),
    }


def crossing_times(f, abs_levels, rel_levels):
    out = {}
    f0 = max(f[0], 1e-300)
    for lv in abs_levels:
        key = f"abs_{lv:g}"
        if f0 >= lv / 10.0:
            out[key] = None
            continue
        idx = int(np.argmax(f >= lv))
        out[key] = idx if f[idx] >= lv else None
    for r in rel_levels:
        key = f"rel_{r:g}x"
        lv = r * f0
        if lv >= 1.0:
            out[key] = None
            continue
        idx = int(np.argmax(f >= lv))
        out[key] = idx if f[idx] >= lv else None
    return out


def early_growth_rate(f):
    """f が min(100 f(0), 1e-2) に達するまでの窓で ln f の線形フィット傾き。"""
    f0 = max(f[0], 1e-300)
    cap = min(100.0 * f0, 1e-2)
    idx = int(np.argmax(f >= cap))
    end = idx if f[idx] >= cap else len(f)
    if end < 8:
        return None, int(end)
    tau = np.arange(end)
    y = np.log(np.maximum(f[:end], 1e-300))
    slope = float(np.polyfit(tau, y, 1)[0])
    return slope, int(end)


def ratio_trajectory(rec):
    return rec["h"] / rec["h_total"][:, None]


def divergence_curve(ratio_a, ratio_b):
    return np.max(np.abs(ratio_a - ratio_b), axis=1)


def divergence_rate(div, floor=1e-16, cap=1e-3):
    """軌道間距離の指数増大率（カオス鋭敏性の指標）。"""
    mask = (div > floor) & (div < cap)
    idx = np.where(mask)[0]
    if len(idx) < 8:
        return None
    y = np.log(div[idx])
    slope = float(np.polyfit(idx, y, 1)[0])
    return slope


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    A = line_graph_adjacency(N_BODIES)
    m = A.shape[0]

    summary = {
        "params": {
            "n_bodies": N_BODIES,
            "m_relation_waves": m,
            "gamma": GAMMA,
            "steps": STEPS,
            "deltas": DELTAS,
            "seeds": SEEDS,
            "self_consistency_iters": SC_ITERS,
            "self_consistency_beta": SC_BETA,
            "active_eps2": ACTIVE_EPS2,
            "f_abs_levels": F_ABS_LEVELS,
            "f_rel_levels": F_REL_LEVELS,
        },
        "runs": [],
        "scale_invariance": {},
        "sensitivity": {},
    }

    f_curves = {}
    rep_records = {}

    for delta in DELTAS:
        f_curves[delta] = []
        for seed in SEEDS:
            rng = np.random.default_rng(20260721 + seed)
            Z0, sc_residual = prepare_initial_state(rng, A, m, delta)
            K_ref = sine_generator(np.angle(Z0), A)
            planes_ref, kernel_ref = plane_decomposition(K_ref)
            n_pl = len(planes_ref)
            h0, _ = registers(Z0, planes_ref)
            h_tot0 = float(np.real(np.conj(Z0) @ Z0))
            init_dorm = 1.0 - h0[0] / h_tot0

            rec_fix = run_dynamics(Z0, A, planes_ref, STEPS, sequential=False)
            rec_seq = run_dynamics(Z0, A, planes_ref, STEPS, sequential=True)

            f_fix = dormant_fraction(rec_fix)
            f_seq = dormant_fraction(rec_seq)
            f_curves[delta].append(f_seq)

            h_fix = ratio_trajectory(rec_fix)
            fix_reg_dev = float(np.max(np.abs(h_fix - h_fix[0])))

            slope, fit_end = early_growth_rate(f_seq)

            run_entry = {
                "delta": delta,
                "seed": seed,
                "n_planes": n_pl,
                "sigmas": [pl["sigma"] for pl in planes_ref],
                "self_consistency_residual": sc_residual,
                "abs_ztz_initial": float(abs(complex(Z0 @ Z0))),
                "initial_dormant_fraction": float(init_dorm),
                "concentration_quality": float(init_dorm / (delta * delta)),
                "fixed": {
                    "max_register_ratio_dev": fix_reg_dev,
                    **conservation_errors(rec_fix),
                    "final_dormant_fraction": float(f_fix[-1]),
                },
                "sequential": {
                    **conservation_errors(rec_seq),
                    "final_dormant_fraction": float(f_seq[-1]),
                    "max_dormant_fraction": float(np.max(f_seq)),
                    "crossing_times": crossing_times(f_seq, F_ABS_LEVELS, F_REL_LEVELS),
                    "early_growth_rate_per_step": slope,
                    "growth_fit_window": fit_end,
                    "n_planes_final": int(rec_seq["n_planes"][-1]),
                    "n_planes_max": int(np.max(rec_seq["n_planes"])),
                    "sigma2_max": float(np.max(rec_seq["sigma2"])),
                    "max_k_dist_from_ref": float(np.max(rec_seq["k_dist"])),
                },
            }
            summary["runs"].append(run_entry)

            if seed == 0:
                rep_records[delta] = {"fix": rec_fix, "seq": rec_seq, "Z0": Z0}

    # ---- δ=0 検証: 固有モード親状態が自己無撞着周期軌道であること ----
    rng0 = np.random.default_rng(20260721)
    Z_orbit, orbit_residual = prepare_initial_state(rng0, A, m, 0.0)
    K_ref0 = sine_generator(np.angle(Z_orbit), A)
    planes_ref0, _ = plane_decomposition(K_ref0)
    rec_orbit = run_dynamics(Z_orbit, A, planes_ref0, STEPS, sequential=True)
    f_orbit = dormant_fraction(rec_orbit)
    period = int(round(math.pi / math.atan(GAMMA * planes_ref0[0]["sigma"])))
    idx_visible = np.where(f_orbit > 1e-10)[0]
    summary["zero_delta_orbit"] = {
        "self_consistency_residual": orbit_residual,
        "f_initial": float(f_orbit[0]),
        "f_after_one_period": float(f_orbit[min(period, STEPS)]),
        "max_k_dist_first_period": float(np.max(rec_orbit["k_dist"][: period + 1])),
        "f_final": float(f_orbit[-1]),
        "tau_f_exceeds_1e-10": int(idx_visible[0]) if len(idx_visible) else None,
        "period_steps": period,
    }

    # ---- スケール不変性・鋭敏性の切り分け（代表条件 δ=1e-3, seed=0） ----
    delta_rep = 1e-3
    Z0 = rep_records[delta_rep]["Z0"]
    K_ref = sine_generator(np.angle(Z0), A)
    planes_ref, _ = plane_decomposition(K_ref)
    base = rep_records[delta_rep]["seq"]
    base_ratio = ratio_trajectory(base)

    # (a) λ = 1024 = 2^10: 浮動小数まるめ誤差ゼロのスケール変換 → bitwise 一致が期待値
    rec = run_dynamics(1024.0 * Z0, A, planes_ref, STEPS, sequential=True)
    summary["scale_invariance"]["lambda_1024_max_ratio_dev"] = float(
        np.max(np.abs(ratio_trajectory(rec) - base_ratio))
    )

    # (b) 非2冪 λ / 共通位相回転 / 微小双子摂動: 軌道発散率を比較し、
    #     発散がカオス鋭敏性（初期まるめ差の指数増幅）であることを確認
    probes = {}
    rec_l = run_dynamics(987.654 * Z0, A, planes_ref, STEPS, sequential=True)
    probes["lambda_987p654"] = divergence_curve(ratio_trajectory(rec_l), base_ratio)
    lam_ph = complex(math.cos(0.7), math.sin(0.7))
    rec_p = run_dynamics(lam_ph * Z0, A, planes_ref, STEPS, sequential=True)
    probes["phase_0p7"] = divergence_curve(ratio_trajectory(rec_p), base_ratio)
    Z_twin = Z0.copy()
    Z_twin[0] *= 1.0 + 1e-13
    rec_t = run_dynamics(Z_twin, A, planes_ref, STEPS, sequential=True)
    probes["twin_1e-13"] = divergence_curve(ratio_trajectory(rec_t), base_ratio)

    for name, div in probes.items():
        summary["sensitivity"][name] = {
            "final_divergence": float(div[-1]),
            "max_divergence": float(np.max(div)),
            "divergence_rate_per_step": divergence_rate(div),
        }

    # 統制側の同じ摂動（固定生成子はカオスでないことの対照）
    rec_t_fix = run_dynamics(Z_twin, A, planes_ref, STEPS, sequential=False)
    base_fix_ratio = ratio_trajectory(rep_records[delta_rep]["fix"])
    div_fix = divergence_curve(ratio_trajectory(rec_t_fix), base_fix_ratio)
    summary["sensitivity"]["twin_1e-13_fixed_control"] = {
        "final_divergence": float(div_fix[-1]),
        "max_divergence": float(np.max(div_fix)),
        "divergence_rate_per_step": divergence_rate(div_fix),
    }

    tau = np.arange(STEPS + 1)

    # ---- 図1: 統制 vs 逐次再構成の参照レジスタ軌道（δ=1e-3, seed=0） ----
    rec_fix = rep_records[delta_rep]["fix"]
    rec_seq = rep_records[delta_rep]["seq"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    for ax, rec, title in [
        (axes[0], rec_fix, "Fixed generator (control, Exp N)"),
        (axes[1], rec_seq, "Sequentially reconstructed generator (Exp O1)"),
    ]:
        ratio = ratio_trajectory(rec)
        ax.semilogy(tau, np.maximum(ratio[:, 0], 1e-18), label="parent plane")
        ax.semilogy(
            tau, np.maximum(1.0 - ratio[:, 0], 1e-18), label="kernel (dormant) fraction"
        )
        ax.set_xlabel(r"$\tau$")
        ax.set_title(title, fontsize=10)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"$h_j / h_{\mathrm{total}}$")
    axes[0].legend(fontsize=8)
    fig.suptitle(
        r"Reference-plane register trajectories ($N=4$, $\delta=10^{-3}$)", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(
        os.path.join(RESULT_DIR, "register_migration_fixed_vs_sequential_v1.png"),
        dpi=160,
    )
    plt.close(fig)

    # ---- 図2: 休眠フラクション f(τ) の δ 系列（seed 中央値） ----
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for delta in DELTAS:
        curves = np.array(f_curves[delta])
        med = np.median(curves, axis=0)
        ax.semilogy(tau, np.maximum(med, 1e-18), label=rf"$\delta={delta:g}$")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"dormant fraction $f(\tau) = 1 - h_1/h_{\mathrm{total}}$")
    ax.set_title(
        "Dormant-register growth under sequential reconstruction (median over seeds)"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(
        os.path.join(RESULT_DIR, "dormant_fraction_growth_delta_series_v1.png"), dpi=160
    )
    plt.close(fig)

    # ---- 図3: 閾値到達時刻 vs log10(1/δ) ----
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    style = {
        f"abs_{F_ABS_LEVELS[0]:g}": ("o", rf"$f \geq {F_ABS_LEVELS[0]:g}$"),
        f"abs_{F_ABS_LEVELS[1]:g}": ("s", rf"$f \geq {F_ABS_LEVELS[1]:g}$"),
        f"rel_{F_REL_LEVELS[0]:g}x": ("^", r"$f \geq 10 f(0)$"),
        f"rel_{F_REL_LEVELS[1]:g}x": ("v", r"$f \geq 100 f(0)$"),
    }
    for key, (mk, lab) in style.items():
        xs, ys = [], []
        for run in summary["runs"]:
            tc = run["sequential"]["crossing_times"].get(key)
            if tc is not None:
                xs.append(math.log10(1.0 / run["delta"]))
                ys.append(tc)
        if xs:
            ax.plot(xs, ys, mk, alpha=0.6, label=lab)
    ax.set_xlabel(r"$\log_{10}(1/\delta)$")
    ax.set_ylabel(r"crossing step $\tau_c$")
    ax.set_title("Threshold-crossing time vs seed amplitude")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "crossing_time_vs_delta_v1.png"), dpi=160)
    plt.close(fig)

    # ---- 図4: 保存量誤差と鋭敏性 ----
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axes[0]
    for rec, ls, name in [(rec_fix, "--", "fixed"), (rec_seq, "-", "sequential")]:
        c = rec["c_total"]
        ax.semilogy(
            tau,
            np.maximum(np.abs(c.real - c.real[0]), 1e-20),
            ls,
            label=f"{name}: |ΔRe(ZᵀZ)|",
        )
        ax.semilogy(
            tau,
            np.maximum(np.abs(c.imag - c.imag[0]), 1e-20),
            ls,
            alpha=0.6,
            label=f"{name}: |ΔIm(ZᵀZ)|",
        )
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("absolute deviation")
    ax.set_title("Exact conservation of the complex square-sum closure", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    for name, div in probes.items():
        ax.semilogy(tau, np.maximum(div, 1e-18), label=name)
    ax.semilogy(tau, np.maximum(div_fix, 1e-18), "--", label="twin 1e-13 (fixed ctrl)")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"trajectory divergence $\max_j |\Delta(h_j/h_{\mathrm{total}})|$")
    ax.set_title("Sensitivity probes (chaos vs scale-symmetry breaking)", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        os.path.join(RESULT_DIR, "closure_conservation_and_sensitivity_v1.png"), dpi=160
    )
    plt.close(fig)

    # ---- 図5: 瞬時生成子の第2σと平面数（分裂の固定/非固定の読出し） ----
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for delta in DELTAS:
        rec = rep_records[delta]["seq"]
        ax.semilogy(tau, np.maximum(rec["sigma2"], 1e-18), alpha=0.85,
                    label=rf"$\delta={delta:g}$ (seq)")
    ax.axhline(0.05, color="k", ls=":", alpha=0.7,
               label="active-plane threshold (0.05, relative)")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"second singular direction $\sigma_2$ of instantaneous $K$")
    ax.set_title("No stable second rotation plane forms "
                 r"($\sigma_2 \ll \sigma_1 = 1$ throughout)", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "instantaneous_plane_count_v1.png"), dpi=160)
    plt.close(fig)

    # ---- CSV: f(τ) 中央値（サブサンプル） ----
    import csv as csvmod

    with open(
        os.path.join(RESULT_DIR, "dormant_fraction_curves_v1.csv"), "w", newline=""
    ) as fh:
        w = csvmod.writer(fh)
        w.writerow(["tau"] + [f"f_median_delta_{d:g}" for d in DELTAS])
        med_all = {d: np.median(np.array(f_curves[d]), axis=0) for d in DELTAS}
        for t in range(0, STEPS + 1, SUBSAMPLE):
            w.writerow([t] + [f"{med_all[d][t]:.12e}" for d in DELTAS])

    with open(os.path.join(RESULT_DIR, "summary_v1.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    # ---- コンソール要約 ----
    print("=== 自発的分裂 予備実験 v1 要約 ===")
    print(f"N={N_BODIES}, M={m}, steps={STEPS}, gamma={GAMMA:.6e}")
    scr = [r["self_consistency_residual"] for r in summary["runs"]]
    print(f"自己無撞着残差: median={np.median(scr):.3e} max={max(scr):.3e}")
    cq = [r["concentration_quality"] for r in summary["runs"]]
    print(f"集中度 init_dorm/δ²: median={np.median(cq):.3f} max={max(cq):.3f}")
    fix_devs = [r["fixed"]["max_register_ratio_dev"] for r in summary["runs"]]
    print(f"[統制/実験N] 参照レジスタ比率の最大変動: max={max(fix_devs):.3e}")
    cons_keys = ["max_abs_dev_re_c", "max_abs_dev_im_c", "max_rel_dev_h_total"]
    for mode in ["fixed", "sequential"]:
        vals = {k: max(r[mode][k] for r in summary["runs"]) for k in cons_keys}
        print(f"[{mode}] 保存量最大誤差: {vals}")
    for delta in DELTAS:
        runs = [r for r in summary["runs"] if r["delta"] == delta]
        finals = [r["sequential"]["final_dormant_fraction"] for r in runs]
        maxima = [r["sequential"]["max_dormant_fraction"] for r in runs]
        rates = [
            r["sequential"]["early_growth_rate_per_step"]
            for r in runs
            if r["sequential"]["early_growth_rate_per_step"] is not None
        ]
        tc_abs = [
            r["sequential"]["crossing_times"][f"abs_{F_ABS_LEVELS[0]:g}"]
            for r in runs
            if r["sequential"]["crossing_times"][f"abs_{F_ABS_LEVELS[0]:g}"] is not None
        ]
        tc_rel = [
            r["sequential"]["crossing_times"][f"rel_{F_REL_LEVELS[0]:g}x"]
            for r in runs
            if r["sequential"]["crossing_times"][f"rel_{F_REL_LEVELS[0]:g}x"] is not None
        ]
        init_d = np.median([r["initial_dormant_fraction"] for r in runs])
        print(
            f"[δ={delta:g}] init_dorm~{init_d:.2e}"
            f" f_final(med)={np.median(finals):.3e} f_max(med)={np.median(maxima):.3e}"
            f" rate(med)={(np.median(rates) if rates else float('nan')):.3e}/step"
            f" tau_c(f>1e-2)(med)={np.median(tc_abs) if tc_abs else 'n/a'}"
            f" tau_c(10x)(med)={np.median(tc_rel) if tc_rel else 'n/a'}"
            f" n_planes(final/max)="
            f"{np.median([r['sequential']['n_planes_final'] for r in runs]):.0f}/"
            f"{np.median([r['sequential']['n_planes_max'] for r in runs]):.0f}"
        )
    zo = summary["zero_delta_orbit"]
    print(
        f"[δ=0 周期軌道検証] 残差={zo['self_consistency_residual']:.3e}"
        f" f0={zo['f_initial']:.3e} f(1周期)={zo['f_after_one_period']:.3e}"
        f" maxKdist(1周期)={zo['max_k_dist_first_period']:.3e}"
        f" f_final={zo['f_final']:.3e} τ(f>1e-10)={zo['tau_f_exceeds_1e-10']}"
    )
    print(f"スケール厳密性 λ=1024 (bitwise期待): "
          f"max dev = {summary['scale_invariance']['lambda_1024_max_ratio_dev']:.3e}")
    for name, entry in summary["sensitivity"].items():
        rate = entry["divergence_rate_per_step"]
        print(
            f"[鋭敏性 {name}] final={entry['final_divergence']:.3e}"
            f" max={entry['max_divergence']:.3e}"
            f" rate={(rate if rate is not None else float('nan')):.3e}/step"
        )
    print(f"結果出力先: {RESULT_DIR}")


if __name__ == "__main__":
    main()
