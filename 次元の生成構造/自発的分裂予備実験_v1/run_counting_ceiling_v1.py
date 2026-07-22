#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
実験U系列: 条件付き計数上限の予備検証

主定理候補: 計数規約 C=(ω0, ε) を固定したとき
    n_read(C) ≤ min( ⌊Ω/ω0⌋, 1/ε², min(N,⌊M/2⌋) )

設計思想: 走行と計数を分離する。走行中は「重み付き振動数スペクトル」
  {(θ_j, h_j/h_total)}（瞬時生成子の各回転平面の回転角と配分比）だけを保存し、
  任意の規約での n_read は事後の再解析で得る。同じ走行データから床の選び方で
  異なる n_read が出ること自体が「数は読出しの属性」の直接デモになる。

n_read の操作的定義:
  振動数軸 (0, Ω]（Ω = 2 arctan γ）を分解能 ω0 = Ω/B でビン分けし、
  ビン内配分比の和が ε² 以上のビンを「読出せる波」と数える。

  U0: 親状態（一本波）→ 全規約で n_read=1 の較正
  U1: 位相のみ力学の拡大過程（N=8,12,16）で n_read(τ) がスロット天井 B・
      ランク天井 N で飽和し、決して超えないこと
  U2: 同一の最終スペクトルに対する (B, ε²) 格子の再解析ヒートマップ。
      硬い上界 n_read ≤ min(B, 1/ε², N) の全セル検証（違反数=0が判定）
      振幅天井 1/ε² は「最大到達値の稜線」として現れる（等分配では
      B > 1/ε² でビン配分が床を割り、計数は飽和でなく崩落する）
  U3: 充填か手前安定化か——熱化状態（充填）vs T2レジスタ凍結状態（手前）
      vs 親状態（=1）の充填率対比

出力: counting_ceiling_result_v1/ に JSON と図。
"""

import json
import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
import run_n_scaling_lowrank_v1 as lr
import run_register_modulation_v2 as t2

RESULT_DIR = os.path.join(BASE_DIR, "counting_ceiling_result_v1")
GAMMA = lr.GAMMA
OMEGA = 2.0 * math.atan(GAMMA)  # 帯域上端（本規約）


# ---------- スペクトル読出し ----------

def weighted_spectrum(sys_lr, Z):
    """瞬時生成子の重み付き振動数スペクトル {(θ_j, share_j)} と核配分。"""
    sys_lr.set_theta(np.angle(Z))
    ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
    idx = np.where(ev.imag > 1e-9)[0]
    sig = ev.imag[idx]
    smax = sig.max()
    h_total = float(np.real(np.conj(Z) @ Z))
    thetas, shares = [], []
    for k, i in enumerate(idx):
        v = sys_lr.w(EV[:, i].astype(complex))
        nv = np.linalg.norm(v)
        if nv < 1e-12:
            continue
        v = v / nv
        h = abs(np.vdot(v, Z)) ** 2 + abs(v @ Z) ** 2
        thetas.append(2.0 * math.atan(GAMMA * sig[k] / smax))
        shares.append(h / h_total)
    return np.array(thetas), np.array(shares)


def weighted_spectrum_dense(K, W):
    """密行列版（拡張系用）。K 反対称、W 複素状態。"""
    mu, V = np.linalg.eigh(1j * K)
    idx = np.where(mu > 1e-9)[0]
    sig = mu[idx]
    smax = sig.max()
    h_total = float(np.real(np.conj(W) @ W))
    thetas, shares = [], []
    for k, i in enumerate(idx):
        v = V[:, i]
        h = abs(np.vdot(v, W)) ** 2 + abs(v @ W) ** 2
        thetas.append(2.0 * math.atan(GAMMA * sig[k] / smax))
        shares.append(h / h_total)
    return np.array(thetas), np.array(shares)


def n_read(thetas, shares, n_bins, eps2):
    """規約 (ω0=Ω/B, ε) での読出し数。"""
    if len(thetas) == 0:
        return 0
    b = np.minimum((thetas / (OMEGA / n_bins)).astype(int), n_bins - 1)
    acc = np.zeros(n_bins)
    np.add.at(acc, b, shares)
    return int(np.sum(acc >= eps2))


# ---------- U0: 較正（親状態） ----------

def u0_calibration():
    out = {}
    for n in [8, 12]:
        sys_lr = lr.LowRankSystem(n)
        rng = np.random.default_rng(60260722)
        v, res, _ = lr.make_parent(sys_lr, rng)
        th, sh = weighted_spectrum(sys_lr, v)
        counts = {f"B={B},eps2={e}": n_read(th, sh, B, e)
                  for B in [2, 8, 32] for e in [1e-4, 1e-2, 0.2]}
        out[n] = {"parent_residual": res, "counts": counts}
    return out


# ---------- U1: 拡大過程の n_read(τ) ----------

def u1_expansion(n, steps=4000, sub=25, seed=0):
    sys_lr = lr.LowRankSystem(n)
    rng = np.random.default_rng(70260722 + seed)
    Z = lr.zero_closure_generic(rng, sys_lr.m)
    spectra = []
    taus = []
    for t in range(steps + 1):
        if t % sub == 0:
            th, sh = weighted_spectrum(sys_lr, Z)
            spectra.append((th, sh))
            taus.append(t)
        sys_lr.set_theta(np.angle(Z))
        sig = sys_lr.sigma_spectrum()[0]
        Z = sys_lr.cayley_step(Z, sig)
    return taus, spectra


# ---------- メイン ----------

def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    summary = {"omega_band": OMEGA}

    # U0
    print("=== U0: 較正（親状態＝一本波） ===", flush=True)
    u0 = u0_calibration()
    for n, d in u0.items():
        vals = set(d["counts"].values())
        print(f"  N={n}: 親残差={d['parent_residual']:.1e} "
              f"n_read 全規約={sorted(vals)}（期待: {{1}}）", flush=True)
    summary["u0"] = {str(k): v["counts"] for k, v in u0.items()}

    # U1
    print("=== U1: 拡大過程の飽和（位相のみ力学） ===", flush=True)
    conventions = [
        ("slot4", 4, 1e-4),     # スロット天井 B=4 が最小
        ("slot16", 16, 1e-4),   # スロット天井 B=16
        ("rank", 64, 1e-4),     # ランク天井 N が最小
    ]
    u1_data = {}
    violations = 0
    for n in [8, 12, 16]:
        taus, spectra = u1_expansion(n)
        curves = {}
        for name, B, e2 in conventions:
            c = [n_read(th, sh, B, e2) for th, sh in spectra]
            ceil = min(B, int(1.0 / e2), n)
            violations += sum(1 for x in c if x > ceil)
            curves[name] = (c, ceil)
            print(f"  N={n} [{name}] 天井={ceil} 終値={c[-1]} 最大={max(c)}"
                  f" {'OK' if max(c) <= ceil else '違反!'}", flush=True)
        u1_data[n] = (taus, spectra, curves)

    # U2: 最終スペクトルの (B, ε²) 格子再解析
    print("=== U2: 床格子の再解析（同一データ・走行なし） ===", flush=True)
    Bs = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    eps2s = np.logspace(-4, math.log10(0.5), 16)
    u2 = {}
    for n in [8, 16]:
        th, sh = u1_data[n][1][-1]
        grid = np.zeros((len(Bs), len(eps2s)), dtype=int)
        pred = np.zeros_like(grid)
        viol = 0
        for i, B in enumerate(Bs):
            for j, e2 in enumerate(eps2s):
                c = n_read(th, sh, B, e2)
                grid[i, j] = c
                pred[i, j] = min(B, int(1.0 / e2), n)
                if c > pred[i, j]:
                    viol += 1
        violations += viol
        ridge = grid.max(axis=0)  # 各 ε² での最大到達値（振幅天井の稜線）
        ridge_bound = np.minimum((1.0 / eps2s).astype(int), n)
        u2[n] = {"grid": grid, "pred": pred, "ridge": ridge,
                 "ridge_bound": ridge_bound, "violations": viol}
        print(f"  N={n}: 格子 {grid.shape} 全{grid.size}セル 違反={viol}"
              f" ｜稜線 max_B n_read ≤ min(1/ε², N): "
              f"{'OK' if np.all(ridge <= ridge_bound) else '違反!'}", flush=True)

    # U3: 充填 vs 手前安定化
    print("=== U3: 充填率の力学依存（規約 B=16, ε²=1e-3 固定） ===", flush=True)
    B3, e3 = 16, 1e-3
    u3 = {}
    for n in [8, 12]:
        ceil = min(B3, int(1.0 / e3), n)
        # (i) 熱化状態（U1終状態）
        th, sh = u1_data[n][1][-1]
        c_therm = n_read(th, sh, B3, e3)
        # (ii) T2 レジスタ凍結状態（β=1, 30000ステップ）
        A_edge, m = t2.build_edge_adjacency(n)
        rng = np.random.default_rng(95260722 + 1000 * n)
        W = t2.make_initial(rng, m)
        for _ in range(30000):
            g = t2.coupling_from_state(W, m, 1.0)
            W = t2.step(W, A_edge, m, g)
        theta = np.angle(W)
        K = np.zeros((m + 1, m + 1))
        K[:m, :m] = A_edge * np.sin(theta[None, :m] - theta[:m, None])
        b = g * np.sin(theta[m] - theta[:m])
        K[:m, m] = b
        K[m, :m] = -b
        th2, sh2 = weighted_spectrum_dense(K, W)
        c_frozen = n_read(th2, sh2, B3, e3)
        # (iii) 親状態
        sys_lr = lr.LowRankSystem(n)
        v, _, _ = lr.make_parent(sys_lr, np.random.default_rng(60260722))
        th3, sh3 = weighted_spectrum(sys_lr, v)
        c_parent = n_read(th3, sh3, B3, e3)
        u3[n] = {"ceiling": ceil, "thermal": c_therm,
                 "frozen": c_frozen, "parent": c_parent}
        print(f"  N={n} 天井={ceil}: 熱化={c_therm}（充填率{c_therm/ceil:.2f}） "
              f"凍結={c_frozen}（{c_frozen/ceil:.2f}） 親={c_parent}", flush=True)

    summary["u1"] = {str(n): {name: {"final": int(c[-1]), "max": int(max(c)),
                                     "ceiling": ceil}
                              for name, (c, ceil) in u1_data[n][2].items()}
                     for n in u1_data}
    summary["u2"] = {str(n): {"violations": int(d["violations"]),
                              "ridge": d["ridge"].tolist(),
                              "ridge_bound": d["ridge_bound"].tolist()}
                     for n, d in u2.items()}
    summary["u3"] = {str(n): d for n, d in u3.items()}
    summary["total_hard_bound_violations"] = int(violations)
    with open(os.path.join(RESULT_DIR, "summary_v1.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    # ---- 図1: U1 n_read(τ) と天井 ----
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
    for ax, n in zip(axes, [8, 12, 16]):
        taus = u1_data[n][0]
        for name, (c, ceil) in u1_data[n][2].items():
            line, = ax.plot(taus, c, lw=1.2, label=f"{name} (ceil {ceil})")
            ax.axhline(ceil, color=line.get_color(), ls=":", alpha=0.6)
        ax.set_xlabel(r"$\tau$")
        ax.set_title(rf"$N={n}$", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].set_ylabel(r"$n_{\mathrm{read}}(\tau;\mathcal{C})$")
    fig.suptitle("U1: readable-wave count under fixed conventions "
                 "(dotted = predicted ceiling, never exceeded)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "u1_saturation_v1.png"), dpi=160)
    plt.close(fig)

    # ---- 図2: U2 ヒートマップ（N=16） ----
    n = 16
    grid = u2[n]["grid"]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    im = axes[0].imshow(grid, aspect="auto", origin="lower", cmap="viridis")
    axes[0].set_xticks(range(0, len(eps2s), 3))
    axes[0].set_xticklabels([f"{e:.0e}" for e in eps2s[::3]], fontsize=7)
    axes[0].set_yticks(range(len(Bs)))
    axes[0].set_yticklabels(Bs, fontsize=7)
    axes[0].set_xlabel(r"$\varepsilon^2$")
    axes[0].set_ylabel(r"$B=\Omega/\omega_0$")
    axes[0].set_title(rf"$N={n}$: measured $n_{{\mathrm{{read}}}}$ "
                      "(same data, reanalysis only)", fontsize=10)
    fig.colorbar(im, ax=axes[0])
    axes[1].plot(eps2s, u2[n]["ridge"], "o-", label=r"$\max_B n_{\mathrm{read}}$")
    axes[1].plot(eps2s, u2[n]["ridge_bound"], "k--",
                 label=r"$\min(1/\varepsilon^2, N)$")
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"$\varepsilon^2$")
    axes[1].set_ylabel("count")
    axes[1].set_title("Amplitude-ceiling ridge", fontsize=10)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "u2_floor_grid_v1.png"), dpi=160)
    plt.close(fig)

    # ---- 図3: U3 充填率 ----
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    xs = np.arange(2)
    width = 0.25
    for k, kind in enumerate(["thermal", "frozen", "parent"]):
        vals = [u3[n][kind] / u3[n]["ceiling"] for n in [8, 12]]
        ax.bar(xs + (k - 1) * width, vals, width,
               label={"thermal": "phase-only (thermalized)",
                      "frozen": "register-feedback (frozen)",
                      "parent": "parent (single wave)"}[kind])
    ax.set_xticks(xs)
    ax.set_xticklabels([r"$N=8$", r"$N=12$"])
    ax.set_ylabel(r"filling ratio $n_{\mathrm{read}}/n_{\max}(\mathcal{C})$")
    ax.set_title("U3: ceiling is set by convention; filling is set by dynamics",
                 fontsize=10)
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "u3_filling_v1.png"), dpi=160)
    plt.close(fig)

    print(f"硬い上界の違反総数 = {violations}（0が判定基準）")
    print(f"出力: {RESULT_DIR}")


if __name__ == "__main__":
    main()
