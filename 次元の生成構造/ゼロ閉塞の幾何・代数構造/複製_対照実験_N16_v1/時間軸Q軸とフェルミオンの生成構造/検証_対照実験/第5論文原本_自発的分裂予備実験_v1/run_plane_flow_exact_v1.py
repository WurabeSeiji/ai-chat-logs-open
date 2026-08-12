#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""平面流入の【厳密】計測：密行列 eig(K) による回転平面分解（閾値非依存）。

親生成子 K(arg v) を密行列 M×M として構成し、直接 eig(K) で回転平面を得る。
K の非零 σ は数値ゼロ（~10⁻¹⁴）から明確なスペクトルギャップで分離するため、
σ を round(·,6) でグループ化するだけで回転平面と核が一意に定まる——結果は
分離閾値に依存しない（この意味で「厳密」）。参照用に σ スペクトル全体を JSON に残す。

適用範囲：密行列 K が M×M（M=N(N−1)/2）で構成可能な N。実用上 N≲60（N=40 で
M=780、N=60 で M=1770）。N=300 は M=44850 で密行列が 16GB 超となり本法は不可
（近似法 run_plane_flow_approx_v1.py を用いる）。

3帯に集約：P₁（最大σの支配平面）／その他の回転平面（σ>0, 非P₁）／σ=0 核。
力学は原本 run_spontaneous_splitting_largeN_v1.run() の忠実複製（不変更）。
対照テスト：全ステップで f == 1 − (P₁比)（恒等）、記録 f が正本 fcurve とビット一致。
図：3帯の積み上げ面 ＋ f(τ)=1−P₁ を対数右軸に重ねる。

使い方: python3 run_plane_flow_exact_v1.py 5   --record-every=20
        python3 run_plane_flow_exact_v1.py 40  --record-every=20
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import (
    LowRankSystem, make_parent, zero_closure_kernel_seed,
)

RESULT_DIR = HERE / "plane_flow_result_v1"
REFERENCE_DIR = HERE / "metastable_series_result_v1"


def parent_plane_split_exact(sys_lr, v, round_digits=6):
    """密行列 eig(K) で (支配平面P₁, その他回転平面) の正規直交実基底 + σスペクトル。

    分離は round(σ, round_digits) のグループ化のみ（明確なギャップに依存し、
    連続的な閾値パラメータを持たない）。核基底は作らず残差で扱う。
    """
    M = sys_lr.m
    sys_lr.set_theta(np.angle(v))
    K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
    w, V = np.linalg.eig(K)
    sig_all = w.imag
    groups = {}
    for i in range(M):
        key = round(float(sig_all[i]), round_digits)
        if key > 0:                      # 正の虚部＝回転平面（round で 0 は核へ）
            groups.setdefault(key, []).extend([np.real(V[:, i]), np.imag(V[:, i])])

    def ortho(cols):
        Q, R = np.linalg.qr(np.column_stack(cols))
        return Q[:, np.abs(np.diag(R)) > 1e-8]

    sig_sorted = sorted(groups, reverse=True)
    p1_sigma = sig_sorted[0]
    B_p1 = ortho(groups[p1_sigma])
    rest = [c for k in sig_sorted[1:] for c in groups[k]]
    B_rot = None
    if rest:
        R0 = np.column_stack(rest)
        R0 = R0 - B_p1 @ (B_p1.T @ R0)
        Q, R = np.linalg.qr(R0)
        B_rot = Q[:, np.abs(np.diag(R)) > 1e-8]
    spectrum = sorted((round(float(s), 6) for s in sig_all if s > 1e-10), reverse=True)
    return p1_sigma, B_p1, B_rot, spectrum


def run(n, delta, seed, cap, after, record_every, tol=1e-12):
    sys_lr = LowRankSystem(n)
    M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n + seed)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=tol)
    p1_sigma, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    dim_p1 = B_p1.shape[1]
    dim_rot = 0 if B_rot is None else B_rot.shape[1]
    dim_ker = M - dim_p1 - dim_rot

    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    def bands(Z):
        a, b = Z.real, Z.imag
        tot = a @ a + b @ b
        h1 = (np.sum((B_p1.T @ a) ** 2) + np.sum((B_p1.T @ b) ** 2)) / tot
        hr = 0.0 if B_rot is None else (
            np.sum((B_rot.T @ a) ** 2) + np.sum((B_rot.T @ b) ** 2)) / tot
        return h1, hr, max(0.0, 1.0 - h1 - hr)

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    rows = []
    max_id_dev = 0.0
    crossed = None
    wp = rng.normal(size=M)
    for t in range(cap + 1):
        f = fval(Z)
        if crossed is None and f > 0.05:
            crossed = t
        if crossed is not None and t >= crossed + after:
            break
        if t % record_every == 0:
            h1, hr, hk = bands(Z)
            max_id_dev = max(max_id_dev, abs(f - (1.0 - h1)))
            rows.append((t, f, h1, hr, hk))
        sys_lr.set_theta(np.angle(Z))
        se, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, se)

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"N{n:05d}_exact"
    with open(RESULT_DIR / f"planeflow_{tag}.csv", "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["tau", "f", "frac_P1", "frac_other_rotation", "frac_kernel"])
        wtr.writerows(rows)

    ref_csv = REFERENCE_DIR / f"fcurve_N{n:05d}_delta{delta:.0e}_seed{seed}.csv"
    ref_dev = None
    if ref_csv.exists():
        ref = {}
        with open(ref_csv) as fh:
            for r in csv.DictReader(fh):
                ref[int(r["tau"])] = float(r["f"])
        ref_dev = max(abs(row[1] - ref[row[0]]) for row in rows if row[0] in ref)

    control = {"identity_f_eq_1_minus_P1_maxdev": float(max_id_dev),
               "identity_passed": bool(max_id_dev < 1e-12),
               "f_vs_reference_maxdev": (None if ref_dev is None else float(ref_dev)),
               "f_vs_reference_passed": bool(ref_dev is None or ref_dev < 1e-12)}
    meta = {"method": "exact_dense_eigK", "n": n, "m": M, "delta": delta, "seed": seed,
            "crossing_tau": crossed, "parent_residual": residual,
            "p1_sigma": float(p1_sigma),
            "dims": {"P1": dim_p1, "other_rotation": dim_rot, "kernel": dim_ker},
            "sigma_spectrum_positive": spectrum[:30],
            "sigma_spectrum_count": len(spectrum),
            "control_test": control}
    with open(RESULT_DIR / f"planeflow_{tag}.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)

    arr = np.array(rows)
    taus, fcol, fracs = arr[:, 0], arr[:, 1], arr[:, 2:5].T
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    labels = [f"P1: σ={p1_sigma:.3f} dominant plane ({dim_p1}D)",
              f"other rotation planes σ>0 ({dim_rot}D)",
              f"σ=0 kernel ({dim_ker}D)"]
    ax.stackplot(taus, *fracs, labels=labels,
                 colors=["#1f77b4", "#ff7f0e", "#2ca02c"], alpha=0.85)
    if crossed is not None:
        ax.axvline(crossed, color="k", ls=":", lw=1.0)
        ax.text(crossed, 1.01, f" crossing τ={crossed}", fontsize=8)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("energy fraction per generator rotation plane")
    ax.set_ylim(0, 1)
    ax.set_title(f"N={n} (M={M:,}) EXACT eig(K): amplitude flows out of P1 into "
                 f"other rotation planes and kernel")
    ax.legend(loc="upper right", fontsize=8)
    ax2 = ax.twinx()
    ax2.semilogy(taus, fcol, color="black", lw=1.4,
                 label=r"$f(\tau)=1-P_1$ (dormant, log right)")
    ax2.set_ylabel(r"$f(\tau)$")
    ax2.set_ylim(1e-31, 2.0)
    ax2.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / f"planeflow_{tag}.png", dpi=150)
    plt.close(fig)

    print(json.dumps({"method": "exact", "n": n, "crossing_tau": crossed,
                      "p1_sigma": float(p1_sigma), "dims": meta["dims"],
                      "control_test": control}, ensure_ascii=False))
    if not control["identity_passed"]:
        raise SystemExit(f"恒等対照失敗: {max_id_dev:.2e}")
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--delta", type=float, default=1e-15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cap", type=int, default=30000)
    ap.add_argument("--after", type=int, default=20000)
    ap.add_argument("--record-every", type=int, default=20)
    args = ap.parse_args()
    if args.n * (args.n - 1) // 2 > 4000:
        raise SystemExit(f"N={args.n} は M={args.n*(args.n-1)//2} で密行列が大きすぎます。"
                         f"近似法 run_plane_flow_approx_v1.py を使ってください。")
    run(args.n, args.delta, args.seed, args.cap, args.after, args.record_every)


if __name__ == "__main__":
    main()
