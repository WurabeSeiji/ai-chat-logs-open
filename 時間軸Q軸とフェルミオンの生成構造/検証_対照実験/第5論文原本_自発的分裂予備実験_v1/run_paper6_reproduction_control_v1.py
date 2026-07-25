#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 論文6追試対照：同一軌道・二観測法の窓別比較。

A3の基準軌道（パラメータ・seed・初期条件・更新則・σ法・精度を一切変更しない）を最初から
crossing+50000 まで一度だけ走らせ、同じ一本の時系列に対して二つの観測法を同時適用する：

 (P6) 論文6型（固定基底）：t=0 の親生成子から確定した固定部分空間 B_P1(親支配平面),
      B_rot(親のその他回転平面) への占有 h1,hr,hk と、固定親平面分裂 f=1-|Π_{P0}Z|²。
      ← run_plane_flow_exact_v1 の parent_plane_split_exact / bands を verbatim 使用。
 (A3) 瞬時型：各時刻の瞬時生成子 K(argZ) を分解した 瞬時支配平面占有・全非支配回転帯占有・
      核占有・有効ランク・瞬時支配平面と固定親平面の重なり。

論文6の値へ合わせ込まない。変更するのは観測終了時刻だけ。窓は事前固定。
使い方: python3 run_paper6_reproduction_control_v1.py 10 40 --after 50000
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_halfband_stability_a2_v1 import plane_decomp, band_bases, occ, proj_overlap
from run_transverse_stability_v1 import DELTA, GUARD, LEARN, VALID

RESULT_DIR = HERE / "paper6_reproduction_control_v1"
REC_P6 = 20


def run_N(n, after):
    sys_lr = LowRankSystem(n)
    M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n)          # A3/論文6と同一seed
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    # --- 論文6の固定基底（親生成子, ループ前に一度だけ確定）---
    p1_sigma, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    # --- A3/論文6と同一の初期状態 ---
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=M)
    rec_a3 = 100 if n <= 10 else 200

    def paper6_bands(Z):
        a, b = Z.real, Z.imag
        tot = a @ a + b @ b
        h1 = (np.sum((B_p1.T @ a) ** 2) + np.sum((B_p1.T @ b) ** 2)) / tot
        hr = 0.0 if B_rot is None else (
            np.sum((B_rot.T @ a) ** 2) + np.sum((B_rot.T @ b) ** 2)) / tot
        return float(h1), float(hr), float(max(0.0, 1 - h1 - hr))

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    rows_p6 = []      # (tau, f, h1, hr, hk)
    rows_a3 = []      # (tau, inst_dom_occ, inst_nondom_band, inst_kernel, dom_parent_overlap)
    a3_frames = []    # (tau, Z) trailing for eff_rank
    crossed = None
    stop_t = None
    id_dev = 0.0
    t = 0
    while True:
        f = fval(Z)
        if crossed is None and f > 0.05:
            crossed = t
            stop_t = crossed + after
        if t % REC_P6 == 0:
            h1, hr, hk = paper6_bands(Z)
            id_dev = max(id_dev, abs(f - (1 - h1)))
            rows_p6.append((t, f, h1, hr, hk))
        if t % rec_a3 == 0:
            planes, smax = plane_decomp(sys_lr, Z)
            Bd, Bh, sb = band_bases(planes, smax)
            dom = occ(Bd, Z); band = occ(Bh, Z)
            ov = proj_overlap(B_p1, Bd)
            a3_frames.append((t, Z.copy()))
            if len(a3_frames) > 20:
                a3_frames.pop(0)
            X = np.column_stack([c for (_, ZZ) in a3_frames for c in (ZZ.real, ZZ.imag)])
            s = np.linalg.svd(X, compute_uv=False); lam = s ** 2
            er = float((lam.sum() ** 2) / np.sum(lam ** 2))
            rows_a3.append((t, float(dom), float(band), float(1 - dom - band), float(ov), er))
        if stop_t is not None and t >= stop_t:
            break
        sys_lr.set_theta(np.angle(Z))
        se, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, se)
        t += 1

    # --- 窓（事前固定, 観測終了時刻だけ変える）---
    windows = {
        "paper6_short[cross+100,cross+20000]": (crossed + 100, crossed + 20000),
        "A1[cross+500,cross+2500]": (crossed + GUARD, crossed + GUARD + LEARN + VALID),
        "end_10000[cross+100,cross+10000]": (crossed + 100, crossed + 10000),
        "end_20000[cross+100,cross+20000]": (crossed + 100, crossed + 20000),
        "end_50000[cross+100,cross+50000]": (crossed + 100, crossed + 50000),
    }
    p6 = np.array(rows_p6)          # tau,f,h1,hr,hk
    a3 = np.array(rows_a3)          # tau,dom,band,ker,ov,er

    def agg(lo, hi):
        mp = (p6[:, 0] >= lo) & (p6[:, 0] <= hi)
        ma = (a3[:, 0] >= lo) & (a3[:, 0] <= hi)
        wp6, wa3 = p6[mp], a3[ma]
        return {
            "P6_f_median": float(np.median(wp6[:, 1])),
            "P6_h1_median(parent_dom)": float(np.median(wp6[:, 2])),
            "P6_hr_median(other_rotation)": float(np.median(wp6[:, 3])),
            "P6_hk_median(kernel)": float(np.median(wp6[:, 4])),
            "A3_inst_dom_occ_median": float(np.median(wa3[:, 1])),
            "A3_inst_nondom_band_median": float(np.median(wa3[:, 2])),
            "A3_inst_nondom_band_end": float(wa3[-1, 2]),
            "A3_kernel_median": float(np.median(wa3[:, 3])),
            "A3_dom_parent_overlap_median": float(np.median(wa3[:, 4])),
            "A3_eff_rank_median": float(np.median(wa3[:, 5])),
            "A3_eff_rank_end": float(wa3[-1, 5]),
        }

    window_stats = {name: agg(lo, hi) for name, (lo, hi) in windows.items()}

    report = {
        "n": n, "m": M, "crossing_tau": crossed, "after": after,
        "paper6_fixed_basis_dims": {"dim_P1": int(B_p1.shape[1]),
                                    "dim_rot": int(0 if B_rot is None else B_rot.shape[1]),
                                    "dim_kernel": int(M - B_p1.shape[1] - (0 if B_rot is None else B_rot.shape[1]))},
        "identity_f_eq_1_minus_h1_maxdev": float(id_dev),
        "windows": window_stats,
        "final_P6": {"f": float(p6[-1, 1]), "h1": float(p6[-1, 2]),
                     "hr": float(p6[-1, 3]), "hk": float(p6[-1, 4])},
        "final_A3": {"inst_dom_occ": float(a3[-1, 1]), "inst_nondom_band": float(a3[-1, 2]),
                     "dom_parent_overlap": float(a3[-1, 4]), "eff_rank": float(a3[-1, 5])},
    }
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULT_DIR / f"control_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    with open(RESULT_DIR / f"p6curve_N{n:05d}.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["tau", "f", "h1_parent_dom", "hr_other_rotation", "hk_kernel"])
        w.writerows(rows_p6)
    with open(RESULT_DIR / f"a3curve_N{n:05d}.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["tau", "inst_dom_occ", "inst_nondom_band", "inst_kernel",
                                        "dom_parent_overlap", "eff_rank"])
        w.writerows(rows_a3)
    make_figure(n, rows_p6, rows_a3, crossed, windows)
    return report


def make_figure(n, rows_p6, rows_a3, crossed, windows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    p6 = np.array(rows_p6); a3 = np.array(rows_a3)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    # Panel A: Paper 6 fixed-basis stackplot + f
    ax1.stackplot(p6[:, 0], p6[:, 2], p6[:, 3], p6[:, 4],
                  labels=["h1: parent dominant plane (fixed)",
                          "hr: other rotation planes (fixed)", "hk: kernel"],
                  colors=["#4C78A8", "#F58518", "#B0B0B0"], alpha=0.85)
    ax1b = ax1.twinx()
    ax1b.semilogy(p6[:, 0], np.maximum(p6[:, 1], 1e-12), "k-", lw=1.2, label="f = 1-h1 (log)")
    ax1.set_ylabel("Paper6 fixed-basis occupation")
    ax1b.set_ylabel("f (log)")
    ax1.set_title(f"N={n}: SAME trajectory, two observation methods")
    ax1.legend(loc="center right", fontsize=8); ax1.set_ylim(0, 1)
    # Panel B: A3 instantaneous
    ax2.semilogy(a3[:, 0], np.maximum(a3[:, 2], 1e-20), color="#E45756", lw=1.4,
                 label="inst. non-dominant band occ (log)")
    ax2.semilogy(a3[:, 0], np.maximum(1 - a3[:, 1], 1e-20), color="#72B7B2", lw=1.0, ls="--",
                 label="1 - inst. dominant occ (log)")
    ax2b = ax2.twinx()
    ax2b.plot(a3[:, 0], a3[:, 5], color="#54A24B", lw=1.2, label="effective rank")
    ax2b.plot(a3[:, 0], a3[:, 4], color="#B279A2", lw=1.0, ls=":", label="dom-parent overlap")
    ax2.set_ylabel("inst. non-dominant (log)"); ax2b.set_ylabel("eff rank / overlap")
    ax2.set_xlabel("tau"); ax2b.set_ylim(0, 3.2)
    ax2.legend(loc="center left", fontsize=8); ax2b.legend(loc="center right", fontsize=8)
    for ax in (ax1, ax2):
        ax.axvline(crossed, color="k", ls=":", lw=0.8)
        for (lo, hi) in [windows["A1[cross+500,cross+2500]"],
                         windows["paper6_short[cross+100,cross+20000]"]]:
            ax.axvspan(lo, hi, color="yellow", alpha=0.06)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / f"control_overlay_N{n:05d}.png", dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", type=int, nargs="+")
    ap.add_argument("--after", type=int, default=50000)
    args = ap.parse_args()
    for n in args.ns:
        r = run_N(n, args.after)
        print(f"\n===== N={n} (M={r['m']}) crossing={r['crossing_tau']}  "
              f"恒等式 |f-(1-h1)|max={r['identity_f_eq_1_minus_h1_maxdev']:.1e}  "
              f"固定基底次元 P1={r['paper6_fixed_basis_dims']['dim_P1']} "
              f"rot={r['paper6_fixed_basis_dims']['dim_rot']} ker={r['paper6_fixed_basis_dims']['dim_kernel']} =====")
        print(f"{'窓':<38}{'P6:f':>7}{'P6:hr':>8}{'P6:hk':>8} | {'A3瞬時支配':>9}{'A3瞬時非支配(中)':>13}{'A3非支配(終)':>11}{'A3有効ランク':>10}{'歳差重なり':>9}")
        for name, s in r["windows"].items():
            print(f"{name:<38}{s['P6_f_median']:>7.3f}{s['P6_hr_median(other_rotation)']:>8.3f}"
                  f"{s['P6_hk_median(kernel)']:>8.1e} | {s['A3_inst_dom_occ_median']:>9.5f}"
                  f"{s['A3_inst_nondom_band_median']:>13.2e}{s['A3_inst_nondom_band_end']:>11.1e}"
                  f"{s['A3_eff_rank_median']:>10.3f}{s['A3_dom_parent_overlap_median']:>9.3f}")
        fa = r["final_A3"]; fp = r["final_P6"]
        print(f"  最終(cross+{args.after}): P6 f={fp['f']:.3f} hr={fp['hr']:.3f} | "
              f"A3 瞬時支配={fa['inst_dom_occ']:.6f} 瞬時非支配={fa['inst_nondom_band']:.1e} "
              f"有効ランク={fa['eff_rank']:.3f} 歳差重なり={fa['dom_parent_overlap']:.3f}")


if __name__ == "__main__":
    main()
