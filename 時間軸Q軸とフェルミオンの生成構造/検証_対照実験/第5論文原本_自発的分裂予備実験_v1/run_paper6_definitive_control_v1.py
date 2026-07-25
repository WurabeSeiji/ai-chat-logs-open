#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 論文6追試 確定対照：N=5,40,300 同一軌道・二観測法比較（図1〜8）。

論文6の三代表例 N=5,40,300 について、長時間実験系をそのまま使い、各Nで一本の軌道だけを
生成する。その同一状態列へ二観測法を同時適用する：
  A 論文6型固定基底（parent_plane_split_exact/approx + bands, verbatim）：h1,hr,hk,f
  B 第7論文型瞬時分解：E_dom,E_nondom,E_ker,E_H,有効ランク,親平面重なり,支配平面歳差
N=5,40 は厳密 eig(K)、N=300 は低ランク JG。N=40 で厳密vs近似を交差検証（事前許容誤差）。

パラメータ調整禁止。原本エンジン不変更。観測窓は事前固定。
使い方:
  python3 run_paper6_definitive_control_v1.py obs 5 40      # 観測（速い）
  python3 run_paper6_definitive_control_v1.py obs 300       # 近似・約13分
  python3 run_paper6_definitive_control_v1.py crossval 40   # N=40 厳密vs近似
  python3 run_paper6_definitive_control_v1.py figs          # 図1〜8
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_plane_flow_approx_v1 import parent_plane_split_approx
from run_halfband_stability_a2_v1 import occ, proj_overlap
from run_transverse_stability_v1 import DELTA, GUARD, LEARN, VALID

RESULT_DIR = HERE / "paper6_definitive_control_v1"
SIGREL = 1e-6
REC_INST = {5: 50, 40: 50, 300: 200}
AFTER = 50000


# ---------- 瞬時分解（exact: dense eig(K) / approx: 低ランク JG） ----------
def inst_bands(sys_lr, Z, method, B_parent):
    sys_lr.set_theta(np.angle(Z))
    if method == "exact":
        M = sys_lr.m
        K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
        w, V = np.linalg.eig(K)
        cols = [(float(w[i].imag), V[:, i].real.copy(), V[:, i].imag.copy()) for i in range(M)]
    else:
        ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
        cols = []
        for i in range(len(ev)):
            lifted = sys_lr.w(EV[:, i].astype(complex))
            cols.append((float(ev[i].imag), np.real(lifted).copy(), np.imag(lifted).copy()))
    smax = max(c[0] for c in cols)
    thr = SIGREL * smax
    groups = defaultdict(list)
    for (s, vr, vi) in cols:
        if s > thr:
            groups[round(s, 4)].append((vr, vi))
    sig_sorted = sorted(groups, reverse=True)

    def ortho(items):
        cm = []
        for (vr, vi) in items:
            cm.append(vr); cm.append(vi)
        Q, R = np.linalg.qr(np.column_stack(cm))
        return Q[:, np.abs(np.diag(R)) > 1e-8]

    def ortho_minus(items, B):
        cm = []
        for (vr, vi) in items:
            cm.append(vr); cm.append(vi)
        R0 = np.column_stack(cm)
        if B is not None:
            R0 = R0 - B @ (B.T @ R0)
        Q, R = np.linalg.qr(R0)
        return Q[:, np.abs(np.diag(R)) > 1e-8]

    B_dom = ortho(groups[sig_sorted[0]])
    rest = [it for k in sig_sorted[1:] for it in groups[k]]
    B_non = ortho_minus(rest, B_dom) if rest else None
    hb = [it for k in sig_sorted if 0.4 * smax <= k <= 0.6 * smax for it in groups[k]]
    B_H = ortho_minus(hb, B_dom) if hb else None
    Edom = occ(B_dom, Z)
    Enon = occ(B_non, Z) if B_non is not None else 0.0
    EH = occ(B_H, Z) if B_H is not None else 0.0
    return {"smax": smax, "Edom": Edom, "Enon": Enon, "Eker": max(0.0, 1 - Edom - Enon),
            "EH": EH, "B_dom": B_dom, "ov_parent": proj_overlap(B_parent, B_dom),
            "n_planes": len(sig_sorted)}


def build_init(n):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    if n <= 40:
        p1s, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    else:
        p1s, B_p1, B_rot, smax, thr = parent_plane_split_approx(sys_lr, v, SIGREL)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=sys_lr.m)
    return sys_lr, v, B_p1, B_rot, Z, p, q, wp


def paper6_bands(Z, B_p1, B_rot):
    a, b = Z.real, Z.imag
    tot = a @ a + b @ b
    h1 = (np.sum((B_p1.T @ a) ** 2) + np.sum((B_p1.T @ b) ** 2)) / tot
    hr = 0.0 if B_rot is None else (np.sum((B_rot.T @ a) ** 2) + np.sum((B_rot.T @ b) ** 2)) / tot
    return float(h1), float(hr), float(max(0.0, 1 - h1 - hr))


def fval(Z, p, q):
    Zp = Z - p * (p @ Z) - q * (q @ Z)
    return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))


def observe(n):
    sys_lr, v, B_p1, B_rot, Z, p, q, wp = build_init(n)
    method = "exact" if n <= 40 else "approx"
    rec = REC_INST[n]
    rows = []
    frames = []
    B_dom0 = None
    crossed = None
    stop_t = None
    iddev = 0.0
    t = 0
    while True:
        f = fval(Z, p, q)
        if crossed is None and f > 0.05:
            crossed = t; stop_t = crossed + AFTER
        if t % rec == 0:
            h1, hr, hk = paper6_bands(Z, B_p1, B_rot)
            iddev = max(iddev, abs(f - (1 - h1)))
            ib = inst_bands(sys_lr, Z, method, B_p1)
            frames.append(Z.copy())
            if len(frames) > 20:
                frames.pop(0)
            X = np.column_stack([c for ZZ in frames for c in (ZZ.real, ZZ.imag)])
            s = np.linalg.svd(X, compute_uv=False); lam = s ** 2
            er = float((lam.sum() ** 2) / np.sum(lam ** 2))
            if B_dom0 is None:
                B_dom0 = ib["B_dom"]
            selfprec = proj_overlap(B_dom0, ib["B_dom"])
            rows.append((t, f, h1, hr, hk, ib["Edom"], ib["Enon"], ib["Eker"], ib["EH"],
                         er, ib["ov_parent"], selfprec, ib["smax"], ib["n_planes"]))
        if stop_t is not None and t >= stop_t:
            break
        sys_lr.set_theta(np.angle(Z))
        se, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, se)
        t += 1

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    hdr = ["tau", "f", "h1", "hr", "hk", "E_dom", "E_nondom", "E_ker", "E_H",
           "eff_rank", "dom_parent_overlap", "dom_self_precession", "sigma_dom", "n_planes"]
    with open(RESULT_DIR / f"obs_N{n:05d}.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(hdr); w.writerows(rows)
    arr = np.array(rows)

    # 窓（事前固定）
    windows = {
        "paper6[cross+100,cross+20000]": (crossed + 100, crossed + 20000),
        "long[cross+100,cross+50000]": (crossed + 100, crossed + 50000),
    }

    def agg(lo, hi):
        m = (arr[:, 0] >= lo) & (arr[:, 0] <= hi)
        w = arr[m]
        return {"f_med": float(np.median(w[:, 1])), "hr_med": float(np.median(w[:, 3])),
                "hk_med": float(np.median(w[:, 4])), "Edom_med": float(np.median(w[:, 5])),
                "Enon_med": float(np.median(w[:, 6])), "Enon_end": float(w[-1, 6]),
                "EH_med": float(np.median(w[:, 8])), "effrank_med": float(np.median(w[:, 9])),
                "effrank_end": float(w[-1, 9]), "ovparent_med": float(np.median(w[:, 10]))}
    # 非支配減衰率
    tt = arr[:, 0]; en = np.maximum(arr[:, 6], 1e-300)
    mdec = (tt > crossed) & (en > 1e-200)
    decay = float(np.polyfit(tt[mdec], np.log(en[mdec]), 1)[0]) if mdec.sum() > 4 else None
    summary = {
        "n": n, "m": int(sys_lr.m), "method_inst": method, "crossing": crossed,
        "identity_maxdev": float(iddev), "n_planes_final": int(arr[-1, 13]),
        "windows": {k: agg(*v) for k, v in windows.items()},
        "Enondom_decay_rate_per_step": decay,
        "final": {"f": float(arr[-1, 1]), "hr": float(arr[-1, 3]), "hk": float(arr[-1, 4]),
                  "E_dom": float(arr[-1, 5]), "E_nondom": float(arr[-1, 6]),
                  "eff_rank": float(arr[-1, 9]), "dom_parent_overlap": float(arr[-1, 10])},
    }
    with open(RESULT_DIR / f"obs_summary_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"[obs] N={n} M={sys_lr.m} 交差={crossed} 恒等式偏差={iddev:.1e} 平面数={int(arr[-1,13])}")
    for k, v in summary["windows"].items():
        print(f"  {k}: f={v['f_med']:.3f} hr={v['hr_med']:.4f} | E_dom={v['Edom_med']:.6f} "
              f"E_nondom(中)={v['Enon_med']:.2e} E_nondom(終)={v['Enon_end']:.1e} "
              f"有効ランク={v['effrank_med']:.3f}→{v['effrank_end']:.3f} 親重なり={v['ovparent_med']:.3f}")
    print(f"  非支配減衰率={decay:.2e}/step  最終: E_dom={summary['final']['E_dom']:.6f} "
          f"E_nondom={summary['final']['E_nondom']:.1e} 有効ランク={summary['final']['eff_rank']:.3f}")
    return summary


# ---------- N=40 厳密vs近似 交差検証 ----------
def crossval(n=40):
    sys_lr, v, B_p1, B_rot, Z, p, q, wp = build_init(n)
    rec = REC_INST[n]
    rows = []
    crossed = None; stop_t = None
    t = 0
    while True:
        f = fval(Z, p, q)
        if crossed is None and f > 0.05:
            crossed = t; stop_t = crossed + AFTER
        if t % rec == 0:
            ex = inst_bands(sys_lr, Z, "exact", B_p1)
            ap = inst_bands(sys_lr, Z, "approx", B_p1)
            dom_overlap = proj_overlap(ex["B_dom"], ap["B_dom"])
            rows.append((t, abs(ex["smax"] - ap["smax"]), dom_overlap,
                         abs(ex["Edom"] - ap["Edom"]), abs(ex["Enon"] - ap["Enon"]),
                         abs(ex["Eker"] - ap["Eker"]), abs(ex["EH"] - ap["EH"]),
                         abs(ex["ov_parent"] - ap["ov_parent"])))
        if stop_t is not None and t >= stop_t:
            break
        sys_lr.set_theta(np.angle(Z))
        se, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, se)
        t += 1
    arr = np.array(rows)
    tol = {"sigma": 1e-6, "dom_overlap_min": 1 - 1e-8, "occ": 1e-6}
    res = {
        "n": n, "n_records": len(rows),
        "max_sigma_dev": float(arr[:, 1].max()),
        "min_dom_plane_overlap": float(arr[:, 2].min()),
        "max_Edom_dev": float(arr[:, 3].max()), "max_Enon_dev": float(arr[:, 4].max()),
        "max_Eker_dev": float(arr[:, 5].max()), "max_EH_dev": float(arr[:, 6].max()),
        "max_ovparent_dev": float(arr[:, 7].max()),
        "tolerance": tol,
        "passed": bool(arr[:, 1].max() < 1e-4 and arr[:, 2].min() > 1 - 1e-6
                       and arr[:, 3].max() < 1e-5 and arr[:, 4].max() < 1e-5
                       and arr[:, 6].max() < 1e-5),
    }
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULT_DIR / f"crossval_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print(f"[crossval] N={n} 厳密vs近似（{len(rows)}時刻）:")
    print(f"  σ偏差max={res['max_sigma_dev']:.2e} 支配平面重なりmin={res['min_dom_plane_overlap']:.8f}")
    print(f"  E_dom偏差max={res['max_Edom_dev']:.2e} E_nondom偏差max={res['max_Enon_dev']:.2e} "
          f"E_H偏差max={res['max_EH_dev']:.2e} 親重なり偏差max={res['max_ovparent_dev']:.2e}")
    print(f"  → 合格={res['passed']}")
    return res


def load(n):
    arr = []
    with open(RESULT_DIR / f"obs_N{n:05d}.csv") as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            arr.append([float(x) for x in row])
    return np.array(arr)


def make_figures(ns=(5, 40, 300)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = {}
    for n in ns:
        try:
            data[n] = load(n)
        except FileNotFoundError:
            print(f"  N={n} の obs CSV なし、図から除外")
    ns = [n for n in ns if n in data]
    C = {"h1": "#4C78A8", "hr": "#F58518", "hk": "#B0B0B0",
         "dom": "#4C78A8", "non": "#E45756", "ker": "#B0B0B0", "rank": "#54A24B",
         "ov": "#B279A2"}
    # 列インデックス: 0tau 1f 2h1 3hr 4hk 5Edom 6Enon 7Eker 8EH 9rank 10ovpar 11selfprec 12sig 13np
    def crossing(a):
        idx = np.argmax(a[:, 1] > 0.05); return a[idx, 0]

    # 図1: 二観測法 直接比較（上段固定基底 stack / 下段瞬時）
    fig, axes = plt.subplots(2, len(ns), figsize=(5 * len(ns), 8), squeeze=False)
    for j, n in enumerate(ns):
        a = data[n]; tau = a[:, 0]
        axes[0][j].stackplot(tau, a[:, 2], a[:, 3], a[:, 4],
                             labels=["h1 parent-dom", "hr other-rot", "hk kernel"],
                             colors=[C["h1"], C["hr"], C["hk"]], alpha=0.85)
        axes[0][j].plot(tau, a[:, 1], "k-", lw=1.0, label="f")
        axes[0][j].set_title(f"N={n} (M={int(a[0,13]) and '' }) — top: Paper6 fixed basis")
        axes[0][j].set_ylim(0, 1); axes[0][j].set_title(f"N={n}: Paper6 fixed basis")
        if j == 0: axes[0][j].legend(fontsize=7, loc="center right")
        axes[1][j].semilogy(tau, np.maximum(a[:, 6], 1e-20), color=C["non"], label="E_nondom (inst)")
        axes[1][j].semilogy(tau, np.maximum(1 - a[:, 5], 1e-20), color=C["dom"], ls="--", lw=0.9,
                            label="1-E_dom (inst)")
        ax2 = axes[1][j].twinx(); ax2.plot(tau, a[:, 9], color=C["rank"], lw=1.2, label="eff rank")
        ax2.set_ylim(0, 3.2)
        axes[1][j].set_title(f"N={n}: inst decomposition"); axes[1][j].set_xlabel("tau")
        if j == 0: axes[1][j].legend(fontsize=7, loc="center left")
        if j == len(ns) - 1: ax2.legend(fontsize=7, loc="center right")
    fig.suptitle("Fig1: SAME trajectory, two observation methods (top: multi-plane appears; bottom: inst non-dominant decays)")
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig1_two_methods.png", dpi=130); plt.close(fig)

    # 図2: f_fixed vs 1-E_dom（通常軸・片対数）
    fig, axes = plt.subplots(2, len(ns), figsize=(5 * len(ns), 7), squeeze=False)
    for j, n in enumerate(ns):
        a = data[n]; tau = a[:, 0]
        for row, logy in [(0, False), (1, True)]:
            ax = axes[row][j]
            if logy:
                ax.semilogy(tau, np.maximum(a[:, 1], 1e-20), "k-", label="f_fixed=1-h1")
                ax.semilogy(tau, np.maximum(1 - a[:, 5], 1e-20), color=C["non"], label="1-E_dom (inst)")
            else:
                ax.plot(tau, a[:, 1], "k-", label="f_fixed"); ax.plot(tau, 1 - a[:, 5], color=C["non"], label="1-E_dom")
            ax.set_title(f"N={n}" + (" (log)" if logy else ""));
            if j == 0 and row == 0: ax.legend(fontsize=7)
    fig.suptitle("Fig2: fixed-basis splitting f -> nonzero plateau, but inst 1-E_dom -> 0")
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig2_f_vs_instdom.png", dpi=130); plt.close(fig)

    # 図3: hr vs E_nondom
    fig, axes = plt.subplots(1, len(ns), figsize=(5 * len(ns), 4), squeeze=False)
    for j, n in enumerate(ns):
        a = data[n]; tau = a[:, 0]
        ax = axes[0][j]
        ax.semilogy(tau, np.maximum(a[:, 3], 1e-20), color=C["hr"], label="hr (fixed other-rot)")
        ax.semilogy(tau, np.maximum(a[:, 6], 1e-20), color=C["non"], label="E_nondom (inst)")
        ax.set_title(f"N={n}"); ax.set_xlabel("tau")
        if j == 0: ax.legend(fontsize=7)
    fig.suptitle("Fig3: hr (basis for 'second plane') plateaus; inst E_nondom decays")
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig3_hr_vs_Enondom.png", dpi=130); plt.close(fig)

    # 図4: 歳差による再構成 h1 vs overlap^? / 実測h1 と 支配平面重なりで予測したh1
    fig, axes = plt.subplots(1, len(ns), figsize=(5 * len(ns), 4), squeeze=False)
    for j, n in enumerate(ns):
        a = data[n]; tau = a[:, 0]
        ax = axes[0][j]
        ax.plot(tau, a[:, 2], color=C["h1"], label="h1 measured (fixed)")
        # 予測: 状態は瞬時支配平面にほぼ全占有 → h1 ≈ overlap(parent,dom)
        ax.plot(tau, a[:, 10], color=C["ov"], ls="--", label="dom-parent overlap (predict h1)")
        ax.set_title(f"N={n}"); ax.set_xlabel("tau"); ax.set_ylim(0, 1.02)
        if j == 0: ax.legend(fontsize=7)
    fig.suptitle("Fig4: fixed h1 reconstructed by dominant-plane reorientation (precession)")
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig4_precession_reconstruct.png", dpi=130); plt.close(fig)

    # 図5: E_nondom 片対数
    fig, ax = plt.subplots(figsize=(8, 5))
    for n in ns:
        a = data[n]
        ax.semilogy(a[:, 0] - crossing(a), np.maximum(a[:, 6], 1e-20), label=f"N={n}")
    ax.set_xlabel("tau - crossing"); ax.set_ylabel("E_nondom (inst, log)")
    ax.set_title("Fig5: inst non-dominant occupation long-time decay"); ax.legend()
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig5_Enondom_decay.png", dpi=130); plt.close(fig)

    # 図6: 有効ランク
    fig, ax = plt.subplots(figsize=(8, 5))
    for n in ns:
        a = data[n]
        ax.plot(a[:, 0] - crossing(a), a[:, 9], label=f"N={n}")
    ax.axhline(2.0, color="k", ls=":", lw=0.8); ax.set_ylim(1.8, 3.0)
    ax.set_xlabel("tau - crossing"); ax.set_ylabel("effective rank")
    ax.set_title("Fig6: effective rank -> 2 (single mode)"); ax.legend()
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig6_eff_rank.png", dpi=130); plt.close(fig)

    # 図7: 短時間窓 vs 長時間窓（N毎に E_nondom と hr、窓を色分け）
    fig, axes = plt.subplots(1, len(ns), figsize=(5 * len(ns), 4), squeeze=False)
    for j, n in enumerate(ns):
        a = data[n]; cr = crossing(a); tau = a[:, 0]
        ax = axes[0][j]
        ax.semilogy(tau, np.maximum(a[:, 6], 1e-20), color=C["non"], label="E_nondom (inst)")
        ax.axvspan(cr + 100, cr + 20000, color="orange", alpha=0.10, label="Paper6 window")
        ax.axvspan(cr + 20000, cr + AFTER, color="green", alpha=0.08, label="long window")
        ax.set_title(f"N={n}"); ax.set_xlabel("tau")
        if j == 0: ax.legend(fontsize=7)
    fig.suptitle("Fig7: short (Paper6) window vs long window — multi-plane appearance is transient")
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig7_windows.png", dpi=130); plt.close(fig)

    # 図8: N依存統合
    fig, ax = plt.subplots(figsize=(8, 5))
    hr_pl = []; en_fin = []; ed_fin = []; rk_fin = []; ov_fin = []
    for n in ns:
        a = data[n]
        hr_pl.append(np.median(a[a[:, 0] > crossing(a) + 20000, 3]))
        en_fin.append(a[-1, 6]); ed_fin.append(a[-1, 5]); rk_fin.append(a[-1, 9]); ov_fin.append(a[-1, 10])
    x = np.array(ns)
    ax.plot(x, hr_pl, "o-", label="hr plateau (fixed)")
    ax.plot(x, ov_fin, "s-", label="final dom-parent overlap")
    ax.plot(x, np.array(rk_fin) - 2, "^-", label="final eff_rank - 2")
    ax.set_xscale("log"); ax.set_xlabel("N"); ax.set_ylabel("value")
    ax.set_title("Fig8: N-dependence (3 points; no continuous law claimed)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RESULT_DIR / "fig8_Ndependence.png", dpi=130); plt.close(fig)
    print(f"[figs] 図1〜8 を {RESULT_DIR} に出力（N={ns}）")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["obs", "crossval", "figs", "all"])
    ap.add_argument("ns", type=int, nargs="*", default=[])
    args = ap.parse_args()
    if args.cmd == "obs":
        for n in args.ns:
            observe(n)
    elif args.cmd == "crossval":
        crossval(args.ns[0] if args.ns else 40)
    elif args.cmd == "figs":
        make_figures(tuple(args.ns) if args.ns else (5, 40, 300))


if __name__ == "__main__":
    main()
