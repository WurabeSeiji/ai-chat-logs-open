#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 第1段階：N=5,40 厳密固有スペクトル観測の図（1A〜8）。解釈なし。

raw/ の CSV を後処理して図を生成する。全枝を残し、小枝を消さない。
使い方: python3 make_figures_exact_lowN_v1.py 5
        python3 make_figures_exact_lowN_v1.py 40
        python3 make_figures_exact_lowN_v1.py comparison 5 40
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE_DIR = Path(__file__).resolve().parent
BASE = CODE_DIR.parent


def load_rows(path):
    with open(path) as fh:
        r = csv.DictReader(fh)
        return list(r), r.fieldnames


def pivot(path, key_col, val_col, time_col="time"):
    rows, _ = load_rows(path)
    d = defaultdict(dict)
    keys = set()
    for row in rows:
        t = int(float(row[time_col])); k = int(float(row[key_col]))
        d[t][k] = float(row[val_col])
        keys.add(k)
    times = sorted(d)
    keys = sorted(keys)
    M = np.full((len(times), len(keys)), np.nan)
    ki = {k: i for i, k in enumerate(keys)}
    for i, t in enumerate(times):
        for k, val in d[t].items():
            M[i, ki[k]] = val
    return np.array(times), keys, M


def crossing_of(n):
    import json
    return json.load(open(BASE / "diagnostics" / f"N{n:05d}.json"))["crossing"]


def figdir(n):
    d = BASE / "figures" / f"N{n:05d}"; d.mkdir(parents=True, exist_ok=True); return d


def make_all(n):
    raw = BASE / "raw" / f"N{n:05d}"
    fd = figdir(n)
    cr = crossing_of(n)
    ev = raw / "eigenvalues.csv"

    # --- 図1A/1B: σ_j/σ_1 全順位 線形/対数 ---
    t, ranks, R = pivot(ev, "rank_index", "sigma_over_sigma1")
    for logy, tag in [(False, "1A_ratio_linear"), (True, "1B_ratio_log")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(R.shape[1]):
            (ax.semilogy if logy else ax.plot)(t, R[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlabel("time (step)"); ax.set_ylabel("sigma_j / sigma_1")
        ax.set_title(f"N={n} Fig{tag}: all rank branches sigma_j/sigma_1" + (" (log)" if logy else ""))
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # --- 図1C: ヒートマップ time×rank, log10(ratio) ---
    fig, ax = plt.subplots(figsize=(11, 5))
    L = np.log10(np.clip(R.T, 1e-20, None))
    im = ax.imshow(L, aspect="auto", origin="lower",
                   extent=[t[0], t[-1], -0.5, R.shape[1] - 0.5], cmap="viridis")
    ax.axvline(cr, color="w", ls=":", lw=0.8)
    ax.set_xlabel("time (step)"); ax.set_ylabel("rank index")
    ax.set_title(f"N={n} Fig1C: log10(sigma_j/sigma_1) heatmap")
    fig.colorbar(im, ax=ax, label="log10(sigma_j/sigma_1)")
    fig.tight_layout(); fig.savefig(fd / "fig1C_ratio_heatmap.png", dpi=130); plt.close(fig)

    # --- 図1D: branch追跡 ---
    tb, branches, Rb = pivot(ev, "branch_id", "sigma_over_sigma1")
    fig, ax = plt.subplots(figsize=(10, 6))
    for j in range(Rb.shape[1]):
        ax.semilogy(tb, Rb[:, j], lw=0.6)
    ax.axvline(cr, color="k", ls=":", lw=0.8)
    ax.set_xlabel("time (step)"); ax.set_ylabel("sigma / sigma_1 (branch-tracked)")
    ax.set_title(f"N={n} Fig1D: branch-tracked sigma/sigma_1 ({Rb.shape[1]} branches)")
    fig.tight_layout(); fig.savefig(fd / "fig1D_branch_ratio.png", dpi=130); plt.close(fig)

    # --- 図2A/2B: N・N² 倍 ---
    for col, tag, lbl in [("N_sigma_over_sigma1", "2A_N_ratio", "N * sigma_j/sigma_1"),
                          ("N2_sigma_over_sigma1", "2B_N2_ratio", "N^2 * sigma_j/sigma_1")]:
        t2, _, R2 = pivot(ev, "rank_index", col)
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(R2.shape[1]):
            ax.plot(t2, R2[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlabel("time (step)"); ax.set_ylabel(lbl)
        ax.set_title(f"N={n} Fig{tag}")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # --- 図3: 全 δ_j 線形/対数 ---
    td, _, D = pivot(raw / "delta.csv", "rank_index", "delta")
    for logy, tag in [(False, "3_delta_linear"), (True, "3_delta_log")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(D.shape[1]):
            (ax.semilogy if logy else ax.plot)(td, np.clip(D[:, j], 1e-20, None) if logy else D[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlabel("time (step)"); ax.set_ylabel("delta_j = sqrt(1 - overlap_with_parent)")
        ax.set_title(f"N={n} Fig{tag}")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # --- 図4: 残差行列 全特異値 線形/対数/N/N² ---
    ts, ks, S = pivot(raw / "residual_svd.csv", "k", "s_delta")
    for logy, tag in [(False, "4_sdelta_linear"), (True, "4_sdelta_log")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(S.shape[1]):
            (ax.semilogy if logy else ax.plot)(ts, np.clip(S[:, j], 1e-20, None) if logy else S[:, j], lw=0.5)
        ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlabel("time (step)"); ax.set_ylabel("s_k (residual matrix R_all singular values)")
        ax.set_title(f"N={n} Fig{tag}: all singular values of (I-Pi0)B_all")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)
    _, _, SN = pivot(raw / "residual_svd.csv", "k", "N_s_delta")
    _, _, SN2 = pivot(raw / "residual_svd.csv", "k", "N2_s_delta")
    for Mat, tag in [(SN, "4_sdelta_N"), (SN2, "4_sdelta_N2")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(Mat.shape[1]):
            ax.plot(ts, Mat[:, j], lw=0.5)
        ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlabel("time (step)"); ax.set_ylabel(tag)
        ax.set_title(f"N={n} Fig{tag}")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # --- 図5: q1..q4 ＋ q3,q4 拡大 ---
    qt = []; Q = []
    for row, _ in [load_rows(raw / "q_svd.csv")]:
        for r in row:
            qt.append(int(float(r["time"])))
            Q.append([float(r["q1"]), float(r["q2"]), float(r["q3"]), float(r["q4"])])
    qt = np.array(qt); Q = np.array(Q)
    fig, ax = plt.subplots(figsize=(10, 6))
    for j, lab in enumerate(["q1", "q2", "q3", "q4"]):
        ax.plot(qt, Q[:, j], lw=0.8, label=lab)
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.legend()
    ax.set_xlabel("time (step)"); ax.set_ylabel("singular values of Q=[B0|B1]")
    ax.set_title(f"N={n} Fig5: four-basis combined matrix singular values")
    fig.tight_layout(); fig.savefig(fd / "fig5_q_svd.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(qt, np.clip(Q[:, 2], 1e-20, None), lw=0.8, label="q3")
    ax.semilogy(qt, np.clip(Q[:, 3], 1e-20, None), lw=0.8, label="q4")
    ax.axvline(cr, color="k", ls=":", lw=0.8); ax.legend()
    ax.set_xlabel("time (step)"); ax.set_ylabel("q3, q4 (log)")
    ax.set_title(f"N={n} Fig5 zoom: q3, q4")
    fig.tight_layout(); fig.savefig(fd / "fig5_q34_zoom.png", dpi=130); plt.close(fig)

    # --- 図6: 全固有平面占有 E_j 線形/対数 ---
    te, _, E = pivot(raw / "occupation.csv", "rank_index", "occupation_fraction")
    for logy, tag in [(False, "6_occ_linear"), (True, "6_occ_log")]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in range(E.shape[1]):
            (ax.semilogy if logy else ax.plot)(te, np.clip(E[:, j], 1e-20, None) if logy else E[:, j], lw=0.6)
        ax.axvline(cr, color="k", ls=":", lw=0.8)
        ax.set_xlabel("time (step)"); ax.set_ylabel("E_j / |Z|^2")
        ax.set_title(f"N={n} Fig{tag}: occupation per eigenplane")
        fig.tight_layout(); fig.savefig(fd / f"fig{tag}.png", dpi=130); plt.close(fig)

    # --- 図7: σ/σ1・δ・E 対応（rank0,1,2 と最小占有枝を代表） ---
    _, _, RR = pivot(ev, "rank_index", "sigma_over_sigma1")
    _, _, DD = pivot(raw / "delta.csv", "rank_index", "delta")
    _, _, EE = pivot(raw / "occupation.csv", "rank_index", "occupation_fraction")
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for j in range(min(RR.shape[1], EE.shape[1])):
        axes[0].plot(t, RR[:, j], lw=0.6)
        axes[1].semilogy(td, np.clip(DD[:, j], 1e-20, None), lw=0.6)
        axes[2].semilogy(te, np.clip(EE[:, j], 1e-20, None), lw=0.6)
    for a in axes: a.axvline(cr, color="k", ls=":", lw=0.8)
    axes[0].set_ylabel("sigma_j/sigma_1"); axes[1].set_ylabel("delta_j (log)")
    axes[2].set_ylabel("E_j fraction (log)"); axes[2].set_xlabel("time (step)")
    axes[0].set_title(f"N={n} Fig7: sigma-ratio / delta / occupation per mode")
    fig.tight_layout(); fig.savefig(fd / "fig7_mode_correspondence.png", dpi=130); plt.close(fig)

    # --- 図8: 代表時刻スペクトル（順位 vs 比）線形/対数/N/N² ---
    import json
    rep = json.load(open(BASE / "diagnostics" / f"N{n:05d}.json"))["representative_times"]
    # 各代表時刻の全順位比を eigenvalues.csv から抽出
    rows, _ = load_rows(ev)
    byt = defaultdict(list)
    for r in rows:
        byt[int(float(r["time"]))].append((int(float(r["rank_index"])),
                                           float(r["sigma_over_sigma1"]),
                                           float(r["N_sigma_over_sigma1"]),
                                           float(r["N2_sigma_over_sigma1"])))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    panels = [(0, 0, 1, False, "sigma_j/sigma_1"), (0, 1, 1, True, "sigma_j/sigma_1 (log)"),
              (1, 0, 2, False, "N*sigma_j/sigma_1"), (1, 1, 3, False, "N^2*sigma_j/sigma_1")]
    for (rr, cc, col, logy, lbl) in panels:
        ax = axes[rr][cc]
        for label, tt in rep.items():
            data = sorted(byt.get(int(tt), []))
            if not data:
                continue
            xr = [d[0] for d in data]; yr = [d[col] for d in data]
            (ax.semilogy if logy else ax.plot)(xr, np.clip(yr, 1e-20, None) if logy else yr,
                                               marker="o", ms=3, lw=0.8, label=label)
        ax.set_xlabel("rank index"); ax.set_ylabel(lbl)
        if rr == 0 and cc == 0: ax.legend(fontsize=7)
    fig.suptitle(f"N={n} Fig8: representative-time spectra (rank vs ratio)")
    fig.tight_layout(); fig.savefig(fd / "fig8_representative_spectra.png", dpi=130); plt.close(fig)
    print(f"[figs] N={n}: 図1A〜8 を {fd} に出力")


def make_comparison(ns):
    cd = BASE / "figures" / "comparison"; cd.mkdir(parents=True, exist_ok=True)
    # 非支配最大 δ と 残差最大特異値の時間推移を N で重ねる
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for n in ns:
        raw = BASE / "raw" / f"N{n:05d}"; cr = crossing_of(n)
        ts, ks, S = pivot(raw / "residual_svd.csv", "k", "s_delta")
        axes[0].semilogy(ts - cr, np.clip(S[:, 0], 1e-20, None), lw=0.8, label=f"N={n}")
        _, _, D = pivot(raw / "delta.csv", "rank_index", "delta")
        # 支配(rank0)以外の最大δ
        axes[1].plot(ts - cr, np.nanmax(D[:, 1:], axis=1) if D.shape[1] > 1 else np.full(len(ts), np.nan),
                     lw=0.8, label=f"N={n}")
    axes[0].set_xlabel("time - crossing"); axes[0].set_ylabel("s_1 (residual matrix top singular value, log)")
    axes[0].set_title("Comparison: top residual singular value"); axes[0].legend()
    axes[1].set_xlabel("time - crossing"); axes[1].set_ylabel("max delta over non-dominant ranks")
    axes[1].set_title("Comparison: max non-dominant delta"); axes[1].legend()
    fig.tight_layout(); fig.savefig(cd / "comparison_N5_N40.png", dpi=130); plt.close(fig)
    print(f"[figs] comparison {ns} を {cd} に出力")


if __name__ == "__main__":
    if sys.argv[1] == "comparison":
        make_comparison([int(x) for x in sys.argv[2:]])
    else:
        for a in sys.argv[1:]:
            make_all(int(a))
