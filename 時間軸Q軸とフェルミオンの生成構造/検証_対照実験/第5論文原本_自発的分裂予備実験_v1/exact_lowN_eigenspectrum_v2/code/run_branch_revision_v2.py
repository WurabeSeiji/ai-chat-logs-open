#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 第1段階 v2 branch修正（後処理のみ）：N=5。解釈なし。N=40/N=300 未着手。

固有分解本体（eigh(iK)・全固有値・クラスタ）は維持（run_exact_lowN_eigenspectrum_v2.decompose）。
後処理を修正：
  §2 branch追跡を Hungarian 一対一（重複割当禁止・同時刻ID一意）。最小重なりは 0.5/0.7/0.9/0.99 診断、主 0.7。
  §3 分裂・合流は lineage edge（continuation/split_candidate/merge_candidate/birth/death/ambiguous）で別保存。
  §4 initial_floor は初期時刻部分空間との直接重なりで定義（旧・前時刻対応失敗＝床 を廃止）。
  §5 併合閾値 1e-10..1e-14 の感度表。
  §6 crossing=1167（真）と代表 step=1165（直前サンプル）を区別。真 crossing 表も追加。
出力 §7、図 §8（A〜F）。

使い方: python3 run_branch_revision_v2.py 5
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
BASE = CODE.parent
sys.path.insert(0, str(CODE))
from run_exact_lowN_eigenspectrum_v2 import (
    build_init, fval, dense_K, decompose, sample_schedule, EPS,
)

ACCEPT_PRIMARY = 0.7
ACCEPT_DIAG = [0.5, 0.7, 0.9, 0.99]
LINEAGE_TOL = 0.3
FLOOR_RHO_PRIMARY = 1e3
FLOOR_RHO_LABELS = [1e2, 1e3, 1e4, 1e6]
MERGE_SWEEP = [1e-10, 1e-11, 1e-12, 1e-13, 1e-14]
ORIGIN_THRESHOLDS = [0.90, 0.95, 0.99, 0.999]
ORIGIN_PRIMARY = 0.99
ORIGIN_CLOSURE_TOL = 1e-12
FMT = "%.17e"


def norm_overlap(Bi, Bj):
    """Tr(Pi Pj)/sqrt(rank_i rank_j) = ||Bi^T Bj||_F^2 / sqrt(dim_i dim_j) ∈ [0,1]."""
    return float(np.sum((Bi.T @ Bj) ** 2) / np.sqrt(Bi.shape[1] * Bj.shape[1]))


def overlap_space(Ucols, B):
    """O(U,B)=||U† B||_F² / dim(B) = Tr(Π_U Π_B)/rank(Π_B)。クラスタ基底の**全列**を使う。"""
    if B.shape[1] == 0 or Ucols.shape[1] == 0:
        return 0.0
    Om = Ucols.conj().T @ B.astype(complex)
    return float(np.real_if_close(np.linalg.norm(Om, "fro") ** 2 / B.shape[1]))


def origin_status(of, onf, sumv, thr, closure_tol=ORIGIN_CLOSURE_TOL):
    if abs(sumv - 1.0) > closure_tol:
        return "undetermined"
    if of >= thr:
        return "initial_floor"
    if onf >= thr:
        return "initial_nonfloor"
    if of > 0.01 and onf > 0.01:
        return "mixed"
    return "undetermined"


def collect(n):
    """軌道を走らせ、サンプル時刻の Z を保存。真 crossing=1167 も明示追加。"""
    sys_lr, v, B0, p, q, Z, wp = build_init(n)
    # crossing 検出
    Zc = Z.copy(); wpc = wp.copy(); crossing = None; t = 0
    while True:
        if fval(Zc, p, q) > 0.05:
            crossing = t; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc)
        Zc = sys_lr.cayley_step(Zc, se); t += 1
    sched, end = sample_schedule(crossing)
    sched.add(crossing)                 # 真 crossing step を追加
    # 軌道保存
    Zs = {}
    Zr = Z.copy(); wpr = wp.copy(); t = 0
    while True:
        if t in sched:
            Zs[t] = Zr.copy()
        if t >= end:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr)
        Zr = sys_lr.cayley_step(Zr, se); t += 1
    return sys_lr, B0, p, q, Zs, crossing, end


def main(n=5):
    assert n in (5, 40), "N=5 検収済み。N=40 は実行指示により許可。N=300 未着手。"
    raw = BASE / "raw" / f"N{n:05d}"; fd = BASE / "figures" / f"N{n:05d}"
    raw.mkdir(parents=True, exist_ok=True); fd.mkdir(parents=True, exist_ok=True)
    sys_lr, B0, p, q, Zs, crossing, end = collect(n)
    times = sorted(Zs)
    rep_crossing = min((x for x in times if x <= crossing), key=lambda x: abs(x - crossing))

    # 主分解（merge_tol=1e-12）を全時刻で
    dec = {}
    for t in times:
        K = dense_K(sys_lr, Zs[t])
        dec[t] = decompose(K, Zs[t], B0, merge_tol=1e-12)

    # 初期時刻 t0 の床/非床部分空間（複素固有ベクトル）
    t0 = times[0]
    d0 = dec[t0]
    rho0 = d0["floor_ratio"]
    floor_idx0 = {lbl: np.where(rho0 < lbl)[0] for lbl in FLOOR_RHO_LABELS}
    U0 = d0["U"]
    Uf = U0[:, floor_idx0[FLOOR_RHO_PRIMARY]]
    Unf = U0[:, np.setdiff1d(np.arange(U0.shape[1]), floor_idx0[FLOOR_RHO_PRIMARY])]
    # 代表時刻（summary 用, 併合感度と共有）
    reps = {"initial": t0, "nearest_before_crossing_step%d" % rep_crossing: rep_crossing,
            "true_crossing_step%d" % crossing: crossing,
            "post_crossing": min(times, key=lambda x: abs(x - (crossing + 200))),
            "plateau_start": min(times, key=lambda x: abs(x - (crossing + 5000))), "final": end}

    # ---- Hungarian 一対一 branch 追跡 ----
    f_bij = open(raw / "branch_tracking_bijective.csv", "w", newline=""); w_bij = csv.writer(f_bij)
    w_bij.writerow(["source_step", "target_step", "source_cluster_id", "target_cluster_id",
                    "source_branch_id", "target_branch_id", "overlap", "assignment_cost",
                    "accepted", "tracking_status"])
    f_lin = open(raw / "cluster_lineage.csv", "w", newline=""); w_lin = csv.writer(f_lin)
    w_lin.writerow(["source_step", "target_step", "source_cluster_id", "target_cluster_id",
                    "source_branch_id", "target_branch_id", "overlap", "relation_type"])
    f_org = open(raw / "initial_space_origin.csv", "w", newline=""); w_org = csv.writer(f_org)
    w_org.writerow(["step", "time", "cluster_id", "branch_id", "cluster_dimension",
                    "sigma_min", "sigma_max", "sigma_representative",
                    "overlap_with_initial_floor_space", "overlap_with_initial_nonfloor_space",
                    "origin_overlap_sum", "origin_overlap_closure_error",
                    "initial_origin_status", "origin_threshold", "floor_threshold_label"])

    branch_of = {}          # t -> {cid: branch_id}
    next_branch = [0]
    dup_count = 0
    unmatched_total = 0
    ambiguous_total = 0
    split_total = 0
    merge_total = 0
    min_overlap_series = []

    # t0 初期化
    branch_of[t0] = {}
    for c in dec[t0]["clusters"]:
        branch_of[t0][c["cid"]] = next_branch[0]; next_branch[0] += 1

    for a, b in zip(times[:-1], times[1:]):
        Pa = dec[a]["clusters"]; Pb = dec[b]["clusters"]
        na, nb = len(Pa), len(Pb)
        C = np.zeros((na, nb))
        for i in range(na):
            for j in range(nb):
                C[i, j] = norm_overlap(Pa[i]["B"], Pb[j]["B"])
        # Hungarian（最大重なり＝ -C の最小コスト）
        ri, ci = linear_sum_assignment(-C)
        matched_target = {}
        assigned = {}
        for i, j in zip(ri, ci):
            ov = C[i, j]
            acc = ov >= ACCEPT_PRIMARY
            assigned[j] = (i, ov, acc)
        branch_of[b] = {}
        min_ov_accepted = 1.0
        for j in range(nb):
            if j in assigned and assigned[j][2]:
                i, ov, acc = assigned[j]
                src_branch = branch_of[a][Pa[i]["cid"]]
                branch_of[b][Pb[j]["cid"]] = src_branch
                status = "continuation"
                min_ov_accepted = min(min_ov_accepted, ov)
                w_bij.writerow([a, b, Pa[i]["cid"], Pb[j]["cid"], src_branch, src_branch,
                                FMT % ov, FMT % (1 - ov), 1, status])
            else:
                branch_of[b][Pb[j]["cid"]] = next_branch[0]
                nb_id = next_branch[0]; next_branch[0] += 1
                ov = assigned[j][1] if j in assigned else 0.0
                unmatched_total += 1
                w_bij.writerow([a, b, (Pa[assigned[j][0]]["cid"] if j in assigned else -1), Pb[j]["cid"],
                                (branch_of[a][Pa[assigned[j][0]]["cid"]] if j in assigned else -1), nb_id,
                                FMT % ov, FMT % (1 - ov), 0, "unmatched"])
        min_overlap_series.append((b, min_ov_accepted if min_ov_accepted <= 1.0 else np.nan))
        # 同時刻 branch ID 一意性検査
        ids = list(branch_of[b].values())
        assert len(ids) == len(set(ids)), f"branch重複 at step {b}"
        if len(ids) != len(set(ids)):
            dup_count += 1

        # ---- lineage（閾値以上の全対応）----
        src_deg = defaultdict(int); tgt_deg = defaultdict(int)
        edges = []
        for i in range(na):
            for j in range(nb):
                if C[i, j] >= LINEAGE_TOL:
                    edges.append((i, j, C[i, j])); src_deg[i] += 1; tgt_deg[j] += 1
        matched_src = set(i for i, _, _ in edges); matched_tgt = set(j for _, j, _ in edges)
        for (i, j, ov) in edges:
            if src_deg[i] > 1 and tgt_deg[j] > 1:
                rel = "ambiguous"; ambiguous_total += 1
            elif src_deg[i] > 1:
                rel = "split_candidate"; split_total += 1
            elif tgt_deg[j] > 1:
                rel = "merge_candidate"; merge_total += 1
            else:
                rel = "continuation"
            w_lin.writerow([a, b, Pa[i]["cid"], Pb[j]["cid"], branch_of[a][Pa[i]["cid"]],
                            branch_of[b][Pb[j]["cid"]], FMT % ov, rel])
        for j in range(nb):
            if j not in matched_tgt:
                w_lin.writerow([a, b, -1, Pb[j]["cid"], -1, branch_of[b][Pb[j]["cid"]], FMT % 0.0, "birth"])
        for i in range(na):
            if i not in matched_src:
                w_lin.writerow([a, b, Pa[i]["cid"], -1, branch_of[a][Pa[i]["cid"]], -1, FMT % 0.0, "death"])

    # ---- initial_space_origin（全時刻・全クラスタ, 全列 Frobenius）----
    origin_rows = []           # (t, cid, branch, dim, of, onf, sumv, err, status_primary)
    closure_errs = []
    status_counts = {thr: defaultdict(int) for thr in ORIGIN_THRESHOLDS}
    for t in times:
        for c in dec[t]["clusters"]:
            of = overlap_space(Uf, c["B"]); onf = overlap_space(Unf, c["B"])
            sumv = of + onf; err = abs(sumv - 1.0)
            closure_errs.append(err)
            st_primary = origin_status(of, onf, sumv, ORIGIN_PRIMARY)
            for thr in ORIGIN_THRESHOLDS:
                status_counts[thr][origin_status(of, onf, sumv, thr)] += 1
            origin_rows.append((t, c["cid"], branch_of[t][c["cid"]], c["dim"], of, onf, sumv, err, st_primary))
            w_org.writerow([t, t, c["cid"], branch_of[t][c["cid"]], c["dim"], FMT % c["sigma_min"],
                            FMT % c["sigma_max"], FMT % c["sigma"], FMT % of, FMT % onf,
                            FMT % sumv, FMT % err, st_primary, ORIGIN_PRIMARY, f"rho<{FLOOR_RHO_PRIMARY:.0e}"])
    for fh in (f_bij, f_lin, f_org):
        fh.close()

    # origin 診断 JSON
    ce = np.array(closure_errs)
    prim = status_counts[ORIGIN_PRIMARY]
    origin_diag = {
        "N": n, "overlap_definition": "||U_dagger B||_F^2 / dim(B) (all columns)",
        "floor_rho_primary": FLOOR_RHO_PRIMARY, "origin_primary_threshold": ORIGIN_PRIMARY,
        "origin_closure_tol": ORIGIN_CLOSURE_TOL,
        "max_origin_overlap_closure_error": float(ce.max()),
        "median_origin_overlap_closure_error": float(np.median(ce)),
        "count_origin_overlap_error_gt_1e-12": int(np.sum(ce > 1e-12)),
        "count_origin_overlap_error_gt_1e-10": int(np.sum(ce > 1e-10)),
        "closure_pass_lt_1e-12": bool(ce.max() < 1e-12),
        "initial_floor_count": prim["initial_floor"], "initial_nonfloor_count": prim["initial_nonfloor"],
        "mixed_count": prim["mixed"], "undetermined_count": prim["undetermined"],
        "status_counts_by_origin_threshold": {str(thr): dict(status_counts[thr]) for thr in ORIGIN_THRESHOLDS},
    }
    with open(BASE / "diagnostics" / f"N{n:05d}_initial_origin_revision.json", "w", encoding="utf-8") as fh:
        json.dump(origin_diag, fh, indent=2, ensure_ascii=False)

    # origin 集計表（代表時刻）
    tab = BASE / "tables" / f"N{n:05d}"; tab.mkdir(parents=True, exist_ok=True)
    with open(tab / "initial_origin_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "label", "cluster_count", "initial_floor_count", "initial_nonfloor_count",
                    "mixed_count", "undetermined_count", "max_origin_overlap_closure_error"])
        for label, t in reps.items():
            rr = [r for r in origin_rows if r[0] == t]
            cnt = defaultdict(int)
            for r in rr:
                cnt[r[8]] += 1
            mxe = max((r[7] for r in rr), default=0.0)
            w.writerow([t, label, len(rr), cnt["initial_floor"], cnt["initial_nonfloor"],
                        cnt["mixed"], cnt["undetermined"], FMT % mxe])

    # ---- 併合閾値感度（代表時刻, reps は上で定義）----
    f_sen = open(raw / "merge_tolerance_sensitivity.csv", "w", newline=""); w_sen = csv.writer(f_sen)
    w_sen.writerow(["step", "label", "merge_tolerance", "cluster_count", "cluster_dimensions",
                    "dominant_cluster_dimension", "dominant_delta", "dominant_occupation",
                    "q1", "q2", "q3", "q4", "max_intercluster_overlap", "closure_error"])
    max_inter_by_tol = defaultdict(float); closure_by_tol = defaultdict(float)
    for label, t in reps.items():
        K = dense_K(sys_lr, Zs[t])
        for mt in MERGE_SWEEP:
            d = decompose(K, Zs[t], B0, merge_tol=mt)
            cls = d["clusters"]
            dom = max(cls, key=lambda c: c["sigma"])
            Bdom = dom["B"][:, :2]
            Q = np.column_stack([B0, Bdom]); qs = np.linalg.svd(Q, compute_uv=False)
            qs = np.pad(qs, (0, max(0, 4 - len(qs))))[:4]
            dims = "|".join(str(c["dim"]) for c in cls)
            w_sen.writerow([t, label, FMT % mt, len(cls), dims, dom["dim"], FMT % dom["delta"],
                            FMT % (dom["occ"] / 1.0), FMT % qs[0], FMT % qs[1], FMT % qs[2], FMT % qs[3],
                            FMT % d["max_cluster_inter"], FMT % d["closure"]])
            max_inter_by_tol[mt] = max(max_inter_by_tol[mt], d["max_cluster_inter"])
            closure_by_tol[mt] = max(closure_by_tol[mt], d["closure"])
    f_sen.close()

    # ---- 診断 JSON ----
    diagj = {
        "N": n, "crossing_step": crossing, "representative_crossing_step": rep_crossing,
        "note_crossing": (f"実際のcrossingはstep={crossing}。v2主表(run_exact_lowN_eigenspectrum_v2)は"
                          f"5step格子のため直前点step=1165を'crossing'表記に使用していた。本branch修正では"
                          f"真crossing step={crossing}をサンプルに追加し、代表crossing step={rep_crossing}"
                          f"（=真crossing）で後処理する。"),
        "acceptance_primary": ACCEPT_PRIMARY, "acceptance_diagnostic_thresholds": ACCEPT_DIAG,
        "lineage_overlap_threshold": LINEAGE_TOL, "floor_rho_primary": FLOOR_RHO_PRIMARY,
        "duplicate_branch_id_count": dup_count, "unmatched_cluster_count": unmatched_total,
        "ambiguous_assignment_count": ambiguous_total, "split_candidate_count": split_total,
        "merge_candidate_count": merge_total,
        "max_intercluster_overlap_by_tolerance": {("%.0e" % k): v for k, v in max_inter_by_tol.items()},
        "closure_error_by_tolerance": {("%.0e" % k): v for k, v in closure_by_tol.items()},
        "final_spectrum_sigma_over_sigma1": sorted(
            [c["sigma"] / dec[end]["sig1"] for c in dec[end]["clusters"]], reverse=True),
        "n_branches_total": next_branch[0],
    }
    (BASE / "diagnostics").mkdir(exist_ok=True)
    with open(BASE / "diagnostics" / f"N{n:05d}_branch_revision.json", "w", encoding="utf-8") as fh:
        json.dump(diagj, fh, indent=2, ensure_ascii=False)

    make_figures(n, dec, times, crossing, branch_of, raw, fd, min_overlap_series, reps)
    print(f"[branch修正 N={n}] crossing真={crossing} 代表={rep_crossing}")
    print(f"  一対一: branch重複={dup_count}(assert通過) unmatched={unmatched_total} "
          f"ambiguous={ambiguous_total} split={split_total} merge={merge_total} 総branch数={next_branch[0]}")
    print(f"  併合感度 max_intercluster_overlap: " +
          " ".join(f"{k:.0e}:{v:.1e}" for k, v in max_inter_by_tol.items()))
    print(f"  最終スペクトル σ/σ1: " + " ".join(f"{x:.6f}" for x in diagj['final_spectrum_sigma_over_sigma1']))
    return diagj


def make_figures(n, dec, times, crossing, branch_of, raw, fd, min_overlap_series, reps):
    # 図A：一対一 branch 固有値推移
    brt = defaultdict(dict)
    for t in times:
        sig1 = dec[t]["sig1"]
        for c in dec[t]["clusters"]:
            brt[branch_of[t][c["cid"]]][t] = c["sigma"] / sig1
    fig, ax = plt.subplots(figsize=(10, 6))
    for b in brt:
        xs = sorted(brt[b]); ax.plot(xs, [brt[b][x] for x in xs], lw=0.7)
    ax.axvline(crossing, color="k", ls=":", lw=0.8)
    ax.set_xlabel("time"); ax.set_ylabel("sigma_branch/sigma_1")
    ax.set_title(f"N={n} FigA: bijective branch sigma/sigma1 ({len(brt)} branches, no duplicate id/time)")
    fig.tight_layout(); fig.savefig(fd / "figA_bijective_branch_ratio.png", dpi=130); plt.close(fig)

    # 図B：lineage（split/merge/birth/death を時刻ごとに件数）
    lin = list(csv.DictReader(open(raw / "cluster_lineage.csv")))
    byt = defaultdict(lambda: defaultdict(int))
    for r in lin:
        byt[int(float(r["target_step"]))][r["relation_type"]] += 1
    ts = sorted(byt)
    fig, ax = plt.subplots(figsize=(10, 6))
    for rel in ["split_candidate", "merge_candidate", "birth", "death", "ambiguous"]:
        ax.plot(ts, [byt[x].get(rel, 0) for x in ts], lw=0.8, label=rel)
    ax.axvline(crossing, color="k", ls=":", lw=0.8); ax.legend(fontsize=8)
    ax.set_xlabel("time"); ax.set_ylabel("edge count"); ax.set_title(f"N={n} FigB: lineage relation counts")
    fig.tight_layout(); fig.savefig(fd / "figB_lineage_counts.png", dpi=130); plt.close(fig)

    # 図C：初期床/非床空間重なり（branch別, 全列 Frobenius）＋ 閉鎖 O_floor+O_nonfloor
    org = list(csv.DictReader(open(raw / "initial_space_origin.csv")))
    of_b = defaultdict(dict); onf_b = defaultdict(dict); sum_b = defaultdict(dict)
    for r in org:
        b = int(r["branch_id"]); t = int(float(r["step"]))
        of_b[b][t] = float(r["overlap_with_initial_floor_space"])
        onf_b[b][t] = float(r["overlap_with_initial_nonfloor_space"])
        sum_b[b][t] = float(r["origin_overlap_sum"])
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for b in of_b:
        xs = sorted(of_b[b]); a1.plot(xs, [of_b[b][x] for x in xs], lw=0.6)
        a2.plot(xs, [onf_b[b][x] for x in xs], lw=0.6)
    for a in (a1, a2): a.axvline(crossing, color="k", ls=":", lw=0.8)
    a1.set_ylabel("overlap w/ initial FLOOR space"); a2.set_ylabel("overlap w/ initial NONFLOOR space")
    a2.set_xlabel("time"); a1.set_title(f"N={n} FigC: overlap with initial-time subspaces (all-column Frobenius, per branch)")
    fig.tight_layout(); fig.savefig(fd / "figC_initial_space_overlap.png", dpi=130); plt.close(fig)
    # 図C 閉鎖: O_floor + O_nonfloor（≈1 の確認）
    fig, ax = plt.subplots(figsize=(10, 5))
    for b in sum_b:
        xs = sorted(sum_b[b]); ax.plot(xs, [sum_b[b][x] for x in xs], lw=0.6)
    ax.axhline(1.0, color="r", ls="--", lw=0.8); ax.axvline(crossing, color="k", ls=":", lw=0.8)
    ax.set_xlabel("time"); ax.set_ylabel("O_floor + O_nonfloor"); ax.set_ylim(0.5, 1.5)
    ax.set_title(f"N={n} FigC(closure): O_floor+O_nonfloor (should be 1)")
    fig.tight_layout(); fig.savefig(fd / "figC_origin_overlap_closure.png", dpi=130); plt.close(fig)

    # 図D：修正後 O_floor が閾値以上の branch の σ/σ1（旧 initial_floor_flag は不使用）
    ORIGIN_THR = 0.99
    floor_branches = set(int(r["branch_id"]) for r in org
                         if float(r["overlap_with_initial_floor_space"]) >= ORIGIN_THR)
    fig, ax = plt.subplots(figsize=(10, 6))
    for b in (floor_branches or []):
        xs = sorted(brt[b]); ax.plot(xs, [brt[b][x] for x in xs], lw=0.7)
    ax.axvline(crossing, color="k", ls=":", lw=0.8); ax.set_xlabel("time"); ax.set_ylabel("sigma/sigma1")
    ax.set_title(f"N={n} FigD: branches with initial-FLOOR overlap>={ORIGIN_THR} "
                 f"(all-column, {len(floor_branches)} branches)")
    fig.tight_layout(); fig.savefig(fd / "figD_initial_floor_branches.png", dpi=130); plt.close(fig)

    # 図E：併合閾値感度（代表時刻ごとに cluster数, dominant δ, dominant占有, q3, q4）
    sen = list(csv.DictReader(open(raw / "merge_tolerance_sensitivity.csv")))
    labels = sorted(set(r["label"] for r in sen))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    metrics = [("cluster_count", "cluster count"), ("dominant_delta", "dominant delta"),
               ("dominant_occupation", "dominant occupation"), ("q3", "q3"), ("q4", "q4"),
               ("max_intercluster_overlap", "max intercluster overlap")]
    for ax, (col, ttl) in zip(axes.flat, metrics):
        for lb in labels:
            rr = [r for r in sen if r["label"] == lb]
            xs = [float(r["merge_tolerance"]) for r in rr]; ys = [float(r[col]) for r in rr]
            o = np.argsort(xs)
            ax.plot(np.array(xs)[o], np.array(ys)[o], marker="o", ms=3, lw=0.8, label=lb)
        ax.set_xscale("log")
        if col in ("dominant_occupation", "max_intercluster_overlap"):
            ax.set_yscale("log")
        ax.set_xlabel("merge_tolerance"); ax.set_ylabel(ttl); ax.set_title(ttl)
    axes.flat[0].legend(fontsize=6)
    fig.suptitle(f"N={n} FigE: merge-tolerance sensitivity")
    fig.tight_layout(); fig.savefig(fd / "figE_merge_sensitivity.png", dpi=130); plt.close(fig)

    # 図F：採用 branch 対応の最小重なり
    mo = np.array([(t, v) for t, v in min_overlap_series if not np.isnan(v)])
    fig, ax = plt.subplots(figsize=(10, 6))
    if len(mo):
        ax.plot(mo[:, 0], mo[:, 1], lw=0.8)
    for thr in ACCEPT_DIAG:
        ax.axhline(thr, ls="--", lw=0.6, color="gray")
    ax.axvline(crossing, color="k", ls=":", lw=0.8); ax.set_xlabel("time")
    ax.set_ylabel("min accepted branch overlap"); ax.set_ylim(0, 1.02)
    ax.set_title(f"N={n} FigF: minimum accepted branch-correspondence overlap")
    fig.tight_layout(); fig.savefig(fd / "figF_min_branch_overlap.png", dpi=130); plt.close(fig)
    print(f"[branch修正 figs] A〜F を {fd} に出力")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
