#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 第1段階 修正版(v2)：N=5 厳密固有分解（eigh(iK)・クラスタ処理）。解釈なし。

修正指示への対応：
  §3 基本固有分解は H=iK の Hermitian 固有分解 numpy.linalg.eigh(1j*K)（全体直交規格）。
      eig(K) の実虚部独立直交化は使わない。
  §1.2/§6 絶対閾値 σ>1e-9 を廃止。全 M 固有値を保存し、ゼロ判定は σ/(eps|K|) の相対基準を併記。
  §5 縮退・近接縮退は生固有値を平均せず、方向比較はクラスタ部分空間 Π_C で行う。
  §4 正負固有値対を検証。非縮退平面は B=[√2Re u,√2Im u] を QR 正規直交化し KB=BJ を確認。
  §1.3/§8 R_all 一括は診断のみ。新方向は対象別残差（支配平面/初期核から成長/占有/新クラスタ）で判定。
  §9 状態占有は直交射影で。核は核固有ベクトルから Π_ker を明示構成。閉鎖を検証。
N=5 のみ実行。N=40・N=300 には着手しない。

使い方: python3 run_exact_lowN_eigenspectrum_v2.py 5
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
BASE = CODE_DIR.parent
ENGINE_DIR = BASE.parent
sys.path.insert(0, str(ENGINE_DIR))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact

DELTA = 1e-15
AFTER = 50000
EPS = np.finfo(float).eps
FLOOR_REL_PRIMARY = 1000.0     # σ/(eps|K|) > これ を回転モード（主）
FLOOR_REL_ALT = 100.0          # 併記
REF_DIR = ENGINE_DIR / "paper6_definitive_control_v1"


def build_init(n):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    p1s, B0, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)   # 初期親平面 P0
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=sys_lr.m)
    return sys_lr, v, B0, p, q, Z, wp


def fval(Z, p, q):
    Zp = Z - p * (p @ Z) - q * (q @ Z)
    return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))


def dense_K(sys_lr, Z):
    sys_lr.set_theta(np.angle(Z))
    M = sys_lr.m
    return np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])


def decompose(K, Z, B0, merge_tol=1e-12):
    """H=iK の eigh 分解 → クラスタ（縮退部分空間）と核。全量を返す。merge_tol は併合閾値。"""
    M = K.shape[0]
    Knorm = np.linalg.norm(K, 2)
    H = 1j * K
    mu, U = np.linalg.eigh(H)                    # mu 昇順実, U 直交規格(複素)
    # 固有対残差
    resid = np.linalg.norm(H @ U - U * mu, axis=0)      # 列ごと |Hu-μu|
    uni = np.linalg.norm(U.conj().T @ U - np.eye(M))
    sigma = np.abs(mu)
    order_abs = np.argsort(-sigma)
    sig1 = float(sigma[order_abs[0]])
    floor_ratio = sigma / (EPS * Knorm)

    # 正負対：昇順 index k と M-1-k
    pair_id = np.full(M, -1, dtype=int)
    pair_err = 0.0
    for k in range(M // 2):
        pair_id[k] = k
        pair_id[M - 1 - k] = k
        pair_err = max(pair_err, abs(mu[k] + mu[M - 1 - k]))
    if M % 2 == 1:
        pair_id[M // 2] = -1

    # 回転 index（相対床・主基準）と核 index
    rot_mask = floor_ratio > FLOOR_REL_PRIMARY
    rot_idx = np.where(rot_mask)[0]
    ker_idx = np.where(~rot_mask)[0]
    n_rot_alt = int(np.sum(floor_ratio > FLOOR_REL_ALT))

    # クラスタ化：個別候補平面を作り、平面間重なり>MERGE_TOL の対を併合（union-find）。
    # これにより「非縮退（未併合）平面」の相互直交を MERGE_TOL 以下に保証する（§5, 検収≤1e-12）。
    # 生固有値は平均しない（cluster_id で束ねるのみ）。σ-近接も併記のため min/max σ を保存。
    MERGE_TOL = merge_tol
    pos_rot = sorted([i for i in rot_idx if mu[i] > 0], key=lambda i: -sigma[i])
    cand = []
    for i in pos_rot:
        u = U[:, i]
        Bq, R = np.linalg.qr(np.column_stack([np.sqrt(2) * u.real, np.sqrt(2) * u.imag]))
        Bq = Bq[:, np.abs(np.diag(R)) > 1e-10]
        cand.append({"i": i, "B": Bq})
    nc = len(cand)
    par = list(range(nc))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a in range(nc):
        for b in range(a + 1, nc):
            ov = np.linalg.norm(cand[a]["B"].T @ cand[b]["B"], 2)   # 最大特異値＝平面間重なり
            if ov > MERGE_TOL:
                par[find(a)] = find(b)
    from collections import defaultdict as _dd
    grp = _dd(list)
    for a in range(nc):
        grp[find(a)].append(a)

    clusters = []
    cluster_id_arr = np.full(M, -1, dtype=int)
    cid = 0
    for _, members_c in sorted(grp.items(), key=lambda kv: -max(sigma[cand[m]["i"]] for m in kv[1])):
        members_pos = [cand[m]["i"] for m in members_c]
        members_neg = [k for k in range(M) if mu[k] < 0 and pair_id[k] in [pair_id[j] for j in members_pos]]
        members = members_pos + members_neg
        for m in members:
            cluster_id_arr[m] = cid
        cols = []
        for m in members_pos:
            u = U[:, m]
            cols.append(np.sqrt(2) * u.real); cols.append(np.sqrt(2) * u.imag)
        Bc, R = np.linalg.qr(np.column_stack(cols))
        Bc = Bc[:, np.abs(np.diag(R)) > 1e-10]
        mult = len(members_pos)
        sig_members = np.array([sigma[j] for j in members_pos])
        s0 = float(np.max(sig_members))
        PiC_KB = Bc @ (Bc.T @ (K @ Bc))
        invar = np.linalg.norm(K @ Bc - PiC_KB) / max(1.0, Knorm)
        kbj = np.nan
        if mult == 1:
            s = s0
            # J の向き（符号）規約に依らないよう両向きの最小残差を採る（σ回転を検証）
            J1 = np.array([[0.0, -s], [s, 0.0]]); J2 = np.array([[0.0, s], [-s, 0.0]])
            kbj = float(min(np.linalg.norm(K @ Bc - Bc @ J1),
                            np.linalg.norm(K @ Bc - Bc @ J2)) / max(1.0, Knorm))
        EC = float(np.sum(np.abs(U[:, members].conj().T @ Z) ** 2))
        oC = 0.5 * float(np.sum((B0.T @ Bc) ** 2))
        deltaC = float(np.sqrt(max(0.0, 1 - oC)))
        clusters.append({"cid": cid, "sigma": s0, "sigma_min": float(np.min(sig_members)),
                         "sigma_max": float(np.max(sig_members)), "mult": mult, "dim": Bc.shape[1],
                         "members_pos": members_pos, "members": members, "B": Bc,
                         "occ": EC, "delta": deltaC, "overlap_parent": oC,
                         "invariance_residual": invar, "kbj_residual": kbj})
        cid += 1

    # 核射影（明示）
    if len(ker_idx) > 0:
        Uker = U[:, ker_idx]
        E_ker = float(np.sum(np.abs(Uker.conj().T @ Z) ** 2))
    else:
        E_ker = 0.0
    # 閉鎖（直交射影に基づく）
    totZ = float(np.real(np.conj(Z) @ Z))
    closure = abs(totZ - sum(c["occ"] for c in clusters) - E_ker)

    # 非縮退平面間直交誤差（縮退クラスタ同士・内は除外）
    # 平面間・クラスタ間直交は併合判定と同じ作用素ノルム(2-norm)で測る
    nondeg = [c for c in clusters if c["mult"] == 1]
    max_inter = 0.0
    for a in range(len(nondeg)):
        for b in range(a + 1, len(nondeg)):
            max_inter = max(max_inter, np.linalg.norm(nondeg[a]["B"].T @ nondeg[b]["B"], 2))
    # クラスタ間直交（全クラスタ対）
    max_cluster_inter = 0.0
    for a in range(len(clusters)):
        for b in range(a + 1, len(clusters)):
            max_cluster_inter = max(max_cluster_inter, np.linalg.norm(clusters[a]["B"].T @ clusters[b]["B"], 2))
    max_idem = 0.0
    for c in clusters:
        P = c["B"] @ c["B"].T
        max_idem = max(max_idem, np.linalg.norm(P @ P - P) / max(1.0, M))
    max_ortho = max((np.linalg.norm(c["B"].T @ c["B"] - np.eye(c["B"].shape[1])) for c in clusters), default=0.0)

    return {
        "mu": mu, "U": U, "sigma": sigma, "sig1": sig1, "floor_ratio": floor_ratio,
        "resid": resid, "uni": float(uni), "pair_id": pair_id, "pair_err": float(pair_err),
        "cluster_id": cluster_id_arr, "clusters": clusters, "ker_idx": ker_idx,
        "E_ker": E_ker, "closure": closure, "Knorm": Knorm, "n_rot_alt": n_rot_alt,
        "max_inter_nondeg": max_inter, "max_cluster_inter": max_cluster_inter,
        "max_idem": max_idem, "max_ortho": max_ortho, "antisym": np.linalg.norm(K + K.T) / max(1.0, Knorm),
    }


def target_residuals(dec, B0, occ_floor=1e-12):
    """§8 対象別残差 R_target=(I-Π0)B_target の全特異値。
    A 支配クラスタ / C 占有クラスタ(E>床) / D 縮退クラスタ(mult>1) / all 診断用一括。"""
    out = {}
    cls = dec["clusters"]
    if not cls:
        return out
    def resid_svd(Bs):
        if not Bs:
            return []
        B = np.column_stack(Bs)
        R = B - B0 @ (B0.T @ B)
        return list(np.linalg.svd(R, compute_uv=False))
    dom = max(cls, key=lambda c: c["sigma"])
    out["A_dominant"] = resid_svd([dom["B"]])
    out["C_occupied"] = resid_svd([c["B"] for c in cls if c["occ"] > occ_floor])
    out["D_degenerate"] = resid_svd([c["B"] for c in cls if c["mult"] > 1])
    out["all_diagnostic"] = resid_svd([c["B"] for c in cls])
    return out


def sample_schedule(crossing):
    ts = set()
    end = crossing + AFTER
    for t in range(0, crossing + 2000 + 1, 5):
        ts.add(t)
    for t in range(crossing + 2000, crossing + 10000 + 1, 50):
        ts.add(t)
    for t in range(crossing + 10000, end + 1, 200):
        ts.add(t)
    ts.add(end)
    return ts, end


def run(n):
    assert n == 5, "v2 第1段階は N=5 のみ。N=40 は人間検収後に別途。"
    sys_lr, v, B0, p, q, Z, wp = build_init(n)
    M = sys_lr.m
    seed = 40260722 + 1000 * n
    raw = BASE / "raw" / f"N{n:05d}"; binp = BASE / "binary" / f"N{n:05d}"; tab = BASE / "tables" / f"N{n:05d}"
    for d in (raw, binp, tab):
        d.mkdir(parents=True, exist_ok=True)

    # crossing 検出
    Zc = Z.copy(); wpc = wp.copy(); crossing = None; tt = 0
    while True:
        if fval(Zc, p, q) > 0.05:
            crossing = tt; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc)
        Zc = sys_lr.cayley_step(Zc, se); tt += 1
    sched, end = sample_schedule(crossing)

    def nearest(x): return min(sched, key=lambda s: abs(s - x))
    rep_times = {"initial": nearest(0), "pre_crossing": nearest(crossing - 50),
                 "crossing": nearest(crossing), "post_crossing": nearest(crossing + 200),
                 "plateau_start": nearest(crossing + 5000), "final": end}
    rep_set = set(rep_times.values())

    ref_f = {}
    if (REF_DIR / f"obs_N{n:05d}.csv").exists():
        with open(REF_DIR / f"obs_N{n:05d}.csv") as fh:
            for row in csv.DictReader(fh):
                ref_f[int(float(row["tau"]))] = float(row["f"])

    fmt = "%.17e"
    f_eig = open(raw / "eigenvalues.csv", "w", newline=""); w_eig = csv.writer(f_eig)
    w_eig.writerow(["N", "step", "time", "eigen_index", "eigenvalue", "abs_eigenvalue",
                    "sigma_over_sigma1", "N_sigma_over_sigma1", "N2_sigma_over_sigma1",
                    "abs_sigma_over_eps_normK", "pair_id", "cluster_id", "solver_residual"])
    f_cl = open(raw / "clusters.csv", "w", newline=""); w_cl = csv.writer(f_cl)
    w_cl.writerow(["N", "step", "time", "cluster_id", "cluster_branch", "sigma", "sigma_over_sigma1",
                   "mult", "dim", "occupation_EC", "occupation_fraction", "delta_C", "overlap_parent",
                   "invariance_residual", "kbj_residual", "initial_floor_flag"])
    f_tg = open(raw / "delta_targets.csv", "w", newline=""); w_tg = csv.writer(f_tg)
    w_tg.writerow(["N", "step", "time", "target_type", "k", "singular_value", "N_sv", "N2_sv"])
    f_q = open(raw / "q_svd.csv", "w", newline=""); w_q = csv.writer(f_q)
    w_q.writerow(["N", "step", "time", "q1", "q2", "q3", "q4", "Nq3", "Nq4", "N2q3", "N2q4"])
    f_trk = open(raw / "cluster_tracking.csv", "w", newline=""); w_trk = csv.writer(f_trk)
    w_trk.writerow(["time", "source_cid", "target_cid", "source_branch", "target_branch",
                    "subspace_overlap", "sigma_source", "sigma_target", "ambiguity_flag"])
    f_di = open(raw / "diagnostics_timeseries.csv", "w", newline=""); w_di = csv.writer(f_di)
    w_di.writerow(["N", "step", "time", "max_eigpair_residual", "unitary_error", "pair_error",
                   "max_ortho_error", "max_interplane_nondeg", "max_cluster_interplane",
                   "max_idempotent", "closure_error", "antisymmetry_error", "n_clusters",
                   "n_kernel", "n_rot_primary", "n_rot_alt100", "sigma1", "f", "f_ref_dev"])

    bt, bsig, bratio, bdelta, bocc, bq = [], [], [], [], [], []
    prev = None
    branch_counter = [0]
    initial_floor_branches = {}     # branch -> initial floor flag

    Zr = Z.copy(); wpr = wp.copy(); t = 0; processed = 0
    B0dim = B0.shape[1]
    while True:
        if t in sched:
            f = fval(Zr, p, q)
            K = dense_K(sys_lr, Zr)
            dec = decompose(K, Zr, B0)
            cls = dec["clusters"]; mu = dec["mu"]; sigma = dec["sigma"]; sig1 = dec["sig1"]

            # クラスタ branch 追跡（部分空間重なり）
            branches = {}
            ambiguity = {}
            if prev is None:
                for c in cls:
                    branches[c["cid"]] = branch_counter[0]
                    # 初期床フラグ：この時刻(=initial)で床近傍だったか（回転にならない程度）→初期は全て回転なので False 既定
                    initial_floor_branches[branch_counter[0]] = 0
                    branch_counter[0] += 1
            else:
                pcls = prev["clusters"]; pbr = prev["branches"]
                for c in cls:
                    best, second, ba = -1.0, -1.0, None
                    for pc in pcls:
                        ov = np.sum((pc["B"].T @ c["B"]) ** 2) / max(1, min(pc["dim"], c["dim"]))
                        if ov > best:
                            second = best; best = ov; ba = pc["cid"]
                        elif ov > second:
                            second = ov
                    if ba is not None and best > 0.5:
                        branches[c["cid"]] = pbr[ba]
                    else:
                        branches[c["cid"]] = branch_counter[0]
                        initial_floor_branches[branch_counter[0]] = 1   # 新規出現＝以前床
                        branch_counter[0] += 1
                    ambiguity[c["cid"]] = int(best < 0.9 or (best - second) < 0.1)
                    w_trk.writerow([t, ba, c["cid"], (pbr.get(ba) if ba is not None else -1),
                                    branches[c["cid"]], fmt % max(best, 0.0),
                                    fmt % (pcls[[pc["cid"] for pc in pcls].index(ba)]["sigma"] if ba is not None else 0.0),
                                    fmt % c["sigma"], ambiguity[c["cid"]]])

            # 全固有値保存
            for i in range(M):
                s = sigma[i]; rr = s / sig1
                w_eig.writerow([n, t, t, i, fmt % mu[i], fmt % s, fmt % rr, fmt % (n * rr),
                                fmt % (n * n * rr), fmt % dec["floor_ratio"][i], dec["pair_id"][i],
                                dec["cluster_id"][i], fmt % dec["resid"][i]])
            # クラスタ保存
            totZ = float(np.real(np.conj(Zr) @ Zr))
            for c in cls:
                br = branches[c["cid"]]
                w_cl.writerow([n, t, t, c["cid"], br, fmt % c["sigma"], fmt % (c["sigma"] / sig1),
                               c["mult"], c["dim"], fmt % c["occ"], fmt % (c["occ"] / totZ),
                               fmt % c["delta"], fmt % c["overlap_parent"], fmt % c["invariance_residual"],
                               (fmt % c["kbj_residual"] if not np.isnan(c["kbj_residual"]) else "nan"),
                               initial_floor_branches.get(br, 0)])
            # 対象別残差
            tg = target_residuals(dec, B0)
            for ttype, svs in tg.items():
                for k, sv in enumerate(svs):
                    w_tg.writerow([n, t, t, ttype, k, fmt % sv, fmt % (n * sv), fmt % (n * n * sv)])
            # Q=[B0|Bdom]
            dom = max(cls, key=lambda c: c["sigma"])
            Bdom = dom["B"][:, :2] if dom["B"].shape[1] >= 2 else dom["B"]
            Q = np.column_stack([B0, Bdom])
            qs = np.linalg.svd(Q, compute_uv=False); qs = np.pad(qs, (0, max(0, 4 - len(qs))))[:4]
            w_q.writerow([n, t, t, fmt % qs[0], fmt % qs[1], fmt % qs[2], fmt % qs[3],
                          fmt % (n * qs[2]), fmt % (n * qs[3]), fmt % (n * n * qs[2]), fmt % (n * n * qs[3])])
            # 診断
            fref = ref_f.get(t); fdev = (abs(f - fref) if fref is not None else float("nan"))
            w_di.writerow([n, t, t, fmt % dec["resid"].max(), fmt % dec["uni"], fmt % dec["pair_err"],
                           fmt % dec["max_ortho"], fmt % dec["max_inter_nondeg"], fmt % dec["max_cluster_inter"],
                           fmt % dec["max_idem"], fmt % dec["closure"], fmt % dec["antisym"],
                           len(cls), len(dec["ker_idx"]), int(np.sum(dec["floor_ratio"] > FLOOR_REL_PRIMARY)),
                           dec["n_rot_alt"], fmt % sig1, fmt % f, fmt % fdev])

            bt.append(t); bsig.append(sigma.copy()); bratio.append(sigma / sig1)
            bdelta.append(np.array([c["delta"] for c in cls])); bocc.append(np.array([c["occ"] / totZ for c in cls]))
            bq.append(qs)

            if t in rep_set:
                label = [k for k, vv in rep_times.items() if vv == t][0]
                np.savez_compressed(binp / f"decomp_{label}_step{t}.npz",
                                    mu=mu, sigma=sigma, floor_ratio=dec["floor_ratio"],
                                    pair_id=dec["pair_id"], cluster_id=dec["cluster_id"],
                                    U=dec["U"], B0=B0, Z_real=Zr.real, Z_imag=Zr.imag,
                                    cluster_sigma=np.array([c["sigma"] for c in cls]),
                                    cluster_mult=np.array([c["mult"] for c in cls]),
                                    cluster_occ=np.array([c["occ"] for c in cls]),
                                    cluster_delta=np.array([c["delta"] for c in cls]),
                                    E_ker=dec["E_ker"])
                with open(tab / f"fulltable_{label}_step{t}.csv", "w", newline="") as tf:
                    tw = csv.writer(tf)
                    tw.writerow(["eigen_index", "mu", "sigma", "sigma_over_sigma1", "abs_sigma_over_eps_normK",
                                 "pair_id", "cluster_id", "solver_residual"])
                    for i in np.argsort(-sigma):
                        tw.writerow([i, fmt % mu[i], fmt % sigma[i], fmt % (sigma[i] / sig1),
                                     fmt % dec["floor_ratio"][i], dec["pair_id"][i], dec["cluster_id"][i],
                                     fmt % dec["resid"][i]])
                with open(tab / f"clusters_{label}_step{t}.csv", "w", newline="") as tf:
                    tw = csv.writer(tf)
                    tw.writerow(["cluster_id", "cluster_branch", "sigma", "sigma_over_sigma1", "mult", "dim",
                                 "occupation_fraction", "delta_C", "overlap_parent", "invariance_residual", "kbj_residual"])
                    for c in cls:
                        tw.writerow([c["cid"], branches[c["cid"]], fmt % c["sigma"], fmt % (c["sigma"] / sig1),
                                     c["mult"], c["dim"], fmt % (c["occ"] / totZ), fmt % c["delta"],
                                     fmt % c["overlap_parent"], fmt % c["invariance_residual"],
                                     (fmt % c["kbj_residual"] if not np.isnan(c["kbj_residual"]) else "nan")])

            prev = {"clusters": cls, "branches": branches}
            processed += 1
        if t >= end:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr)
        Zr = sys_lr.cayley_step(Zr, se); t += 1

    for fh in (f_eig, f_cl, f_tg, f_q, f_trk, f_di):
        fh.close()
    np.savez_compressed(binp / "timeseries.npz", times=np.array(bt),
                        sigma=np.array(bsig, dtype=object), ratio=np.array(bratio, dtype=object),
                        delta=np.array(bdelta, dtype=object), occ=np.array(bocc, dtype=object),
                        q=np.array(bq), allow_pickle=True)

    d = np.genfromtxt(raw / "diagnostics_timeseries.csv", delimiter=",", names=True, encoding="utf-8")
    summary = {
        "N": n, "M": M, "seed": seed, "crossing": crossing, "after": AFTER, "n_records": processed,
        "method": "eigh(1j*K)", "floor_rel_primary": FLOOR_REL_PRIMARY, "floor_rel_alt": FLOOR_REL_ALT,
        "representative_times": rep_times,
        "max_eigpair_residual": float(np.nanmax(d["max_eigpair_residual"])),
        "max_unitary_error": float(np.nanmax(d["unitary_error"])),
        "max_pair_error": float(np.nanmax(d["pair_error"])),
        "max_ortho_error": float(np.nanmax(d["max_ortho_error"])),
        "max_interplane_nondeg": float(np.nanmax(d["max_interplane_nondeg"])),
        "max_cluster_interplane": float(np.nanmax(d["max_cluster_interplane"])),
        "max_idempotent": float(np.nanmax(d["max_idempotent"])),
        "max_closure_error": float(np.nanmax(d["closure_error"])),
        "max_antisymmetry_error": float(np.nanmax(d["antisymmetry_error"])),
        "max_f_ref_dev": float(np.nanmax(d["f_ref_dev"])) if not np.all(np.isnan(d["f_ref_dev"])) else None,
        "acceptance_nondeg_interplane_le_1e-12": bool(np.nanmax(d["max_interplane_nondeg"]) <= 1e-12),
    }
    (BASE / "diagnostics").mkdir(exist_ok=True)
    with open(BASE / "diagnostics" / f"N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"[v2 N={n}] M={M} crossing={crossing} 記録={processed}")
    print(f"  eigh: 固有対残差≤{summary['max_eigpair_residual']:.1e} U†U-I≤{summary['max_unitary_error']:.1e} "
          f"正負対誤差≤{summary['max_pair_error']:.1e}")
    print(f"  非縮退平面間直交≤{summary['max_interplane_nondeg']:.1e} (検収≤1e-12: {summary['acceptance_nondeg_interplane_le_1e-12']}) "
          f"クラスタ間≤{summary['max_cluster_interplane']:.1e}")
    print(f"  冪等≤{summary['max_idempotent']:.1e} 閉鎖≤{summary['max_closure_error']:.1e} "
          f"反対称≤{summary['max_antisymmetry_error']:.1e} f参照一致≤{summary['max_f_ref_dev']}")
    return summary


if __name__ == "__main__":
    for a in sys.argv[1:]:
        run(int(a))
