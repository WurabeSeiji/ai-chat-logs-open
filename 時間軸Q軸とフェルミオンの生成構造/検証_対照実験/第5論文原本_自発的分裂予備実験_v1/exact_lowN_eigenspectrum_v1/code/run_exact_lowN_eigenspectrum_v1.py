#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 第1段階：N=5,40 厳密全固有値・全固有方向 観測実験（解釈なし）。

論文6と同一軌道（同 seed・同 DELTA・同更新則・同 crossing 定義・crossing 後 50000 step）を
後処理し、各記録時刻の瞬時生成子 K(t) を**完全厳密固有分解**して全量を保存する：
  全 σ_j, σ_j/σ_1, Nσ_j/σ_1, N²σ_j/σ_1（全枝, 小枝も削除しない）
  各固有回転平面の実基底 B_j と射影 Π_j
  親平面外重なり o_j, δ_j（全モード）
  親平面外残差行列 R_all の全特異値 s_k^δ
  四基底結合行列 Q=[B_0|B_1] の特異値 q_1..q_4
  状態の各固有平面占有 E_j, 核占有, 閉鎖誤差
  順位別＋branch追跡の両データ
  数値診断（固有対残差・直交・冪等・平面間直交・閉鎖・反対称）

近似禁止：全行列 dense eig(K) のみ。低ランク/反復/上位限定/閾値削除/合算/縮退平均を用いない。
N=300 には着手しない。解釈は書かない。

使い方:
  python3 run_exact_lowN_eigenspectrum_v1.py 5
  python3 run_exact_lowN_eigenspectrum_v1.py 40
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
BASE = CODE_DIR.parent                     # exact_lowN_eigenspectrum_v1/
ENGINE_DIR = BASE.parent                    # 第5論文原本_..._v1/
sys.path.insert(0, str(ENGINE_DIR))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact

DELTA = 1e-15
GUARD = 500
LEARN = 1000
VALID = 1000
AFTER = 50000
SIG_TOL = 1e-9        # 正の回転率とみなす下限（削除でなく核判定のためだけ）
REF_DIR = ENGINE_DIR / "paper6_definitive_control_v1"


# ---------------- 軌道（同一条件で再生成） ----------------
def build_init(n):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    p1s, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)   # 初期親平面 P0
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=sys_lr.m)
    return sys_lr, v, B_p1, p, q, Z, wp, residual


def fval(Z, p, q):
    Zp = Z - p * (p @ Z) - q * (q @ Z)
    return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))


def dense_K(sys_lr, Z):
    sys_lr.set_theta(np.angle(Z))
    M = sys_lr.m
    return np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])


# ---------------- 厳密全固有分解 ----------------
def exact_decompose(K):
    """実反対称 K の完全固有分解。正の σ を降順、各モードの実平面基底 B_j(M×2) を返す。"""
    M = K.shape[0]
    w, V = np.linalg.eig(K)
    Knorm = np.linalg.norm(K, 2)
    order = np.argsort(-w.imag)             # imag 降順
    modes = []
    sig1 = None
    for i in order:
        s = float(w[i].imag)
        if s <= SIG_TOL:
            continue
        v = V[:, i]
        # 固有対残差
        r = np.linalg.norm(K @ v - w[i] * v) / (max(1.0, Knorm) * np.linalg.norm(v))
        # 実平面：Re v, Im v を正規直交化
        a, b = v.real.copy(), v.imag.copy()
        na = np.linalg.norm(a)
        e1 = a / na
        b = b - (b @ e1) * e1
        nb = np.linalg.norm(b)
        e2 = b / nb
        B = np.column_stack([e1, e2])
        modes.append({"sigma": s, "eig": complex(w[i]), "B": B, "res": float(r)})
        if sig1 is None:
            sig1 = s
    n_ker = M - 2 * len(modes)
    return modes, sig1, n_ker, Knorm


def plane_overlap(Bi, Bj):
    return 0.5 * float(np.sum((Bi.T @ Bj) ** 2))     # ½Tr(Π_i Π_j)


def occ_plane(B, Z):
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))


# ---------------- サンプリング時刻 ----------------
def sample_schedule(crossing):
    """初期・crossing近傍を高密度、長時間を粗く。全て事前固定。"""
    ts = set()
    end = crossing + AFTER
    for t in range(0, crossing + 2000 + 1, 5):
        ts.add(t)
    for t in range(crossing + 2000, crossing + 10000 + 1, 50):
        ts.add(t)
    for t in range(crossing + 10000, end + 1, 200):
        ts.add(t)
    ts.add(end)
    return ts


# ---------------- 本体 ----------------
def run(n):
    sys_lr, v, B0, p, q, Z, wp, parent_res = build_init(n)
    M = sys_lr.m
    seed = 40260722 + 1000 * n
    raw_dir = BASE / "raw" / f"N{n:05d}"
    bin_dir = BASE / "binary" / f"N{n:05d}"
    tab_dir = BASE / "tables" / f"N{n:05d}"
    for d in (raw_dir, bin_dir, tab_dir):
        d.mkdir(parents=True, exist_ok=True)

    # crossing 検出のため一度 f を走らせる（軌道は決定論的）
    # 効率のため：本走行で crossing 検出後にスケジュールを確定し、全 step を回して記録
    # 先に crossing を求める軽量パス
    Zc = Z.copy(); wpc = wp.copy(); crossing = None
    tt = 0
    while True:
        if fval(Zc, p, q) > 0.05:
            crossing = tt; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc)
        Zc = sys_lr.cayley_step(Zc, se); tt += 1
        if tt > 200000:
            raise RuntimeError("crossing 未検出")
    sched = sample_schedule(crossing)
    end = crossing + AFTER

    # 代表時刻（fig8・全モード表・固有ベクトルbinary用）
    def nearest(target):
        return min(sched, key=lambda x: abs(x - target))
    rep_times = {
        "initial": nearest(0),
        "pre_crossing": nearest(crossing - 50),
        "crossing": nearest(crossing),
        "post_crossing": nearest(crossing + 200),
        "plateau_start": nearest(crossing + 5000),
        "final": end,
    }
    rep_set = set(rep_times.values())

    # 参照 f（既存 obs CSV）との一致検証
    ref_f = {}
    ref_csv = REF_DIR / f"obs_N{n:05d}.csv"
    if ref_csv.exists():
        with open(ref_csv) as fh:
            r = csv.DictReader(fh)
            for row in r:
                ref_f[int(float(row["tau"]))] = float(row["f"])

    # 出力ファイル
    f_eig = open(raw_dir / "eigenvalues.csv", "w", newline="")
    w_eig = csv.writer(f_eig)
    w_eig.writerow(["N", "M", "seed", "step", "time", "crossing_relative_step", "rank_index",
                    "branch_id", "eigenvalue_real", "eigenvalue_imag", "sigma", "sigma_over_sigma1",
                    "N_sigma_over_sigma1", "N2_sigma_over_sigma1", "log10_abs_sigma",
                    "log10_sigma_over_sigma1", "solver_residual"])
    f_delta = open(raw_dir / "delta.csv", "w", newline="")
    w_delta = csv.writer(f_delta)
    w_delta.writerow(["N", "step", "time", "rank_index", "branch_id", "sigma", "sigma_over_sigma1",
                      "overlap_with_parent", "delta2", "delta"])
    f_occ = open(raw_dir / "occupation.csv", "w", newline="")
    w_occ = csv.writer(f_occ)
    w_occ.writerow(["N", "step", "time", "rank_index", "branch_id", "sigma", "sigma_over_sigma1",
                    "delta", "occupation_Ej", "occupation_fraction"])
    f_svd = open(raw_dir / "residual_svd.csv", "w", newline="")
    w_svd = csv.writer(f_svd)
    w_svd.writerow(["N", "step", "time", "k", "s_delta", "N_s_delta", "N2_s_delta"])
    f_q = open(raw_dir / "q_svd.csv", "w", newline="")
    w_q = csv.writer(f_q)
    w_q.writerow(["N", "step", "time", "q1", "q2", "q3", "q4",
                  "Nq3", "Nq4", "N2q3", "N2q4"])
    f_branch = open(raw_dir / "branch_tracking.csv", "w", newline="")
    w_branch = csv.writer(f_branch)
    w_branch.writerow(["time", "source_rank", "target_rank", "source_branch", "target_branch",
                       "plane_overlap", "sigma_source", "sigma_target", "tracking_ambiguity_flag"])
    f_diag = open(raw_dir / "diagnostics_timeseries.csv", "w", newline="")
    w_diag = csv.writer(f_diag)
    w_diag.writerow(["N", "step", "time", "max_eigpair_residual", "max_ortho_error",
                     "max_idempotent_error", "max_interplane_error", "closure_error",
                     "antisymmetry_error", "n_modes", "n_kernel", "sigma1", "f", "f_ref_dev"])

    # 時系列バイナリ蓄積
    bin_times, bin_sigma, bin_ratio, bin_delta, bin_E, bin_sdelta, bin_q = [], [], [], [], [], [], []

    prev = None                # (times, modes list) for branch tracking
    branch_counter = [0]
    fmt = "%.17e"

    # 全 step 走行、sched 時刻で厳密解析
    Zr = Z.copy(); wpr = wp.copy()
    t = 0
    processed = 0
    while True:
        if t in sched:
            f = fval(Zr, p, q)
            K = dense_K(sys_lr, Zr)
            modes, sig1, n_ker, Knorm = exact_decompose(K)
            nm = len(modes)
            # 派生量
            sigmas = np.array([mm["sigma"] for mm in modes])
            ratios = sigmas / sig1
            Bs = [mm["B"] for mm in modes]
            # δ, 占有
            o = np.array([plane_overlap(B0, B) for B in Bs])
            delta2 = np.maximum(0.0, 1 - o)
            delta = np.sqrt(delta2)
            E = np.array([occ_plane(B, Zr) for B in Bs])
            totZ = float(np.real(np.conj(Zr) @ Zr))
            E_ker = totZ - float(E.sum())
            closure = abs(totZ - float(E.sum()) - E_ker)   # =0 by construction; 別途厳密核占有下記
            # 厳密核占有：I - Σ Π_j を Z に作用
            PZr = Zr.real.copy(); PZi = Zr.imag.copy()
            acc_r = np.zeros(M); acc_i = np.zeros(M)
            for B in Bs:
                acc_r += B @ (B.T @ Zr.real); acc_i += B @ (B.T @ Zr.imag)
            ker_r = Zr.real - acc_r; ker_i = Zr.imag - acc_i
            E_ker_exact = float(ker_r @ ker_r + ker_i @ ker_i)
            closure_exact = abs(totZ - float(E.sum()) - E_ker_exact)
            # 残差行列 R_all の全特異値
            Rcols = []
            for B in Bs:
                Rj = B - B0 @ (B0.T @ B)      # (I-Π0)B_j, M×2
                Rcols.append(Rj)
            Rall = np.column_stack(Rcols) if Rcols else np.zeros((M, 1))
            s_delta = np.linalg.svd(Rall, compute_uv=False)
            # 四基底結合 Q=[B0|B1]
            Q = np.column_stack([B0, Bs[0]])
            qs = np.linalg.svd(Q, compute_uv=False)
            qs = np.pad(qs, (0, max(0, 4 - len(qs))))[:4]
            # 数値診断
            max_res = max((mm["res"] for mm in modes), default=0.0)
            max_ortho = max((np.linalg.norm(B.T @ B - np.eye(2)) for B in Bs), default=0.0)
            max_idem = 0.0
            for B in Bs:
                P = B @ B.T
                max_idem = max(max_idem, np.linalg.norm(P @ P - P) / max(1.0, M))
            # 平面間直交（サンプル：隣接順位対のみ全対は O(nm²)、nm=40で可）
            max_inter = 0.0
            for i in range(nm):
                for jj in range(i + 1, nm):
                    max_inter = max(max_inter, np.linalg.norm(Bs[i].T @ Bs[jj]))
            antisym = np.linalg.norm(K + K.T) / max(1.0, Knorm)

            # branch tracking（隣接記録時刻）
            branches = [None] * nm
            ambiguity = [0] * nm
            if prev is None:
                for j in range(nm):
                    branches[j] = branch_counter[0]; branch_counter[0] += 1
            else:
                pmodes = prev["modes"]; pbr = prev["branches"]
                pn = len(pmodes)
                ov = np.zeros((pn, nm))
                for a in range(pn):
                    for bb in range(nm):
                        ov[a, bb] = plane_overlap(pmodes[a]["B"], Bs[bb])
                # 各 target に対し最大 overlap の source
                for bb in range(nm):
                    a = int(np.argmax(ov[:, bb]))
                    best = ov[a, bb]
                    second = np.partition(ov[:, bb], -2)[-2] if pn > 1 else 0.0
                    branches[bb] = pbr[a]
                    ambiguity[bb] = int(best < 0.9 or (best - second) < 0.1)
                    w_branch.writerow([t, a, bb, pbr[a], branches[bb], fmt % best,
                                       fmt % pmodes[a]["sigma"], fmt % modes[bb]["sigma"], ambiguity[bb]])
                # 同一 source が複数 target に付いた場合の新規枝付与（衝突時）
                seen = {}
                for bb in range(nm):
                    if branches[bb] in seen:
                        branches[bb] = branch_counter[0]; branch_counter[0] += 1
                    seen[branches[bb]] = bb

            # --- 書き出し ---
            for j in range(nm):
                s = sigmas[j]; rr = ratios[j]
                w_eig.writerow([n, M, seed, t, t, t - crossing, j, branches[j],
                                fmt % modes[j]["eig"].real, fmt % modes[j]["eig"].imag,
                                fmt % s, fmt % rr, fmt % (n * rr), fmt % (n * n * rr),
                                fmt % np.log10(abs(s)), fmt % np.log10(rr), fmt % modes[j]["res"]])
                w_delta.writerow([n, t, t, j, branches[j], fmt % s, fmt % rr,
                                  fmt % o[j], fmt % delta2[j], fmt % delta[j]])
                w_occ.writerow([n, t, t, j, branches[j], fmt % s, fmt % rr, fmt % delta[j],
                                fmt % E[j], fmt % (E[j] / totZ)])
            for k, sv in enumerate(s_delta):
                w_svd.writerow([n, t, t, k, fmt % sv, fmt % (n * sv), fmt % (n * n * sv)])
            w_q.writerow([n, t, t, fmt % qs[0], fmt % qs[1], fmt % qs[2], fmt % qs[3],
                          fmt % (n * qs[2]), fmt % (n * qs[3]), fmt % (n * n * qs[2]), fmt % (n * n * qs[3])])
            fref = ref_f.get(t)
            fdev = (abs(f - fref) if fref is not None else float("nan"))
            w_diag.writerow([n, t, t, fmt % max_res, fmt % max_ortho, fmt % max_idem,
                             fmt % max_inter, fmt % closure_exact, fmt % antisym,
                             nm, n_ker, fmt % sig1, fmt % f, fmt % fdev])

            # バイナリ蓄積
            bin_times.append(t)
            bin_sigma.append(sigmas); bin_ratio.append(ratios)
            bin_delta.append(delta); bin_E.append(E / totZ)
            bin_sdelta.append(s_delta); bin_q.append(qs)

            # 代表時刻：固有ベクトル・平面基底・射影をバイナリ保存＋全モード表
            if t in rep_set:
                label = [k for k, vv in rep_times.items() if vv == t][0]
                np.savez_compressed(bin_dir / f"planes_{label}_step{t}.npz",
                                    sigma=sigmas, ratio=ratios, delta=delta, E=E / totZ,
                                    B_stack=np.stack(Bs, axis=0), B0=B0,
                                    s_delta=s_delta, q=qs, Z_real=Zr.real, Z_imag=Zr.imag,
                                    eig_real=np.array([mm["eig"].real for mm in modes]),
                                    eig_imag=np.array([mm["eig"].imag for mm in modes]))
                with open(tab_dir / f"fulltable_{label}_step{t}.csv", "w", newline="") as tf:
                    tw = csv.writer(tf)
                    tw.writerow(["j", "branch_id", "sigma", "sigma_over_sigma1", "N_sigma_over_sigma1",
                                 "N2_sigma_over_sigma1", "delta", "E_j", "occupation_fraction",
                                 "solver_residual"])
                    for j in range(nm):
                        tw.writerow([j, branches[j], fmt % sigmas[j], fmt % ratios[j],
                                     fmt % (n * ratios[j]), fmt % (n * n * ratios[j]),
                                     fmt % delta[j], fmt % E[j], fmt % (E[j] / totZ),
                                     fmt % modes[j]["res"]])

            prev = {"modes": modes, "branches": branches}
            processed += 1

        if t >= end:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr)
        Zr = sys_lr.cayley_step(Zr, se); t += 1

    for fh in (f_eig, f_delta, f_occ, f_svd, f_q, f_branch, f_diag):
        fh.close()

    # 時系列バイナリ（不揃い長は object 配列）
    np.savez_compressed(bin_dir / "timeseries.npz",
                        times=np.array(bin_times),
                        sigma=np.array(bin_sigma, dtype=object),
                        ratio=np.array(bin_ratio, dtype=object),
                        delta=np.array(bin_delta, dtype=object),
                        E=np.array(bin_E, dtype=object),
                        s_delta=np.array(bin_sdelta, dtype=object),
                        q=np.array(bin_q, dtype=object), allow_pickle=True)

    # 参照一致・診断集約
    with open(REF_DIR / f"obs_N{n:05d}.csv") as _:
        pass
    diag = np.genfromtxt(raw_dir / "diagnostics_timeseries.csv", delimiter=",",
                         names=True, dtype=float, encoding="utf-8")
    summary = {
        "N": n, "M": M, "seed": seed, "crossing": crossing, "after": AFTER,
        "n_records": processed, "n_modes": int(diag["n_modes"][0]),
        "n_kernel": int(diag["n_kernel"][0]),
        "representative_times": rep_times,
        "max_eigpair_residual": float(np.nanmax(diag["max_eigpair_residual"])),
        "max_ortho_error": float(np.nanmax(diag["max_ortho_error"])),
        "max_idempotent_error": float(np.nanmax(diag["max_idempotent_error"])),
        "max_interplane_error": float(np.nanmax(diag["max_interplane_error"])),
        "max_closure_error": float(np.nanmax(diag["closure_error"])),
        "max_antisymmetry_error": float(np.nanmax(diag["antisymmetry_error"])),
        "max_f_ref_dev": float(np.nanmax(diag["f_ref_dev"])) if not np.all(np.isnan(diag["f_ref_dev"])) else None,
    }
    (BASE / "diagnostics").mkdir(exist_ok=True)
    with open(BASE / "diagnostics" / f"N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"[N={n}] M={M} crossing={crossing} 記録時刻数={processed} モード数={summary['n_modes']} 核={summary['n_kernel']}")
    print(f"  数値診断: 固有対残差≤{summary['max_eigpair_residual']:.1e} 直交≤{summary['max_ortho_error']:.1e} "
          f"冪等≤{summary['max_idempotent_error']:.1e} 平面間≤{summary['max_interplane_error']:.1e}")
    print(f"           閉鎖≤{summary['max_closure_error']:.1e} 反対称≤{summary['max_antisymmetry_error']:.1e} "
          f"f参照一致≤{summary['max_f_ref_dev']}")
    return summary


if __name__ == "__main__":
    for a in sys.argv[1:]:
        run(int(a))
