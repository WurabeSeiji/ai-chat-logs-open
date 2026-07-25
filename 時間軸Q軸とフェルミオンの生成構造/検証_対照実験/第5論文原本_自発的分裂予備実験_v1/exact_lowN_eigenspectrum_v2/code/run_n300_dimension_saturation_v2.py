#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 N=300 有効4方向飽和検査（Gram縮約 G=WᵀW、計算量最適化版）。解釈なし。

巨大密行列 K∈ℝ^{44850×44850} を作らず、K=WJWᵀ の低ランク構造を用い、G=WᵀW（2N×2N≤600）の
固有分解で縮約反対称生成子 K_r=Λ^{1/2}V_rᵀ J V_r Λ^{1/2}（≤600次元）を作り eigh(iK_r) する。
支配平面2列のみ辺空間へ持ち上げる。数値精度内で密行列法と同値（N=40で検証済）。

検査A：[B0|Bdom(t)] の q1..q4 と結合階数。
検査B：全時間の支配平面が単一4次元へ閉じるか（閉鎖残差、5本目のストリーミング探索）。
検査C：間引き結合の全時間特異値 s1..s10、s5/s1。
同一コードで N=5,40,300。N=300 は B_dom を間引き保存（メモリ節約）。

使い方: python3 run_n300_dimension_saturation_v2.py 300
        python3 run_n300_dimension_saturation_v2.py 40   # 検収・比較用
"""

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
BASE = CODE.parent
ENGINE = BASE.parent
sys.path.insert(0, str(ENGINE))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed

DELTA = 1e-15
AFTER = 50000
TAU_G_PRIMARY = 1e-12
TAU_G_SWEEP = [1e-10, 1e-11, 1e-12, 1e-13, 1e-14]
TAU_STREAM_SWEEP = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
TAU_STREAM_PRIMARY = 1e-10
Q_REL_TAU = 1e-8            # q_j > τ q1 判定
FMT = "%.17e"
EPS = np.finfo(float).eps
MAX_STORE = 220            # 間引き保存する支配平面の最大時刻数（閉鎖・大域SVD用）


def gram_reduce(sys_lr, Z, tau_G=TAU_G_PRIMARY):
    sys_lr.set_theta(np.angle(Z))
    G = 0.5 * (sys_lr.G + sys_lr.G.T)
    Gnorm = np.linalg.norm(G) + 1e-300
    lam, V = np.linalg.eigh(G)
    order = np.argsort(-lam)
    lam = lam[order]; V = V[:, order]
    lmax = lam[0]
    keep = lam > tau_G * lmax
    lam_r = lam[keep]; V_r = V[:, keep]; sq = np.sqrt(lam_r)
    S = V_r * sq[None, :]
    Kr = S.T @ sys_lr.J @ S
    Kr = 0.5 * (Kr - Kr.T)
    Krnorm = np.linalg.norm(Kr) + 1e-300
    mu, Ur = np.linalg.eigh(1j * Kr)
    diag = {
        "gram_symmetry_error": float(np.linalg.norm(sys_lr.G - sys_lr.G.T) / max(1.0, Gnorm)),
        "gram_reconstruction_error": float(
            np.linalg.norm(G - (V_r * lam_r[None, :]) @ V_r.T) / max(1.0, Gnorm)),
        "reduced_skew_error": float(np.linalg.norm(Kr + Kr.T) / max(1.0, Krnorm)),
        "r_G": int(keep.sum()),
    }
    return {"lam_all": lam, "lam_r": lam_r, "V_r": V_r, "sq": sq, "Kr": Kr,
            "mu": mu, "Ur": Ur, "lmax": lmax, "diag": diag}


def dominant_plane(sys_lr, gr):
    mu, Ur, V_r, sq = gr["mu"], gr["Ur"], gr["V_r"], gr["sq"]
    idx = int(np.argmax(mu))                 # 最大正 μ = σ_dom
    sig = float(mu[idx]); ur = Ur[:, idx]
    coef = V_r @ (ur / sq)                    # 2N 複素
    u = sys_lr.w(coef.real) + 1j * sys_lr.w(coef.imag)   # 辺空間 M
    b = np.column_stack([np.sqrt(2) * u.real, np.sqrt(2) * u.imag])
    Q, R = np.linalg.qr(b)
    B = Q[:, :2]
    # 持ち上げ固有対残差（密行列を作らず Ku=WJ(Wᵀu)）
    Wt_u = sys_lr.wt(u)                       # 2N 複素
    Ku = sys_lr.w((sys_lr.J @ Wt_u.real)) + 1j * sys_lr.w((sys_lr.J @ Wt_u.imag))
    lifted_res = float(np.linalg.norm(Ku + 1j * sig * u) / max(1.0, abs(sig)))
    orthB = float(np.linalg.norm(B.T @ B - np.eye(2)))
    return sig, B, u, lifted_res, orthB


def qsv4(B0, Bd):
    Q4 = np.column_stack([B0, Bd])
    C4 = Q4.T @ Q4
    ev = np.clip(np.linalg.eigvalsh(C4)[::-1], 0, None)
    q = np.sqrt(ev)
    orthC = None
    return q, C4


def principal_angles(A, B):
    s = np.clip(np.linalg.svd(A.T @ B, compute_uv=False), -1, 1)
    return np.degrees(np.arccos(s))


def sample_schedule(cr):
    ts = set()
    for t in range(0, cr + 2000 + 1, 5):
        ts.add(t)
    for t in range(cr + 2000, cr + 10000 + 1, 50):
        ts.add(t)
    for t in range(cr + 10000, cr + AFTER + 1, 200):
        ts.add(t)
    ts.add(cr + AFTER)
    return ts


def run(n):
    assert n in (5, 40, 300)
    sys_lr = LowRankSystem(n); M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    # B0 = 親 v の支配平面（Gram法）
    gr0 = gram_reduce(sys_lr, v)
    _, B0, _, _, _ = dominant_plane(sys_lr, gr0)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g; Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p; q = q / np.linalg.norm(q)
    wp = rng.normal(size=M)

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    # crossing
    Zc = Z.copy(); wpc = wp.copy(); crossing = None; t = 0
    while True:
        if fval(Zc) > 0.05:
            crossing = t; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc); Zc = sys_lr.cayley_step(Zc, se); t += 1
    sched = sample_schedule(crossing); end = crossing + AFTER
    times = sorted(sched)
    store_every = max(1, len(times) // MAX_STORE)

    outdir = BASE / "raw" / f"N{n:05d}_dimension_saturation_v2"
    figdir = BASE / "figures" / f"N{n:05d}_dimension_saturation_v2"
    repdir = BASE / "reports"
    for d in (outdir, figdir, repdir):
        d.mkdir(parents=True, exist_ok=True)

    # 閉鎖候補 C4 を作る代表時刻（crossing直後/準安定前半/中央/最終）
    rep_for_C4 = {min(times, key=lambda x: abs(x - (crossing + 200))): "post_crossing",
                  min(times, key=lambda x: abs(x - (crossing + 5000))): "meta_early",
                  min(times, key=lambda x: abs(x - (crossing + 25000))): "meta_mid",
                  end: "final"}

    # 出力ファイル
    f_q = open(outdir / f"q_svd_N{n:05d}.csv", "w", newline=""); w_q = csv.writer(f_q)
    w_q.writerow(["step", "time", "relative_time", "gram_rank", "dominant_eigenvalue",
                  "q1", "q2", "q3", "q4", "rank_q", "theta1_deg", "theta2_deg", "orth_B0", "orth_Bdom"])
    f_di = open(outdir / f"diagnostics_N{n:05d}.csv", "w", newline=""); w_di = csv.writer(f_di)
    w_di.writerow(["step", "time", "gram_symmetry_error", "gram_reconstruction_error",
                   "reduced_skew_error", "reduced_eigen_residual", "lifted_eigen_residual",
                   "orth_Bdom", "orth_C4", "conservation_error"])
    f_bm = open(outdir / f"benchmark_N{n:05d}.csv", "w", newline=""); w_bm = csv.writer(f_bm)
    w_bm.writerow(["step", "time_build_G", "time_eigh_G", "time_build_Kr", "time_eigh_Kr",
                   "time_lift_dominant", "time_q", "time_closure", "time_total", "peak_memory_mb"])

    # ストリーミング基底（各 τ_stream）
    streams = {tau: B0.copy() for tau in TAU_STREAM_SWEEP}
    stream_hist = {tau: [] for tau in TAU_STREAM_SWEEP}
    # gram rank sweep 記録
    grank_sweep = []
    # 間引き保存
    stored_times = []; stored_B = []
    rep_B = {}

    STREAM_CAP = 2 * sys_lr.n + 20         # span(W)=2N が上限。安全余裕付き

    def stream_update(C, cols, tau):
        added = 0
        for col in cols.T:
            if C.shape[1] >= STREAM_CAP:
                break
            vv = col.copy()
            vv -= C @ (C.T @ vv); vv -= C @ (C.T @ vv)
            nv = np.linalg.norm(vv)
            if nv > tau:
                C = np.column_stack([C, vv / nv]); added += 1
        return C, added

    Zr = Z.copy(); wpr = wp.copy(); t = 0; ridx = 0
    Zprev_norm = None
    while True:
        if t in sched:
            tt0 = time.time()
            tb = time.time(); gr = gram_reduce(sys_lr, Zr, TAU_G_PRIMARY); t_gr = time.time() - tb
            sig, Bd, u, lifted_res, orthB = dominant_plane(sys_lr, gr)
            # reduced eigen residual
            idx = int(np.argmax(gr["mu"])); ur = gr["Ur"][:, idx]
            red_res = float(np.linalg.norm(1j * gr["Kr"] @ ur - gr["mu"][idx] * ur) / max(1.0, abs(gr["mu"][idx])))
            tq = time.time(); q4, C4mat = qsv4(B0, Bd); t_q = time.time() - tq
            rank_q = int(np.sum(q4 > Q_REL_TAU * q4[0]))
            ang = principal_angles(B0, Bd)
            th1 = float(ang[0]); th2 = float(ang[1]) if len(ang) > 1 else 0.0
            conservation = abs(float(np.real(np.conj(Zr) @ Zr)) - 1.0)
            # gram rank sweep
            gr_ranks = {tau: int(np.sum(gr["lam_all"] > tau * gr["lmax"])) for tau in TAU_G_SWEEP}
            grank_sweep.append((t, gr_ranks))
            # streaming
            for tau in TAU_STREAM_SWEEP:
                streams[tau], added = stream_update(streams[tau], Bd, tau)
                stream_hist[tau].append((t, streams[tau].shape[1], added))
            # store subsample + rep
            if ridx % store_every == 0 or t in rep_for_C4:
                stored_times.append(t); stored_B.append(Bd.copy())
            if t in rep_for_C4:
                rep_B[t] = Bd.copy()
            t_tot = time.time() - tt0
            w_q.writerow([t, t, t - crossing, gr["diag"]["r_G"], FMT % sig,
                          FMT % q4[0], FMT % q4[1], FMT % q4[2], FMT % q4[3], rank_q,
                          FMT % th1, FMT % th2, FMT % float(np.linalg.norm(B0.T @ B0 - np.eye(2))), FMT % orthB])
            w_di.writerow([t, t, FMT % gr["diag"]["gram_symmetry_error"],
                           FMT % gr["diag"]["gram_reconstruction_error"],
                           FMT % gr["diag"]["reduced_skew_error"], FMT % red_res, FMT % lifted_res,
                           FMT % orthB, "nan", FMT % conservation])
            if ridx < 10:
                w_bm.writerow([t, FMT % 0, FMT % 0, FMT % 0, FMT % 0, FMT % 0, FMT % t_q, FMT % 0,
                               FMT % t_tot, FMT % 0])
            ridx += 1
        if t >= end:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr); Zr = sys_lr.cayley_step(Zr, se); t += 1
    for fh in (f_q, f_di, f_bm):
        fh.close()

    # ---- 検査B：閉鎖残差（4候補C4）----
    with open(outdir / f"closure4_N{n:05d}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["candidate_step", "candidate_label", "step", "time", "relative_time",
                    "eps_fro", "eps_op", "eta4"])
        closure_summary = {}
        for cstep, clabel in rep_for_C4.items():
            Craw = np.column_stack([B0, rep_B[cstep]])
            C4, _ = np.linalg.qr(Craw)
            C4 = C4[:, :4]
            epsF = []
            for st, Bd in zip(stored_times, stored_B):
                R4 = Bd - C4 @ (C4.T @ Bd)
                ef = float(np.linalg.norm(R4)); eo = float(np.linalg.norm(R4, 2))
                eta = ef ** 2 / float(np.linalg.norm(Bd) ** 2)
                epsF.append(ef)
                w.writerow([cstep, clabel, st, st, st - crossing, FMT % ef, FMT % eo, FMT % eta])
            closure_summary[clabel] = {"candidate_step": cstep,
                                       "eps_fro_max": float(np.max(epsF)),
                                       "eta4_max": float(max((float(np.linalg.norm(Bd - C4 @ (C4.T @ Bd)) ** 2 /
                                                                     np.linalg.norm(Bd) ** 2)) for Bd in stored_B))}

    # ---- 検査C：間引き結合の全時間特異値 ----
    Bsub = np.column_stack([B0] + stored_B)
    GB = Bsub.T @ Bsub
    sB = np.sqrt(np.clip(np.linalg.eigvalsh(GB)[::-1], 0, None))
    with open(outdir / f"global_svd_N{n:05d}.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["index", "singular_value", "relative_to_s1"])
        for i, s in enumerate(sB[:20]):
            w.writerow([i + 1, FMT % s, FMT % (s / sB[0])])
    s5_s1 = float(sB[4] / sB[0]) if len(sB) > 4 else 0.0

    # ---- 集約 ----
    endrow = None
    with open(outdir / f"q_svd_N{n:05d}.csv") as fh:
        rows = list(csv.DictReader(fh))
        endrow = rows[-1]
        meta_row = min(rows, key=lambda r: abs(int(float(r["time"])) - (crossing + 5000)))
    stream_final = {("%.0e" % tau): stream_hist[tau][-1][1] for tau in TAU_STREAM_SWEEP}
    summary = {
        "N": n, "M": M, "crossing": crossing, "n_records": ridx, "store_count": len(stored_times),
        "gram_rank_final": int(endrow["gram_rank"]),
        "q3_final": float(endrow["q3"]), "q4_final": float(endrow["q4"]),
        "q3_meta": float(meta_row["q3"]), "q4_meta": float(meta_row["q4"]),
        "rank_q_final": int(endrow["rank_q"]),
        "global_s": [float(x) for x in sB[:10]], "s5_over_s1": s5_s1,
        "closure": closure_summary,
        "stream_final_dim_by_tau": stream_final,
        "gram_rank_sweep_final": {("%.0e" % tau): int(np.sum(gr["lam_all"] > tau * gr["lmax"])) for tau in TAU_G_SWEEP},
    }
    with open(BASE / "diagnostics" / f"N{n:05d}_saturation.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    make_figures(n, outdir, figdir, crossing)
    print(f"[飽和検査 N={n}] M={M} crossing={crossing} 記録={ridx} 間引き保存={len(stored_times)} gram_rank={summary['gram_rank_final']}")
    print(f"  q3_final={summary['q3_final']:.4f} q4_final={summary['q4_final']:.4f} rank_q={summary['rank_q_final']} "
          f"(準安定 q3={summary['q3_meta']:.4f} q4={summary['q4_meta']:.4f})")
    print(f"  大域特異値 s1..s6={[f'{x:.3e}' for x in summary['global_s'][:6]]}  s5/s1={s5_s1:.3e}")
    print(f"  閉鎖 eps_fro_max(候補別): " + " ".join(f"{k}={v['eps_fro_max']:.2e}" for k, v in closure_summary.items()))
    print(f"  ストリーミング基底次数(τ別, 終): {stream_final}")
    return summary


def make_figures(n, outdir, figdir, crossing):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def load(name):
        with open(outdir / name) as fh:
            return list(csv.DictReader(fh))
    q = load(f"q_svd_N{n:05d}.csv")
    t = np.array([int(float(r["time"])) for r in q]); rt = t - crossing
    q1 = np.array([float(r["q1"]) for r in q]); q2 = np.array([float(r["q2"]) for r in q])
    q3 = np.array([float(r["q3"]) for r in q]); q4 = np.array([float(r["q4"]) for r in q])
    gr = np.array([float(r["gram_rank"]) for r in q])
    # 図1 q1-4
    fig, ax = plt.subplots(figsize=(10, 6))
    for y, l in [(q1, "q1"), (q2, "q2"), (q3, "q3"), (q4, "q4")]:
        ax.plot(rt, y, lw=0.7, label=l)
    ax.axvline(0, color="k", ls=":", lw=0.8); ax.legend(); ax.set_xlabel("time - crossing"); ax.set_ylabel("q")
    ax.set_title(f"N={n} sat-Fig1: q1..q4 of [B0|Bdom]"); fig.tight_layout()
    fig.savefig(figdir / "fig1_q1234.png", dpi=130); plt.close(fig)
    # 図2 q3,q4 拡大
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rt, q3, lw=0.7, label="q3"); ax.plot(rt, q4, lw=0.7, label="q4")
    ax.axvline(0, color="k", ls=":", lw=0.8); ax.legend(); ax.set_xlabel("time - crossing"); ax.set_ylabel("q3,q4")
    ax.set_title(f"N={n} sat-Fig2: q3,q4 zoom"); fig.tight_layout()
    fig.savefig(figdir / "fig2_q3q4.png", dpi=130); plt.close(fig)
    # 図3 gram rank
    fig, ax = plt.subplots(figsize=(10, 5)); ax.plot(rt, gr, lw=0.8)
    ax.axvline(0, color="k", ls=":", lw=0.8); ax.set_xlabel("time - crossing"); ax.set_ylabel("gram_rank r_G")
    ax.set_title(f"N={n} sat-Fig3: Gram rank"); fig.tight_layout()
    fig.savefig(figdir / "fig3_gram_rank.png", dpi=130); plt.close(fig)
    # 図4/5 閉鎖残差
    try:
        cl = load(f"closure4_N{n:05d}.csv")
        labels = sorted(set(r["candidate_label"] for r in cl))
        for logy, tag in [(False, "fig4_closure_lin"), (True, "fig5_closure_log")]:
            fig, ax = plt.subplots(figsize=(10, 6))
            for lb in labels:
                rr = [r for r in cl if r["candidate_label"] == lb]
                x = np.array([int(float(r["step"])) for r in rr]) - crossing
                y = np.array([float(r["eps_fro"]) for r in rr])
                o = np.argsort(x)
                (ax.semilogy if logy else ax.plot)(x[o], np.clip(y[o], 1e-20, None) if logy else y[o], lw=0.7, label=lb)
            ax.axvline(0, color="k", ls=":", lw=0.8); ax.legend(fontsize=8)
            ax.set_xlabel("time - crossing"); ax.set_ylabel("eps_fro (4D closure residual)")
            ax.set_title(f"N={n} sat-{tag}"); fig.tight_layout(); fig.savefig(figdir / f"{tag}.png", dpi=130); plt.close(fig)
    except Exception:
        pass
    # 図8/9 大域特異値
    try:
        gs = load(f"global_svd_N{n:05d}.csv")
        idx = np.array([int(r["index"]) for r in gs]); sv = np.array([float(r["singular_value"]) for r in gs])
        rel = np.array([float(r["relative_to_s1"]) for r in gs])
        fig, ax = plt.subplots(figsize=(9, 5)); ax.plot(idx, sv, "o-", ms=4)
        ax.set_xlabel("index"); ax.set_ylabel("singular value"); ax.set_title(f"N={n} sat-Fig8: global combined SVD")
        fig.tight_layout(); fig.savefig(figdir / "fig8_global_svd.png", dpi=130); plt.close(fig)
        fig, ax = plt.subplots(figsize=(9, 5)); ax.semilogy(idx, np.clip(rel, 1e-20, None), "o-", ms=4)
        ax.axhline(1e-8, color="r", ls="--", lw=0.8); ax.set_xlabel("index"); ax.set_ylabel("s_j/s_1 (log)")
        ax.set_title(f"N={n} sat-Fig9: s_j/s_1"); fig.tight_layout(); fig.savefig(figdir / "fig9_sj_s1.png", dpi=130); plt.close(fig)
    except Exception:
        pass


if __name__ == "__main__":
    for a in sys.argv[1:]:
        run(int(a))
