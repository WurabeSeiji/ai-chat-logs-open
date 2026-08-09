#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N 掃引 1→20（同一シード・同一条件）— 何が生まれ、何が生まれないかを全 N で読む

目的（木原指示 2026-08-08）:
  **「作れないこと」を調べるのも目的である。** 電荷を持たないフェルミオン場の
  真空が作れない N が存在すること自体が結果であり、避けて通ってはならない。
  N=1 には何もなく、N=2 から何かは測れるが空間はまだ生まれない——**何が
  生まれないか**を全 N について記録する。これにより、どこまでが原子・分子
  オーダの系で、どこからが上位構造のオーダの系かの境界を見る。
  テストベッド N=12 がその階梯のどこに位置するのかを判定する。

  **万能関数が NaN を返してもアボートしない。** 構成が失敗した N も、その事実を
  結果として記録し、任意の N について図化する。NaN の扱いは本スクリプト
  （上位側）の責務である。

  【是正記録】資格審査（run_qualification_dimension_v1.py）では make_parent が
  N=8 で失敗した際にシードを退避させて無理に構成した。あれは「その条件では
  その N の真空が作れない」という結果を消す操作だった。本掃引では**シードを
  固定し、失敗は失敗としてそのまま記録する**。

使用する万能関数（2026-08-08 分岐版の3本＋選択層）:
  F = ../統一万能関数_v1/unified_interaction_v2.py
  D = ../統一万能関数_v1/unified_dimension_v1.py
  G = ../統一万能関数_v1/unified_readout_v3.py
  S = ../統一万能関数_v1/selection_v1.py

条件（全 N 共通・すべて宣言値）:
  帯 Nn=16・毛 Nη=8・T=4000・ニュートリノ型シード（帯 k=1・巻き η=0・δ=10⁻²・
  親 seed=2 固定）・ポンプ（帯 k=2・巻き η=0）・D の宣言引数 cell=(2,0)・order=6・
  準安定窓 τ∈[2000,4000]・選択 s_clock_acquirable（既定床）
  対照として各 N で真空宇宙（δ=0）も走らせる。

記録する「生まれた／生まれない」（各 N）:
  構成      : 宇宙そのものが構成できたか（失敗理由も記録）
  空間      : crossing（f₂>0.05）が起きるか
  次元      : 第3次元の確定度 align・実効平面数 n_eff・直交化残差
  物質      : f_seed が立つか
  時間      : 時計の担い手が現れるか・ω̂ が π/72 に定着するか
  凝縮体    : ゼロ閉塞ブロックができるか（セル別閉塞残差）

図:
  各 N: fig_nsweep_N{n}_v1.png（4段: 空間 f₂／次元 n_eff・align／時間 担い手・ω̂／閉塞）
  総括: fig_nsweep_summary_v1.png（N 依存性: 各量の準安定窓中央値と誕生時刻）
        fig_nsweep_birth_matrix_v1.png（N × 何が生まれたかの可否行列）

使い方: python3 run_tb_nsweep_1to20_v1.py [Nmin Nmax]（省略時 1 20）
"""
from __future__ import annotations
import importlib.util, json, sys, time, traceback
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
UF = HERE.parent / "統一万能関数_v1"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F = load("f_ns", UF / "unified_interaction_v2.py")
D = load("d_ns", UF / "unified_dimension_v1.py")
G = load("g_ns", UF / "unified_readout_v3.py")
S = load("s_ns", UF / "selection_v1.py")

NN, NETA = 16, 8
T = 4000
DELTA = 1e-2
SEED = 2                 # 固定（退避しない）
CELL = (2, 0)
ORDER = 6
WIN = (2000, 4000)
OMEGA_REF = np.pi / 72.0
KEYS = ("f2", "f_seed", "align", "n_eff", "nonunif", "closure",
        "carrier_power", "coherence", "omega_hat", "phi_weight")


def first_true(mask):
    idx = np.flatnonzero(mask)
    return int(idx[0]) + 1 if len(idx) else None


def run_one(n, delta):
    """1 条件の走行。例外は呼び出し側で捕捉する（ここでは握り潰さない）。"""
    eng, p2, q2 = F.build_standard_universe(n, delta, Nn=NN, Neta=NETA, seed=SEED)
    syslr = F.abl.LowRankSystem(n)
    H = {k: np.full(T, np.nan) for k in KEYS}
    H_resid = np.full((T, 3), np.nan)
    H_acq = np.zeros(T, bool)
    cond_closure = np.full(T, np.nan)      # 凝縮体セル (2,0) の閉塞残差
    seed_closure = np.full(T, np.nan)      # シードセル (1,0) の閉塞残差
    prev_fr = None
    carry = {"C_flat": None, "c_gen": None}
    for t in range(T):
        eng.step()
        C2 = eng.C2()
        th = np.angle(np.sum(C2.reshape(C2.shape[0], -1), axis=1))
        syslr.set_theta(th)
        dp = D.d_panel(C2, syslr.kmatvec, p2, q2, cell=CELL, order=ORDER,
                       frame_prev=prev_fr)
        prev_fr = dp["_frame"]
        gp = G.g_panel(C2, p2, q2, carry["C_flat"], carry["c_gen"])
        carry = gp["_carry"]
        acq = S.s_clock_acquirable(gp)
        H["f2"][t] = gp["f2"]; H["f_seed"][t] = gp["f_seed"]
        H["align"][t] = dp["frame_align"]; H["n_eff"][t] = dp["ladder_n_eff"]
        H["nonunif"][t] = dp["gauge_nonunif"]; H["closure"][t] = dp["total_closure"]
        H["carrier_power"][t] = gp["carrier_power"]; H["coherence"][t] = gp["coherence"]
        H["omega_hat"][t] = gp["phase"]; H["phi_weight"][t] = dp["clock_phi_weight"]
        H_resid[t] = dp["frame_resid"]
        H_acq[t] = bool(acq["acquirable"])
        cond_closure[t] = dp["cell_closure"][CELL[0], CELL[1]]
        seed_closure[t] = dp["cell_closure"][1, 0]
    return H, H_resid, H_acq, cond_closure, seed_closure


def summarize(n, Hm, Rm, Am, Ccm, Csm, Hv, Av):
    w = slice(*WIN)
    def med(a):
        x = a[w]; x = x[np.isfinite(x)]
        return float(np.median(x)) if len(x) else float("nan")
    tau_space = first_true(np.nan_to_num(Hm["f2"], nan=-1) > 0.05)
    tau_space_v = first_true(np.nan_to_num(Hv["f2"], nan=-1) > 0.05)
    tau_time = first_true(Am)
    lock = np.abs(Hm["omega_hat"] - OMEGA_REF) / OMEGA_REF < 0.01
    tau_lock = first_true(np.nan_to_num(lock, nan=False).astype(bool))
    return {"N": n, "M": n * (n - 1) // 2, "built": True,
            "tau_space": tau_space, "tau_space_vacuum": tau_space_v,
            "tau_time": tau_time, "tau_lock": tau_lock,
            "space_born": tau_space is not None,
            "matter_born": bool(np.nanmax(Hm["f_seed"]) > 1e-30),
            "time_born": tau_time is not None,
            "vacuum_time_born": bool(Av.any()),
            "align_med": med(Hm["align"]), "n_eff_med": med(Hm["n_eff"]),
            "nonunif_med": med(Hm["nonunif"]), "closure_med": med(Hm["closure"]),
            "f_seed_med": med(Hm["f_seed"]), "carrier_med": med(Hm["carrier_power"]),
            "coherence_med": med(Hm["coherence"]), "omega_med": med(Hm["omega_hat"]),
            "cond_closure_med": med(Ccm), "seed_closure_med": med(Csm),
            "resid_med": [float(np.nanmedian(Rm[w, j])) for j in range(3)],
            "align_med_vacuum": med(Hv["align"]), "n_eff_med_vacuum": med(Hv["n_eff"]),
            "clock_rate": float(Am[w].mean()), "clock_rate_vacuum": float(Av[w].mean())}


def fig_one(n, Hm, Hv, Am, Ccm, Csm, rec):
    ts = np.arange(1, T + 1)
    fig, ax = plt.subplots(4, 1, figsize=(8.5, 10), sharex=True)
    ax[0].semilogy(ts, np.maximum(Hm["f2"], 1e-18), lw=0.8, label="物質")
    ax[0].semilogy(ts, np.maximum(Hv["f2"], 1e-18), "k--", lw=0.8, label="真空")
    ax[0].axhline(0.05, color="tab:blue", lw=0.7, ls=":")
    ax[0].set_ylabel("空間 f₂"); ax[0].legend(fontsize=7)
    ax[0].set_title(f"N={n}（M={n*(n-1)//2}）  空間 τ={rec['tau_space']}  "
                    f"時間 τ={rec['tau_time']}  時計定着 τ={rec['tau_lock']}")
    ax[1].plot(ts, Hm["n_eff"], lw=0.8, label="n_eff 物質")
    ax[1].plot(ts, Hv["n_eff"], "k--", lw=0.8, label="n_eff 真空")
    ax[1].plot(ts, Hm["align"], lw=0.8, color="tab:orange", label="align 物質")
    ax[1].set_ylabel("次元"); ax[1].legend(fontsize=7)
    ax[2].semilogy(ts, np.maximum(Hm["carrier_power"], 1e-30), lw=0.9,
                   color="tab:red", label="担い手 物質")
    ax[2].semilogy(ts, np.maximum(Hv["carrier_power"], 1e-30), "k--", lw=0.9,
                   label="担い手 真空")
    ax[2].fill_between(ts, 1e-30, 1e0, where=Am, color="tab:red", alpha=0.10,
                       label="時間が読める")
    ax[2].set_ylabel("時間"); ax[2].legend(fontsize=7)
    ax[3].semilogy(ts, np.maximum(Hm["closure"], 1e-20), lw=0.8, label="全体")
    ax[3].semilogy(ts, np.maximum(Ccm, 1e-20), lw=0.8, label="凝縮体セル(2,0)")
    ax[3].semilogy(ts, np.maximum(Csm, 1e-20), lw=0.8, label="シードセル(1,0)")
    ax[3].set_ylabel("閉塞残差"); ax[3].set_xlabel("τ（step）"); ax[3].legend(fontsize=7)
    for a in ax:
        if rec["tau_space"]:
            a.axvline(rec["tau_space"], color="tab:blue", lw=0.8, ls="-.")
        if rec["tau_time"]:
            a.axvline(rec["tau_time"], color="tab:red", lw=0.8, ls="-.")
        a.axvspan(WIN[0], WIN[1], color="green", alpha=0.06)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_nsweep_N{n}_v1.png", dpi=110)
    plt.close(fig)


def fig_summary(recs, fails):
    ns = [r["N"] for r in recs]
    fig, ax = plt.subplots(3, 2, figsize=(12, 11))
    a = ax[0, 0]
    a.plot(ns, [r["tau_space"] or np.nan for r in recs], "o-", label="空間 τ_space")
    a.plot(ns, [r["tau_time"] or np.nan for r in recs], "s-", label="時間 τ_time")
    a.plot(ns, [r["tau_lock"] or np.nan for r in recs], "^-", label="時計定着 τ_lock")
    a.set_xlabel("N"); a.set_ylabel("誕生時刻 τ（step）"); a.legend(fontsize=8)
    a.set_title("何がいつ生まれるか")
    a = ax[0, 1]
    a.plot(ns, [r["n_eff_med"] for r in recs], "o-", label="物質")
    a.plot(ns, [r["n_eff_med_vacuum"] for r in recs], "k--s", label="真空")
    a.axhline(1.0, color="red", lw=0.8, ls=":")
    a.set_xlabel("N"); a.set_ylabel("n_eff（準安定窓中央値）"); a.legend(fontsize=8)
    a.set_title("凝縮体の数＝平面の縮退")
    a = ax[1, 0]
    a.plot(ns, [r["align_med"] for r in recs], "o-", label="物質")
    a.plot(ns, [r["align_med_vacuum"] for r in recs], "k--s", label="真空")
    a.set_xlabel("N"); a.set_ylabel("第3次元の確定度"); a.legend(fontsize=8)
    a.set_title("次元は結晶化するか")
    a = ax[1, 1]
    ms = [r["M"] for r in recs]
    a.loglog(ms, [max(1 - r["align_med"], 1e-6) for r in recs], "o-", label="1−align")
    mm = np.array([m for m in ms if m > 0], float)
    if len(mm):
        a.loglog(mm, 1.0 / mm, "k:", label="1/M 参照線")
    a.set_xlabel("M = N(N−1)/2"); a.set_ylabel("曖昧さ 1−align"); a.legend(fontsize=8)
    a.set_title("曖昧さの M スケーリング")
    a = ax[2, 0]
    a.semilogy(ns, [max(r["cond_closure_med"], 1e-20) for r in recs], "o-",
               label="凝縮体セル(2,0)")
    a.semilogy(ns, [max(r["seed_closure_med"], 1e-20) for r in recs], "s-",
               label="シードセル(1,0)")
    a.semilogy(ns, [max(r["closure_med"], 1e-20) for r in recs], "^-", label="全体")
    a.set_xlabel("N"); a.set_ylabel("閉塞残差（中央値）"); a.legend(fontsize=8)
    a.set_title("ゼロ閉塞ブロック（凝縮体）はできるか")
    a = ax[2, 1]
    a.semilogy(ns, [max(r["carrier_med"], 1e-32) for r in recs], "o-", label="担い手")
    a.semilogy(ns, [max(r["f_seed_med"], 1e-32) for r in recs], "s-", label="物質分率")
    a.set_xlabel("N"); a.set_ylabel("量（中央値）"); a.legend(fontsize=8)
    a.set_title("物質と時計の担い手")
    if fails:
        fig.suptitle("構成できなかった N: " + ", ".join(str(f) for f in fails),
                     fontsize=10, y=0.995)
    fig.tight_layout()
    fig.savefig(HERE / "fig_nsweep_summary_v1.png", dpi=130)
    plt.close(fig)


def fig_matrix(recs, fails, nmin, nmax):
    labels = ["構成できた", "空間が生まれた", "物質が生まれた", "時間が生まれた",
              "3次元が確定(align>0.8)", "凝縮体(閉塞<1e-5)"]
    ns = list(range(nmin, nmax + 1))
    Mx = np.zeros((len(labels), len(ns)))
    byN = {r["N"]: r for r in recs}
    for j, n in enumerate(ns):
        r = byN.get(n)
        if r is None:
            Mx[:, j] = 0.0
            continue
        Mx[0, j] = 1.0
        Mx[1, j] = 1.0 if r["space_born"] else 0.0
        Mx[2, j] = 1.0 if r["matter_born"] else 0.0
        Mx[3, j] = 1.0 if r["time_born"] else 0.0
        Mx[4, j] = 1.0 if (np.isfinite(r["align_med"]) and r["align_med"] > 0.8) else 0.0
        Mx[5, j] = 1.0 if (np.isfinite(r["cond_closure_med"])
                           and r["cond_closure_med"] < 1e-5) else 0.0
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.imshow(Mx, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(ns))); ax.set_xticklabels(ns)
    ax.set_yticks(np.arange(len(labels))); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("N")
    ax.set_title("N ごとに何が生まれ、何が生まれないか（緑=生まれた・赤=生まれない）")
    for j, n in enumerate(ns):
        for i in range(len(labels)):
            ax.text(j, i, "○" if Mx[i, j] > 0.5 else "×", ha="center",
                    va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "fig_nsweep_birth_matrix_v1.png", dpi=130)
    plt.close(fig)


def main():
    t0 = time.time()
    nmin = int(sys.argv[1]) if len(sys.argv) > 2 else 1
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    print(f"=== N 掃引 {nmin}→{nmax}（Nn={NN}・Nη={NETA}・δ={DELTA}・"
          f"seed={SEED} 固定・T={T}）===")
    recs, fails, out = [], [], {"env": {
        "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA, "seed": SEED,
        "cell": list(CELL), "order": ORDER, "window": list(WIN),
        "functions": ["unified_interaction_v2", "unified_dimension_v1",
                      "unified_readout_v3", "selection_v1"]}, "N": {}, "failed": {}}
    for n in range(nmin, nmax + 1):
        t1 = time.time()
        try:
            Hm, Rm, Am, Ccm, Csm = run_one(n, DELTA)
            Hv, Rv, Av, Ccv, Csv = run_one(n, 0.0)
        except Exception as ex:
            msg = f"{type(ex).__name__}: {ex}"
            fails.append(n)
            out["failed"][n] = msg
            out["N"][n] = {"N": n, "M": n * (n - 1) // 2, "built": False,
                           "error": msg}
            print(f"N={n:3d} M={n*(n-1)//2:4d}: **構成不能** {msg[:70]}")
            (HERE / "result_tb_nsweep_1to20_v1.json").write_text(
                json.dumps(out, indent=1, ensure_ascii=False, default=float))
            continue
        rec = summarize(n, Hm, Rm, Am, Ccm, Csm, Hv, Av)
        recs.append(rec); out["N"][n] = rec
        fig_one(n, Hm, Hv, Am, Ccm, Csm, rec)
        np.savez_compressed(HERE / f"tb_nsweep_N{n}_v1.npz",
                            **{f"m_{k}": Hm[k] for k in KEYS},
                            **{f"v_{k}": Hv[k] for k in KEYS},
                            m_resid=Rm, m_acq=Am, v_acq=Av,
                            m_cond_closure=Ccm, m_seed_closure=Csm)
        print(f"N={n:3d} M={rec['M']:4d}: 空間τ={str(rec['tau_space']):>5} "
              f"物質={str(rec['matter_born']):>5} 時間τ={str(rec['tau_time']):>5} "
              f"時計定着τ={str(rec['tau_lock']):>5} "
              f"align={rec['align_med']:.4f} n_eff={rec['n_eff_med']:.3f} "
              f"凝縮体閉塞={rec['cond_closure_med']:.2e} [{time.time()-t1:.0f}s]")
        (HERE / "result_tb_nsweep_1to20_v1.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False, default=float))
    if recs:
        fig_summary(recs, fails)
    fig_matrix(recs, fails, nmin, nmax)
    out["failed_N"] = fails
    out["runtime_sec"] = time.time() - t0
    (HERE / "result_tb_nsweep_1to20_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"\n構成できなかった N: {fails if fails else 'なし'}")
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
