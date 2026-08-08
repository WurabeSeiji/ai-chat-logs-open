#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""万能次元読出関数 D の単体テスト v2（N=12・τ=0→準安定・空間と時間の誕生を図化）

v1 からの変更（2026-08-08・木原指摘）: v1 は**空間側の量しか記録・図化して
おらず、「時間がいつから読めるようになるか」がどの図でも読めなかった**。
時計の担い手（質量）・時計の鋭さ・取得可否・ω̂ の定着を記録し、
**空間の誕生と時間の誕生を同じ τ 軸で対比する**パネルを追加した。
v1 は再現性のため残置する。

使用する万能関数（2026-08-08 分岐版の3本）:
  F = ../統一万能関数_v1/unified_interaction_v2.py
  D = ../統一万能関数_v1/unified_dimension_v1.py
  G = ../統一万能関数_v1/unified_readout_v3.py
  S = ../統一万能関数_v1/selection_v1.py（取得可否の判定＝宣言した床を適用）

条件（すべて宣言値）:
  N=12（M=66）・帯 Nn=16・毛 Nη=8・T=4000
  ニュートリノ型シード: 帯 k=1・巻き η=0・振幅 δ=10⁻²・親 seed=2
  ポンプ（真空）: 帯 k=2・巻き η=0
  D の宣言引数: cell=(2,0)・order=6
  準安定窓: τ∈[2000,4000]
  選択の宣言: s_clock_acquirable（既定床 FLOOR_OVERLAP=1e-30・FLOOR_CARRIER=1e-12）
  対照: 真空宇宙 δ=0

事前登録した判定（実行前固定・閾値はすべて本スクリプト側）:
 (U1) 準安定窓で D の主要出力（frame_resid・align・n_eff・gauge_nonunif）が有限。
 (U2) 準安定窓で n_eff の変動係数 < 0.10。
 (U3) 凝縮体セル(2,0) の閉塞残差がシードセル(1,0) より 3 桁以上小さい。
 (U4) 準安定窓の align 中央値 > 潜伏期の中央値。
 (U5) 真空対照: 空間は形成される（f₂>0.05 に到達）が物質は生じない。
 (U6) **時間の誕生**: 物質宇宙で s_clock_acquirable が True になる時刻 τ_time が
      存在し、真空宇宙では全 τ で False であること。
 (U7) **時計の定着**: 物質宇宙で |ω̂ − π/72|/(π/72) < 0.01 となる時刻が存在する。
 記録（判定でなく取得）: 空間の誕生 τ_space（f₂>0.05）と τ_time の前後関係。

図（横軸はすべて τ（step））:
 図A fig_dim_unittest_N12_v2.png        D の全出力＋**時間の誕生**（6段）
 図B fig_dim_condensate_birth_N12_v2.png 凝縮体はどこで・いくつ生まれるか
 図C fig_dim_spacetime_birth_N12_v2.png  **空間の誕生 vs 時間の誕生**（対比専用）

使い方: python3 run_tb_dimension_unittest_N12_v2.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
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


F = load("f_dut2", UF / "unified_interaction_v2.py")
D = load("d_dut2", UF / "unified_dimension_v1.py")
G = load("g_dut2", UF / "unified_readout_v3.py")
S = load("s_dut2", UF / "selection_v1.py")

N, NN, NETA = 12, 16, 8
T = 4000
DELTA = 1e-2
SEED = 2
CELL = (2, 0)
ORDER = 6
WIN = (2000, 4000)
OMEGA_REF = np.pi / 72.0
TRACK = [(2, 0), (1, 0), (3, 0), (0, 0), (4, 0)]
KEYS = ("align", "n_eff", "nonunif", "closure", "f2", "f_seed",
        "phi_weight", "axis_persist", "omega_gen",
        "carrier_power", "coherence", "omega_hat")


def run(delta):
    eng, p2, q2 = F.build_standard_universe(N, delta, Nn=NN, Neta=NETA, seed=SEED)
    syslr = F.abl.LowRankSystem(N)
    H = {k: np.zeros(T) for k in KEYS}
    H_resid = np.zeros((T, 3))
    H_cell = np.zeros((T, len(TRACK)))
    H_wts = np.zeros((T, ORDER))
    H_acq = np.zeros(T, bool)
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
        # 【宣言】選択 S = s_clock_acquirable（既定床）で時間の取得可否を判定する
        acq = S.s_clock_acquirable(gp)
        H["align"][t] = dp["frame_align"]
        H["n_eff"][t] = dp["ladder_n_eff"]
        H["nonunif"][t] = dp["gauge_nonunif"]
        H["closure"][t] = dp["total_closure"]
        H["phi_weight"][t] = dp["clock_phi_weight"]
        H["axis_persist"][t] = dp["pers_axis_persist"]
        H["omega_gen"][t] = dp["frame_omega_gen"]
        H["f2"][t] = gp["f2"]
        H["f_seed"][t] = gp["f_seed"]
        H["carrier_power"][t] = gp["carrier_power"]
        H["coherence"][t] = gp["coherence"]
        H["omega_hat"][t] = gp["phase"]
        H_acq[t] = bool(acq["acquirable"])
        H_resid[t] = dp["frame_resid"]
        H_cell[t] = [dp["cell_closure"][k, e] for (k, e) in TRACK]
        w = dp["ladder_weights"]
        sw = w.sum() if len(w) else 0.0
        H_wts[t, :min(ORDER, len(w))] = (w[:ORDER] / sw) if sw > 0 else 0.0
    return H, H_resid, H_cell, H_wts, H_acq


def first_true(mask):
    idx = np.flatnonzero(mask)
    return int(idx[0]) + 1 if len(idx) else None


def main():
    t0 = time.time()
    print(f"=== D 単体テスト v2（N={N}・Nn={NN}・δ={DELTA}・seed={SEED}・"
          f"cell={CELL}・order={ORDER}）===")
    Hm, Rm, Cm, Wm, Am = run(DELTA)
    print(f"  物質宇宙 走行完了 [{time.time()-t0:.0f}s]")
    Hv, Rv, Cv, Wv, Av = run(0.0)
    print(f"  真空宇宙 走行完了 [{time.time()-t0:.0f}s]")

    w = slice(*WIN)
    tau_space = first_true(Hm["f2"] > 0.05)
    tau_space_v = first_true(Hv["f2"] > 0.05)
    lat = slice(0, max((tau_space or 2) - 1, 1))

    fin = {k: bool(np.all(np.isfinite(Hm[k][w]))) for k in
           ("align", "n_eff", "nonunif")}
    fin["resid"] = bool(np.all(np.isfinite(Rm[w])))
    U1 = all(fin.values())
    ne = Hm["n_eff"][w]
    cv = float(np.std(ne) / np.mean(ne))
    U2 = cv < 0.10
    c_cond = float(np.median(Cm[w, 0])); c_seed = float(np.median(Cm[w, 1]))
    ratio = c_seed / c_cond
    U3 = ratio > 1e3
    a_lat = float(np.median(Hm["align"][lat])); a_win = float(np.median(Hm["align"][w]))
    U4 = a_win > a_lat
    U5 = bool(tau_space_v is not None and float(Hv["f_seed"][-1]) == 0.0)
    tau_time = first_true(Am)
    U6 = bool(tau_time is not None and not Av.any())
    lock = np.abs(Hm["omega_hat"] - OMEGA_REF) / OMEGA_REF < 0.01
    tau_lock = first_true(lock)
    U7 = tau_lock is not None

    print(f"\n(U1) 準安定窓での可読性: {fin} → {'通過' if U1 else '不成立'}")
    print(f"(U2) n_eff の収束: 中央={np.median(ne):.4f} 変動係数={cv:.4f} → "
          f"{'通過' if U2 else '不成立'}")
    print(f"(U3) 凝縮体の検出: (2,0)={c_cond:.3e} / (1,0)={c_seed:.3e} "
          f"比={ratio:.2e} → {'通過' if U3 else '不成立'}")
    print(f"(U4) 結晶化: align 潜伏期={a_lat:.4f} → 準安定={a_win:.4f} → "
          f"{'通過' if U4 else '不成立'}")
    print(f"(U5) 真空対照: crossing={tau_space_v}・f_seed 終端"
          f"={Hv['f_seed'][-1]:.3e} → {'通過' if U5 else '不成立'}")
    print(f"(U6) **時間の誕生**: 物質 τ_time={tau_time}（時計が取得可能になる時刻）"
          f"・真空は全 τ で取得不能={not Av.any()} → {'通過' if U6 else '不成立'}")
    print(f"(U7) 時計の定着: |ω̂−π/72|/(π/72)<0.01 の初回 τ={tau_lock} → "
          f"{'通過' if U7 else '不成立'}")
    print(f"\n--- 空間と時間の誕生 ---")
    print(f"  空間の誕生 τ_space（f₂>0.05）: 物質={tau_space}  真空={tau_space_v}")
    print(f"  時間の誕生 τ_time（時計取得可）: 物質={tau_time}  真空=なし")
    print(f"  時計の定着 τ_lock（ω̂→π/72）  : 物質={tau_lock}")
    print(f"  n_eff が 1.5 を下回る時刻（縮退解消）: 物質="
          f"{first_true(Hm['n_eff'] < 1.5)}  真空={first_true(Hv['n_eff'] < 1.5)}")
    print(f"\n--- 準安定窓の中央値 ---")
    for k in ("align", "n_eff", "nonunif", "closure", "carrier_power",
              "coherence", "omega_hat", "phi_weight"):
        print(f"  {k:14s} 物質={np.median(Hm[k][w]):.6e}  "
              f"真空={np.median(Hv[k][w]):.6e}")
    print(f"  時計取得率（準安定窓）: 物質={Am[w].mean():.3f}  真空={Av[w].mean():.3f}")

    out = {"env": {"N": N, "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA,
                   "seed": SEED, "cell": list(CELL), "order": ORDER,
                   "window": list(WIN), "omega_ref": OMEGA_REF,
                   "functions": ["unified_interaction_v2", "unified_dimension_v1",
                                 "unified_readout_v3", "selection_v1"],
                   "selection": "s_clock_acquirable（既定床）"},
           "tau_space_matter": tau_space, "tau_space_vacuum": tau_space_v,
           "tau_time_matter": tau_time, "tau_time_vacuum": None,
           "tau_clock_lock": tau_lock,
           "tau_neff_below_1p5": {"matter": first_true(Hm["n_eff"] < 1.5),
                                  "vacuum": first_true(Hv["n_eff"] < 1.5)},
           "U1": bool(U1), "U2": bool(U2), "U3": bool(U3), "U4": bool(U4),
           "U5": bool(U5), "U6": bool(U6), "U7": bool(U7),
           "detail": {"fin": fin, "n_eff_cv": cv, "closure_ratio": ratio,
                      "align_lat": a_lat, "align_win": a_win},
           "window_median": {k: float(np.median(Hm[k][w])) for k in KEYS},
           "window_median_vacuum": {k: float(np.median(Hv[k][w])) for k in KEYS},
           "clock_rate_window": {"matter": float(Am[w].mean()),
                                 "vacuum": float(Av[w].mean())},
           "frame_resid_median": [float(x) for x in np.median(Rm[w], axis=0)],
           "plane_weights_median": [float(x) for x in np.median(Wm[w], axis=0)],
           "cell_closure_median": {f"{k}_{e}": float(np.median(Cm[w, i]))
                                   for i, (k, e) in enumerate(TRACK)},
           "runtime_sec": time.time() - t0}
    (HERE / "result_tb_dimension_unittest_N12_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    np.savez_compressed(HERE / "tb_dimension_unittest_N12_v2.npz",
                        **{f"m_{k}": Hm[k] for k in KEYS},
                        **{f"v_{k}": Hv[k] for k in KEYS},
                        m_resid=Rm, v_resid=Rv, m_cell=Cm, v_cell=Cv,
                        m_wts=Wm, v_wts=Wv, m_acq=Am, v_acq=Av)

    ts = np.arange(1, T + 1)

    def marks(a):
        a.axvline(tau_space, color="tab:blue", lw=0.9, ls="-.",
                  label=f"空間の誕生 τ={tau_space}")
        if tau_time:
            a.axvline(tau_time, color="tab:red", lw=0.9, ls="-.",
                      label=f"時間の誕生 τ={tau_time}")
        a.axvspan(WIN[0], WIN[1], color="green", alpha=0.06)

    # ---- 図A: D の全出力＋時間の誕生
    fig, ax = plt.subplots(6, 1, figsize=(9, 15), sharex=True)
    ax[0].semilogy(ts, np.maximum(Hm["closure"], 1e-20), lw=0.8, label="物質 全体")
    ax[0].semilogy(ts, np.maximum(Hv["closure"], 1e-20), "k--", lw=0.8, label="真空 全体")
    for i, (k, e) in enumerate(TRACK):
        ax[0].semilogy(ts, np.maximum(Cm[:, i], 1e-20), lw=0.6, alpha=0.7,
                       label=f"物質 セル(k={k})")
    ax[0].set_ylabel("閉塞残差"); ax[0].legend(fontsize=6, ncol=3)
    ax[0].set_title(f"D 単体テスト v2  N={N}・Nn={NN}・δ={DELTA}（F v2 / D v1 / G v3）")
    for j in range(3):
        ax[1].semilogy(ts, np.maximum(Rm[:, j], 1e-20), lw=0.8,
                       label=["位置", "速度", "加速度"][j])
    ax[1].set_ylabel("直交化残差"); ax[1].legend(fontsize=8)
    ax[2].plot(ts, Hm["align"], lw=0.8, label="物質")
    ax[2].plot(ts, Hv["align"], "k--", lw=0.8, label="真空")
    ax[2].set_ylabel("第3次元の確定度"); ax[2].legend(fontsize=7)
    ax[3].plot(ts, Hm["n_eff"], lw=0.8, label="物質")
    ax[3].plot(ts, Hv["n_eff"], "k--", lw=0.8, label="真空")
    ax[3].axhline(1.0, color="red", lw=0.8, ls=":")
    ax[3].set_ylabel("実効平面数 n_eff"); ax[3].legend(fontsize=7)
    ax[4].plot(ts, Hm["nonunif"], lw=0.8, label="物質")
    ax[4].plot(ts, Hv["nonunif"], "k--", lw=0.8, label="真空")
    ax[4].set_ylabel("ゲージ非一様度（＝重力）"); ax[4].legend(fontsize=7)
    ax[5].semilogy(ts, np.maximum(Hm["carrier_power"], 1e-30), lw=0.9,
                   color="tab:red", label="物質 時計の担い手パワー")
    ax[5].semilogy(ts, np.maximum(Hv["carrier_power"], 1e-30), "k--", lw=0.9,
                   label="真空 時計の担い手パワー（不在）")
    ax[5].fill_between(ts, 1e-30, 1e0, where=Am, color="tab:red", alpha=0.10,
                       label="時間が読める区間（物質）")
    ax[5].set_ylabel("時間の誕生"); ax[5].set_xlabel("τ（step）")
    ax[5].legend(fontsize=7)
    for a in ax:
        marks(a)
    ax[0].legend(fontsize=6, ncol=3)
    fig.tight_layout()
    fig.savefig(HERE / "fig_dim_unittest_N12_v2.png", dpi=130)
    plt.close(fig)

    # ---- 図B: 凝縮体はどこで・いくつ
    fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    ax[0].plot(ts, Hm["n_eff"], lw=0.9, label="物質")
    ax[0].plot(ts, Hv["n_eff"], "k--", lw=0.9, label="真空")
    ax[0].axhline(1.0, color="red", lw=0.8, ls=":")
    ax[0].set_ylabel("いくつ: 実効平面数 n_eff")
    ax[0].set_title(f"凝縮体はどこで・いくつ生まれるか（N={N}・Nn={NN}）")
    im = ax[1].pcolormesh(ts, np.arange(len(TRACK)),
                          np.log10(np.maximum(Cm.T, 1e-20)),
                          shading="auto", cmap="viridis")
    ax[1].set_yticks(np.arange(len(TRACK)))
    ax[1].set_yticklabels([f"k={k},η={e}" for (k, e) in TRACK], fontsize=8)
    ax[1].set_ylabel("どこで: セル別の閉塞残差")
    fig.colorbar(im, ax=ax[1], label="log₁₀ 閉塞残差（小＝ゼロ閉塞ブロック）")
    for j in range(min(4, ORDER)):
        ax[2].plot(ts, Wm[:, j], lw=0.8, label=f"平面{j+1}")
    ax[2].set_ylabel("平面の重み"); ax[2].set_xlabel("τ（step）")
    ax[2].legend(fontsize=8)
    for a in ax:
        marks(a)
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "fig_dim_condensate_birth_N12_v2.png", dpi=130)
    plt.close(fig)

    # ---- 図C: 空間の誕生 vs 時間の誕生（対比専用）
    fig, ax = plt.subplots(4, 1, figsize=(9, 11), sharex=True)
    ax[0].semilogy(ts, np.maximum(Hm["f2"], 1e-18), lw=0.9, label="物質")
    ax[0].semilogy(ts, np.maximum(Hv["f2"], 1e-18), "k--", lw=0.9, label="真空")
    ax[0].axhline(0.05, color="tab:blue", lw=0.8, ls=":", label="crossing 判定 0.05")
    ax[0].set_ylabel("空間: f₂（空間形成史）"); ax[0].legend(fontsize=7)
    ax[0].set_title("空間の誕生と時間の誕生（同一 τ 軸・N=12）")
    ax[1].plot(ts, Hm["n_eff"], lw=0.9, label="物質")
    ax[1].plot(ts, Hv["n_eff"], "k--", lw=0.9, label="真空")
    ax[1].set_ylabel("空間: 縮退の解消 n_eff"); ax[1].legend(fontsize=7)
    ax[2].semilogy(ts, np.maximum(Hm["carrier_power"], 1e-30), lw=0.9,
                   color="tab:red", label="物質")
    ax[2].semilogy(ts, np.maximum(Hv["carrier_power"], 1e-30), "k--", lw=0.9,
                   label="真空（担い手なし）")
    ax[2].fill_between(ts, 1e-30, 1e0, where=Am, color="tab:red", alpha=0.10,
                       label="時間が読める区間")
    ax[2].set_ylabel("時間: 時計の担い手パワー"); ax[2].legend(fontsize=7)
    ax[3].plot(ts, Hm["omega_hat"], ".", ms=1, color="tab:red", label="物質 ω̂")
    ax[3].plot(ts, Hv["omega_hat"], ".", ms=1, color="gray", label="真空 ω̂")
    ax[3].axhline(OMEGA_REF, color="black", lw=0.9, ls="--", label="π/72")
    if tau_lock:
        ax[3].axvline(tau_lock, color="tab:green", lw=0.9, ls="-.",
                      label=f"時計の定着 τ={tau_lock}")
    ax[3].set_ylabel("時間: 物質時計 ω̂"); ax[3].set_xlabel("τ（step）")
    ax[3].legend(fontsize=7)
    for a in ax:
        marks(a)
    fig.tight_layout()
    fig.savefig(HERE / "fig_dim_spacetime_birth_N12_v2.png", dpi=130)
    plt.close(fig)
    print(f"\n完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
