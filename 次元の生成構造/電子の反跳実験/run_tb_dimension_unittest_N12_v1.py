#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""万能次元読出関数 D の単体テスト（N=12・τ=0→準安定・図化つき）

目的: 特定条件（N=12・ニュートリノ型シード）で通し創世させ、インフレーション後の
準安定状況で **D が正しく次元を読み出せるか**を単体で確かめる。条件スキャンは
目的ではない。

使用する万能関数（2026-08-08 分岐版の3本）:
  F = ../統一万能関数_v1/unified_interaction_v2.py
  D = ../統一万能関数_v1/unified_dimension_v1.py
  G = ../統一万能関数_v1/unified_readout_v3.py
  （閾値・IF による代用は3本とも撤廃済み。判定と停止条件は本スクリプト側の責務）

条件（すべて宣言値）:
  N=12（M=66）・帯 Nn=16・毛 Nη=8・T=4000
  ニュートリノ型シード: 帯 k=1・巻き η=0・振幅 δ=10⁻²・親 seed=2
  ポンプ（真空）: 帯 k=2・巻き η=0
  D の宣言引数: cell=(2,0)（凝縮体セルからフレームを作る）・order=6（Krylov 打切り）
  準安定窓: τ∈[2000,4000]
  対照: 真空宇宙 δ=0（F v2 では頂点スキップ分岐を撤廃したため v1 と軌道が
        一致しない見込み——本ランで実測する）

事前登録した判定（実行前固定）:
 (U1) 可読性: 準安定窓の全ステップで D の主要出力（frame_resid・align・n_eff・
      gauge_nonunif）が有限であること（NaN は「読めない」の正当な表現だが、
      準安定窓では読めていなければならない）。
 (U2) 収束: 準安定窓で n_eff の変動係数 std/mean < 0.10。
 (U3) 凝縮体の検出: 凝縮体セル (k=2,η=0) の閉塞残差が、シードセル (k=1,η=0) の
      それより 3 桁以上小さいこと（ゼロ閉塞ブロックとして分離できるか）。
 (U4) 結晶化: 準安定窓の align 中央値が潜伏期（τ<crossing）の中央値より大きい。
 (U5) 真空対照: 真空宇宙でも空間は形成される（f₂ が立つ）が、物質は生じない。
      F v2 と v1 の軌道差の有無を記録（判定ではなく取得）。
判定に用いる閾値は本スクリプト側に置く（D には持たせない・規約 R9）。

図（横軸はすべて τ（step）——処理ステップであり座標時間 t ではない）:
 図A fig_dim_unittest_N12_v1.png     D の全出力の推移（5段）
 図B fig_dim_condensate_birth_N12_v1.png  凝縮体は「どこで・いくつ」生まれるか

使い方: python3 run_tb_dimension_unittest_N12_v1.py
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


F = load("f_dut", UF / "unified_interaction_v2.py")
D = load("d_dut", UF / "unified_dimension_v1.py")
G = load("g_dut", UF / "unified_readout_v3.py")

N, NN, NETA = 12, 16, 8
T = 4000
DELTA = 1e-2
SEED = 2
CELL = (2, 0)          # フレームを作るセル（宣言）
ORDER = 6              # Krylov 打切り次数（宣言）
WIN = (2000, 4000)     # 準安定窓（宣言）
TRACK = [(2, 0), (1, 0), (3, 0), (0, 0), (4, 0)]
KEYS = ("align", "n_eff", "nonunif", "closure", "f2", "f_seed",
        "phi_weight", "axis_persist", "omega_gen")


def run(delta):
    eng, p2, q2 = F.build_standard_universe(N, delta, Nn=NN, Neta=NETA, seed=SEED)
    syslr = F.abl.LowRankSystem(N)
    H = {k: np.zeros(T) for k in KEYS}
    H_resid = np.zeros((T, 3))
    H_cell = np.zeros((T, len(TRACK)))
    H_wts = np.zeros((T, ORDER))
    prev_fr = None
    carry = {"C_flat": None, "c_gen": None}
    for t in range(T):
        eng.step()
        C2 = eng.C2()
        # 生成子（正本 LowRankSystem・read-only）を現ステップの位相で構成
        th = np.angle(np.sum(C2.reshape(C2.shape[0], -1), axis=1))
        syslr.set_theta(th)
        dp = D.d_panel(C2, syslr.kmatvec, p2, q2, cell=CELL, order=ORDER,
                       frame_prev=prev_fr)
        prev_fr = dp["_frame"]
        gp = G.g_panel(C2, p2, q2, carry["C_flat"], carry["c_gen"])
        carry = gp["_carry"]
        H["align"][t] = dp["frame_align"]
        H["n_eff"][t] = dp["ladder_n_eff"]
        H["nonunif"][t] = dp["gauge_nonunif"]
        H["closure"][t] = dp["total_closure"]
        H["phi_weight"][t] = dp["clock_phi_weight"]
        H["axis_persist"][t] = dp["pers_axis_persist"]
        H["omega_gen"][t] = dp["frame_omega_gen"]
        H["f2"][t] = gp["f2"]
        H["f_seed"][t] = gp["f_seed"]
        H_resid[t] = dp["frame_resid"]
        H_cell[t] = [dp["cell_closure"][k, e] for (k, e) in TRACK]
        w = dp["ladder_weights"]
        H_wts[t, :min(ORDER, len(w))] = w[:ORDER] / (w.sum() if w.sum() > 0 else 1)
    return H, H_resid, H_cell, H_wts, eng


def main():
    t0 = time.time()
    print(f"=== D 単体テスト（N={N}・Nn={NN}・δ={DELTA}・seed={SEED}・"
          f"cell={CELL}・order={ORDER}）===")
    Hm, Rm, Cm, Wm, engm = run(DELTA)
    print(f"  物質宇宙 走行完了 [{time.time()-t0:.0f}s]")
    Hv, Rv, Cv, Wv, engv = run(0.0)
    print(f"  真空宇宙 走行完了 [{time.time()-t0:.0f}s]")

    w = slice(*WIN)
    cross = int(np.argmax(Hm["f2"] > 0.05)) + 1
    lat = slice(0, max(cross - 1, 1))

    # (U1) 準安定窓での可読性
    fin = {k: bool(np.all(np.isfinite(Hm[k][w]))) for k in
           ("align", "n_eff", "nonunif")}
    fin["resid"] = bool(np.all(np.isfinite(Rm[w])))
    U1 = all(fin.values())
    # (U2) n_eff の収束
    ne = Hm["n_eff"][w]
    cv = float(np.std(ne) / np.mean(ne))
    U2 = cv < 0.10
    # (U3) 凝縮体セルの分離
    c_cond = float(np.median(Cm[w, 0]))
    c_seed = float(np.median(Cm[w, 1]))
    ratio = c_seed / c_cond
    U3 = ratio > 1e3
    # (U4) 結晶化
    a_lat = float(np.median(Hm["align"][lat]))
    a_win = float(np.median(Hm["align"][w]))
    U4 = a_win > a_lat
    # (U5) 真空対照
    cross_v = int(np.argmax(Hv["f2"] > 0.05)) + 1
    fs_v = float(Hv["f_seed"][-1])

    print(f"\n(U1) 準安定窓での可読性: {fin} → {'通過' if U1 else '不成立'}")
    print(f"(U2) n_eff の収束: 中央={np.median(ne):.4f} 変動係数={cv:.4f}<0.10 → "
          f"{'通過' if U2 else '不成立'}")
    print(f"(U3) 凝縮体の検出: 凝縮体セル(2,0)={c_cond:.3e} / "
          f"シードセル(1,0)={c_seed:.3e} 比={ratio:.2e}>1e3 → "
          f"{'通過' if U3 else '不成立'}")
    print(f"(U4) 結晶化: align 潜伏期={a_lat:.4f} → 準安定={a_win:.4f} → "
          f"{'通過' if U4 else '不成立'}")
    print(f"(U5) 真空対照: crossing 物質={cross} 真空={cross_v}・"
          f"真空の f_seed 終端={fs_v:.3e}（物質は {Hm['f_seed'][-1]:.3e}）")

    print(f"\n--- 準安定窓の D 出力（中央値）---")
    for k in ("align", "n_eff", "nonunif", "closure", "omega_gen",
              "axis_persist", "phi_weight"):
        print(f"  {k:14s} 物質={np.median(Hm[k][w]):.6f}  "
              f"真空={np.median(Hv[k][w]):.6f}")
    print(f"  frame_resid    物質={np.round(np.median(Rm[w], axis=0), 4)}  "
          f"真空={np.round(np.median(Rv[w], axis=0), 4)}")
    print(f"  平面の重み(上位3・物質)={np.round(np.median(Wm[w], axis=0)[:3], 4)}")

    out = {"env": {"N": N, "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA,
                   "seed": SEED, "cell": list(CELL), "order": ORDER,
                   "window": list(WIN),
                   "functions": ["unified_interaction_v2", "unified_dimension_v1",
                                 "unified_readout_v3"]},
           "crossing_matter": cross, "crossing_vacuum": cross_v,
           "U1_readable": {"detail": fin, "ok": bool(U1)},
           "U2_converged": {"n_eff_med": float(np.median(ne)), "cv": cv,
                            "ok": bool(U2)},
           "U3_condensate": {"cell_2_0": c_cond, "cell_1_0": c_seed,
                             "ratio": ratio, "ok": bool(U3)},
           "U4_crystallize": {"align_latency": a_lat, "align_window": a_win,
                              "ok": bool(U4)},
           "U5_vacuum": {"crossing": cross_v, "f_seed_final": fs_v},
           "window_median": {k: float(np.median(Hm[k][w])) for k in KEYS},
           "window_median_vacuum": {k: float(np.median(Hv[k][w])) for k in KEYS},
           "frame_resid_median": [float(x) for x in np.median(Rm[w], axis=0)],
           "plane_weights_median": [float(x) for x in np.median(Wm[w], axis=0)],
           "cell_closure_median": {f"{k}_{e}": float(np.median(Cm[w, i]))
                                   for i, (k, e) in enumerate(TRACK)},
           "runtime_sec": time.time() - t0}
    (HERE / "result_tb_dimension_unittest_N12_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    np.savez_compressed(HERE / "tb_dimension_unittest_N12_v1.npz",
                        **{f"m_{k}": Hm[k] for k in KEYS},
                        **{f"v_{k}": Hv[k] for k in KEYS},
                        m_resid=Rm, v_resid=Rv, m_cell=Cm, v_cell=Cv,
                        m_wts=Wm, v_wts=Wv)

    ts = np.arange(1, T + 1)
    # ---- 図A: D の全出力
    fig, ax = plt.subplots(5, 1, figsize=(9, 13), sharex=True)
    ax[0].semilogy(ts, np.maximum(Hm["closure"], 1e-20), lw=0.8, label="物質 全体")
    ax[0].semilogy(ts, np.maximum(Hv["closure"], 1e-20), "k--", lw=0.8,
                   label="真空 全体")
    for i, (k, e) in enumerate(TRACK):
        ax[0].semilogy(ts, np.maximum(Cm[:, i], 1e-20), lw=0.6, alpha=0.7,
                       label=f"物質 セル(k={k})")
    ax[0].set_ylabel("閉塞残差 |Σxₙ²|/Σ|xₙ|²")
    ax[0].legend(fontsize=6, ncol=3)
    ax[0].set_title(f"D 単体テスト N={N}・Nn={NN}・δ={DELTA}（F v2 / D v1 / G v3）")
    for j in range(3):
        ax[1].semilogy(ts, np.maximum(Rm[:, j], 1e-20), lw=0.8,
                       label=["位置", "速度", "加速度"][j])
    ax[1].set_ylabel("直交化残差 resid"); ax[1].legend(fontsize=8)
    ax[2].plot(ts, Hm["align"], lw=0.8, label="物質")
    ax[2].plot(ts, Hv["align"], "k--", lw=0.8, label="真空")
    ax[2].set_ylabel("第3次元の確定度 |n̂·â|"); ax[2].legend(fontsize=8)
    ax[3].plot(ts, Hm["n_eff"], lw=0.8, label="物質")
    ax[3].plot(ts, Hv["n_eff"], "k--", lw=0.8, label="真空")
    ax[3].axhline(1.0, color="red", lw=0.8, ls=":", label="1枚（縮退なし）")
    ax[3].set_ylabel("実効平面数 n_eff（凝縮体の数）"); ax[3].legend(fontsize=8)
    ax[4].plot(ts, Hm["nonunif"], lw=0.8, label="物質")
    ax[4].plot(ts, Hv["nonunif"], "k--", lw=0.8, label="真空")
    ax[4].set_ylabel("ゲージ非一様度（＝重力）"); ax[4].set_xlabel("τ（step）")
    ax[4].legend(fontsize=8)
    for a in ax:
        a.axvline(cross, color="gray", lw=0.8)
        a.axvspan(WIN[0], WIN[1], color="green", alpha=0.06)
    fig.tight_layout()
    fig.savefig(HERE / "fig_dim_unittest_N12_v1.png", dpi=130)
    plt.close(fig)

    # ---- 図B: 凝縮体は「どこで・いくつ」生まれるか
    fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    ax[0].plot(ts, Hm["n_eff"], lw=0.9, label="物質")
    ax[0].plot(ts, Hv["n_eff"], "k--", lw=0.9, label="真空")
    ax[0].axhline(1.0, color="red", lw=0.8, ls=":")
    ax[0].set_ylabel("いくつ: 実効平面数 n_eff"); ax[0].legend(fontsize=8)
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
    ax[2].set_ylabel("平面の重み（規格化）"); ax[2].set_xlabel("τ（step）")
    ax[2].legend(fontsize=8)
    for a in ax:
        a.axvline(cross, color="gray", lw=0.8)
        a.axvspan(WIN[0], WIN[1], color="green", alpha=0.06)
    fig.tight_layout()
    fig.savefig(HERE / "fig_dim_condensate_birth_N12_v1.png", dpi=130)
    plt.close(fig)
    print(f"\n完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
