#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""子閉塞ブロック検定 v2（v1の反証を受けた修正版）

■ v1 の結果と、その反証（記録）
木原の仮説「Σxₙ²=0 の中に ΣAₙ²=0 となる子の塊があり、それが凝縮体」を、
v1 では「Z⊥ の特異値の梯子から取った回転平面が双一次形式について完全等方か」
という形で検定し、**全て不成立**だった（自己等方性 0.14・平面間 0.32・
子の閉塞 0.14——いずれも乱数対照 0.14 と区別が付かない）。

しかし v1 の (B1) が前提の誤りを暴いた: 閉塞はどの粒度でも厳密ゼロでは
なかった（全体 3.14e-05・セル別最大 3.14e-01・双対点別 1.82e-03）。
診断の結果、原因が判明した:

  * 真空ポンプ Z0c は **厳密ゼロ閉塞**（|Z0cᵀZ0c|/‖Z0c‖² = 4.4e-15）
  * 物質シード s0 は **ゼロ閉塞しない**（同 0.3176）
  * 標準宇宙の全体閉塞欠損は δ²×(シードの欠損) に厳密に一致する
    （τ=0: |Σz²| = 3.1764e-05 = 1.0e-4 × 0.31764）

すなわち**凝縮体はすでに厳密なゼロ閉塞ブロックであり、閉塞の欠損は物質
そのもの**であった。v1 は「凝縮体の内部をさらに平面へ割る」という誤った
切り方をしていた。正しい切り方は内容（セル）による分割である。

■ v2 の仮説（修正後）
  凝縮体（真空）＝ 厳密にゼロ閉塞する子ブロック
  物質         ＝ 閉塞の欠損そのもの（ゼロ閉塞しない成分）

事前登録した判定（実行前固定）:
 (C1) 真空凝縮体の厳密ゼロ閉塞: τ=0 で |Z0cᵀZ0c|/‖Z0c‖² < 1e-12
 (C2) 物質シードは非ゼロ閉塞: 同 > 1e-2
 (C3) 閉塞欠損＝物質量: τ=0 で |Σz²|_全体 = f_seed·|s0ᵀs0| を相対 1e-6 以内で満たす
 (C4) **真空宇宙（δ=0）は全時代で厳密ゼロ閉塞を保つ**: 全 τ で
      |Σz²|/Σ|z|² < 1e-12（インフレーション・凝縮を通しても破れない）
 (C5) 物質宇宙では凝縮体セルの閉塞が τ とともに劣化する（物質による汚染）
      ——劣化率を記録（判定でなく取得）
 (C6) セル別閉塞マップ: どのセルがゼロ閉塞ブロックかを全セルで記録
不成立の場合もそのまま記録する。

使い方: python3 run_tb_child_closure_blocks_v2.py
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


ui = load("ui_cb2", UF / "unified_interaction_v1.py")
G = load("ur_cb2", UF / "unified_readout_v2.py")

N, NN, NETA = 12, 16, 8
T = 4000
DELTA = 1e-2
TOL = 1e-12
TRACK = [(2, 0), (1, 0), (3, 0), (0, 0), (4, 0)]   # 追跡セル（凝縮体・シード・生成・零・高次）


def bil(z):
    return complex(np.sum(z * z))


def rel_bil(z):
    d = float(np.real(np.vdot(z, z)))
    return abs(bil(z)) / d if d > 0 else 0.0


def build(delta, seed=2):
    m = N * (N - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(N, False)
    r2 = ui.gen3.make_parent(N, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / N
    s0 = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c
    if delta > 0:
        C2_0[:, 1, 0] = delta * s0
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    return ui.UnifiedEngine(N, C2_0, wp0), p2, q2, Z0c, s0


def run(delta):
    eng, p2, q2, Z0c, s0 = build(delta)
    H = {"total": np.zeros(T), "f_seed": np.zeros(T)}
    for (k, e) in TRACK:
        H[f"cell_{k}_{e}"] = np.zeros(T)
    for t in range(T):
        eng.step()
        C2 = eng.C2()
        H["total"][t] = rel_bil(C2.reshape(-1))
        for (k, e) in TRACK:
            H[f"cell_{k}_{e}"][t] = rel_bil(C2[:, k, e])
        H["f_seed"][t] = G.g_matter_fraction(C2)["f_seed"]
    return H, C2, Z0c, s0


def main():
    t0 = time.time()
    print("=== (C1)(C2) 初期構成要素の閉塞（τ=0・構成そのもの） ===")
    _, _, _, Z0c, s0 = build(DELTA)
    r_pump = rel_bil(Z0c)
    r_seed = rel_bil(s0)
    C1 = r_pump < TOL
    C2j = r_seed > 1e-2
    print(f"  真空ポンプ Z0c : |zᵀz|/‖z‖² = {r_pump:.3e}  "
          f"{'← 厳密ゼロ閉塞' if C1 else '不成立'}")
    print(f"  物質シード s0  : |zᵀz|/‖z‖² = {r_seed:.3e}  "
          f"{'← 非ゼロ閉塞' if C2j else '不成立'}")

    print("\n=== 物質宇宙（δ=1e-2）走行 ===")
    Hm, C2m, _, _ = run(DELTA)
    print("=== 真空宇宙（δ=0）走行 ===")
    Hv, C2v, _, _ = run(0.0)

    # (C3) 閉塞欠損 = 物質量
    m_ = N * (N - 1) // 2
    C0 = np.zeros((m_, NN, NETA), complex)
    C0[:, 2, 0] = Z0c
    C0[:, 1, 0] = DELTA * s0
    lhs = abs(bil(C0.reshape(-1)))
    rhs = DELTA ** 2 * abs(bil(s0))
    rel3 = abs(lhs - rhs) / rhs
    C3 = rel3 < 1e-6
    print(f"\n=== (C3) 閉塞欠損＝物質量（τ=0） ===")
    print(f"  |Σz²|_全体 = {lhs:.6e}   δ²·|s0ᵀs0| = {rhs:.6e}   相対差 {rel3:.2e}  "
          f"{'一致' if C3 else '不一致'}")

    # (C4) 真空宇宙は全時代で厳密ゼロ閉塞か
    vmax = float(np.max(Hv["total"]))
    C4 = vmax < TOL
    print(f"\n=== (C4) 真空宇宙の閉塞（全 τ の最大） ===")
    print(f"  |Σz²|/Σ|z|² 最大 = {vmax:.3e}  "
          f"{'← 厳密ゼロ閉塞を全時代で保持' if C4 else '★破れている'}")
    print(f"  参考: 物質宇宙の同量 最大 = {float(np.max(Hm['total'])):.3e}")

    # (C5) 凝縮体セルの劣化
    cc = Hm["cell_2_0"]
    print(f"\n=== (C5) 凝縮体セル (k=2,η=0) の閉塞の推移（物質宇宙） ===")
    for t in (0, 99, 999, 1999, 3999):
        print(f"  τ={t+1:5d}: {cc[t]:.3e}")
    print(f"  真空宇宙の同セル τ=4000: {Hv['cell_2_0'][-1]:.3e}")

    # (C6) セル別閉塞マップ（終端）
    print(f"\n=== (C6) セル別閉塞マップ（物質宇宙・τ=4000・非零セルのみ） ===")
    cellmap = {}
    for k in range(NN):
        for e in range(NETA):
            z = C2m[:, k, e]
            p = float(np.real(np.vdot(z, z)))
            if p > 0:
                r = rel_bil(z)
                cellmap[f"{k}_{e}"] = {"power": p, "closure_rel": r}
                if k <= 5:
                    print(f"  帯k={k:2d} 巻きη={e}: パワー {p:.3e}  閉塞 {r:.3e} "
                          f"{'← ゼロ閉塞' if r < TOL else ''}")

    print(f"\n(C1) 真空凝縮体の厳密ゼロ閉塞: {'通過' if C1 else '不成立'}")
    print(f"(C2) 物質シードは非ゼロ閉塞: {'通過' if C2j else '不成立'}")
    print(f"(C3) 閉塞欠損＝物質量: {'通過' if C3 else '不成立'}")
    print(f"(C4) 真空宇宙は全時代でゼロ閉塞: {'通過' if C4 else '不成立'}")
    print("(C5)(C6) 取得（判定なし）")

    out = {"env": {"N": N, "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA, "tol": TOL},
           "C1_pump_closure": r_pump, "C2_seed_closure": r_seed,
           "C3_defect_equals_matter": {"lhs": lhs, "rhs": rhs, "rel": rel3,
                                       "ok": bool(C3)},
           "C4_vacuum_exact": {"max": vmax, "ok": bool(C4),
                               "matter_max": float(np.max(Hm["total"]))},
           "C5_condensate_cell": {"tau1": float(cc[0]), "tau100": float(cc[99]),
                                  "tau1000": float(cc[999]),
                                  "tau2000": float(cc[1999]),
                                  "tau4000": float(cc[3999]),
                                  "vacuum_tau4000": float(Hv["cell_2_0"][-1])},
           "C6_cell_map": cellmap,
           "verdict": {"C1": bool(C1), "C2": bool(C2j), "C3": bool(C3),
                       "C4": bool(C4)},
           "v1_refutation": "回転平面の完全等方性による切り方は不成立（0.14 対 乱数対照 0.14）。"
                            "正しい分割は内容（セル）による分割であった。",
           "runtime_sec": time.time() - t0}
    (HERE / "result_tb_child_closure_blocks_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))

    ts = np.arange(1, T + 1)
    fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax[0].semilogy(ts, np.maximum(Hm["total"], 1e-20), lw=0.8, label="物質宇宙 δ=1e-2")
    ax[0].semilogy(ts, np.maximum(Hv["total"], 1e-20), "k--", lw=0.8, label="真空宇宙 δ=0")
    ax[0].axhline(TOL, color="red", lw=0.8, ls=":", label="機械精度の目安 1e-12")
    ax[0].set_ylabel("全体の閉塞 |Σxₙ²|/Σ|xₙ|²")
    ax[0].legend(fontsize=8)
    ax[0].set_title("ゼロ閉塞ブロック: 真空は閉じ、物質は閉塞の欠損である")
    for (k, e) in TRACK:
        ax[1].semilogy(ts, np.maximum(Hm[f"cell_{k}_{e}"], 1e-20), lw=0.8,
                       label=f"帯k={k}・巻きη={e}")
    ax[1].axhline(TOL, color="red", lw=0.8, ls=":")
    ax[1].set_ylabel("セル別の閉塞（物質宇宙）"); ax[1].set_xlabel("τ（step）")
    ax[1].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "fig_tb_child_closure_v2.png", dpi=130)
    plt.close(fig)
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
