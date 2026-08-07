#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ε較正実験 — 「個数」を測定量として定義する（テストベッド仕様書 §8.1a の再構築）

問題（監査 2026-08-08）: 旧仕様の「波の本数」は count = #(振幅²>0) の
台集合濃度であり、パワー 1.0 の凝縮体も 9.1e-22 の数値塵も同じ M 本になる
——個数の量として機能していなかった。本実験は個数を較正で定義し直す。

定義（本実験で検証する）:
  ε      := 単一局在ドメインあたりの物質パワー（較正値）
  n̂_P    := P_matter / ε                      （パワー基準の個数）
  n̂_PR   := pr_n / pr_n(単一ドメイン)          （局在ドメイン基準の個数）
  pr_n は統一 G v2 の位置スペクトルが返す双対占有の実効セル数（連続量）。
  奇数帯（フェルミオン型）内容は対蹠二点構造（周期表柱7）を持つため
  単一ドメインの pr_n は 2 になる見込み——これを較正で確定する。

構成: 標準環境（N=12・Nn=16・Nη=8・中性レシピ）に、同一の局在パケットを
j=1,2,3 個、双対レジスタ上の互いに素な位置に置く。位置は差が Nn/2=8 の
倍数にならないよう選ぶ（対蹠像の重なりを避ける）。1個あたりの振幅は同一。

事前登録した判定（実行前固定）:
 (E1) 加法性・パワー: P_matter(j) = j·ε を相対誤差 <1e-12 で満たす（t=0）。
 (E2) 加法性・ドメイン: pr_n(j) = j·pr_n(1) を相対誤差 <1e-12 で満たす（t=0）。
 (E3) 較正の耐久性: 力学発展後（T=2000 の安定窓中央値）も
      n̂_PR が整数 j から ±0.15 以内に留まる。
 (E4) 記録: 位置スペクトルの巻きモーメント重み分布（選択前の束）を保存する。
判定は「較正が成立するか」であり、成立しない場合もそのまま記録する。

使い方: python3 run_tb_epsilon_calibration_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UF = HERE.parent / "統一万能関数_v1"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ui = load("ui_eps", UF / "unified_interaction_v1.py")
G = load("ur_eps", UF / "unified_readout_v2.py")
S = load("sel_eps", UF / "selection_v1.py")

N, NN, NETA = 12, 16, 8
T = 2000
DELTA = 1e-2
SITES = [4, 9, 14]          # 差 5, 5, 10 — いずれも Nn/2=8 の倍数でない
WIN = slice(1000, 2000)


def build(n_packets, seed=2):
    m = N * (N - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(N, False)
    r2 = ui.gen3.make_parent(N, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / N
    s0 = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c                       # 真空（光子型凝縮体・偶数帯）
    odd = list(range(1, NN, 2))
    amp = DELTA / np.sqrt(len(odd))           # 1パケットあたりの振幅（同一）
    for i in range(n_packets):
        n0 = SITES[i]
        for k in odd:
            C2_0[:, k, 0] += amp * s0 * np.exp(-2j * np.pi * k * n0 / NN)
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    return ui.UnifiedEngine(N, C2_0, wp0), p2, q2


def measure(C2, p2, q2, carry=None):
    pan = G.g_panel(C2, p2, q2,
                    carry["C_flat"] if carry else None,
                    carry["c_gen"] if carry else None)
    odd = G._odd_mask(NN)
    P_matter = float(pan["cell_power"][odd, :].sum())
    return pan, P_matter


def main():
    t0 = time.time()
    out = {"env": {"N": N, "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA,
                   "sites": SITES}, "j": {}}
    print(f"=== ε較正（N={N}・Nn={NN}・δ={DELTA}・位置{SITES}） ===")
    for j in (1, 2, 3):
        eng, p2, q2 = build(j)
        pan0, P0 = measure(eng.C2(), p2, q2)
        rec = {"P_matter_t0": P0, "pr_n_t0": pan0["pr_n"],
               "mix_mean_t0": pan0["mix_mean"],
               "pos_weight_t0": list(map(float, pan0["pos_weight"])),
               "pos_x_t0": list(map(float, pan0["pos_x"]))}
        carry = {"C_flat": None, "c_gen": None}
        prs, Pms = [], []
        for t in range(T):
            eng.step()
            pan, Pm = measure(eng.C2(), p2, q2, carry)
            carry = pan["_carry"]
            prs.append(pan["pr_n"]); Pms.append(Pm)
        prs = np.array(prs); Pms = np.array(Pms)
        rec.update({"pr_n_win_median": float(np.median(prs[WIN])),
                    "pr_n_win_std": float(np.std(prs[WIN])),
                    "P_matter_win_median": float(np.median(Pms[WIN])),
                    "pos_weight_end": list(map(float, pan["pos_weight"])),
                    "pos_x_end": list(map(float, pan["pos_x"])),
                    "mix_mean_end": float(pan["mix_mean"])})
        out["j"][j] = rec
        print(f"  j={j}: t=0 P_matter={P0:.6e} pr_n={pan0['pr_n']:.6f}  "
              f"→ 安定窓 P={rec['P_matter_win_median']:.6e} "
              f"pr_n={rec['pr_n_win_median']:.4f}")

    eps = out["j"][1]["P_matter_t0"]
    pr1 = out["j"][1]["pr_n_t0"]
    e1 = max(abs(out["j"][j]["P_matter_t0"] / (j * eps) - 1) for j in (1, 2, 3))
    e2 = max(abs(out["j"][j]["pr_n_t0"] / (j * pr1) - 1) for j in (1, 2, 3))
    n_hat = {j: out["j"][j]["pr_n_win_median"] / pr1 for j in (1, 2, 3)}
    e3 = max(abs(n_hat[j] - j) for j in (1, 2, 3))
    E1, E2, E3 = e1 < 1e-12, e2 < 1e-12, e3 < 0.15
    print(f"\n(E1) パワー加法性 P(j)=j·ε: 最大相対誤差 {e1:.2e} <1e-12 → "
          f"{'通過' if E1 else '不成立'}")
    print(f"(E2) ドメイン加法性 pr_n(j)=j·pr_n(1): 最大相対誤差 {e2:.2e} <1e-12 → "
          f"{'通過' if E2 else '不成立'}")
    print(f"(E3) 較正の耐久性 n̂_PR={{{', '.join(f'{j}:{n_hat[j]:.3f}' for j in n_hat)}}} "
          f"最大偏差 {e3:.3f} <0.15 → {'通過' if E3 else '不成立'}")
    print(f"\n較正値: ε = {eps:.9e}（単一局在ドメインの物質パワー）")
    print(f"        pr_n(単一ドメイン) = {pr1:.6f}"
          f"（対蹠二点構造なら 2）")
    print(f"個数の定義: n̂_P = P_matter/ε ・ n̂_PR = pr_n/{pr1:.4f}")
    out.update({"epsilon": eps, "pr_n_single": pr1,
                "E1_power_additivity": {"max_rel_err": e1, "ok": bool(E1)},
                "E2_domain_additivity": {"max_rel_err": e2, "ok": bool(E2)},
                "E3_durability": {"n_hat": n_hat, "max_dev": e3, "ok": bool(E3)},
                "all_pass": bool(E1 and E2 and E3),
                "runtime_sec": time.time() - t0})
    (HERE / "result_tb_epsilon_calibration_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
