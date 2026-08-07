#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ゼロ点対照ラン（標準環境 Nn=16・シード集団）— 基準値と基準粒子構成表の取得

問題（監査 2026-08-08）: 旧 §8 の「基準値」は単一シード・単一終端の一点値で
ゆらぎ帯がなく、しかも時代・census は Nn=5、位置は Nn=16 と条件が割れていた。
実験（ABCD投入・反跳・EPR）が走るのは Nn=16 側なのに、そこには基準がない。
本ランは条件を Nn=16 に一本化し、シード集団で平均とゆらぎ帯を取得する。

構成: N=12・Nn=16・Nη=8・中性レシピ・T=4000・統一 G v2 を第0步から毎ステップ。
  物質宇宙: δ=1e-2、シード集団 = gen3.make_parent(N, seed=s), s∈{2,3,4,5,6}
  真空宇宙: δ=0（構成は完全決定論——真空側にシード自由度は無いので1本）
読出しは統一 G v2（曖昧さ保存）。確定値が要る箇所は選択層 selection_v1 を
宣言して適用する（S∘G）。

事前登録した判定（実行前固定）:
 (Z1) 全時代定義性: 全ステップ・全シードでパネルに NaN/Inf が無い
      （位置の x は不在時 NaN を正規の表現とするため content_power で判定）。
 (Z2) 受動性: パネル併走の有無で終状態ビット同一。
 (Z3) ゼロ点帯の取得義務: 安定窓の平均・std・5/50/95分位を全メンバーで記録。
 (Z4) 定常性の分類: |窓内ドリフト/平均| < 0.05 を満たす量と満たさない量を
      分類して記録する（合否判定ではなく事実の分類）。
 (Z5) 真空対照: 物質パワー厳密 0・選択層 s_clock_acquirable が False。
 (Z6) 個数のゼロ点: 較正値 ε（result_tb_epsilon_calibration_v1.json）により
      n̂_P = P_matter/ε、n̂_PR = pr_n/pr_n(単一) を記録する。
      非局在内容では pr_n = Nn（最大）＝「可算な局在ドメインなし」。

使い方: python3 run_tb_zeropoint_Nn16_v1.py
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


ui = load("ui_zp", UF / "unified_interaction_v1.py")
G = load("ur_zp", UF / "unified_readout_v2.py")
S = load("sel_zp", UF / "selection_v1.py")

N, NN, NETA = 12, 16, 8
T = 4000
DELTA = 1e-2
SEEDS = [2, 3, 4, 5, 6]
WIN = (2000, 4000)
KEYS = ("f2", "f_seed", "r", "pr_n", "content_power", "mix_mean",
        "coherence", "carrier_power", "P_tot", "phase")


def build(delta, seed):
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
    return ui.UnifiedEngine(N, C2_0, wp0), p2, q2


def run(delta, seed, panel_on=True):
    eng, p2, q2 = build(delta, seed)
    H = {k: np.zeros(T) for k in KEYS}
    carry = {"C_flat": None, "c_gen": None}
    nonfinite = 0
    acq = np.zeros(T, bool)
    for t in range(T):
        eng.step()
        if not panel_on:
            continue
        pan = G.g_panel(eng.C2(), p2, q2, carry["C_flat"], carry["c_gen"])
        carry = pan["_carry"]
        for k in KEYS:
            H[k][t] = pan.get(k, np.nan) if k != "r" else pan.get("r", 0.0)
        acq[t] = S.s_clock_acquirable(pan)["acquirable"]
        vals = [pan["f2"], pan["f_seed"], pan["P_tot"], pan["pr_n"],
                pan["content_power"], pan["mix_mean"]]
        if not all(np.isfinite(v) for v in vals):
            nonfinite += 1
    cell = G.g_cell_ledger(eng.C2())
    return eng, H, acq, nonfinite, cell


def band(x):
    w = x[WIN[0]:WIN[1]]
    w = w[np.isfinite(w)]
    if len(w) == 0:
        return None
    drift = float(np.polyfit(np.arange(len(w)), w, 1)[0]) * len(w)
    mean = float(np.mean(w))
    return {"mean": mean, "std": float(np.std(w)),
            "p05": float(np.percentile(w, 5)), "p50": float(np.percentile(w, 50)),
            "p95": float(np.percentile(w, 95)), "drift": drift,
            "drift_rel": abs(drift / mean) if mean != 0 else float("inf")}


def main():
    t0 = time.time()
    cal = json.loads((HERE / "result_tb_epsilon_calibration_v1.json").read_text())
    eps, pr1 = cal["epsilon"], cal["pr_n_single"]
    print(f"=== ゼロ点対照ラン（N={N}・Nn={NN}・T={T}・シード{SEEDS}） ===")
    print(f"較正値: ε={eps:.6e}  pr_n(単一ドメイン)={pr1:.4f}")

    runs, cells, nf_tot = {}, {}, 0
    for s in SEEDS:
        _, H, acq, nf, cell = run(DELTA, s)
        runs[s] = (H, acq)
        cells[s] = cell
        nf_tot += nf
        b = band(H["f_seed"])
        print(f"  seed={s}: f_seed={b['mean']:.6e}±{b['std']:.1e} "
              f"pr_n={band(H['pr_n'])['mean']:.3f} "
              f"時計取得率={acq[WIN[0]:WIN[1]].mean():.3f} [{time.time()-t0:.0f}s]")
    _, Hv, acqv, nfv, cellv = run(0.0, SEEDS[0])
    nf_tot += nfv

    # --- Z3: メンバーごとのゼロ点帯（集団統合）
    zero_point = {}
    for k in KEYS:
        per_seed = [band(runs[s][0][k]) for s in SEEDS]
        per_seed = [b for b in per_seed if b]
        if not per_seed:
            continue
        means = np.array([b["mean"] for b in per_seed])
        zero_point[k] = {
            "ensemble_mean": float(means.mean()),
            "ensemble_std": float(means.std()),
            "within_std_mean": float(np.mean([b["std"] for b in per_seed])),
            "p05_min": float(min(b["p05"] for b in per_seed)),
            "p95_max": float(max(b["p95"] for b in per_seed)),
            "drift_rel_max": float(max(b["drift_rel"] for b in per_seed)),
        }
    stationary = {k: v["drift_rel_max"] < 0.05 for k, v in zero_point.items()}

    # --- Z5: 真空対照
    oddm = G._odd_mask(NN)
    P_mat_vac = float(cellv["cell_power"][oddm, :].sum())
    acq_vac = float(acqv[WIN[0]:WIN[1]].mean())
    Z5 = bool(P_mat_vac == 0.0 and acq_vac == 0.0)

    # --- Z2: 受動性
    ea, _, _, _, _ = run(DELTA, SEEDS[0], panel_on=True)
    eb, _, _, _, _ = run(DELTA, SEEDS[0], panel_on=False)
    Z2 = bool(np.array_equal(ea.C, eb.C))
    Z1 = bool(nf_tot == 0)

    # --- Z6: 基準粒子構成表（セル別・集団平均±ゆらぎ・個数）
    Pstack = np.stack([cells[s]["cell_power"] for s in SEEDS])   # S×Nn×Nη
    PRstack = np.stack([cells[s]["cell_pr_m"] for s in SEEDS])
    Pm, Ps = Pstack.mean(axis=0), Pstack.std(axis=0)
    table = []
    for k in range(NN):
        for e in range(NETA):
            if Pm[k, e] == 0.0 and Ps[k, e] == 0.0:
                continue
            table.append({"k": int(k), "eta": int(e),
                          "parity": "奇" if k % 2 == 1 else ("零" if k == 0 else "偶"),
                          "power_mean": float(Pm[k, e]), "power_std": float(Ps[k, e]),
                          "pr_m_mean": float(PRstack.mean(axis=0)[k, e]),
                          "n_hat_P": float(Pm[k, e] / eps),
                          "exact_zero": False})
    n_zero_cells = int((Pstack == 0.0).all(axis=0).sum())
    prn_mean = zero_point["pr_n"]["ensemble_mean"]
    n_hat_PR = prn_mean / pr1
    # 局在の判定は正本のドメイン基準 pr_n < 0.5·Nn（位置検収 P1 と同一）。
    # 【自己是正】初版は |pr_n − Nn| > 1e-6 という過剰に厳しい判定で、
    # 一様（pr_n=15.99983）を「局在」と誤分類した。
    localized = bool(prn_mean < 0.5 * NN)
    # 【自己是正】個数は物質パワー（帳簿の奇数帯パワー）で数える。初版は
    # content_power（双対プロファイルの総和＝Nn·Nη 倍の規格化）を使っており
    # 128 倍ずれていた。ε 較正と同じ量（cell_power[odd] の和）に統一する。
    P_matter_mean = float(Pstack[:, oddm, :].sum(axis=(1, 2)).mean())
    P_matter_std = float(Pstack[:, oddm, :].sum(axis=(1, 2)).std())
    n_hat_P = P_matter_mean / eps

    print(f"\n(Z1) 全時代定義性（非有限 {nf_tot} 件）: {'通過' if Z1 else '不成立'}")
    print(f"(Z2) 受動性（パネル有無でビット同一）: {'通過' if Z2 else '不成立'}")
    print(f"(Z3) ゼロ点帯: {len(zero_point)} メンバーで平均・std・分位を取得")
    print(f"(Z4) 定常性の分類（|ドリフト/平均|<0.05）:")
    for k, v in zero_point.items():
        print(f"      {k:15s} 帯平均={v['ensemble_mean']:.6e} "
              f"シード間std={v['ensemble_std']:.2e} "
              f"窓内std={v['within_std_mean']:.2e} "
              f"ドリフト比={v['drift_rel_max']:.3f} "
              f"{'定常' if stationary[k] else '★非定常'}")
    print(f"(Z5) 真空対照（物質パワー={P_mat_vac:.1e}・時計取得率={acq_vac:.2f}）: "
          f"{'通過' if Z5 else '不成立'}")
    print(f"(Z6) 個数のゼロ点: 物質パワー {P_matter_mean:.6e}±{P_matter_std:.1e} "
          f"→ n̂_P = {n_hat_P:.4f} 個相当（ε={eps:.3e}）／ "
          f"pr_n={prn_mean:.4f}（Nn={NN}＝一様の署名・基準 0.5·Nn={0.5*NN}）→ "
          f"局在ドメイン数 = "
          f"{'{:.3f}'.format(n_hat_PR) if localized else '0（非局在＝可算ドメインなし）'}")
    print(f"      厳密0セル: {n_zero_cells}/{NN*NETA}   非零セル: {len(table)}")

    out = {"env": {"N": N, "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA,
                   "seeds": SEEDS, "window": WIN,
                   "readout": "unified_readout_v2 (曖昧さ保存)",
                   "selection": "selection_v1: s_clock_acquirable(既定床)"},
           "calibration": {"epsilon": eps, "pr_n_single": pr1,
                           "source": "result_tb_epsilon_calibration_v1.json"},
           "Z1_defined": Z1, "Z2_passive": Z2, "zero_point": zero_point,
           "Z4_stationary": stationary,
           "Z5_vacuum": {"P_matter": P_mat_vac, "clock_rate": acq_vac, "ok": Z5},
           "Z6_counting": {"P_matter_mean": P_matter_mean,
                           "P_matter_std": P_matter_std,
                           "n_hat_P": n_hat_P, "pr_n": prn_mean,
                           "n_hat_PR_if_localized": n_hat_PR,
                           "localized": localized,
                           "domain_count": (n_hat_PR if localized else 0.0)},
           "baseline_table": table, "exact_zero_cells": n_zero_cells,
           "nonfinite": nf_tot, "runtime_sec": time.time() - t0}
    (HERE / "result_tb_zeropoint_Nn16_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))

    # --- 図
    ts = np.arange(1, T + 1)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for s in SEEDS:
        H = runs[s][0]
        axes[0].semilogy(ts, np.maximum(H["f2"], 1e-16), lw=0.7, alpha=0.8,
                         label=f"seed {s}")
        axes[1].semilogy(ts, np.maximum(H["f_seed"], 1e-16), lw=0.7, alpha=0.8)
        axes[2].plot(ts, H["pr_n"], lw=0.7, alpha=0.8)
    axes[0].semilogy(ts, np.maximum(Hv["f2"], 1e-16), "k--", lw=0.9, label="真空 δ=0")
    axes[1].semilogy(ts, np.maximum(Hv["f_seed"], 1e-16), "k--", lw=0.9)
    axes[0].set_ylabel("f₂（空間形成史）"); axes[0].legend(fontsize=8, ncol=3)
    axes[0].set_title(f"ゼロ点対照ラン N={N}・Nn={NN}・シード集団（統一G v2・常時実行）")
    axes[1].set_ylabel("f_seed（物質分率）")
    axes[2].axhline(NN, color="red", lw=0.8, label=f"Nn={NN}（一様＝ドメインなし）")
    axes[2].set_ylabel("pr_n（双対占有）"); axes[2].set_xlabel("t（step）")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_tb_zeropoint_Nn{NN}_v1.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(np.log10(np.maximum(Pm, 1e-300)), aspect="auto",
                   origin="lower", cmap="viridis")
    ax.set_xlabel("巻き η（電荷）"); ax.set_ylabel("帯 k")
    ax.set_title(f"基準粒子構成（集団平均 log₁₀P・N={N}・Nn={NN}・T={T}）")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_tb_zeropoint_census_Nn{NN}_v1.png", dpi=130)
    plt.close(fig)
    np.savez_compressed(HERE / f"tb_zeropoint_Nn{NN}_v1.npz",
                        **{f"s{s}_{k}": runs[s][0][k] for s in SEEDS for k in KEYS},
                        **{f"vac_{k}": Hv[k] for k in KEYS},
                        cell_power_mean=Pm, cell_power_std=Ps)
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
