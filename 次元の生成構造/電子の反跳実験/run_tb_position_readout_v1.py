#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TB位置読出し実演: 空間が読めることの検収（仕様書§7 xyzt節のデータ取得）

統一関数のみ使用（木原ルール）。標準物質宇宙（N=12・中性レシピ）で、
ニュートリノ型物質（周期表[4]同定: 巻き0×奇数帯=中性フェルミオン）を
レジスタ局在パケットとして生成し、統一読出し g_position_1d が
位置を読めること・生誕位置に留まること（粒子は位置を持つ・[2]）を検収。

構成:
 (a) 局在シード: 奇数帯 k=1,3 の同振幅・双対点 n0=2 同位相パケット（中性）
 (b) 非局在シード: 単一帯 k=1（標準・双対空間で一様）
 (c) 真空: シードなし
判定（v2・正本基準に整合）:
 (P1) 局在（正本 genesis v3 の V3-P2 と同形のドメイン基準）:
      present=True 全窓・安定窓の PR_n < 0.5·Nn（局在ドメインの維持）・
      circular 平均位置が生誕点 n0=2 と ±1 セル以内（生誕域残留）。
      （初版の「±0.5セル重心静止」は正本に無い過剰基準だった——
      2帯パケットの呼吸による重心振動はドメイン維持と両立する。修正記録。）
 (P2) 非局在: present=True・PR_n が高い（>0.6·Nn＝一様の署名）を記録。
 (P3) 真空: 物質不在で present=False（空間は「空として」正しく読める）。
図: fig_tb_position_Nn%d_v1.png" % NN if False else "fig_tb_position_v1.png（x(t) 軌跡・局在/非局在/真空）。
使い方: python3 run_tb_position_readout_v1.py
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


ui = load("ui_pos", UF / "unified_interaction_v1.py")
ur = load("ur_pos", UF / "unified_readout_v2.py")     # 曖昧さ保存版（現行）
S = load("sel_pos", UF / "selection_v1.py")           # 選択層（宣言して使う）

N = 12
NN = int(sys.argv[1]) if len(sys.argv) > 1 else 5   # 空間分解能（正本v3=16）
NETA = 8
T = 4000
DELTA = 1e-2
N0 = 2 if NN == 5 else NN // 4   # 生誕点（双対レジスタ）


def build(seed_mode):
    m = N * (N - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(N, False)
    r2 = ui.gen3.make_parent(N, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / N
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c
    if seed_mode == "localized":
        # 全奇数帯の同振幅・双対点 n0 同位相パケット（中性・毛0）——
        # 正本 genesis v3 のレジスタ局在構成。総パワーは非局在と同一。
        odd_ks = list(range(1, NN, 2))
        for k in odd_ks:
            C2_0[:, k, 0] = (DELTA / np.sqrt(len(odd_ks))) * seed_state * \
                np.exp(-2j * np.pi * k * N0 / NN)
    elif seed_mode == "delocalized":
        C2_0[:, 1, 0] = DELTA * seed_state
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    nq = np.linalg.norm(q2)
    q2 = q2 / nq if nq > 1e-12 else np.zeros_like(p2)
    return ui.UnifiedEngine(N, C2_0, wp0), p2, q2


def run(seed_mode):
    eng, p2, q2 = build(seed_mode)
    xs, prs, present = [], [], []
    carry_f, carry_g = None, None
    for t in range(T):
        eng.step()
        pan = ur.g_panel(eng.C2(), p2, q2, carry_f, carry_g)
        carry_f = pan["_carry"]["C_flat"]; carry_g = pan["_carry"]["c_gen"]
        # 【宣言】選択 S = s_position_maxmoment（巻き m=1,2 のうち重み最大を採る）
        # ＋ s_present（既定床 1e-280）。G は束を返すだけで選択しない（R7）。
        pos = S.s_position_maxmoment(pan)
        present.append(pos["present"])
        xs.append(pos["x"] if pos["present"] else np.nan)
        # PR_n（双対占有）は統一 G の読出し値（是正記録: 初版は「補助統計」と
        # 称してインラインFFTで内部状態を直接読んでいた＝ルール違反——
        # 統一 G の pr_n 出力に置換）
        prs.append(pan["pr_n"] if pos["present"] else np.nan)
    return np.array(xs), np.array(prs), np.array(present)


def circ_stats(x, Nn):
    ang = x * 2 * np.pi / Nn
    z = np.nanmean(np.exp(1j * ang))
    mean = (np.angle(z) * Nn / (2 * np.pi)) % Nn
    R = abs(z)
    std = np.sqrt(max(-2 * np.log(max(R, 1e-300)), 0.0)) * Nn / (2 * np.pi)
    return float(mean), float(std)


def main():
    t0 = time.time()
    xs_l, pr_l, ps_l = run("localized")
    xs_d, pr_d, ps_d = run("delocalized")
    xs_v, pr_v, ps_v = run("vacuum")
    W = slice(2000, 4000)
    half = NN / 2.0   # 奇数帯内容は被覆2 → 位置は mod Nn/2（統一読出しの cover 出力）
    mean_l, std_l = circ_stats(xs_l[W] % half, half)
    dev = abs(((mean_l - (N0 % half)) + half / 2) % half - half / 2)
    pr_l_med = float(np.nanmedian(pr_l[W]))
    p1 = bool(ps_l[W].all() and pr_l_med < 0.5 * NN and dev < 1.0)
    pr_d_med = float(np.nanmedian(pr_d[W]))
    p2j = bool(ps_d[W].all() and pr_d_med > 0.6 * NN)
    p3 = bool(~ps_v.any())
    print(f"(P1) 局在ニュートリノ型: 位置={mean_l:.3f}（生誕点{N0}・偏差{dev:.3f}<1.0）"
          f" PR_n中央値={pr_l_med:.2f}（<{0.5*NN:.1f}＝ドメイン維持） "
          f"重心std={std_l:.3f}セル（記録） → {'通過' if p1 else '不成立'}")
    print(f"(P2) 非局在: present全True・PR_n中央値={pr_d_med:.2f}"
          f"（>{0.6*NN:.1f}） → {'通過' if p2j else '不成立'}")
    print(f"(P3) 真空: present={ps_v.any()}（False=空として読める） → "
          f"{'通過' if p3 else '不成立'}")
    ok = p1 and p2j and p3
    print("位置読出し検収: " + ("PASS" if ok else "不成立あり"))

    ts = np.arange(1, T + 1)
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax = axes[0]
    ax.plot(ts, xs_l, ".", ms=1, label="局在ニュートリノ型（生誕点 n₀=2）")
    ax.plot(ts, xs_d, ".", ms=1, alpha=0.4, label="非局在（一様）")
    ax.axhline(N0, color="red", lw=0.8, label="生誕点")
    ax.set_ylabel("位置読出し x（双対レジスタ）")
    ax.legend(); ax.set_title(f"N={N} 位置読出しの検収（真空は present=False・非表示）")
    ax = axes[1]
    ax.plot(ts, pr_l, ".", ms=1, label="局在 PR_n")
    ax.plot(ts, pr_d, ".", ms=1, alpha=0.4, label="非局在 PR_n")
    ax.set_ylabel("双対占有 PR_n"); ax.set_xlabel("t（step）"); ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / f"fig_tb_position_Nn{NN}_v1.png", dpi=130)
    plt.close(fig)

    out = {"N": N, "n0": N0, "T": T,
           "P1": {"mean": mean_l, "std": std_l, "dev": dev,
                  "PR_median": pr_l_med, "ok": p1},
           "P2": {"PR_median": pr_d_med, "ok": p2j},
           "P3": {"present_any": bool(ps_v.any()), "ok": p3},
           "all_pass": bool(ok), "runtime_sec": time.time() - t0}
    (HERE / f"result_tb_position_readout_Nn{NN}_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
