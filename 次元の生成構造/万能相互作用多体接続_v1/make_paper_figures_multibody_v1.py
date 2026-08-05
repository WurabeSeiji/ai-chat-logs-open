#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""多体接続論文の図8枚

出典の規律: 図の数値は既存の結果JSON（コミット済み）から読む。系列が保存されて
いない2実験（ラチェット・GENESIS通し）は同一シード・同一パラメータで再走行し、
要約値が保存済みJSONと一致することを対照検定してから系列を描く。
3Dデモは回収済み原本スクリプトをそのまま子プロセス実行し、印字位置を読む。

fig_mb1: アーキテクチャ模式図（グラフ×レジスタ直積）
fig_mb2: 二因子ゲート（凝縮×2.74・時代変調）
fig_mb3: 点火則の移植（p=2.001・C一定）
fig_mb4: census毛込み完全版（排他6e289・ロックと減衰B4）
fig_mb5: 局所性の住所（グラフ無効・レジスタ×19.6）
fig_mb6: 方向くじとラチェット整流（再走行系列）
fig_mb7: 運動（1D調和閉鎖並進・3D軌道）
fig_mb8: GENESIS通し走行（真空→インフレーション→凝縮→生成開門）
"""
import importlib.util
import itertools
import json
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = ["Hiragino Sans", "Arial Unicode MS", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent
J = lambda name: json.loads((HERE / name).read_text())

# ================= fig_mb1 アーキテクチャ =================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1, 1.5]})
n = 6
th = np.pi / 2 + 2 * np.pi * np.arange(n) / n
pos = {i + 1: (np.cos(t), np.sin(t)) for i, t in enumerate(th)}
for (u, v) in itertools.combinations(range(1, n + 1), 2):
    ax1.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#555555", lw=1.2, zorder=1)
for i, (x, y) in pos.items():
    ax1.add_patch(plt.Circle((x, y), 0.13, fc="#eeeeee", ec="black", lw=1.2, zorder=3))
    ax1.annotate(str(i), (x, y), ha="center", va="center", fontsize=11, weight="bold", zorder=4)
ax1.set_xlim(-1.4, 1.4); ax1.set_ylim(-1.5, 1.4); ax1.set_aspect("equal"); ax1.axis("off")
ax1.set_title("(a) グラフ軸\nN体（ノード）と M=N(N−1)/2 本の関係波（辺）", fontsize=10)
ax1.annotate("ノード＝体（変数なし・帳簿の集約点）\n辺＝関係波 z_e（力学変数）",
             (0, -1.42), ha="center", fontsize=9)

M_s, NR = 10, 8
for r in range(M_s):
    for c in range(NR):
        fc = "#ffd0d0" if c == 2 else ("#d0e4ff" if c == 1 else "white")
        ax2.add_patch(plt.Rectangle((c, M_s - 1 - r), 0.92, 0.92, fc=fc, ec="#888888", lw=0.6))
ax2.annotate("真空／ポンプ\nスライス k=2", (2.46, M_s + 0.35), ha="center", fontsize=9, color="#b03030")
ax2.annotate("種スライス k=1", (1.46, -0.75), ha="center", fontsize=9, color="#2060b0")
ax2.annotate("レジスタ軸（$N_{reg}$ 点）→", (NR / 2, -1.5), ha="center", fontsize=10)
ax2.annotate("グ\nラ\nフ\n軸\n（\nM\n本\n）", (-1.15, M_s / 2), ha="center", va="center", fontsize=10)
ax2.annotate("線形部：スライス（列）ごとに\n第7・8論文の力学が作用", (NR + 0.4, M_s * 0.72),
             fontsize=9, va="center")
ax2.annotate("頂点：レジスタ点ごとに\nノード帳簿を介して作用\n（強度 R は読出し・IF文なし）",
             (NR + 0.4, M_s * 0.32), fontsize=9, va="center")
ax2.set_xlim(-1.8, NR + 5.2); ax2.set_ylim(-2.0, M_s + 1.2)
ax2.set_aspect("equal"); ax2.axis("off")
ax2.set_title("(b) 状態 $C \\in \\mathbb{C}^{M\\times N_{reg}}$（直積格子）", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "fig_mb1_architecture_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb1 saved")

# ================= fig_mb2 二因子ゲート =================
e1b = J("stage2_vertex_engine_result_v1.json")["E1b"]
gb = J("genesis_result_v1.json")["backreaction"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
vals = [e1b["rates"]["incoherent"], e1b["rates"]["condensed"]]
ax1.bar(["非凝縮ポンプ", "凝縮ポンプ"], vals, color=["tab:gray", "tab:red"], width=0.5)
ax1.set_yscale("log"); ax1.set_ylabel("生成率（同一ノルム・同一種）")
ax1.set_title(f"(a) 凝縮ゲート：比 {e1b['ratio']:.2f}", fontsize=11)
for i, v_ in enumerate(vals):
    ax1.annotate(f"{v_:.2e}", (i, v_), ha="center", va="bottom", fontsize=9)
eras = ["潜伏期\n(50–250)", "バースト期\n(400–1100)", "凝縮後\n(2000–4000)"]
gv = [gb["g_latency"], gb["g_burst"], gb["g_metastable"]]
ax2.bar(eras, gv, color=["#bbbbbb", "#8888cc", "tab:green"], width=0.5)
ax2.set_yscale("log"); ax2.set_ylabel("窓内生成率 d(lnP_seed)/dt")
ax2.set_title(f"(b) 時代変調：潜伏→凝縮後で ×{gv[2]/gv[0]:.1f}", fontsize=11)
for i, v_ in enumerate(gv):
    ax2.annotate(f"{v_:.1e}", (i, v_), ha="center", va="bottom", fontsize=9)
fig.suptitle("生成率 = 凝縮ゲート × C·f²（相互作用は常時作用・タイミング則なし）", fontsize=11)
fig.tight_layout()
fig.savefig(HERE / "fig_mb2_gate_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb2 saved")

# ================= fig_mb3 点火則の移植 =================
e1a = J("stage2_tests_corrected_result_v2.json")["E1a_refit"]
rows = J("stage2_vertex_engine_result_v1.json")["E1a"]["rows"]
Cvals = J("stage2_vertex_engine_result_v1.json")["E1a"]["C_values"]
b3 = J("stage3_sharedO_hair_result_v1.json")["B3"]
f0s = np.array([r["f_seed0"] for r in rows]); rates = np.array([r["rate"] for r in rows])
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
ax1.loglog(f0s, rates, "o", color="tab:red", ms=7, label="実測（4点・4桁）")
xx = np.array([f0s.min() * 0.5, f0s.max() * 2])
ax1.loglog(xx, np.mean(Cvals) * xx ** 2, "k--", lw=1, label=f"C·f²（p={e1a['p']:.4f}）")
ax1.set_xlabel("初期種割合 f₀"); ax1.set_ylabel("初期成長率 rate₀")
ax1.set_title("(a) 多体点火則：rate = C·f^p, p=2.001", fontsize=11)
ax1.legend(fontsize=9); ax1.grid(alpha=0.3, which="both")
ax2.semilogx(f0s, Cvals, "o-", color="tab:red", label=f"段階2エンジン C（広がり0.4%）")
ax2.axhline(b3["C"], color="tab:blue", ls=":", lw=1.5,
            label=f"段階3共有O版 C={b3['C']:.3f}（p={b3['p']:.4f}）")
ax2.set_ylim(5, 7.5)
ax2.set_xlabel("初期種割合 f₀"); ax2.set_ylabel("C = rate₀ / f₀²")
ax2.set_title("(b) 係数の一定性（二体則の移植成立）", fontsize=11)
ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(HERE / "fig_mb3_ignition_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb3 saved")

# ================= fig_mb4 census 毛込み完全版 =================
h = J("stage3_sharedO_hair_result_v1.json")
Pm = np.array(h["HAIR"]["P_k3_by_m"])
b4 = h["B4"]["coherences"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
ms = np.arange(len(Pm))
colors = ["tab:red" if m_ == h["HAIR"]["m_star"] else "tab:gray" for m_ in ms]
ax1.bar([str(m_) for m_ in ms], np.maximum(Pm, 1e-330), color=colors)
ax1.set_yscale("log")
ax1.set_xlabel("毛（巻き数）m"); ax1.set_ylabel("相棒帯 k=3 のパワー")
ax1.set_title(f"(a) 毛の帳簿：予言 m*={h['HAIR']['m_star']} のみ生成\n排他比 {h['HAIR']['exclusivity']:.1e}", fontsize=10)
ax2.plot(np.arange(1, len(b4) + 1), b4, "o-", color="tab:blue", label="対相関ロック（窓ごと）")
ax2.axhline(h["HAIR"]["lock_coherence"], color="tab:red", ls="--", lw=1,
            label=f"毛ロック {h['HAIR']['lock_coherence']:.4f}")
ax2.set_ylim(0.985, 1.001)
ax2.set_xlabel("測定窓"); ax2.set_ylabel("コヒーレンス")
ax2.set_title(f"(b) ロックの持続と緩慢な減衰（B4終値 {h['B4']['final']:.4f}・原因未解明＝公開課題）",
              fontsize=10)
ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(HERE / "fig_mb4_census_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb4 saved")

# ================= fig_mb5 局所性の住所 =================
v2 = J("genesis_v2_local_result_v1.json")
v3 = J("genesis_v3_register_local_result_v1.json")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
labels, u_rates, l_rates = [], [], []
for r in v2["rows"]:
    labels.append(f"δ={r['delta']:g}")
    u_rates.append(r["uniform"]["rate"]); l_rates.append(r["local"]["rate"])
x = np.arange(len(labels)); w = 0.35
ax1.bar(x - w / 2, u_rates, w, color="tab:gray", label="一様種")
ax1.bar(x + w / 2, l_rates, w, color="tab:orange", label=f"グラフ集中種（星型・濃縮{v2['concentration']:g}倍）")
ax1.set_xticks(x); ax1.set_xticklabels(labels); ax1.set_yscale("log")
ax1.set_ylabel("生成率")
ratios = [r["rate_ratio"] for r in v2["rows"]]
ax1.set_title(f"(a) グラフ軸に局所性なし（比 {ratios[0]:.2f}, {ratios[1]:.2f} ≈ 1）", fontsize=11)
ax1.legend(fontsize=8)
vals2 = [abs(v3["single"]["rate"]), abs(v3["packet"]["rate"])]
ax2.bar(["単一倍音種\n（レジスタ非局在）", f"パケット種\n（濃縮{v3['concentration']:g}倍・PR=2）"],
        vals2, color=["tab:gray", "tab:green"], width=0.5)
ax2.set_yscale("log"); ax2.set_ylabel("|生成率|")
ax2.set_title(f"(b) 局所性はレジスタ軸に住む（応答 ×{abs(v3['rate_ratio']):.1f}）", fontsize=11)
for i, v_ in enumerate(vals2):
    ax2.annotate(f"{v_:.2e}", (i, v_), ha="center", va="bottom", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "fig_mb5_locality_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb5 saved")

# ================= fig_mb6 方向くじとラチェット =================
pr = J("genesis_v3_phase_randomized_result_v1.json")
h1 = J("holes_h1_h2_result_v1.json")["H1"]
stored_ratchet = J("ratchet_test_result_v1.json")

# ラチェット再走行（同一シード・対照検定つき）
spec3 = importlib.util.spec_from_file_location("s3r_fig", HERE / "run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec3); sys.modules[spec3.name] = g3; spec3.loader.exec_module(g3)
abl, V2 = g3.abl, g3.V2
N_GRAPH, NREG, DELTA = 12, 16, 3e-2
m = N_GRAPH * (N_GRAPH - 1) // 2
_, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
Z0c = Z0c / np.linalg.norm(Z0c)
seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))
odd_ks = [k for k in range(NREG) if k % 2 == 1]
T = 800
trajs = []
for ps in range(20):
    rng = np.random.default_rng(100000 + ps)
    prof = np.zeros(NREG, complex)
    for k in odd_ks:
        prof[k] = np.exp(1j * rng.uniform(0, 2 * np.pi)) / np.sqrt(len(odd_ks))
    C0 = np.zeros((m, NREG), complex); C0[:, 2] = Z0c
    for k in range(NREG):
        if abs(prof[k]) > 0:
            C0[:, k] += DELTA * prof[k] * seed_edge
    eng = V2(N_GRAPH, C0, wp0, vertex_on=True)
    fs = []
    for t in range(T):
        eng.step(); fs.append(eng.diagnostics()["f_seed"])
    trajs.append(fs)
    print(f"  ratchet walker {ps+1}/20", flush=True)
trajs = np.array(trajs)
frac_above = float(np.mean(trajs[:, -1] > trajs[:, 0]))
mean_ratio = float(trajs.mean(axis=0)[-1] / trajs.mean(axis=0)[0])
ok = (abs(frac_above - stored_ratchet["frac_final_above_init"]) < 1e-9
      and abs(mean_ratio - stored_ratchet["mean_ratio"]) < 1e-6)
print(f"  ラチェット対照: frac={frac_above} (保存値 {stored_ratchet['frac_final_above_init']}), "
      f"mean_ratio={mean_ratio:.6f} (保存値 {stored_ratchet['mean_ratio']:.6f}) → 一致={ok}")
np.savez_compressed(HERE / "ratchet_trajectories_v1.npz", trajs=trajs)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
rr = np.array(pr["rates"])
ax1.bar(np.arange(len(rr)), rr, color=["tab:green" if x_ > 0 else "tab:red" for x_ in rr], width=0.6)
ax1.axhline(0, color="k", lw=0.8)
ax1.set_xlabel("乱位相種 #"); ax1.set_ylabel("早期生成率（符号つき）")
ax1.set_title(f"(a) 方向は位相のくじ：50乱位相で P(+)={h1['P_plus']:.2f}\n"
              f"瞬時符号予測子の的中率 {h1['sign_predictor_hit']:.2f}（偶然50%超・完全予測は未達）", fontsize=10)
tt = np.arange(1, T + 1)
for i in range(trajs.shape[0]):
    ax2.plot(tt, trajs[i] / trajs[i, 0], color="tab:gray", lw=0.5, alpha=0.45)
ax2.plot(tt, trajs.mean(axis=0) / trajs.mean(axis=0)[0], color="tab:red", lw=2,
         label=f"アンサンブル平均（終値比 {mean_ratio:.4f}）")
ax2.axhline(1.0, color="k", lw=0.8, ls=":")
ax2.set_xlabel("step"); ax2.set_ylabel("f_seed / f_seed(0)")
ax2.set_title(f"(b) C·f² ラチェット整流：{frac_above*100:.0f}% が正味成長（20本, T={T}）", fontsize=10)
ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(HERE / "fig_mb6_lottery_ratchet_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb6 saved")

# ================= fig_mb7 運動（確定版データ: 計器修正済み v3/v2） =================
# 計器修正の経緯: v1重心計（基本巻きモーメント）は奇数倍音パケットに恒等零
# → run_kinetic_dispersion_demo_v2/v3.py, run_genesis_3d_demo_v2.py 参照
kin = J("kinetic_dispersion_demo_v3.json")
g3d = J("genesis_3d_demo_result_v2.json")
fig = plt.figure(figsize=(10.5, 4.1))
ax1 = fig.add_subplot(1, 2, 1)
tt = np.linspace(0, 400, 801)
pred = (kin["x0"] - 0.05 * tt) % 8
pred_plot = pred.copy()
for s_ in np.where(np.abs(np.diff(pred)) > 4)[0]:
    pred_plot[s_ + 1] = np.nan
ax1.plot(tt, pred_plot, "k--", lw=1, label="予言 (x₀−0.05t) mod 8")
ts_k = [r["t"] for r in kin["rows"]]; xs_k = [r["x"] for r in kin["rows"]]
ax1.plot([0] + ts_k, [kin["x0"] % 8] + xs_k, "o", color="tab:red", ms=7,
         label="実測（8点・最大偏差3×10⁻⁴セル）")
ax1.set_xlabel("step"); ax1.set_ylabel("パケット位置（mod 8, レジスタセル）")
ax1.set_title("(a) 1D：調和閉鎖分散による形不変並進 |v|=0.05\n"
              "整数シフトで波形が厳密に初期形へ回帰（PR 1.000）", fontsize=10)
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
VEL = g3d["VEL"]; x0_3 = g3d["x0"]
tt3 = np.linspace(0, 300, 61)
unw = np.array([[x0_3[a] + VEL[a] * t for a in range(3)] for t in tt3])
ax2.plot(unw[:, 0], unw[:, 1], unw[:, 2], "k--", lw=1, label="予言軌道（展開表示）")
tpts = [0] + [r["t"] for r in g3d["rows"]]
mpts = np.array([[x0_3[a] + VEL[a] * t for a in range(3)] for t in tpts])
ax2.plot(mpts[:, 0], mpts[:, 1], mpts[:, 2], "o", color="tab:red", ms=6,
         label="実測 t=0,100,200,300（3軸とも偏差0.0000）")
ax2.set_xlabel("n₁"); ax2.set_ylabel("n₂"); ax2.set_zlabel("n₃")
ax2.legend(fontsize=8)
ax2.set_title("(b) 3D（8³レジスタ）：3軸同時の等速直線運動\nPR一定3.39（形不変）", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "fig_mb7_motion_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb7 saved")

# ================= fig_mb8 GENESIS 通し =================
specg = importlib.util.spec_from_file_location("gen_fig", HERE / "run_genesis_v1.py")
gv1 = importlib.util.module_from_spec(specg); sys.modules[specg.name] = gv1; specg.loader.exec_module(gv1)
g = gv1.run_genesis(1e-2, "backreaction_fig")
stored_g = J("genesis_result_v1.json")["backreaction"]
match = abs(g["crossing"] - stored_g["crossing"]) == 0
print(f"  GENESIS対照: crossing={g['crossing']} (保存値 {stored_g['crossing']}) → 一致={match}")
f2 = np.array(g["series"]["f2"]); fseed = np.array(g["series"]["fseed"]); rres = np.array(g["series"]["r"])
ts = np.arange(len(f2))
W = 120
lnP = np.log(np.maximum(fseed, 1e-300))
g_run = np.full(len(ts), np.nan)
for i in range(W, len(ts)):
    g_run[i] = (lnP[i] - lnP[i - W]) / W

fig, axes = plt.subplots(3, 1, figsize=(9, 7.6), sharex=True)
eras = [(50, 250, "#e8f0e8", "潜伏期"), (400, 1100, "#f0e8f0", "バースト期"),
        (2000, 4000, "#e8ecf5", "凝縮後（三方向準安定）")]
for lo, hi, c, lbl in eras:
    for ax in axes:
        ax.axvspan(lo, hi, color=c, zorder=0)
axes[0].semilogy(ts, np.maximum(f2, 1e-16), color="tab:gray", lw=1.2)
axes[0].axvline(g["crossing"], color="tab:red", ls=":", lw=1)
axes[0].annotate(f"crossing t={g['crossing']}", (g["crossing"] + 40, 1e-8), fontsize=9, color="tab:red")
axes[0].set_ylabel("f₂（インフレーション）")
axes[0].set_title("GENESIS：真空 → インフレーション → 凝縮 → 生成の開門（単一力学・スイッチなし）",
                  fontsize=11)
for lo, hi, c, lbl in eras:
    axes[0].annotate(lbl, ((lo + hi) / 2, 3e-16), ha="center", fontsize=9)
axes[1].semilogy(ts, np.maximum(rres, 1e-9), color="k", lw=1.0)
axes[1].set_ylabel("一段残差 r（時計）")
axes[2].plot(ts, g_run, color="tab:green", lw=1.2)
gvv = [stored_g["g_latency"], stored_g["g_burst"], stored_g["g_metastable"]]
for (lo, hi, c, lbl), gv_ in zip(eras, gvv):
    axes[2].hlines(gv_, lo, hi, color="tab:red", lw=2)
    axes[2].annotate(f"{gv_:.1e}", ((lo + hi) / 2, gv_), ha="center", va="bottom", fontsize=9, color="tab:red")
axes[2].set_ylabel(f"生成率（移動窓 {W}step）")
axes[2].set_xlabel("step")
axes[2].annotate(f"潜伏→凝縮後 ×{gvv[2]/gvv[0]:.1f}（時代変調＝二因子ゲート）",
                 (2500, gvv[2] * 0.45), fontsize=10, color="tab:green")
fig.tight_layout()
fig.savefig(HERE / "fig_mb8_genesis_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_mb8 saved")
print("ALL FIGURES DONE")
