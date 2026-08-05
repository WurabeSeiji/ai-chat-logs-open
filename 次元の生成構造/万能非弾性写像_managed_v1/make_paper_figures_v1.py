#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""二体論文用の図3枚＋census RK4再測定

図1 点火則: RK4掃引の lnP_f 軌跡と C=rate/f0²（RK4 vs RK2較正）
図2 運命: f(j) 往復型→統計的平衡（RK4、RK2重ね描き）
図3 対構造census: bin別成長の排他性＋毛コヒーレンス q走査（RK4再測定）

census再測定の判定（実行前固定）: P1/P2/P3 が RK2 実測と同一の合否・
q*=+4・コヒーレンス±0.05 で一致すること（積分器非依存の確認）。
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = ["Hiragino Sans", "Arial Unicode MS", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


rk4 = load("rk4_fig", HERE / "run_ignition_fate_rk4_v2.py")
pc3 = load("pc3_fig", HERE / "run_pair_structure_census_v3.py")
pc = sys.modules["pc2_for3"]
v3, v1, toy, base = rk4.v3, rk4.v1, rk4.toy, rk4.base

params = base.Params(high_n=63, recursive_collision_count=200)
sp = base.build_source_params(params)

# ================= 図1 点火則 =================
S = 8.0
series = {}
for seed_amp in (1e-3, 1e-2, 1e-1):
    a = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    a = a + v1.make_bundle(sp, v1.ODD_KS, "A", scale=seed_amp * S)
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
    fs = []
    for _ in range(200):
        a, b, _ = rk4.collision_step_rk4(a, b, sp)
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        fs.append(v1.fermionic_power_raw(a, b, sp) / tot)
    series[seed_amp] = np.array(fs)

rk4_sweep = json.loads((HERE / "ignition_fate_rk4_result_v2.json").read_text())["C_sweep"]
rk2_sweep = json.loads((HERE / "seed_fraction_sweep_result_v3.json").read_text())["rows"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
for seed_amp, fs in series.items():
    y = (np.log(fs) - np.log(fs[0])) / fs[0] ** 2
    ax1.plot(np.arange(1, 201), y,
             label=f"種振幅比 {seed_amp:g}（f₀={fs[0]:.1e}）")
jj = np.arange(1, 201)
ax1.plot(jj, 11.45 * jj, "k--", lw=1, label="C·j（C=11.45）")
ax1.set_xlabel("衝突回数 j")
ax1.set_ylabel("(ln f − ln f₀) / f₀²")
ax1.set_title("(a) 成長軌跡の収束（4桁の f₀ が一本に重なる）")
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

f0s_rk4 = [r["f0"] for r in rk4_sweep]
Cs_rk4 = [r["C"] for r in rk4_sweep]
f0s_rk2 = [r["f0"] for r in rk2_sweep]
Cs_rk2 = [r["C"] for r in rk2_sweep]
ax2.semilogx(f0s_rk4, Cs_rk4, "o-", color="tab:red", label="RK4: C=11.45（採用値）")
ax2.semilogx(f0s_rk2, Cs_rk2, "s--", mfc="none", color="tab:gray",
             label="RK2: C=10.4（積分バイアス約10%）")
ax2.set_ylim(9, 13)
ax2.set_xlabel("初期種割合 f₀")
ax2.set_ylabel("C = rate₀ / f₀²")
ax2.set_title("(b) 点火則 rate = C·f² の係数一定性")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(HERE / "fig1_ignition_law_v1.png", dpi=160)
plt.close(fig)
print("fig1 saved")

# ================= 図2 運命 =================
d4 = json.loads((HERE / "ignition_fate_rk4_result_v2.json").read_text())
d2 = json.loads((HERE / "ignition_fate_result_v1.json").read_text())
j4 = np.arange(len(d4["f_series_every10"])) * 10
j2 = np.arange(len(d2["f_series_every10"])) * 10

fig, ax = plt.subplots(figsize=(8, 3.8))
ax.plot(j2, d2["f_series_every10"], color="tab:gray", lw=1, alpha=0.6,
        label="RK2（v1・ドリフト2.9e-3）")
ax.plot(j4, d4["f_series_every10"], color="tab:blue", lw=1.5,
        label="RK4（本稿・ドリフト2.8e-9）")
ax.axhline(0.4690, color="tab:red", ls="--", lw=1,
           label="統計的平衡 f* = 0.4690（後半500平均）")
ax.axhline(0.494140625, color="tab:green", ls=":", lw=1,
           label="マスク位相空間割合 0.494（等分配候補値）")
ax.set_xlabel("衝突回数 j")
ax.set_ylabel("フェルミオン関係量比 f")
ax.set_title("点火後の運命：往復型 → 統計的平衡（S=8, 種振幅比0.1, 3000衝突）")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(HERE / "fig2_fate_v1.png", dpi=160)
plt.close(fig)
print("fig2 saved")

# ================= 図3 census（RK4再測定） =================
n, ne = sp.chi_grid_n, sp.eta_grid_n
Sc, PUMP_KS, SEED_K, SEED_AMP, J_WINDOW = 8.0, (30, 32, 34), 21, 0.2, 40
a0 = pc.single_winding(v1.make_bundle(sp, PUMP_KS, "A", scale=1.0), sp) * Sc
b0 = pc.single_winding(v1.make_bundle(sp, PUMP_KS, "B", scale=1.0), sp) * Sc
seed = pc.single_winding(v1.make_bundle(sp, (SEED_K,), "A", scale=1.0), sp) * (SEED_AMP * Sc)
p_pump = pc.power_bins(a0, b0, sp)
pump_bins = sorted(int(k) for k in range(n) if p_pump[k] > 1e-6)
r_s = int(np.argmax(pc.power_bins(seed, np.zeros_like(seed), sp)))
pump_sums = sorted(set((k1 + k2) % n for k1 in pump_bins for k2 in pump_bins))
partner_band = sorted(set((sm - r_s) % n for sm in pump_sums))
pump_diffs = sorted(set(abs(k1 - k2) for k1 in pump_bins for k2 in pump_bins if k1 != k2))
sideband = sorted(set((r_s + d) % n for d in pump_diffs) | set((r_s - d) % n for d in pump_diffs))

a, b = a0 + seed, b0
p0 = pc.power_bins(a, b, sp)
partner_p0 = float(sum(p0[k] for k in partner_band))
for _ in range(J_WINDOW):
    a, b, _ = rk4.collision_step_rk4(a, b, sp)
p1 = pc.power_bins(a, b, sp)
growth = np.maximum(p1 - p0, 0.0)

pred_all = set(partner_band) | set(sideband) | set(pump_sums) | set(pump_bins) \
    | set(range(max(0, r_s - 2), min(n, r_s + 3)))
others = [k for k in range(n) if k not in pred_all and p0[k] < 1e-12]
med_other = float(np.median(growth[others]))
mean_partner = float(np.mean([growth[k] for k in partner_band]))
fa, fb = pc.chi_spectra(a, sp), pc.chi_spectra(b, sp)
r_p = int(partner_band[int(np.argmax([p1[k] for k in partner_band]))])
cohs = {int(q): pc3.hair_coherence(fa, fb, ne, r_s, r_p, q) for q in range(-6, 7)}
q_star = max(cohs, key=cohs.get)
g_partner = float(sum(growth[k] for k in partner_band))
g_side = float(sum(growth[k] for k in sideband))
ratio = g_partner / max(g_side, 1e-300)
p1_pass = med_other < 1e-15 * mean_partner
p2_pass = cohs[q_star] > 0.5 and q_star == 4 and cohs[0] < 0.1
p3_pass = partner_p0 < 1e-25 and 0.1 <= ratio <= 10.0
print(f"census RK4: P1 中央値={med_other:.1e} 相棒平均={mean_partner:.1e} → {p1_pass}")
print(f"            P2 q*={q_star:+d} coh={cohs[q_star]:.4f} q0={cohs[0]:.1e} → {p2_pass}")
print(f"            P3 初期={partner_p0:.1e} 比={ratio:.3f} → {p3_pass}")
json.dump({"P1": {"median_other": med_other, "mean_partner": mean_partner,
                    "n_other": len(others), "pass": bool(p1_pass)},
           "P2": {"q_star": int(q_star), "coh_max": float(cohs[q_star]),
                    "coh_q0": float(cohs[0]), "pass": bool(p2_pass),
                    "coherence_by_q": {str(q): float(c) for q, c in cohs.items()}},
           "P3": {"partner_initial": partner_p0, "ratio": ratio, "pass": bool(p3_pass)},
           "integrator": "RK4"},
          open(HERE / "pair_structure_census_rk4_result_v1.json", "w"),
          ensure_ascii=False, indent=2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
ks = np.arange(n)
floor = 1e-34
gplot = np.maximum(growth, floor)
ax1.semilogy(ks, gplot, color="lightgray", lw=0.6, zorder=1)
sets = [(others, "lightgray", "予言外の空bin（485個・機械ゼロ）", "."),
        (sideband, "tab:blue", "XPM側帯（予言）", "o"),
        (pump_sums, "tab:green", "ポンプ和帯（予言）", "^"),
        (partner_band, "tab:red", "相棒帯（対生成の予言bin）", "s")]
for bins, c, lbl, m in sets:
    ax1.semilogy(bins, gplot[list(bins)], m, color=c, ms=4, label=lbl, zorder=3)
ax1.axvline(r_s, color="k", lw=0.8, ls=":", label=f"シード bin {r_s}")
ax1.set_xlabel("生スペクトル bin k")
ax1.set_ylabel("パワー成長（40衝突, 下限1e-34）")
ax1.set_title("(a) 和則の排他性：予言binのみが成長")
ax1.set_xlim(0, 130)
ax1.legend(fontsize=7, loc="upper right")

qs = sorted(cohs)
ax2.bar([str(q) for q in qs], [max(cohs[q], 1e-16) for q in qs],
        color=["tab:red" if q == q_star else "tab:gray" for q in qs])
ax2.set_yscale("log")
ax2.set_ylim(1e-16, 3)
ax2.set_xlabel("毛巻き数オフセット q")
ax2.set_ylabel("異常対相関コヒーレンス")
ax2.set_title(f"(b) 毛の帳簿：q*={q_star:+d} で {cohs[q_star]:.2f}、他は選択則で消滅")
fig.tight_layout()
fig.savefig(HERE / "fig3_census_v1.png", dpi=160)
plt.close(fig)
print("fig3 saved")
