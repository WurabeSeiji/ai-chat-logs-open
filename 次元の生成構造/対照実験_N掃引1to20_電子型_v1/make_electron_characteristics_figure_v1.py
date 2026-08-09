#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図E: 電子の特性値（E1 巻き集中度・E2 電荷）の推移——時間推移と N 依存性

入力（走行済みデータのみ・再走行なし・決定論）:
  result_nsweep_electron_v1.json（準安定窓の中央値）
  nsweep_electron_N{n}_v1.npz（全步時系列 hair_m_* / hair_v_*）
  参照_中性_result_tb_nsweep_1to20_v1.json（比較用・空間側のみ）

パネル:
 (a) E1 集中度 P(m*=-3)/ΣP の時間推移（全 N 重ね描き・色=N）
 (b) E2 優勢巻き m̂ の時間推移（全 N 重ね描き）
 (c) 巻きが確定する時刻 τ_wind（集中度が初めて 0.9 を超えた步）の N 依存性
     ——空間 τ_space・時間 τ_time と同じ軸で比較
 (d) E1 の N 依存性（準安定窓中央値・k=3 帯と全帯）
 (e) E2 の N 依存性（Q̂ と可読率）
 (f) 相棒帯 k=3 のパワー（物質 vs 真空対照）の N 依存性

使い方: python3 make_electron_characteristics_figure_v1.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, font_manager

for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
RES = json.loads((HERE / "result_nsweep_electron_v1.json").read_text())
NEU = json.loads((HERE / "参照_中性_result_tb_nsweep_1to20_v1.json").read_text())
WIN = RES["env"]["window"]
M_STAR = RES["env"]["recipe"]["m_star"]
PK = RES["env"]["recipe"]["partner_band_k"]
THRESH = 0.9

built = sorted(int(k) for k, v in RES["N"].items() if v.get("built"))
series = {}
for n in built:
    p = HERE / f"nsweep_electron_N{n}_v1.npz"
    if p.exists():
        z = np.load(p)
        series[n] = {"conc": z["hair_m_conc_k3"], "dom": z["hair_m_dom_m"],
                     "pw": z["hair_m_k3_power"], "pw_v": z["hair_v_k3_power"]}


def first_over(x, thr):
    idx = np.flatnonzero(np.nan_to_num(x, nan=-1.0) > thr)
    return int(idx[0]) + 1 if len(idx) else None


tau_wind = {n: first_over(s["conc"], THRESH) for n, s in series.items()}
colors = {n: cm.viridis((i + 0.5) / max(len(series), 1))
          for i, n in enumerate(sorted(series))}

fig, ax = plt.subplots(3, 2, figsize=(13, 12))

# (a) 集中度の時間推移
a = ax[0, 0]
for n, s in sorted(series.items()):
    a.plot(np.arange(1, len(s["conc"]) + 1), s["conc"], lw=0.7,
           color=colors[n], label=f"N={n}")
a.axhline(THRESH, color="tab:red", lw=0.9, ls=":", label=f"{THRESH}（V2a 閾値）")
a.axvspan(WIN[0], WIN[1], color="green", alpha=0.06)
a.set_xlabel("τ（step）"); a.set_ylabel(f"E1 集中度 P(m*={M_STAR})/ΣP  [k={PK}帯]")
a.set_title("(a) E1 巻き集中度の時間推移（全N）", fontsize=11)
a.legend(fontsize=6, ncol=3)

# (b) 優勢巻きの時間推移
a = ax[0, 1]
for n, s in sorted(series.items()):
    a.plot(np.arange(1, len(s["dom"]) + 1), s["dom"], lw=0.7, color=colors[n])
a.axhline(M_STAR, color="tab:red", lw=0.9, ls=":", label=f"m*={M_STAR}")
a.axvspan(WIN[0], WIN[1], color="green", alpha=0.06)
a.set_xlabel("τ（step）"); a.set_ylabel("E2 優勢巻き m̂（符号付き）")
a.set_title("(b) E2 優勢巻きの時間推移（全N）", fontsize=11)
a.legend(fontsize=8)

# (c) 巻き確定時刻の N 依存性
a = ax[1, 0]
ns_ = sorted(series)
a.plot(ns_, [tau_wind[n] or np.nan for n in ns_], "o-", color="tab:purple",
       label=f"τ_wind（集中度>{THRESH}）")
a.plot(ns_, [RES["N"][str(n)]["tau_space"] or np.nan for n in ns_], "s--",
       color="tab:blue", label="τ_space（f₂>0.05）")
a.plot(ns_, [RES["N"][str(n)]["tau_time"] or np.nan for n in ns_], "^--",
       color="tab:red", label="τ_time（時計取得）")
a.set_xlabel("N"); a.set_ylabel("誕生時刻 τ（step）")
a.set_title("(c) 巻きはいつ確定するか（空間・時間との比較）", fontsize=11)
a.legend(fontsize=8)

# (d) E1 の N 依存性
a = ax[1, 1]
a.plot(ns_, [RES["N"][str(n)]["conc_k3_med"] for n in ns_], "o-",
       label=f"k={PK}帯（準安定窓中央値）")
a.plot(ns_, [RES["N"][str(n)]["conc_all_med"] for n in ns_], "s--",
       color="tab:orange", label="全帯")
a.axhline(THRESH, color="tab:red", lw=0.9, ls=":")
a.set_ylim(-0.02, 1.05)
a.set_xlabel("N"); a.set_ylabel("E1 集中度")
a.set_title("(d) E1 の N 依存性", fontsize=11); a.legend(fontsize=8)

# (e) E2 の N 依存性
a = ax[2, 0]
a.plot(ns_, [RES["N"][str(n)]["q_hat_med"] for n in ns_], "o-",
       color="tab:green", label="Q̂ = m̂/3（窓中央値）")
a.axhline(M_STAR / 3.0, color="tab:red", lw=0.9, ls=":",
          label=f"Q = {M_STAR}/3 = {M_STAR/3:+.0f}")
a2 = a.twinx()
a2.plot(ns_, [RES["N"][str(n)]["readable_rate"] for n in ns_], "^--",
        color="tab:gray", label="可読率（3|m̂）")
a2.set_ylabel("可読率"); a2.set_ylim(-0.05, 1.1)
a.set_xlabel("N"); a.set_ylabel("Q̂")
a.set_title("(e) E2 電荷と単独可読性の N 依存性", fontsize=11)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")

# (f) 相棒帯パワー（物質 vs 真空）
a = ax[2, 1]
a.semilogy(ns_, [max(RES["N"][str(n)]["k3_power_med"], 1e-40) for n in ns_],
           "o-", label="物質 δ=1e−2")
a.semilogy(ns_, [max(RES["N"][str(n)]["k3_power_med_vacuum"], 1e-40) for n in ns_],
           "k--s", label="真空対照 δ=0")
a.set_xlabel("N"); a.set_ylabel(f"k={PK}帯 パワー（窓中央値）")
a.set_title("(f) 相棒帯パワー：物質と真空対照", fontsize=11); a.legend(fontsize=8)

fig.suptitle(f"図E  電子型フェルミオン（m*={M_STAR}・Q={M_STAR/3:+.0f}）の特性値の推移"
             f"——E1 巻き集中度・E2 電荷（N={min(ns_)}〜{max(ns_)}・"
             f"構成成功 {len(ns_)} 点）", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.965))
out = HERE / "fig_electron_characteristics_v1.png"
fig.savefig(out, dpi=140)
plt.close(fig)

summary = {"図": out.name, "構成成功N": ns_,
           "tau_wind": {str(n): tau_wind[n] for n in ns_},
           "conc_k3_med": {str(n): RES["N"][str(n)]["conc_k3_med"] for n in ns_},
           "q_hat_med": {str(n): RES["N"][str(n)]["q_hat_med"] for n in ns_},
           "readable_rate": {str(n): RES["N"][str(n)]["readable_rate"] for n in ns_},
           "閾値": THRESH}
(HERE / "fig_electron_characteristics_summary_v1.json").write_text(
    json.dumps(summary, indent=1, ensure_ascii=False))
print(f"→ {out.name}")
print("N        :", ns_)
print("τ_wind   :", [tau_wind[n] for n in ns_])
print("集中度   :", [round(RES["N"][str(n)]["conc_k3_med"], 4) for n in ns_])
print("Q̂        :", [round(RES["N"][str(n)]["q_hat_med"], 3) for n in ns_])
print("可読率   :", [round(RES["N"][str(n)]["readable_rate"], 2) for n in ns_])
