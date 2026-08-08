#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図U1生成: 番地走査に存在するが公開論文に未記載の観測量（対照実行データのみから）

入力は対照テスト ALL PASS 済みの走査結果
`検証_対照実験/次元の生成構造/番地走査_v1/periodic_address_scan_result_v1.json`
（公開ベースラインと 1522/1522 完全一致・差分0）。走行は行わない・決定論。

4パネル:
 (a) 見かけ次元 n_sig（>5% 特異値の本数）——N=4 で 2 に落ちる
 (b) 親平面外パワー比 perp_ratio——低N域と高N域の二領域構造
 (c) per-step 時計比 ω/(π/72)＝clock_over_pi72 / SAMPLE_EVERY——論文§3の 0.988±0.028 の素データ
 (d) 質量度（Gram 非コヒーレンス 4detΓ/trΓ²）——全域で実質ゼロ

使い方: python3 make_unrecorded_observables_figure_v1.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
SRC = (HERE / "検証_対照実験" / "次元の生成構造" / "番地走査_v1"
       / "periodic_address_scan_result_v1.json")
OUT_PNG = HERE / "fig_u1_unrecorded_observables_v1.png"
OUT_JSON = HERE / "fig_u1_unrecorded_observables_summary_v1.json"

d = json.loads(SRC.read_text())
S = d["scan"]
SE = d["SAMPLE_EVERY"]
N = np.array([r["N"] for r in S])
n_sig = np.array([r["n_sig"] for r in S])
perp = np.array([r["partitions"]["perp_ratio"] for r in S])
axis_sh = np.array([r["partitions"]["axis_share"] for r in S])
clk = np.array([r["clock_over_pi72"] / SE for r in S])
mass = np.array([r["mass_deg"] for r in S])

LOW = N <= 16          # 低N域（少数体・バリオン領域を含む）
HIGH = N >= 40         # 高N域（鋭い凝縮）

fig, axes = plt.subplots(2, 2, figsize=(11, 7.4))

# (a) 見かけ次元
ax = axes[0, 0]
ax.plot(N, n_sig, "o-", ms=5, lw=1.2, color="tab:blue")
ax.axhline(3, color="k", ls="--", lw=0.8, label="3（空間3方向）")
ax.annotate("N=4: n_sig=2\n第3方向が立たない", xy=(4, 2), xytext=(6, 4.2),
            fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xscale("log"); ax.set_xlabel("N（関係波数＝実効エネルギー）")
ax.set_ylabel("n_sig（>5% 特異値の本数）")
ax.set_title("(a) 見かけ次元：N=4 だけ 2 に落ちる", fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# (b) 親平面外パワー比＝二領域構造
ax = axes[0, 1]
ax.plot(N[LOW], perp[LOW], "o", ms=6, color="tab:red", label=f"N 16以下（平均 {perp[LOW].mean():.3f}）")
ax.plot(N[~LOW & ~HIGH], perp[~LOW & ~HIGH], "^", ms=6, color="tab:gray", label="N 20〜36（遷移）")
ax.plot(N[HIGH], perp[HIGH], "s", ms=6, color="tab:blue", label=f"N 40以上（平均 {perp[HIGH].mean():.3f}）")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("N"); ax.set_ylabel("perp_ratio ＝ 親平面外パワー / 全パワー")
ax.set_title("(b) 二領域構造：高N で凝縮が鋭くなる（0.630 → 0.093）", fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

# (c) per-step 時計比
ax = axes[1, 0]
ax.axhline(1.0, color="k", ls="--", lw=0.8, label="ω=π/72（一周=144ステップ）")
ax.plot(N, clk, "o", ms=5, color="tab:green")
ax.axhspan(clk.mean() - clk.std(), clk.mean() + clk.std(), color="tab:green", alpha=0.12,
           label=f"平均±標準偏差 {clk.mean():.4f}±{clk.std():.4f}")
ax.set_xscale("log"); ax.set_xlabel("N"); ax.set_ylabel("ω_clock / (π/72)  [per step]")
ax.set_title("(c) per-step 時計比：÷SAMPLE_EVERY(=5) が正規化の実体", fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# (d) 質量度
ax = axes[1, 1]
ax.plot(N, mass, "o-", ms=5, lw=1.0, color="tab:purple")
ax.axhline(0.0, color="k", ls="--", lw=0.8)
ax.set_xscale("log"); ax.set_xlabel("N")
ax.set_ylabel("質量度 4detΓ/(trΓ)²")
ax.set_title(f"(d) 質量度は全域で実質ゼロ（最大 {mass.max():.2e}）", fontsize=11)
ax.grid(alpha=0.3)

fig.suptitle("図U1  番地走査に存在するが公開論文に未記載の観測量（対照実行データ・28点）",
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT_PNG, dpi=150)
plt.close(fig)

summary = {
    "入力": str(SRC.relative_to(HERE)),
    "点数": len(N),
    "SAMPLE_EVERY": SE,
    "per_step時計比": {"平均": float(clk.mean()), "標準偏差": float(clk.std()),
                       "最小": float(clk.min()), "最大": float(clk.max())},
    "n_sig": {"N=4": int(n_sig[0]), "N=5": int(n_sig[1]),
              "N>=6_最小": int(n_sig[2:].min()), "N>=6_最大": int(n_sig[2:].max())},
    "perp_ratio": {"N<=16_平均": float(perp[LOW].mean()),
                   "N>=40_平均": float(perp[HIGH].mean()),
                   "N=144": float(perp[-1])},
    "axis_share": {"最小": float(axis_sh.min()), "最大": float(axis_sh.max())},
    "質量度": {"最小": float(mass.min()), "最大": float(mass.max())},
}
OUT_JSON.write_text(json.dumps(summary, indent=1, ensure_ascii=False))
print(json.dumps(summary, indent=1, ensure_ascii=False))
print(f"→ {OUT_PNG.name} / {OUT_JSON.name}")
