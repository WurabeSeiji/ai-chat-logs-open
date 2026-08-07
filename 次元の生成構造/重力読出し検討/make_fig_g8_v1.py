#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図G8生成（日英）: ゲージ/重力の同時分離（result_g9のJSONのみから）"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

JA_FONT = None
for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        JA_FONT = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent
d = json.loads((HERE / "result_g9_gauge_gravity_separation_v1.json").read_text())

T = {"ja": dict(
      title="図G8  決定実験: 同一読出しからのゲージ応答と重力応答の分離",
      y="位相平均結合 Ē", note=("全対で引力（重力=電荷盲目）／E(++)−E(−−)=7e-10（共役厳密）／\n"
                             "E(++)−E(+−)=−1.43（ゲージ=電荷構造依存チャネル）／中性プローブも重力"),
      sfx="_ja"),
     "en": dict(
      title="Fig. G8  Decisive test: separation of gauge and gravity responses from one readout",
      y="phase-averaged coupling Ē",
      note=("attraction for all pairs (gravity = charge-blind) / E(++)−E(−−)=7e-10 (exact conjugation) /\n"
            "E(++)−E(+−)=−1.43 (gauge = charge-structure channel) / neutral probe still gravitates"),
      sfx="_en")}
keys = ["0,0", "+1,+1", "-1,-1", "+1,-1", "+1,0"]
labels = ["(0,0)", "(+1,+1)", "(−1,−1)", "(+1,−1)", "(+1,0)"]
Eb = [d["pairs"][k]["E_bar"] for k in keys]
Er = [d["pairs"][k]["spread"] for k in keys]
for lang, t in T.items():
    plt.rcParams["font.family"] = JA_FONT if lang == "ja" else "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    cols = ["tab:green", "tab:blue", "tab:blue", "tab:orange", "tab:purple"]
    ax.bar(range(5), Eb, yerr=Er, color=cols, alpha=0.85, capsize=4)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(5)); ax.set_xticklabels(labels)
    ax.set_ylabel(t["y"]); ax.set_title(t["title"], fontsize=10.5)
    ax.text(0.02, 0.05, t["note"], transform=ax.transAxes, fontsize=8.2,
            bbox=dict(boxstyle="round", fc="#f2f2f2"))
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(HERE / f"fig_g8_gauge_gravity{t['sfx']}_v1.png", dpi=150)
    plt.close(fig)
print("fig G8 ja/en done")
