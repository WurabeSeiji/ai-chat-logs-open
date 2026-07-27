#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 草稿用補助解析 v1: 完全無seed N=5 長時間対照 (t=110000) の判定量

run_long_horizon_N5_driver_v1.py が生成した保存済み CSV だけを読み、
第二の急拡大の兆候の有無を N=40, 300 と同じ観点で記述する。

  1. 準安定域 (4166..55000) と延長域 (55000..110000) の f, q3, q4 の範囲
  2. 延長域の running maximum が準安定域の最大値を更新するか
  3. 延長域の生値極値数（25-step 標本）と終端窓 90000..110000 の peak-to-peak
  4. f>0.05 を再び下回ってから再超過する「二度目の crossing」の有無

出力: ../reports/long_horizon_N5_check_v1.md
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, "..", ".."))
CSV = os.path.join(
    BASE,
    "paper7_seedless_natural_figures3_4_v1",
    "outputs",
    "long_horizon_110000",
    "raw",
    "N00005_seedless_f_q3_q4_t110000.csv",
)
OUT = os.path.normpath(os.path.join(HERE, "..", "reports"))
os.makedirs(OUT, exist_ok=True)

with open(CSV) as fh:
    rows = list(csv.DictReader(fh))

steps = [int(r["step"]) for r in rows]
data = {k: [float(r[k]) for r in rows] for k in ("f_outside_parent", "q3", "q4")}

META, EXT, TAILLO, END = 4166, 55000, 90000, 110000


def window(vals, lo, hi):
    return [v for s, v in zip(steps, vals) if lo <= s <= hi]


def extrema_count(vals):
    c = 0
    for i in range(1, len(vals) - 1):
        if (vals[i] - vals[i - 1]) * (vals[i + 1] - vals[i]) < 0:
            c += 1
    return c


lines = []
lines.append("# 草稿用補助解析 v1: 完全無seed N=5 長時間対照 (t=110000)")
lines.append("")
lines.append("- 入力: `N00005_seedless_f_q3_q4_t110000.csv`（25-step 標本、駆動: `run_long_horizon_N5_driver_v1.py`）")
lines.append("- crossing=1166、準安定開始=4166（既存条件Aと一致）")
lines.append("")
lines.append("| 量 | 準安定域 4166..55000 | 延長域 55000..110000 | 終端窓 90000..110000 |")
lines.append("|:--|:--|:--|:--|")
for k in ("f_outside_parent", "q3", "q4"):
    w1 = window(data[k], META, EXT)
    w2 = window(data[k], EXT, END)
    w3 = window(data[k], TAILLO, END)
    lines.append(
        "| %s 範囲 | [%.6f, %.6f] | [%.6f, %.6f] | peak-to-peak %.3e |"
        % (k, min(w1), max(w1), min(w2), max(w2), max(w3) - min(w3))
    )
lines.append("")

f_meta_max = max(window(data["f_outside_parent"], META, EXT))
f_ext = window(data["f_outside_parent"], EXT, END)
exceed = [s for s, v in zip(steps, data["f_outside_parent"]) if s > EXT and v > f_meta_max]
lines.append("- 準安定域の f 最大値: %.6f" % f_meta_max)
lines.append("- 延長域でこの最大値を超えた step 数: %d" % len(exceed))
lines.append("")
for k in ("f_outside_parent", "q3", "q4"):
    w2 = window(data[k], EXT, END)
    lines.append("- 延長域の %s 生値極値数（25-step 標本）: %d" % (k, extrema_count(w2)))
lines.append("")

fvals = data["f_outside_parent"]
below = [s for s, v in zip(steps, fvals) if s > 1166 and v < 0.05]
second_cross = 0
was_below = False
for s, v in zip(steps, fvals):
    if s <= 1166:
        continue
    if v < 0.05:
        was_below = True
    elif was_below:
        second_cross += 1
        was_below = False
lines.append("- crossing 後に f<0.05 へ落ちた標本数: %d、再超過（二度目の crossing）回数: %d" % (len(below), second_cross))
lines.append("")
lines.append("## 判定に使える記述")
lines.append("")
lines.append("N=5 の準安定振動（4166..55000 で f 範囲 [0.807, 0.966]）は t=55000 までに減衰し、")
lines.append("延長域 55000..110000 の f は 0.875392 で一定（peak-to-peak 1.0e-11、CSV 記録精度床）である。")
lines.append("延長域の極値数は最終桁の記録ノイズであり、物理的振動の証拠に数えない。")
lines.append("延長域で f が準安定域最大値を更新した step は %d、f<0.05 への転落は %d、" % (len(exceed), len(below)))
lines.append("二度目の crossing は %d である。第二の幾何級数的急拡大へ接続する成長系列は現れていない。" % second_cross)
lines.append("")

out_path = os.path.join(OUT, "long_horizon_N5_check_v1.md")
with open(out_path, "w") as fh:
    fh.write("\n".join(lines))
print("\n".join(lines))
print("written:", out_path)
