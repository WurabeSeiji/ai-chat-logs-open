#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 草稿用補助解析 v1

Stage A2a の保存済み正本CSV（完全無seed N=5 float64, exec 1）だけを読み、
論文本文で引用する以下の量を再計算して固定する。新しい時間発展は行わない。

  1. f の毎step差分の符号分布（step 0 から f>=1e-2 初回到達の直前まで）
     - 「成長→停止→下位準安定棚→再成長」区間の有無の直接検査
  2. log10(f) の線形回帰（f が数値床を離れた 1e-23 初回到達step以降、
     1e-2 初回到達stepまで）: 傾き lambda（自然対数/step）と R^2
  3. 無seed軌道での rank_q=4 の初出step（q_timeseries 保存行, 5 step間隔）と、
     その時点の f・direction 3/4 占有（占有は挟む実保存行）

入力（読み取り専用・Stage A2a 正本）:
  ../../paper8_stage_A2a_seedless_N5/raw/A2a_N5_seedless_f64_e1/f_timeseries.csv
  ../../paper8_stage_A2a_seedless_N5/raw/A2a_N5_seedless_f64_e1/q_timeseries.csv
  ../../paper8_stage_A2a_seedless_N5/raw/A2a_N5_seedless_f64_e1/occupation_timeseries.csv
  ../../paper8_stage_A2a_seedless_N5/processed/seedless_first_passage_levels.csv

出力: ../reports/monotonicity_and_regression_check_v1.md
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(BASE, "paper8_stage_A2a_seedless_N5", "raw", "A2a_N5_seedless_f64_e1")
PROC = os.path.join(BASE, "paper8_stage_A2a_seedless_N5", "processed")
OUT = os.path.normpath(os.path.join(HERE, "..", "reports"))
os.makedirs(OUT, exist_ok=True)


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


f_rows = read_csv(os.path.join(RAW, "f_timeseries.csv"))
q_rows = read_csv(os.path.join(RAW, "q_timeseries.csv"))
occ_rows = read_csv(os.path.join(RAW, "occupation_timeseries.csv"))
fp_rows = read_csv(os.path.join(PROC, "seedless_first_passage_levels.csv"))

steps = [int(r["step"]) for r in f_rows]
fvals = [float(r["f"]) for r in f_rows]
assert steps == list(range(len(steps))), "f_timeseries must be 1-step contiguous"

fp = {r["level_label"]: int(r["first_passage_step"]) for r in fp_rows if r["status"] == "found"}
step_1e2 = fp["1e-02"]
step_1e23 = fp["1e-23"]

# 1. 毎step差分の符号（0 .. step_1e2-1 の差分）
pos = neg = zero = 0
min_diff = None
for t in range(step_1e2):
    d = fvals[t + 1] - fvals[t]
    if d > 0:
        pos += 1
        min_diff = d if min_diff is None else min(min_diff, d)
    elif d < 0:
        neg += 1
    else:
        zero += 1

# 2. log(f) 回帰（step_1e23 .. step_1e2）
xs = list(range(step_1e23, step_1e2 + 1))
ys = [math.log(fvals[t]) for t in xs]
n = len(xs)
mx = sum(xs) / n
my = sum(ys) / n
sxx = sum((x - mx) ** 2 for x in xs)
sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
lam = sxy / sxx
b0 = my - lam * mx
ss_res = sum((y - (b0 + lam * x)) ** 2 for x, y in zip(xs, ys))
ss_tot = sum((y - my) ** 2 for y in ys)
r2 = 1.0 - ss_res / ss_tot

# 3. rank_q=4 初出（保存行, 5 step間隔）
first_r4 = None
for r in q_rows:
    if int(float(r["rank_q"])) >= 4:
        first_r4 = r
        break
r4_step = int(first_r4["step"])
r4_f = fvals[r4_step]
occ_before = max((r for r in occ_rows if int(r["step"]) <= r4_step), key=lambda r: int(r["step"]))
occ_after = min((r for r in occ_rows if int(r["step"]) >= r4_step), key=lambda r: int(r["step"]))

lines = []
lines.append("# 草稿用補助解析 v1: 単調性・回帰・rank_q初出（完全無seed N=5, exec 1）")
lines.append("")
lines.append("入力は Stage A2a 保存済み正本CSVのみ。新しい時間発展・seed注入は行っていない。")
lines.append("")
lines.append("## 1. f の毎step差分の符号（step 0 〜 %d, 差分 %d 本）" % (step_1e2, step_1e2))
lines.append("")
lines.append("- 正の差分: %d" % pos)
lines.append("- 負の差分: %d" % neg)
lines.append("- ゼロ差分: %d" % zero)
lines.append("- 最小の正差分: %.6e" % min_diff)
lines.append("- f(0)=%.17e, f(%d)=%.17e" % (fvals[0], step_1e2, fvals[step_1e2]))
lines.append("")
lines.append("## 2. log(f) 線形回帰（step %d 〜 %d, %d 点）" % (step_1e23, step_1e2, n))
lines.append("")
lines.append("- lambda（自然対数/step）: %.8f" % lam)
lines.append("- R^2: %.8f" % r2)
lines.append("- 対応する1 stepあたり振幅比 exp(lambda/2): %.7f" % math.exp(lam / 2))
lines.append("- 対応する f 比 exp(lambda): %.7f" % math.exp(lam))
lines.append("")
lines.append("## 3. rank_q=4 の初出（q保存行, 5 step間隔）")
lines.append("")
lines.append("- 初出保存step: %d" % r4_step)
lines.append("- その時点の f: %.17e" % r4_f)
lines.append("- q3/q1: %s, q4/q1: %s" % (first_r4["q3_over_q1"], first_r4["q4_over_q1"]))
lines.append(
    "- 挟む占有実保存行 step %s: d3=%s, d4=%s"
    % (occ_before["step"], occ_before["direction_3_occupation"], occ_before["direction_4_occupation"])
)
lines.append(
    "- 挟む占有実保存行 step %s: d3=%s, d4=%s"
    % (occ_after["step"], occ_after["direction_3_occupation"], occ_after["direction_4_occupation"])
)
lines.append("")

out_path = os.path.join(OUT, "monotonicity_and_regression_check_v1.md")
with open(out_path, "w") as fh:
    fh.write("\n".join(lines))
print("\n".join(lines))
print("written:", out_path)
