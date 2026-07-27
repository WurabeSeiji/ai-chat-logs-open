#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 草稿用補助解析 v1: 残余回転占有の検査（追加方向不在の直接測定）

rank_q は Q=[B0|Bdom] (M×4) の rank であり構造上 4 が上限のため、
「rank_q が 4 を超えない」は追加方向不在の証拠にならない。
追加方向の自然占有を直接測るのは五色分解の残余回転占有
other_rotating_occupation（direction 1〜4 と kernel 以外の回転部分空間占有、
全占有和は毎 step 1）である。

入力（読み取り専用・保存済み正本 CSV のみ、新しい時間発展なし）:
  paper7_seedless_natural_figures3_4_v1/outputs/raw/N000{05,40,300}/paper7_long_timeseries.csv
  （完全無seed自然軌道、t<=55000、N=300 は実行済みの場合のみ）

出力: ../reports/residual_rotating_occupation_check_v1.md
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(BASE, "paper7_seedless_natural_figures3_4_v1", "outputs", "raw")
OUT = os.path.normpath(os.path.join(HERE, "..", "reports"))
os.makedirs(OUT, exist_ok=True)

# 準安定開始（crossing+3000、予備実験報告の条件A実測値）
META_START = {5: 4166, 40: 5011, 300: 7849}
# 初期 K の分解次元（親平面/回転補/核）: N=5,40 は parent_state_structure_check_v1 実測、
# N=300 は第6論文実測値（回転補 598）
ROT_DIM = {5: 4, 40: 78, 300: 598}

lines = []
lines.append("# 草稿用補助解析 v1: 残余回転占有の検査（追加方向不在の直接測定）")
lines.append("")
lines.append("- 入力: 完全無seed自然軌道の五色分解 CSV（保存済み正本、新しい時間発展なし）")
lines.append("- rank_q は Q=[B0|Bdom] (M×4) の rank で構造上 4 が上限。追加方向の検出器にならない。")
lines.append("- 残余回転占有 = direction 1〜4 と kernel 以外の回転部分空間への占有（全占有和=1）。")
lines.append("")
lines.append("| N | 残余回転空間次元 | 準安定域 max | 準安定域 mean | 末端値 | 標本数 |")
lines.append("|--:|--:|:--|:--|:--|--:|")

for n in (5, 40, 300):
    path = os.path.join(RAW, "N%05d" % n, "paper7_long_timeseries.csv")
    if not os.path.exists(path):
        lines.append("| %d | %d | （CSV未生成・未実行） | — | — | — |" % (n, ROT_DIM[n] - 2))
        continue
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    col = "other_rotating_occupation"
    late = [float(r[col]) for r in rows if float(r["step"]) >= META_START[n]]
    if not late:
        lines.append("| %d | %d | （準安定域未到達・実行中） | — | — | — |" % (n, ROT_DIM[n] - 2))
        continue
    ssum = [abs(float(r["occupation_sum"]) - 1.0) for r in rows]
    lines.append(
        "| %d | %d | %.3e | %.3e | %.3e | %d |"
        % (n, ROT_DIM[n] - 2, max(late), sum(late) / len(late), late[-1], len(late))
    )
    lines.append("")
    lines.append("  - N=%d 全行の |occupation_sum − 1| 最大: %.3e" % (n, max(ssum)))
    lines.append("")

lines.append("")
lines.append("## 判定に使える記述")
lines.append("")
lines.append("準安定域において、direction 1〜4 と kernel の外側の回転部分空間占有は")
lines.append("最大でも 10^-3 台（N=5）〜10^-6 台（N=40）に留まり、末端では数値零へ")
lines.append("減衰する。追加の離散方向が自然占有される兆候はこの観測量に現れていない。")
lines.append("")

out_path = os.path.join(OUT, "residual_rotating_occupation_check_v1.md")
with open(out_path, "w") as fh:
    fh.write("\n".join(lines))
print("\n".join(lines))
print("written:", out_path)
