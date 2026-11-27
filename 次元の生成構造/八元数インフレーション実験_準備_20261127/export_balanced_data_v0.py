#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
balanced_init_pm1_v3.py が作る初期データを、行の入れ替えなしで、そのままファイルに書き出す。
検算プログラム (verify_balanced_data_v0.py) は、このファイルだけを読んで検算する。
  出力: results/balanced_data/balanced_N{N:03d}.txt  (1 行目は説明、2 行目以降が元 1 個ずつ。成分は +1 / -1 を空白区切り)
        results/balanced_data/SHA256SUMS.txt          (書き出したファイルの SHA256)
対象の N: 4 から 100 までの 4 の倍数のうち、作れるもの全部。
使い方: python3 export_balanced_data_v0.py
"""
import hashlib
import os

from balanced_init_pm1_v3 import balanced_init

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "balanced_data")
os.makedirs(OUT, exist_ok=True)
sums = []
for N in range(4, 101, 4):
    try:
        T, how = balanced_init(N)                 # 行の入れ替えなし
    except ValueError as e:
        print("N = %3d: 書き出さない (%s)" % (N, e))
        continue
    name = "balanced_N%03d.txt" % N
    body = "# N = %d、元 %d 個 × 成分 %d 個、作り方: %s、行の入れ替えなし、生成: balanced_init_pm1_v3.py\n" % (N, len(T), len(T[0]), how)
    body += "".join(" ".join("+1" if v > 0 else "-1" for v in r) + "\n" for r in T)
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as f:
        f.write(body)
    sums.append((hashlib.sha256(body.encode("utf-8")).hexdigest(), name))
    print("N = %3d: %s  (元 %d 個 × 成分 %d 個)" % (N, name, len(T), len(T[0])))
with open(os.path.join(OUT, "SHA256SUMS.txt"), "w", encoding="utf-8", newline="\n") as f:
    f.write("".join("%s  %s\n" % s for s in sums))
print("%d 個のファイルを書き出した: results/balanced_data/" % len(sums))
