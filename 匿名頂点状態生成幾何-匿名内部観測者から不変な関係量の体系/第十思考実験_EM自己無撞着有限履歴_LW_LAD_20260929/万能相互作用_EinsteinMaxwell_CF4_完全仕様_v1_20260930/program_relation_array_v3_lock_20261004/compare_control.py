#!/usr/bin/env python3
"""対照走行の照合：v2（--lock none）の csv の先頭 20 列が v1 の csv と文字列として一致するかを数える。"""
import csv, sys

a, b = sys.argv[1], sys.argv[2]
with open(a) as fa, open(b) as fb:
    ra = list(csv.reader(fa))
    rb = list(csv.reader(fb))
ncol = 20
mism = 0
rows = min(len(ra), len(rb))
for i in range(rows):
    if ra[i][:ncol] != rb[i][:ncol]:
        mism += 1
        if mism <= 3:
            print(f"row {i}: differs\n  v1: {ra[i][:ncol]}\n  v2: {rb[i][:ncol]}")
print(f"rows compared: {rows} (v1 {len(ra)}, v2 {len(rb)}), columns: {ncol}, rows differing: {mism}")
sys.exit(0 if (mism == 0 and len(ra) == len(rb)) else 1)
