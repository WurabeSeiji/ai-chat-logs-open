#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
条件A/B の毎ステップ版（0-5000）。駆け上がりを時間分解する。

  python3 run_ablation_everystep_v1.py 40 A B

`run_ablation_dual_f_v1.py` の run() をそのまま使い、se_ev=1 / xmax=5000 を渡すだけ。
記録間隔を変えても、公開CSV の step にある行（25 の倍数、0..5000 の 201 行）は
既存27列が全て一致しなければならない。run() 内でそれを検査する。

第7論文の5色時系列とは違い、本系列には align_2d による記録間隔依存の列が無い。
q1..q4 は qsv4(B0, Bdom) で毎回その場の Bdom から計算され、前回記録に依存しない。
したがって**27列すべてを照合できる**。
"""
import sys
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "dual", os.path.join(HERE, "run_ablation_dual_f_v1.py"))
dual = importlib.util.module_from_spec(spec); spec.loader.exec_module(dual)

def main():
    args = sys.argv[1:]
    n = int(args[0]) if args else 40
    conds = [a.upper() for a in args[1:]] or ["A", "B"]
    print("=== 原本の照合 ===")
    mod = dual.load_original()
    print()
    for cond in conds:
        dual.run(mod, n, cond, se_ev=1, xmax=5000, tag="_everystep")

if __name__ == "__main__":
    main()
