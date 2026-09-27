#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の対数半径軌道図プログラムを、入力と保存先だけを変えて実行する（ラッパー）。

  python3 fig50_log_radial_orbits.py <pattern_id>

第六のプログラムは CSV を読む。結合系の row を第六の行の形式に並べ替え、図に使う範囲
（先頭から N_ORBITS 周ぶんの行）を CSV に書き出して渡す。値の計算や間引きはしない。
"""
from __future__ import annotations
import csv
import sys
from fig_common import PAPER6_COLUMNS, figure_dir, load_program, paper6_rows, read_unified_rows, relation_conditions


def main() -> None:
    pattern_id = sys.argv[1]
    base = figure_dir(pattern_id)
    data, columns = read_unified_rows(pattern_id)
    plot = load_program("plot_first_10_orbits_log_radial_scale.py")
    n_rows = plot.N_ORBITS * plot.STEPS_PER_ORBIT + 1
    integer_columns = {"step", "q_index"}
    for rel, cid, *_r in relation_conditions(pattern_id):
        g = paper6_rows(data, columns, rel)[:n_rows]
        work = base / cid / "log_radial_input"
        work.mkdir(parents=True, exist_ok=True)
        plot.RAW = work / "raw_macro_trajectory_first_orbits.csv"
        plot.SVG_OUT = base / cid / "figure07_first_10_orbits_log_radial_scale.svg"
        plot.PNG_OUT = base / cid / "figure07_first_10_orbits_log_radial_scale.png"
        with plot.RAW.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(PAPER6_COLUMNS)
            for row in g:
                w.writerow([int(v) if name in integer_columns else float(v) for name, v in zip(PAPER6_COLUMNS, row)])
        print(cid, flush=True)
        plot.main()


if __name__ == "__main__":
    main()
