#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の対数半径軌道図プログラムを、入力と保存先だけを変えて実行する。

第六のプログラムは CSV を読む。今回の生成器の出力は HDF5 なので、図に使う範囲
（先頭から N_ORBITS 周ぶんの行）を、列名も値もそのまま CSV に書き出してから渡す。
値の計算や間引きはしない。
"""
from __future__ import annotations
import csv
import json
import h5py
from wrapper_common import HERE, SELF_CASES, load_program


def export_rows(case_dir, n_rows, csv_path):
    manifest = json.loads((case_dir / "part_manifest.json").read_text(encoding="utf-8"))
    first = case_dir / manifest["parts"][0]["file"]
    with h5py.File(first, "r") as hf:
        ds = hf["raw_macro_trajectory"]
        columns = json.loads(hf.attrs["columns_json"])
        if ds.shape[0] < n_rows:
            raise RuntimeError(f"{first.name}: rows {ds.shape[0]} < {n_rows}")
        block = ds[:n_rows]
    integer_columns = {"step", "q_index"}
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for row in block:
            w.writerow([int(v) if name in integer_columns else float(v) for name, v in zip(columns, row)])


def main() -> None:
    plot = load_program("plot_first_10_orbits_log_radial_scale.py")
    n_rows = plot.N_ORBITS * plot.STEPS_PER_ORBIT + 1
    for cid, *_rest in SELF_CASES:
        case_dir = HERE / "cases" / cid
        work = case_dir / "log_radial_input"
        work.mkdir(exist_ok=True)
        plot.RAW = work / "raw_macro_trajectory_first_orbits.csv"
        plot.SVG_OUT = case_dir / "figures" / "figure07_first_10_orbits_log_radial_scale.svg"
        plot.PNG_OUT = case_dir / "figures" / "figure07_first_10_orbits_log_radial_scale.png"
        export_rows(case_dir, n_rows, plot.RAW)
        print(cid, flush=True)
        plot.main()


if __name__ == "__main__":
    main()
