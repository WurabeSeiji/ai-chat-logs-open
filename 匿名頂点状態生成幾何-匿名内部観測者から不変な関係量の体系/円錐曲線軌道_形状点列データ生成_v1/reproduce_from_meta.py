#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
保存した生成パラメータだけから点列を再生成し、入力ファイルと一致することを確かめる（数値実験01）。

    P_n = O + R(φ) · p / (1 + e cos θ_n) · (cos θ_n, sin θ_n)

p, e, O, φ は data/truth/orbit_Dxx_meta.json から、θ_n は data/truth/orbit_Dxx_truth.csv の
theta_rad 列から読む。結果は data/truth/reproduction_check.json に保存する。

使い方
    python reproduce_from_meta.py
    python reproduce_from_meta.py --data データフォルダ
"""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f))
    header, body = rows[0], rows[1:]
    return {h: np.array([float(r[i]) for r in body]) for i, h in enumerate(header)}


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="保存パラメータから点列を再生成して入力と照合する")
    parser.add_argument("--data", default=None, help="データフォルダ（既定: このスクリプトと同じフォルダの data）")
    args = parser.parse_args()
    data = Path(args.data) if args.data else Path(__file__).resolve().parent / "data"

    manifest = load_json(data / "truth" / "manifest.json")
    results = {}
    worst = 0.0
    for d in manifest["datasets"]:
        meta = load_json(data / d["meta_file"])
        inp = read_csv(data / d["input_file"])
        tru = read_csv(data / d["truth_file"])

        p = meta["shape_parameters"]["p"]
        e = meta["shape_parameters"]["e"]
        O = meta["pose"]["focus_offset_O_xy"]
        phi = meta["pose"]["rotation_phi_rad"]
        th = tru["theta_rad"]

        r = p / (1.0 + e * np.cos(th))
        xc, yc = r * np.cos(th), r * np.sin(th)
        c, s = math.cos(phi), math.sin(phi)
        x = O[0] + c * xc - s * yc
        y = O[1] + s * xc + c * yc

        scale = max(1.0, float(np.max(np.abs(np.concatenate((inp["x"], inp["y"]))))))
        dev = float(max(np.max(np.abs(x - inp["x"])), np.max(np.abs(y - inp["y"]))) / scale)
        rows_match = bool(np.array_equal(inp["n"], tru["n"]))
        worst = max(worst, dev)
        results[d["dataset_id"]] = {"n_points": int(th.size),
                                    "max_deviation_relative_to_coordinate_scale": dev,
                                    "input_and_truth_rows_aligned": rows_match}
        print(f'{d["dataset_id"]}  N={th.size:4d}  max relative deviation = {dev:.3e}  rows aligned = {rows_match}')

    passed = worst < 1e-12 and all(v["input_and_truth_rows_aligned"] for v in results.values())
    out = {"formula": "P_n = O + R(phi) * p / (1 + e cos theta_n) * (cos theta_n, sin theta_n)",
           "max_deviation_over_all_datasets": worst,
           "passed": passed,
           "datasets": results}
    with open(data / "truth" / "reproduction_check.json", "w", encoding="utf-8", newline="\n") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print(f"max relative deviation over all datasets = {worst:.3e}  passed = {passed}")


if __name__ == "__main__":
    main()
