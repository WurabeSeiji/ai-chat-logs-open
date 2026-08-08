#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""本来の実験系列 N=5,40,300：指数増幅から準安定域までを観察する。

O6原本対照（run_original_O6_control_v1.py）とは目的が異なる。
  - O6原本対照：論文の厳密再現。N=40,300,1000、交差後1500ステップ。
  - 本ラッパー：準安定状態の観察。N=5,40,300、交差後を大きく取り
    （既定 after=20000, 上限 cap=30000）指数増幅→飽和→準安定プラトーまで追う。

原本三ファイル（run_n_scaling_lowrank_v1.py /
run_spontaneous_splitting_largeN_v1.py / make_largeN_figure_v1.py）は
一切変更しない。SHA-256 を検査したうえで import し、走行器の出力先だけを
metastable_series_result_v1/ に切り替える。正本 largeN_splitting_result_v1/
と O6_原本対照結果_v1/ は破壊しない。

使い方:
    python3 run_metastable_series_v1.py run          # N=5,40,300 を順に走行
    python3 run_metastable_series_v1.py run 5        # 個別
    python3 run_metastable_series_v1.py plot         # 集約図
    python3 run_metastable_series_v1.py all          # 走行して作図
"""

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE.parents[2] / "次元の生成構造" / "自発的分裂予備実験_v1"
RESULT_DIR = HERE / "metastable_series_result_v1"

CORE_HASHES = {
    "run_n_scaling_lowrank_v1.py":
        "ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d",
    "run_spontaneous_splitting_largeN_v1.py":
        "13baf6f5158c53ee92d1a08a0a5b60832424a83222b0601563a2676a392515ac",
    "make_largeN_figure_v1.py":
        "f4cbc4efe2af2ffbffa45eec8c71ecb0340a1a3c940b78a2d83847fb03ea5804",
}

N_VALUES = (5, 40, 300)
DELTA = 1e-15
SEED = 0
AFTER = 20000   # 交差後にプラトーを観察するための追跡ステップ数
CAP = 30000     # 上限（交差時刻 + AFTER がこれを超えたら打ち切り）
TOL = 1e-12


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_core_hashes():
    failed = []
    for filename, expected in CORE_HASHES.items():
        local_hash = sha256(HERE / filename)
        source_hash = sha256(SOURCE_DIR / filename)
        ok = local_hash == expected and source_hash == expected
        print(f"[{'OK' if ok else 'NG'}] {filename} {local_hash}")
        if not ok:
            failed.append(filename)
    if failed:
        raise RuntimeError("原本ハッシュ不一致: " + ", ".join(failed))


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tag_for(n):
    return f"N{n:05d}_delta{DELTA:.0e}_seed{SEED}"


def run_one(n):
    if n not in N_VALUES:
        raise ValueError(f"本系列のNは {N_VALUES} のいずれかです: {n}")
    verify_core_hashes()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(HERE))
    runner = load_module(
        "metastable_runner", HERE / "run_spontaneous_splitting_largeN_v1.py"
    )
    runner.RESULT_DIR = str(RESULT_DIR)
    print(
        f"=== 準安定系列 N={n}, delta={DELTA}, seed={SEED}, "
        f"cap={CAP}, after={AFTER}, tol={TOL} ===",
        flush=True,
    )
    summary = runner.run(n, DELTA, SEED, CAP, AFTER, tol=TOL)
    # 準安定域の要約（後期窓 f の中央値・振れ幅）を追記して別ファイルに保存
    tag = tag_for(n)
    taus, fs = [], []
    with open(RESULT_DIR / f"fcurve_{tag}.csv") as fh:
        import csv
        for row in csv.DictReader(fh):
            taus.append(int(row["tau"]))
            fs.append(float(row["f"]))
    fs = np.array(fs)
    cross = summary.get("crossing_tau")
    if cross is not None and len(fs) > cross + 100:
        tail = fs[cross + 100:]
        plateau = {
            "tail_start_tau": cross + 100,
            "tail_median_f": float(np.median(tail)),
            "tail_min_f": float(np.min(tail)),
            "tail_max_f": float(np.max(tail)),
            "tail_q05_f": float(np.quantile(tail, 0.05)),
            "tail_q95_f": float(np.quantile(tail, 0.95)),
        }
    else:
        plateau = {"note": "交差後の準安定窓が不足（プラトー未到達）"}
    meta = {"experiment": "metastable_series", "after": AFTER, "cap": CAP,
            "plateau": plateau, "base_summary": summary}
    with open(RESULT_DIR / f"metastable_{tag}.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)
    print(f"[準安定要約] N={n}: {json.dumps(plateau, ensure_ascii=False)}")
    return meta


def make_plot():
    plotter = load_module(
        "metastable_plotter", HERE / "make_metastable_series_figure_v1.py"
    )
    plotter.RESULT_DIR = str(RESULT_DIR)
    plotter.N_VALUES = N_VALUES
    plotter.DELTA = DELTA
    plotter.SEED = SEED
    plotter.main()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "plot", "all"))
    parser.add_argument("n", nargs="?", type=int)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "run":
        if args.n is not None:
            run_one(args.n)
        else:
            for n in N_VALUES:
                run_one(n)
    elif args.command == "plot":
        make_plot()
    elif args.command == "all":
        for n in N_VALUES:
            run_one(n)
        make_plot()


if __name__ == "__main__":
    main()
