#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図化サンプル4枚のための対照走行（条件固定ランナー）v1

**このスクリプトは走行条件を固定して記録するだけのラッパーである。**
正本 `run_nsweep_three_series_v2.py` は一行も書き換えず、
コマンドライン引数を与えて起動するだけ。物理・記録項目・判定・図の
生成ロジックには一切触れない。

----------------------------------------------------------------------
条件と、その選択根拠（既存データからの実測）
----------------------------------------------------------------------
  mode  = electron   （単一種。図の解釈が一意になる）
  N     = 12         （長時間走行が存在する唯一の N。N=16 の長時間データは無い）
  delta = 0.1        （シード強度）
  T     = 4000       （後述の通り、この条件では 4000 步で完結する）

**なぜ delta = 0.1 か**
既存 `nsweep_electron_T42000_d*_N12_v2.npz` 8 水準の先頭 4000 步を実測した:

    delta      f_seed@0     f_seed@4000    倍率     T=42000 最終値への到達度
    1e-15      1.004e-30    1.696e-30      1.69     59.3 %
    1e-08      1.000e-16    1.000e-16      1.00    100.0 %（何も起きない）
    1e-04      1.000e-08    1.000e-08      1.00    100.0 %（何も起きない）
    1e-03      1.000e-06    1.000e-06      1.00    100.0 %（何も起きない）
    1e-02      9.999e-05    1.000e-04      1.001    94.1 %（物質が動かない）
    3.16e-02   9.990e-04    1.053e-03      1.054     0.2 %（立ち上がりが遅すぎる）
    4.36e-02   1.895e-03    2.308e-03      1.218     0.5 %（同上）
    1.00e-01   9.901e-03    5.054e-01     51.05    102.0 %  ← 採用

  delta <= 1e-2 では物質がまったく育たない。
  delta = 3.16e-2 は T=42000 では 504 倍まで育つが、T=4000 では 0.2 % しか
  進まない（閾値を超えていても立ち上がりが遅い）。
  delta = 0.1 のみが T=4000 で最終値に到達する。T を伸ばす必要がない。

**この走行で得たい 4 時点（図化サンプル）**
  1. 開始                          tau = 0
  2. インフレーション開始           tau = 立ち上がりの開始（本走行で確定する）
  3. インフレーション終了＝準安定開始 tau = 飽和点（同上）
  4. 終了                          tau = 4000

**本走行では新しい記録項目を追加しない。**
まず正本のままで対照データを取り、既存データと同じものが得られることを
確認する。データ抽出の追加はその後の別ステップとする。

----------------------------------------------------------------------
出力（正本の命名規則にそのまま従う。TAG = "_d0.1"）
----------------------------------------------------------------------
  result_nsweep_electron_d0.1_v2.json
  nsweep_electron_d0.1_N12_v2.npz
  fig_electron_d0.1_4panel_N12_v2.png
  fig_electron_d0.1_mix_N12_v2.png
  fig_electron_d0.1_ledger_N12_v2.png
  fig_electron_d0.1_summary_v2.png
  fig_electron_d0.1_birth_matrix_v2.png
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGET = HERE / "run_nsweep_three_series_v2.py"

MODE = "electron"
NMIN = "12"
NMAX = "12"
T = "4000"
DELTA = "0.1"

ARGV = [sys.executable, str(TARGET), MODE, NMIN, NMAX, T, DELTA]

if __name__ == "__main__":
    assert TARGET.exists(), f"正本が見つからない: {TARGET}"
    print("=" * 78)
    print("図化サンプル4枚のための対照走行")
    print(f"  正本      : {TARGET.name}（無改変）")
    print(f"  mode      : {MODE}")
    print(f"  N         : {NMIN}")
    print(f"  T         : {T}")
    print(f"  delta     : {DELTA}")
    print(f"  実行       : {' '.join(ARGV[1:])}")
    print("=" * 78)
    sys.exit(subprocess.call(ARGV, cwd=str(HERE)))
