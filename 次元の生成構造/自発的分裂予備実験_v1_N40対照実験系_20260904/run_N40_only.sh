#!/bin/bash
# N=40 対照実験系ラッパー（元プログラムは一切無変更。物理・初期データ変更なし）
# 引数は正本 summary から同定済みの元走行条件: 40 1e-15 --after=1500 --tol=1e-12
set -e
cd "$(dirname "$0")"
CANON='../自発的分裂予備実験_v1/largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv'

echo '===== [1/4] 対照走行（元プログラム無変更） ====='
python3 run_spontaneous_splitting_largeN_v1.py 40 1e-15 --after=1500 --tol=1e-12
diff largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv "$CANON" \
  && echo 'GATE1 PASS: fcurve は 2026-07-22 正本と全行 bit 一致'

echo '===== [2/4] 状態保存走行（差分は保存追記のみ・力学無変更の検証込み） ====='
python3 run_spontaneous_splitting_largeN_v1_savestate.py 40 1e-15 --after=1500 --tol=1e-12
diff largeN_splitting_result_v1/fcurve_N00040_delta1e-15_seed0.csv "$CANON" \
  && echo 'GATE2 PASS: 状態保存版でも fcurve は正本と全行 bit 一致（物理不変）'

echo '===== [3/4] インフレーション図（元の図化プログラム無変更・N=40のみ） ====='
python3 make_largeN_figure_v1.py

echo '===== [4/4] 追加読出し図: step0・最大step複素図・凝縮部拡大図 ====='
python3 plot_complex_plane_N40_v1.py

echo 'N40 CONTROL SYSTEM ALL DONE'
