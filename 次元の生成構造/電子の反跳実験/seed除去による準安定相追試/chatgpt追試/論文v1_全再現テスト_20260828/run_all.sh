#!/bin/bash
# 全工程：zip 収集 → original/rerun 展開 → パス変更 4 箇所適用 → 再実行 → 突合
set -e; cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg
mkdir -p zips original rerun results/logs
for z in ../*_20260826.zip; do n=$(basename "$z"); case "$n" in K_sigma*|N3_N4*|N5_complex_simplex_complete*|N6_N7*|N8_N9*|N10_N11*|N12_N13*|N14_N15*|N16_complex*|N3_N16_partial_zero_closure_analysis*|N3_N16_nontrivial*|N14_N16_complete*|N5_dynamics_followup*|complex_simplex_decompactification*) cp -n "$z" zips/ ;; esac; done
cp -n "../../../../ゼロ閉塞の幾何・代数構造/複製_ダンプ版_v1/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/A2a_N5_ab_probe_20260825/N5_gamma_continuum_test_bundle_20260825.zip" zips/ 2>/dev/null || true
for z in zips/*.zip; do n=$(basename "$z" .zip); [ -d "original/$n" ] || { mkdir -p "original/$n"; unzip -q -o "$z" -d "original/$n"; }; [ -d "rerun/$n" ] || cp -R "original/$n" "rerun/$n"; done
# パス変更（results/path_patches.diff と同内容）
( cd rerun && patch -p1 -N < ../results/path_patches.diff ) || true
mkdir -p rerun/N14_N16_complete_nontrivial_zero_closure_search_20260826/compat/bits; cp compat_bits_stdc++.h rerun/N14_N16_complete_nontrivial_zero_closure_search_20260826/compat/bits/stdc++.h
( cd rerun/N5_dynamics_followup_theorems_and_stability_20260826 && ln -sfn ../K_sigma_normalization_artifact_test_N4_N5_20260826 N5_sigma_normalization_artifact_test && ln -sfn ../N5_complex_simplex_complete_analysis_20260826 N5_complex_simplex_complete_analysis_20260826 )
date +%s > results/run_start_epoch.txt
./run_all_reruns.sh
python3 compare_all.py; python3 compare_physics.py > /dev/null
