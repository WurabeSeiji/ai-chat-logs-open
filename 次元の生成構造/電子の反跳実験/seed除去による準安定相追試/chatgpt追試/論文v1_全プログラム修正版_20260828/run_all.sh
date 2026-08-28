#!/bin/bash
# 全工程：original（../論文v1_全再現テスト_20260828/original）から fixed/ を作り、修正を適用・検証、再実行、突合、比較図、報告 md
set -e; cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg
rm -rf fixed && mkdir fixed && cp -R ../論文v1_全再現テスト_20260828/original/* fixed/ && rm -rf fixed/*/__pycache__ fixed/*/*/__pycache__
( cd fixed && patch -p1 -N -s < ../../論文v1_全再現テスト_20260828/results/path_patches.diff )
python3 apply_fixes.py
python3 - <<'PY'
p="fixed/complex_simplex_decompactification_N5_N16_20260826/run_complex_simplex_decompactification.py"; s=open(p,encoding="utf-8").read()
old="f\"- R_perp early log growth rate: {s['R_perp_log_growth_rate_per_step']:.6f}/step\","
s=s.replace(old,"(f\"- R_perp early log growth rate: {s['R_perp_log_growth_rate_per_step']:.6f}/step\" if s['R_perp_log_growth_rate_per_step'] is not None else \"- R_perp early log growth rate: (no exponential regime found; fit window empty)\"),  # ROBUSTNESS PATCH"); open(p,"w",encoding="utf-8").write(s)
PY
mkdir -p fixed/N14_N16_complete_nontrivial_zero_closure_search_20260826/compat/bits && cp ../論文v1_全再現テスト_20260828/compat_bits_stdc++.h fixed/N14_N16_complete_nontrivial_zero_closure_search_20260826/compat/bits/stdc++.h
( cd fixed/N5_dynamics_followup_theorems_and_stability_20260826 && ln -sfn ../K_sigma_normalization_artifact_test_N4_N5_20260826 N5_sigma_normalization_artifact_test && ln -sfn ../N5_complex_simplex_complete_analysis_20260826 N5_complex_simplex_complete_analysis_20260826 )
python3 verify_fixes_all.py | tee results/verify_fixes_all.log
./run_all_fixed.sh
( cd fixed/N3_N16_partial_zero_closure_analysis_20260826 && for N in 3 4 6 7 8 9 10 11 12 13 14 15 16; do cp ../N3_N16_nontrivial_zero_closure_analysis_20260826/SOURCE_N${N}_step5000_final_edges.csv .; done )
( cd fixed/N14_N16_complete_nontrivial_zero_closure_search_20260826 && for N in 14 15 16; do python3 mitm56.py $N rerun_mitm56_N$N.csv; done )
python3 compare_all_fixed.py; python3 compare_physics_fixed.py > /dev/null; python3 make_comparison_figures.py; python3 write_reports.py
