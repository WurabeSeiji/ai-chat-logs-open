#!/bin/bash
# 復元生成器（patched）を、修正の影響を受ける 3 つの root に適用し、下流解析（followup の entropy、partial の N5）を更新する。
# usage: apply_reconstruction_to_roots.sh <ROOT> [KMODE]   ROOT = …/rerun | …/fixed | …/fixed_baseline
set -e; cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg CPLUS_INCLUDE_PATH="$PWD/compat"; [ -n "$2" ] && export KMODE="$2"
P="$PWD/extracted_patched/論文v1_全プログラム修正版_20260828"; R="$(cd "$1" && pwd)"; name=$(basename "$R"); L="$PWD/results/logs_apply_$name"; mkdir -p "$L"
N5="$R/N5_complex_simplex_complete_analysis_20260826"; N16="$R/N16_complex_simplex_complete_analysis_20260826"; P34="$R/N3_N4_complex_simplex_complete_analysis_20260826"
NON="$R/N3_N16_nontrivial_zero_closure_analysis_20260826"; PART="$R/N3_N16_partial_zero_closure_analysis_20260826"
python3 "$P/B1_N5/reconstruct_N5_missing_outputs.py" --trajectory-csv "$N5/N5_phase_by_edge_5000steps.csv" --raw-observables "$R/K_sigma_normalization_artifact_test_N4_N5_20260826/N5_raw_K_raw_observables.csv" --outdir "$N5" > "$L/B1.log" 2>&1; echo "B1 $?"
python3 "$P/copy_A_aliases.py" --root "$R" --decompact-results "$R/complex_simplex_decompactification_N5_N16_20260826/results" --n5-raw-k-source "$R/K_sigma_normalization_artifact_test_N4_N5_20260826/N5_raw_K_raw_observables.csv" > "$L/copyA.log" 2>&1; echo "copyA $?"
python3 "$P/B2_N16/reconstruct_N16_missing_outputs.py" --selected-snapshots "$N16/N16_selected_snapshots_long.csv" --global-summary "$N16/N16_global_summary.csv" --geometry-summary "$N16/decompact_N16_geometry_summary.csv" --outdir "$N16" > "$L/B2.log" 2>&1; echo "B2 $?"
python3 "$P/B3_N3_N4/reconstruct_N3_N4_phase_outputs.py" --n3-all-steps "$P34/N3_all_steps_long.csv" --n4-all-steps "$P34/N4_all_steps_long.csv" --outdir "$P34" > "$L/B3.log" 2>&1; echo "B3 $?"
for pr in "3 4" "6 7" "8 9" "10 11" "12 13" "14 15"; do set -- $pr; pd="$R/N${1}_N${2}_complex_simplex_complete_analysis_20260826"; python3 "$P/B4_comparison/build_pair_comparison_summary.py" --pair-dir "$pd" --n1 $1 --n2 $2 --out "$pd/N${1}_N${2}_comparison_summary.csv" >> "$L/B4.log" 2>&1; done; echo "B4 $?"
python3 "$P/B5_nontrivial/run_mitm_surveys.py" --source-dir "$NON" --outdir "$NON" > "$L/B5_mitm.log" 2>&1; echo "B5mitm $?"
python3 "$P/B5_nontrivial/rebuild_nontrivial_zero_closure_analysis.py" --source-dir "$NON" --mitm-small "$NON/MITM_N6_N13_k5_k6.csv" --mitm-targeted "$NON/targeted_MITM_N14_N16.csv" --n5-all-steps "$N5/N5_all_steps_a_b_a2_b2_ab.csv" --n14-all-steps "$R/N14_N15_complex_simplex_complete_analysis_20260826/N14_all_steps_long.csv" --outdir "$NON" > "$L/B5_rebuild.log" 2>&1; echo "B5rebuild $?"
# 下流：partial（N5 は N5_all_steps を読む）と followup analyze（entropy は N5_all_steps を読む）
( cd "$PART" && python3 analyze_partial_zero_closures_N3_N16.py > "$L/partial.log" 2>&1 ); echo "partial $?"
python3 "$P/B7_partial/write_zero_subset_alias.py" --source "$PART/out/N3_N16_zero_triple_exact_covers.csv" --out "$PART/out/N3_N16_zero_subset_exact_covers.csv" > "$L/B7.log" 2>&1; echo "B7 $?"
( cd "$R/N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826" && python3 analyze_followup.py > "$L/followup_analyze.log" 2>&1 ); echo "followup $?"
echo "DONE $name"
