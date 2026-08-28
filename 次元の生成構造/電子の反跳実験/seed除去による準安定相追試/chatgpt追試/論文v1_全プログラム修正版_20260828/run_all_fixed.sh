#!/bin/bash
# 修正版の全パッケージ再実行。usage: run_all_fixed.sh <ROOT(fixed|fixed_baseline)> <KMODE(amplitude|phase)>
# 物理走行（独立）は並列 5、その後 依存順に：decompact→N16 解析、followup、SOURCE コピー→partial→nontrivial→closure search
cd "$(dirname "$0")"; ROOT="$1"; export KMODE="$2" LC_ALL=en_US.UTF-8 MPLBACKEND=Agg; L="results/logs/$ROOT"; mkdir -p "$L"; : > "$L/_status.txt"
run_job(){ d="${1%%|*}"; c="${1#*|}"; n="$(echo "$d" | tr '/' '_')"; s=$(date +%s); ( cd "$ROOT/$d" && eval "$c" ) > "$L/$n.log" 2>&1; echo "$? $(( $(date +%s)-s ))s $d" >> "$L/_status.txt"; }
export -f run_job; export ROOT L
cat > "$L/_jobs.txt" <<'J'
K_sigma_normalization_artifact_test_N4_N5_20260826|python3 run_artifact_comparison_N4_N5.py
N3_N4_complex_simplex_complete_analysis_20260826|python3 run_N3_N4_complete_analysis.py
N5_complex_simplex_complete_analysis_20260826|python3 run_N5_physical_phase_step_test.py && python3 analyze_N5_complex_simplex_structure.py && python3 plot_N5_inflation_vs_ordering.py
N6_N7_complex_simplex_complete_analysis_20260826|python3 run_N6_N7_complete_analysis.py
N8_N9_complex_simplex_complete_analysis_20260826|python3 run_N8_N9_complete_analysis.py
N10_N11_complex_simplex_complete_analysis_20260826|python3 run_N10_N11_complete_analysis.py
N12_N13_complex_simplex_complete_analysis_20260826|python3 run_N12_N13_complete_analysis.py
N14_N15_complex_simplex_complete_analysis_20260826|python3 run_N14_N15_complete_analysis.py
N16_complex_simplex_complete_analysis_20260826|python3 run_N16_complex_simplex_physics.py
complex_simplex_decompactification_N5_N16_20260826|python3 run_complex_simplex_decompactification.py
N5_gamma_continuum_test_bundle_20260825|python3 run_N5_gamma_continuum_test.py --outdir .
J
xargs -P 5 -I{} bash -c 'run_job "$@"' _ {} < "$L/_jobs.txt"
# 依存解析
run_job "N16_complex_simplex_complete_analysis_20260826|for f in geometry_summary perp_axis_growth_rates takagi_axes; do cp ../complex_simplex_decompactification_N5_N16_20260826/results/N16_\$f.csv decompact_N16_\$f.csv; done && python3 analyze_and_plot_N16.py"
( cd "$ROOT/N5_dynamics_followup_theorems_and_stability_20260826" && ln -sfn ../K_sigma_normalization_artifact_test_N4_N5_20260826 N5_sigma_normalization_artifact_test && ln -sfn ../N5_complex_simplex_complete_analysis_20260826 N5_complex_simplex_complete_analysis_20260826 )
run_job "N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826|python3 run_followup_experiments.py && cp N5_moduli_seed_sweep.csv N5_moduli_seed_sweep_from_run_followup_20seeds.csv && python3 run_moduli_sweep_fast.py && python3 analyze_followup.py"
# SOURCE コピー（A 群：修正後の最終状態で差し替える）
( cd "$ROOT" && for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do case $N in 3|4) P=N3_N4;; 5) P=N5;; 6|7) P=N6_N7;; 8|9) P=N8_N9;; 10|11) P=N10_N11;; 12|13) P=N12_N13;; 14|15) P=N14_N15;; 16) P=N16;; esac; src="${P}_complex_simplex_complete_analysis_20260826/N${N}_step5000_final_edges.csv"; for dst in N3_N16_partial_zero_closure_analysis_20260826 N3_N16_nontrivial_zero_closure_analysis_20260826; do cp "$src" "$dst/SOURCE_N${N}_step5000_final_edges.csv"; done; [ $N -ge 14 ] && cp "$src" "N14_N16_complete_nontrivial_zero_closure_search_20260826/SOURCE_N${N}_step5000_final_edges.csv"; done; cp K_sigma_normalization_artifact_test_N4_N5_20260826/N5_raw_K_raw_observables.csv N5_complex_simplex_complete_analysis_20260826/N5_raw_K_raw_observables.csv )
run_job "N3_N16_partial_zero_closure_analysis_20260826|python3 analyze_partial_zero_closures_N3_N16.py"
run_job "N3_N16_nontrivial_zero_closure_analysis_20260826|python3 analyze_nontrivial_zero_closures.py"
mkdir -p "$ROOT/N14_N16_complete_nontrivial_zero_closure_search_20260826/compat/bits"; cp compat_bits_stdc++.h "$ROOT/N14_N16_complete_nontrivial_zero_closure_search_20260826/compat/bits/stdc++.h"
run_job "N14_N16_complete_nontrivial_zero_closure_search_20260826|c++ -O2 -std=c++17 -I compat -o exact_k234 exact_k234.cpp && c++ -O2 -std=c++17 -I compat -o search_subsets search_subsets.cpp && for N in 14 15 16; do ./exact_k234 \$N SOURCE_N\${N}_step5000_final_edges.csv rerun_exact_k234_N\$N.csv; ./search_subsets \$N SOURCE_N\${N}_step5000_final_edges.csv rerun_search_subsets_N\$N.csv; python3 mitm56.py \$N rerun_mitm56_N\$N.csv; done"
echo DONE >> "$L/_status.txt"
