#!/usr/bin/env python3
"""Orchestrate reconstructed missing-output generators against an extracted package root.

By default B6 stochastic k>=7 search is NOT rerun, because the historical
steps/seeds were not recovered. Pass --run-new-stochastic-search to use the
bundled new deterministic reproduction matrix (explicitly not historical).
"""
from __future__ import annotations
import argparse,subprocess,sys,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent

def run(cmd): print('+',' '.join(map(str,cmd))); subprocess.run(list(map(str,cmd)),check=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--decompact-results',type=Path,default=None);ap.add_argument('--n5-raw-k-source',type=Path,default=None);ap.add_argument('--run-new-stochastic-search',action='store_true');a=ap.parse_args();r=a.root
 # Copy-only inputs first.
 cmd=[sys.executable,HERE/'copy_A_aliases.py','--root',r]
 if a.decompact_results:cmd += ['--decompact-results',a.decompact_results]
 if a.n5_raw_k_source:cmd += ['--n5-raw-k-source',a.n5_raw_k_source]
 run(cmd)
 eng=r/'N5_complex_simplex_complete_analysis_20260826'/'run_n_scaling_lowrank_v1_no_sigma_norm.py'; n5=r/'N5_complex_simplex_complete_analysis_20260826'
 run([sys.executable,HERE/'B1_N5/reconstruct_N5_missing_outputs.py','--engine',eng,'--raw-observables',n5/'N5_raw_K_raw_observables.csv','--outdir',n5])
 n16=r/'N16_complex_simplex_complete_analysis_20260826'
 run([sys.executable,HERE/'B2_N16/reconstruct_N16_missing_outputs.py','--selected-snapshots',n16/'N16_selected_snapshots_long.csv','--global-summary',n16/'N16_global_summary.csv','--geometry-summary',n16/'decompact_N16_geometry_summary.csv','--outdir',n16])
 p34=r/'N3_N4_complex_simplex_complete_analysis_20260826';run([sys.executable,HERE/'B3_N3_N4/reconstruct_N3_N4_phase_outputs.py','--n3-all-steps',p34/'N3_all_steps_long.csv','--n4-all-steps',p34/'N4_all_steps_long.csv','--outdir',p34])
 pairs=[(3,4),(6,7),(8,9),(10,11),(12,13),(14,15)]
 for n1,n2 in pairs:
  pd=r/f'N{n1}_N{n2}_complex_simplex_complete_analysis_20260826';run([sys.executable,HERE/'B4_comparison/build_pair_comparison_summary.py','--pair-dir',pd,'--n1',n1,'--n2',n2,'--out',pd/f'N{n1}_N{n2}_comparison_summary.csv'])
 non=r/'N3_N16_nontrivial_zero_closure_analysis_20260826';run([sys.executable,HERE/'B5_nontrivial/run_mitm_surveys.py','--source-dir',non,'--outdir',non]);run([sys.executable,HERE/'B5_nontrivial/rebuild_nontrivial_zero_closure_analysis.py','--source-dir',non,'--mitm-small',non/'MITM_N6_N13_k5_k6.csv','--mitm-targeted',non/'targeted_MITM_N14_N16.csv','--n5-all-steps',n5/'N5_all_steps_a_b_a2_b2_ab.csv','--n14-all-steps',r/'N14_N15_complex_simplex_complete_analysis_20260826'/'N14_all_steps_long.csv','--outdir',non])
 partial=r/'N3_N16_partial_zero_closure_analysis_20260826';run([sys.executable,HERE/'B7_partial/write_zero_subset_alias.py','--source',partial/'N3_N16_zero_triple_exact_covers.csv','--out',partial/'N3_N16_zero_subset_exact_covers.csv'])
 if a.run_new_stochastic_search:
  srch=r/'N14_N16_complete_nontrivial_zero_closure_search_20260826';work=srch/'method_results_reproduction';run([sys.executable,HERE/'B6_all_cardinalities/run_search_methods.py','--source-dir',srch,'--anneal-args',HERE/'B6_all_cardinalities/reproduction_anneal_args.csv','--workdir',work]);run([sys.executable,HERE/'B6_all_cardinalities/aggregate_search_results.py','--workdir',work,'--outdir',srch]);run([sys.executable,HERE/'B6_all_cardinalities/build_time_evolution_and_n5_covers.py','--search-results',srch/'N14_N16_all_cardinalities_best_results.csv','--n14-all-steps',r/'N14_N15_complex_simplex_complete_analysis_20260826'/'N14_all_steps_long.csv','--n5-all-steps',n5/'N5_all_steps_a_b_a2_b2_ab.csv','--n5-final',non/'SOURCE_N5_step5000_final_edges.csv','--outdir',srch])
  # install shared B6 products into B5
  for f in ['N14_best_k6_candidate_time_evolution.csv','N5_best_nontrivial_pair_time_evolution.csv','N5_nontrivial_pair_exact_covers.csv']:
   shutil.copy2(srch/f,non/f)
 else:
  print('B6 stochastic search skipped: historical steps/seeds were not recovered. Use --run-new-stochastic-search for the new explicit reproduction matrix.')
if __name__=='__main__':main()
