#!/usr/bin/env python3
"""Install the copy/alias-only files identified as group A in the audit.

The script intentionally performs no physics calculation.  It copies canonical
base outputs into the analysis packages that consume them. N=5 is special:
SOURCE_N5_step5000_final_edges.csv is derived by selecting step 5000 from the
canonical all-steps table and appending i,j in K5 edge order.
"""
from __future__ import annotations
import argparse,shutil
from pathlib import Path
import numpy as np,pandas as pd

def pkg_for_n(root,n):
 if n in (3,4): name='N3_N4_complex_simplex_complete_analysis_20260826'
 elif n==5: name='N5_complex_simplex_complete_analysis_20260826'
 elif n in (6,7): name='N6_N7_complex_simplex_complete_analysis_20260826'
 elif n in (8,9): name='N8_N9_complex_simplex_complete_analysis_20260826'
 elif n in (10,11): name='N10_N11_complex_simplex_complete_analysis_20260826'
 elif n in (12,13): name='N12_N13_complex_simplex_complete_analysis_20260826'
 elif n in (14,15): name='N14_N15_complex_simplex_complete_analysis_20260826'
 elif n==16: name='N16_complex_simplex_complete_analysis_20260826'
 return root/name

def source_edges(root,n,tmp):
 p=pkg_for_n(root,n)
 f=p/f'N{n}_step5000_final_edges.csv'
 if f.exists(): return f
 if n==5:
  d=pd.read_csv(p/'N5_all_steps_a_b_a2_b2_ab.csv'); d=d[d.step==5000].sort_values('edge_index').copy(); ea,eb=np.triu_indices(5,k=1); d['i']=ea+1;d['j']=eb+1;out=tmp/'SOURCE_N5_step5000_final_edges.csv';d.to_csv(out,index=False);return out
 raise FileNotFoundError(f)

def cp(src,dst): dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);print(src,'->',dst)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True,help='directory containing the original analysis package folders');ap.add_argument('--decompact-results',type=Path,default=None);ap.add_argument('--n5-raw-k-source',type=Path,default=None);a=ap.parse_args();tmp=a.root/'.copy_alias_tmp';tmp.mkdir(exist_ok=True)
 partial=a.root/'N3_N16_partial_zero_closure_analysis_20260826';nontr=a.root/'N3_N16_nontrivial_zero_closure_analysis_20260826';search=a.root/'N14_N16_complete_nontrivial_zero_closure_search_20260826'
 for n in range(3,17):
  src=source_edges(a.root,n,tmp);cp(src,partial/f'SOURCE_N{n}_step5000_final_edges.csv');cp(src,nontr/f'SOURCE_N{n}_step5000_final_edges.csv')
  if n>=14:cp(src,search/f'SOURCE_N{n}_step5000_final_edges.csv')
 if a.decompact_results:
  n16=pkg_for_n(a.root,16)
  for stem in ['geometry_summary','perp_axis_growth_rates','takagi_axes']:
   src=a.decompact_results/f'N16_{stem}.csv'; cp(src,n16/f'decompact_N16_{stem}.csv')
 if a.n5_raw_k_source: cp(a.n5_raw_k_source,pkg_for_n(a.root,5)/'N5_raw_K_raw_observables.csv')
 # B6 is canonical for these shared outputs once it has run.
 shared=['N14_best_k6_candidate_time_evolution.csv','N5_best_nontrivial_pair_time_evolution.csv','N5_nontrivial_pair_exact_covers.csv']
 for f in shared:
  src=search/f
  if src.exists(): cp(src,nontr/f)
 try: tmp.rmdir()
 except OSError: pass
if __name__=='__main__':main()
