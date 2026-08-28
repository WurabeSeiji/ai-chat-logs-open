#!/usr/bin/env python3
"""Reconstruct the four missing N=16 post-processing outputs.

Inputs are base outputs of run_N16_complex_simplex_physics.py and the preserved
decompactification geometry summary; none of the four target files is read.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

SELECTED=[0,50,100,150,200,250,300,350,400,450,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000]
TOLS=[1e-2,1e-3,1e-4,1e-6,1e-8]

def greedy_fixed_centers(vals,tol):
    centers=[]
    for v in vals:
        hit=False
        for c in centers:
            if np.max(np.abs(v-c)) < tol:
                hit=True; break
        if not hit: centers.append(v.copy())
    return len(centers)

def forever_step(steps,values,threshold):
    values=np.asarray(values,float); good=values < threshold
    suffix=np.logical_and.accumulate(good[::-1])[::-1]
    idx=np.flatnonzero(suffix)
    return int(steps.iloc[idx[0]]) if len(idx) else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--selected-snapshots',type=Path,required=True,help='N16_selected_snapshots_long.csv')
    ap.add_argument('--global-summary',type=Path,required=True,help='N16_global_summary.csv')
    ap.add_argument('--geometry-summary',type=Path,required=True,help='decompact_N16_geometry_summary.csv')
    ap.add_argument('--outdir',type=Path,required=True)
    a=ap.parse_args(); out=a.outdir; out.mkdir(parents=True,exist_ok=True)
    snap=pd.read_csv(a.selected_snapshots)
    glob=pd.read_csv(a.global_summary)
    geo=pd.read_csv(a.geometry_summary)
    final=snap[snap.step==5000].sort_values('edge_index').copy()
    final['a2n']=final.a2/final.r2; final['b2n']=final.b2/final.r2; final['abn']=final.ab/final.r2
    final['z2n_re']=final.z2_re/final.r2; final['z2n_im']=final.z2_im/final.r2
    outcols=list(snap.columns)+['a2n','b2n','abn','z2n_re','z2n_im']
    final[outcols].to_csv(out/'N16_step5000_final_edges_with_normalized_components.csv',index=False)
    # cluster counts at selected snapshots. The historical implementation used fixed first-member centers.
    crows=[]
    for st in SELECTED:
        g=snap[snap.step==st].sort_values('edge_index')
        vals=np.column_stack([g.a2/g.r2,g.b2/g.r2,g.ab/g.r2])
        for tol in TOLS: crows.append((st,tol,greedy_fixed_centers(vals,tol)))
    cdf=pd.DataFrame(crows,columns=['step','tol','triplet_cluster_count'])
    cdf.to_csv(out/'N16_triplet_cluster_counts_selected_steps.csv',index=False)
    # final phase gaps modulo pi
    th=np.mod(final.theta.to_numpy(float),np.pi)/np.pi
    s=np.sort(th); gaps=np.diff(np.r_[s,s[0]+1.0])
    def cfinal(tol):
        vals=np.column_stack([final.a2n,final.b2n,final.abn])
        return greedy_fixed_centers(vals,tol)
    ps={
      'final_phase_mod_pi_min_gap':float(gaps.min()),
      'final_phase_mod_pi_max_gap':float(gaps.max()),
      'final_phase_mod_pi_mean_gap':float(gaps.mean()),
      'final_phase_mod_pi_gap_std':float(gaps.std()),
      'final_triplet_clusters_tol_1e-2':int(cfinal(1e-2)),
      'final_triplet_clusters_tol_1e-4':int(cfinal(1e-4)),
      'final_triplet_clusters_tol_1e-6':int(cfinal(1e-6)),
    }
    (out/'N16_phase_structure_summary.json').write_text(json.dumps(ps,indent=2),encoding='utf-8')
    # time-separation milestones
    g=geo[['step','H_perp']].merge(glob[['step','r2_min','r2_max']],on='step',how='inner').sort_values('step')
    hf=float(g.loc[g.step==5000,'H_perp'].iloc[0])
    rows=[]
    for p in [.5,.9,.95,.99]:
        hit=g.loc[g.H_perp >= p*hf,'step']; rows.append((f'H_perp >= {int(p*100)}% final',int(hit.iloc[0])))
    mean=.5*(g.r2_min+g.r2_max); rel=(g.r2_max-g.r2_min)/mean
    for th in [1e-1,1e-2,1e-3,1e-4,1e-6,1e-8,1e-10]:
        rows.append((f'r2 relative spread < {th:g} forever',forever_step(g.step,rel,th)))
    pd.DataFrame(rows,columns=['metric','step']).to_csv(out/'N16_time_separation_milestones.csv',index=False)
    print('N16 post-processing reconstructed:',out)
if __name__=='__main__': main()
