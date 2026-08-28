#!/usr/bin/env python3
"""Generate N*_N*_comparison_summary.csv for the six paired packages."""
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd

def wrap_pi(x): return ((x+np.pi/2)%np.pi)-np.pi/2

def phase_dist(a,b): return abs(wrap_pi(a-b))/np.pi

def tol_count(path:Path,tol=1e-2):
    d=pd.read_csv(path)
    q=d[(d.step==5000)&np.isclose(d.tol,tol)]
    if q.empty: raise RuntimeError(f'no step=5000 tol={tol} in {path}')
    return int(q.iloc[0].triplet_cluster_count), str(q.iloc[0].get('cluster_sizes_desc',''))

def summary(path): return json.loads(path.read_text())

def build_34(d:Path,out:Path):
    rows=[]
    for n in (3,4):
        df=pd.read_csv(d/f'N{n}_all_steps_long.csv'); f=df[df.step==5000].sort_values('edge_index'); th=f.theta.to_numpy(float)
        if n==3:
            dd=np.array([phase_dist(th[i],th[j]) for i in range(3) for j in range(i+1,3)]); err=float(np.max(np.abs(dd-1/3))); sizes='1+1+1'
        else:
            pairs=[(0,5),(1,4),(2,3)]
            def cm(x): return float((np.angle(np.mean(np.exp(2j*x)))/2)%np.pi)/np.pi
            ctr=np.sort(np.array([cm(th[[a,b]]) for a,b in pairs])); gaps=np.diff(np.r_[ctr,ctr[0]+1]); err=float(np.max(np.abs(gaps-1/3))); sizes='2+2+2'
        rows.append((f'N={n}',n,n*(n-1)//2,n-1,float(f.r2.min()),float(f.r2.max()),3,sizes,err))
    pd.DataFrame(rows,columns=['case','N','M','simplex_rank','final_r2_min','final_r2_max','final_phase_distance_classes_tol1e-8','class_sizes','phase_thirds_error']).to_csv(out,index=False)

def build_generic(d:Path,n1:int,n2:int,out:Path):
    vals=[]
    for n in (n1,n2):
        s=summary(d/f'N{n}_summary.json'); c2,sizes=tol_count(d/f'N{n}_triplet_cluster_counts.csv',1e-2)
        vals.append([f'N={n}',n,int(s['M']),int(s['final_rank']),float(s['final_H_perp']),float(s['final_r2_min']),float(s['final_r2_max']),int(s['final_class_count_tol1e-6']),c2,sizes])
    if (n1,n2)==(6,7):
        coarse=['2x6 + 1x3 at tol1e-2','all distinct at tol1e-2']
        rows=[v[:-1]+[coarse[i]] for i,v in enumerate(vals)]
        cols=['case','N','M','rank','final_H_perp','r2_min','r2_max','classes_tol1e-6','classes_tol1e-2','coarse_pattern']
    elif (n1,n2)==(8,9):
        coarse=['3+2+2+2+2 + 17 singles at tol1e-2','2+2 + 32 singles at tol1e-2']
        rows=[v[:-1]+[coarse[i]] for i,v in enumerate(vals)]
        cols=['case','N','M','rank','final_H_perp','r2_min','r2_max','classes_tol1e-6','classes_tol1e-2','coarse_pattern']
    else:
        rows=[v[:-1] for v in vals]
        cols=['case','N','M','rank','final_H_perp','r2_min','r2_max','classes_tol1e-6','classes_tol1e-2']
    pd.DataFrame(rows,columns=cols).to_csv(out,index=False)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--pair-dir',type=Path,required=True); ap.add_argument('--n1',type=int,required=True); ap.add_argument('--n2',type=int,required=True); ap.add_argument('--out',type=Path,default=None); a=ap.parse_args()
    out=a.out or a.pair_dir/f'N{a.n1}_N{a.n2}_comparison_summary.csv'; out.parent.mkdir(parents=True,exist_ok=True)
    if (a.n1,a.n2)==(3,4): build_34(a.pair_dir,out)
    else: build_generic(a.pair_dir,a.n1,a.n2,out)
    print(out)
if __name__=='__main__': main()
