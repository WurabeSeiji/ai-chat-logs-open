#!/usr/bin/env python3
"""Aggregate exact/MITM/annealing method outputs into the two published tables + plots."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workdir',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True);a=ap.parse_args();a.outdir.mkdir(parents=True,exist_ok=True);rows=[]
 for N in (14,15,16):
  e=pd.read_csv(a.workdir/f'exact_N{N}_k234.csv')
  for _,r in e.iterrows():rows.append([N,int(r.M),int(r.k),float(r.best_residual),str(r.edges).rstrip(';').replace(';',','),'EXACT exhaustive',True,float(r.best_residual)<1e-6])
 m=pd.read_csv(a.workdir/'mitm_N14_N16_k56.csv')
 for _,r in m.iterrows():rows.append([int(r.N),int(r.M),int(r.k),float(r.best_residual),str(r.edges),'MITM nearest-neighbor',False,float(r.best_residual)<1e-6])
 an=pd.read_csv(a.workdir/'anneal_results.csv')
 for _,r in an.iterrows():rows.append([int(r.N),int(r.M),int(r.k),float(r.best_residual),str(r.edges),'multistart simulated annealing + 1-swap descent',False,float(r.best_residual)<1e-6])
 d=pd.DataFrame(rows,columns=['N','M','k','best_residual_found','best_edges','method','globally_certified_for_this_k','below_1e-6']).sort_values(['N','k']);d.to_csv(a.outdir/'N14_N16_all_cardinalities_best_results.csv',index=False)
 sr=[]
 for N,g in d.groupby('N'):
  b=g.loc[g.best_residual_found.idxmin()]; good=g[g['below_1e-6']];sr.append([N,int(b.M),int(g.k.min()),int(g.k.max()),float(b.best_residual_found),int(b.k),b.best_edges,'' if good.empty else ','.join(map(str,good.k.astype(int))),len(good)])
 s=pd.DataFrame(sr,columns=['N','M','k_min_searched','k_max_searched','best_residual_found','best_k','best_edges','k_below_1e-6','count_k_below_1e-6']);s.to_csv(a.outdir/'N14_N16_summary.csv',index=False)
 fig,ax=plt.subplots(figsize=(8,5));
 for N,g in d.groupby('N'):ax.semilogy(g.k,g.best_residual_found,marker='o',label=f'N={N}')
 ax.axhline(1e-6,ls=':',lw=.8);ax.set_xlabel('subset cardinality k');ax.set_ylabel('best residual');ax.set_title('N=14..16 best closure residual vs k');ax.legend();fig.tight_layout();fig.savefig(a.outdir/'N14_N16_best_residual_vs_k.png',dpi=180);plt.close(fig)
 piv=d.pivot(index='N',columns='k',values='best_residual_found');fig,ax=plt.subplots(figsize=(9,3.5));im=ax.imshow(np.log10(piv.to_numpy()),aspect='auto');ax.set_yticks(range(len(piv.index)),piv.index);ax.set_xticks(range(len(piv.columns)),piv.columns);ax.set_xlabel('k');ax.set_ylabel('N');ax.set_title('log10 best closure residual');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(a.outdir/'N14_N16_residual_heatmap.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
