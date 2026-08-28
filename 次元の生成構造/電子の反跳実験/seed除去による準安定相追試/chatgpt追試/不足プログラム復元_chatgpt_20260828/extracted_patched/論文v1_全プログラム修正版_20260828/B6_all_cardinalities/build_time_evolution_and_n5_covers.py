#!/usr/bin/env python3
"""Generate the B-6 N14/N5 time-series, N5 exact covers, and comparison plot."""
from __future__ import annotations
import argparse,itertools
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt

def labels_for_df(d):
 if 'i' in d and 'j' in d:return np.array([f'{int(i)}-{int(j)}' for i,j in zip(d.i,d.j)],object)
 m=int(d.groupby('step').size().iloc[0]);n=int((1+np.sqrt(1+8*m))/2);ea,eb=np.triu_indices(n,k=1);mp=np.array([f'{i+1}-{j+1}' for i,j in zip(ea,eb)],object);return mp[d.edge_index.to_numpy(int)]
def evolution(path,edges):
 d=pd.read_csv(path);wanted=set(edges.split(','));laball=labels_for_df(d);d=d.copy();d['_lab']=laball;rows=[]
 for st,g in d.groupby('step',sort=True):z=g[g._lab.isin(wanted)].z2_re.to_numpy()+1j*g[g._lab.isin(wanted)].z2_im.to_numpy();rows.append([int(st),float(abs(z.sum())/np.abs(z).sum())])
 return pd.DataFrame(rows,columns=['step','residual'])
def n5_pairs_and_covers(final):
 d=pd.read_csv(final).reset_index(drop=True);q=d.z2_re.to_numpy()+1j*d.z2_im.to_numpy();pairs=[]
 for a,b in itertools.combinations(range(len(d)),2):
  r=float(abs(q[[a,b]].sum())/np.abs(q[[a,b]]).sum())
  if r<1e-6:pairs.append(((a,b),r,f'{int(d.iloc[a].i)}-{int(d.iloc[a].j)}+{int(d.iloc[b].i)}-{int(d.iloc[b].j)}'))
 covers=[]
 def rec(used,ch):
  if len(used)==10:covers.append(ch.copy());return
  first=min(set(range(10))-used)
  for p in pairs:
   if first in p[0] and not any(x in used for x in p[0]):rec(used|set(p[0]),ch+[p])
 rec(set(),[]);rows=[]
 for i,c in enumerate(covers):rows.append([i,5,max(p[1] for p in c),' | '.join(p[2] for p in c)])
 return pd.DataFrame(rows,columns=['cover_id','blocks','max_block_residual','pair_blocks'])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--search-results',type=Path,required=True);ap.add_argument('--n14-all-steps',type=Path,required=True);ap.add_argument('--n5-all-steps',type=Path,required=True);ap.add_argument('--n5-final',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True);a=ap.parse_args();a.outdir.mkdir(parents=True,exist_ok=True);s=pd.read_csv(a.search_results);r=s[(s.N==14)&(s.k==6)].iloc[0];t14=evolution(a.n14_all_steps,str(r.best_edges));t14.to_csv(a.outdir/'N14_best_k6_candidate_time_evolution.csv',index=False)
 # historical N5 best nontrivial pair is found directly by minimum 2-edge residual below 1e-6
 f=pd.read_csv(a.n5_final);q=f.z2_re.to_numpy()+1j*f.z2_im.to_numpy();best=(np.inf,None)
 for x in itertools.combinations(range(len(f)),2):
  rr=float(abs(q[list(x)].sum())/np.abs(q[list(x)]).sum())
  if rr<best[0]:best=(rr,x)
 ed=','.join(f'{int(f.iloc[i].i)}-{int(f.iloc[i].j)}' for i in best[1]);t5=evolution(a.n5_all_steps,ed);t5.to_csv(a.outdir/'N5_best_nontrivial_pair_time_evolution.csv',index=False);n5_pairs_and_covers(a.n5_final).to_csv(a.outdir/'N5_nontrivial_pair_exact_covers.csv',index=False)
 fig,ax=plt.subplots(figsize=(8,5));ax.semilogy(t5.step,np.maximum(t5.residual,1e-18),label='N=5 best pair');ax.semilogy(t14.step,np.maximum(t14.residual,1e-18),label='N=14 k=6 candidate');ax.axhline(1e-6,ls=':',lw=.8);ax.set_xlabel('step');ax.set_ylabel('closure residual');ax.set_title('N=5 benchmark vs N=14 quasi-closure');ax.legend();fig.tight_layout();fig.savefig(a.outdir/'N14_candidate_vs_N5_benchmark_time_evolution.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
