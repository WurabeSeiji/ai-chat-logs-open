#!/usr/bin/env python3
"""Run exact k=2..4, MITM k=5..6, and deterministic annealing k>=7 for N=14..16.

The annealing args come from an explicit CSV. If that CSV is the bundled
reproduction_anneal_args.csv, the run is reproducible but NOT historical.
"""
from __future__ import annotations
import argparse,subprocess,time,sys
from pathlib import Path
import pandas as pd
from mitm_repro import best

def parse_anneal_line(s):
 r,e=s.strip().split('\t',1);return float(r),e.rstrip(';').replace(';',',')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,required=True);ap.add_argument('--anneal-args',type=Path,required=True);ap.add_argument('--workdir',type=Path,required=True);a=ap.parse_args();w=a.workdir;w.mkdir(parents=True,exist_ok=True);here=Path(__file__).resolve().parent
 exact=w/'exact_k234_bin';ann=w/'anneal_subsets_bin';subprocess.run(['g++','-O3','-std=c++17',str(here/'exact_k234.cpp'),'-o',str(exact)],check=True);subprocess.run(['g++','-O3','-std=c++17',str(here/'anneal_subsets.cpp'),'-o',str(ann)],check=True)
 mitmrows=[]
 for N in (14,15,16):
  src=a.source_dir/f'SOURCE_N{N}_step5000_final_edges.csv'; subprocess.run([str(exact),str(N),str(src),str(w/f'exact_N{N}_k234.csv')],check=True)
  df=pd.read_csv(src)
  for k in (5,6):
   t=time.time();r,idx=best(df,k);rt=time.time()-t;edges=','.join(f'{int(df.iloc[i].i)}-{int(df.iloc[i].j)}' for i in idx);mitmrows.append([N,len(df),k,r,edges,rt])
 pd.DataFrame(mitmrows,columns=['N','M','k','best_residual','edges','runtime_s']).to_csv(w/'mitm_N14_N16_k56.csv',index=False)
 args=pd.read_csv(a.anneal_args);arows=[]
 for _,r in args.iterrows():
  N=int(r.N);k=int(r.k);steps=int(r.steps);seed=int(r.seed);src=a.source_dir/f'SOURCE_N{N}_step5000_final_edges.csv';t=time.time();p=subprocess.run([str(ann),str(N),str(k),str(src),str(steps),str(seed)],capture_output=True,text=True,check=True);rt=time.time()-t;br,edges=parse_anneal_line(p.stdout); rp=subprocess.run([sys.executable,str(here/'refine_subset_repro.py'),'--input',str(src),'--edges',edges],capture_output=True,text=True,check=True); br,edges=parse_anneal_line(rp.stdout); arws=[N,N*(N-1)//2,k,br,edges,steps,seed,rt,str(r.get('status',''))];arows.append(arws);print(N,k,br)
 pd.DataFrame(arows,columns=['N','M','k','best_residual','edges','steps','seed','runtime_s','args_status']).to_csv(w/'anneal_results.csv',index=False)
 print('method results:',w)
if __name__=='__main__':main()
