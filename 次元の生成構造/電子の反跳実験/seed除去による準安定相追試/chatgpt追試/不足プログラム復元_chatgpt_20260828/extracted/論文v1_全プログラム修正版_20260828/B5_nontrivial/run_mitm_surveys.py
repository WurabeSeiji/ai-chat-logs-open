#!/usr/bin/env python3
"""Regenerate MITM_N6_N13_k5_k6.csv and targeted_MITM_N14_N16.csv."""
from __future__ import annotations
import argparse,subprocess,tempfile,time,sys
from pathlib import Path
import pandas as pd
from mitm_subset_search import best

def row_mitm(N,k,src):
    df=pd.read_csv(src); t=time.time(); r,idx=best(df,k); rt=time.time()-t
    edges=','.join(f'{int(df.iloc[e].i)}-{int(df.iloc[e].j)}' for e in idx)
    return [N,len(df),k,r,edges,rt]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--source-dir',type=Path,required=True); ap.add_argument('--outdir',type=Path,required=True); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for N in range(6,14):
        ks=[6] if N<=9 else [5,6]
        for k in ks: rows.append(row_mitm(N,k,a.source_dir/f'SOURCE_N{N}_step5000_final_edges.csv'))
    cols=['N','M','k','best_residual','edges','runtime_s']; pd.DataFrame(rows,columns=cols).to_csv(a.outdir/'MITM_N6_N13_k5_k6.csv',index=False)
    rows=[row_mitm(14,6,a.source_dir/'SOURCE_N14_step5000_final_edges.csv'),row_mitm(15,6,a.source_dir/'SOURCE_N15_step5000_final_edges.csv')]
    # N16 k4: the original all-cardinality package used exact enumeration for k<=4.
    cpp=Path(__file__).with_name('exact_k234.cpp'); exe=a.outdir/'exact_k234_bin'
    subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
    tmp=a.outdir/'N16_exact_k234_tmp.csv'; t=time.time(); subprocess.run([str(exe),'16',str(a.source_dir/'SOURCE_N16_step5000_final_edges.csv'),str(tmp)],check=True,stdout=subprocess.DEVNULL); rt=time.time()-t
    e=pd.read_csv(tmp); r=e[e.k==4].iloc[0]; edges=str(r.edges).rstrip(';').replace(';',','); rows.append([16,int(r.M),4,float(r.best_residual),edges,rt])
    rows += [row_mitm(16,5,a.source_dir/'SOURCE_N16_step5000_final_edges.csv'),row_mitm(16,6,a.source_dir/'SOURCE_N16_step5000_final_edges.csv')]
    pd.DataFrame(rows,columns=cols).to_csv(a.outdir/'targeted_MITM_N14_N16.csv',index=False)
    try: exe.unlink(); tmp.unlink()
    except OSError: pass
if __name__=='__main__': main()
