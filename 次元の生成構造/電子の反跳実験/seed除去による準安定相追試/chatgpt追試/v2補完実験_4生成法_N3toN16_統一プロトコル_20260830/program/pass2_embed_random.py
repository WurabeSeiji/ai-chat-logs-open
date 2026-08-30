#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス2：主張 2「複素シンプレックスは制約を与えない」の数値実証。
N=3〜16 の各 N について乱数状態（自己無撞着でない、閉塞もしていない）100 個を作り、
 (a) 複素：d²=z² から B=−½JD²J（複素対称）を作り Takagi 分解で座標 x を構成、距離再現誤差 max|(x_i−x_j)·(x_i−x_j)−d²_ij| と rank を測る。
 (b) 実の対照：同じ状態を実シンプレックスとして読む（d²=|z|²、実対称 B）と B の負固有値の有無（半正定値でなければ R^{N−1} に埋め込めない）。
     また、実の d²≥0 で閉塞 Σd²=0 が全零以外で不可能であることの数値確認（Σ|z|²>0）。
出力：results/embed_random.csv、results/embed_random_summary.csv"""
import os, sys, csv
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import edges, gram_takagi
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng=np.random.default_rng(2026)
rows=[]; summ=[]
for N in range(3,17):
    E=edges(N); M=len(E); errs=[]; ranks=[]; negfrac=[]; mineig=[]
    for t in range(100):
        v=rng.standard_normal(M)+1j*rng.standard_normal(M)
        gt=gram_takagi(N,v); errs.append(gt['embed_err']/np.abs(v*v).max()); ranks.append(gt['rank'])
        D=np.zeros((N,N))
        for val,(i,j) in zip(np.abs(v)**2,E): D[i,j]=D[j,i]=val
        J=np.eye(N)-np.ones((N,N))/N; Br=-0.5*J@D@J; w=np.linalg.eigvalsh(Br)
        mineig.append(float(w.min()/abs(w).max())); negfrac.append(float((w<-1e-12*abs(w).max()).sum()))
        rows.append(dict(N=N,trial=t,complex_embed_err_rel=f'{errs[-1]:.2e}',complex_rank=ranks[-1],real_B_min_eig_rel=f'{mineig[-1]:.3e}',real_B_n_negative=int(negfrac[-1]),sum_abs_z2=f'{(np.abs(v)**2).sum():.4f}'))
    psd=sum(1 for x in negfrac if x==0)
    summ.append(dict(N=N,M=M,n_trials=100,complex_embed_err_rel_max=f'{max(errs):.2e}',complex_rank_all_Nminus1=int(all(r==N-1 for r in ranks)),real_PSD_count=psd,real_not_embeddable_count=100-psd,real_min_eig_rel_median=f'{np.median(mineig):.3e}'))
    print(f"N={N:2d} M={M:3d} 複素：埋め込み誤差 max={max(errs):.1e} rank=N−1 全て={all(r==N-1 for r in ranks)} | 実（d²=|z|²）：半正定値（埋め込み可）{psd}/100、負固有値あり {100-psd}/100、最小固有値中央値 {np.median(mineig):+.3f}")
with open(os.path.join(ROOT,'results','embed_random.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(os.path.join(ROOT,'results','embed_random_summary.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
print("PASS2 OK")
