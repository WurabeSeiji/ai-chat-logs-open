#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス2：自己無撞着写像 F(v)=iK(v)v−μ(v)v（μ(v)=Re⟨v,iKv⟩/⟨v,v⟩）の実ヤコビアン（2M×2M）の核次元を、
等モジュラー点で N=3〜16 について数値計算。核には自明方向 2（全体スケール・全体位相）が必ず含まれる。
核次元−2 = 等モジュラー点を通る自己無撞着解集合の接空間の次元（上界：解集合が滑らかでない交差点なら過大評価）。
クラス一様仮定を置かないので、パス1 の族次元 q−3 と比較して「クラス一様族が全てか」を判定する。
F は v の 3 次多項式なので中心差分の誤差は O(h²)。h=1e-6、特異値の閾値は最大特異値×1e-7。
出力：results/jacobian_nullity.csv"""
import os, csv, sys
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R2=1.0/15.0
def F(N,x,A):
    M=len(x)//2; v=x[:M]+1j*x[M:]; K=K_of(N,v,A); hv=1j*(K@v); mu=(np.vdot(v,hv)/np.vdot(v,v)).real
    r=hv-mu*v; return np.r_[r.real,r.imag]
rows=[]
for N in range(3,17):
    th,cls,q=phases(N); A=adjacency(N); M=len(cls)
    v=np.sqrt(R2)*np.exp(1j*th); x=np.r_[v.real,v.imag]; h=1e-6
    Jm=np.zeros((2*M,2*M))
    for j in range(2*M):
        e=np.zeros(2*M); e[j]=h; Jm[:,j]=(F(N,x+e,A)-F(N,x-e,A))/(2*h)
    s=np.linalg.svd(Jm,compute_uv=False); thr=s.max()*1e-7; null=int((s<thr).sum())
    # 自明方向の確認：スケール方向 x、位相方向 (−Im v, Re v) が核にあるか
    sc=np.linalg.norm(Jm@x)/np.linalg.norm(Jm,2)/np.linalg.norm(x); ph=np.r_[-v.imag,v.real]; phn=np.linalg.norm(Jm@ph)/np.linalg.norm(Jm,2)/np.linalg.norm(ph)
    gap=s[::-1][null]/s.max() if null<2*M else float('nan')  # 核の直上の特異値（相対）
    fd=family_dim(q) if N!=3 else 0
    rows.append(dict(N=N,M=M,q=q,nullity=null,nontrivial=null-2,class_family_dim=fd,scale_dir_resid=f'{sc:.1e}',phase_dir_resid=f'{phn:.1e}',
                     smallest_nonnull_rel=f'{gap:.3e}',largest_null_rel=f'{(s[::-1][null-1]/s.max()) if null>0 else 0:.1e}'))
    print(f"N={N:2d} M={M:3d} q={q:2d} 核次元={null} 非自明={null-2} クラス一様族 q−3={fd} | 核直上 σ/σmax={gap:.2e} 核内最大={rows[-1]['largest_null_rel']} スケール/位相方向残差={sc:.0e}/{phn:.0e}")
with open(os.path.join(ROOT,'results','jacobian_nullity.csv'),'w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print("PASS2 OK")
