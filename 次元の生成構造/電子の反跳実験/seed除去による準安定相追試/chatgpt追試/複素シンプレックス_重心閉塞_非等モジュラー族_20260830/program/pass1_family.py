#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス1：クラス重み付き非等モジュラー族の実測（N=3〜16）。
(a) 各 N で q、族の次元 q−3（理論）、Σa_cω^c=0 の実線形条件の rank を数値で確認（理論との照合）。
(b) 等モジュラー点と、cos モード k=2..q−2 の族メンバー（ε=0.6）で自己無撞着残差・μ・閉塞・rank・Takagi 軸・|z| の幅を実測。
(c) 対照：Σa_cω^c≠0（k=1）では自己無撞着が破れることを実測。
(d) 1 パラメータ路 ε∈[−0.9,0.9]（k=2）に沿う Takagi 軸を data/family_path_N*.csv に保存（q≥4 の N のみ）。
出力：results/family_table.csv, results/family_members.csv, data/family_path_N*.csv, results/pass1.log（標準出力）"""
import os, csv, sys
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R2=1.0/15.0   # 基礎文書と同じスケール規約（等モジュラー点で r²=1/15）
rows=[]; members=[]
for N in range(3,17):
    th,cls,q=phases(N); A=adjacency(N); M=len(cls)
    # (a) 線形条件の rank：写像 a∈R^q → Σ a_c ω^c ∈ C ≅ R^2
    w=np.exp(2j*np.pi*np.arange(q)/q); L=np.vstack([w.real,w.imag]) if N!=3 else np.zeros((2,3))
    rank_cond=int(np.linalg.matrix_rank(L,tol=1e-12)) if N!=3 else 0
    dim_theory=family_dim(q) if N!=3 else 0
    dim_measured=(q-rank_cond-1) if N!=3 else 0   # −1 は全体スケール
    # 等モジュラー点
    a_eq=np.full(q,R2) if N!=3 else np.full(3,R2)
    v=state(N,a_eq) if N!=3 else np.sqrt(R2)*np.exp(1j*th)
    sc=selfconsistency(N,v,A); gt=gram_takagi(N,v); rnd=is_round(gt['axes'],gt['rank'])
    uax=sorted(set(np.round(gt['axes'][:gt['rank']],9)),reverse=True)
    rows.append(dict(N=N,M=M,q=q,cond_rank=rank_cond,dim_theory=dim_theory,dim_measured=dim_measured,
                     eq_residual=sc['residual'],eq_mu=sc['mu'],eq_mu_over_mean_amp2=sc['mu_over_mean_amp2'],
                     eq_global=sc['global_closure'],eq_local=sc['local_closure'],eq_rank=gt['rank'],eq_round=rnd,
                     eq_axes_distinct=' '.join(f'{x:.9f}' for x in uax),eq_vertex_null=gt['vertex_null_dev']))
    print(f"N={N:2d} M={M:3d} q={q:2d} 条件rank={rank_cond} 族次元 理論={dim_theory} 実測={dim_measured} | 等モジュラー: 残差={sc['residual']:.1e} μ/⟨|z|²⟩={sc['mu_over_mean_amp2']:+.6f} 大域={sc['global_closure']:.1e} 局所={sc['local_closure']:.1e} rank={gt['rank']} 丸い={rnd} 軸={uax}")
    members.append(dict(N=N,q=q,kind='equimodular',k=0,eps=0.0,T=0.0,residual=sc['residual'],mu=sc['mu'],
                        global_closure=sc['global_closure'],local_closure=sc['local_closure'],rank=gt['rank'],
                        amp_min=float(abs(v).min()),amp_max=float(abs(v).max()),amp_spread=float(np.ptp(abs(v))),
                        round=rnd,axes=' '.join(f'{x:.6f}' for x in gt['axes'][:gt['rank']]),vertex_null=gt['vertex_null_dev'],embed_err=gt['embed_err']))
    if N==3: continue
    # (b) 族メンバー
    for k in range(2,q-1):
        a=R2*cos_mode(q,k,0.6); T=float(abs((a*w).sum()))
        v=state(N,a); sc=selfconsistency(N,v,A); gt=gram_takagi(N,v); rnd=is_round(gt['axes'],gt['rank'])
        members.append(dict(N=N,q=q,kind='member',k=k,eps=0.6,T=T,residual=sc['residual'],mu=sc['mu'],
                            global_closure=sc['global_closure'],local_closure=sc['local_closure'],rank=gt['rank'],
                            amp_min=float(abs(v).min()),amp_max=float(abs(v).max()),amp_spread=float(np.ptp(abs(v))),
                            round=rnd,axes=' '.join(f'{x:.6f}' for x in gt['axes'][:gt['rank']]),vertex_null=gt['vertex_null_dev'],embed_err=gt['embed_err']))
        print(f"      族 k={k} ε=0.6: |Σaω^c|={T:.1e} 残差={sc['residual']:.1e} μ={sc['mu']:+.6f} 大域={sc['global_closure']:.1e} 局所={sc['local_closure']:.1e} rank={gt['rank']} |z|∈[{abs(v).min():.4f},{abs(v).max():.4f}] 丸い={rnd} 頂点ヌル={gt['vertex_null_dev']:.1e}")
    # (c) 対照 k=1
    a=R2*cos_mode(q,1,0.6); T=float(abs((a*w).sum())); v=state(N,a); sc=selfconsistency(N,v,A); gt=gram_takagi(N,v)
    members.append(dict(N=N,q=q,kind='control_k1',k=1,eps=0.6,T=T,residual=sc['residual'],mu=sc['mu'],
                        global_closure=sc['global_closure'],local_closure=sc['local_closure'],rank=gt['rank'],
                        amp_min=float(abs(v).min()),amp_max=float(abs(v).max()),amp_spread=float(np.ptp(abs(v))),
                        round=False,axes='',vertex_null=gt['vertex_null_dev'],embed_err=gt['embed_err']))
    print(f"      対照 k=1: |Σaω^c|={T:.2e} 残差={sc['residual']:.2e} 大域閉塞={sc['global_closure']:.2e}（自己無撞着でない）")
    # (d) 1 パラメータ路
    if q>=4:
        with open(os.path.join(ROOT,'data',f'family_path_N{N}.csv'),'w',newline='') as f:
            wr=csv.writer(f); wr.writerow(['eps','residual','mu','global_closure','local_closure','rank','amp_min','amp_max']+[f'axis{i+1}' for i in range(N-1)])
            for eps in np.linspace(-0.9,0.9,37):
                a=R2*cos_mode(q,2,eps); v=state(N,a); sc=selfconsistency(N,v,A); gt=gram_takagi(N,v)
                wr.writerow([f'{eps:.3f}',f"{sc['residual']:.3e}",f"{sc['mu']:.9f}",f"{sc['global_closure']:.3e}",f"{sc['local_closure']:.3e}",gt['rank'],f'{abs(v).min():.6f}',f'{abs(v).max():.6f}']+[f'{x:.9f}' for x in gt['axes'][:N-1]])
with open(os.path.join(ROOT,'results','family_table.csv'),'w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
with open(os.path.join(ROOT,'results','family_members.csv'),'w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(members[0].keys())); wr.writeheader(); wr.writerows(members)
print("PASS1 OK")
