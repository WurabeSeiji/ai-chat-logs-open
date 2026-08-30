#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス2b：ヤコビアン核方向の摂動が本物の自己無撞着解へつながるか（族か特異点か）を Newton 精密化で判別。
等モジュラー点 v0 から、核（スケール・位相を除く）内のランダム方向 d に δ=1e-3 だけ動かし、
F(v)=iK(v)v−μ(v)v=0 を Gauss–Newton（最小ノルム更新）で解く。‖v‖ は v0 に固定、全体位相は ⟨v0,v⟩ 実で固定。
判定：収束した解 v* について (i) 残差 (ii) 親からの距離 ‖v*−v0‖（δ と同程度なら族、δ² 程度に戻れば特異点）
(iii) |z| のクラス内ばらつき（クラス一様族の外か）(iv) 位相のクラス内ばらつき (v) 大域・局所閉塞・rank・丸さ（軸の相対幅）・μ（全収束解の中央値と最大）。
(vi) generic 点（収束解の 1 つ）でのヤコビアン核次元——等モジュラー点の核次元が対称点の特異性でなく多様体の次元であることの確認。
さらに収束解 100 個の差分 v*−v0 の PCA で「実際に到達した解集合」の次元を推定（特異値の折れ）。
出力：results/nullspace_probe.csv, results/nullspace_probe_pca_N*.csv"""
import os, csv, sys
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R2=1.0/15.0; DELTA=1e-3; NS=100
def Fvec(N,v,A):
    K=K_of(N,v,A); hv=1j*(K@v); mu=(np.vdot(v,hv)/np.vdot(v,v)).real; r=hv-mu*v; return np.r_[r.real,r.imag]
def jac(N,v,A,h=1e-6):
    x=np.r_[v.real,v.imag]; M=len(v); J=np.zeros((2*M,2*M))
    for j in range(2*M):
        e=np.zeros(2*M); e[j]=h
        J[:,j]=(Fvec(N,(x+e)[:M]+1j*(x+e)[M:],A)-Fvec(N,(x-e)[:M]+1j*(x-e)[M:],A))/(2*h)
    return J
def newton(N,v,v0,A,iters=60):
    n0=np.linalg.norm(v0)
    for it in range(iters):
        r=Fvec(N,v,A)
        if np.linalg.norm(r)<1e-14*n0: break
        J=jac(N,v,A)
        # 拘束：‖v‖=‖v0‖、Im⟨v0,v⟩=0 を行として追加
        x=np.r_[v.real,v.imag]; g1=np.r_[v.real,v.imag]/np.linalg.norm(x); g2=np.r_[-v0.imag,v0.real]
        Ja=np.vstack([J,g1,g2]); ra=np.r_[r,np.linalg.norm(x)-n0,float(np.vdot(v0,v).imag)]
        dx=-np.linalg.lstsq(Ja,ra,rcond=1e-10)[0]; x=x+dx; v=x[:len(v)]+1j*x[len(v):]
    return v,float(np.linalg.norm(Fvec(N,v,A))/n0),it
rows=[]
for N in [4,5,6,7,8]:
    th,cls,q=phases(N); A=adjacency(N); M=len(cls); v0=np.sqrt(R2)*np.exp(1j*th)
    J=jac(N,v0,A); U,s,Vt=np.linalg.svd(J); null=int((s<s.max()*1e-7).sum()); Nsp=Vt[-null:].T   # 核の正規直交基底
    # スケール・位相方向を核から除く
    triv=np.c_[np.r_[v0.real,v0.imag],np.r_[-v0.imag,v0.real]]; triv,_=np.linalg.qr(triv)
    P=Nsp-triv@(triv.T@Nsp); Q,_=np.linalg.qr(P); Q=Q[:,:null-2]
    rng=np.random.default_rng(0); sols=[]; ok=0; dists=[]; ampcls=[]; phcls=[]; back=[]; locs=[]; globs=[]; ranks=[]; rounds=[]; mus=[]
    for t in range(NS):
        d=Q@rng.standard_normal(Q.shape[1]); d/=np.linalg.norm(d); x=np.r_[v0.real,v0.imag]+DELTA*np.linalg.norm(v0)*d
        v,res,it=newton(N,x[:M]+1j*x[M:],v0,A)
        dist=np.linalg.norm(v-v0)/np.linalg.norm(v0)
        amp=abs(v); ph=np.angle(v*np.exp(-1j*th))  # 親位相からのずれ
        ac=max(np.ptp(amp[cls==c])/amp.mean() for c in range(q)); pc=max(np.ptp(ph[cls==c]) for c in range(q))
        conv=res<1e-12
        if conv:
            ok+=1; sols.append(v-v0); dists.append(dist); ampcls.append(ac); phcls.append(pc)
            s2=selfconsistency(N,v,A); g2=gram_takagi(N,v); ax_=g2['axes'][:g2['rank']]
            locs.append(s2['local_closure']); globs.append(s2['global_closure']); ranks.append(g2['rank']); rounds.append(float((ax_.max()-ax_.min())/ax_.max())); mus.append(s2['mu'])
    S=np.array([np.r_[s_.real,s_.imag] for s_ in sols]); sv=np.linalg.svd(S,compute_uv=False) if len(S) else np.array([])
    svn=sv/sv[0] if len(sv) else sv
    dim_est=int((svn>1e-3).sum()) if len(sv) else 0
    sc=selfconsistency(N,v0+sols[0],A) if sols else {}
    gt=gram_takagi(N,v0+sols[0]) if sols else {}
    # generic 点での核次元
    Jg=jac(N,v0+sols[0],A); sg=np.linalg.svd(Jg,compute_uv=False); null_g=int((sg<sg.max()*1e-7).sum())
    ax0=gram_takagi(N,v0)['axes'][:N-1]; round0=float((ax0.max()-ax0.min())/ax0.max())
    rows.append(dict(N=N,M=M,q=q,nullity=null,nontrivial_null=null-2,class_family_dim=family_dim(q),n_tried=NS,n_converged=ok,
                     dist_med=float(np.median(dists)) if dists else float('nan'),dist_min=float(np.min(dists)) if dists else float('nan'),dist_max=float(np.max(dists)) if dists else float('nan'),
                     amp_within_class_spread_med=float(np.median(ampcls)) if ampcls else float('nan'),
                     phase_within_class_spread_med=float(np.median(phcls)) if phcls else float('nan'),
                     pca_dim_est=dim_est,nullity_generic_point=null_g,
                     global_closure_med=float(np.median(globs)),global_closure_max=float(np.max(globs)),
                     local_closure_med=float(np.median(locs)),local_closure_max=float(np.max(locs)),local_closure_min=float(np.min(locs)),
                     n_local_closed_lt_1e10=int(sum(l<1e-10 for l in locs)),rank_min=int(min(ranks)),rank_max=int(max(ranks)),
                     roundness_parent=round0,roundness_med=float(np.median(rounds)),roundness_max=float(np.max(rounds)),
                     mu_parent=selfconsistency(N,v0,A)['mu'],mu_dev_max=float(np.max(np.abs(np.array(mus)-selfconsistency(N,v0,A)['mu'])))))
    with open(os.path.join(ROOT,'results',f'nullspace_probe_pca_N{N}.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['index','singular_value_rel']); [w.writerow([i,f'{x:.3e}']) for i,x in enumerate(svn)]
    r_=rows[-1]
    print(f"N={N} 核={null}(非自明 {null-2}) generic点の核={null_g} クラス族={family_dim(q)} 収束={ok}/{NS} 距離/δ 中央={np.median(dists)/DELTA:.2f} クラス内|z|幅={np.median(ampcls):.1e} クラス内位相幅={np.median(phcls):.1e} PCA次元={dim_est}")
    print(f"     大域閉塞 max={r_['global_closure_max']:.1e} | 局所閉塞 中央={r_['local_closure_med']:.1e} min={r_['local_closure_min']:.1e} max={r_['local_closure_max']:.1e} 局所閉塞保持(<1e-10)={r_['n_local_closed_lt_1e10']}/{ok} | rank∈[{r_['rank_min']},{r_['rank_max']}] | 丸さ(軸相対幅) 親={round0:.1e} 中央={r_['roundness_med']:.2e} max={r_['roundness_max']:.2e} | μ ずれ max={r_['mu_dev_max']:.1e}")
    print("   PCA 特異値(相対) 先頭:", ' '.join(f'{x:.2e}' for x in svn[:16]))
with open(os.path.join(ROOT,'results','nullspace_probe.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS2b OK")
