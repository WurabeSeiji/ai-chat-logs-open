#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス1：4 生成法 × N=3〜16 の親を生成・受け入れ検査・保存し、走行前に共回転モノドロミー予測を固定する。
生成法（タグ）：
  mp = make_parent 等モジュラー（3 段階、original_engine.make_parent、rng=default_rng(40260721+1000N+0)、N=3 は pre_steps=40000）
  hm = 手作り等モジュラー（基礎文書の構成：偶数 1-因子分解／奇数 距離クラス／N=3 Z3。state_provider.equimodular）
  ne = 非等モジュラー（q≥4 クラス重み付き族 a_c=r̄²(1+0.6cos(4πc/q))、q≤3 は多様体上の代表点。state_provider.state）
  rb = 乱数均衡親 {S_i=0, W_i=W0}（対称性なし、seed=100+N、N≥5）
スケール規約：各 N で ‖v‖ を mp 親の ‖v‖ に揃える（スケールは回転速度にしか効かない）。
受け入れ：残差 <1e-10、|Σz²|/H <1e-12、μ≠0。局所閉塞は記録のみ（ne の N=3,4 は破れる）。
予測：Φ(z)=exp(−i(2π/L)·iK(z))z の実ヤコビアンを中心差分で求め、G=R(+φ)DΦ（φ=(2π/L)μ）の最大固有値絶対値 ρ。
判定規則（較正済み）：ρ−1>1e-3 → inflating、それ未満 → neutral（床）。λ_f=2lnρ、t50 予測=(ln0.5−ln f0)/λ_f、f0=3e-32。
出力：data/<tag>/parent_v.{npz,csv}, parent_checks.json；results/parents_predictions.csv；results/closure_step0_4methods.csv"""
import os, sys, csv, json, math
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
import original_engine as eng
from common import edges, adjacency, K_of, selfconsistency, gram_takagi
import state_provider as sp
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L=124; F0=3e-32
def step_map(N,z,A):
    K=K_of(N,z,A); w,V=np.linalg.eigh(1j*K); return V@(np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@z))
def monodromy(N,v,A,h=1e-7):
    M=len(v); x=np.r_[v.real,v.imag]; D=np.zeros((2*M,2*M))
    for j in range(2*M):
        e=np.zeros(2*M); e[j]=h; zp=step_map(N,(x+e)[:M]+1j*(x+e)[M:],A); zm=step_map(N,(x-e)[:M]+1j*(x-e)[M:],A); D[:,j]=np.r_[(zp-zm).real,(zp-zm).imag]/(2*h)
    mu=selfconsistency(N,v,A)['mu']; phi=(2*math.pi/L)*mu; c,s=math.cos(phi),math.sin(phi)
    R=np.block([[c*np.eye(M),-s*np.eye(M)],[s*np.eye(M),c*np.eye(M)]]); ev=np.linalg.eigvals(R@D)
    return float(np.abs(ev).max()), int((np.abs(ev)>1+1e-9).sum())
def SW(N,v):
    E=edges(N); S=np.zeros(N,complex); W=np.zeros(N)
    for k,(i,j) in enumerate(E): S[i]+=v[k]**2; S[j]+=v[k]**2; W[i]+=abs(v[k])**2; W[j]+=abs(v[k])**2
    return S,W
def solve_balanced(N,rng,W0,iters=100):
    E=edges(N); M=len(E); v=rng.standard_normal(M)+1j*rng.standard_normal(M); v*=np.sqrt(N*W0/2)/np.linalg.norm(v)
    def G(x):
        vv=x[:M]+1j*x[M:]; S,W=SW(N,vv); return np.r_[S.real,S.imag,W-W0]
    for it in range(iters):
        x=np.r_[v.real,v.imag]; g=G(x)
        if np.linalg.norm(g)<1e-15: break
        J=np.zeros((3*N,2*M)); h=1e-7
        for j in range(2*M):
            e=np.zeros(2*M); e[j]=h; J[:,j]=(G(x+e)-G(x-e))/(2*h)
        dx=-np.linalg.lstsq(J,g,rcond=1e-12)[0]; x=x+dx; v=x[:M]+1j*x[M:]
    return v
def save(tag,N,v,design,extra):
    A=adjacency(N); sc=selfconsistency(N,v,A); gt=gram_takagi(N,v); S,W=SW(N,v); E=edges(N)
    ok=dict(residual=sc['residual']<1e-10,global_closure=sc['global_closure']<1e-12,mu_nonzero=abs(sc['mu'])>1e-6)
    rho,nun=monodromy(N,v,A); lam=2*math.log(rho); kind='inflating' if rho-1>1e-3 else 'neutral'
    dd=os.path.join(ROOT,'data',tag); os.makedirs(dd,exist_ok=True)
    sig=np.sort(np.linalg.eigvalsh(1j*K_of(N,v,A)))[::-1]
    np.savez_compressed(os.path.join(dd,'parent_v.npz'),v=v,edges=np.array(E),color=np.zeros(len(E),int),theta=np.angle(v),design=design,r=float(np.sqrt(sc['mean_amp2'])),sigma=sig,mu=sc['mu'],residual=sc['residual'])
    with open(os.path.join(dd,'parent_v.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['edge_index','i','j','theta_deg','a_Re','b_Im','abs_z'])
        for k,e in enumerate(E): w.writerow([k,e[0],e[1],np.degrees(np.angle(v[k]))%360,v[k].real,v[k].imag,abs(v[k])])
    ax=gt['axes'][:gt['rank']]
    rep=dict(tag=tag,N=N,method=tag.split('_')[0],design=design,M=len(E),norm=float(np.linalg.norm(v)),mean_amp2=sc['mean_amp2'],amp_min=float(abs(v).min()),amp_max=float(abs(v).max()),amp_spread_rel=float(np.ptp(abs(v)**2)/(abs(v)**2).mean()),
             residual=sc['residual'],mu=sc['mu'],mu_over_mean_amp2=sc['mu_over_mean_amp2'],global_closure=sc['global_closure'],local_closure=sc['local_closure'],local_closed=bool(sc['local_closure']<1e-10),
             W_spread=float(np.ptp(W)),rank=gt['rank'],roundness=float((ax.max()-ax.min())/ax.max()),axes=[float(x) for x in ax],vertex_null=gt['vertex_null_dev'],
             pred_rho=rho,pred_rho_minus_1=rho-1,pred_lambda_f=lam,pred_n_unstable=nun,pred_kind=kind,pred_t50=(float((math.log(0.5)-math.log(F0))/lam) if lam>1e-6 else None),ok={k:bool(x) for k,x in ok.items()},**extra)
    with open(os.path.join(dd,'parent_checks.json'),'w') as f: json.dump(rep,f,indent=2,ensure_ascii=False)
    print(f"{tag}: ‖v‖={rep['norm']:.6f} |z|∈[{rep['amp_min']:.4f},{rep['amp_max']:.4f}] 残差={sc['residual']:.1e} μ={sc['mu']:+.6f} 閉塞={sc['global_closure']:.1e} 局所={sc['local_closure']:.1e} rank={gt['rank']} 丸さ={rep['roundness']:.3f} | ρ−1={rho-1:.2e} λ_f={lam:.5f} → {kind} 判定={ok}")
    if not all(ok.values()): raise SystemExit(f"ABORT {tag}: {ok}")
    return rep
rows=[]; NORM={}
for N in range(3,17):
    # --- mp ---
    s0=eng.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N+0)
    v,res,sig=eng.make_parent(s0,rng,pre_steps=(40000 if N==3 else 20000))
    NORM[N]=float(np.linalg.norm(v))
    rows.append(save(f'mp_N{N}',N,v,'make_parent_equimodular_3stage',dict(seed_rng=40260721+1000*N)))
    # --- hm ---
    v=sp.equimodular(N); v=v*NORM[N]/np.linalg.norm(v)
    rows.append(save(f'hm_N{N}',N,v,'handmade_equimodular_'+('Z3' if N==3 else '1factor' if N%2==0 else 'distance_classes'),dict()))
    # --- ne ---
    v,kind,col,q,step=sp.state(N); v=v*NORM[N]/np.linalg.norm(v)
    rows.append(save(f'ne_N{N}',N,v,'nonequimodular_'+('class_family_k2_eps0.6' if kind=='class' else 'manifold_point'),dict(q=q)))
    # --- rb ---
    if N>=5:
        rng=np.random.default_rng(100+N); W0=2*NORM[N]**2/N; v=solve_balanced(N,rng,W0)
        rows.append(save(f'rb_N{N}',N,v,'random_balanced_S0_Wconst',dict(seed_rng=100+N)))
keys=[k for k in rows[0].keys() if k not in ('axes','ok')]
with open(os.path.join(ROOT,'results','parents_predictions.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
with open(os.path.join(ROOT,'results','closure_step0_4methods.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['tag','N','method','norm','amp_spread_rel','global_closure_abs_sum_z2_over_H','local_closure_max_over_H','residual','mu'])
    for r in rows: w.writerow([r['tag'],r['N'],r['method'],f"{r['norm']:.9f}",f"{r['amp_spread_rel']:.3e}",f"{r['global_closure']:.3e}",f"{r['local_closure']:.3e}",f"{r['residual']:.3e}",f"{r['mu']:.9f}"])
json.dump(NORM,open(os.path.join(ROOT,'results','norm_by_N.json'),'w'),indent=1)
print("PASS1 OK（予測は走行前に固定）")
