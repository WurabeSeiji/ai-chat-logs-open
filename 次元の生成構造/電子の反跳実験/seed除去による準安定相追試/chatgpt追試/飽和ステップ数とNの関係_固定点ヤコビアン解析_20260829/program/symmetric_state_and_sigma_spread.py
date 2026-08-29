# -*- coding: utf-8 -*-
"""(a) 各親の σ スペクトル（iK_amp の正固有値）の「対称性からのずれ」：非最大 σ の相対幅 spread=(max−min)/mean と流れの成長率 a の関係。
(b) Z_N 対称候補 z_ab = ω^{a+b}（ω=e^{2πi/N}）：等モジュラー・ヌル単体（N≠2,4 の約数）で、自己無撞着残差・σ スペクトル・線形安定性を N=5..20 で計算。
    もし全 N で厳密な相対平衡かつ安定なら、「大 N で親がこの対称解に収束する」を示せば証明に近づく。
出力: data/sigma_spread_vs_a.csv, data/symmetric_candidate.csv"""
import os, math, importlib.util, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
def flow_a(A,v):
    M=len(v); f=lambda z: eng._K_amplitude_aware(A,z)@z; fv=f(v); sg=float((np.vdot(v,fv)/np.vdot(v,v)).imag); to_r=lambda z: np.concatenate([z.real,z.imag]); h=1e-6; L=np.zeros((2*M,2*M))
    for k in range(2*M):
        e=np.zeros(2*M); e[k]=h; dd=e[:M]+1j*e[M:]; L[:,k]=(to_r(f(v+dd))-to_r(f(v-dd)))/(2*h)
    Jm=np.block([[np.zeros((M,M)),-np.eye(M)],[np.eye(M),np.zeros((M,M))]]); ev=np.linalg.eigvals(L-sg*Jm); return 2*float(ev.real.max()), sg
def sigma_spec(A,v):
    K=eng._K_amplitude_aware(A,v); w=np.sort(np.linalg.eigvalsh(1j*K))[::-1]; return w[w>1e-9]
# (a)
rows=[]
for N in range(4,15):
    for seed in range(25 if 9<=N<=14 else 5):
        s=eng.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N+seed)
        try: v,res,sig=eng.make_parent(s,rng)
        except RuntimeError: continue
        A=eng._adjacency(s); w=sigma_spec(A,v); rest=w[1:]; spread=float((rest.max()-rest.min())/rest.mean()) if len(rest)>1 else 0.0
        a,sg=flow_a(A,v); rows.append(dict(N=N,seed=seed,sigma_max=float(w[0]),n_sigma=len(w),sigma_rest_mean=float(rest.mean()) if len(rest) else float("nan"),sigma_rest_spread=spread,sigma_2nd=float(w[1]) if len(w)>1 else float("nan"),sigma_min=float(w[-1]),a_flow_per_tau=a))
D=pd.DataFrame(rows); D.to_csv(os.path.join(ROOT,"data","sigma_spread_vs_a.csv"),index=False)
print("=== σ 非最大部の相対幅 vs a（N ごと：不安定親 / 安定親）")
for N,g in D.groupby("N"):
    u=g[g.a_flow_per_tau>1e-4]; st=g[g.a_flow_per_tau<=1e-4]
    print(f"N={N}: 不安定 {len(u)} 件 spread={np.round(u.sigma_rest_spread.to_numpy(),3)} σmax={np.round(u.sigma_max.to_numpy(),4)} | 安定 {len(st)} 件 spread 範囲 [{st.sigma_rest_spread.min():.3f},{st.sigma_rest_spread.max():.3f}] σmax 範囲 [{st.sigma_max.min():.4f},{st.sigma_max.max():.4f}]")
# (b)
print("=== Z_N 対称候補 z_ab = ω^{a+b}")
rows=[]
for N in range(5,21):
    s=eng.LowRankSystem(N); A=eng._adjacency(s); om=np.exp(2j*math.pi/N); v=np.array([om**(int(a)+int(b)) for a,b in zip(s.ea,s.eb)])*math.sqrt(N-1)/math.sqrt(len(s.ea))*0+np.array([om**(int(a)+int(b)) for a,b in zip(s.ea,s.eb)])
    null=abs(np.sum(v*v))/np.sum(abs(v)**2); res=eng._selfconsistency_residual(eng._K_amplitude_aware(A,v),v); w=sigma_spec(A,v)
    a,sg=flow_a(A,v) if res<1e-6 else (float("nan"),float("nan"))
    rows.append(dict(N=N,null_residual=null,selfconsistency_residual=res,sigma_max=float(w[0]),sigma_rest_spread=float((w[1:].max()-w[1:].min())/w[1:].mean()) if len(w)>2 else 0.0,n_sigma=len(w),a_flow_per_tau=a)); print(f"N={N}: |Σz²|/Σ|z|²={null:.1e} 残差={res:.1e} σ={np.round(w[:4],4)} n_σ={len(w)} a={a}")
pd.DataFrame(rows).to_csv(os.path.join(ROOT,"data","symmetric_candidate.csv"),index=False)
