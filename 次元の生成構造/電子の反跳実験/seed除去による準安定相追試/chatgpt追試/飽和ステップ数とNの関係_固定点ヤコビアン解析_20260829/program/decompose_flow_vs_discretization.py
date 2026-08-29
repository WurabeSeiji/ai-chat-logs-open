# -*- coding: utf-8 -*-
"""成長率の分解：λ_f(ANGLE) = a·ANGLE + b·ANGLE²。
a = 連続時間の流れ自体の不安定性（τ あたりの成長率、刻みに依らない）、b·ANGLE² = 線形回転の刻み由来（N=16 で b=1/8 を確認）。
N=4..20 × 親 5 実現について L=124 と L=248 の 2 点から a, b を解く。出力: data/flow_vs_discretization.csv, figures/flow_rate_a_vs_N.png"""
import os, math, importlib.util, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"]=["Hiragino Sans","Arial Unicode MS","DejaVu Sans"]
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
def lam_f(N,v,ang,h=1e-6):
    s=eng.LowRankSystem(N); A=eng._adjacency(s); M=len(v); to_r=lambda z: np.concatenate([z.real,z.imag])
    F=lambda z: eng._exp_step(eng._K_amplitude_aware(A,z),z,ang); phi=float(np.angle(np.vdot(v,F(v)))); J=np.zeros((2*M,2*M))
    for k in range(2*M):
        e=np.zeros(2*M); e[k]=h; dd=e[:M]+1j*e[M:]; J[:,k]=(to_r(F(v+dd))-to_r(F(v-dd)))/(2*h)
    c,sn=np.cos(-phi),np.sin(-phi); R=np.block([[c*np.eye(M),-sn*np.eye(M)],[sn*np.eye(M),c*np.eye(M)]]); return 2*math.log(np.abs(np.linalg.eigvals(R@J)).max())
rows=[]; a1,a2=2*math.pi/124,2*math.pi/248
for N in range(4,21):
    for seed in range(5):
        s=eng.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N+seed)
        try: v,res,sig=eng.make_parent(s,rng)
        except RuntimeError: continue
        l1,l2=lam_f(N,v,a1),lam_f(N,v,a2)
        # l1 = a·a1 + b·a1², l2 = a·a2 + b·a2²  (a2=a1/2) → b = (l1 − 2 l2)/(a1² − 2 a2²) = (l1−2l2)/(a1²/2),  a = (l1 − b a1²)/a1
        b=(l1-2*l2)/(a1**2-2*a2**2); a=(l1-b*a1**2)/a1
        rows.append(dict(N=N,seed=seed,parent_residual=float(res),sigma_max=float(sig[0]),lambda_f_L124=l1,lambda_f_L248=l2,a_flow_per_tau=a,b_discretization=b,artifact_share_L124=b*a1**2/l1))
    sub=pd.DataFrame([r for r in rows if r["N"]==N]); print(f"N={N}: a(流れ,/τ)={np.round(sub.a_flow_per_tau.to_numpy(),4)}  b={np.round(sub.b_discretization.to_numpy(),3)}",flush=True)
D=pd.DataFrame(rows); D.to_csv(os.path.join(ROOT,"data","flow_vs_discretization.csv"),index=False)
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for N,g in D.groupby("N"): ax[0].plot([N]*len(g),np.maximum(g.a_flow_per_tau,1e-6),"o",c="#d62728",ms=5,alpha=.7); ax[1].plot([N]*len(g),g.b_discretization,"o",c="#1f77b4",ms=5,alpha=.7)
ax[0].set_yscale("log"); ax[0].set_xlabel("N"); ax[0].set_ylabel("a：流れの成長率（H⊥/H、/τ）"); ax[0].set_title("連続時間の不安定性 a（親 5 実現）"); ax[0].grid(alpha=.3)
ax[1].axhline(0.125,ls="--",c="gray",lw=.8); ax[1].text(4.2,0.13,"1/8",color="gray"); ax[1].set_xlabel("N"); ax[1].set_ylabel("b：刻み由来 λ_f = b·ANGLE²"); ax[1].set_title("刻み由来の係数 b"); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(ROOT,"figures","flow_rate_a_vs_N.png"),dpi=150); print(D.groupby("N").agg(a_min=("a_flow_per_tau","min"),a_med=("a_flow_per_tau","median"),a_max=("a_flow_per_tau","max"),b_med=("b_discretization","median"),b_min=("b_discretization","min"),b_max=("b_discretization","max")).to_string())
