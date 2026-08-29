# -*- coding: utf-8 -*-
"""連続時間の流れ ż = K_amp(z)·z の相対平衡 v における線形化スペクトル（共回転系）。
v は K(v)v = iσ_max v（相対平衡）。線形化 L δ = K(v)δ + (DK(v)[δ]) v を実 2M で中心差分し、共回転系 L_rot = L − σ_max·J（J は i 倍の実表現）。
ハミルトン系なので固有値は ±λ, ±λ̄ の四つ組。max Re λ > 0 ⇔ 流れの不安定性（a = 2·max Re λ が H⊥/H の τ あたり成長率）。
周波数（Im λ）の分布と、不安定固有値がどの周波数で生じるか（Krein 衝突の位置）を記録。
usage: python3 flow_linearization_spectrum.py Nmin Nmax seed0 seed1  → data/flow_spectrum_N{Nmin}_{Nmax}_s{seed0}_{seed1}.csv"""
import os, sys, math, importlib.util, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
Nmin,Nmax,S0,S1=map(int,sys.argv[1:5])
rows=[]
for N in range(Nmin,Nmax+1):
    for seed in range(S0,S1):
        s=eng.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N+seed)
        try: v,res,sig=eng.make_parent(s,rng)
        except RuntimeError: rows.append(dict(N=N,seed=seed,parent_ok=False)); continue
        A=eng._adjacency(s); M=len(v); f=lambda z: eng._K_amplitude_aware(A,z)@z
        fv=f(v); sig_re=float((np.vdot(v,fv)/np.vdot(v,v)).imag); rel_def=float(np.linalg.norm(fv-1j*sig_re*v)/np.linalg.norm(fv))
        to_r=lambda z: np.concatenate([z.real,z.imag]); h=1e-6; L=np.zeros((2*M,2*M))
        for k in range(2*M):
            e=np.zeros(2*M); e[k]=h; dd=e[:M]+1j*e[M:]; L[:,k]=(to_r(f(v+dd))-to_r(f(v-dd)))/(2*h)
        Jm=np.block([[np.zeros((M,M)),-np.eye(M)],[np.eye(M),np.zeros((M,M))]])  # i 倍
        Lrot=L-sig_re*Jm; ev=np.linalg.eigvals(Lrot); re=ev.real; im=np.abs(ev.imag)
        mx=float(re.max()); k=int(np.argmax(re)); unstable=mx>1e-6
        freqs=np.sort(im[im>1e-9]); dedup=[]; 
        for w in freqs:
            if not dedup or w-dedup[-1]>1e-6: dedup.append(float(w))
        rows.append(dict(N=N,seed=seed,parent_ok=True,M=M,parent_residual=float(res),sigma_max=float(sig[0]),sigma_from_flow=sig_re,rel_equilibrium_defect=rel_def,max_Re=mx,a_flow_per_tau=2*mx,n_unstable=int((re>1e-6).sum()),unstable_freq=float(im[k]) if unstable else float("nan"),n_zero_freq=int((im<=1e-9).sum()),n_distinct_freq=len(dedup),freq_min=dedup[0] if dedup else float("nan"),freq_max=dedup[-1] if dedup else float("nan"),freqs=";".join(f"{w:.6f}" for w in dedup)))
        print(f"N={N} seed={seed}: σ={sig_re:.5f} def={rel_def:.1e} maxRe={mx:.3e} a={2*mx:.4f} n_unst={int((re>1e-6).sum())} unst_freq={rows[-1]['unstable_freq']:.4f} distinct freqs={len(dedup)} [{dedup[0] if dedup else 0:.4f}..{dedup[-1] if dedup else 0:.4f}]",flush=True)
D=pd.DataFrame(rows); out=os.path.join(ROOT,"data",f"flow_spectrum_N{Nmin}_{Nmax}_s{S0}_{S1}.csv"); D.to_csv(out,index=False); print("saved",out)
ok=D[D.parent_ok==True]; print(ok.groupby("N").agg(n=("seed","count"),n_unstable=("a_flow_per_tau",lambda x:(x>1e-4).sum()),a_max=("a_flow_per_tau","max"),n_distinct_freq=("n_distinct_freq","median"),freq_min=("freq_min","min")).to_string())
