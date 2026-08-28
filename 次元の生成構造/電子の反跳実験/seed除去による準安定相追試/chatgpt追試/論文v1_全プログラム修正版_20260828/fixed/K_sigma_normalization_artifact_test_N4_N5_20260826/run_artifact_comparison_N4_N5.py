#!/usr/bin/env python3
# Reproducibility runner for K/sigma normalization artifact test.
# Requires run_n_scaling_lowrank_v1_NORMALIZED_ORIGINAL.py and run_n_scaling_lowrank_v1_RAW_K.py
# Conditions: same parent, same zero-closure seed, same Z0, same GAMMA, same wp0.
# Branch difference: only K/sigma normalization inside the linear rotation exp((ANGLE/sigma) K) vs exp(ANGLE K).

import importlib.util, math, numpy as np, pandas as pd
from pathlib import Path

HERE=Path(__file__).resolve().parent
N_VALUES=(4,5)
STEPS=5000
DELTA=1e-12

def loadmod(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

norm=loadmod(HERE/"run_n_scaling_lowrank_v1_NORMALIZED_ORIGINAL.py","norm")
raw=loadmod(HERE/"run_n_scaling_lowrank_v1_RAW_K.py","raw")

def obs(Z,p,q):
    Zp=Z-p*(p@Z)-q*(q@Z)
    Ht=float(np.real(np.conj(Z)@Z))
    Ho=float(np.real(np.conj(Zp)@Zp))
    Hp=Ht-Ho
    return Ht,Hp,Ho,math.sqrt(max(Ho,0.0)),Ho/Ht,complex(Z@Z)

def run_branch(mod,N,Z0,p,q,wp0):
    sys=mod.LowRankSystem(N); Z=Z0.copy(); wp=wp0.copy(); rows=[]
    for t in range(STEPS+1):
        Ht,Hp,Ho,Ao,f,ztz=obs(Z,p,q)
        sys.set_state(Z)  # FIX4
        se,wpn=sys.sigma_max_power(wp)
        sex=float(sys.sigma_spectrum()[0])
        rows.append((t,Ht,Hp,Ho,Ao,f,ztz.real,ztz.imag,abs(ztz),se,sex))
        if t<STEPS:
            Z=sys.linear_rotation_step(Z,se); wp=wpn
    return pd.DataFrame(rows,columns=["step","H_total","H_parallel","H_perp","A_perp","f","ztz_re","ztz_im","abs_ztz","sigma_est","sigma_exact"])

for N in N_VALUES:
    seed=40260722+1000*N
    sys0=norm.LowRankSystem(N)
    rng=np.random.default_rng(seed)
    v,residual,sig=norm.make_parent(sys0,rng,iters=1200,tol=1e-12)
    Z0=v.copy()  # FIX2: 外部 seed と正規化を除去（DELTA は不使用）
    p=v.real/np.linalg.norm(v.real)
    q=v.imag-(v.imag@p)*p
    q=q/np.linalg.norm(q)
    wp0=rng.normal(size=sys0.m)
    run_branch(norm,N,Z0,p,q,wp0).to_csv(HERE/f"N{N}_normalized_K_raw_observables.csv",index=False)
    run_branch(raw,N,Z0,p,q,wp0).to_csv(HERE/f"N{N}_raw_K_raw_observables.csv",index=False)
