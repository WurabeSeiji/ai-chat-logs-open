#!/usr/bin/env python3
# Reproducibility runner for K/sigma normalization artifact test.
# Requires run_n_scaling_lowrank_v1_NORMALIZED_ORIGINAL.py and run_n_scaling_lowrank_v1_RAW_K.py
# Conditions: same parent, same zero-closure seed, same Z0, same GAMMA, same wp0.
# R2 abolished: raw branch only (linear rotation exp(ANGLE K), amplitude-aware K).

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

# R2 廃止: K/σ 正規化枝（NORMALIZED_ORIGINAL）は実行しない
raw=loadmod(HERE/"run_n_scaling_lowrank_v1_RAW_K.py","raw")

def obs(Z,p,q):
    Zp=Z-p*(p@Z)-q*(q@Z)
    Ht=float(np.real(np.conj(Z)@Z))
    Ho=float(np.real(np.conj(Zp)@Zp))
    Hp=Ht-Ho
    return Ht,Hp,Ho,math.sqrt(max(Ho,0.0)),Ho/Ht,complex(Z@Z)

def run_branch(mod,N,Z0,p,q,wp0):
    sys=mod.LowRankSystem(N); Z=Z0.copy(); rows=[]  # S2: wp は使わない
    for t in range(STEPS+1):
        Ht,Hp,Ho,Ao,f,ztz=obs(Z,p,q)
        sys.set_state(Z)  # A4
        sex=float(sys.sigma_spectrum()[0])  # A6(b): 実際の生成子の σ₁（sigma_est 列も同値）
        rows.append((t,Ht,Hp,Ho,Ao,f,ztz.real,ztz.imag,abs(ztz),sex,sex))
        if t<STEPS:
            Z=sys.linear_rotation_step(Z)  # R1
    return pd.DataFrame(rows,columns=["step","H_total","H_parallel","H_perp","A_perp","f","ztz_re","ztz_im","abs_ztz","sigma_est","sigma_exact"])

for N in N_VALUES:
    seed=40260722+1000*N
    sys0=raw.LowRankSystem(N)
    rng=np.random.default_rng(seed)
    v,residual,sig=raw.make_parent(sys0,rng,iters=1200,tol=1e-12)
    Z0=v.copy()  # A2/A3/S1: 外部 seed も正規化も無し
    p=v.real/np.linalg.norm(v.real)
    q=v.imag-(v.imag@p)*p
    q=q/np.linalg.norm(q)
    wp0=None  # S2
    run_branch(raw,N,Z0,p,q,wp0).to_csv(HERE/f"N{N}_raw_K_raw_observables.csv",index=False)
