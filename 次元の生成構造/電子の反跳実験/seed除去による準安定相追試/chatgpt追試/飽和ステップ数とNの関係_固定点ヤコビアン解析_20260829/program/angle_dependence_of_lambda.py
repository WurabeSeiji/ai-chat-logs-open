# -*- coding: utf-8 -*-
"""λ_G の刻み角依存性。走行の親 Z[0]（N=6, 16）を固定し、写像 F_L(z)=exp((2π/L)K_amp(z))z の共回転モノドロミー最大乗数を L=124..1984 で計算。
λ_f = 2λ_G が ANGLE² に比例するなら、τ あたりの成長率 λ_f/ANGLE は ANGLE→0 で 0 になる（刻み由来）。出力: data/angle_dependence.csv"""
import os, math, importlib.util, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); UP=os.path.dirname(ROOT)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
RUNS={6:"N6_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260829",16:"N16_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828"}
rows=[]
for N,d in RUNS.items():
    v=np.load(os.path.join(UP,d,"data","states_treatment.npz"))["Z"][0]; s=eng.LowRankSystem(N); A=eng._adjacency(s); M=len(v); to_r=lambda z: np.concatenate([z.real,z.imag])
    for L in (62,124,144,248,496,992,1984):
        ang=2*math.pi/L; F=lambda z: eng._exp_step(eng._K_amplitude_aware(A,z),z,ang); phi=float(np.angle(np.vdot(v,F(v))))
        J=np.zeros((2*M,2*M)); h=1e-6
        for k in range(2*M):
            e=np.zeros(2*M); e[k]=h; dd=e[:M]+1j*e[M:]; J[:,k]=(to_r(F(v+dd))-to_r(F(v-dd)))/(2*h)
        c,sn=np.cos(-phi),np.sin(-phi); R=np.block([[c*np.eye(M),-sn*np.eye(M)],[sn*np.eye(M),c*np.eye(M)]]); a=np.abs(np.linalg.eigvals(R@J)).max()
        lf=2*math.log(a); rows.append(dict(N=N,L=L,angle=ang,lambda_f_per_step=lf,lambda_f_over_angle2=lf/ang**2,lambda_f_per_tau=lf/ang)); print(f"N={N} L={L}: λ_f={lf:.4e}/step  λ_f/ANGLE²={lf/ang**2:.4f}  λ_f/ANGLE(τあたり)={lf/ang:.4e}",flush=True)
pd.DataFrame(rows).to_csv(os.path.join(ROOT,"data","angle_dependence.csv"),index=False)
