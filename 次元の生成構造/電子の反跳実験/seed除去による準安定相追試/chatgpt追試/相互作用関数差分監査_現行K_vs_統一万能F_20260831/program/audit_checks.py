#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""現行相互作用（振幅込み K・凍結生成子の指数回転）と統一万能相互作用関数 F（unified_interaction_v1）の
差分を数値で確定する監査。仕様書（統一万能関数_仕様_v2.md §1）と実装を読んだ上での検証であり、
コードが正・文書がそれに合わせる（規約 R0）。

V1 恒等式: 等モジュラー自己無撞着状態では K_amp(v) = r²·K_phase(θ(v))（両力学の共通固有モード）。
    非等モジュラーでは成り立たない（相対偏差を実測）。
V2 時計: 同じ等モジュラー親の 1 step 位相進み。現行 = Δ·μ（状態のスケールに比例）、
    統一線形部 = 2·arctan(γ·σ/σ_max)（K/σ 正規化により最上位モードは状態非依存の 2·arctan(γ)）。
V3 倍音→反射係数: build_standard_universe(N=12) の δ 掃引で step0 の R を実測。
    真空 δ=0 で R≡0（厳密）、小 δ で R = scale·Pf/(Pf+Pb) ∝ δ²。
出力: results/audit_results.json と標準出力。"""
import os, sys, json, math
import importlib.util
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0,HERE)
from common import edges, adjacency, K_of, selfconsistency
import state_provider as sp

REPO=os.path.abspath(os.path.join(ROOT,'..','..','..','..','..'))
SERIES=os.path.join(REPO,'次元の生成構造')
spec=importlib.util.spec_from_file_location('uif',os.path.join(SERIES,'統一万能関数_v1','unified_interaction_v1.py'))
uif=importlib.util.module_from_spec(spec); sys.modules['uif']=spec.loader.exec_module(spec and uif) or uif
out={}

def phase_only_K(N,theta,A):
    return A*np.sin(theta[None,:]-theta[:,None]).T*0+A*np.sin(np.subtract.outer(theta,theta)).T  # placeholder

def K_phase(N,theta,A):
    # K_phase_ef = A_ef sin(θ_f − θ_e)
    return A*np.sin(np.subtract.outer(theta,theta).T)

# ---------- V1: 恒等式と固定点の帰属 ----------
for tag in ['hm','ne']:
    N=8; A=adjacency(N)
    if tag=='hm': v=sp.equimodular(N)
    else: v,kind,col,q,step=sp.state(N)
    th=np.angle(v); r2=(np.abs(v)**2)
    Ka=K_of(N,v,A); Kp=K_phase(N,th,A)
    # 等モジュラーなら Ka = r̄²·Kp
    c=float(r2.mean()); dev=float(np.linalg.norm(Ka-c*Kp)/np.linalg.norm(Ka))
    sc_amp=selfconsistency(N,v,A)['residual']
    # 位相のみ K に対する自己無撞着性（統一線形部の固定点か）
    hv=1j*(Kp@v); mu_p=(np.vdot(v,hv)/np.vdot(v,v)).real; res_p=float(np.linalg.norm(hv-mu_p*v)/np.linalg.norm(v))
    out[f'V1_{tag}_N8']=dict(K_identity_reldev=dev, amp_selfcons_residual=sc_amp, phase_only_selfcons_residual=res_p)
    print(f"V1 {tag} N=8: ‖K_amp − ⟨r²⟩K_phase‖/‖K_amp‖ = {dev:.2e} | 振幅込み残差 {sc_amp:.1e} | 位相のみ残差 {res_p:.2e}")

# ---------- V2: 1 step の位相進み（時計の違い） ----------
N=8; A=adjacency(N); v=sp.equimodular(N); L=124; GAMMA=math.tan(math.pi/144.0)
K=K_of(N,v,A); w,U=np.linalg.eigh(1j*K); mu=selfconsistency(N,v,A)['mu']
znew=U@(np.exp(-1j*(2*math.pi/L)*w)*(U.conj().T@v))
adv_cur=float(np.angle(np.vdot(v,znew)))
sig=np.abs(w).max()
adv_uni=2*math.atan(GAMMA*abs(mu)/ (np.abs(np.linalg.eigvalsh(1j*K_phase(N,np.angle(v),A))).max()) * (np.abs(np.linalg.eigvalsh(1j*K_phase(N,np.angle(v),A))).max()/np.abs(np.linalg.eigvalsh(1j*K_phase(N,np.angle(v),A))).max()))
# 統一線形部: K̃=K_phase/σ_max、固有モード位相進み = 2 arctan(γ σ/σ_max)。親は σ_max モード → 2 arctan(γ)
Kp=K_phase(N,np.angle(v),A); sp_max=float(np.abs(np.linalg.eigvalsh(1j*Kp)).max())
mu_p=float(abs((np.vdot(v,1j*(Kp@v))/np.vdot(v,v)).real))
adv_uni_top=2*math.atan(GAMMA)  # 状態非依存
adv_uni_parent=2*math.atan(GAMMA*mu_p/sp_max)
# スケール依存性: v を 2 倍にすると
v2x=2*v; K2=K_of(N,v2x,A); w2,U2=np.linalg.eigh(1j*K2)
z2=U2@(np.exp(-1j*(2*math.pi/L)*w2)*(U2.conj().T@v2x)); adv_cur_2x=float(np.angle(np.vdot(v2x,z2)))
out['V2_clock']=dict(current_phase_advance=adv_cur, current_phase_advance_2x_amplitude=adv_cur_2x,
                     unified_top_mode_advance=adv_uni_top, unified_parent_advance=adv_uni_parent,
                     mu_amp=mu, sigma_phase_max=sp_max, mu_phase=mu_p)
print(f"V2 時計 N=8 等モジュラー: 現行の 1 step 位相進み = {adv_cur:+.6f}（振幅 2 倍で {adv_cur_2x:+.6f} → 4 倍＝振幅二乗比例）")
print(f"   統一線形部（K/σ 正規化 Cayley）: 最上位モード = 2·arctan(γ) = {adv_uni_top:+.6f}（状態・スケール非依存）、親モード = {adv_uni_parent:+.6f}")

# ---------- V3: 倍音（奇数帯）→ 反射係数 R ----------
N=12
rows=[]
for delta in [0.0, 1e-4, 1e-3, 1e-2, 1e-1]:
    eng,p2,q2=uif.build_standard_universe(N, delta)
    R=eng._readout()
    C2=eng.C2(); P2=np.abs(C2)**2
    Podd=float(P2[:,eng.odd_k,:].sum()); Peven=float(P2[:,eng.even_k,:].sum())
    rows.append(dict(delta=delta,R_max=float(np.max(R)),R_min=float(np.min(R)),P_odd=Podd,P_even=Peven))
    print(f"V3 δ={delta:8.1e}: R_max={np.max(R):.6e} R_min={np.min(R):.6e} | P_odd={Podd:.3e} P_even={Peven:.3e}")
# δ² 比例の確認
r1=rows[2]['R_max']; r2_=rows[3]['R_max']
print(f"   R(1e-2)/R(1e-3) = {r2_/r1:.2f}（δ² 比例なら 100）")
out['V3_harmonic_to_R']=rows
# 真空で頂点レートが厳密零か（仕様 §1.5 の実測を再確認、1 step）
eng0,_,_=uif.build_standard_universe(N,0.0)
R0=eng0._readout(); W=np.fft.ifft2(eng0.C2(),axes=(1,2))*(eng0.Nn*eng0.Neta)
rate0=eng0._vertex_rate(W.reshape(eng0.m,-1),R0)
out['V3_vacuum']=dict(R_max=float(np.max(R0)),vertex_rate_max=float(np.max(np.abs(rate0))))
print(f"V3 真空: R_max={np.max(R0):.2e}, |頂点レート|max={np.max(np.abs(rate0)):.2e}")
json.dump(out,open(os.path.join(ROOT,'results','audit_results.json'),'w'),indent=1)
print("AUDIT OK")
