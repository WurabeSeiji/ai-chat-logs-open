# -*- coding: utf-8 -*-
"""fixed/ の全エンジンについて 4 修正を検証。参照 N5_linear124_all3fix の dense_K_amplitude_aware / exp 回転と数値一致を確認。出力 results/verify_fixes_all.json"""
import os, glob, json, math, importlib.util, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(HERE,"fixed"); out={}
def load(p):
    spec=importlib.util.spec_from_file_location("m"+str(abs(hash(p))),p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
def adjacency(s):
    A=np.zeros((s.m,s.m))
    for i in range(s.m):
        share=(s.ea==s.ea[i])|(s.ea==s.eb[i])|(s.eb==s.ea[i])|(s.eb==s.eb[i]); share[i]=False; A[i,share]=1.0
    return A
files=sorted(glob.glob(os.path.join(FX,"**","run_n_scaling_lowrank_v1_*.py"),recursive=True)+glob.glob(os.path.join(FX,"**","original_RAW_K_reference.py"),recursive=True)+glob.glob(os.path.join(FX,"**","run_followup_experiments.py"),recursive=True))
for p in files:
    rel=os.path.relpath(p,FX); src=open(p,encoding="utf-8").read(); r={}
    r["static_no_cayley"]="cayley_step" not in src; r["static_linear_rotation"]="linear_rotation_step" in src and "np.exp(-1j" in src
    r["static_set_state_amplitude"]="def set_state" in src and "np.real(z)" in src
    r["static_no_parent_norm"]=("v = v / np.linalg.norm(v)" not in src) and ("v=v/np.linalg.norm(v)" not in src)
    r["static_no_seed_init"]=("Z = v + delta * g" not in src) and ("Z = Z / np.linalg.norm(Z)" not in src)
    m=load(p); m.progress=lambda *a,**k: None
    normalized="ANGLE / sigma" in src; N=5; s=m.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N)
    res=m.make_parent(s,rng,iters=400,tol=1e-8,restarts=1); v=res[0]; c=float(np.linalg.norm(v)); r["runtime_parent_norm"]=c; r["runtime_parent_norm_ne_1"]=abs(c-1)>1e-3
    z=v*np.exp(1j*rng.uniform(0,2*np.pi,s.m))*rng.uniform(0.5,1.5,s.m)  # 振幅も位相も一般の状態
    s.set_state(z); K=s.dense_K(); Kref=adjacency(s)*np.imag(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(Kref,0.0)
    r["runtime_K_vs_reference_all3fix"]=float(np.linalg.norm(K-Kref)); r["runtime_K_antisym"]=float(np.linalg.norm(K+K.T))
    s.set_state(2*z); r["runtime_K(2z)-4K(z)"]=float(np.linalg.norm(s.dense_K()-4*K)); s.set_state(z)
    sig=float(s.sigma_spectrum()[0]); U=np.column_stack([s.linear_rotation_step(np.eye(s.m)[:,j].astype(complex),sig) for j in range(s.m)])
    w,V=np.linalg.eig(U); logU=V@np.diag(np.log(w))@np.linalg.inv(V); ang=(m.ANGLE/sig) if normalized else m.ANGLE
    r["runtime_||U^H U-I||"]=float(np.linalg.norm(U.conj().T@U-np.eye(s.m))); r["runtime_||logU-angle*K||"]=float(np.linalg.norm(logU-ang*K)); r["runtime_U_is_real"]=float(np.abs(U.imag).max())
    r["angle"]=float(m.ANGLE); r["normalized_branch"]=normalized
    r["ok"]=all([r["static_no_cayley"],r["static_linear_rotation"],r["static_set_state_amplitude"],r["static_no_parent_norm"],r["static_no_seed_init"],r["runtime_parent_norm_ne_1"],r["runtime_K_vs_reference_all3fix"]<1e-12,r["runtime_K(2z)-4K(z)"]<1e-12,r["runtime_||logU-angle*K||"]<1e-9,r["runtime_U_is_real"]<1e-12])
    out[rel]=r; print(("OK  " if r["ok"] else "NG  ")+rel, {k:(f"{v_:.2e}" if isinstance(v_,float) else v_) for k,v_ in r.items() if k.startswith("runtime")})
out["all_ok"]=all(v["ok"] for k,v in out.items() if k!="all_ok"); print("ALL OK:",out["all_ok"])
json.dump(out,open(os.path.join(HERE,"results","verify_fixes_all.json"),"w"),indent=1)
