# -*- coding: utf-8 -*-
"""fixed/ の全エンジン複製と followup エンジンについて、判断 1〜6 の実装を静的・動的に検証。出力 results/verify_fixes_all.json"""
import os, glob, json, math, importlib.util, numpy as np
H=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(H,"fixed"); out={}
def load(p):
    spec=importlib.util.spec_from_file_location("eng_"+str(abs(hash(p))),p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); m.progress=lambda s:None; return m
def adjacency(s):
    A=np.zeros((s.m,s.m))
    for i in range(s.m):
        sh=(s.ea==s.ea[i])|(s.ea==s.eb[i])|(s.eb==s.ea[i])|(s.eb==s.eb[i]); sh[i]=False; A[i,sh]=1
    return A
engines=sorted(p for p in glob.glob(os.path.join(FX,"**","run_n_scaling_lowrank_v1_*.py"),recursive=True)+glob.glob(os.path.join(FX,"**","original_RAW_K_reference.py"),recursive=True) if "NORMALIZED_ORIGINAL" not in p and "ORIGINAL.py" not in p and "make_parent_path" not in p)
for p in engines:
    src=open(p,encoding="utf-8").read(); r={}
    r["static_A1_no_parent_norm"]="v = v / np.linalg.norm(v)" not in src
    r["static_A2A3_onset_seedless"]=("Z = v + delta * g" not in src) and ("Z = Z / np.linalg.norm(Z)" not in src)
    r["static_A5_generic_no_norm"]="return Z / np.linalg.norm(Z)" not in src
    r["static_R1_no_cayley_call"]=".cayley_step(" not in src and "def linear_rotation_step" in src
    r["static_S1_no_seed_call"]="zero_closure_kernel_seed(sys_lr, rng)" not in src.replace("def zero_closure_kernel_seed(sys_lr, rng)","")
    r["static_S2_no_power_iteration_call"]="sigma_max_power(wp)" not in src
    for mode in ("amplitude","phase"):
        os.environ["KMODE"]=mode; e=load(p); s=e.LowRankSystem(5); rng=np.random.default_rng(40265721)
        v,res,sig=e.make_parent(s,rng); c=float(np.linalg.norm(v)); A=adjacency(s)
        s.set_state(v); K=s.dense_K()
        zz=v/np.abs(v) if mode=="phase" else v
        Kref=A*np.imag(np.conj(zz)[:,None]*zz[None,:])
        s.set_state(2*v); K2=s.dense_K()
        U=np.column_stack([s.linear_rotation_step(np.eye(s.m)[:,j].astype(complex)) for j in range(s.m)])  # state set to 2v: exp(ANGLE K(2v))
        w,V=np.linalg.eig(U); logU=V@np.diag(np.log(w))@np.linalg.inv(V)
        rr={"parent_norm_c":c,"A1_ok":abs(c-1)>1e-3,"A4_K_matches_Im(conj z z)":float(np.linalg.norm(K-Kref)),
            "K(2z)_vs_K(z)":("4K" if np.linalg.norm(K2-4*K)<1e-12 else ("K" if np.linalg.norm(K2-K)<1e-12 else "other")),
            "R1_||logU-ANGLE*K(2v)||":float(np.linalg.norm(logU-e.ANGLE*K2)),"R1_unitarity":float(np.linalg.norm(U.conj().T@U-np.eye(s.m))),
            "A6_sigma1_actual_generator":float(np.max(np.linalg.eigvalsh(1j*K))),"validate":e.validate_against_dense(5,0,steps=100)}
        rr["mode_ok"]=rr["A4_K_matches_Im(conj z z)"]<1e-12 and rr["R1_||logU-ANGLE*K(2v)||"]<1e-9 and (rr["K(2z)_vs_K(z)"]==("4K" if mode=="amplitude" else "K")) and rr["validate"]["max_traj_dev"]<1e-9
        r[mode]=rr
    r["all_ok"]=all(v for k,v in r.items() if k.startswith("static_")) and r["amplitude"]["mode_ok"] and r["phase"]["mode_ok"] and r["amplitude"]["A1_ok"]
    out[os.path.relpath(p,FX)]=r; print(("OK  " if r["all_ok"] else "NG  ")+os.path.relpath(p,FX), "c=%.4f"%r["amplitude"]["parent_norm_c"], "σ1(amp)=%.4f σ1(phase)=%.4f"%(r["amplitude"]["A6_sigma1_actual_generator"],r["phase"]["A6_sigma1_actual_generator"]), "traj_dev=%.1e"%r["amplitude"]["validate"]["max_traj_dev"])
# followup 自前エンジン
fp=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","run_followup_experiments.py")
src=open(fp,encoding="utf-8").read(); r={"static_A1":"v=v/np.linalg.norm(v)" not in src,"static_R1":".cayley_step(" not in src and "def linear_rotation_step" in src}
for mode in ("amplitude","phase"):
    os.environ["KMODE"]=mode; e=load(fp); s=e.LowRankSystem(5); rng=np.random.default_rng(40260721+5000); v,res,sig,nit=e.make_parent(s,rng,iters=800,tol=1e-10,restarts=3)
    s.set_state(v); K=s.dense_K(); A=adjacency(s); zz=v/np.abs(v) if mode=="phase" else v; Kref=A*np.imag(np.conj(zz)[:,None]*zz[None,:])
    r[mode]={"parent_norm_c":float(np.linalg.norm(v)),"A4_K_matches":float(np.linalg.norm(K-Kref)),"sigma1":float(np.max(np.linalg.eigvalsh(1j*K)))}
r["all_ok"]=r["static_A1"] and r["static_R1"] and r["amplitude"]["A4_K_matches"]<1e-12 and r["phase"]["A4_K_matches"]<1e-12
out["followup/run_followup_experiments.py"]=r; print(("OK  " if r["all_ok"] else "NG  ")+"followup engine", r)
json.dump(out,open(os.path.join(H,"results","verify_fixes_all.json"),"w"),indent=1); print("ALL OK:", all(v["all_ok"] for v in out.values()))
