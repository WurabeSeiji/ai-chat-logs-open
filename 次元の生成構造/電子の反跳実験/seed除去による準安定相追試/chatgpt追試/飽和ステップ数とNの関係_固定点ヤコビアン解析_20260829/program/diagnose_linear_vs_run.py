# -*- coding: utf-8 -*-
"""ヤコビアン最大乗数（λ_J ~ 1e-5）と実測成長率（λ_emp ~ 1e-2）の食い違いの診断。
(a) 非線形写像 F の反復: z0 = v（丸め誤差のみ）と z0 = v + 1e-10·g（陽な摂動）から H⊥/H を追跡し成長率を測る
(b) 線形化 DF（中心差分）の反復（べき乗法）: 乱数ベクトルに DF を t 回掛けたノルム成長 → 非正規性・Jordan 構造込みの線形成長率
(c) DF の最大特異値、および DF のスペクトルの分布
usage: python3 diagnose_linear_vs_run.py N steps"""
import os, sys, json, math, importlib.util, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
N=int(sys.argv[1]); T=int(sys.argv[2]); ANGLE=eng.ANGLE
sys_lr=eng.LowRankSystem(N); rng=np.random.default_rng(40260722+1000*N); v,res,sig=eng.make_parent(sys_lr,rng); A=eng._adjacency(sys_lr); M=len(v)
F=lambda z: eng._exp_step(eng._K_amplitude_aware(A,z),z,ANGLE)
p=v.real/np.linalg.norm(v.real); q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)
def fperp(z): zp=z-p*(p@z)-q*(q@z); return float(np.vdot(zp,zp).real/np.vdot(z,z).real)
out={"N":N,"M":M,"parent_residual":res,"steps":T}
# (a)
g=rng.standard_normal(M)+1j*rng.standard_normal(M); g/=np.linalg.norm(g)
for lab,z0 in (("roundoff_only",v.copy()),("perturbed_1e-10",v+1e-10*np.linalg.norm(v)*g)):
    z=z0.copy(); fs=[]
    for t in range(T): z=F(z); fs.append(fperp(z))
    fs=np.array(fs); s=np.arange(1,T+1); m=(fs>1e-16)&(fs<1e-3)
    lam=float(np.polyfit(s[m],np.log(fs[m]),1)[0]) if m.sum()>20 else None
    out[lab]={"f_first":fs[0],"f_final":fs[-1],"lambda_fit":lam,"n_fit":int(m.sum())}; print(lab,out[lab],flush=True)
# (b),(c)
h=1e-6; J=np.zeros((2*M,2*M)); to_r=lambda z: np.concatenate([z.real,z.imag])
for k in range(2*M):
    e=np.zeros(2*M); e[k]=h; d=e[:M]+1j*e[M:]; J[:,k]=(to_r(F(v+d))-to_r(F(v-d)))/(2*h)
mu=np.linalg.eigvals(J); sv=np.linalg.svd(J,compute_uv=False)
out["DF_max_abs_eig"]=float(np.abs(mu).max()); out["DF_max_singular"]=float(sv.max()); out["DF_min_singular"]=float(sv.min())
x=rng.standard_normal(2*M); x/=np.linalg.norm(x); norms=[]
for t in range(T): x=J@x; n=np.linalg.norm(x); norms.append(n); x/=n
ln=np.cumsum(np.log(norms)); s=np.arange(1,T+1); m=s>T//2
out["power_iter_lambda_second_half"]=float(np.polyfit(s[m],ln[m],1)[0]); out["power_iter_total_ln_growth"]=float(ln[-1])
# 線形化 DF を H⊥ 成分で追跡（摂動ベクトルの成長）: 実 2M → 複素
x=to_r(g); fs=[]
for t in range(T): x=J@x; zc=x[:M]+1j*x[M:]; fs.append(float(np.linalg.norm(zc-p*(p@zc)-q*(q@zc))**2))
fs=np.array(fs); m=s>T//2; out["linear_DF_Hperp_lambda_second_half"]=float(np.polyfit(s[m],np.log(fs[m]),1)[0])
print(json.dumps({k:(v if not isinstance(v,float) else float(f"{v:.6g}")) for k,v in out.items() if not isinstance(v,dict)},ensure_ascii=False))
json.dump(out,open(os.path.join(ROOT,"results",f"diagnose_N{N}.json"),"w"),indent=1,ensure_ascii=False)
