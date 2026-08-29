# -*- coding: utf-8 -*-
"""エネルギー・カシミール判定：親 v における H_int の球面上ヘッセ行列の慣性指数。
H_int(z) = ½ Σ_{隣接辺対} (Im z̄_i z_j)²、保存量 ‖z‖²。v は ∇H_int(v) = μ v（球面上の臨界点）。
Hess_L = Hess(H_int) − μ I を接空間 T = {δ ⊥ v} ∩ {δ ⊥ iv}（球面接空間から全体位相方向を除く）に制限し固有値の符号を数える。
 定符号（n+ = 0 または n− = 0、零モードは連続族＝moduli 分）なら軌道は非線形安定（インフレーション不在の証明）。
 不定なら線形不安定の可能性があり、流れの成長率 a（2 刻み分解）と突合する。
usage: python3 hessian_inertia_energy_casimir.py Nmin Nmax n_seeds  → data/hessian_inertia.csv"""
import os, sys, math, importlib.util, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
Nmin,Nmax,NS=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
def grad(A,x,y):
    K=A*(np.outer(x,y)-np.outer(y,x))  # K_ij = x_i y_j − y_i x_j
    return np.concatenate([(K*y[None,:]).sum(1), -(K*x[None,:]).sum(1)])
def Hint(A,x,y): K=A*(np.outer(x,y)-np.outer(y,x)); return 0.25*np.sum(K**2)  # ½Σ_{i<j} = ¼Σ_{i,j}
def lam_f(A,v,ang,h=1e-6):
    M=len(v); to_r=lambda z: np.concatenate([z.real,z.imag]); F=lambda z: eng._exp_step(eng._K_amplitude_aware(A,z),z,ang); phi=float(np.angle(np.vdot(v,F(v)))); J=np.zeros((2*M,2*M))
    for k in range(2*M):
        e=np.zeros(2*M); e[k]=h; dd=e[:M]+1j*e[M:]; J[:,k]=(to_r(F(v+dd))-to_r(F(v-dd)))/(2*h)
    c,sn=np.cos(-phi),np.sin(-phi); R=np.block([[c*np.eye(M),-sn*np.eye(M)],[sn*np.eye(M),c*np.eye(M)]]); return 2*math.log(np.abs(np.linalg.eigvals(R@J)).max())
rows=[]
for N in range(Nmin,Nmax+1):
    for seed in range(NS):
        s=eng.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N+seed)
        try: v,res,sig=eng.make_parent(s,rng)
        except RuntimeError: continue
        A=eng._adjacency(s); M=len(v); x,y=v.real.copy(),v.imag.copy(); w=np.concatenate([x,y])
        g=grad(A,x,y); mu=float(g@w/(w@w)); crit=float(np.linalg.norm(g-mu*w)/np.linalg.norm(g))
        # ヘッセ行列（勾配の中心差分）
        h=1e-5; Hs=np.zeros((2*M,2*M))
        for k in range(2*M):
            e=np.zeros(2*M); e[k]=h; Hs[:,k]=(grad(A,(w+e)[:M],(w+e)[M:])-grad(A,(w-e)[:M],(w-e)[M:]))/(2*h)
        Hs=0.5*(Hs+Hs.T); HL=Hs-mu*np.eye(2*M)
        # 接空間: v と iv を除く
        iv=np.concatenate([-y,x]); B=np.stack([w/np.linalg.norm(w),iv/np.linalg.norm(iv)],1); Q,_=np.linalg.qr(np.concatenate([B,np.eye(2*M)],1)); T=Q[:,2:2*M]
        ev=np.linalg.eigvalsh(T.T@HL@T); scale=np.abs(ev).max(); tol=1e-7*scale
        npos,nneg,nzero=int((ev>tol).sum()),int((ev<-tol).sum()),int((np.abs(ev)<=tol).sum())
        a1,a2=2*math.pi/124,2*math.pi/248; l1,l2=lam_f(A,v,a1),lam_f(A,v,a2); b=(l1-2*l2)/(a1**2-2*a2**2); a=(l1-b*a1**2)/a1
        rows.append(dict(N=N,seed=seed,M=M,parent_residual=float(res),sigma_max=float(sig[0]),mu=mu,critical_point_residual=crit,H_int=Hint(A,x,y),n_pos=npos,n_neg=nneg,n_zero=nzero,ev_min=float(ev.min()),ev_max=float(ev.max()),definite=(npos==0 or nneg==0),a_flow_per_tau=a,b_discretization=b))
        print(f"N={N} seed={seed}: crit={crit:.1e} inertia (+,−,0)=({npos},{nneg},{nzero}) ev[{ev.min():.3e},{ev.max():.3e}] definite={npos==0 or nneg==0}  a={a:.4f} b={b:.3f}",flush=True)
D=pd.DataFrame(rows); out=os.path.join(ROOT,"data",f"hessian_inertia_N{Nmin}_{Nmax}_s{NS}.csv"); D.to_csv(out,index=False); print("saved",out)
print(D.groupby("N").agg(n_parents=("seed","count"),n_definite=("definite","sum"),n_unstable_a=("a_flow_per_tau",lambda x:(x>1e-4).sum()),n_neg_min=("n_neg","min"),n_neg_max=("n_neg","max"),n_pos_min=("n_pos","min"),n_pos_max=("n_pos","max"),n_zero_max=("n_zero","max")).to_string())
