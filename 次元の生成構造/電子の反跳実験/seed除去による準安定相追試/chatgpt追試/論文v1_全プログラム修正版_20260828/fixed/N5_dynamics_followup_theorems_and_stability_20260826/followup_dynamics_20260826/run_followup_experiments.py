#!/usr/bin/env python3
import math, json, itertools, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
GAMMA=math.tan(math.pi/144.0)  # 旧 Cayley 刻み（記録用）
ANGLE=2.0*math.pi/144.0        # FIX3: 線形回転の刻み角

# Physics core copied verbatim in structure from run_n_scaling_lowrank_v1_RAW_K.py.
def build_edges(n):
    ea, eb=np.triu_indices(n,k=1); return ea.astype(np.int64), eb.astype(np.int64)
class LowRankSystem:
    def __init__(self,n):
        self.n=n; self.ea,self.eb=build_edges(n); self.m=len(self.ea)
        self.J=np.zeros((2*n,2*n)); self.J[:n,n:]=np.eye(n); self.J[n:,:n]=-np.eye(n)
    def set_theta(self,theta): self.set_state(np.exp(1j*np.asarray(theta,dtype=float)))  # 位相のみ（親の固有モード反復でのみ使用）
    def set_state(self,z):
        # FIX4: 振幅込み生成子 K_ij = Im(conj(z_i) z_j)。c=Re z, s=Im z
        n=self.n; z=np.asarray(z,dtype=complex); self.c=np.real(z).copy(); self.s=np.imag(z).copy()
        CT=np.zeros((n,n)); ST=np.zeros((n,n)); CT[self.ea,self.eb]=self.c; CT[self.eb,self.ea]=self.c; ST[self.ea,self.eb]=self.s; ST[self.eb,self.ea]=self.s
        Gcc=CT*CT; Gcs=CT*ST; Gss=ST*ST
        np.fill_diagonal(Gcc,Gcc.sum(axis=1)); np.fill_diagonal(Gcs,Gcs.sum(axis=1)); np.fill_diagonal(Gss,Gss.sum(axis=1))
        G=np.empty((2*n,2*n)); G[:n,:n]=Gcc; G[:n,n:]=Gcs; G[n:,:n]=Gcs; G[n:,n:]=Gss; self.G=G
    def vsum(self,vals):
        n=self.n
        if np.iscomplexobj(vals):
            re=np.bincount(self.ea,vals.real,n)+np.bincount(self.eb,vals.real,n)
            im=np.bincount(self.ea,vals.imag,n)+np.bincount(self.eb,vals.imag,n); return re+1j*im
        return np.bincount(self.ea,vals,n)+np.bincount(self.eb,vals,n)
    def wt(self,z): return np.concatenate([self.vsum(self.c*z), self.vsum(self.s*z)])
    def w(self,y):
        n=self.n; yc,ys=y[:n],y[n:]
        return self.c*(yc[self.ea]+yc[self.eb])+self.s*(ys[self.ea]+ys[self.eb])
    def kmatvec(self,z):
        vs=self.vsum(self.s*z); vc=self.vsum(self.c*z)
        return self.c*(vs[self.ea]+vs[self.eb])-self.s*(vc[self.ea]+vc[self.eb])
    def sigma_spectrum(self):
        ev=np.linalg.eigvals(self.J@self.G); return np.sort(ev.imag[ev.imag>1e-12])[::-1]
    def dense_K(self): return np.column_stack([self.kmatvec(e) for e in np.eye(self.m)])
    def linear_rotation_step(self,z,sigma_unused=None):
        # FIX3: z ← exp(ANGLE·K) z（厳密な線形回転、正規化なし）
        K=self.dense_K(); w,V=np.linalg.eigh(1j*K); return V@(np.exp(-1j*ANGLE*w)*(V.conj().T@z))

def eigenmode_residual(sys,v):
    kv=sys.kmatvec(v); mu=float(np.real(np.conj(v)@(1j*kv))); return float(np.linalg.norm(1j*kv-mu*v))
def make_parent(sys,rng,iters=400,beta=0.5,tol=1e-8,restarts=3):
    best=(None,np.inf,None,0)
    for _ in range(restarts):
        theta=rng.uniform(0,2*np.pi,sys.m); v=None
        for it in range(iters):
            sys.set_theta(theta); ev,EV=np.linalg.eig(sys.J@sys.G); idx=int(np.argmin(ev.imag))
            v=sys.w(EV[:,idx].astype(complex)); theta_new=np.angle(v)  # FIX1: 振幅正規化を除去
            mix=(1-beta)*np.exp(1j*theta)+beta*np.exp(1j*theta_new); theta=np.angle(mix)
            if it%10==9:
                sys.set_theta(np.angle(v)); rr=eigenmode_residual(sys,v)
                if rr<tol: break
        sys.set_theta(np.angle(v)); rr=eigenmode_residual(sys,v)
        if rr<best[1]: best=(v.copy(),rr,sys.sigma_spectrum(),it+1)
        if rr<tol: break
    v,rr,sig,nit=best; sys.set_theta(np.angle(v)); return v,rr,sig,nit

def f_relative_to_parent(Z,v):
    p=v.real/np.linalg.norm(v.real); q=v.imag-(v.imag@p)*p; q=q/np.linalg.norm(q)
    hpar=abs(p@Z)**2+abs(q@Z)**2; htot=float(np.real(np.vdot(Z,Z))); return float(max(0,1-hpar/htot)), float(hpar), htot

def seedless_tol_sweep():
    rows=[]; series=[]
    for tol in [1e-6,1e-8,1e-10,1e-12]:
        sys=LowRankSystem(5); rng=np.random.default_rng(40260721+5000)
        v,res,sig,nit=make_parent(sys,rng,iters=2000,beta=0.5,tol=tol,restarts=3)
        Z=v.copy(); onset=None; vals=[]
        for t in range(5001):
            f,hpar,ht=f_relative_to_parent(Z,v); vals.append((t,f))
            if onset is None and f>=1e-8: onset=t
            if t<5000:
                sys.set_state(Z); Z=sys.linear_rotation_step(Z)  # FIX3/4
        arr=np.array(vals); mask=(arr[:,1]>1e-12)&(arr[:,1]<1e-4)
        rate=float(np.polyfit(arr[mask,0],np.log(arr[mask,1]),1)[0]) if mask.sum()>=5 else np.nan
        rows.append((tol,res,nit,onset,rate,arr[-1,1])); series.append(pd.DataFrame({'tol':tol,'step':arr[:,0].astype(int),'f':arr[:,1]}))
    pd.DataFrame(rows,columns=['requested_tol','parent_residual','parent_iterations','onset_f_ge_1e-8','log_f_growth_rate','f_step5000']).to_csv(HERE/'tol_sweep_summary.csv',index=False)
    pd.concat(series).to_csv(HERE/'tol_sweep_timeseries.csv',index=False)

def floquet_jacobian():
    sys=LowRankSystem(5); rng=np.random.default_rng(40260721+5000); v,res,sig,nit=make_parent(sys,rng,iters=2000,tol=1e-12,restarts=3)
    sys.set_state(v); Fv=sys.linear_rotation_step(v); phase=np.vdot(v,Fv); phase=phase/abs(phase)  # FIX3/4
    def R(z):
        sys.set_state(z); return sys.linear_rotation_step(z)/phase  # FIX3/4
    defect=np.linalg.norm(R(v)-v)
    x0=np.r_[v.real,v.imag]; m=len(v); results=[]
    for eps in [3e-6,1e-6,3e-7,1e-7]:
        J=np.zeros((2*m,2*m))
        for j in range(2*m):
            dx=np.zeros(2*m); dx[j]=eps
            zp=(x0+dx)[:m]+1j*(x0+dx)[m:]; zm=(x0-dx)[:m]+1j*(x0-dx)[m:]
            yp=R(zp); ym=R(zm); dy=(yp-ym)/(2*eps); J[:m,j]=dy.real; J[m:,j]=dy.imag
        ev=np.linalg.eigvals(J); mods=np.abs(ev); order=np.argsort(mods)[::-1]
        for rank,ii in enumerate(order): results.append((eps,rank+1,ev[ii].real,ev[ii].imag,mods[ii],math.log(mods[ii]) if mods[ii]>0 else -np.inf))
    pd.DataFrame(results,columns=['fd_eps','rank','eig_re','eig_im','modulus','log_modulus_per_step']).to_csv(HERE/'floquet_spectrum.csv',index=False)
    meta={'parent_residual':res,'parent_iterations':nit,'relative_fixed_point_defect':float(defect),'phase_re':float(phase.real),'phase_im':float(phase.imag),'sigma':sig.tolist()}
    (HERE/'floquet_meta.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')

def n5_moduli_seed_sweep():
    rows=[]
    for seed in range(20):
        sys=LowRankSystem(5); rng=np.random.default_rng(40260721+5000+seed); v,res,sig,nit=make_parent(sys,rng,iters=800,tol=1e-10,restarts=3); Z=v.copy()
        for t in range(5000): sys.set_state(Z); Z=sys.linear_rotation_step(Z)  # FIX3/4
        q=Z**2
        # Enumerate partitions into sizes 3,3,2,2 minimizing within-cluster SSE; 10 points => manageable.
        inds=set(range(10)); best=None
        for A in itertools.combinations(range(10),3):
            rem1=inds-set(A)
            for B in itertools.combinations(sorted(rem1),3):
                rem2=rem1-set(B)
                # avoid A/B permutation by mean phase lexicographic not needed; all considered
                rem2s=sorted(rem2)
                for C in itertools.combinations(rem2s,2):
                    D=tuple(sorted(rem2-set(C))); groups=[tuple(A),tuple(B),tuple(C),D]
                    centers=[q[list(g)].mean() for g in groups]
                    sse=sum(float(np.sum(np.abs(q[list(g)]-c)**2)) for g,c in zip(groups,centers))
                    if best is None or sse<best[0]: best=(sse,groups,centers)
        sse,groups,centers=best
        # Pair centers into opposite pairs minimizing |c_i+c_j|.
        pairings=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
        pairing=min(pairings,key=lambda P: sum(abs(centers[a]+centers[b]) for a,b in P))
        fam=[]
        for a,b in pairing: fam.append((centers[a]-centers[b])/2)
        rel=(np.angle(fam[1]/fam[0])+np.pi)%(2*np.pi)-np.pi
        # modulo pi because each family has sign ambiguity
        rel_mod_pi=((rel+np.pi/2)%np.pi)-np.pi/2
        rows.append((seed,res,sse,abs(fam[0]),abs(fam[1]),rel,rel_mod_pi,groups,pairing))
    df=pd.DataFrame(rows,columns=['seed','parent_residual','cluster_sse','family_A_modulus','family_B_modulus','relative_phase_rad','relative_phase_mod_pi_rad','groups','opposite_pairing'])
    df.to_csv(HERE/'N5_moduli_seed_sweep.csv',index=False)

def main():
    seedless_tol_sweep(); floquet_jacobian(); n5_moduli_seed_sweep()
if __name__=='__main__': main()
