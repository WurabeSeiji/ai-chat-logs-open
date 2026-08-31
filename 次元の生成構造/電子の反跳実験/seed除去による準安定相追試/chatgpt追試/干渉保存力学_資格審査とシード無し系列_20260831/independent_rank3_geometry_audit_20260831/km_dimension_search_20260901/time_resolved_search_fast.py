import numpy as np, math, csv, json, os
from scipy.optimize import least_squares
from math import gcd
from functools import reduce

DATA={5:'/mnt/data/km_dimension_search/hm_N5_states_treatment.npz',6:'/mnt/data/km_dimension_search/hm_N6_states_treatment.npz'}
OUT='/mnt/data/km_dimension_search'

def edges(N): return [(i,j) for i in range(N) for j in range(i+1,N)]

def extract_lambda(fn,W=4096):
    Z=np.load(fn)['Z'][-W:]
    ph=np.unwrap(np.angle(Z),axis=0)
    t=np.arange(W,dtype=float); tc=t-t.mean(); den=tc@tc
    slope=(tc[:,None]*(ph-ph.mean(0))).sum(0)/den
    lam=2*np.pi/np.maximum(np.abs(slope),1e-300)
    return lam/lam.min()

def gram_from_lengths(N,E,L):
    D2=np.zeros((N,N))
    for val,(i,j) in zip(L,E): D2[i,j]=D2[j,i]=val*val
    J=np.eye(N)-np.ones((N,N))/N
    B=-0.5*J@D2@J
    ev=np.linalg.eigvalsh(B)[::-1]
    tol=max(ev[0],1.0)*1e-8
    rank=int((ev>tol).sum())
    psd=ev[-1]>=-tol
    return B,ev,rank,psd

def primitive(k):
    g=reduce(gcd,[int(x) for x in k])
    return k//g,g

def residual_X(x,N,d,E,target):
    X=x.reshape(N,d); X=X-X.mean(0)
    ds=np.array([np.linalg.norm(X[i]-X[j]) for i,j in E])
    # relative residual; center penalty unnecessary after centering
    return (ds-target)/np.maximum(target,1e-12)

def search(N,d,lam,kmax=100,starts=250,outer=15,seed=20260901):
    E=edges(N); M=len(E); rng=np.random.default_rng(seed+100*N+d)
    best=None; solutions={}
    # scales cover low to high k, log-uniform
    for s in range(starts):
        scale=np.exp(rng.uniform(np.log(1.0),np.log(kmax*0.7)))
        X=rng.normal(size=(N,d)); X-=X.mean(0)
        # normalize typical edge, then scale relative to lambdas
        ds=np.array([np.linalg.norm(X[i]-X[j]) for i,j in E]); X*=scale/max(np.median(2*ds/lam),1e-12)
        prev=None
        for it in range(outer):
            ds=np.array([np.linalg.norm(X[i]-X[j]) for i,j in E])
            q=2*ds/lam
            k=np.clip(np.rint(q),1,kmax).astype(int)
            kp,g=primitive(k)
            if kp.max()<=kmax and kp.min()>=1:
                k=kp
            target=0.5*k*lam
            res=least_squares(residual_X,X.ravel(),args=(N,d,E,target),max_nfev=250,ftol=1e-11,xtol=1e-11,gtol=1e-11)
            X=res.x.reshape(N,d); X-=X.mean(0)
            if prev is not None and np.array_equal(k,prev): break
            prev=k.copy()
        # exact check from integer-derived lengths, not optimized X
        k,_=primitive(k)
        L=.5*k*lam
        B,ev,rank,psd=gram_from_lengths(N,E,L)
        tailneg=max(0.0,-ev[-1]/max(ev[0],1e-300))
        tailrank=(ev[d]/max(ev[0],1e-300) if d < len(ev) else 0.0)
        # distance mismatch of best d-dimensional realization via positive eig reconstruction
        w,V=np.linalg.eigh(B); idx=np.argsort(w)[::-1]
        wp=np.maximum(w[idx[:d]],0); Xb=V[:,idx[:d]]*np.sqrt(wp)
        Lb=np.array([np.linalg.norm(Xb[i]-Xb[j]) for i,j in E])
        rms=float(np.sqrt(np.mean(((Lb-L)/np.maximum(L,1e-12))**2)))
        score=rms + 10*tailneg + max(0,rank-d)*1e-3
        rec=dict(N=N,d=d,k=k.copy(),ev=ev.copy(),rank=rank,psd=psd,rms=rms,score=score,kmax=int(k.max()))
        if best is None or score<best['score']: best=rec
        if psd and rank==d and rms<1e-7:
            key=tuple(k.tolist()); solutions[key]=rec
    return best,list(solutions.values())


def lam_window(fn,end,W=2048):
    Z=np.load(fn)['Z']; zz=Z[end-W+1:end+1]; ph=np.unwrap(np.angle(zz),axis=0); t=np.arange(W,dtype=float); tc=t-t.mean(); sl=(tc[:,None]*(ph-ph.mean(0))).sum(0)/(tc@tc); lam=2*np.pi/np.abs(sl); return lam/lam.min()
rows=[]
for N in [5,6]:
  fn=DATA[N]
  for end in [33856,35904,37952,40000]:
    lam=lam_window(fn,end)
    for d in range(3,N):
      st=10 if N==5 else 8
      best,sol=search(N,d,lam,kmax=100,starts=st,outer=8,seed=20260910+end)
      print('N',N,'end',end,'d',d,'sol',len(sol),'best_rank',best['rank'],'psd',best['psd'],'rms',best['rms'],'kmax',best['kmax'], flush=True)
      rows.append(dict(N=N,end=end,target_d=d,n_exact=len(sol),best_rank=best['rank'],best_psd=best['psd'],best_rms=best['rms'],best_kmax=best['kmax'],best_k=' '.join(map(str,best['k']))))
with open(os.path.join(OUT,'hm_N5_N6_time_resolved_dimension_search.csv'),'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
