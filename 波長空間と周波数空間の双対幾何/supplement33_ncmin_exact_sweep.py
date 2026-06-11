#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺33 v0.2: 最小証明帯域 nc_min(s) のスイープ (厳密連続フーリエ、格子ジッタなし)
# c_n = (1/P^4) * prod_i sinc(n_i/P) * S(n),  S(n) = sum_cells exp(-2пi n.k/P)
import numpy as np, itertools, math, time

OUT=[]
def log(s):
    print(s, flush=True); OUT.append(s)

def centers4(s):
    R=math.sqrt(s); K=int(math.floor(R))+1
    rng=np.arange(-K,K+1)
    g=np.stack(np.meshgrid(rng,rng,rng,rng,indexing='ij'),axis=-1).reshape(-1,4)
    val=((np.abs(g)+0.5)**2).sum(axis=1)
    return g[val<=s+1e-9]

def band_axis(nc):
    vals=[0]
    for n in range(1,nc+1,2):
        vals+= [n,-n]
    return np.array(sorted(vals))

def eye_exact(s, nc_max=9):
    R=math.sqrt(s); P=2*R
    cells=centers4(s).astype(np.float64); N=len(cells)
    K=int(np.abs(cells).max())
    ax=band_axis(nc_max); A=len(ax)
    n4=np.array(list(itertools.product(ax,repeat=4)),dtype=np.float64)
    M=len(n4)
    S=np.zeros(M,dtype=np.complex128)
    chunk=max(1, int(4e7//max(N,1)))
    w=-2j*np.pi/P
    for i in range(0,M,chunk):
        nn=n4[i:i+chunk]
        ph=cells@nn.T
        S[i:i+chunk]=np.exp(w*ph).sum(axis=0)
    def phi(n):
        out=np.ones_like(n)
        nz=n!=0
        out[nz]=np.sin(np.pi*n[nz]/P)/(np.pi*n[nz]/P)
        return out
    PHI=phi(n4[:,0])*phi(n4[:,1])*phi(n4[:,2])*phi(n4[:,3])
    C=(S*PHI/(P**4)).reshape(A,A,A,A)
    slots=np.arange(-(K+1),K+2)
    E=np.exp(2j*np.pi*np.outer(slots,ax)/P)
    Ls=len(slots)
    occmask=np.zeros((Ls,)*4,dtype=bool)
    for k in cells.astype(int):
        occmask[tuple(t+K+1 for t in k)]=True
    res={}
    for nc in range(1,nc_max+1,2):
        keep=np.abs(ax)<=nc
        Ck=C*(keep[:,None,None,None]&keep[None,:,None,None]&keep[None,None,:,None]&keep[None,None,None,:])
        f=np.einsum('abcd,pa,qb,rc,sd->pqrs',Ck,E,E,E,E,optimize=True).real
        vo=f[occmask].min(); ve=f[~occmask].max()
        res[nc]=vo-ve
    return dict(s=s,R=R,N=N,res=res)

log("最小証明帯域 nc_min(s): 厳密連続フーリエによるアイ開口 (奇数帯域 DC+|n|奇数<=nc)")
log("  s    R      セル数   eye(1)    eye(3)    eye(5)    eye(7)    eye(9)   nc_min")
ladder=[1,3,5,7,9,11,13,17,21,25,33,49,81,121]
summary={}
for s in ladder:
    t0=time.time()
    r=eye_exact(s, nc_max=9)
    e=[r['res'][nc] for nc in (1,3,5,7,9)]
    ncmin=next((nc for nc in (1,3,5,7,9) if r['res'][nc]>0), None)
    summary[s]=(r,ncmin)
    nm=str(ncmin) if ncmin else ">9"
    log(f"  {s:3d} {r['R']:6.3f} {r['N']:7d}  "+"  ".join(f"{v:+.4f}" for v in e)+f"   {nm}   [{time.time()-t0:.0f}s]")

# s=7 など >9 のケースを拡張帯域で再判定
for s,(r,ncmin) in summary.items():
    if ncmin is None:
        r2=eye_exact(s, nc_max=21)
        e2={nc:r2['res'][nc] for nc in range(1,22,2)}
        nm=next((nc for nc in range(1,22,2) if e2[nc]>0), None)
        log(f"  [拡張] s={s}: nc=11..21 eye = "+", ".join(f"{e2[nc]:+.4f}" for nc in range(11,22,2))
            +f"  -> nc_min={'>21' if nm is None else nm}")

with open('/tmp/ncmin_results.txt','w') as f: f.write("\n".join(OUT))
log("done")
