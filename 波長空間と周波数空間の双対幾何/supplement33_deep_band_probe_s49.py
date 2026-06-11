#!/usr/bin/env python3
# s=49 奇数セクター深掘り: nc=23..31 でも開かないか
import numpy as np, itertools, math
def centers4(s):
    R=math.sqrt(s); K=int(math.floor(R))+1
    rng=np.arange(-K,K+1)
    g=np.stack(np.meshgrid(rng,rng,rng,rng,indexing='ij'),axis=-1).reshape(-1,4)
    return g[((np.abs(g)+0.5)**2).sum(axis=1)<=s+1e-9]
s=49; R=7.0; P=14.0; nc_max=31
cells=centers4(s).astype(np.float64); N=len(cells)
K=int(np.abs(cells).max())
ax=np.array(sorted([0]+[v for n in range(1,nc_max+1,2) for v in (n,-n)]))
A=len(ax)
n4=np.array(list(itertools.product(ax,repeat=4)),dtype=np.float64)
M=len(n4); S=np.zeros(M,dtype=np.complex128)
chunk=max(1,int(1e7//N)); w=-2j*np.pi/P
for i in range(0,M,chunk):
    nn=n4[i:i+chunk]; ph=cells@nn.T
    S[i:i+chunk]=np.exp(w*ph).sum(axis=0)
def phi(n):
    out=np.ones_like(n); nz=n!=0
    out[nz]=np.sin(np.pi*n[nz]/P)/(np.pi*n[nz]/P); return out
PHI=phi(n4[:,0])*phi(n4[:,1])*phi(n4[:,2])*phi(n4[:,3])
C=(S*PHI/(P**4)).reshape(A,A,A,A)
slots=np.arange(-(K+1),K+2)
E=np.exp(2j*np.pi*np.outer(slots,ax)/P)
occmask=np.zeros((len(slots),)*4,dtype=bool)
for k in cells.astype(int): occmask[tuple(t+K+1 for t in k)]=True
print(f"s=49 (N={N}, K={K}) 奇数帯域深掘り:")
for nc in range(21,nc_max+1,2):
    keep=np.abs(ax)<=nc
    Ck=C*(keep[:,None,None,None]&keep[None,:,None,None]&keep[None,None,:,None]&keep[None,None,None,:])
    f=np.einsum('abcd,pa,qb,rc,sd->pqrs',Ck,E,E,E,E,optimize=True).real
    eye=f[occmask].min()-f[~occmask].max()
    lam=2*R/nc
    print(f"  nc={nc}: eye={eye:+.4f}  (lambda={lam:.2f} cell)  {'OPEN' if eye>0 else 'closed'}")
