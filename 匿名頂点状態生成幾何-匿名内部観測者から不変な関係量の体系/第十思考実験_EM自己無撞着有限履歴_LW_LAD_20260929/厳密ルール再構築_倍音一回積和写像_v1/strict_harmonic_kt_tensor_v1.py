#!/usr/bin/env python3
"""Strict harmonic-only Kastor-Traschen state generator + external exact-field readout.

Persistent state consists only of harmonic coefficients for anonymous scalar channels:
unit, scale, scale_ratio, source_strength_1, source_strength_2, source_x_1, source_x_2.
One fixed sparse quadratic interaction array performs every update synchronously.
No time t, H, exact KT state, geometry, or constraints feed the transition.
KT geometry and constraints are observation-only readouts from the generated state.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy import sparse

VARS=['unit','scale','ratio','s1','s2','x1','x2']
NV=len(VARS)

def encode(s): return np.fft.fft(s,axis=1)/s.shape[1]
def decode(c): return np.fft.ifft(c*c.shape[1],axis=1)

def idx(v,m,K): return v*K+m

def compile_map(K=4):
    # sample bilinear rules (vo,ho,v1,h1,v2,h2,c)
    last=K-1; rules=[]; unit=0
    for v in range(NV):
        for h in range(K-1): rules.append((v,h,v,h+1,unit,last,1+0j))
    # newest samples: all through same quadratic map
    # unit'=unit*unit; scale'=scale*ratio; everything else preserved via *unit
    rules.append((0,last,0,last,0,last,1+0j))
    rules.append((1,last,1,last,2,last,1+0j))
    for v in [2,3,4,5,6]: rules.append((v,last,v,last,0,last,1+0j))
    roots=np.array([1+0j,1j,-1+0j,-1j]) if K==4 else np.exp(2j*np.pi*np.arange(K)/K)
    pair_to_col={}; rows=[];cols=[];vals=[]
    for vo,ho,v1,h1,v2,h2,c0 in rules:
        for m1 in range(K):
            for m2 in range(K):
                pair=(idx(v1,m1,K),idx(v2,m2,K)); col=pair_to_col.setdefault(pair,len(pair_to_col))
                e=roots[m1]**h1 * roots[m2]**h2
                for mo in range(K):
                    c=c0*e*np.exp(-2j*np.pi*mo*ho/K)/K
                    if abs(c)>1e-15:
                        rows.append(idx(vo,mo,K));cols.append(col);vals.append(c)
    pairs=[None]*len(pair_to_col)
    for p,c in pair_to_col.items():pairs[c]=p
    pi=np.array([p[0] for p in pairs],np.int32); pj=np.array([p[1] for p in pairs],np.int32)
    M=sparse.coo_matrix((np.array(vals,complex),(rows,cols)),shape=(NV*K,len(pairs))).tocsr();M.sum_duplicates()
    M.data[np.abs(M.data)<1e-13]=0;M.eliminate_zeros()
    return M,pi,pj,len(rules)

def transition(z,M,pi,pj): return M @ (z[pi]*z[pj])

def init_state(K,H,dt,sep=8,m1=1,m2=1,a0=1):
    s=np.zeros((NV,K),complex)
    vals=[1,a0,math.exp(H*dt),m1,m2,-sep/2,sep/2]
    for v,val in enumerate(vals): s[v,:]=val
    return encode(s).reshape(-1)

def state_values(z,K):
    s=decode(z.reshape(NV,K)); return np.real(s[:,-1])

def centered_grad(f,h):
    return ((f[2:,1:-1,1:-1]-f[:-2,1:-1,1:-1])/(2*h),
            (f[1:-1,2:,1:-1]-f[1:-1,:-2,1:-1])/(2*h),
            (f[1:-1,1:-1,2:]-f[1:-1,1:-1,:-2])/(2*h))
def centered_lap(f,h):
    c=f[1:-1,1:-1,1:-1]
    return (f[2:,1:-1,1:-1]+f[:-2,1:-1,1:-1]+f[1:-1,2:,1:-1]+f[1:-1,:-2,1:-1]+f[1:-1,1:-1,2:]+f[1:-1,1:-1,:-2]-6*c)/(h*h)

def readout(z,K,N=49,L=12,dt=0.25,exclusion=1.2):
    unit,a,ratio,m1,m2,x1,x2=state_values(z,K)
    # observation-only H from state ratio
    H=math.log(ratio)/dt; Lam=3*H*H
    x=np.linspace(-L,L,N); h=x[1]-x[0]; X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
    r1=np.sqrt((X-x1)**2+Y*Y+Z*Z);r2=np.sqrt((X-x2)**2+Y*Y+Z*Z);eps=1e-15
    u1=m1/np.maximum(r1,eps);u2=m2/np.maximum(r2,eps);W=a+u1+u2
    gx,gy,gz=centered_grad(W,h);lap=centered_lap(W,h);Wi=W[1:-1,1:-1,1:-1];grad2=gx*gx+gy*gy+gz*gz
    R3=-4*lap/Wi**3+2*grad2/Wi**4;rho=grad2/(8*math.pi*Wi**4)
    Ham=R3+6*H*H-16*math.pi*rho-2*Lam;Gauss=-lap/Wi**3
    rr1=r1[1:-1,1:-1,1:-1];rr2=r2[1:-1,1:-1,1:-1];mask=(rr1>exclusion)&(rr2>exclusion)
    def stats(v):
        q=np.abs(v[mask]);return {'max_abs':float(q.max()),'rms':float(np.sqrt(np.mean(q*q))),'median_abs':float(np.median(q))}
    i=N//2;ua=float(u1[i,i,i]);ub=float(u2[i,i,i]);W0=float(W[i,i,i])
    decomp={'a2':a*a,'linear_a':2*a*ua,'linear_b':2*a*ub,'aa':ua*ua,'ab':2*ua*ub,'bb':ub*ub,
            'sum':a*a+2*a*ua+2*a*ub+ua*ua+2*ua*ub+ub*ub,'W2':W0*W0}
    return {'state':{'unit':unit,'scale':a,'ratio':ratio,'s1':m1,'s2':m2,'x1':x1,'x2':x2,'H_readout':H},
            'background_physical_separation':a*(x2-x1),'constraints':{'Hamiltonian':stats(Ham),'Gauss':stats(Gauss)},
            'midpoint':{'W':W0,'decomposition':decomp,'decomposition_error':abs(decomp['sum']-decomp['W2'])}}

def run(H=-0.05,dt=0.25,steps=32,K=4,Ns=(33,49,65)):
    M,pi,pj,nrules=compile_map(K); z=init_state(K,H,dt); checkpoints={0:z.copy()}
    for n in range(1,steps+1):
        z=transition(z,M,pi,pj)
        if n in (steps//2,steps):checkpoints[n]=z.copy()
    out={'model':'strict harmonic-only quadratic KT state generator','persistent_state':'harmonic coefficients only',
         'transition':'one fixed sparse quadratic product-sum interaction array','K':K,'dt':dt,'steps':steps,
         'interaction':{'shape':list(M.shape),'nnz':int(M.nnz),'pairs':len(pi),'sample_rules':nrules},'results':[]}
    ratio=math.exp(H*dt)
    for n,zv in checkpoints.items():
        val=state_values(zv,K); exact_scale=ratio**n
        rr={'step':n,'state_scale_error_vs_external_exact':abs(val[1]-exact_scale),'grids':{}}
        for N in Ns: rr['grids'][str(N)]=readout(zv,K,N=N,dt=dt)
        out['results'].append(rr)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);ap.add_argument('--H',type=float,default=-0.05);ap.add_argument('--steps',type=int,default=32)
    a=ap.parse_args();res=run(H=a.H,steps=a.steps);Path(a.out).write_text(json.dumps(res,indent=2,ensure_ascii=False));print(json.dumps(res,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
