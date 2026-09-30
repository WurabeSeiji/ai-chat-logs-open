#!/usr/bin/env python3
"""Strict one-map harmonic spin hierarchy check.
Persistent state: harmonic coefficients only.
Single fixed sparse quadratic product-sum update.
No spin labels are used by transition; interpretations are readout only.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy import sparse

# anonymous scalar channels: unit, rhalf, rone, v, a, b
NV=6

def enc(s): return np.fft.fft(s,axis=1)/s.shape[1]
def dec(c): return np.fft.ifft(c*c.shape[1],axis=1)
def idx(v,m,K): return v*K+m

def compile_map(K=4):
    last=K-1; rules=[];unit=0
    for v in range(NV):
        for h in range(K-1): rules.append((v,h,v,h+1,unit,last,1+0j))
    # all newest values through same quadratic product-sum map
    rules += [
      (0,last,0,last,0,last,1+0j),
      (1,last,1,last,0,last,1+0j),
      (2,last,2,last,0,last,1+0j),
      (3,last,3,last,1,last,1+0j), # v *= rhalf
      (4,last,4,last,2,last,1+0j), # a *= rone
      (5,last,5,last,2,last,1+0j), # b *= rone
    ]
    roots=np.array([1+0j,1j,-1+0j,-1j]) if K==4 else np.exp(2j*np.pi*np.arange(K)/K)
    pc={};rows=[];cols=[];vals=[]
    for vo,ho,v1,h1,v2,h2,c0 in rules:
      for m1 in range(K):
       for m2 in range(K):
        pair=(idx(v1,m1,K),idx(v2,m2,K));col=pc.setdefault(pair,len(pc));e=roots[m1]**h1*roots[m2]**h2
        for mo in range(K):
         c=c0*e*np.exp(-2j*np.pi*mo*ho/K)/K
         if abs(c)>1e-15:rows.append(idx(vo,mo,K));cols.append(col);vals.append(c)
    pairs=[None]*len(pc)
    for p,c in pc.items():pairs[c]=p
    pi=np.array([p[0] for p in pairs],np.int32);pj=np.array([p[1] for p in pairs],np.int32)
    M=sparse.coo_matrix((np.array(vals,complex),(rows,cols)),shape=(NV*K,len(pairs))).tocsr();M.sum_duplicates();M.data[np.abs(M.data)<1e-13]=0;M.eliminate_zeros()
    return M,pi,pj,len(rules)

def transition(z,M,pi,pj): return M@(z[pi]*z[pj])
def init(K,N):
    s=np.zeros((NV,K),complex);rh=np.exp(1j*np.pi/N);r1=np.exp(2j*np.pi/N)
    vals=[1+0j,rh,r1,1+0j,1+0j,np.exp(1j*0.37)]
    for v,x in enumerate(vals):s[v,:]=x
    return enc(s).reshape(-1)
def newest(z,K):return dec(z.reshape(NV,K))[:,-1]

def run(N=64,K=4):
    M,pi,pj,nrules=compile_map(K);z=init(K,N);z0=newest(z,K); rec=[]
    for n in range(1,2*N+1):
      z=transition(z,M,pi,pj);x=newest(z,K)
      if n in [1,N,2*N]:
        unit,rh,r1,v,a,b=x; aa=a*a;ab=2*a*b;bb=b*b
        rec.append({'step':n,'V': [v.real,v.imag],'a':[a.real,a.imag],'b':[b.real,b.imag],
                    'aa':[aa.real,aa.imag],'ab':[ab.real,ab.imag],'bb':[bb.real,bb.imag],
                    'rhalf2_minus_rone':float(abs(rh*rh-r1))})
    x=newest(z,K)
    # independent readout tests at exactly N and 2N need rerun captures
    def at(k):
      zz=init(K,N)
      for _ in range(k):zz=transition(zz,M,pi,pj)
      return newest(zz,K)
    xN=at(N);x2=at(2*N)
    # weight-2 covariance in one step from readouts
    xbefore=at(0);xafter=at(1); aa0=xbefore[4]**2;aa1=xafter[4]**2
    expected=xbefore[2]**2*aa0
    return {'model':'strict harmonic-only spin hierarchy','persistent_state':'harmonic coefficients only','transition':'one fixed sparse quadratic interaction array',
            'N':N,'K':K,'interaction':{'shape':list(M.shape),'nnz':int(M.nnz),'pairs':len(pi),'sample_rules':nrules},
            'tests':{'half_2pi_signflip_error':float(abs(xN[3]+xbefore[3])),
                     'half_4pi_return_error':float(abs(x2[3]-xbefore[3])),
                     'weight2_one_step_covariance_error':float(abs(aa1-expected)),
                     'generator_consistency_error':float(abs(xbefore[1]**2-xbefore[2]))},'records':rec}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);ap.add_argument('--N',type=int,default=64);a=ap.parse_args();d=run(a.N);Path(a.out).write_text(json.dumps(d,indent=2));print(json.dumps(d,indent=2))
if __name__=='__main__':main()
