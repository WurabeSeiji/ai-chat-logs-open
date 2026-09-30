#!/usr/bin/env python3
import json, hashlib
from pathlib import Path
import numpy as np
from scipy import sparse
import strict_harmonic_rn_tensor_v1 as rn
import strict_harmonic_kt_tensor_v1 as kt
import strict_harmonic_spin_tensor_v1 as sp

def hash_sparse(M,pi,pj):
    h=hashlib.sha256()
    for a in [M.data,M.indices,M.indptr,pi,pj]: h.update(np.ascontiguousarray(a).view(np.uint8))
    return h.hexdigest()

def save_sparse(path,M,pi,pj):
    np.savez_compressed(path,data=M.data,indices=M.indices,indptr=M.indptr,shape=np.array(M.shape),pi=pi,pj=pj)

def direct_rn_reference(q=0.6,channel=0,nx=31,K=4,rmin=2.4,rmax=32.4,dt=0.08,steps=120,x0=10,sigma=3,amp=1e-4):
    r=np.linspace(rmin,rmax,nx);L=rn.Layout(nx,K);A=rn.radial_operator_bases(r,dt)
    # sample histories directly, external validator only
    z=rn.init_state(L,r,q,x0,sigma,amp,channel); sh=rn.decode_state(z,L)
    qpow=np.array([1,q,q*q,q*q*q,q*q*q*q])
    M,pi,pj,_=rn.compile_sparse_quadratic(L,r,dt); zh=z.copy()
    for _ in range(steps):
        # direct sample-space recurrence
        old=sh.copy(); new=np.empty_like(old); new[:,:-1]=old[:,1:]
        # preserve q power channels
        for p in range(5): new[L.q_var(p),-1]=old[L.q_var(p),-1]
        for co in (0,1):
            for i in range(nx):
                y=0j
                if 0<i<nx-1: y-=old[L.field_var(co,i),K-2]
                for p in range(5):
                    for ci in (0,1):
                        row=A[p,co,ci,i]
                        y += qpow[p]*np.dot(row,[old[L.field_var(ci,j),K-1] for j in range(nx)])
                new[L.field_var(co,i),-1]=y
        sh=new
        zh=rn.transition(zh,M,pi,pj)
    dh=rn.decode_state(zh,L)
    return float(np.max(np.abs(dh-sh))),hash_sparse(M,pi,pj)

def main():
    out={}
    err,h=direct_rn_reference(); out['rn_harmonic_vs_direct_sample_max_error']=err;out['rn_interaction_sha256']=h
    r=np.linspace(2.4,62.4,61);L=rn.Layout(61,4);M,pi,pj,_=rn.compile_sparse_quadratic(L,r,0.12);save_sparse('/mnt/data/rebuild_harmonic_strict_v1/strict_rn_interaction_v1.npz',M,pi,pj)
    out['rn_final_interaction_sha256']=hash_sparse(M,pi,pj)
    M2,pi2,pj2,_=kt.compile_map(4);save_sparse('/mnt/data/rebuild_harmonic_strict_v1/strict_kt_interaction_v1.npz',M2,pi2,pj2);out['kt_interaction_sha256']=hash_sparse(M2,pi2,pj2)
    M3,pi3,pj3,_=sp.compile_map(4);save_sparse('/mnt/data/rebuild_harmonic_strict_v1/strict_spin_interaction_v1.npz',M3,pi3,pj3);out['spin_interaction_sha256']=hash_sparse(M3,pi3,pj3)
    Path('/mnt/data/rebuild_harmonic_strict_v1/strict_map_validation_v1.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
