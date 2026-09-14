#!/usr/bin/env python3
import numpy as np, itertools

def elementary(z):
    m=len(z); e=np.zeros(m+1,dtype=complex); e[0]=1
    for v in z:
        for k in range(m,0,-1): e[k]+=e[k-1]*v
    return e

def J_invariants(z):
    e=elementary(z)
    if abs(e[1])<1e-12:
        raise ValueError("chart e1 != 0 required")
    return np.array([e[k]/e[1]**k for k in range(2,len(z)+1)])

def kappa(z):
    z=np.asarray(z,dtype=complex)
    return abs(np.sum(z*z))/np.sum(np.abs(z)**2)

def reconstruct(J):
    m=len(J)+1
    coeff=[1,-1]+[((-1)**k)*Jk for k,Jk in enumerate(J,start=2)]
    return np.roots(np.array(coeff,dtype=complex))

# Theory:
# Bob-star has m=N-1 unordered complex relations.
# Gauge is C* x S_m.
# On chart e1 != 0:
#   J_k=e_k/e_1^k, k=2..m
# are invariant and locally complete.
# Normalized unordered relations y_j=z_j/e1 are roots of:
#   y^m-y^(m-1)+J2 y^(m-2)-J3 y^(m-3)+...+(-1)^m Jm=0.
