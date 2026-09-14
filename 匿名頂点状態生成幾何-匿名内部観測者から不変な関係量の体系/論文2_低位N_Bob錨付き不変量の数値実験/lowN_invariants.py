#!/usr/bin/env python3
import numpy as np, itertools

def elementary(z):
    m=len(z)
    e=np.zeros(m+1,dtype=complex)
    e[0]=1
    for v in z:
        for k in range(m,0,-1):
            e[k]+=e[k-1]*v
    return e

def J_invariants(z):
    e=elementary(z)
    if abs(e[1])<1e-12:
        raise ValueError("chart e1 != 0 required")
    return np.array([e[k]/e[1]**k for k in range(2,len(z)+1)])

def kappa(z):
    z=np.asarray(z,dtype=complex)
    return abs(np.sum(z*z))/np.sum(np.abs(z)**2)
