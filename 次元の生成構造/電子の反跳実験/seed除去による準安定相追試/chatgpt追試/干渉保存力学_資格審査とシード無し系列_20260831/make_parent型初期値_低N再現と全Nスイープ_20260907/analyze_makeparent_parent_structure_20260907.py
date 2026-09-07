#!/usr/bin/env python3
import os,csv,math
import numpy as np
BASE=os.path.dirname(os.path.abspath(__file__)); P=os.path.join(BASE,'parents_actual_N3_N7')

def adjacency(N):
    ea,eb=np.triu_indices(N,k=1); M=len(ea); A=np.zeros((M,M))
    for e in range(M):
        sh=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); sh[e]=False; A[e,sh]=1
    return A
rows=[]
for N in range(3,8):
    z=np.asarray(np.load(os.path.join(P,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],complex)
    th=np.angle(z); r=np.abs(z); A=adjacency(N)
    W=A*np.sin(th[None,:]-th[:,None])**2
    lam=float(np.vdot(r,W@r).real/np.vdot(r,r).real)
    res=float(np.linalg.norm(W@r-lam*r)/np.linalg.norm(r))
    # global phase modulo pi/2: nearest Z4 residual
    phi=th[0]; d=((th-phi+np.pi/4)%(np.pi/2))-np.pi/4
    q=np.rint((th-phi-d)/(np.pi/2)).astype(int)%4
    p=[float(np.sum(r[q==k]**2)) for k in range(4)]
    rows.append([N,len(z),abs(z@z),abs(np.sum(z)),len(np.unique(np.round(r,12))),lam,res,max(abs(d)),*p,p[0]+p[2],p[1]+p[3]])
with open(os.path.join(BASE,'parent_structure_summary_N3_N7.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','M','abs_sum_z2','abs_sum_z','amplitude_species_1e-12','W_perron_lambda','W_perron_residual','max_Z4_phase_residual_rad','P0','P1','P2','P3','P_even','P_odd']); w.writerows(rows)
print('STRUCTURE ANALYSIS DONE')
