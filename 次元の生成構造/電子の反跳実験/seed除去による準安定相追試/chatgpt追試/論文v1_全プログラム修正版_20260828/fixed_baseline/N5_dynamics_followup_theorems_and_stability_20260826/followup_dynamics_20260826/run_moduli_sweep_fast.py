#!/usr/bin/env python3
from pathlib import Path
import sys, json
import numpy as np, pandas as pd
from scipy.cluster.vq import kmeans2
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_followup_experiments import LowRankSystem, make_parent
rows=[]
for seed in range(8):
    sysm=LowRankSystem(5); rng=np.random.default_rng(40260721+5000+seed)
    v,res,sig,nit=make_parent(sysm,rng,iters=800,tol=1e-10,restarts=3); Z=v.copy()
    for t in range(5000):
        sysm.set_state(Z); Z=sysm.linear_rotation_step(Z)  # A4/R1
    q=Z**2
    data=np.c_[q.real,q.imag]
    # deterministic initial centers: four points with farthest-first selection
    chosen=[int(np.argmax(np.abs(q)))]
    while len(chosen)<4:
        dist=np.min([np.abs(q-q[c]) for c in chosen],axis=0)
        for c in chosen: dist[c]=-1
        chosen.append(int(np.argmax(dist)))
    centers,labels=kmeans2(data,data[chosen],iter=100,minit='matrix')
    cs=centers[:,0]+1j*centers[:,1]
    counts=np.bincount(labels,minlength=4)
    pairings=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
    pairing=min(pairings,key=lambda P: sum(abs(cs[a]+cs[b]) for a,b in P))
    fam=[]
    for a,b in pairing: fam.append((cs[a]-cs[b])/2)
    rel=np.angle(fam[1]/fam[0]); rel_mod_pi=((rel+np.pi/2)%np.pi)-np.pi/2
    opposite_error=sum(abs(cs[a]+cs[b]) for a,b in pairing)
    sse=sum(abs(q[i]-cs[labels[i]])**2 for i in range(10))
    rows.append(dict(seed=seed,parent_residual=res,parent_iterations=nit,counts=';'.join(map(str,sorted(counts,reverse=True))),cluster_sse=sse,opposite_error=opposite_error,family1_mod=abs(fam[0]),family2_mod=abs(fam[1]),relative_phase_rad=rel,relative_phase_mod_pi_rad=rel_mod_pi))
pd.DataFrame(rows).to_csv(HERE/'N5_moduli_seed_sweep.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
