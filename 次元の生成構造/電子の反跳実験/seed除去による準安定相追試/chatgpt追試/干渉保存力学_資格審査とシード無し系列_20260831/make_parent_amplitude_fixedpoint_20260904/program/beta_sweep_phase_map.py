import sys, numpy as np, csv
from pathlib import Path
ROOT=Path('/mnt/data/make_parent_amplitude_fixedpoint_20260904')
sys.path.insert(0,str(ROOT/'runtime'))
import original_engine as eng
OUT=ROOT/'results'/'beta_sweep_phase_map.csv'

def proj_spread(x):
    # x phase differences; projective mod pi represented by exp(2ix)
    z=np.exp(2j*x); m=np.mean(z)
    return float(np.max(np.abs(z-m)))

def run(N,seed,beta,iters=1200):
    s=eng.LowRankSystem(N); rng=np.random.default_rng(seed); th=rng.uniform(0,2*np.pi,s.m)
    history=[]
    for it in range(iters):
        s.set_theta(th); ev,Y=np.linalg.eig(s.J@s.G); k=int(np.argmin(ev.imag)); sig=float(-ev[k].imag); v=s.w(Y[:,k]); ph=np.angle(v)
        mix=(1-beta)*np.exp(1j*th)+beta*np.exp(1j*ph); thn=np.angle(mix)
        if it>=iters-4: history.append((th.copy(),ph.copy(),sig,v.copy()))
        th=thn
    # recompute states t0->t1->t2 using exact beta
    def F(theta):
        s.set_theta(theta); ev,Y=np.linalg.eig(s.J@s.G); k=int(np.argmin(ev.imag)); return np.angle(s.w(Y[:,k])), float(-ev[k].imag)
    t0=th.copy(); p0,sg0=F(t0); t1=np.angle((1-beta)*np.exp(1j*t0)+beta*np.exp(1j*p0)); p1,sg1=F(t1); t2=np.angle((1-beta)*np.exp(1j*t1)+beta*np.exp(1j*p1));
    cyc=proj_spread(t2-t0); pout=proj_spread(p1-p0)
    delta=p0-t0; cube=float(np.max(np.abs(np.exp(3j*delta)-np.mean(np.exp(3j*delta)))))
    # K ratio best scalar p minimizing ||Kp-aKt||
    s.set_theta(t0); K0=np.zeros((s.m,s.m));
    # construct via action on identity
    I=np.eye(s.m); K0=np.column_stack([s.kmatvec(I[:,j]) for j in range(s.m)])
    s.set_theta(p0); Kp=np.column_stack([s.kmatvec(I[:,j]) for j in range(s.m)])
    a=float(np.vdot(K0,Kp).real/np.vdot(K0,K0).real)
    kerr=float(np.linalg.norm(Kp-a*K0)/np.linalg.norm(Kp))
    return [N,seed,beta,sg0,sg1,cyc,pout,cube,a,kerr]

rows=[]
for beta in [0.25,1/3,0.4,0.45,0.5,0.55,0.6,2/3,0.75]:
  for N in range(3,9):
    for j in range(3):
      rows.append(run(N,40260721+1000*N+j,beta))
with open(OUT,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','seed','beta','sigma0','sigma1','proj_cycle_spread','output_projective_spread','cube_spread','K_best_scalar','K_relerr']); w.writerows(rows)
print(OUT)
for beta in sorted(set(r[2] for r in rows)):
    rr=[r for r in rows if r[2]==beta]
    print(beta,'median cycle',np.median([x[5] for x in rr]),'median cube',np.median([x[7] for x in rr]),'median Kscalar',np.median([x[8] for x in rr]),'median Kerr',np.median([x[9] for x in rr]))
