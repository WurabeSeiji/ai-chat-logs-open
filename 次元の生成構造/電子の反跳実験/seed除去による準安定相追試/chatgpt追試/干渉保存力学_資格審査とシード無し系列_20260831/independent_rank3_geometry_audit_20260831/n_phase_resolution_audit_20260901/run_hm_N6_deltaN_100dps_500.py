#!/usr/bin/env python3
import mpmath as mp
import csv, os, time

mp.mp.dps = 100
N=6
M=15
STEPS=500
delta=2*mp.pi/N
OUT='/mnt/data/hm100'
os.makedirs(OUT,exist_ok=True)

# Edge order is the original hm_N6 parent_v.npz order.
edges=[(i,j) for i in range(N) for j in range(i+1,N)]
# Exact phase classes inferred from the original handmade_equimodular_1factor parent:
# theta = class * pi/5 (0,36,72,108,144 deg)
cls=[3,1,4,2,0, 4,2,0,1, 0,3,2, 1,3,4]
r=mp.sqrt(mp.mpf(2))/5  # r^2 = 0.08, H_total = 15*0.08 = 1.2
z0=[r*mp.e**(1j*mp.pi*mp.mpf(c)/5) for c in cls]

# adjacency list in edge space: two K6 edges interact iff they share a vertex
nbr=[]
for e,(i,j) in enumerate(edges):
    ns=[]
    for f,(k,l) in enumerate(edges):
        if e!=f and (i==k or i==l or j==k or j==l): ns.append(f)
    nbr.append(ns)

# rank-2 real plane basis from initial Re/Im, exactly as the existing readout logic.
def dot_real_complex(a,z):
    return mp.fsum([a[i]*z[i] for i in range(M)])

def norm_real(a):
    return mp.sqrt(mp.fsum([x*x for x in a]))

pr=[mp.re(x) for x in z0]
pn=norm_real(pr); p=[x/pn for x in pr]
qi=[mp.im(x) for x in z0]
proj=mp.fsum([qi[i]*p[i] for i in range(M)])
qraw=[qi[i]-proj*p[i] for i in range(M)]
qn=norm_real(qraw); q=[x/qn for x in qraw]

def H_apply_fixed(z, x):
    # H_ef = A_ef conj(z_e) z_f ; z is frozen during one macro step
    y=[0j]*M
    for e in range(M):
        s=mp.mpc(0)
        for f in nbr[e]:
            s += z[f]*x[f]
        y[e]=mp.conj(z[e])*s
    return y

def exp_action(z):
    # exp(-i delta H(z)) z by adaptive Taylor action, all at 100 decimal digits.
    out=[mp.mpc(x) for x in z]
    term=[mp.mpc(x) for x in z]
    # tolerance safely below 100-digit observable floor; relative to vector norm.
    tol=mp.mpf('1e-108')
    for k in range(1,600):
        ht=H_apply_fixed(z,term)
        fac=(-1j*delta)/k
        term=[fac*x for x in ht]
        out=[out[i]+term[i] for i in range(M)]
        tn=mp.sqrt(mp.fsum([abs(x)**2 for x in term]))
        on=mp.sqrt(mp.fsum([abs(x)**2 for x in out]))
        if tn <= tol*max(on,mp.mpf(1)):
            return out,k
    raise RuntimeError('Taylor series did not converge')

def metrics(z):
    Htot=mp.fsum([abs(x)**2 for x in z])
    ap=dot_real_complex(p,z); aq=dot_real_complex(q,z)
    zp=[z[i]-p[i]*ap-q[i]*aq for i in range(M)]
    Hp=mp.fsum([abs(x)**2 for x in zp])
    frac=Hp/Htot
    clos=abs(mp.fsum([x*x for x in z]))/Htot
    return Htot,frac,clos

csv_path=os.path.join(OUT,'hm_N6_deltaN_100dps_500.csv')
rows=[]
z=[mp.mpc(x) for x in z0]
t0=time.time(); ks=[]
for step in range(STEPS+1):
    Htot,frac,clos=metrics(z)
    rows.append((step, mp.nstr(Htot,105), mp.nstr(frac,105), mp.nstr(mp.log10(frac),60) if frac>0 else '-inf', mp.nstr(clos,105)))
    if step in (0,1,2,5,10,20,50,100,200,300,400,500):
        print(step, 'log10Hperp=', rows[-1][3], 'closure=', mp.nstr(clos,8), flush=True)
    if step<STEPS:
        z,k=exp_action(z); ks.append(k)
print('runtime_sec',time.time()-t0,'taylor_terms min/max',min(ks),max(ks),flush=True)
with open(csv_path,'w',newline='') as f:
    w=csv.writer(f);w.writerow(['step','H_total','Hperp_frac','log10_Hperp_frac','global_closure']);w.writerows(rows)

# Save final state as decimal text for full reproducibility.
with open(os.path.join(OUT,'hm_N6_deltaN_100dps_final_state.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['edge','re','im'])
    for e,x in zip(edges,z): w.writerow([f'{e[0]}-{e[1]}',mp.nstr(mp.re(x),105),mp.nstr(mp.im(x),105)])
