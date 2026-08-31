import numpy as np, mpmath as mp, csv, json, os, time
mp.mp.dps=100
IN='/mnt/data/parent_v.npz'; OUT='/mnt/data/hm100_redo'; os.makedirs(OUT,exist_ok=True)
D=np.load(IN,allow_pickle=True); v0=np.array(D['v'],dtype=np.complex128); edges=np.array(D['edges'],dtype=int)
M=len(v0); N=6; STEPS=500

def mf(x):
    n,d=float(x).as_integer_ratio(); return mp.mpf(n)/mp.mpf(d)
def mc(z): return mp.mpc(mf(z.real),mf(z.imag))
z=mp.matrix([mc(x) for x in v0])
# Original float64 observation plane, then exact promotion.
p64=v0.real/np.linalg.norm(v0.real)
q64=v0.imag-(v0.imag@p64)*p64; q64/=np.linalg.norm(q64)
p=[mf(x) for x in p64]; q=[mf(x) for x in q64]
A=np.zeros((M,M),dtype=np.int8)
for e,(a,b) in enumerate(edges):
    for f,(c,d) in enumerate(edges):
        if e!=f and (a==c or a==d or b==c or b==d): A[e,f]=1
DELTA=2*mp.pi/N

def metrics(z):
    Htot=mp.fsum([abs(z[i])**2 for i in range(M)])
    cp=mp.fsum([p[i]*z[i] for i in range(M)]); cq=mp.fsum([q[i]*z[i] for i in range(M)])
    hp=mp.mpf('0')
    for i in range(M):
        u=z[i]-p[i]*cp-q[i]*cq; hp += abs(u)**2
    cl=abs(mp.fsum([z[i]*z[i] for i in range(M)]))/Htot
    return Htot,hp,hp/Htot,cl

def Hmat(z):
    H=mp.matrix(M,M)
    for e in range(M):
        ce=mp.conj(z[e])
        for f in range(M):
            if A[e,f]: H[e,f]=ce*z[f]
    return H

def step(z):
    H=Hmat(z)
    ev,Q=mp.eighe(H)
    c=Q.H*z
    for i in range(M): c[i] *= mp.exp(-1j*DELTA*ev[i])
    return Q*c

rows=[]; t0=time.time()
for t in range(STEPS+1):
    met=metrics(z); rows.append((t,*met))
    if t in (0,1,10,50,100,200,300,400,500):
        print(t, mp.nstr(met[2],12), 'sec',round(time.time()-t0,2), flush=True)
    if t==STEPS: break
    z=step(z)

csvp=os.path.join(OUT,'hm_N6_deltaN_100dps_preserved_float64_seed_500.csv')
with open(csvp,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['step','H_total','H_perp','Hperp_frac','closure_frac'])
    for r in rows: w.writerow([r[0]]+[mp.nstr(x,110) for x in r[1:]])
meta={'N':6,'steps':500,'mp_dps':100,'delta':'2*pi/6','initial_source':'original hm_N6 parent_v.npz complex128 values promoted exactly; original float64 numerical noise preserved','projection_basis':'p,q computed in float64 exactly as pass2_run.py then promoted exactly','interaction':'H_ef=A_ef*conj(z_e)*z_f','update':'mpmath.eighe Hermitian diagonalization; z_next=exp(-i delta H(z))z','initial_Hperp_frac':mp.nstr(rows[0][3],40),'final_Hperp_frac':mp.nstr(rows[-1][3],40),'elapsed_seconds':time.time()-t0}
with open(os.path.join(OUT,'manifest.json'),'w') as f: json.dump(meta,f,indent=2)
print(json.dumps(meta,indent=2))
