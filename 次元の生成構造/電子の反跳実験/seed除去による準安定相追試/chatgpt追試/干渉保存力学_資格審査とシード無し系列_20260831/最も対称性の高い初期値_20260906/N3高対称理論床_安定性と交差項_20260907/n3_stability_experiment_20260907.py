#!/usr/bin/env python3
import csv, math, os
import numpy as np
from scipy.linalg import expm

OUT = os.path.dirname(os.path.abspath(__file__))
N = 3
DT = 2.0 * math.pi / N
STEPS = 80
EPS_LIST = [1e-1,1e-2,1e-3,1e-4]

def edge_data(N=3):
    ea, eb = np.triu_indices(N, k=1)
    M=len(ea)
    A=np.zeros((M,M),dtype=float)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e])
        share[e]=False
        A[e,share]=1.0
    return ea,eb,A
EA,EB,A=edge_data(N)

def K_of(z):
    th=np.angle(z)
    return A*np.sin(th[None,:]-th[:,None])

def one_step(z):
    return expm(DT*K_of(z)) @ z

def floor_state():
    # Same cyclic-Fourier ansatz used in theoretical-floor generator.
    theta=2.0*math.pi*(EA+EB)/N
    K=A*np.sin(theta[None,:]-theta[:,None])
    H=1j*K
    # N=3 has one distance class, hence equal amplitude.
    z=np.exp(1j*theta)/math.sqrt(3.0)
    # phase convention does not matter physically; normalize exactly.
    z=z/np.linalg.norm(z)
    lam=np.vdot(z,H@z).real/np.vdot(z,z).real
    return z.astype(np.complex128), float(lam)

Z0,LAM=floor_state()

def orbit_distance(z, ref=Z0):
    # distance modulo one global U(1) phase
    c=np.vdot(ref,z)
    ph=1.0+0j if abs(c)==0 else c/abs(c)
    return float(np.linalg.norm(z-ph*ref)/np.linalg.norm(ref))

def cross_sum(z):
    return z[0]*z[1]+z[0]*z[2]+z[1]*z[2]

def metrics(z):
    return {
        'norm2': float(np.vdot(z,z).real),
        'centroid_abs': float(abs(np.sum(z))),
        'square_closure_abs': float(abs(np.sum(z*z))),
        'cross_sum_abs': float(abs(cross_sum(z))),
        'orbit_distance': orbit_distance(z),
    }

def perturb_closure(eps):
    # Deterministic nearby perturbation: phase-kick z0[0], hold z0[1],
    # solve z2^2=-(z0^2+z1^2), choose square-root branch nearest original z0[2].
    z=Z0.copy()
    z[0]*=np.exp(1j*eps)
    target=-(z[0]*z[0]+z[1]*z[1])
    roots=[np.sqrt(target),-np.sqrt(target)]
    z[2]=min(roots,key=lambda r: abs(r-Z0[2]))
    z/=np.linalg.norm(z)
    return z

def perturb_nonclosure(eps):
    z=Z0.copy()
    z[0]*=np.exp(1j*eps)
    z/=np.linalg.norm(z)
    return z

def run_case(kind,eps):
    z=perturb_closure(eps) if kind=='closure_preserving' else perturb_nonclosure(eps)
    rows=[]
    for t in range(STEPS+1):
        m=metrics(z)
        rows.append(dict(case=kind,epsilon=eps,step=t,**m))
        if t<STEPS:
            z=one_step(z)
    return rows

allrows=[]
for kind in ('closure_preserving','nonclosure_control'):
    for eps in EPS_LIST:
        allrows.extend(run_case(kind,eps))

csv_path=os.path.join(OUT,'n3_stability_timeseries_20260907.csv')
fields=['case','epsilon','step','orbit_distance','norm2','centroid_abs','square_closure_abs','cross_sum_abs']
with open(csv_path,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=fields)
    w.writeheader(); w.writerows(allrows)

# Summary at selected steps
summary_path=os.path.join(OUT,'n3_stability_summary_20260907.csv')
sel={0,1,2,5,10,20,40,80}
with open(summary_path,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=fields)
    w.writeheader()
    for r in allrows:
        if r['step'] in sel: w.writerow(r)

# Jacobian of one-step map in R^6 at the floor state; report singular/eigen values.
def pack(z): return np.r_[z.real,z.imag]
def unpack(x): return x[:3]+1j*x[3:]
x0=pack(Z0)
h=1e-7
J=np.zeros((6,6))
for j in range(6):
    dx=np.zeros(6); dx[j]=h
    J[:,j]=(pack(one_step(unpack(x0+dx)))-pack(one_step(unpack(x0-dx))))/(2*h)
np.savetxt(os.path.join(OUT,'n3_one_step_jacobian_20260907.csv'),J,delimiter=',')
vals=np.linalg.eigvals(J)
with open(os.path.join(OUT,'n3_jacobian_eigenvalues_20260907.csv'),'w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['real','imag','abs'])
    for v in sorted(vals,key=lambda x:abs(x),reverse=True): w.writerow([v.real,v.imag,abs(v)])

base=metrics(Z0)
with open(os.path.join(OUT,'n3_base_state_20260907.csv'),'w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['edge','real','imag','abs','phase_rad'])
    for k,z in enumerate(Z0): w.writerow([f'{EA[k]}-{EB[k]}',z.real,z.imag,abs(z),np.angle(z)])

print('lambda',LAM)
print('base',base)
for eps in EPS_LIST:
    rr=[r for r in allrows if r['case']=='closure_preserving' and r['epsilon']==eps]
    print('closure',eps,'d0,d1,d10,d80',rr[0]['orbit_distance'],rr[1]['orbit_distance'],rr[10]['orbit_distance'],rr[80]['orbit_distance'], 's2max',max(r['square_closure_abs'] for r in rr))
for eps in EPS_LIST:
    rr=[r for r in allrows if r['case']=='nonclosure_control' and r['epsilon']==eps]
    print('control',eps,'d0,d1,d10,d80',rr[0]['orbit_distance'],rr[1]['orbit_distance'],rr[10]['orbit_distance'],rr[80]['orbit_distance'], 's2',rr[0]['square_closure_abs'],rr[80]['square_closure_abs'])

# Linearization restricted to tangent of norm=1 and sum(z^2)=0 (3 real constraints).
# Constraint Jacobian C has rows d(norm2), d Re(sum z^2), d Im(sum z^2).
a=Z0.real; b=Z0.imag
C=np.vstack([
    np.r_[2*a,2*b],
    np.r_[2*a,-2*b],
    np.r_[-2*b,2*a],
])
U,S,Vh=np.linalg.svd(C)
Btan=Vh[3:].T  # 6x3 orthonormal nullspace
Jtan=Btan.T@J@Btan
vals_t=np.linalg.eigvals(Jtan)
with open(os.path.join(OUT,'n3_closure_tangent_jacobian_20260907.csv'),'w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['i','j','value'])
    for i in range(3):
        for j in range(3): w.writerow([i,j,Jtan[i,j]])
with open(os.path.join(OUT,'n3_closure_tangent_eigenvalues_20260907.csv'),'w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['real','imag','abs'])
    for v in sorted(vals_t,key=lambda x:abs(x),reverse=True): w.writerow([v.real,v.imag,abs(v)])
print('closure tangent eigs',vals_t)
