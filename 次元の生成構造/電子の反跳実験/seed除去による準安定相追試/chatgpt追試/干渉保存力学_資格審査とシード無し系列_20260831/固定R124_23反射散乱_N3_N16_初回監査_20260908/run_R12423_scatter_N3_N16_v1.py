#!/usr/bin/env python3
import math, csv, json, hashlib, platform
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
STEPS=500
PHI=2*math.pi*101/124
Q=np.exp(1j*PHI)

def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(int),b.astype(int)

def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),float)
    for e in range(M):
        sh=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); sh[e]=False; A[e,sh]=1.0
    return A

def build_floor(N):
    ea,eb=edges(N); M=len(ea); d=np.minimum((eb-ea)%N,(ea-eb)%N).astype(int); D=N//2
    theta=2*math.pi*(ea+eb)/N; A=adjacency(N)
    K=A*np.sin(theta[None,:]-theta[:,None]); H=1j*K
    nd=np.array([np.sum(d==dd+1) for dd in range(D)],float)
    B=np.zeros((M,D),complex)
    for dd in range(D):
        mask=d==dd+1; B[mask,dd]=np.exp(1j*theta[mask])/math.sqrt(nd[dd])
    w,V=np.linalg.eigh(B.conj().T@H@B); chosen=None
    for k in range(len(w)):
        v=V[:,k]; j=int(np.argmax(np.abs(v))); v*=np.exp(-1j*np.angle(v[j]))
        if np.max(np.abs(v.imag))<1e-10:
            c=v.real
            if np.all(c>1e-12) or np.all(c<-1e-12): chosen=(float(w[k]),np.abs(c)); break
    if chosen is None: raise RuntimeError(f'N={N}: floor branch not found')
    lam,c=chosen; rd=c/np.sqrt(nd); z=(rd[d-1]*np.exp(1j*theta)).astype(np.complex128)
    return z,A,dict(lam=lam,floor_eigen_residual=float(np.linalg.norm(H@z-lam*z)/np.linalg.norm(z)),floor_sumz_abs=float(abs(z.sum())),floor_sumz2_abs=float(abs(z@z)))

def current_step(z,A,N):
    u=np.exp(1j*np.angle(z)); H=A*(np.conj(u)[:,None]*u[None,:]); np.fill_diagonal(H,0)
    H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    w,V=np.linalg.eigh(H); ph=np.exp(-1j*(2*math.pi/N)*w)
    return (V@(ph*(V.conj().T@z))).astype(np.complex128)

def scatter_operator(A):
    # G=sum over unordered adjacent edge-pairs P_-^(ef), P_-=1/2(e-f)(e-f)^T = 1/2 graph Laplacian
    deg=A.sum(axis=1); G=0.5*(np.diag(deg)-A)
    w,V=np.linalg.eigh(G)
    half=np.exp(0.5j*PHI*w)
    full=np.exp(1j*PHI*w)
    Uh=(V*half)@V.T; U=(V*full)@V.T
    return G,w,Uh.astype(complex),U.astype(complex)

def plane(z0):
    p=z0.real.copy(); p/=np.linalg.norm(p); q=z0.imag.copy(); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q

def metrics(z,p,q):
    h=float(np.vdot(z,z).real); zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=float(np.vdot(zp,zp).real)
    return hp/h,h,float(abs(z@z)/h),float(abs(z.sum())/math.sqrt(h))

def pair_coverage(z,tol=1e-10):
    # maximum matching of algebraic ±i candidates; local import so experiment remains reproducible if networkx exists
    import networkx as nx
    M=len(z); G=nx.Graph(); G.add_nodes_from(range(M))
    scale=max(float(np.mean(np.abs(z)**2)),1e-300)
    for i in range(M):
        for j in range(i+1,M):
            if abs(z[i]*z[i]+z[j]*z[j])/scale < tol:
                G.add_edge(i,j)
    m=nx.algorithms.matching.max_weight_matching(G,maxcardinality=True)
    return len(m),2*len(m)/M,G.number_of_edges()

summary=[]; ts=[]; spectra=[]
for N in range(3,17):
    z0,A,fc=build_floor(N); M=len(z0); p,q=plane(z0)
    G,gw,Uh,U=scatter_operator(A)
    I=np.eye(M,dtype=complex)
    order124=float(np.linalg.norm(np.linalg.matrix_power(U,124)-I,2))
    order248=float(np.linalg.norm(np.linalg.matrix_power(U,248)-I,2))
    # line-graph Laplacian spectrum / globalization diagnostics
    uniq=[]
    for x in gw:
        if not uniq or abs(x-uniq[-1][0])>1e-9: uniq.append([float(x),1])
        else: uniq[-1][1]+=1
    for val,mul in uniq: spectra.append((N,M,val,mul))
    # isolated pair exact gate
    Upair=0.5*np.array([[1+Q,1-Q],[1-Q,1+Q]],complex)
    pair124=float(np.linalg.norm(np.linalg.matrix_power(Upair,124)-np.eye(2),2))
    c0=pair_coverage(z0)
    runs={}
    for mode in ('control','scatter'):
        z=z0.copy(); arr=np.empty((STEPS+1,M),complex)
        for t in range(STEPS+1):
            arr[t]=z; hp,h,cl,cent=metrics(z,p,q); ts.append((N,mode,t,hp,h,cl,cent))
            if t<STEPS:
                if mode=='scatter': z=Uh@z
                z=current_step(z,A,N)
                if mode=='scatter': z=Uh@z
        np.savez_compressed(OUT/f'N{N:02d}_{mode}_states_500.npz',Z=arr,N=N,steps=STEPS,R=np.cos(23*math.pi/124)**2,phi=PHI)
        runs[mode]=(arr,metrics(arr[-1],p,q),pair_coverage(arr[-1]))
    sc=runs['scatter']; co=runs['control']
    summary.append((N,M,fc['floor_eigen_residual'],fc['floor_sumz_abs'],fc['floor_sumz2_abs'],pair124,order124,order248,
                    c0[0],c0[1],c0[2],co[1][0],sc[1][0],co[1][2],sc[1][2],co[2][0],co[2][1],sc[2][0],sc[2][1]))
    print('done',N,'M',M,'U124',order124,'paircov0',c0[1],'Hperp scat',sc[1][0],flush=True)

with open(OUT/'summary_N3_N16.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','M','floor_eigen_residual','floor_sumz_abs','floor_sumz2_abs','isolated_pair_U124_error','global_scatter_U124_error','global_scatter_U248_error','step0_pair_matching_pairs','step0_pair_coverage','step0_pair_candidate_edges','control_step500_Hperp_frac','scatter_step500_Hperp_frac','control_step500_closure','scatter_step500_closure','control_step500_pair_pairs','control_step500_pair_coverage','scatter_step500_pair_pairs','scatter_step500_pair_coverage']); w.writerows(summary)
with open(OUT/'timeseries_N3_N16.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','mode','step','Hperp_frac','H_total','global_closure','centroid_rel']); w.writerows(ts)
with open(OUT/'global_scatter_generator_spectrum.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','M','G_eigenvalue','multiplicity']); w.writerows(spectra)
meta={'created':'2026-09-08','N_range':[3,16],'steps':STEPS,'R':float(np.cos(23*math.pi/124)**2),'phi_over_2pi':'101/124','composition':'Uscatter_half @ Ucurrent(z) @ Uscatter_half','scatter_G':'sum P_minus over all unordered edge-pairs sharing one original vertex = 0.5*L(line graph K_N)','control_dtau':'2pi/N','python':platform.python_version(),'numpy':np.__version__}
(OUT/'RUN_METADATA.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
print('ALL DONE')
