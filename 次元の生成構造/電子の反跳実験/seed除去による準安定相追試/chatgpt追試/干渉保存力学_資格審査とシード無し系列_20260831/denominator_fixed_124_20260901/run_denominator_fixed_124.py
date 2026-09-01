import os, math, csv, json, platform
import numpy as np
import matplotlib.pyplot as plt
IN=os.path.join(os.path.dirname(__file__), '..', 'data')
OUT=os.path.dirname(__file__)
STEPS=500
FIXED_DEN=124
assert np.dtype(np.float64).itemsize==8 and np.dtype(np.complex128).itemsize==16

def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A
def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
def one_step(z,A,den):
    H=H_of(z,A); w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
def plane(v):
    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p); q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q
def metrics(z,p,q):
    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real; return float(hp/h),float(h),float(abs(z@z)/h)

N_RANGE=list(range(3, 17))
rows=[]; summaries=[]
for N in N_RANGE:
    TAG=f'hm_N{N}'
    data_path=os.path.join(IN, TAG)
    if not os.path.exists(data_path):
        print(f'WARNING: {data_path} not found, skipping N={N}')
        continue
    npz_file=os.path.join(data_path, 'states_treatment.npz')
    if not os.path.exists(npz_file):
        print(f'WARNING: {npz_file} not found, skipping N={N}')
        continue
    Zsrc=np.load(npz_file)['Z']
    z0=np.array(Zsrc[0],dtype=np.complex128,copy=True)
    A=adjacency(N); p,q=plane(z0)

    z=z0.copy(); vals=np.empty(STEPS+1,np.float64); states=np.empty((STEPS+1,z.size),np.complex128); closures=np.empty(STEPS+1,np.float64); htot=np.empty(STEPS+1,np.float64)
    for t in range(STEPS+1):
        states[t]=z; vals[t],htot[t],closures[t]=metrics(z,p,q)
        if t<STEPS: z=one_step(z,A,FIXED_DEN)
    np.savez_compressed(os.path.join(OUT,f'hm_N{N}_den_124_states_500.npz'),Z=states,N=np.int64(N),denominator=np.int64(FIXED_DEN),steps=np.int64(STEPS))
    rows.extend((N,FIXED_DEN,t,vals[t],htot[t],closures[t]) for t in range(STEPS+1))
    ix=np.flatnonzero(vals>0.05)
    summaries.append((N,FIXED_DEN,int(ix[0]) if ix.size else -1,float(vals[0]),float(vals[1]),float(vals[-1]),float(vals.max())))
    print('done N',N,flush=True)

with open(os.path.join(OUT,'timeseries_fixed_124.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','denominator','step','Hperp_frac','H_total','global_closure']); w.writerows(rows)
with open(os.path.join(OUT,'summary_fixed_124.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','denominator','onset_gt_0.05','initial','step1','final','max']); w.writerows(summaries)

by={}
for N,den,t,h,hTot,cl in rows: by.setdefault(N,([],[]))[0].append(t); by[N][1].append(h)
fig,axs=plt.subplots(4,4,figsize=(16,12)); axs=axs.ravel()
for k,N in enumerate(range(3,17)):
    ax=axs[k]
    if N in by:
        x,y=by[N]; ax.semilogy(x,y,label=f'Δτ=2π/124',linewidth=2.0,color='red')
    ax.set_xlim(0,500); ax.set_ylim(1e-34,3); ax.set_title(f'N={N}'); ax.grid(alpha=.25)
    if k//4==3: ax.set_xlabel('step')
    if k%4==0: ax.set_ylabel('Hperp/H')
    ax.legend(fontsize=8,loc='lower right')
for k in range(14,16): axs[k].axis('off')
fig.suptitle('Hperp/H: Fixed Δτ=2π/124 across all N (N=3..16)',y=.995); fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_Hperp_fixed_denominator_124.png'),dpi=180); plt.close(fig)

with open(os.path.join(OUT,'RUN_METADATA.json'),'w') as f:
    json.dump({'dtype_state':'complex128','dtype_real':'float64','steps':STEPS,'N_range':[3,16],'denominator_fixed':FIXED_DEN,'description':'Single-denominator run with Δτ=2π/124 (fixed across all N)','numpy':np.__version__,'python':platform.python_version(),'input':'original saved states_treatment.npz Z[0] for each hm_N'},f,indent=2)
print('ALL DONE (fixed denominator 124)')
