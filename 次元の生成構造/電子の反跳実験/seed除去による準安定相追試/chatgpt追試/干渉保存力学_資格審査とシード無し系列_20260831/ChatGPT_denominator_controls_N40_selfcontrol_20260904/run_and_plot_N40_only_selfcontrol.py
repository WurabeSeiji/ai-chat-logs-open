import os, math, csv, json, platform
import numpy as np
import matplotlib.pyplot as plt

ROOT='/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831'
OLD_IN=os.path.join(ROOT,'data')
NEW_IN=os.path.join(ROOT,'hm_mp_free_N3_N40_20260901','data')
OUT=os.path.join(ROOT,'ChatGPT_denominator_controls_N40_selfcontrol_20260904','results')
STEPS=500; OFFSETS=(-2,-1,0,1,2)
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
rows=[]; summaries=[]
for N in range(40,41):
    if N<=16:
        Zsrc=np.load(os.path.join(OLD_IN,f'hm_N{N}','states_treatment.npz'))['Z']; z0=np.array(Zsrc[0],dtype=np.complex128,copy=True)
    else:
        vsrc=np.load(os.path.join(NEW_IN,f'hm_N{N}','parent_v.npz'))['v']; z0=np.array(vsrc,dtype=np.complex128,copy=True)
    A=adjacency(N); p,q=plane(z0)
    pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]
    for den,label in pairs:
        z=z0.copy(); vals=np.empty(STEPS+1,np.float64); states=np.empty((STEPS+1,z.size),np.complex128); closures=np.empty(STEPS+1,np.float64); htot=np.empty(STEPS+1,np.float64)
        for t in range(STEPS+1):
            states[t]=z; vals[t],htot[t],closures[t]=metrics(z,p,q)
            if t<STEPS: z=one_step(z,A,den)
        np.savez_compressed(os.path.join(OUT,f'hm_N{N}_den_{den}_states_500.npz'),Z=states,N=np.int64(N),denominator=np.int64(den),steps=np.int64(STEPS))
        rows.extend((N,label,den,t,vals[t],htot[t],closures[t]) for t in range(STEPS+1)); ix=np.flatnonzero(vals>0.05)
        summaries.append((N,label,den,int(ix[0]) if ix.size else -1,float(vals[0]),float(vals[1]),float(vals[-1]),float(vals.max())))
    print('done N',N,flush=True)
with open(os.path.join(OUT,'timeseries_64bit_with124_N3_N40.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','series','denominator','step','Hperp_frac','H_total','global_closure']); w.writerows(rows)
with open(os.path.join(OUT,'summary_64bit_with124_N3_N40.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','series','denominator','onset_gt_0.05','initial','step1','final','max']); w.writerows(summaries)
by={}
for N,label,den,t,h,hTot,cl in rows: by.setdefault((N,label),([],[]))[0].append(t); by[(N,label)][1].append(h)
fig,axs=plt.subplots(8,5,figsize=(20,24)); axs=axs.ravel(); order=['N-2','N-1','N','N+1','N+2','124']
for k,N in enumerate(range(3,41)):
    ax=axs[k]
    for label in order:
        if (N,label) in by:
            x,y=by[(N,label)]; ax.semilogy(x,y,label=('2pi/124' if label=='124' else label),linewidth=1.0)
    ax.set_xlim(0,500); ax.set_ylim(1e-34,3); ax.set_title(f'N={N}'); ax.grid(alpha=.25)
    if k//5==7: ax.set_xlabel('step')
    if k%5==0: ax.set_ylabel('Hperp/H')
    ax.legend(fontsize=7,loc='lower right')
for k in range(38,40): axs[k].axis('off')
fig.suptitle('Hperp/H denominator control (float64/complex128): 2pi/(N-2), 2pi/(N-1), 2pi/N, 2pi/(N+1), 2pi/(N+2), 2pi/124; N=3..40',y=.995); fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_Hperp_denominator_controls_with_124_N3_N40.png'),dpi=180); plt.close(fig)
with open(os.path.join(OUT,'RUN_METADATA_N3_N40.json'),'w') as f: json.dump({'dtype_state':'complex128','dtype_real':'float64','steps':STEPS,'N_range':[3,40],'denominators':'N-2,N-1,N,N+1,N+2,124','numpy':np.__version__,'python':platform.python_version(),'input_N3_N16':'original saved data/hm_N*/states_treatment.npz Z[0]','input_N17_N40':'hm_mp_free_N3_N40_20260901/data/hm_N*/parent_v.npz v'},f,indent=2)
print('ALL DONE')
