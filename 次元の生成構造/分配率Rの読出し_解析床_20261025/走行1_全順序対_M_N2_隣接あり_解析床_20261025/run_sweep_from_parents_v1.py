import os, math, csv, json, platform, time, sys
import numpy as np
import matplotlib.pyplot as plt

ROOT=os.path.dirname(os.path.abspath(__file__))
# 引数: PARENT_DIR PARENT_PATTERN OUT TITLE N1 N2 ...（親 npz は parents_dense_v1.py の形式: Z0, ea, eb を持つ）
PARENT_DIR=sys.argv[1]; PARENT_PATTERN=sys.argv[2]; OUT=sys.argv[3]; TITLE=sys.argv[4]
os.makedirs(OUT,exist_ok=True)
STEPS=500; OFFSETS=(-2,-1,0,1,2)
NS=tuple(int(x) for x in sys.argv[5:])   # N は引数で与える
assert np.dtype(np.float64).itemsize==8 and np.dtype(np.complex128).itemsize==16

def adjacency(ea,eb):   # 辺集合は親 npz の ea, eb から読む（構造を親と走行で共有する）。規則は原本と同一
    M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A
def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
def one_step(z,A,den):
    # 段3の最小変更（唯一の力学変更点）: 位相のみ生成子 Ĥ の虚部だけを取る H=i·K（K=sin(Δθ) 実反対称）。
    # exp(-iΔτ·iK)=exp(Δτ·K) の実直交回転となり、Z^T Z（零閉塞）と ‖Z‖ を厳密保存する。
    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
def plane(v):
    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p); q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q
def metrics(z,p,q):
    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real; return float(hp/h),float(h),float(abs(z@z)/h)
rows=[]; summaries=[]
for N in NS:   # 案B: range(3,41) → NS
    # 初期データ: 各 N の静的親ファイルの Z0 を使用
    _d=np.load(os.path.join(PARENT_DIR,PARENT_PATTERN.format(N=N)))
    z0=np.array(_d['Z0'],dtype=np.complex128,copy=True); ea,eb=_d['ea'],_d['eb']; assert len(ea)==z0.size
    A=adjacency(ea,eb); p,q=plane(z0)
    pairs=[(N,'N'),(124,'124')]   # 案B: 6 系列 → den=N と 124 の 2 系列
    for den,label in pairs:
        t_start=time.time(); z=z0.copy(); vals=np.empty(STEPS+1,np.float64); states=np.empty((STEPS+1,z.size),np.complex128); closures=np.empty(STEPS+1,np.float64); htot=np.empty(STEPS+1,np.float64)
        for t in range(STEPS+1):
            states[t]=z; vals[t],htot[t],closures[t]=metrics(z,p,q)
            if t<STEPS: z=one_step(z,A,den)
        np.savez_compressed(os.path.join(OUT,f'hm_N{N}_den_{den}_states_500.npz'),Z=states,N=np.int64(N),denominator=np.int64(den),steps=np.int64(STEPS))
        rows.extend((N,label,den,t,vals[t],htot[t],closures[t]) for t in range(STEPS+1)); ix=np.flatnonzero(vals>0.05)
        summaries.append((N,label,den,int(ix[0]) if ix.size else -1,float(vals[0]),float(vals[1]),float(vals[-1]),float(vals.max())))
        print(f'  N={N} den={den}: {time.time()-t_start:.1f} s',flush=True)   # 案B: 組ごとの経過時間
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
fig.suptitle(TITLE,y=.995)
with open(os.path.join(OUT,'RUN_METADATA_N3_N40_stage123.json'),'w') as f: json.dump({'dtype_state':'complex128','dtype_real':'float64','steps':STEPS,'N_list':list(NS),'denominators':'N,124','parent_dir':PARENT_DIR,'parent_pattern':PARENT_PATTERN,'title':TITLE,'numpy':np.__version__,'python':platform.python_version(),'input':'per-N static make_parent parents/parent_static_N*_makeparent_20260905.npz Z0 (N=40 bit-identical to canonical)','generator_normalization':'stage1+2+3: phase-only imaginary part H = 1j*Im(H_of(exp(1j*angle(z)))) = i*K, map = exp(dtau*K) real orthogonal, fixed dtau=2pi/den'},f,indent=2)
print('ALL DONE')
