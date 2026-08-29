#!/usr/bin/env python3
"""手作り自己無撞着親（3 分類）の構成・和則・スケール補正した共回転安定性・直接走行。
入力（read-only）: ../論文v1_全プログラム修正版_20260828/fixed_equimodular/ の step 0 親、
                  ../N{5,6,7,8,10,16,20}_linear124_equimodular_selfconsistent_directHperp_treatment_only_2026082{8,9}/data/states_treatment.npz
出力: results/*.csv, results/run.log
流れ: ż = K(z) z、K_ij = A_ij Im(conj(z_i) z_j)（隣接辺）。自己無撞着 iK(v)v = μv。共回転線形化 L_rot = DF(v) + μ J（J z = i z）。"""
import csv, glob, gzip, os, sys, warnings, itertools
import numpy as np
warnings.filterwarnings('ignore')
HERE=os.path.dirname(os.path.abspath(__file__)); UP=os.path.dirname(HERE); OUT=os.path.join(HERE,'results'); os.makedirs(OUT,exist_ok=True)
FE=os.path.join(UP,'論文v1_全プログラム修正版_20260828','fixed_equimodular')
PK={3:'N3_N4',4:'N3_N4',5:'N5',6:'N6_N7',7:'N6_N7',8:'N8_N9',9:'N8_N9',10:'N10_N11',11:'N10_N11',12:'N12_N13',13:'N12_N13',14:'N14_N15',15:'N14_N15',16:'N16'}
NPZ={5:'N5_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828',6:'N6_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260829',7:'N7_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260829',8:'N8_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828',10:'N10_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828',16:'N16_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828',20:'N20_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828'}
def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def adjacency(N):
    E=edges(N); M=len(E); A=np.zeros((M,M))
    for i in range(M):
        for j in range(M):
            if i!=j and set(E[i])&set(E[j]): A[i,j]=1.0
    return A
def Kof(A,z): return A*np.imag(np.conj(z)[:,None]*z[None,:])
def selfcons(A,v):
    hv=1j*(Kof(A,v)@v); mu=(np.vdot(v,hv)/np.vdot(v,v)).real; return float(np.linalg.norm(hv-mu*v)/np.linalg.norm(v)), float(mu)
def jac(A,z):
    a,b=z.real,z.imag; Sab=A@(a*b); Saa=A@(a*a); Sbb=A@(b*b)
    return np.block([[np.diag(Sab)+A*(np.outer(a,b)-2*np.outer(b,a)), np.diag(-Saa)+A*np.outer(a,a)],[np.diag(Sbb)-A*np.outer(b,b), np.diag(-Sab)+A*(2*np.outer(a,b)-np.outer(b,a))]])
def corot(A,v):
    res,mu=selfcons(A,v); M=len(v); J=np.zeros((2*M,2*M)); J[:M,M:]=-np.eye(M); J[M:,:M]=np.eye(M)
    ev=np.linalg.eigvals(jac(A,v)+mu*J); i=int(np.argmax(ev.real)); r2=float((abs(v)**2).mean())
    return dict(residual=res,mu=mu,mu_over_r2=mu/r2,a=float(ev.real.max()),a_over_r2=float(ev.real.max()/r2),a_over_mu=float(ev.real.max()/abs(mu)),freq_over_mu=float(abs(ev[i].imag)/abs(mu)),n_unstable=int((ev.real>1e-7*max(1,abs(mu))).sum()),r2=r2,
                closure=float(abs((v*v).sum())/(abs(v)**2).sum()),nonzero=int((abs(v)>1e-12).sum()),spread=float(((abs(v)**2).max()-(abs(v)**2).min())/r2))
def step0(N):
    fs=sorted(glob.glob(os.path.join(FE,f'{PK[N]}_complex_simplex_complete_analysis_20260826',f'N{N}_all_steps*'))); f=[x for x in fs if not x.endswith('.gz')] or fs; f=f[0]
    fh=gzip.open(f,'rt') if f.endswith('.gz') else open(f); rd=csv.DictReader(fh)
    if 'edge_index' in rd.fieldnames:
        rows=[]
        for x in rd:
            if x['step']!='0': break
            rows.append(x)
        rows.sort(key=lambda x:int(x['edge_index'])); return np.array([float(x['a'])+1j*float(x['b']) for x in rows])
    x=next(rd); ac=[c for c in rd.fieldnames if c.endswith('_a')]; return np.array([float(x[c]) for c in ac])+1j*np.array([float(x[c[:-2]+'_b']) for c in ac])
def one_factorization(N):
    n=N-1; col={}
    for r in range(n):
        col[tuple(sorted((r,N-1)))]=r
        for k in range(1,N//2): col[tuple(sorted(((r-k)%n,(r+k)%n)))]=r
    return col
def designA(N,signs=None):
    E=edges(N); col=one_factorization(N); th=np.array([col[e]*np.pi/(N-1) for e in E]); 
    if signs is not None: th=th+signs
    return np.exp(1j*th)
def designB(N,G): E=edges(N); return np.exp(1j*np.array([np.pi/2 if (((j-i)%N in G) or ((i-j)%N in G)) else 0.0 for (i,j) in E]))
def colorN(N): E=edges(N); return np.exp(1j*np.array([((i+j)%N)*np.pi/N for (i,j) in E]))
def star(N,equal=True,seed=0):
    E=edges(N); M=len(E); adj=[f for f in range(1,M) if set(E[0])&set(E[f])]; v=np.zeros(M,complex)
    a=np.ones(len(adj)) if equal else np.random.default_rng(seed).normal(size=len(adj)); v[adj]=a; v[0]=1j*np.sqrt((a*a).sum()); return v
log=open(os.path.join(OUT,'run.log'),'w')
def P(*s):
    t=' '.join(str(x) for x in s); print(t); log.write(t+'\n'); log.flush()
# ---- 1. 和則と隣接位相差の分類（make_parent 親, fixed_equimodular step 0）
P('## 1. 和則 Σsin²φ=N−1, Σsin2φ=0 と隣接位相差の分類（fixed_equimodular step 0）')
w=csv.writer(open(os.path.join(OUT,'sum_rules_and_neighbor_angles.csv'),'w',newline='')); w.writerow(['N','M','max_abs_sum_sin2_minus_Nm1','max_abs_sum_sin2phi','frac_90','frac_0_180','frac_mid','per_edge_90','per_edge_0_180','per_edge_mid','need_90','need_0_180','N_mod_4'])
for N in range(3,17):
    z=step0(N); th=np.angle(z); E=edges(N); M=len(z); adj=[[f for f in range(M) if f!=e and set(E[e])&set(E[f])] for e in range(M)]
    S=[];S2=[];c90=[];c0=[];cm=[]
    for e in range(M):
        ph=np.array([th[f]-th[e] for f in adj[e]]); S.append((np.sin(ph)**2).sum()); S2.append(np.sin(2*ph).sum()); d=np.abs(np.degrees((ph+np.pi)%(2*np.pi)-np.pi))
        c90.append(((d>80)&(d<100)).sum()); c0.append(((d<10)|(d>170)).sum()); cm.append(len(d)-c90[-1]-c0[-1])
    tot=sum(len(a) for a in adj); row=[N,M,abs(np.array(S)-(N-1)).max(),np.abs(S2).max(),sum(c90)/tot,sum(c0)/tot,sum(cm)/tot,np.mean(c90),np.mean(c0),np.mean(cm),N-1,N-3,N%4]
    w.writerow(row); P(f"N={N:2d} 和則誤差 {row[2]:.1e} {row[3]:.1e} | 90°/0-180°/中間 = {row[4]:.2f}/{row[5]:.2f}/{row[6]:.2f} 1辺あたり [{row[7]:.1f},{row[8]:.1f},{row[9]:.1f}] 純設計必要 [{N-1},{N-3},0] N mod 4={N%4}")
    if N==3: P('   N=3 位相(deg):',np.round(np.degrees(th),2).tolist(),' 隣接|Δθ|:',sorted(set(round(abs(np.degrees((th[b]-th[a]+np.pi)%(2*np.pi)-np.pi)),1) for a,b in itertools.combinations(range(3),2))))
# ---- 2. 手作り設計の構成と残差 ＋ 3. スケール補正した共回転成長率
P('\n## 2-3. 手作り設計の残差と、スケール補正した共回転成長率 a/r²（振幅レベル、H⊥率は2倍）')
w=csv.writer(open(os.path.join(OUT,'stability_by_parent_type.csv'),'w',newline='')); keys=['residual','mu','mu_over_r2','a','a_over_r2','a_over_mu','freq_over_mu','n_unstable','r2','closure','nonzero','spread']; w.writerow(['type','N','label']+keys)
def rec(typ,N,v,label):
    d=corot(adjacency(N),v); w.writerow([typ,N,label]+[d[k] for k in keys])
    P(f"  {typ:12s} N={N:2d} {label:22s} 残差={d['residual']:.0e} μ/r²={d['mu_over_r2']:+7.3f} 非零辺={d['nonzero']:3d} |z|²幅={d['spread']:5.2f} 閉塞={d['closure']:.0e} | a/r²={d['a_over_r2']:+.3e} a/|μ|={d['a_over_mu']:+.2e} 不安定数={d['n_unstable']} 周波数/|μ|={d['freq_over_mu']:.3f}")
rng=np.random.default_rng(1)
for N in [5,6,7,8,10,16,20]: rec('make_parent',N,np.load(os.path.join(UP,NPZ[N],'data','states_treatment.npz'))['Z'][0],'npz Z[0]')
for N in [4,6,8,10,12,14,16,20]:
    rec('designA',N,designA(N),'1-因子分解 符号+'); rec('designA',N,designA(N,rng.choice([0,np.pi],size=N*(N-1)//2)),'1-因子分解 符号ランダム')
for N,G in [(5,{1}),(9,{1,2}),(9,{1,3}),(9,{1,4}),(13,{1,3,4}),(13,{1,2,3}),(17,{1,2,4,8})]: rec('designB',N,designB(N,G),f'直交2族 G={sorted(G)}')
for N in [5,7,9,11]: rec('colorN_fail',N,colorN(N),'N色彩色 π/N (不成立)')
for N in [7,8]: rec('designB_fail',N,designB(N,{1,2}),'直交 circulant(1,2) (不成立)')
for N in [6,10,16]: rec('star',N,star(N),'1虚+隣接実 等振幅'); rec('star',N,star(N,False),'1虚+隣接実 乱振幅')
# ---- 4. 直接走行
P('\n## 4. 直接走行（種なし、exp(Δ·K) 厳密回転、直接読出し H⊥/H）')
w=csv.writer(open(os.path.join(OUT,'direct_runs.csv'),'w',newline='')); w.writerow(['label','N','tau','Hperp_over_H'])
def run(N,v,tau_max,dt,label):
    A=adjacency(N); Z=v.copy(); p=v.real/np.linalg.norm(v.real); q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)
    steps=int(round(tau_max/dt)); chk={int(round(t/dt)) for t in [0,5,10,20,30,40,60,80,100] if t<=tau_max}; out=[]
    for t in range(steps+1):
        if t in chk:
            Zp=Z-p*(p@Z)-q*(q@Z); f=float(np.vdot(Zp,Zp).real/np.vdot(Z,Z).real); out.append((t*dt,f)); w.writerow([label,N,t*dt,f])
        K=Kof(A,Z); ww,V=np.linalg.eigh(1j*K); Z=V@(np.exp(-1j*dt*ww)*(V.conj().T@Z))
    P(f"  {label:34s} "+'  '.join(f"τ={t:>3.0f}:{f:.1e}" for t,f in out))
z10=np.load(os.path.join(UP,NPZ[10],'data','states_treatment.npz'))['Z'][0]
run(10,designA(10),40,0.002,'N=10 設計A r²=1'); run(10,z10/np.sqrt((abs(z10)**2).mean()),40,0.002,'N=10 make_parent 親 r²=1 に拡大')
run(6,designA(6),60,0.002,'N=6 設計A r²=1'); s6=star(6); run(6,s6/np.sqrt((abs(s6)**2).mean()),60,0.002,'N=6 星型 r²=1 に規格')
run(5,designB(5,{1}),100,0.002,'N=5 設計B δ=0 r²=1')
P('\nresults →',OUT); log.close()
