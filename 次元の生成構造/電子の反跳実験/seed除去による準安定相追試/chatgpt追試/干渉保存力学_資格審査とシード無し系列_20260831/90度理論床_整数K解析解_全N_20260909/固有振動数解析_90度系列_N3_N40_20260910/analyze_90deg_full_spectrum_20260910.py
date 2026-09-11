import os, math, csv, json
from collections import Counter
import numpy as np

CANON='/mnt/data/run_N3_N40_stage123_v1.py'
OUT='/mnt/data/90deg_full_eigenspectrum_20260910'
os.makedirs(OUT,exist_ok=True)

src=open(CANON,encoding='utf-8').read()
g={'__name__':'canonical','__file__':CANON}
exec(compile(src.split('rows=[]; summaries=[]')[0],CANON,'exec'),g)
adjacency,one_step=g['adjacency'],g['one_step']

def integer_K(labels,A):
    diff=(labels[None,:]-labels[:,None])%4
    s=np.zeros_like(diff,dtype=float)
    s[diff==1]=1.0; s[diff==3]=-1.0
    return A*s

def releq_residual(z,A,N):
    w=one_step(z,A,N)
    ip=np.vdot(z,w)
    return float(np.linalg.norm(w-np.exp(1j*np.angle(ip))*z)/np.linalg.norm(z))

def build_analytic_floor(N,max_iter=500):
    ea,eb=np.triu_indices(N,k=1); A=adjacency(N)
    labels=((ea+eb)%4).astype(int)
    best=None; seen={}
    for it in range(max_iter):
        K=integer_K(labels,A)
        ew,V=np.linalg.eigh(1j*K)
        v=V[:,-1]; z=v/np.linalg.norm(v)
        res=releq_residual(z,A,N)
        if best is None or res<best[0]: best=(res,z.copy(),labels.copy())
        newlab=(np.round(np.degrees(np.angle(v))%360/90).astype(int))%4
        cand=None; agmax=-1
        for sh in range(4):
            c=(newlab+sh)%4; ag=int(np.sum(c==labels))
            if ag>agmax: agmax=ag; cand=c
        key=tuple(cand.tolist())
        if np.array_equal(cand,labels) or key in seen: break
        seen[key]=it; labels=cand
    return best[1],best[2],best[0]

def cluster(vals,tol=1e-8):
    vals=sorted(vals)
    if not vals: return []
    groups=[[vals[0]]]
    for x in vals[1:]:
        ref=sum(groups[-1])/len(groups[-1])
        if abs(x-ref) <= tol*max(1.0,abs(ref),abs(x)):
            groups[-1].append(x)
        else: groups.append([x])
    return [(float(sum(gr)/len(gr)),len(gr)) for gr in groups]

all_rows=[]; summary=[]; positive_unique=[]
for N in range(3,41):
    z,labels,res=build_analytic_floor(N)
    A=adjacency(N); K=integer_K(labels,A)
    vals=np.linalg.eigvalsh(1j*K)
    vals=np.real_if_close(vals).astype(float)
    pos=[x for x in vals if x>1e-9]
    cls=cluster(pos,1e-8)
    w0=cls[0][0] if cls else float('nan')
    n_int=0; max_int_err=0.0
    for j,(w,mult) in enumerate(cls,1):
        ratio=w/w0
        near=round(ratio)
        err=abs(ratio-near)
        if err<1e-8: n_int+=1
        max_int_err=max(max_int_err,err)
        row=dict(N=N,M=len(z),mode_index=j,omega=w,multiplicity=mult,
                 omega_over_min=ratio,nearest_integer=int(near),integer_error=err,
                 sigma_max=float(vals[-1]),releq_residual=res)
        all_rows.append(row); positive_unique.append((N,w,mult))
    nzero=int(np.sum(np.abs(vals)<1e-9))
    summary.append(dict(N=N,M=len(z),n_positive_unique=len(cls),n_zero=nzero,
                        omega_min=w0,sigma_max=float(vals[-1]),
                        sigma_over_omega_min=float(vals[-1]/w0),
                        n_integer_multiple_of_min=n_int,
                        all_positive_integer_multiples=(n_int==len(cls)),
                        max_integer_error=max_int_err,
                        releq_residual=res))
    print('done',N,'M',len(z),'unique+',len(cls),'zero',nzero,'w0',w0,'sig',vals[-1],flush=True)

with open(os.path.join(OUT,'spectrum_unique_positive_N3_N40.csv'),'w',newline='',encoding='utf-8') as f:
    wr=csv.DictWriter(f,fieldnames=list(all_rows[0].keys())); wr.writeheader(); wr.writerows(all_rows)
with open(os.path.join(OUT,'spectrum_summary_N3_N40.csv'),'w',newline='',encoding='utf-8') as f:
    wr=csv.DictWriter(f,fieldnames=list(summary[0].keys())); wr.writeheader(); wr.writerows(summary)

# Cross-N exact-ish harmonic pairs: ratios near small integers 2..8
cross=[]
for ia,(N1,w1,m1) in enumerate(positive_unique):
    for N2,w2,m2 in positive_unique[ia+1:]:
        if N2==N1: continue
        loN,lo,ml,hiN,hi,mh=(N1,w1,m1,N2,w2,m2) if w1<w2 else (N2,w2,m2,N1,w1,m1)
        ratio=hi/lo
        k=round(ratio)
        if 2<=k<=8:
            relerr=abs(ratio-k)
            if relerr<1e-8:
                cross.append(dict(N_low=loN,omega_low=lo,mult_low=ml,N_high=hiN,omega_high=hi,mult_high=mh,
                                  harmonic=k,ratio=ratio,error=relerr))
with open(os.path.join(OUT,'crossN_exact_integer_harmonic_pairs.csv'),'w',newline='',encoding='utf-8') as f:
    if cross:
        wr=csv.DictWriter(f,fieldnames=list(cross[0].keys())); wr.writeheader(); wr.writerows(cross)
    else:
        f.write('N_low,omega_low,mult_low,N_high,omega_high,mult_high,harmonic,ratio,error\n')

# Compact markdown report
ints=[r for r in summary if r['all_positive_integer_multiples']]
report=[]
report.append('# 90度自己無撞着床 N=3..40 全固有振動スペクトル解析')
report.append('')
report.append('生成子は各床の Z4 ラベルから K_ef=A_ef sin(90deg(k_f-k_e)) を構成し、Hermitian 行列 iK の全固有値を直接対角化した。')
report.append('')
report.append('## 結果要約')
report.append(f'- 解析N: 3..40 (38床)')
report.append(f'- 各N内で「最小正固有値の整数倍だけ」で全正スペクトルが閉じるN: {[r["N"] for r in ints]}')
report.append(f'- N横断で 2..8 倍の機械精度整数比を持つ固有値ペア数: {len(cross)}')
report.append('')
report.append('|N|M|正の異なる固有値数|zero多重度|omega_min|sigma_max|sigma/omega_min|全て整数倍?|')
report.append('|---:|---:|---:|---:|---:|---:|---:|:---:|')
for r in summary:
    report.append(f"|{r['N']}|{r['M']}|{r['n_positive_unique']}|{r['n_zero']}|{r['omega_min']:.9g}|{r['sigma_max']:.9g}|{r['sigma_over_omega_min']:.9g}|{'YES' if r['all_positive_integer_multiples'] else 'no'}|")
report.append('')
report.append('## 注意')
report.append('- ここでの omega は iK の正固有値を固有振動数候補として読んだスペクトル量。one_step の物理時間への換算は dt=2pi/N を別途掛けて位相進みとして扱う。')
report.append('- 「整数倍でない」は共通基音が最小正固有値ではない可能性を排除しない。次段で全N横断の代数的共通因子・平方根整数構造を調べる。')
open(os.path.join(OUT,'REPORT.md'),'w',encoding='utf-8').write('\n'.join(report)+'\n')

print('SUMMARY_INT_N', [r['N'] for r in ints])
print('CROSS_COUNT',len(cross))
print('CROSS_FIRST20')
for r in cross[:20]: print(r)
