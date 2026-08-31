#!/usr/bin/env python3
import os, csv, math
import numpy as np

ROOT='/mnt/data/hm_rank3_inputs'
OUT='/mnt/data/amplitude2_geometry_audit'
os.makedirs(OUT, exist_ok=True)

KEY_STEPS=[0,4096,8192,12288,16384,20480,24576,28672,32768,35905,40000]
WINDOW=4096


def edges(N):
    return [(i,j) for i in range(N) for j in range(i+1,N)]

def geom_metrics(N, L):
    E=edges(N)
    L=np.asarray(L,float)
    # scale normalization only for numerical comparability
    pos=L[L>0]
    if len(pos)==0:
        return dict(valid=False, reason='all_zero')
    scale=np.mean(pos)
    Ln=L/scale
    D2=np.zeros((N,N),float)
    for x,(i,j) in zip(Ln,E):
        D2[i,j]=D2[j,i]=x*x
    # triangle margin: minimum (Lik+Lkj-Lij), normalized by mean positive L
    tri_margin=np.inf
    tri_ok=True
    for i in range(N):
        for j in range(i+1,N):
            lij=math.sqrt(D2[i,j])
            for k in range(N):
                if k==i or k==j: continue
                m=math.sqrt(D2[i,k])+math.sqrt(D2[k,j])-lij
                tri_margin=min(tri_margin,m)
                if m < -1e-10:
                    tri_ok=False
    J=np.eye(N)-np.ones((N,N))/N
    B=-0.5*J@D2@J
    ev=np.linalg.eigvalsh(B)
    mu_max=max(float(np.max(np.abs(ev))),1e-300)
    neg_ratio=max(0.0,-float(ev[0]))/mu_max
    psd=ev[0] >= -1e-9*mu_max
    rank=int(np.sum(ev > 1e-8*mu_max)) if psd else None
    axis_ratio=None; centroid_cv=None; eigpos=[]
    if psd:
        ep=ev[ev>1e-8*mu_max]
        eigpos=ep.tolist()
        if len(ep):
            axis_ratio=float(math.sqrt(ep.max()/ep.min()))
        # reconstruct centered coordinates B=X X^T from positive eigenspace
        w,V=np.linalg.eigh(B)
        keep=w>1e-8*mu_max
        X=V[:,keep]*np.sqrt(w[keep])[None,:]
        r=np.linalg.norm(X,axis=1)
        centroid_cv=float(np.std(r)/np.mean(r)) if np.mean(r)>0 else 0.0
    return dict(valid=True, tri_ok=tri_ok, tri_margin=float(tri_margin), psd=psd,
                neg_ratio=float(neg_ratio), rank=rank, axis_ratio=axis_ratio,
                centroid_cv=centroid_cv, lmin=float(L.min()),lmax=float(L.max()),
                lratio=float(L.max()/max(L.min(),1e-300)), eigpos=eigpos)

rows=[]
edge_rows=[]
for N in range(3,17):
    path=os.path.join(ROOT,f'hm_N{N}_states_treatment.npz')
    Z=np.load(path)['Z']
    M=Z.shape[1]
    E=edges(N)
    assert M==len(E)
    # instantaneous key steps
    for t in KEY_STEPS:
        if t>=len(Z): continue
        L=np.abs(Z[t])**2
        g=geom_metrics(N,L)
        row=dict(N=N,kind='snapshot',start=t,end=t+1)
        row.update({k:v for k,v in g.items() if k!='eigpos'})
        rows.append(row)
        for e,(ij,x) in enumerate(zip(E,L)):
            edge_rows.append(dict(N=N,kind='snapshot',start=t,end=t+1,edge=e,i=ij[0],j=ij[1],L_amp2=float(x)))
    # non-overlapping windows plus final window
    starts=list(range(0, len(Z)-WINDOW+1, WINDOW))
    fs=len(Z)-WINDOW
    if fs not in starts: starts.append(fs)
    starts=sorted(set(starts))
    for s in starts:
        e=s+WINDOW
        L=np.mean(np.abs(Z[s:e])**2,axis=0)
        g=geom_metrics(N,L)
        row=dict(N=N,kind='window_mean',start=s,end=e)
        row.update({k:v for k,v in g.items() if k!='eigpos'})
        rows.append(row)
        for idx,(ij,x) in enumerate(zip(E,L)):
            edge_rows.append(dict(N=N,kind='window_mean',start=s,end=e,edge=idx,i=ij[0],j=ij[1],L_amp2=float(x)))

# write
fields=['N','kind','start','end','valid','reason','tri_ok','tri_margin','psd','neg_ratio','rank','axis_ratio','centroid_cv','lmin','lmax','lratio']
with open(os.path.join(OUT,'hm_amplitude2_geometry_time.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT,'hm_amplitude2_edges.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['N','kind','start','end','edge','i','j','L_amp2']); w.writeheader(); w.writerows(edge_rows)

# summary initial/final and first failure
summary=[]
for N in range(3,17):
    rr=[r for r in rows if r['N']==N and r['kind']=='snapshot']
    rr=sorted(rr,key=lambda x:x['start'])
    init=rr[0]; fin=rr[-1]
    first_bad=next((r for r in rr if not (r.get('tri_ok') and r.get('psd'))),None)
    # last good snapshot
    good=[r for r in rr if r.get('tri_ok') and r.get('psd')]
    lg=good[-1] if good else None
    summary.append(dict(N=N,
        initial_rank=init.get('rank'), initial_axis_ratio=init.get('axis_ratio'), initial_centroid_cv=init.get('centroid_cv'),
        final_tri_ok=fin.get('tri_ok'), final_psd=fin.get('psd'), final_neg_ratio=fin.get('neg_ratio'), final_lratio=fin.get('lratio'),
        final_rank=fin.get('rank'), final_axis_ratio=fin.get('axis_ratio'), final_centroid_cv=fin.get('centroid_cv'),
        first_bad_step=None if first_bad is None else first_bad['start'],
        last_good_step=None if lg is None else lg['start'],
        last_good_rank=None if lg is None else lg.get('rank'),
        last_good_axis_ratio=None if lg is None else lg.get('axis_ratio'),
        last_good_centroid_cv=None if lg is None else lg.get('centroid_cv'),
        last_good_lratio=None if lg is None else lg.get('lratio')))
with open(os.path.join(OUT,'hm_amplitude2_geometry_summary.csv'),'w',newline='') as f:
    fields2=list(summary[0].keys()); w=csv.DictWriter(f,fieldnames=fields2); w.writeheader(); w.writerows(summary)

# md
lines=['# hm系列: 振幅二乗 |z_ij|^2 をそのまま辺長とした実ユークリッド幾何監査','',
       '距離候補: L_ij(t)=|z_ij(t)|^2。k_m, λ, 位相は使わない。','',
       '## snapshot 要約','',
       '|N|初期rank|初期axis比|最終Lmax/Lmin|最終実ユークリッド?|最終負固有値比|最初の不成立step|最後の成立step|最後の成立axis比|',
       '|---:|---:|---:|---:|:---:|---:|---:|---:|---:|']
for s in summary:
    eu=bool(s['final_tri_ok'] and s['final_psd'])
    def fmt(x, spec='.4g'):
        return '—' if x is None else format(x,spec)
    lines.append(f"|{s['N']}|{s['initial_rank']}|{fmt(s['initial_axis_ratio'])}|{fmt(s['final_lratio'])}|{'yes' if eu else 'no'}|{fmt(s['final_neg_ratio'],'.3g')}|{s['first_bad_step'] if s['first_bad_step'] is not None else '—'}|{s['last_good_step'] if s['last_good_step'] is not None else '—'}|{fmt(s['last_good_axis_ratio'])}|")
lines += ['', '## 判定規約','',
          '- 距離二乗から中心化 Gram B=-1/2 J D^2 J。',
          '- B が PSD かつ三角不等式を満たすときだけ実ユークリッド配置として rank/主軸を採用。',
          '- axis ratio = sqrt(mu_max/mu_min_positive)。',
          '- centroid CV = 重心から各頂点までの半径の変動係数。',
          '- window_mean では各4096-step窓で mean(|z|^2) を距離候補にした。']
with open(os.path.join(OUT,'analysis_amplitude2_geometry.md'),'w') as f: f.write('\n'.join(lines)+'\n')

print('wrote',OUT)
for s in summary: print(s)
