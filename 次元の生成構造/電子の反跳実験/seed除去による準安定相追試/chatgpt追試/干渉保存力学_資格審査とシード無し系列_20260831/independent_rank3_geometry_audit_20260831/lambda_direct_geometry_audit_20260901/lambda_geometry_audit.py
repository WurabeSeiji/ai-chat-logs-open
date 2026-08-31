import numpy as np, pandas as pd, math, os, json
from pathlib import Path

ROOT=Path('/mnt/data/hm_rank3_inputs')
OUT=Path('/mnt/data/lambda_geometry_audit')
OUT.mkdir(exist_ok=True)
WIN=4096; PAD=4

def edges(N): return [(i,j) for i in range(N) for j in range(i+1,N)]

def dominant_freq_matrix(Z):
    # Z: [T,M] complex. Hann + padded FFT, choose strongest nonzero-frequency peak.
    n,M=Z.shape
    h=np.hanning(n)[:,None]
    nfft=n*PAD
    F=np.fft.fft(Z*h,n=nfft,axis=0)
    P=np.abs(F)**2
    fr=np.fft.fftfreq(nfft)*2*np.pi
    # suppress very near zero bin to avoid DC contamination
    mask=np.abs(fr)>2*np.pi/nfft
    P2=P.copy(); P2[~mask,:]=-1
    idx=np.argmax(P2,axis=0)
    om=fr[idx].astype(float)
    # parabolic interpolation in power around peak, circular indexing
    df=2*np.pi/nfft
    for m,k in enumerate(idx):
        km=(k-1)%nfft; kp=(k+1)%nfft
        y1,y2,y3=P[km,m],P[k,m],P[kp,m]
        den=y1-2*y2+y3
        if den!=0 and np.isfinite(den):
            delta=0.5*(y1-y3)/den
            if abs(delta)<=1:
                om[m]=fr[k]+delta*df
    # crude monochromaticity: main peak / sum of local maxima above 1% main, separated 4 original bins
    mono=[]
    sep=4*PAD
    for m,k in enumerate(idx):
        p=P[:,m]; pmax=p[k]; cand=[]
        # top candidates only
        ids=np.argpartition(p,-min(40,len(p)))[-min(40,len(p)):]
        ids=ids[np.argsort(p[ids])[::-1]]
        for q in ids:
            if p[q] < .01*pmax: break
            if all(min((q-r)%nfft,(r-q)%nfft)>sep for r in cand): cand.append(int(q))
            if len(cand)>=6: break
        tot=sum(p[q] for q in cand) if cand else pmax
        mono.append(float(pmax/tot) if tot>0 else np.nan)
    return om,np.array(mono)

def geometry_from_lambdas(N,lam):
    E=edges(N); D2=np.zeros((N,N))
    for L,(i,j) in zip(lam,E): D2[i,j]=D2[j,i]=L*L
    J=np.eye(N)-np.ones((N,N))/N
    B=-0.5*J@D2@J
    ev,U=np.linalg.eigh(B)
    order=np.argsort(ev)[::-1]; ev=ev[order]; U=U[:,order]
    scale=max(np.max(np.abs(ev)),1e-300)
    # exact/numerical diagnostics
    neg_ratio=max(0.0,-ev[-1]/scale)
    rank_1e8=int(np.sum(ev/scale>1e-8))
    rank_1e6=int(np.sum(ev/scale>1e-6))
    # strict triangle margin
    tri_margin=1e99
    for i in range(N):
      for j in range(i+1,N):
       for k in range(N):
        if k in (i,j): continue
        tri_margin=min(tri_margin, D2[i,k]**0.5+D2[k,j]**0.5-D2[i,j]**0.5)
    # positive axes (numerical 1e-8)
    pos=ev[ev/scale>1e-8]
    axes=np.sqrt(pos) if len(pos) else np.array([])
    axis_ratio=float(axes[0]/axes[-1]) if len(axes)>1 else 1.0
    # centered coordinates from positive part. Diagonal B gives squared radii if PSD.
    radii=np.sqrt(np.clip(np.diag(B),0,None))
    rmean=float(radii.mean()) if len(radii) else np.nan
    rcv=float(radii.std()/rmean) if rmean>0 else np.nan
    rratio=float(radii.max()/radii.min()) if np.all(radii>0) else np.inf
    return ev, neg_ratio, rank_1e8,rank_1e6,tri_margin,axes,axis_ratio,rcv,rratio,B,U

rows=[]; edge_rows=[]
for N in range(3,17):
    fp=ROOT/f'hm_N{N}_states_treatment.npz'
    Zall=np.load(fp)['Z']
    T=len(Zall); M=Zall.shape[1]
    windows=[(s,min(s+WIN,T)) for s in range(0,T-WIN+1,WIN)]
    # exact final 4096 window too
    fw=(T-WIN,T)
    if fw not in windows: windows.append(fw)
    for wi,(a,b) in enumerate(windows):
        Z=Zall[a:b]
        amp=np.abs(Z).mean(axis=0); amax=amp.max()
        live=amp>=1e-6*amax
        om,mono=dominant_freq_matrix(Z)
        valid=live & (np.abs(om)>1e-8)
        if not np.all(valid):
            rows.append(dict(N=N,window_start=a,window_end=b,M=M,all_edges_wave=False,n_live=int(valid.sum())))
            continue
        lam=1/np.abs(om); lam/=lam.min()
        # split-half frequency variation as uncertainty proxy
        om1,_=dominant_freq_matrix(Z[:len(Z)//2]); om2,_=dominant_freq_matrix(Z[len(Z)//2:])
        relerr=np.abs(np.abs(om1)-np.abs(om2))/np.maximum(np.abs(om),1e-300)
        ev,neg,r8,r6,trim,axes,ar,rcv,rratio,B,U=geometry_from_lambdas(N,lam)
        row=dict(N=N,window_start=a,window_end=b,M=M,all_edges_wave=True,n_live=M,
                 median_mono=float(np.nanmedian(mono)),min_mono=float(np.nanmin(mono)),
                 median_rel_freq_split=float(np.nanmedian(relerr)),max_rel_freq_split=float(np.nanmax(relerr)),
                 lambda_min=float(lam.min()),lambda_max=float(lam.max()),lambda_ratio=float(lam.max()/lam.min()),
                 triangle_min_margin=float(trim),triangle_ok=bool(trim>=-1e-10),
                 gram_min_eig_ratio=float(-neg if neg>0 else ev[-1]/max(np.max(np.abs(ev)),1e-300)),
                 gram_neg_ratio=float(neg),rank_1e8=r8,rank_1e6=r6,
                 axis_ratio_max_min=float(ar),centroid_radius_cv=float(rcv),centroid_radius_ratio=float(rratio),
                 eig_ratios=';'.join(f'{x/max(ev[0],1e-300):.9g}' for x in ev),
                 axis_lengths=';'.join(f'{x:.9g}' for x in axes))
        rows.append(row)
        for m,((i,j),L,w,mo,re) in enumerate(zip(edges(N),lam,om,mono,relerr)):
            edge_rows.append(dict(N=N,window_start=a,window_end=b,edge_index=m,i=i,j=j,omega=w,lambda_norm=L,mono=mo,rel_freq_split=re))

rdf=pd.DataFrame(rows); edf=pd.DataFrame(edge_rows)
rdf.to_csv(OUT/'hm_lambda_geometry_time_windows.csv',index=False)
edf.to_csv(OUT/'hm_lambda_edges_time_windows.csv',index=False)

# summary: early, final, and best PSD-like window per N; exact PSD tolerance 1e-8 plus uncertainty-aware flag
summ=[]
for N,g in rdf.groupby('N'):
    g=g[g.get('all_edges_wave',False)==True].copy()
    if len(g)==0: continue
    g['psd_num_1e8']=g['gram_neg_ratio']<=1e-8
    # uncertainty-aware rough flag: negative eig ratio <= 4*median relative lambda error (D^2 ~2 err, safety 2)
    g['psd_within_split_unc']=g['gram_neg_ratio'] <= 4*np.maximum(g['median_rel_freq_split'],1e-12)
    early=g.sort_values('window_start').iloc[0]
    final=g.sort_values('window_end').iloc[-1]
    best=g.sort_values('gram_neg_ratio').iloc[0]
    for kind,x in [('early',early),('final',final),('best_psd',best)]:
      summ.append(dict(N=N,kind=kind,window_start=int(x.window_start),window_end=int(x.window_end),
                       lambda_ratio=x.lambda_ratio,median_mono=x.median_mono,median_rel_freq_split=x.median_rel_freq_split,
                       gram_neg_ratio=x.gram_neg_ratio,psd_num_1e8=bool(x.psd_num_1e8),psd_within_split_unc=bool(x.psd_within_split_unc),
                       rank_1e8=int(x.rank_1e8),rank_1e6=int(x.rank_1e6),axis_ratio_max_min=x.axis_ratio_max_min,
                       centroid_radius_cv=x.centroid_radius_cv,centroid_radius_ratio=x.centroid_radius_ratio,
                       eig_ratios=x.eig_ratios,axis_lengths=x.axis_lengths))
sdf=pd.DataFrame(summ); sdf.to_csv(OUT/'hm_lambda_geometry_summary.csv',index=False)
print(sdf.to_string(index=False, max_colwidth=40))
