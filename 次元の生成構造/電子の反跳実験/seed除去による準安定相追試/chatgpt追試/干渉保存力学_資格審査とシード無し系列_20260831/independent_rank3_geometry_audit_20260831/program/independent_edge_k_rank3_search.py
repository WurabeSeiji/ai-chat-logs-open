#!/usr/bin/env python3
import argparse, csv, json, math, os, glob
import numpy as np
from scipy.optimize import least_squares

TWOPI=2*np.pi

def edges(N): return [(i,j) for i in range(N) for j in range(i+1,N)]

def phase_freq(x):
    ph=np.unwrap(np.angle(x)); t=np.arange(len(x),dtype=float)
    t-=t.mean(); ph=ph-ph.mean()
    den=np.dot(t,t)
    return abs(np.dot(t,ph)/den)/TWOPI if den>0 else 0.0

def lambdas_window(Z,start,end):
    f=np.array([phase_freq(Z[start:end,j]) for j in range(Z.shape[1])])
    good=f>1e-14
    if not np.all(good):
        raise ValueError('zero/invalid phase frequency')
    return f.max()/f, f

def choose_stable_tail(Z, candidates=(8192,4096,2048,1024), med_tol=.02, p90_tol=.05):
    T=len(Z)
    rows=[]
    chosen=None
    for W in candidates:
        if T < 2*W: continue
        l1,_=lambdas_window(Z,T-W,T)
        l0,_=lambdas_window(Z,T-2*W,T-W)
        rel=np.abs(l1-l0)/np.maximum(l1,l0)
        rec={'W':W,'median_rel':float(np.median(rel)),'p90_rel':float(np.quantile(rel,.9)),'max_rel':float(rel.max())}
        rows.append(rec)
        if chosen is None and rec['median_rel']<=med_tol and rec['p90_rel']<=p90_tol:
            chosen=rec
    if chosen is None:
        chosen=min(rows,key=lambda r:(r['median_rel']+r['p90_rel']))
        chosen=dict(chosen); chosen['fallback']=True
    else:
        chosen=dict(chosen); chosen['fallback']=False
    W=chosen['W']; lam,f=lambdas_window(Z,T-W,T)
    return lam,f,chosen,rows

def centered_gram_from_lengths(N,E,L):
    D2=np.zeros((N,N))
    for (i,j),v in zip(E,L): D2[i,j]=D2[j,i]=v*v
    J=np.eye(N)-np.ones((N,N))/N
    return -0.5*J@D2@J

def gram_metrics(B, reltol=1e-8):
    ev=np.linalg.eigvalsh(B)[::-1]
    scale=max(abs(ev[0]),1e-30)
    rank=int(np.sum(ev>reltol*scale))
    min_rel=float(ev[-1]/scale)
    tail4=float(ev[3]/scale) if len(ev)>3 else 0.0
    return ev,rank,min_rel,tail4

def unpack(v,N):
    X=v.reshape(N-1,3)
    X=np.vstack([np.zeros((1,3)),X])
    return X

def distances(X,E): return np.array([np.linalg.norm(X[i]-X[j]) for i,j in E])

def search_one(lam,N,kmax,starts,max_outer,seed):
    E=edges(N); M=len(E); rng=np.random.default_rng(seed)
    target_rank=min(3,N-1)
    best=None
    # normalize lambda only affects overall scale; min already 1
    for s in range(starts):
        X=rng.normal(size=(N,3)); X-=X.mean(axis=0)
        # initial scale so typical k not huge
        X*=rng.uniform(0.7,4.0)*np.median(lam)/max(np.median(distances(X,E)),1e-9)
        prevk=None
        for it in range(max_outer):
            d=distances(X,E)
            k=np.clip(np.rint(2*d/lam),1,kmax).astype(int)
            # gauge normalization min k=1 is required; if not, scale X toward primitive range and recompute
            if k.min()>1:
                X/=k.min(); d=distances(X,E); k=np.clip(np.rint(2*d/lam),1,kmax).astype(int)
            L=k*lam/2.0
            def fun(v):
                Y=unpack(v,N); return (distances(Y,E)-L)/np.maximum(L,1e-12)
            # anchor vertex 0; initialize relative to vertex0
            v0=(X[1:]-X[0]).reshape(-1)
            r=least_squares(fun,v0,method='trf',max_nfev=500,ftol=1e-11,xtol=1e-11,gtol=1e-11)
            X=unpack(r.x,N); X-=X.mean(axis=0)
            if prevk is not None and np.array_equal(k,prevk): break
            prevk=k.copy()
        d=distances(X,E); L=k*lam/2.0
        rr=(d-L)/np.maximum(L,1e-12)
        rms=float(np.sqrt(np.mean(rr*rr))); mx=float(np.max(np.abs(rr)))
        B=centered_gram_from_lengths(N,E,L); ev,rank,minrel,tail4=gram_metrics(B)
        neg=float(max(0,-ev[-1]/max(abs(ev[0]),1e-30)))
        exact_rank_ok=(rank==target_rank and neg<1e-8)
        score=rms + 10*neg + (0 if rank==target_rank else 0.1*abs(rank-target_rank))
        rec={'score':score,'rms_rel':rms,'max_rel':mx,'rank':rank,'target_rank':target_rank,'neg_rel':neg,'tail4_rel':tail4,'k':k.copy(),'L':L.copy(),'ev':ev.copy(),'X':X.copy(),'start':s,'exact_rank_ok':exact_rank_ok}
        if best is None or rec['score']<best['score']: best=rec
        if exact_rank_ok and rms<1e-7 and mx<1e-6: break
    return best

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',required=True); ap.add_argument('--out-dir',required=True)
    ap.add_argument('--kmax',type=int,default=100); ap.add_argument('--starts',type=int,default=40); ap.add_argument('--max-outer',type=int,default=12); ap.add_argument('--seed',type=int,default=20260831)
    a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True)
    wav=[]; summary=[]; cand=[]; stability=[]
    files=sorted(glob.glob(os.path.join(a.input_dir,'hm_N*_states_treatment.npz')),key=lambda p:int(os.path.basename(p).split('_')[1][1:]))
    for path in files:
        N=int(os.path.basename(path).split('_')[1][1:]); Z=np.load(path)['Z']; E=edges(N)
        if Z.shape[1]!=len(E): raise ValueError((N,Z.shape,len(E)))
        lam,f,ch,rows=choose_stable_tail(Z)
        for r in rows: stability.append({'N':N,**r,'chosen':int(r['W']==ch['W'])})
        for m,((i,j),la,nu) in enumerate(zip(E,lam,f)):
            wav.append({'N':N,'edge':m,'i':i,'j':j,'window':ch['W'],'window_fallback':ch['fallback'],'lambda':la,'frequency':nu})
        starts=a.starts
        best=search_one(lam,N,a.kmax,starts,a.max_outer,a.seed+N)
        pass_exact=bool(best['exact_rank_ok'] and best['rms_rel']<1e-7 and best['max_rel']<1e-6 and best['k'].min()==1)
        summary.append({'N':N,'M':len(E),'window':ch['W'],'window_fallback':ch['fallback'],'starts':starts,'kmax':a.kmax,'rank':best['rank'],'target_rank':best['target_rank'],'rms_rel':best['rms_rel'],'max_rel':best['max_rel'],'neg_rel':best['neg_rel'],'tail4_rel':best['tail4_rel'],'k_min':int(best['k'].min()),'k_max_found':int(best['k'].max()),'pass_exact':pass_exact,'best_start':best['start']})
        for m,((i,j),la,k,L) in enumerate(zip(E,lam,best['k'],best['L'])):
            cand.append({'N':N,'edge':m,'i':i,'j':j,'lambda':la,'k':int(k),'L':L})
        print(f'N={N:2d} W={ch["W"]:4d} rank={best["rank"]} rms={best["rms_rel"]:.3e} max={best["max_rel"]:.3e} neg={best["neg_rel"]:.2e} k=[{best["k"].min()},{best["k"].max()}] PASS={pass_exact}',flush=True)
    def write(name,rows):
        if not rows:return
        with open(os.path.join(a.out_dir,name),'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    write('wavelength_master.csv',wav);write('stability_windows.csv',stability);write('rank3_search_summary.csv',summary);write('rank3_candidate_edges.csv',cand)
    with open(os.path.join(a.out_dir,'run_parameters.json'),'w') as f: json.dump(vars(a),f,indent=2)
if __name__=='__main__': main()
