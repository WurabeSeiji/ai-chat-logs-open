#!/usr/bin/env python3
import argparse, csv, glob, hashlib, os
import numpy as np

THRESHOLDS=(1e-8,1e-10,1e-12)

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def matrices(z,N):
    ea,eb=np.triu_indices(N,k=1)
    if len(z)!=len(ea): raise ValueError((N,len(z),len(ea)))
    D=np.zeros((N,N),dtype=np.complex128)
    q=z*z
    D[ea,eb]=q; D[eb,ea]=q
    idx=np.arange(1,N)
    d0=D[0,idx]
    G=0.5*(d0[:,None]+d0[None,:]-D[np.ix_(idx,idx)])
    J=np.eye(N)-np.ones((N,N))/N
    Gc=-0.5*J@D@J
    CM=np.ones((N+1,N+1),dtype=np.complex128)
    CM[0,0]=0; CM[1:,1:]=D
    return D,G,Gc,CM

def svd_metrics(A):
    s=np.linalg.svd(A,compute_uv=False)
    smax=float(s[0]) if len(s) else 0.0
    rel=s/smax if smax else np.zeros_like(s)
    ranks={thr:int(np.sum(rel>thr)) for thr in THRESHOLDS}
    return s,rel,ranks

def row_for(z,N,step):
    D,G,Gc,CM=matrices(z,N)
    sg,rg,rkg=svd_metrics(G)
    sc,rc,rkc=svd_metrics(Gc)
    sm,rm,rkm=svd_metrics(CM)
    # Centered Gram always has one translational null direction; second-smallest tracks extra metric radical.
    center_extra_rel=float(rc[-2]) if len(rc)>=2 else np.nan
    return {
      'N':N,'step':step,'M':len(z),
      'closure_abs':float(abs(np.sum(z*z))),
      'centroid_abs':float(abs(np.sum(z))),
      'norm2':float(np.vdot(z,z).real),
      'gram_sigma_max':float(sg[0]),'gram_sigma_min':float(sg[-1]),'gram_sigma_min_rel':float(rg[-1]),
      'gram_rank_rel1e-8':rkg[1e-8],'gram_rank_rel1e-10':rkg[1e-10],'gram_rank_rel1e-12':rkg[1e-12],
      'centered_second_min_rel':center_extra_rel,
      'centered_rank_rel1e-10':rkc[1e-10],
      'cm_sigma_min_rel':float(rm[-1]),'cm_rank_rel1e-10':rkm[1e-10],
    }, sg, sm

def first_cross(rows,key,threshold):
    for r in rows:
        if r[key]>threshold: return int(r['step'])
    return ''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('input_dir')
    ap.add_argument('output_dir')
    args=ap.parse_args()
    os.makedirs(args.output_dir,exist_ok=True)
    files=sorted(glob.glob(os.path.join(args.input_dir,'hm_N*_den_*_states_500.npz')))
    if not files: raise SystemExit('no input npz')
    allrows=[]; spectral=[]; summary=[]; manifest=[]
    for path in files:
        a=np.load(path)
        N=int(a['N']); den=int(a['denominator']); steps=int(a['steps']); Z=a['Z']
        if den!=N or steps!=500 or Z.shape[0]!=501: raise ValueError((path,N,den,steps,Z.shape))
        manifest.append({'N':N,'filename':os.path.basename(path),'sha256':sha256(path),'bytes':os.path.getsize(path)})
        rows=[]
        for t,z in enumerate(Z):
            r,sg,sm=row_for(z,N,t); rows.append(r); allrows.append(r)
            if t in (0,500):
                for i,x in enumerate(sg,1): spectral.append({'N':N,'step':t,'matrix':'reference_gram','index':i,'singular_value':float(x),'relative':float(x/sg[0])})
                for i,x in enumerate(sm,1): spectral.append({'N':N,'step':t,'matrix':'cayley_menger','index':i,'singular_value':float(x),'relative':float(x/sm[0])})
        r0=rows[0]; r5=rows[-1]
        mask=np.array([(r['gram_sigma_min_rel']>1e-14 and r['centroid_abs']>1e-14) for r in rows])
        if mask.sum()>=4:
            x=np.log10([rows[i]['gram_sigma_min_rel'] for i in range(len(rows)) if mask[i]])
            y=np.log10([rows[i]['centroid_abs'] for i in range(len(rows)) if mask[i]])
            corr=float(np.corrcoef(x,y)[0,1])
        else: corr=np.nan
        summary.append({
          'N':N,'M':len(Z[0]),'max_simplex_dim':N-1,
          'rank_step0_rel1e-10':r0['gram_rank_rel1e-10'],'rank_step500_rel1e-10':r5['gram_rank_rel1e-10'],
          'cm_rank_step0_rel1e-10':r0['cm_rank_rel1e-10'],'cm_rank_step500_rel1e-10':r5['cm_rank_rel1e-10'],
          'sigma_min_rel_step0':r0['gram_sigma_min_rel'],'sigma_min_rel_step500':r5['gram_sigma_min_rel'],
          'centroid_abs_step0':r0['centroid_abs'],'centroid_abs_step500':r5['centroid_abs'],
          'closure_abs_step0':r0['closure_abs'],'closure_abs_step500':r5['closure_abs'],
          'first_step_sigma_rel_gt_1e-12':first_cross(rows,'gram_sigma_min_rel',1e-12),
          'first_step_sigma_rel_gt_1e-10':first_cross(rows,'gram_sigma_min_rel',1e-10),
          'first_step_sigma_rel_gt_1e-8':first_cross(rows,'gram_sigma_min_rel',1e-8),
          'log10_corr_sigmaMin_centroid':corr,
        })
    def write_csv(name,rows):
        p=os.path.join(args.output_dir,name)
        with open(p,'w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    write_csv('simplex_dimension_timeseries.csv',allrows)
    write_csv('simplex_dimension_summary.csv',summary)
    write_csv('simplex_singular_values_step0_step500.csv',spectral)
    write_csv('input_manifest.csv',manifest)
    print('WROTE',len(allrows),'timeseries rows')
    for s in summary: print(s)

if __name__=='__main__': main()
