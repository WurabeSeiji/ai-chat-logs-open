#!/usr/bin/env python3
"""Read-only, all-run exploratory analysis of the frozen 99-run survey.
No trajectories are generated or modified. All numerical state files must match
both the original SHA256SUMS and the original QA records before analysis.
Windows are fixed at four non-overlapping 1024-state blocks; t=4096 is retained
for endpoint diagnostics. The full window contains all 4097 states.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, sys, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.linalg import eigh
from scipy import sparse

ROOT=Path('/mnt/data/survey_analysis_20260915')
INPUT=Path('/mnt/data/survey_audit_inputs')
RAW=Path('/mnt/data/survey_raw')
(ROOT/'per_run').mkdir(exist_ok=True)
(ROOT/'timeseries').mkdir(exist_ok=True)
WINDOWS={'full':(0,4097),'W1':(0,1024),'W2':(1024,2048),'W3':(2048,3072),'W4':(3072,4096)}
EPS_AMP=1e-10
QA=json.loads((INPUT/'qa_summary.json').read_text())
QMAP={r['run_id']:r for r in QA['runs']}
MAN=list(csv.DictReader((INPUT/'master_manifest.csv').open()))
MMAP={r['run_id']:r for r in MAN}
EXPECTED={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (INPUT/'SHA256SUMS.txt').read_text().splitlines()}
HASHRUN={r['sha256_states']:r['run_id'] for r in QA['runs']}

def hashfile(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def ranks(evals):
    e=np.maximum(np.asarray(evals,dtype=float),0.0)
    total=float(e.sum())
    if total<=0: return {'r90':None,'r99':None,'r999':None,'Reff':None,'p1':None,'E2':None}, e
    p=e/total; cs=np.cumsum(p)
    out={f'r{k}':int(min(len(p),np.searchsorted(cs,a)+1)) for k,a in [('90',.90),('99',.99),('999',.999)]}
    out.update(Reff=float(1/np.dot(p,p)),p1=float(p[0]),E2=float(p[:2].sum()))
    return out,p

def normrows(a): return np.sqrt(np.sum(np.abs(a)**2,axis=1))
def safe_angle_stats(a,b,scale):
    v=np.abs(a)*np.abs(b)
    ok=(np.abs(a)>EPS_AMP*scale)&(np.abs(b)>EPS_AMP*scale)
    ph=np.angle(a*np.conj(b)); pair=ok[1:]&ok[:-1]
    d=np.angle((a[1:]*np.conj(b[1:]))*np.conj(a[:-1]*np.conj(b[:-1])))
    dv=d[pair]
    if not len(dv): return {'valid_fraction':float(np.mean(ok)),'mean_inc':None,'std_inc':None,'net_principal_turns':None,'directionality':None,'near_pi_fraction':None}
    return {'valid_fraction':float(np.mean(ok)), 'mean_inc':float(dv.mean()),'std_inc':float(dv.std()),
            'net_principal_turns':float(dv.sum()/(2*np.pi)),
            'directionality':float(abs(dv.sum())/max(np.abs(dv).sum(),1e-300)),
            'near_pi_fraction':float(np.mean(np.abs(dv)>.9*np.pi))}

def spectral(Y):
    Yc=Y-Y.mean(axis=0)
    energy=float(np.sum(np.abs(Yc)**2))
    if energy <= 1e-20*float(np.sum(np.abs(Y)**2)):
        return {'dynamic_spectral_valid':False,'peak_bin':None,'peak_period_steps':None,'peak_power_fraction':None,'top3_power_fraction':None,'bins90':None}
    fft=np.fft.fft(Yc,axis=0)
    power=np.sum(np.abs(fft)**2,axis=1); power[0]=0
    power/=power.sum(); ix=int(np.argmax(power)); fs=np.fft.fftfreq(len(Y))
    order=np.sort(power)[::-1]
    return {'dynamic_spectral_valid':True,'peak_bin':int(ix if ix<=len(Y)//2 else ix-len(Y)),
            'peak_period_steps':float(1/abs(fs[ix])), 'peak_power_fraction':float(power[ix]),
            'top3_power_fraction':float(order[:3].sum()),'bins90':int(np.searchsorted(np.cumsum(order),.9)+1)}

def analyze(p,run_id,hsh):
    row=MMAP[run_id]
    with np.load(p,allow_pickle=False) as d:
        Z=d['Z']; pars={k:int(d[k]) for k in ['L','ma','mb','den','T']}
    L,ma,mb,den,T=(pars[k] for k in ['L','ma','mb','den','T']); M=L*(L-1)//2
    assert Z.shape==(4097,M) and Z.dtype==np.complex128,(run_id,Z.shape,Z.dtype)
    assert all(pars[k]==int(row[k]) for k in pars)
    assert np.isfinite(Z).all(),run_id
    ea,eb=np.triu_indices(L,1); delta=eb-ea
    U=np.exp(2j*np.pi/L); ca=np.exp(-1j*np.pi/(2*L)); cb=ca.conjugate()
    B=np.stack([U**(ma*delta),U**(mb*delta)],axis=1)
    expected_z0=ca*B[:,0]+cb*B[:,1]
    initial_diff=float(np.max(np.abs(Z[0]-expected_z0)/np.maximum(1,np.abs(expected_z0))))
    H=np.sum(np.abs(Z)**2,axis=1); Q=np.sum(Z**2,axis=1); H0=float(H[0]); scale=math.sqrt(H0)
    amp=np.abs(Z); weights=amp**2/H[:,None]
    amp_pr=1/np.sum(weights**2,axis=1)
    amp_cv=amp.std(axis=1)/amp.mean(axis=1)
    amp_disp=normrows(amp-amp[0])/scale
    abs_ret=normrows(Z-Z[0])/scale
    overlap=Z@Z[0].conj(); phi=np.angle(overlap)
    proj_ret=normrows(Z-np.exp(1j*phi)[:,None]*Z[0])/scale
    step_move=normrows(np.diff(Z,axis=0))/scale
    step_ip=np.sum(Z[1:]*Z[:-1].conj(),axis=1)
    step_proj=normrows(Z[1:]-np.exp(1j*np.angle(step_ip))[:,None]*Z[:-1])/scale
    coeff=Z@np.linalg.pinv(B).T
    recon=coeff@B.T
    eta=normrows(Z-recon)/np.sqrt(H)
    # Matrix-free line-graph action: A y = B_inc B_inc^T y - 2y.
    Inc=sparse.csr_matrix((np.ones(2*M),(np.r_[np.arange(M),np.arange(M)],np.r_[ea,eb])),shape=(M,L))
    def A(Y): return (Inc@(Inc.T@Y.T)).T-2*Y
    u=np.exp(1j*np.angle(Z))
    KZ=(u.conj()*A(u*Z)-u*A(u.conj()*Z))/(2j)
    gen=normrows(KZ)/np.sqrt(H)
    sin2=(2*(L-2)*M-np.real(np.sum((u*u).conj()*A(u*u),axis=1)))/(4*(L-2)*M)
    quarter_ranks={}; spectra={}; late_vectors=None; train_vectors=None
    for name,(lo,hi) in WINDOWS.items():
        Y=Z[lo:hi]; n=len(Y); mean=Y.mean(axis=0)
        G=Y.conj().T@Y/n; G=(G+G.conj().T)/2
        ev,V=eigh(G,check_finite=False,driver='evr'); ev=ev[::-1]; V=V[:,::-1]
        un,pun=ranks(ev)
        Yc=Y-mean
        Gc=Yc.conj().T@Yc/n; Gc=(Gc+Gc.conj().T)/2
        ec=eigh(Gc,eigvals_only=True,check_finite=False,driver='evr')[::-1]
        dynamic_fraction=float(max(0,np.real(np.trace(Gc)))/np.mean(H[lo:hi]))
        ce,pce=ranks(ec)
        if dynamic_fraction<1e-20: ce={k:None for k in ce}
        qr,_=np.linalg.qr(B); sv=np.linalg.svd(qr.conj().T@V[:,:2].conj(),compute_uv=False)
        angle=np.arccos(np.clip(sv,0,1))
        ac=Y@V[:,:2]
        qrout=dict(un=un,centered=ce,mean_fraction=float(np.vdot(mean,mean).real/np.mean(H[lo:hi])),
                   dynamic_fraction=dynamic_fraction,initial_mode_eta_mean=float(np.mean(eta[lo:hi])),
                   initial_mode_eta_max=float(np.max(eta[lo:hi])),amp_cv_mean=float(np.mean(amp_cv[lo:hi])),
                   amp_pr_fraction_mean=float(np.mean(amp_pr[lo:hi]/M)),generator_norm_mean=float(np.mean(gen[lo:hi])),
                   sin2_mean=float(np.mean(sin2[lo:hi])),principal_angles_initial=angle.tolist(),
                   svd_phase=safe_angle_stats(ac[:,0],ac[:,1],scale),
                   source_phase=safe_angle_stats(coeff[lo:hi,0],coeff[lo:hi,1],1.0),
                   eigenvalues=ev.tolist(),centered_eigenvalues=ec.tolist())
        quarter_ranks[name]=qrout
        if name!='full': spectra[name]=spectral(Y)
        if name=='W3': train_vectors=V[:,:2].copy()
        if name=='W4': late_vectors=V[:,:2].copy()
    # Out-of-window mode generalization (no refitting of the training basis).
    Y=Z[3072:4096]; ac_train=Y@train_vectors
    capture_train=float(np.sum(np.abs(ac_train)**2)/np.sum(np.abs(Y)**2))
    pa_train=np.arccos(np.clip(np.linalg.svd(train_vectors.conj().T@late_vectors,compute_uv=False),0,1)).tolist()
    np.savez_compressed(ROOT/'timeseries'/f'{run_id}.npz',t=np.arange(4097),H=H,Q2=Q,amp_cv=amp_cv,
          amp_pr=amp_pr,amp_disp=amp_disp,return_abs=abs_ret,return_projective=proj_ret,step_move=step_move,
          step_projective=step_proj,eta_initial=eta,source_coeff=coeff,generator_norm=gen,phase_sin2=sin2)
    base={k:(int(v) if k in ['L','M','ma','mb','den','T','ord_a','ord_b','gcd_L_ma','gcd_L_mb','gcd_ma_mb'] else v) for k,v in row.items()}
    base.update(raw_path=str(p),sha256=hsh,n_states=len(Z),initial_diff=initial_diff,
        H0=H0,Q20_re=float(Q[0].real),Q20_im=float(Q[0].imag),chi=float(abs(Q[0])/H0),
        dH_rel_max=float(np.max(np.abs(H-H0))/H0),dQ2_rel_max=float(np.max(np.abs(Q-Q[0]))/max(H0,abs(Q[0]))),
        analytic_fixed=(ma+mb)%L==0,mode_sum=ma+mb,mode_difference=mb-ma,
        min_amp_relative=float(amp.min()/scale),exact_zero_count=int(np.sum(amp==0)),
        orbit_abs_max=float(abs_ret.max()),orbit_amp_max=float(amp_disp.max()),
        step_move_last=float(step_move[-1]),step_move_W4_mean=float(step_move[3072:4096].mean()),
        step_projective_W4_mean=float(step_proj[3072:4096].mean()),
        initial_return_min=float(abs_ret[1:].min()),initial_return_q=int(np.argmin(abs_ret[1:])+1),
        projective_return_min=float(proj_ret[1:].min()),projective_return_q=int(np.argmin(proj_ret[1:])+1),
        generator_norm_initial=float(gen[0]),generator_norm_final=float(gen[-1]),
        amp_cv_initial=float(amp_cv[0]),amp_cv_final=float(amp_cv[-1]),
        W3_basis_capture_W4=capture_train,W3_W4_principal_angles=pa_train,
        W3_basis_phase_W4=safe_angle_stats(ac_train[:,0],ac_train[:,1],scale),
        windows=quarter_ranks,spectra=spectra)
    (ROOT/'per_run'/f'{run_id}.json').write_text(json.dumps(base,indent=2,allow_nan=False))
    return base

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--partial',action='store_true'); args=ap.parse_args()
    paths=list(RAW.glob('*.npz')); found={}; unknown=[]
    for p in paths:
        h=hashfile(p); rid=HASHRUN.get(h)
        if rid is None: unknown.append(str(p)); continue
        assert EXPECTED[f'runs/{rid}/states.npz']==h
        if rid in found: raise RuntimeError(f'Duplicate state data: {rid}')
        found[rid]=(p,h)
    audit={'expected_runs':99,'found_runs':len(found),'missing':sorted(set(MMAP)-set(found)),
           'unrecognized_files':unknown,'all_available_hashes_match':not unknown,
           'source_commit':QA['commit_sha'],'analysis_script_sha256':hashfile(Path(__file__))}
    (ROOT/'raw_hash_audit.json').write_text(json.dumps(audit,indent=2))
    if unknown: raise RuntimeError(f'Unrecognized/corrupt state files: {unknown}')
    if len(found)!=99 and not args.partial: raise RuntimeError(f'{len(found)}/99 available; refusing all-run analysis')
    for j,row in enumerate(MAN):
        rid=row['run_id']
        if rid not in found: continue
        out=ROOT/'per_run'/f'{rid}.json'
        if out.exists(): continue
        start=time.time(); p,h=found[rid]; res=analyze(p,rid,h)
        print(f'ANALYZED {rid} ({time.time()-start:.2f}s)',flush=True)
    if len(found)==99:
        allresults=[json.loads((ROOT/'per_run'/f'{r["run_id"]}.json').read_text()) for r in MAN]
        (ROOT/'all_results.json').write_text(json.dumps(allresults,indent=2))
        flat=[]
        for r in allresults:
            d={k:v for k,v in r.items() if isinstance(v,(str,int,float,bool))}
            for w,s in r['windows'].items():
                for k,v in s['un'].items(): d[f'{w}_{k}']=v
                for k,v in s['centered'].items(): d[f'{w}_centered_{k}']=v
                for k in ['mean_fraction','dynamic_fraction','initial_mode_eta_mean','amp_cv_mean','amp_pr_fraction_mean','generator_norm_mean','sin2_mean']: d[f'{w}_{k}']=s[k]
                for k,v in s['svd_phase'].items(): d[f'{w}_svd_phase_{k}']=v
            for w,s in r['spectra'].items():
                for k,v in s.items(): d[f'{w}_{k}']=v
            flat.append(d)
        with (ROOT/'all_runs_metrics.csv').open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=list(flat[0]));wr.writeheader();wr.writerows(flat)
        print('ALL_99_ANALYSIS_COMPLETE',flush=True)
    else: print(f'PARTIAL_PROCESSING_ONLY {len(found)}/99; no global conclusion',flush=True)
if __name__=='__main__': main()
