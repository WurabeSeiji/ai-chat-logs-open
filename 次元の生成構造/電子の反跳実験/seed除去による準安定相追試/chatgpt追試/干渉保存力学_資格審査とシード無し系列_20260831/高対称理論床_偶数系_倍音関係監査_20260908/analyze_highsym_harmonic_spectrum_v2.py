#!/usr/bin/env python3
import csv, json, math
from pathlib import Path
from fractions import Fraction
import numpy as np

HERE=Path(__file__).resolve().parent
EVEN_NS=list(range(4,41,2))
ABS_TOL=1e-9
REL_TOL=1e-9
RATIONAL_MAX_DEN=16


def edges_for_n(N):
    return [(i,j) for i in range(N) for j in range(i+1,N)]

def cyc_dist(i,j,N):
    return min((j-i)%N,(i-j)%N)

def build_A(edges):
    m=len(edges); A=np.zeros((m,m),float); S=[set(e) for e in edges]
    for a in range(m):
        for b in range(a+1,m):
            if S[a]&S[b]: A[a,b]=A[b,a]=1.0
    return A

def build_theoretical_floor(N):
    edges=edges_for_n(N); M=len(edges); D=N//2
    d=np.array([cyc_dist(i,j,N) for i,j in edges],int)
    theta=np.array([2*np.pi*(i+j)/N for i,j in edges])
    A=build_A(edges)
    K=A*np.sin(theta[None,:]-theta[:,None]); H=1j*K
    B=np.zeros((M,D),complex); nd=[]
    for dd in range(1,D+1):
        idx=np.where(d==dd)[0]; nd.append(len(idx))
        B[idx,dd-1]=np.exp(1j*theta[idx])/np.sqrt(len(idx))
    Q=B.conj().T@H@B
    w,V=np.linalg.eigh(Q)
    chosen=None
    for kk in np.argsort(w):
        c=V[:,kk]; p=int(np.argmax(np.abs(c))); c*=np.exp(-1j*np.angle(c[p]))
        if np.max(np.abs(c.imag))<1e-9:
            cr=c.real
            if np.all(cr>=-1e-10) or np.all(cr<=1e-10):
                if np.sum(cr)<0: cr=-cr
                chosen=(float(w[kk]),cr); break
    if chosen is None:
        raise RuntimeError(f'No same-sign reduced eigenvector N={N}')
    _,c=chosen
    r=np.array([c[dd-1]/np.sqrt(nd[dd-1]) for dd in d])
    z=r*np.exp(1j*theta); z/=np.linalg.norm(z)
    th=np.angle(z); K2=A*np.sin(th[None,:]-th[:,None]); H2=1j*K2
    lam=float(np.real(np.vdot(z,H2@z)/np.vdot(z,z)))
    eigres=float(np.linalg.norm(H2@z-lam*z)/np.linalg.norm(z))
    return z,H2,lam,eigres

def cluster(vals, atol=1e-8, rtol=1e-8):
    vals=np.sort(np.asarray(vals,float))
    out=[]
    for x in vals:
        if not out or abs(x-out[-1][-1]) > atol + rtol*max(abs(x),abs(out[-1][-1])):
            out.append([x])
        else:
            out[-1].append(x)
    return [(float(np.mean(g)),len(g),float(max(g)-min(g))) for g in out]

def close(a,b):
    return abs(a-b) <= ABS_TOL + REL_TOL*max(abs(a),abs(b),1.0)

EMPTY_HEADERS={
    'integer_harmonic_relations.csv':['N','base_index','harmonic_index','multiple','omega_base','omega_harmonic','abs_error','rel_error'],
    'low_denominator_rational_relations.csv':['N','low_index','high_index','omega_low','omega_high','ratio','p','q','ratio_error'],
    'three_wave_sum_resonances.csv':['N','i','j','k','omega_i','omega_j','omega_k','abs_error','rel_error'],
}

def write_csv(path,rows):
    fields=list(rows[0].keys()) if rows else EMPTY_HEADERS.get(path.name)
    if not fields:
        raise RuntimeError(f'No field definition for empty CSV: {path.name}')
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def main():
    summary=[]; modes=[]; harmonics=[]; rationals=[]; threewave=[]
    for N in EVEN_NS:
        z,H,lam,eigres=build_theoretical_floor(N)
        ev=np.linalg.eigvalsh(H).real
        zero_count=int(np.sum(np.abs(ev)<1e-10))
        pos=ev[ev>1e-9]
        cls=cluster(pos)
        uniq=np.array([x[0] for x in cls])
        for idx,(v,mult,spread) in enumerate(cls):
            modes.append({'N':N,'mode_index':idx,'omega':v,'multiplicity':mult,'cluster_spread':spread})
        # integer harmonics among unique positive frequencies
        hcount=0
        for i,a in enumerate(uniq):
            for j,b in enumerate(uniq):
                if j<=i: continue
                for k in (2,3,4):
                    err=abs(b-k*a)
                    rel=err/max(abs(b),abs(k*a),1.0)
                    if err <= 1e-8 + 1e-8*max(abs(b),abs(k*a),1.0):
                        harmonics.append({'N':N,'base_index':i,'harmonic_index':j,'multiple':k,'omega_base':a,'omega_harmonic':b,'abs_error':err,'rel_error':rel})
                        hcount+=1
        # low-denominator rational ratios, excluding 1/1 and exact integer already retained separately
        rcount=0
        for i,a in enumerate(uniq):
            for j,b in enumerate(uniq):
                if j<=i: continue
                ratio=b/a
                fr=Fraction(float(ratio)).limit_denominator(RATIONAL_MAX_DEN)
                approx=fr.numerator/fr.denominator
                err=abs(ratio-approx)
                if err <= 1e-8 and not (fr.denominator==1 and fr.numerator in (2,3,4)):
                    rationals.append({'N':N,'low_index':i,'high_index':j,'omega_low':a,'omega_high':b,'ratio':ratio,'p':fr.numerator,'q':fr.denominator,'ratio_error':err})
                    rcount+=1
        # three-wave positive-frequency sum resonance a+b=c; i<=j<k
        twcount=0
        for i,a in enumerate(uniq):
            for j in range(i,len(uniq)):
                b=uniq[j]
                target=a+b
                kidx=int(np.searchsorted(uniq,target))
                for k in (kidx-1,kidx,kidx+1):
                    if k>j and 0<=k<len(uniq):
                        c=uniq[k]; err=abs(c-target)
                        if err <= 1e-8 + 1e-8*max(abs(c),abs(target),1.0):
                            threewave.append({'N':N,'i':i,'j':j,'k':k,'omega_i':a,'omega_j':b,'omega_k':c,'abs_error':err,'rel_error':err/max(abs(c),1.0)})
                            twcount+=1
        # floor carrier relation: every component belongs to one relative-equilibrium carrier lambda
        summary.append({
            'N':N,'M':len(z),'floor_lambda':lam,'floor_eigen_residual':eigres,
            'full_positive_mode_count':len(pos),'unique_positive_frequency_count':len(uniq),'zero_mode_count':zero_count,
            'integer_harmonic_relations_2to4':hcount,'low_den_rational_relations_q_le_16':rcount,
            'three_wave_sum_resonances':twcount,'photon_pair_complete_expected_N8k':int(N%8==0)
        })
        print(f'N={N:2d} M={len(z):3d} uniq+={len(uniq):2d} harm={hcount:3d} rat={rcount:3d} 3w={twcount:3d} lambda={lam:.12g}')
    write_csv(HERE/'even_N_harmonic_summary.csv',summary)
    write_csv(HERE/'positive_frequency_modes.csv',modes)
    write_csv(HERE/'integer_harmonic_relations.csv',harmonics)
    write_csv(HERE/'low_denominator_rational_relations.csv',rationals)
    write_csv(HERE/'three_wave_sum_resonances.csv',threewave)
    meta={'even_N':EVEN_NS,'abs_tol':ABS_TOL,'rel_tol':REL_TOL,'cluster_tol':1e-8,'rational_max_den':RATIONAL_MAX_DEN,
          'spectrum':'positive eigenvalues of H=iK(theta_floor), with theta_floor from high-symmetry theoretical-floor ansatz',
          'important_limitation':'This audits the frozen-generator spectrum at the relative-equilibrium floor. It is not yet the Jacobian/Floquet spectrum of the full state-dependent stage1+2+3 one-step map.'}
    (HERE/'RUN_METADATA.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__': main()
