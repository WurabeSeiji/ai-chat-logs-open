#!/usr/bin/env python3
import json, math
from pathlib import Path
import numpy as np

# Geometric units G=c=4*pi*eps0=1.


def rn_exact(M,Q,r):
    f=1.0-2.0*M/r+(Q*Q)/(r*r)
    E=Q/(r*r)
    disc=M*M-Q*Q
    if disc>0:
        horizons=[M+math.sqrt(disc),M-math.sqrt(disc)]
        cls='subextremal_black_hole'
    elif abs(disc)<1e-14:
        horizons=[M,M]
        cls='extremal_black_hole'
    else:
        horizons=[]
        cls='overextremal_naked_singularity'
    return dict(f=f,g_tt=-f,g_rr=1.0/f,E_r=E,class_name=cls,horizons=horizons)


def mp_two_center(ma,mb,xa,xb,x):
    ra=float(np.linalg.norm(x-xa)); rb=float(np.linalg.norm(x-xb))
    ua=ma/ra; ub=mb/rb
    U=1.0+ua+ub
    # Exact Majumdar-Papapetrou metric convention ds^2=-U^-2 dt^2+U^2 dx^2
    spatial=U*U
    gtt=-1.0/(U*U)
    At=1.0/U  # gauge choice up to sign/additive constant
    # Exact decomposition of spatial conformal factor U^2.
    parts={
      'background':1.0,
      'linear_a':2.0*ua,
      'linear_b':2.0*ub,
      'aa':ua*ua,
      'ab':2.0*ua*ub,
      'bb':ub*ub,
    }
    return dict(ra=ra,rb=rb,ua=ua,ub=ub,U=U,g_tt=gtt,spatial_factor=spatial,A_t=At,
                decomposition=parts,reconstruction=sum(parts.values()),
                reconstruction_error=abs(sum(parts.values())-spatial))


def spin2_readout(Aa,Ab,phia,phib,angles):
    a=Aa*np.exp(1j*phia); b=Ab*np.exp(1j*phib)
    aa=a*a; ab=2*a*b; bb=b*b; total=(a+b)**2
    rows=[]; maxerr=0.0
    for psi in angles:
        # common transverse-frame phase weight +1 carrier rotation
        ar=a*np.exp(1j*psi); br=b*np.exp(1j*psi)
        terms=np.array([ar*ar,2*ar*br,br*br,(ar+br)**2])
        base=np.array([aa,ab,bb,total])*np.exp(2j*psi)
        err=float(np.max(np.abs(terms-base)))
        maxerr=max(maxerr,err)
        rows.append(dict(psi=float(psi),max_error=err))
    return dict(aa=[aa.real,aa.imag],ab=[ab.real,ab.imag],bb=[bb.real,bb.imag],
                total=[total.real,total.imag],max_spin2_covariance_error=maxerr,tests=rows)


def half_angle_readout(phi):
    V=np.exp(0.5j*phi)
    V2=np.exp(0.5j*(phi+2*np.pi))
    V4=np.exp(0.5j*(phi+4*np.pi))
    return dict(V=[V.real,V.imag],after_2pi=[V2.real,V2.imag],after_4pi=[V4.real,V4.imag],
                err_2pi_sign=float(abs(V2+V)),err_4pi_return=float(abs(V4-V)),
                weight2_from_V4=[(V**4).real,(V**4).imag])


def main(out):
    ratios=[0.0,0.1,0.5,1.0,2.0,10.0]
    rn=[]
    M=1.0; r=20.0
    for lam in ratios:
        Q=lam*M
        rn.append(dict(q_over_m=lam,**rn_exact(M,Q,r)))

    xa=np.array([-1.5,0.0,0.0]); xb=np.array([1.5,0.0,0.0])
    points=[np.array([0.0,2.0,0.0]),np.array([0.5,3.0,0.0]),np.array([4.0,1.0,0.0])]
    mp=[dict(point=p.tolist(),**mp_two_center(1.0,0.7,xa,xb,p)) for p in points]

    spin=spin2_readout(1.2,0.8,0.37,-0.61,np.linspace(0,2*np.pi,33))
    half=half_angle_readout(0.73)

    result={
      'units':'G=c=4*pi*eps0=1',
      'rn_exact_sweep':rn,
      'mp_two_center_exact':mp,
      'spin2_symmetric_square_test':spin,
      'spin_half_double_cover_test':half,
      'interpretation':{
        'rn':'Exact Einstein-Maxwell single-center family; Q/M spans neutral, mass-dominated, extremal, and charge-dominated regimes. Q/M>1 is exact but over-extremal.',
        'mp':'Exact two-center Einstein-Maxwell solution at extremal same-sign charge/mass balance. aa, ab, bb appear together inside U^2, not as separately added forces.',
        'spin2':'If a and b each carry phase weight 1 under a transverse-frame rotation, aa, 2ab, bb and their sum all carry weight 2 exactly.',
        'spin_half':'V=exp(i phi/2) changes sign under 2pi and returns under 4pi; V^4 carries phase weight 2.'
      }
    }
    Path(out).write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps({
      'rn':[(x['q_over_m'],x['class_name'],x['f']) for x in rn],
      'mp_max_reconstruction_error':max(x['reconstruction_error'] for x in mp),
      'spin2_covariance_error':spin['max_spin2_covariance_error'],
      'spin_half_2pi_sign_error':half['err_2pi_sign'],
      'spin_half_4pi_return_error':half['err_4pi_return'],
      'sample_mp_cross_terms':[{k:v for k,v in x['decomposition'].items() if k in ('aa','ab','bb')} for x in mp]
    },ensure_ascii=False,indent=2))

if __name__=='__main__':
    main('/mnt/data/em_selfconsistent_v1/exact_einstein_maxwell_anchor_v1.json')
