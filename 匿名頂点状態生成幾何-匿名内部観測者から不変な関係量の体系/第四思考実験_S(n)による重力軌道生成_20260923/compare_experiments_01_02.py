#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare Paper-3 PN reference with Paper-4 experiments 01 and 02."""
import argparse,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read(path):
    return np.genfromtxt(path,delimiter=',',names=True)


def precession(phi,r):
    i=np.where((r[1:-1]<r[:-2])&(r[1:-1]<=r[2:]))[0]+1
    vals=np.r_[phi[0],phi[i]]
    prec=float(np.mean(np.diff(vals)-2*math.pi)) if len(vals)>1 else float('nan')
    return prec, vals, i


def radial_rms_by_phi(ref,mod):
    refphi=np.unwrap(np.arctan2(ref['y'],ref['x'])); refr=np.hypot(ref['x'],ref['y'])
    modphi=np.unwrap(np.arctan2(mod['y'],mod['x'])); modr=np.hypot(mod['x'],mod['y'])
    max_common=min(float(refphi[-1]),float(modphi[-1]))
    q=np.linspace(0.0,max_common,16000)
    rr=np.interp(q,refphi,refr); rm=np.interp(q,modphi,modr)
    rms=float(np.sqrt(np.mean((rr-rm)**2)))
    return rms,rms/30.0


def peri_radii(data):
    phi=np.unwrap(np.arctan2(data['y'],data['x'])); r=np.hypot(data['x'],data['y'])
    i=np.where((r[1:-1]<r[:-2])&(r[1:-1]<=r[2:]))[0]+1
    vals=np.r_[r[0],r[i]]
    ph=np.r_[phi[0],phi[i]]
    return vals,ph


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--order',type=int,default=12); args=ap.parse_args()
    here=Path(__file__).resolve().parent; res=here/'results'; figs=here/'figures'; figs.mkdir(exist_ok=True)
    pn=read(res/'reference_pn_C1.csv')
    e1=read(res/f'experiment01_sn_conservative_C1_N{args.order}.csv')
    e2=read(res/f'experiment02_sn_2p5pn_rr_C1_N{args.order}.csv')

    pnr=np.hypot(pn['x'],pn['y']); pnphi=np.unwrap(np.arctan2(pn['y'],pn['x']))
    e1r=np.hypot(e1['x'],e1['y']); e1phi=np.unwrap(np.arctan2(e1['y'],e1['x']))
    e2r=np.hypot(e2['x'],e2['y']); e2phi=np.unwrap(np.arctan2(e2['y'],e2['x']))
    pp,_,_=precession(pnphi,pnr); p1,_,_=precession(e1phi,e1r); p2,_,_=precession(e2phi,e2r)
    rms1,rel1=radial_rms_by_phi(pn,e1); rms2,rel2=radial_rms_by_phi(pn,e2)
    pr_pn,ph_pn=peri_radii(pn); pr1,ph1=peri_radii(e1); pr2,ph2=peri_radii(e2)

    n=min(len(pr_pn),len(pr1),len(pr2))
    peri_table=[]
    for k in range(n):
        peri_table.append({
          'index':k,
          'PN_r_peri':float(pr_pn[k]),'Exp01_r_peri':float(pr1[k]),'Exp02_r_peri':float(pr2[k]),
          'PN_phi':float(ph_pn[k]),'Exp01_phi':float(ph1[k]),'Exp02_phi':float(ph2[k])})

    summary={
      'case':'C1','a0':30.0,'b0':27.0,'S_order_N':args.order,
      'experiment_01':'conservative S(n)',
      'experiment_02':'same S(n) + local 2.5PN RR recurrence of osculating p,e',
      'PN_precession_rad_per_orbit':pp,
      'Exp01_precession_rad_per_orbit':p1,
      'Exp02_precession_rad_per_orbit':p2,
      'Exp01_precession_error_rad':float(p1-pp),
      'Exp02_precession_error_rad':float(p2-pp),
      'Exp01_radial_rms_M':rms1,'Exp01_radial_rms_over_a':rel1,
      'Exp02_radial_rms_M':rms2,'Exp02_radial_rms_over_a':rel2,
      'radial_rms_improvement_fraction':float(1.0-rms2/rms1),
      'periapsis_comparison':peri_table,
    }
    outjson=res/f'comparison_experiments_01_02_C1_N{args.order}.json'
    with open(outjson,'w',encoding='utf-8') as f: json.dump(summary,f,ensure_ascii=False,indent=2)

    fig,ax=plt.subplots(figsize=(8.6,8.6))
    ax.plot(pn['x'],pn['y'],lw=1.35,label='Paper 3 reference: 2PN + 2.5PN')
    ax.plot(e1['x'],e1['y'],lw=1.0,ls='--',label=f'Experiment 01: conservative S(n), N={args.order}')
    ax.plot(e2['x'],e2['y'],lw=1.05,ls='-.',label=f'Experiment 02: S(n) + 2.5PN RR recurrence, N={args.order}')
    ax.plot([0],[0],marker='+',ms=9,label='origin')
    ax.set_aspect('equal',adjustable='datalim'); ax.set_xlabel('x [M]'); ax.set_ylabel('y [M]')
    ax.set_title('Paper 4: reference PN orbit vs Experiment 01 and Experiment 02')
    ax.legend(fontsize=8); ax.grid(True,alpha=.25)
    txt=(f"radial RMS/a: Exp01={rel1:.3e}, Exp02={rel2:.3e}\n"
         f"precession [rad/orbit]: PN={pp:.5f}, Exp01={p1:.5f}, Exp02={p2:.5f}")
    ax.text(.02,.02,txt,transform=ax.transAxes,fontsize=8,va='bottom',bbox=dict(boxstyle='round',fc='white',alpha=.85))
    fig.tight_layout()
    stem=f'compare_C1_PN_vs_Exp01_Exp02_N{args.order}'
    fig.savefig(figs/f'{stem}.png',dpi=180)
    fig.savefig(figs/f'{stem}.svg')
    plt.close(fig)
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
