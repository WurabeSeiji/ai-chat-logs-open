#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse,csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(path):
    a=np.genfromtxt(path,delimiter=',',names=True); return a

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--order',type=int,default=12); args=ap.parse_args()
    here=Path(__file__).resolve().parent; res=here/'results'; figs=here/'figures'; figs.mkdir(exist_ok=True)
    pn=read(res/'reference_pn_C1.csv'); sn=read(res/f'sn_gravity_C1_N{args.order}.csv')
    # compare shape after normalizing each by its initial periapsis; parametrizations differ.
    # Nearest polar-angle comparison over first radial cycle.
    pnphi=np.unwrap(np.arctan2(pn['y'],pn['x'])); pnr=np.hypot(pn['x'],pn['y'])
    snphi=np.unwrap(np.arctan2(sn['y'],sn['x'])); snr=np.hypot(sn['x'],sn['y'])
    max_common=min(pnphi[-1],snphi[-1])
    q=np.linspace(0,max_common,12000)
    rp=np.interp(q,pnphi,pnr); rs=np.interp(q,snphi,snr)
    rms=float(np.sqrt(np.mean((rp-rs)**2))); rms_rel=float(rms/30.0)
    # precession estimates from local minima in r vs sequence
    def prec(phi,r):
        i=np.where((r[1:-1]<r[:-2])&(r[1:-1]<=r[2:]))[0]+1
        vals=np.r_[phi[0],phi[i]]
        return float(np.mean(np.diff(vals)-2*math.pi)) if len(vals)>1 else float('nan'), vals
    pp,pvals=prec(pnphi,pnr); sp,svals=prec(snphi,snr)
    summary={'S_order_N':args.order,'radial_rms_difference_M':rms,'radial_rms_difference_over_a':rms_rel,
             'PN_precession_rad_per_orbit':pp,'S_precession_rad_per_orbit':sp,
             'precession_difference_rad':float(sp-pp)}
    with open(res/f'comparison_C1_N{args.order}.json','w',encoding='utf-8') as f: json.dump(summary,f,ensure_ascii=False,indent=2)

    fig,ax=plt.subplots(figsize=(8.2,8.2))
    ax.plot(pn['x'],pn['y'],lw=1.2,label='Paper 3 reference: 2PN + 2.5PN (C1)')
    ax.plot(sn['x'],sn['y'],lw=1.0,ls='--',label=f'S(n) gravity model, N={args.order}')
    ax.plot([0],[0],marker='+',ms=9,label='origin')
    ax.set_aspect('equal',adjustable='datalim'); ax.set_xlabel('x [M]'); ax.set_ylabel('y [M]')
    ax.set_title('Paper 4 numerical experiment: Paper 3 PN orbit vs S(n) orbit')
    ax.legend(fontsize=8); ax.grid(True,alpha=.25)
    txt=(f"radial RMS / a = {rms_rel:.3e}\n"
         f"PN precession = {pp:.5f} rad/orbit\nS(n) precession = {sp:.5f} rad/orbit")
    ax.text(.02,.02,txt,transform=ax.transAxes,fontsize=8,va='bottom',bbox=dict(boxstyle='round',fc='white',alpha=.8))
    fig.tight_layout()
    fig.savefig(figs/f'compare_C1_PN_vs_Sn_N{args.order}.png',dpi=180)
    fig.savefig(figs/f'compare_C1_PN_vs_Sn_N{args.order}.svg')
    plt.close(fig)
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
