#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproducible trajectory/readout visualization for Thought Experiment 10.

The plotting conventions follow the earlier Paper 6/7 orbit figures:
explicit readout coordinates, equal x/y scale for spatial plots, start/end markers,
and simultaneous SVG + PNG output.

Interpretation boundary:
  figure01: actual particle positions in the retarded EM/LAD reference model.
  figure02: complex phase/readout trajectories, not spatial orbits.
  figure03: KT physical center worldlines x_phys=a(t)x_coord, not angular orbits.
  figure04: cumulative radiation readout histories, not orbits.
  figure05: strict-RN wave-packet intensity centroids, observation only.
  figure06: dual-balance wave-packet centroids, observation only.
  figure07: dual-balance spacetime intensity maps, observation only.

No plotted readout is fed back into any state transition.
"""
from __future__ import annotations
import csv, hashlib, json, math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
import em_retarded_harmonic_lad_v3 as em
import strict_harmonic_rn_tensor_v1 as srn
import strict_harmonic_kt_tensor_v1 as skt
import strict_harmonic_spin_tensor_v1 as sspin
import search_dual_balance_strict_v2 as dual

Q_CASES=[0.0,0.3,0.6,0.9,0.99]
DUAL_FILES={0.3:'dual_q0p3_v2.json',0.6:'dual_q0p6_v2.json',0.9:'dual_q0p9_v2.json',0.99:'dual_q0p99_v2.json'}

def save_both(fig,stem):
    fig.savefig(HERE/f'{stem}.svg',bbox_inches='tight')
    fig.savefig(HERE/f'{stem}.png',dpi=220,bbox_inches='tight')
    plt.close(fig)

def write_csv(name,header,rows):
    with (HERE/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(header); w.writerows(rows)

def record_em_worldline(params):
    args=SimpleNamespace(**params)
    p=[dict(q=args.q1,m=args.m1),dict(q=args.q2,m=args.m2)]
    coeffs,_=em.initial_coeffs(args,p); t_end=0.0; rows=[]
    def sample(step,t):
        rb=np.stack([em.read_history(coeffs[i]) for i in range(2)],axis=0)
        for i in range(2):
            r=rb[i,-1,:3]; rows.append((step,t,i,float(r[0]),float(r[1]),float(r[2])))
    sample(0,0.0)
    for step in range(args.steps):
        times,rh,uh,a4h,vh,a3h,_=em.derive_readout(coeffs,t_end,args.dt)
        states=[]
        for i in range(2):
            st=(rh[i,-1].copy(),uh[i,-1].copy(),a4h[i,-1].copy())
            ns,_=em.rk4_one(st,i,p,args.dt,t_end,times,rh,vh,a3h); states.append(ns)
        for i,(ri,ui,ai) in enumerate(states): coeffs[i]=em.harmonic_shift_append(coeffs[i],np.r_[ri,ui,ai])
        t_end+=args.dt; sample(step+1,t_end)
    return rows

def plot_em_worldlines():
    saved_r=json.loads((HERE/'em_harmonic_lad_radial_v3.json').read_text())
    saved_t=json.loads((HERE/'em_harmonic_lad_transverse_v3.json').read_text())
    pr=dict(saved_r['parameters']); pt=dict(saved_t['parameters']); pr['out']='unused.json'; pt['out']='unused.json'
    rr=record_em_worldline(pr); rt=record_em_worldline(pt)
    write_csv('trajectory_em_harmonic_particles.csv',['case','step','time','particle','x','y','z'],
              [('radial',)+r for r in rr]+[('transverse',)+r for r in rt])
    fig,axes=plt.subplots(2,2,figsize=(12,9.2))
    for col,(rows,title) in enumerate([(rr,'Radial control'),(rt,'Transverse initial velocity')]):
        a=np.asarray(rows,float); ax=axes[0,col]
        for i in (0,1):
            s=a[a[:,2]==i]; ax.plot(s[:,3],s[:,4],lw=1.2,label=f'particle {i+1}')
            ax.scatter([s[0,3]],[s[0,4]],s=28,marker='o'); ax.scatter([s[-1,3]],[s[-1,4]],s=34,marker='x')
        ax.set_aspect('equal','box'); ax.set_xlim(-.55,.55); ax.set_ylim(-.55,.55); ax.grid(True,alpha=.25)
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_title(title+' — physical x-y trajectory (equal scale)'); ax.legend(fontsize=8)
        ax=axes[1,col]
        for i in (0,1):
            s=a[a[:,2]==i]; t=s[:,1]; ax.plot(t,s[:,3]-s[0,3],lw=1.1,label=f'particle {i+1} Δx')
            ax.plot(t,s[:,4]-s[0,4],lw=1.1,ls='--',label=f'particle {i+1} Δy')
        ax.grid(True,alpha=.25); ax.set_xlabel('t'); ax.set_ylabel('displacement from initial position')
        ax.set_title(title+' — component readout'); ax.legend(fontsize=7,ncol=2)
    fig.suptitle('Harmonic-state retarded EM/LAD: particle trajectories'); fig.tight_layout()
    save_both(fig,'figure01_em_harmonic_particle_trajectories')

def cp(v): return complex(float(v[0]),float(v[1]))
def plot_spin_orbits():
    anchor=json.loads((HERE/'exact_einstein_maxwell_anchor_v1.json').read_text())['spin2_symmetric_square_test']
    psi=np.linspace(0,2*math.pi,721); vals={k:cp(anchor[k]) for k in ['aa','ab','bb','total']}
    N=64;K=4;M,pi,pj,_=sspin.compile_map(K);z=sspin.init(K,N); rec=[]
    for n in range(2*N+1):
        x=sspin.newest(z,K);V,a,b=x[3],x[4],x[5];rec.append((n,V,a*a,2*a*b,b*b))
        if n<2*N:z=sspin.transition(z,M,pi,pj)
    write_csv('trajectory_spin_complex_orbits.csv',['step','V_re','V_im','aa_re','aa_im','ab_re','ab_im','bb_re','bb_im'],
              [(n,V.real,V.imag,aa.real,aa.imag,ab.real,ab.imag,bb.real,bb.imag) for n,V,aa,ab,bb in rec])
    fig,axes=plt.subplots(1,2,figsize=(12,5.6));ax=axes[0]
    for k,z0 in vals.items():
        zz=z0*np.exp(2j*psi);ax.plot(zz.real,zz.imag,lw=1,label=k);ax.scatter([zz[0].real],[zz[0].imag],s=18)
    ax.set_aspect('equal','box');ax.grid(True,alpha=.25);ax.legend(fontsize=8);ax.set_xlabel('Re');ax.set_ylabel('Im');ax.set_title('Exact anchor: weight-2 aa / ab / bb / total')
    ax=axes[1];V=np.array([r[1] for r in rec]);aa=np.array([r[2] for r in rec]);ax.plot(V.real,V.imag,lw=1,label='V (half-angle state)');ax.plot(aa.real,aa.imag,lw=1,label='aa (weight 2 readout)')
    ax.scatter([V[0].real,V[N].real,V[-1].real],[V[0].imag,V[N].imag,V[-1].imag],s=24);ax.set_aspect('equal','box');ax.grid(True,alpha=.25);ax.legend(fontsize=8);ax.set_xlabel('Re');ax.set_ylabel('Im');ax.set_title('Strict harmonic map: phase trajectories, N=64')
    fig.suptitle('Complex phase-orbit readouts in the spin hierarchy');fig.tight_layout();save_both(fig,'figure02_spin_complex_phase_orbits')

def plot_kt_worldlines():
    H=-.05;dt=.25;steps=32;K=4;sep=8.;tt=np.linspace(0,steps*dt,257);a=np.exp(H*tt)
    M,pi,pj,_=skt.compile_map(K);z=skt.init_state(K,H,dt,sep=sep,m1=1,m2=1,a0=1);rows=[]
    for n in range(steps+1):
        v=skt.state_values(z,K);scale=float(v[1]);x1=float(v[5]);x2=float(v[6]);rows.append((n,n*dt,scale,scale*x1,scale*x2,scale*(x2-x1)))
        if n<steps:z=skt.transition(z,M,pi,pj)
    write_csv('trajectory_kt_centers.csv',['step','time','scale','x1_physical','x2_physical','physical_separation'],rows);s=np.array(rows,float)
    fig,axes=plt.subplots(1,2,figsize=(12,5.4));ax=axes[0];ax.plot(tt,-sep*a/2,lw=1.2,label='exact KT center 1');ax.plot(tt,sep*a/2,lw=1.2,label='exact KT center 2');ax.scatter(s[:,1],s[:,3],s=12,label='strict-map center 1');ax.scatter(s[:,1],s[:,4],s=12,label='strict-map center 2');ax.grid(True,alpha=.25);ax.set_xlabel('t');ax.set_ylabel('physical x = a(t) x_coord');ax.set_title('Two-center physical worldlines');ax.legend(fontsize=7)
    ax=axes[1];ax.plot(tt,sep*a,lw=1.2,label='exact KT');ax.scatter(s[:,1],s[:,5],s=14,label='strict harmonic map');ax.grid(True,alpha=.25);ax.set_xlabel('t');ax.set_ylabel('physical separation');ax.set_title('Contracting background separation');ax.legend(fontsize=8)
    fig.suptitle('Kastor-Traschen center trajectories: exact reference vs strict generator');fig.tight_layout();save_both(fig,'figure03_kt_center_worldlines_exact_vs_strict')

def plot_radiation_histories():
    rw=json.loads((HERE/'rw_odd_l2.json').read_text());rn=json.loads((HERE/'q0p90_grav.json').read_text());fig,axes=plt.subplots(1,2,figsize=(12,5.2));ax=axes[0]
    t=np.array([x['t'] for x in rw['samples']]);eh=np.array([x['E_horizon_cum'] for x in rw['samples']])/rw['E0'];ei=np.array([x['E_infinity_cum'] for x in rw['samples']])/rw['E0'];ax.plot(t,eh,label='horizon');ax.plot(t,ei,label='infinity');ax.plot(t,eh+ei,label='total boundary');ax.grid(True,alpha=.25);ax.set_xlabel('t');ax.set_ylabel('cumulative fraction of initial energy');ax.set_title('RWZ spin-2 reference');ax.legend(fontsize=8)
    ax=axes[1];t=np.array([x['t'] for x in rn['samples']]);Hg=np.array([x['H_grav'] for x in rn['samples']])/rn['E0'];He=np.array([x['H_em'] for x in rn['samples']])/rn['E0'];Ig=np.array([x['I_grav'] for x in rn['samples']])/rn['E0'];Ie=np.array([x['I_em'] for x in rn['samples']])/rn['E0'];ax.plot(t,Hg+Ig,label='GW total boundary');ax.plot(t,He+Ie,label='EM total boundary');ax.plot(t,Hg+He+Ig+Ie,label='all boundary');ax.grid(True,alpha=.25);ax.set_xlabel('t');ax.set_ylabel('cumulative fraction of initial energy');ax.set_title('RN coupled, Q/M=0.9, grav initial');ax.legend(fontsize=8)
    fig.suptitle('Radiation readout histories (reference implementations)');fig.tight_layout();save_both(fig,'figure04_rwz_rn_radiation_readout_histories')

def centroid(field,r):
    w=np.asarray(field,float)**2;den=float(w.sum());return float((r*w).sum()/den) if den>1e-30 else math.nan

def record_strict_rn(q,initial_channel,steps=700,stride=4):
    nx=61;K=4;dt=.12;r=np.linspace(2.4,62.4,nx);layout=srn.Layout(nx=nx,K=K);M,pi,pj,_=srn.compile_sparse_quadratic(layout,r,dt);z=srn.init_state(layout,r,q,x0=13.,sigma=4.,amp=1e-4,channel=initial_channel);rows=[]
    for n in range(steps+1):
        if n%stride==0 or n==steps:
            s=srn.decode_state(z,layout);last=K-1;ff=[np.array([s[layout.field_var(ch,i),last].real for i in range(nx)]) for ch in (0,1)];rows.append((n,n*dt,centroid(ff[0],r),centroid(ff[1],r),float((ff[0]**2).sum()),float((ff[1]**2).sum())))
        if n<steps:z=srn.transition(z,M,pi,pj)
    return rows

def plot_strict_rn_centroids():
    allrows=[];fig,axes=plt.subplots(2,5,figsize=(18,7.4),sharex=True,sharey=True)
    for row,ch0 in enumerate((0,1)):
        for col,q in enumerate(Q_CASES):
            rec=record_strict_rn(q,ch0);allrows += [(q,ch0)+tuple(x) for x in rec];a=np.array(rec,float);ax=axes[row,col];ax.plot(a[:,1],a[:,2],lw=1,label='ch0 centroid');ax.plot(a[:,1],a[:,3],lw=1,label='ch1 centroid');ax.set_title(f'init ch{ch0}, Q/M={q:g}');ax.grid(True,alpha=.2)
            if row==1:ax.set_xlabel('t readout')
            if col==0:ax.set_ylabel('radial intensity centroid')
            if row==0 and col==4:ax.legend(fontsize=7)
    write_csv('trajectory_strict_rn_single_centroids.csv',['Q_over_M','initial_channel','step','time','centroid_ch0','centroid_ch1','norm2_ch0','norm2_ch1'],allrows);fig.suptitle('Strict harmonic RN: observation-only wave-packet centroid trajectories');fig.tight_layout();save_both(fig,'figure05_strict_rn_wavepacket_centroid_trajectories')

def record_dual(q,meta,steps=700,stride=4):
    nx=meta['parameters']['nx'];K=meta['parameters']['K'];dt=meta['parameters']['dt'];r=np.linspace(2.4,62.4,nx);layout=srn.Layout(nx=nx,K=K);M,pi,pj,_=srn.compile_sparse_quadratic(layout,r,dt);f=meta['final'];z=dual.init_two_pulses(layout,r,q,f['x0'],f['x1'],f['eta'],sigma=meta['parameters']['sigma'],amp=meta['parameters']['amp_ch0']);rows=[];times=[];h0=[];h1=[]
    for n in range(steps+1):
        if n%stride==0 or n==steps:
            s=srn.decode_state(z,layout);last=K-1;ff=[np.array([s[layout.field_var(ch,i),last].real for i in range(nx)]) for ch in (0,1)];rows.append((n,n*dt,centroid(ff[0],r),centroid(ff[1],r),float((ff[0]**2).sum()),float((ff[1]**2).sum())));times.append(n*dt);h0.append(ff[0]**2);h1.append(ff[1]**2)
        if n<steps:z=srn.transition(z,M,pi,pj)
    return r,np.array(times),np.array(h0),np.array(h1),rows

def plot_dual_balance():
    cases=[];allrows=[]
    for q,fn in DUAL_FILES.items():
        meta=json.loads((HERE/fn).read_text());r,t,h0,h1,rec=record_dual(q,meta,steps=meta['parameters']['steps']);allrows += [(q,meta['final']['x0'],meta['final']['x1'],meta['final']['eta'])+tuple(x) for x in rec];cases.append((q,meta,r,t,h0,h1,np.array(rec,float)))
    write_csv('trajectory_dual_balance_centroids.csv',['Q_over_M','x0_ch0','x0_ch1','eta','step','time','centroid_ch0','centroid_ch1','norm2_ch0','norm2_ch1'],allrows)
    fig,axes=plt.subplots(2,2,figsize=(12,9),sharex=True,sharey=True)
    for ax,(q,meta,r,t,h0,h1,a) in zip(axes.ravel(),cases):
        ax.plot(a[:,1],a[:,2],lw=1,label='ch0 centroid');ax.plot(a[:,1],a[:,3],lw=1,label='ch1 centroid');ax.axhline(meta['final']['x0'],lw=.7,ls='--',label='initial ch0 position');ax.axhline(meta['final']['x1'],lw=.7,ls=':',label='initial ch1 position');d0,d1=meta['final']['delta'];ax.set_title(f'Q/M={q:g}: delta=({d0:.2e}, {d1:.2e})');ax.grid(True,alpha=.2);ax.set_xlabel('t readout');ax.set_ylabel('radial intensity centroid')
    axes[0,0].legend(fontsize=7);fig.suptitle('Dual GW/EM balance: wave-packet centroid trajectories at saved optima');fig.tight_layout();save_both(fig,'figure06_dual_balance_centroid_trajectories')
    fig,axes=plt.subplots(2,2,figsize=(12,9),sharex=True,sharey=True)
    for ax,(q,meta,r,t,h0,h1,a) in zip(axes.ravel(),cases):
        inten=h0+h1;scale=max(float(inten.max()),1e-300);img=np.log10(inten/scale+1e-12);ax.imshow(img,origin='lower',aspect='auto',extent=[r[0],r[-1],t[0],t[-1]],vmin=-12,vmax=0);ax.axvline(meta['final']['x0'],lw=.8,ls='--');ax.axvline(meta['final']['x1'],lw=.8,ls=':');ax.set_title(f'Q/M={q:g}');ax.set_xlabel('radial readout coordinate');ax.set_ylabel('t readout')
    fig.suptitle('Dual-balance strict map: total wave intensity spacetime trajectories (log10 normalized)');fig.tight_layout();save_both(fig,'figure07_dual_balance_wavepacket_spacetime')

def write_readme():
    (HERE/'README_TRAJECTORY_VISUALIZATION_ja_v1.md').write_text('''# 第十思考実験系列：軌道・世界線・波束軌跡 図化 v1\n\n過去の Paper 6 / Paper 7 の軌道図規約（明示的な読出し座標、空間図は同一x/yスケール、開始/終了点、SVG+PNG同時出力）を継承した再現パッケージ。\n\n## 図の意味\n- figure01: retarded EM/LAD 参照実装の粒子位置。\n- figure02: exact anchor / strict spin の複素位相軌道（空間軌道ではない）。\n- figure03: KT 二中心の物理世界線 `x_phys=a(t)x_coord`。\n- figure04: RWZ / RN coupled の累積放射読出し履歴。\n- figure05: strict RN の波束強度重心。\n- figure06: GW/EM 同時バランス最適初期状態の波束強度重心。\n- figure07: 同じ dual-balance run の時空強度図。\n\n波束重心は `r_c = sum_i r_i |psi_i|^2 / sum_i |psi_i|^2`。readout は transition へ一切フィードバックしない。\n\n## 再現\n`python plot_experiment_trajectory_suite_v1.py`\n\n依存: numpy, scipy, matplotlib。元生成コード、入力JSON、過去のPaper 6/7図化コードも同梱。\n''',encoding='utf-8')

def manifest():
    rows=[]
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.name!='MANIFEST_SHA256_v1.txt':rows.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}')
    (HERE/'MANIFEST_SHA256_v1.txt').write_text('\n'.join(rows)+'\n',encoding='utf-8')

def main():
    plot_em_worldlines();plot_spin_orbits();plot_kt_worldlines();plot_radiation_histories();plot_strict_rn_centroids();plot_dual_balance();write_readme();manifest();print(HERE)
if __name__=='__main__':main()
