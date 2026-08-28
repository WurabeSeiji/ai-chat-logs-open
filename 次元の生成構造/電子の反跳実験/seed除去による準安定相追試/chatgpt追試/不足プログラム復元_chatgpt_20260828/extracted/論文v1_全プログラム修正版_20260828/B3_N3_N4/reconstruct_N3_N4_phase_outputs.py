#!/usr/bin/env python3
"""Reconstruct the N=3 and N=4 special phase-ordering outputs/figures."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def wrap_pi(x): return ((x+np.pi/2)%np.pi)-np.pi/2

def circular_diff_pi(a,b): return abs(wrap_pi(a-b))/np.pi

def circ_mean_pi(theta):
    u=np.mean(np.exp(2j*np.asarray(theta,float)))
    return float((np.angle(u)/2)%np.pi)

def edge_map(n):
    i,j=np.triu_indices(n,k=1); return i,j,[f'{a+1}-{b+1}' for a,b in zip(i,j)]

def greedy_1d(vals,tol=1e-8):
    centers=[]; mem=[]
    for v in vals:
        hit=None
        for k,c in enumerate(centers):
            if abs(v-c) < tol: hit=k; break
        if hit is None: centers.append(v); mem.append([v])
        else: mem[hit].append(v)
    rows=[(float(np.mean(x)),len(x)) for x in mem]
    return sorted(rows,key=lambda z:z[0])

def n3_outputs(df,out):
    ea,eb,labels=edge_map(3); final=df[df.step==5000].sort_values('edge_index').reset_index(drop=True)
    th=final.theta.to_numpy(float); mod=np.mod(th,np.pi)/np.pi
    rows=[(float(x),lab) for x,lab in zip(mod,labels)]
    pd.DataFrame([(i,1,x,lab) for i,(x,lab) in enumerate(rows)],columns=['class_id','count','theta_over_pi_mod1','edges']).to_csv(out/'N3_final_phase_classes.csv',index=False)
    ds=sorted(circular_diff_pi(th[i],th[j]) for i in range(3) for j in range(i+1,3))
    pd.DataFrame([(x,1) for x in ds],columns=['abs_delta_theta_over_pi','count']).to_csv(out/'N3_pairwise_phase_difference_classes.csv',index=False)
    erows=[]
    for st,g in df.groupby('step',sort=True):
        t=g.sort_values('edge_index').theta.to_numpy(float)
        d=np.array([circular_diff_pi(t[i],t[j]) for i in range(3) for j in range(i+1,3)])
        e=np.abs(d-1/3)
        erows.append((int(st),float(e.max()),float(e.mean())))
    edf=pd.DataFrame(erows,columns=['step','max_abs_error_from_one_third','mean_abs_error_from_one_third'])
    edf.to_csv(out/'N3_phase_difference_convergence_to_one_third.csv',index=False)
    # figures
    fig,ax=plt.subplots(figsize=(8,5)); ax.semilogy(edf.step,np.maximum(edf.max_abs_error_from_one_third,1e-18),label='max'); ax.semilogy(edf.step,np.maximum(edf.mean_abs_error_from_one_third,1e-18),label='mean'); ax.set_xlabel('step'); ax.set_ylabel('error from 1/3'); ax.set_title('N=3: phase-difference convergence toward one third'); ax.legend(); fig.tight_layout(); fig.savefig(out/'N3_phase_convergence_to_one_third.png',dpi=180); plt.close(fig)
    pos=np.array([[0,1],[-.9,-.6],[.9,-.6]])
    fig,ax=plt.subplots(figsize=(6,6))
    for m,(i,j) in enumerate(zip(ea,eb)):
        ax.plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],lw=2); mid=(pos[i]+pos[j])/2; ax.text(*mid,f'{labels[m]}\n{mod[m]:.3f} pi',ha='center',fontsize=9)
    ax.scatter(pos[:,0],pos[:,1],s=120); [ax.text(pos[i,0]+.05,pos[i,1]+.05,str(i+1)) for i in range(3)]; ax.axis('off'); ax.set_title('N=3 triangle: final phase classes'); fig.tight_layout(); fig.savefig(out/'N3_triangle_phase_classes.png',dpi=180); plt.close(fig)
    z=(final.z2_re.to_numpy()+1j*final.z2_im.to_numpy())/final.r2.to_numpy()
    phase_circle(z,labels,'N=3 final normalized z^2 phases',out/'N3_final_z2_phase_circle.png')
    return float(edf.iloc[-1].max_abs_error_from_one_third)

def phase_circle(z,labels,title,path):
    fig,ax=plt.subplots(figsize=(6,6)); ang=np.linspace(0,2*np.pi,400); ax.plot(np.cos(ang),np.sin(ang),lw=.8); ax.scatter(z.real,z.imag,s=48)
    for x,y,s in zip(z.real,z.imag,labels): ax.annotate(s,(x,y),xytext=(4,4),textcoords='offset points',fontsize=8)
    ax.axhline(0,lw=.5); ax.axvline(0,lw=.5); ax.set_aspect('equal'); ax.set_xlim(-1.15,1.15); ax.set_ylim(-1.15,1.15); ax.set_xlabel('Re(z^2)/|z|^2'); ax.set_ylabel('Im(z^2)/|z|^2'); ax.set_title(title); fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)

def n4_outputs(df,out):
    ea,eb,labels=edge_map(4); pair_idx=[(0,5),(1,4),(2,3)]
    final=df[df.step==5000].sort_values('edge_index').reset_index(drop=True); th=final.theta.to_numpy(float)
    rows=[]
    for cid,(a,b) in enumerate(pair_idx): rows.append((cid,2,circ_mean_pi(th[[a,b]])/np.pi,f'{labels[a]},{labels[b]}'))
    pd.DataFrame(rows,columns=['class_id','count','theta_over_pi_mod1','edges']).to_csv(out/'N4_final_phase_classes.csv',index=False)
    orows=[]
    for st,g in df.groupby('step',sort=True):
        t=g.sort_values('edge_index').theta.to_numpy(float)
        within=[circular_diff_pi(t[a],t[b]) for a,b in pair_idx]
        ctr=np.sort(np.array([circ_mean_pi(t[[a,b]])/np.pi for a,b in pair_idx]))
        gaps=np.diff(np.r_[ctr,ctr[0]+1.0])
        orows.append((int(st),float(max(within)),float(np.max(np.abs(gaps-1/3))),*map(float,gaps)))
    odf=pd.DataFrame(orows,columns=['step','max_within_opposite_pair_phase_diff','max_class_gap_error_from_one_third','gap1','gap2','gap3'])
    odf.to_csv(out/'N4_opposite_edge_phase_ordering.csv',index=False)
    all_d=[circular_diff_pi(th[i],th[j]) for i in range(6) for j in range(i+1,6)]
    cls=greedy_1d(sorted(all_d),1e-8)
    pd.DataFrame(cls,columns=['abs_delta_theta_over_pi','count']).to_csv(out/'N4_pairwise_phase_difference_classes.csv',index=False)
    fig,ax=plt.subplots(figsize=(8,5)); ax.semilogy(odf.step,np.maximum(odf.max_within_opposite_pair_phase_diff,1e-18),label='within opposite pair'); ax.semilogy(odf.step,np.maximum(odf.max_class_gap_error_from_one_third,1e-18),label='three-class gap error'); ax.set_xlabel('step'); ax.set_ylabel('phase error / pi'); ax.set_title('N=4: opposite-edge phase ordering'); ax.legend(); fig.tight_layout(); fig.savefig(out/'N4_phase_ordering.png',dpi=180); plt.close(fig)
    # tetrahedral schematic (2D projection), opposite pair styles shared
    pos=np.array([[0,1.1],[-1,-.6],[1,-.6],[0,.1]])
    fig,ax=plt.subplots(figsize=(7,6)); styles=['-','--',':']
    pclass={0:0,5:0,1:1,4:1,2:2,3:2}
    for m,(i,j) in enumerate(zip(ea,eb)):
        ax.plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],styles[pclass[m]],lw=2); mid=(pos[i]+pos[j])/2; ax.text(*mid,labels[m],fontsize=8)
    ax.scatter(pos[:,0],pos[:,1],s=120); [ax.text(pos[i,0]+.04,pos[i,1]+.04,str(i+1)) for i in range(4)]; ax.axis('off'); ax.set_title('N=4 tetrahedron: three opposite-edge classes'); fig.tight_layout(); fig.savefig(out/'N4_tetrahedron_opposite_edge_classes.png',dpi=180); plt.close(fig)
    z=(final.z2_re.to_numpy()+1j*final.z2_im.to_numpy())/final.r2.to_numpy(); phase_circle(z,labels,'N=4 final normalized z^2 phases',out/'N4_final_z2_phase_circle.png')
    return float(odf.iloc[-1].max_class_gap_error_from_one_third)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--n3-all-steps',type=Path,required=True); ap.add_argument('--n4-all-steps',type=Path,required=True); ap.add_argument('--outdir',type=Path,required=True); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    d3=pd.read_csv(a.n3_all_steps); d4=pd.read_csv(a.n4_all_steps)
    e3=n3_outputs(d3,a.outdir); e4=n4_outputs(d4,a.outdir)
    # comparison summary belongs to B4 too, but emitting it here reproduces the original N3/N4 package completely.
    rows=[]
    for n,d,err,classes,sizes in [(3,d3,e3,3,'1+1+1'),(4,d4,e4,3,'2+2+2')]:
        f=d[d.step==5000]
        rows.append((f'N={n}',n,n*(n-1)//2,n-1,float(f.r2.min()),float(f.r2.max()),classes,sizes,err))
    pd.DataFrame(rows,columns=['case','N','M','simplex_rank','final_r2_min','final_r2_max','final_phase_distance_classes_tol1e-8','class_sizes','phase_thirds_error']).to_csv(a.outdir/'N3_N4_comparison_summary.csv',index=False)
    print('N3/N4 special phase outputs reconstructed:',a.outdir)
if __name__=='__main__': main()
