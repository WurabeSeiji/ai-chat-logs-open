#!/usr/bin/env python3
import os,csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT='/mnt/data/denominator_controls_20260901'; IN='/mnt/data/hm_rank3_inputs'; FIG=os.path.join(ROOT,'figures'); os.makedirs(FIG,exist_ok=True)
OFFSETS=[-2,-1,0,1,2]

def plane(v):
 p=v.real/np.linalg.norm(v.real); q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q); return p,q
allrows=[]; sums=[]
for off in OFFSETS:
 for N in range(3,17):
  den=N+off; fn=os.path.join(ROOT,f'hm_N{N}_den_{den}_offset_{off:+d}_states_500.npz')
  Z=np.load(fn)['Z']; v=Z[0]; p,q=plane(v)
  vals=[]
  h0=float(np.vdot(v,v).real)
  for t,z in enumerate(Z):
   h=float(np.vdot(z,z).real); zp=z-p*(p@z)-q*(q@z); hp=float(np.vdot(zp,zp).real)/h; clo=float(abs(z@z)/h)
   allrows.append(dict(offset=off,N=N,denominator=den,delta=2*math.pi/den,step=t,H_total=h,Hperp_frac=hp,global_closure=clo)); vals.append(hp)
  vals=np.array(vals); ix=np.where(vals>0.05)[0]
  sums.append(dict(offset=off,N=N,denominator=den,delta=2*math.pi/den,onset_Hperp_005=(int(ix[0]) if len(ix) else ''),max_Hperp_frac=float(vals.max()),final_Hperp_frac=float(vals[-1]),H_drift_abs=abs(float(np.vdot(Z[-1],Z[-1]).real)-h0)))
with open(os.path.join(ROOT,'timeseries_metrics_denominator_controls_500.csv'),'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(allrows[0].keys())); w.writeheader(); w.writerows(allrows)
with open(os.path.join(ROOT,'summary_denominator_controls_500.csv'),'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(sums[0].keys())); w.writeheader(); w.writerows(sums)

D={(o,N):[] for o in OFFSETS for N in range(3,17)}
for r in allrows: D[(int(r['offset']),int(r['N']))].append(r)
for off in OFFSETS:
 label=f'N{off:+d}' if off else 'N'
 # full scale
 for ymin,suffix in [(1e-34,''),(1e-5,'_floor1e-5')]:
  fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
  for i,N in enumerate(range(3,17)):
   ax=axs[i]; rr=D[(off,N)]; x=np.array([r['step'] for r in rr]); y=np.array([r['Hperp_frac'] for r in rr])
   ax.semilogy(x,np.maximum(y,1e-40),color='#1f77b4',lw=1.1,label='hm')
   ax.set_ylim(ymin,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
   if i%4==0: ax.set_ylabel('Hperp/H')
   if i>=10: ax.set_xlabel('step')
  for j in range(14,16): axs[j].axis('off')
  fig.suptitle(f'Hperp/H over 500 steps (delta=2pi/({label}), seedless, direct readout, interference-preserving frame): handmade equimodular',fontsize=11)
  plt.tight_layout(rect=(0,0,1,0.97)); out=os.path.join(FIG,f'fig3_Hperp_grid_N3_N16_den_{label}{suffix}.png'); plt.savefig(out,dpi=130); plt.close()
 # closure same original style
 fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
 for i,N in enumerate(range(3,17)):
  ax=axs[i]; rr=D[(off,N)]; x=np.array([r['step'] for r in rr]); y=np.array([r['global_closure'] for r in rr])
  ax.semilogy(x,np.maximum(y,1e-18),color='#1f77b4',lw=1.1,label='hm'); ax.set_ylim(1e-18,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
  if i%4==0: ax.set_ylabel('|sum Z^2| / H')
  if i>=10: ax.set_xlabel('step')
 for j in range(14,16): axs[j].axis('off')
 fig.suptitle(f'closure |sum Z^2|/H over 500 steps: delta=2pi/({label}), interference-preserving frame',fontsize=11)
 plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,f'fig8_closure_grid_N3_N16_den_{label}.png'),dpi=130); plt.close()

# one overlay figure, same per-N panel layout, for direct control comparison
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
 ax=axs[i]
 for off in OFFSETS:
  rr=D[(off,N)]; x=np.array([r['step'] for r in rr]); y=np.array([r['Hperp_frac'] for r in rr]); label=f'N{off:+d}' if off else 'N'
  ax.semilogy(x,np.maximum(y,1e-40),lw=1.0,label=label)
 ax.set_ylim(1e-34,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
 if i%4==0: ax.set_ylabel('Hperp/H')
 if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('Hperp/H denominator control: 2pi/(N-2), 2pi/(N-1), 2pi/N, 2pi/(N+1), 2pi/(N+2)',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig3_Hperp_denominator_control_overlay.png'),dpi=130); plt.close()
print('DONE')
