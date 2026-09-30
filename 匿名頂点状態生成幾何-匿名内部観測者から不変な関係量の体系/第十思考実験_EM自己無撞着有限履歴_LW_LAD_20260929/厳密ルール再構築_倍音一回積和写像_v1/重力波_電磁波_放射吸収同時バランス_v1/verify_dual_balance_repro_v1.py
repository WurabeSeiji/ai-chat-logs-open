#!/usr/bin/env python3
import json, glob, hashlib
from pathlib import Path
import numpy as np
import strict_harmonic_rn_tensor_v1 as rn
from search_dual_balance_strict_v2 import run_precompiled

base=Path('.')
items=[]
for fn in ['dual_q0p3_v2.json','dual_q0p6_v2.json','dual_q0p9_v2.json','dual_q0p99_v2.json']:
    d=json.loads((base/fn).read_text())
    f=d['final']; p=d['parameters']; q=p['q_initial_state']
    nx=p['nx']; K=p['K']; dt=p['dt']; steps=p['steps']
    r=np.linspace(2.4,62.4,nx); layout=rn.Layout(nx=nx,K=K)
    M,pi,pj,nrules=rn.compile_sparse_quadratic(layout,r,dt)
    rerun=run_precompiled(q,f['x0'],f['x1'],f['eta'],M,pi,pj,layout,r,dt=dt,steps=steps)
    keys=['ch0_inner','ch0_outer','ch1_inner','ch1_outer']
    absdiff={k:abs(rerun['flux'][k]-f['flux'][k]) for k in keys}
    items.append({'source':fn,'q':q,'saved_final':f,'rerun_final':rerun,
                  'max_abs_flux_difference':max(absdiff.values()),'flux_abs_difference':absdiff,
                  'delta_abs_difference':[abs(rerun['delta'][i]-f['delta'][i]) for i in (0,1)]})
out={'verification':'independent direct rerun of saved optimum initial states','transition_source':'strict_harmonic_rn_tensor_v1.py, same compile_sparse_quadratic and transition','cases':items,
     'all_exact_reproduction':all(x['max_abs_flux_difference']==0.0 and max(x['delta_abs_difference'])==0.0 for x in items)}
Path('dual_balance_repro_verification_v1.json').write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps(out,indent=2,ensure_ascii=False))
