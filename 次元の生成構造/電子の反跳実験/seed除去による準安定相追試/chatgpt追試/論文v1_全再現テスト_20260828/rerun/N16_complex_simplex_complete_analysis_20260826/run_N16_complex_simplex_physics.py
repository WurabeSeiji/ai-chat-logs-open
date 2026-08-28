
from pathlib import Path
import importlib.util, numpy as np, pandas as pd, math, json
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ENGINE=HERE/"run_n_scaling_lowrank_v1_no_sigma_norm.py"

spec=importlib.util.spec_from_file_location("eng",ENGINE)
eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)

N=16
M=N*(N-1)//2
STEPS=5000
SEED=40260722+1000*N

def centered_gram(z,n,edges):
    D2=np.zeros((n,n),dtype=complex)
    for val,(i,j) in zip(z*z,edges):
        D2[i,j]=D2[j,i]=val
    J=np.eye(n)-np.ones((n,n))/n
    return -0.5*J@D2@J

sys=eng.LowRankSystem(N)
rng=np.random.default_rng(SEED)
v,res,sig=eng.make_parent(sys,rng,iters=1200,beta=0.5,tol=1e-12,restarts=3)
Z=v.copy()
wp=rng.normal(size=M)
edges=[(i,j) for i in range(N) for j in range(i+1,N)]

raw_rows=[]
summary_rows=[]
axes_rows=[]
for t in range(STEPS+1):
    a=Z.real.copy(); b=Z.imag.copy()
    a2=a*a; b2=b*b; ab=a*b; r2=a2+b2
    z2=Z*Z
    B=centered_gram(Z,N,edges)
    s=np.linalg.svd(B,compute_uv=False)
    rank=int(np.sum(s>max(float(s[0]),1.0)*1e-10))
    tak=np.sqrt(np.maximum(s,0.0))
    # global invariants
    summary_rows.append([
        t,float(np.vdot(Z,Z).real),float(abs(Z@Z)),
        float(a2.sum()),float(b2.sum()),float(ab.sum()),
        float(r2.min()),float(r2.max()),rank
    ])
    axes_rows.append([t]+list(map(float,tak)))
    row=[t]
    for k,(i,j) in enumerate(edges):
        row += [i+1,j+1,float(a[k]),float(b[k]),float(a2[k]),float(b2[k]),float(ab[k]),
                float(r2[k]),float(np.angle(Z[k])),float(z2[k].real),float(z2[k].imag),float(np.angle(z2[k])/np.pi)]
    raw_rows.append(row)
    if t<STEPS:
        sys.set_theta(np.angle(Z))
        se,wp=sys.sigma_max_power(wp)
        Z=sys.cayley_step(Z,se)

# Wide raw data
cols=["step"]
for k,(i,j) in enumerate(edges):
    prefix=f"e{i+1}_{j+1}"
    cols += [f"{prefix}_i",f"{prefix}_j",f"{prefix}_a",f"{prefix}_b",f"{prefix}_a2",f"{prefix}_b2",
             f"{prefix}_ab",f"{prefix}_r2",f"{prefix}_theta",f"{prefix}_z2_re",f"{prefix}_z2_im",f"{prefix}_z2_phase_pi"]
pd.DataFrame(raw_rows,columns=cols).to_csv(HERE/"N16_all_steps_wide.csv",index=False)

pd.DataFrame(summary_rows,columns=["step","H_total","abs_ZtZ","sum_a2","sum_b2","sum_ab","r2_min","r2_max","simplex_rank"]).to_csv(HERE/"N16_global_summary.csv",index=False)
pd.DataFrame(axes_rows,columns=["step"]+[f"takagi_r{k+1}" for k in range(N)]).to_csv(HERE/"N16_takagi_axes.csv",index=False)

# Long final state and selected snapshots
snap_steps=[0,50,100,150,200,250,300,350,400,450,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000]
long_rows=[]
snapshot_rows=[]
for row in raw_rows:
    t=row[0]
    for k,(i,j) in enumerate(edges):
        off=1+12*k
        vals=row[off:off+12]
        rec=[t,k,i+1,j+1]+vals[2:]
        if t in snap_steps: snapshot_rows.append(rec)
        if t==STEPS: long_rows.append(rec)
long_cols=["step","edge_index","i","j","a","b","a2","b2","ab","r2","theta","z2_re","z2_im","z2_phase_pi"]
pd.DataFrame(long_rows,columns=long_cols).to_csv(HERE/"N16_step5000_final_edges.csv",index=False)
pd.DataFrame(snapshot_rows,columns=long_cols).to_csv(HERE/"N16_selected_snapshots_long.csv",index=False)

print(json.dumps({"N":N,"M":M,"steps":STEPS,"seed":SEED,"parent_residual":float(res),"final_rank":int(summary_rows[-1][-1]),
                  "final_r2_min":summary_rows[-1][6],"final_r2_max":summary_rows[-1][7],
                  "final_sum_a2":summary_rows[-1][3],"final_sum_b2":summary_rows[-1][4],"final_sum_ab":summary_rows[-1][5]},indent=2))
