
from pathlib import Path
import importlib.util, numpy as np, pandas as pd, math, json, hashlib, zipfile, shutil
from fractions import Fraction
import matplotlib.pyplot as plt

base=Path(__file__).resolve().parent  # PATH PATCH: was Path("/mnt/data")
engine_path=base/"run_n_scaling_lowrank_v1_no_sigma_norm.py"
OUT=base  # PATH PATCH: was base/"N5_physical_phase_step_test" (original outputs sit at package root)
OUT.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location("engraw",engine_path)
eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)

N=5; M=N*(N-1)//2; STEPS=5000; SEED=40260722+1000*N; MAX_DEN=256
def wrap_pi(x): return ((x+np.pi/2)%np.pi)-np.pi/2
def wrap_2pi(x): return ((x+np.pi)%(2*np.pi))-np.pi

def global_q(vals,tol=1e-10,maxq=MAX_DEN):
    vals=np.asarray(vals,float)
    for q in range(1,maxq+1):
        if float(np.max(np.abs(vals-np.round(vals*q)/q)))<tol: return q
    return None

def metrics(theta):
    D=wrap_pi(theta[:,None]-theta[None,:])/np.pi
    vals=D[np.triu_indices(len(theta),1)]
    errs=[]; dens=[]
    for x in vals:
        fr=Fraction(float(x)).limit_denominator(MAX_DEN)
        errs.append(abs(float(x)-float(fr))); dens.append(fr.denominator)
    errs=np.asarray(errs)
    return {
        "phase_pair_max_rat_err_q_le_256":float(errs.max()),
        "phase_pair_median_rat_err_q_le_256":float(np.median(errs)),
        "phase_pair_frac_err_lt_1e10":float(np.mean(errs<1e-10)),
        "phase_pair_smallest_global_q_tol_1e10":global_q(vals),
        "phase_pair_best_denominators":";".join(map(str,sorted(set(dens))))
    }

sys=eng.LowRankSystem(N)
rng=np.random.default_rng(SEED)
v,res,sig=eng.make_parent(sys,rng,iters=1200,beta=0.5,tol=1e-12,restarts=3)
Z=v.copy()

states=[]; increments=[]; phase_rows=[]
prev_theta=None
for t in range(STEPS+1):
    theta=np.angle(Z)
    met=metrics(theta)
    states.append([t,float(np.vdot(Z,Z).real),float(abs(Z@Z)),*[
        met["phase_pair_max_rat_err_q_le_256"],met["phase_pair_median_rat_err_q_le_256"],
        met["phase_pair_frac_err_lt_1e10"],
        np.nan if met["phase_pair_smallest_global_q_tol_1e10"] is None else met["phase_pair_smallest_global_q_tol_1e10"]
    ]])
    for m in range(M):
        phase_rows.append([t,m,float(theta[m]),float((wrap_pi(theta[m]-theta[0]))/np.pi),float(abs(Z[m]))])
    if prev_theta is not None:
        dth=wrap_2pi(theta-prev_theta)/np.pi
        # increment common denominator directly on component increments / pi
        q=global_q(dth,tol=1e-10,maxq=MAX_DEN)
        errs=[]; dens=[]
        for x in dth:
            fr=Fraction(float(x)).limit_denominator(MAX_DEN)
            errs.append(abs(float(x)-float(fr))); dens.append(fr.denominator)
        increments.append([t,float(np.max(np.abs(dth))),float(np.mean(np.abs(dth))),
                           float(max(errs)),float(np.median(errs)),
                           np.nan if q is None else q,
                           ";".join(map(str,sorted(set(dens))))] + list(map(float,dth)))
    if t<STEPS:
        sys.set_state(Z)  # A4
        Z=sys.linear_rotation_step(Z)  # R1
    prev_theta=theta.copy()

states_df=pd.DataFrame(states,columns=["step","H_total","abs_ZtZ","phase_pair_max_rat_err_q_le_256","phase_pair_median_rat_err_q_le_256","phase_pair_frac_err_lt_1e10","phase_pair_smallest_global_q_tol_1e10"])
states_df.to_csv(OUT/"N5_physical_phase_metrics_5000steps.csv",index=False)

inc_cols=["step","max_abs_phase_increment_over_pi","mean_abs_phase_increment_over_pi","max_rat_err_q_le_256","median_rat_err_q_le_256","smallest_global_q_tol_1e10","best_denominators"]+[f"dtheta_edge_{m}_over_pi" for m in range(M)]
inc_df=pd.DataFrame(increments,columns=inc_cols)
inc_df.to_csv(OUT/"N5_phase_increments_5000steps.csv",index=False)
pd.DataFrame(phase_rows,columns=["step","edge_index","theta","relative_theta_over_pi_mod_pi","amplitude"]).to_csv(OUT/"N5_phase_by_edge_5000steps.csv",index=False)

# plots
plt.figure(figsize=(8,5))
plt.semilogy(states_df.step,np.maximum(states_df.phase_pair_max_rat_err_q_le_256,1e-18))
plt.xlabel("physical step"); plt.ylabel("max rational-lock error (q <= 256)")
plt.title("N=5 raw-K physics: relative-phase rational-lock error")
plt.tight_layout(); plt.savefig(OUT/"N5_physical_relative_phase_lock_error.png",dpi=180); plt.close()

plt.figure(figsize=(8,5))
plt.semilogy(inc_df.step,np.maximum(inc_df.max_rat_err_q_le_256,1e-18))
plt.xlabel("physical step"); plt.ylabel("max increment rational-lock error (q <= 256)")
plt.title("N=5 raw-K physics: phase-increment rational-lock error")
plt.tight_layout(); plt.savefig(OUT/"N5_physical_phase_increment_lock_error.png",dpi=180); plt.close()

plt.figure(figsize=(9,5))
phase_mat=np.empty((STEPS+1,M))
for m in range(M):
    vals=pd.DataFrame(phase_rows,columns=["step","edge_index","theta","relative","amp"])
    break
# reconstruct efficiently from saved rows
pdf=pd.DataFrame(phase_rows,columns=["step","edge_index","theta","relative","amp"])
for m in range(M):
    phase_mat[:,m]=pdf[pdf.edge_index==m].relative.to_numpy()
plt.imshow(phase_mat.T,aspect="auto",origin="lower",extent=[0,STEPS,0,M-1])
plt.xlabel("physical step"); plt.ylabel("edge index"); plt.title("N=5 raw-K physics: relative phase / pi (mod pi)")
plt.colorbar(label="relative phase / pi"); plt.tight_layout(); plt.savefig(OUT/"N5_physical_relative_phase_heatmap.png",dpi=180); plt.close()

plt.figure(figsize=(9,5))
dmat=np.array([[row[7+m] for m in range(M)] for row in increments],float)
plt.imshow(dmat.T,aspect="auto",origin="lower",extent=[1,STEPS,0,M-1])
plt.xlabel("physical step"); plt.ylabel("edge index"); plt.title("N=5 raw-K physics: step phase increments / pi")
plt.colorbar(label="delta theta / pi"); plt.tight_layout(); plt.savefig(OUT/"N5_physical_phase_increment_heatmap.png",dpi=180); plt.close()

# summary distribution of q
q_phase=states_df.phase_pair_smallest_global_q_tol_1e10.dropna().astype(int)
q_inc=inc_df.smallest_global_q_tol_1e10.dropna().astype(int)
summary={
"N":N,"M":M,"steps":STEPS,"seed":SEED,"physics":"raw K Cayley; K/sigma normalization removed; seedless parent",
"parent_residual":float(res),
"H_total_min":float(states_df.H_total.min()),"H_total_max":float(states_df.H_total.max()),
"max_abs_ZtZ":float(states_df.abs_ZtZ.max()),
"relative_phase_steps_with_global_q_tol_1e10":int(len(q_phase)),
"relative_phase_global_q_counts":{str(k):int(v) for k,v in q_phase.value_counts().sort_index().items()},
"phase_increment_steps_with_global_q_tol_1e10":int(len(q_inc)),
"phase_increment_global_q_counts":{str(k):int(v) for k,v in q_inc.value_counts().sort_index().items()},
"initial_relative_phase_q":None if np.isnan(states_df.iloc[0].phase_pair_smallest_global_q_tol_1e10) else int(states_df.iloc[0].phase_pair_smallest_global_q_tol_1e10),
"final_relative_phase_q":None if np.isnan(states_df.iloc[-1].phase_pair_smallest_global_q_tol_1e10) else int(states_df.iloc[-1].phase_pair_smallest_global_q_tol_1e10),
"initial_relative_phase_max_error":float(states_df.iloc[0].phase_pair_max_rat_err_q_le_256),
"final_relative_phase_max_error":float(states_df.iloc[-1].phase_pair_max_rat_err_q_le_256),
"phase_increment_max_error_overall":float(inc_df.max_rat_err_q_le_256.max()),
"phase_increment_median_error_overall":float(inc_df.median_rat_err_q_le_256.median())
}
(OUT/"N5_physical_phase_analysis_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
print(json.dumps(summary,indent=2,ensure_ascii=False))
