#!/usr/bin/env python3
"""Reconstruct the missing N=5 complex-simplex tables and figures.

This is an equivalent reconstruction from the preserved engine/output definitions.
It does not read the missing CSVs as inputs.  The N=5 trajectory is rerun from
run_n_scaling_lowrank_v1_no_sigma_norm.py with the original N-dependent seed.

Required:
  --engine run_n_scaling_lowrank_v1_no_sigma_norm.py
Optional but required for the H_perp milestones / inflation-vs-ordering table:
  --raw-observables N5_raw_K_raw_observables.csv
"""
from __future__ import annotations
import argparse, importlib.util, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

N=5
M=10
STEPS=5000
SEED=40260722+1000*N
KEY_STEPS=[0,1,10,50,81,82,100,200,500,1000,2000,3000,4000,5000]
GROUPS={
    "A_plus":[0,1,9],
    "A_minus":[2,3,4],
    "B_plus":[5,8],
    "B_minus":[6,7],
}

def wrap_pi(x):
    return ((x+np.pi/2)%np.pi)-np.pi/2

def edge_label(i,j):
    return f"{i+1}-{j+1}"

def load_engine(path: Path):
    spec=importlib.util.spec_from_file_location("n5_engine",path)
    if spec is None or spec.loader is None: raise RuntimeError(path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod

def generate_trajectory(engine, trajectory_csv=None):
    """PATCH(Claude, 2026-08-28): trajectory_csv が与えられたら engine を再実行せず、
    その root の力学が書いた N5_phase_by_edge_5000steps.csv（step, edge_index, theta, amplitude）から Z(t) を復元する。
    以降の表・図の計算は無変更。"""
    if trajectory_csv is not None:
        tj=pd.read_csv(trajectory_csv); ea,eb=np.triu_indices(N,k=1); res=float("nan")
        Zs={int(t):(g.sort_values("edge_index").amplitude.to_numpy()*np.exp(1j*g.sort_values("edge_index").theta.to_numpy())) for t,g in tj.groupby("step")}
        class _S: pass
        sys=_S(); sys.ea=ea; sys.eb=eb
    else:
        sys=engine.LowRankSystem(N)
        rng=np.random.default_rng(SEED)
        v,res,_=engine.make_parent(sys,rng,iters=1200,beta=0.5,tol=1e-12,restarts=3)
        Z=v.copy(); wp=rng.normal(size=M)
        ea,eb=sys.ea,sys.eb
    if trajectory_csv is not None: Z=Zs[0]
    rows=[]
    snapshots={}
    for t in range(STEPS+1):
        theta=np.angle(Z); a=Z.real; b=Z.imag
        a2=a*a; b2=b*b; ab=a*b; r2=a2+b2
        a2n=a2/r2; b2n=b2/r2; abn=ab/r2
        zre=a2-b2; zim=2*ab
        zph=np.angle(zre+1j*zim)/np.pi
        rel=wrap_pi(theta-theta[0])/np.pi
        relclass=np.where(a2 < b2,"a2<b2","a2>b2")+np.where(ab < 0,";ab<0",";ab>0")
        c1=a2n+b2n-1.0
        c2=abn*abn-a2n*b2n
        for m in range(M):
            rows.append((t,m,theta[m],rel[m],abs(Z[m]),a[m],b[m],a2[m],b2[m],ab[m],r2[m],
                         a2n[m],b2n[m],abn[m],zre[m],zim[m],zph[m],relclass[m],c1[m],c2[m]))
        if t in KEY_STEPS or t==STEPS:
            snapshots[t]=(Z.copy(),theta.copy())
        if t<STEPS:
            if trajectory_csv is not None: Z=Zs[t+1]
            else: sys.set_theta(theta); se,wp=sys.sigma_max_power(wp); Z=sys.cayley_step(Z,se)
    cols=["step","edge_index","theta","relative_theta_over_pi_mod_pi","amplitude","a","b","a2","b2","ab","r2",
          "a2_norm","b2_norm","ab_norm","z2_re","z2_im","z2_phase_over_pi","relation_class","constraint1","constraint2"]
    return pd.DataFrame(rows,columns=cols), ea, eb, float(res)

def group_metrics(step_df):
    X=step_df[["a2_norm","b2_norm","ab_norm"]].to_numpy(float)
    centers={g:X[idx].mean(axis=0) for g,idx in GROUPS.items()}
    within={g:float(np.max(np.abs(X[idx]-centers[g]))) for g,idx in GROUPS.items()}
    T=lambda x: np.array([x[1],x[0],-x[2]])
    compA=float(np.max(np.abs(T(centers["A_plus"])-centers["A_minus"])))
    compB=float(np.max(np.abs(T(centers["B_plus"])-centers["B_minus"])))
    abdist=float(np.linalg.norm(centers["A_plus"]-centers["B_plus"]))
    return within,centers,compA,compB,abdist

def fixed_first_greedy(values,tol=1e-10):
    centers=[]; members=[]
    for i,v in enumerate(values):
        found=None
        for c,ctr in enumerate(centers):
            if np.max(np.abs(v-ctr)) < tol:
                found=c; break
        if found is None:
            centers.append(v.copy()); members.append([i])
        else: members[found].append(i)
    return members

def forever_step(series, threshold):
    x=np.asarray(series,float)
    good=x < threshold
    suffix=np.logical_and.accumulate(good[::-1])[::-1]
    idx=np.flatnonzero(suffix)
    return int(idx[0]) if len(idx) else None

def circular_group_phase(theta,idx):
    u=np.mean(np.exp(2j*theta[idx]))
    return float(np.angle(u)/2)

def plot_graph(ax, labels=True, square_pyramid=False):
    if square_pyramid:
        pos=np.array([[0,1.3],[-1,-0.4],[1,-0.4],[0.75,-1.2],[-0.75,-1.2]])
    else:
        ang=np.linspace(np.pi/2,np.pi/2+2*np.pi,5,endpoint=False)
        pos=np.c_[np.cos(ang),np.sin(ang)]
    ea=np.array([0,0,0,0,1,1,1,2,2,3]); eb=np.array([1,2,3,4,2,3,4,3,4,4])
    group_for={m:g for g,idx in GROUPS.items() for m in idx}
    styles={"A_plus":"-","A_minus":"--","B_plus":"-.","B_minus":":"}
    for m,(i,j) in enumerate(zip(ea,eb)):
        ax.plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],styles[group_for[m]],lw=1.8,alpha=.8)
        if labels:
            mid=(pos[i]+pos[j])/2; ax.text(mid[0],mid[1],edge_label(i,j),fontsize=7)
    ax.scatter(pos[:,0],pos[:,1],s=100,zorder=5)
    for i,p in enumerate(pos): ax.text(p[0]+.04,p[1]+.04,str(i+1),fontsize=10,weight="bold")
    ax.axis("equal"); ax.axis("off")

def make_figures(out:Path, all_df, pat, summary, ea, eb):
    final=all_df[all_df.step==STEPS].copy()
    labels=[edge_label(i,j) for i,j in zip(ea,eb)]
    # 1 final z2 classes
    fig,ax=plt.subplots(figsize=(6,6))
    ax.scatter(final.z2_re,final.z2_im,s=42)
    for x,y,s in zip(final.z2_re,final.z2_im,labels): ax.annotate(s,(x,y),xytext=(4,4),textcoords="offset points",fontsize=8)
    ax.axhline(0,lw=.6); ax.axvline(0,lw=.6); ax.set_aspect("equal",adjustable="box")
    ax.set_xlabel("Re(z^2) = a^2 - b^2"); ax.set_ylabel("Im(z^2) = 2ab"); ax.set_title("N=5 step 5000: complex squared-distance classes")
    fig.tight_layout(); fig.savefig(out/"N5_final_z2_classes.png",dpi=180); plt.close(fig)
    # 2 four-group convergence
    errcols=["within_A_plus","within_A_minus","within_B_plus","within_B_minus","complement_error_A","complement_error_B"]
    fig,ax=plt.subplots(figsize=(8,5));
    for c in errcols: ax.semilogy(pat.step,np.maximum(pat[c],1e-18),label=c)
    for y in [1e-4,1e-6,1e-8]: ax.axhline(y,lw=.6,ls=":")
    ax.set_xlabel("step"); ax.set_ylabel("four-group component error"); ax.set_title("N=5: convergence to the 3+3+2+2 four-group structure"); ax.legend(fontsize=7,ncol=2)
    fig.tight_layout(); fig.savefig(out/"N5_four_group_convergence.png",dpi=180); plt.close(fig)
    # 3 two-family components
    fig,ax=plt.subplots(figsize=(8,5))
    for prefix in ["Aplus","Bplus"]:
        for c in ["a2n","b2n","abn"]: ax.plot(pat.step,pat[f"{prefix}_{c}"],label=f"{prefix} {c}")
    ax.set_xlabel("step"); ax.set_ylabel("normalized component"); ax.set_title("N=5: A/B family component evolution"); ax.legend(fontsize=8,ncol=2)
    fig.tight_layout(); fig.savefig(out/"N5_two_family_components.png",dpi=180); plt.close(fig)
    # 4 common vs differential phase rotation (period-pi circular centers)
    phase_rows=[]
    for t,g in all_df.groupby("step",sort=True):
        th=g.theta.to_numpy(float)
        ap=circular_group_phase(th,GROUPS["A_plus"]); bp=circular_group_phase(th,GROUPS["B_plus"])
        common=.5*(ap+bp); diff=wrap_pi(ap-bp)
        phase_rows.append((t,common/np.pi,diff/np.pi))
    pdf=pd.DataFrame(phase_rows,columns=["step","common","differential"])
    fig,ax=plt.subplots(figsize=(8,5)); ax.plot(pdf.step,pdf.common,label="common phase / pi"); ax.plot(pdf.step,pdf.differential,label="A+ - B+ phase / pi")
    ax.set_xlabel("step"); ax.set_ylabel("phase / pi"); ax.set_title("N=5: common vs differential phase rotation"); ax.legend()
    fig.tight_layout(); fig.savefig(out/"N5_common_vs_differential_phase_rotation.png",dpi=180); plt.close(fig)
    # 5 square pyramid interpretation
    fig,ax=plt.subplots(figsize=(9,6)); plot_graph(ax,square_pyramid=True); ax.set_title("N=5: 3D square-pyramid readout of the ten relations\n(line style = four-group class)")
    fig.tight_layout(); fig.savefig(out/"N5_square_pyramid_interpretation.png",dpi=170); plt.close(fig)
    # 6 complete infographic
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,3)
    ax=fig.add_subplot(gs[:,0]); plot_graph(ax); ax.set_title("K5 / 10 relations")
    ax=fig.add_subplot(gs[0,1]); ax.scatter(final.z2_re,final.z2_im,s=28); ax.set_title("four z^2 families"); ax.set_xlabel("Re z^2"); ax.set_ylabel("Im z^2")
    ax=fig.add_subplot(gs[1,1]); ferr=pat[errcols].max(axis=1); ax.semilogy(pat.step,np.maximum(ferr,1e-18)); ax.set_title("four-group ordering"); ax.set_xlabel("step"); ax.set_ylabel("max error")
    ax=fig.add_subplot(gs[:,2]); ax.axis("off")
    txt=("N=5 complete analysis\n\n"
         "10 relations -> rank-4 complex simplex\n"
         "final phase/distance groups: 3+3+2+2\n"
         f"four-group <1e-4 forever: step {forever_step(ferr,1e-4)}\n"
         f"four-group <1e-6 forever: step {forever_step(ferr,1e-6)}\n"
         f"four-group <1e-8 forever: step {forever_step(ferr,1e-8)}\n\n"
         "same 10 relations admit a square-pyramid readout\nand nontrivial 2-edge closure in the separate closure analysis")
    ax.text(.03,.95,txt,va="top",fontsize=11)
    fig.tight_layout(); fig.savefig(out/"N5_complete_analysis_infographic.png",dpi=128); plt.close(fig)
    # 7 simplex structure infographic
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[:,0]); plot_graph(ax); ax.set_title("N=5 complex simplex relation graph")
    ax=fig.add_subplot(gs[0,1]);
    for _,r in summary.iterrows(): ax.scatter(r.mean_z2_re,r.mean_z2_im,s=90,label=r.group)
    ax.axhline(0,lw=.5); ax.axvline(0,lw=.5); ax.set_aspect("equal",adjustable="box"); ax.set_title("four final squared-distance centers"); ax.legend(fontsize=8)
    ax=fig.add_subplot(gs[1,1]); ax.axis("off")
    ax.text(.02,.95,"Distance geometry: d_ij^2 = z_ij^2\nN=5 -> M=10 complex distances\nCentered complex Gram rank = N-1 = 4\n\nA-/B- are the 90-degree component complements\nof A+/B+; group sizes are 3,3,2,2.",va="top",fontsize=11)
    fig.tight_layout(); fig.savefig(out/"N5_complex_simplex_structure_infographic.png",dpi=128); plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--engine",type=Path,required=False,default=None)
    ap.add_argument("--trajectory-csv",type=Path,default=None,help="PATCH: N5_phase_by_edge_5000steps.csv から軌道を読む（engine 再実行なし）")
    ap.add_argument("--raw-observables",type=Path,default=None)
    ap.add_argument("--outdir",type=Path,required=True)
    args=ap.parse_args(); out=args.outdir; out.mkdir(parents=True,exist_ok=True)
    engine=load_engine(args.engine) if args.engine else None
    all_df,ea,eb,parent_res=generate_trajectory(engine, trajectory_csv=args.trajectory_csv)
    all_df.to_csv(out/"N5_all_steps_a_b_a2_b2_ab.csv",index=False)
    # four group pattern
    prows=[]
    for t,g in all_df.groupby("step",sort=True):
        within,c,ca,cb,dist=group_metrics(g)
        prows.append((int(t),within["A_plus"],within["A_minus"],within["B_plus"],within["B_minus"],ca,cb,dist,
                      *c["A_plus"],*c["B_plus"]))
    pcols=["step","within_A_plus","within_A_minus","within_B_plus","within_B_minus","complement_error_A","complement_error_B","A_B_center_distance",
           "Aplus_a2n","Aplus_b2n","Aplus_abn","Bplus_a2n","Bplus_b2n","Bplus_abn"]
    pat=pd.DataFrame(prows,columns=pcols); pat.to_csv(out/"N5_final_four_group_pattern_by_step.csv",index=False)
    # relation counts
    cnt=all_df.groupby(["step","relation_class"]).size().unstack(fill_value=0)
    order=["a2<b2;ab<0","a2<b2;ab>0","a2>b2;ab<0","a2>b2;ab>0"]
    for c in order:
        if c not in cnt: cnt[c]=0
    cnt=cnt[order].reset_index(); cnt.to_csv(out/"N5_relation_class_counts_by_step.csv",index=False)
    # selected clustering
    krows=[]
    for t in KEY_STEPS:
        g=all_df[all_df.step==t].sort_values("edge_index")
        vals=g[["a2_norm","b2_norm","ab_norm"]].to_numpy(float)
        for cid,mem in enumerate(fixed_first_greedy(vals,1e-10)):
            ctr=vals[mem].mean(axis=0); krows.append((t,cid,len(mem),*ctr))
    pd.DataFrame(krows,columns=["step","class_id","count","a2_norm","b2_norm","ab_norm"]).to_csv(out/"N5_key_step_ab2_ab_classes.csv",index=False)
    # final group summary
    final=all_df[all_df.step==STEPS].sort_values("edge_index").reset_index(drop=True)
    srows=[]
    for gn,idx in GROUPS.items():
        q=final.iloc[idx]; zre=q.z2_re.mean(); zim=q.z2_im.mean()
        srows.append((gn,",".join(edge_label(ea[i],eb[i]) for i in idx),len(idx),q.r2.mean(),q.a2.mean(),q.b2.mean(),q.ab.mean(),zre,zim,np.angle(zre+1j*zim)/np.pi))
    summary=pd.DataFrame(srows,columns=["group","edges","count","mean_r2","mean_a2","mean_b2","mean_ab","mean_z2_re","mean_z2_im","mean_z2_phase_over_pi"])
    summary.to_csv(out/"N5_step5000_four_group_summary.csv",index=False)
    ferr=pat[["within_A_plus","within_A_minus","within_B_plus","within_B_minus","complement_error_A","complement_error_B"]].max(axis=1)
    milestone=[]
    raw=None
    if args.raw_observables is not None and args.raw_observables.exists():
        raw=pd.read_csv(args.raw_observables)
        hf=float(raw.loc[raw.step==STEPS,"H_perp"].iloc[0])
        for p in [.5,.9,.95,.99]:
            hit=raw.loc[raw.H_perp >= p*hf,"step"]
            milestone.append((f"H_perp >= {int(p*100)}% final",int(hit.iloc[0])))
    for th in [1e-1,1e-2,1e-3,1e-4,1e-6,1e-8]: milestone.append((f"four-group error < {th:g} forever",forever_step(ferr,th)))
    pd.DataFrame(milestone,columns=["metric","step"]).to_csv(out/"N5_time_separation_milestones.csv",index=False)
    if raw is not None:
        x=raw[["step","H_perp","A_perp","H_total"]].merge(pat[["step","A_B_center_distance"]],on="step",how="inner")
        x["four_group_error"]=ferr.to_numpy()[:len(x)]
        x=x[["step","H_perp","A_perp","H_total","four_group_error","A_B_center_distance"]]
        x.to_csv(out/"N5_inflation_vs_ordering_timeseries.csv",index=False)
    else:
        pd.DataFrame({"step":pat.step,"H_perp":np.nan,"A_perp":np.nan,"H_total":np.nan,"four_group_error":ferr,"A_B_center_distance":pat.A_B_center_distance}).to_csv(out/"N5_inflation_vs_ordering_timeseries.csv",index=False)
    make_figures(out,all_df,pat,summary,ea,eb)
    print(f"N5 reconstructed in {out}")
    print(f"parent residual={parent_res:.3e}")
    print("forever thresholds:", {"1e-4":forever_step(ferr,1e-4),"1e-6":forever_step(ferr,1e-6),"1e-8":forever_step(ferr,1e-8)})
if __name__=="__main__": main()
