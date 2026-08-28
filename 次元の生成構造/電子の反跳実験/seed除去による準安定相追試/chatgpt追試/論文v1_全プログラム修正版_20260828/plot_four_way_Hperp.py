# -*- coding: utf-8 -*-
"""四者（原本 / 修正版 / baseline / 等モジュラー）の H⊥/H_total（global_summary の差し引き読出し、床 ~1e-16）を N=3..16 で重ね描き。
出力: results/figures/four_way_Hperp_frac_N3_N16.png, four_way_Hperp_frac_N5_N16.png"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__)); O=os.path.join(os.path.dirname(H),"論文v1_全再現テスト_20260828","original")
ROOTS={"原本 (Cayley, 位相のみK, 親正規化)":O,"修正版 (exp, 振幅込みK)":os.path.join(H,"fixed"),"baseline (exp, 位相のみK)":os.path.join(H,"fixed_baseline"),"等モジュラー (exp, 振幅込みK, 等モジュラー親)":os.path.join(H,"fixed_equimodular")}
COL=["#888888","#d62728","#1f77b4","#2ca02c"]
pk={3:"N3_N4",4:"N3_N4",5:"N5",6:"N6_N7",7:"N6_N7",8:"N8_N9",9:"N8_N9",10:"N10_N11",11:"N10_N11",12:"N12_N13",13:"N12_N13",14:"N14_N15",15:"N14_N15",16:"N16"}
plt.rcParams["font.family"]=["Hiragino Sans","Arial Unicode MS","DejaVu Sans"]
def load(root,N):
    pkg=f"{pk[N]}_complex_simplex_complete_analysis_20260826"
    for f in (f"N{N}_global_summary.csv",f"N{N}_physical_phase_metrics_5000steps.csv"):
        p=os.path.join(root,pkg,f)
        if os.path.exists(p):
            d=pd.read_csv(p)
            if "H_perp" in d and "H_total" in d: return d.step.to_numpy(), (d.H_perp/d.H_total).to_numpy()
            if "H_perp_frac" in d: return d.step.to_numpy(), d.H_perp_frac.to_numpy()
    # N=5, N=16 の global_summary には H_perp が無いので、N5 は K_sigma raw 枝（N5_raw_K_raw_observables.csv）、N16 は decompactification の geometry_summary を用いる
    if N==5:
        p=os.path.join(root,"K_sigma_normalization_artifact_test_N4_N5_20260826","N5_raw_K_raw_observables.csv")
        if os.path.exists(p): d=pd.read_csv(p); return d.step.to_numpy(), (d.H_perp/d.H_total).to_numpy()
    if N==16:
        p=os.path.join(root,"complex_simplex_decompactification_N5_N16_20260826","results","N16_geometry_summary.csv")
        if os.path.exists(p): d=pd.read_csv(p); return d.step.to_numpy(), (d.H_perp/d.H_total).to_numpy()
    return None
for name,Ns in (("N3_N16",list(range(3,17))),("N5_N16",[5,16])):
    nc=4 if len(Ns)>2 else 2; nr=int(np.ceil(len(Ns)/nc)); fig,axs=plt.subplots(nr,nc,figsize=(4.2*nc,3.2*nr),squeeze=False)
    for ax,N in zip(axs.flat,Ns):
        for (lab,root),c in zip(ROOTS.items(),COL):
            r=load(root,N)
            if r is None: ax.plot([],[],color=c,label=lab+"（データなし）"); continue
            s,f=r; ax.semilogy(s,np.maximum(f,1e-19),color=c,lw=1.2,label=lab)
        ax.set_title(f"N={N}"+(" (K_sigma raw 枝)" if N==5 else " (decompactification)" if N==16 else "")); ax.set_xlabel("step (2π/144)"); ax.set_ylabel("H⊥/H_total"); ax.set_ylim(3e-20,2); ax.grid(alpha=.3)
    for ax in list(axs.flat)[len(Ns):]: ax.axis("off")
    axs.flat[0].legend(fontsize=7,loc="lower right"); fig.suptitle("H⊥/H_total（差し引き読出し、床≈1e-16）：原本 / 修正版 / baseline / 等モジュラー親"); fig.tight_layout()
    out=os.path.join(H,"results","figures",f"four_way_Hperp_frac_{name}.png"); fig.savefig(out,dpi=130); plt.close(fig); print(out)
