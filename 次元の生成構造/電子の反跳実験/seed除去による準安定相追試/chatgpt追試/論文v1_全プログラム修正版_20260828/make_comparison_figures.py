# -*- coding: utf-8 -*-
"""原本（旧エンジン）と修正版（FIX1-4）の主要量を重ね描きする比較図。出力 results/figures/*.png"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Hiragino Sans","font.size":11})
H=os.path.dirname(os.path.abspath(__file__)); O=os.path.join(os.path.dirname(H),"論文v1_全再現テスト_20260828","original"); F=os.path.join(H,"fixed"); OUT=os.path.join(H,"results","figures")
def rd(base,p): return pd.read_csv(os.path.join(base,p))
# 1) decompactification: H_perp（対数）N5/N16
fig,ax=plt.subplots(1,2,figsize=(13,5))
for i,N in enumerate((5,16)):
    a=rd(O,f"complex_simplex_decompactification_N5_N16_20260826/results/N{N}_geometry_summary.csv"); b=rd(F,f"complex_simplex_decompactification_N5_N16_20260826/results/N{N}_geometry_summary.csv")
    ax[i].semilogy(a.step,np.maximum(a.H_perp,1e-40),color="#c0392b",label="原本（Cayley・位相のみK・正規化）H⊥"); ax[i].semilogy(a.step,a.H_parallel,color="#c0392b",ls="--",alpha=.6,label="原本 H∥")
    ax[i].semilogy(b.step,np.maximum(b.H_perp,1e-40),color="#1f4e79",label="修正版（線形回転・振幅込みK・正規化なし）H⊥"); ax[i].semilogy(b.step,b.H_parallel,color="#1f4e79",ls="--",alpha=.6,label="修正版 H∥")
    ax[i].set_title(f"N={N}：親平面外成分 H⊥（対数）"); ax[i].set_xlabel("step"); ax[i].set_ylim(1e-34,10); ax[i].legend(fontsize=8,loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"cmp1_Hperp_log_N5_N16.png"),dpi=160); plt.close(fig)
# 2) pump depletion（K_sigma raw N5）
fig,ax=plt.subplots(figsize=(9,4.8)); a=rd(O,"K_sigma_normalization_artifact_test_N4_N5_20260826/N5_raw_K_raw_observables.csv"); b=rd(F,"K_sigma_normalization_artifact_test_N4_N5_20260826/N5_raw_K_raw_observables.csv")
ax.plot(a.step,a.H_parallel,color="#c0392b",ls="--",label="原本 H∥"); ax.plot(a.step,a.H_perp,color="#c0392b",label="原本 H⊥"); ax.plot(b.step,b.H_parallel,color="#1f4e79",ls="--",label="修正版 H∥"); ax.plot(b.step,b.H_perp,color="#1f4e79",label="修正版 H⊥"); ax.plot(b.step,b.H_total,color="k",lw=.8,label="修正版 H_total（=c²）")
ax.set_title("N=5 raw K：親平面 ↔ 直交方向の移送"); ax.set_xlabel("step"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(OUT,"cmp2_pump_depletion_N5.png"),dpi=160); plt.close(fig)
# 3) tol sweep
fig,ax=plt.subplots(1,2,figsize=(12,4.6)); 
for base,col,lab in ((O,"#c0392b","原本"),(F,"#1f4e79","修正版")):
    t=rd(base,"N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826/tol_sweep_timeseries.csv")
    for k,tol in enumerate(sorted(t.tol.unique())): d=t[t.tol==tol]; ax[0].semilogy(d.step,np.maximum(d.f,1e-40),color=col,alpha=.4+.15*k,lw=1,label=f"{lab} tol={tol:g}")
    s=rd(base,"N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826/tol_sweep_summary.csv"); ax[1].plot(-np.log(s.parent_residual),s["onset_f_ge_1e-8"],"o",color=col,ms=9,label=lab)
ax[0].set_title("seedless：f=H⊥/H の時間発展（tol 4 段階）"); ax[0].set_xlabel("step"); ax[0].set_ylim(1e-34,2); ax[0].legend(fontsize=7,ncol=2)
ax[1].set_title("onset step vs −ln(親残差)"); ax[1].set_xlabel("−ln ε"); ax[1].set_ylabel("onset step"); ax[1].legend(); fig.tight_layout(); fig.savefig(os.path.join(OUT,"cmp3_tol_sweep_onset.png"),dpi=160); plt.close(fig)
# 4) Floquet circle
fig,ax=plt.subplots(1,2,figsize=(11,5.2))
for i,(base,lab) in enumerate(((O,"原本"),(F,"修正版"))):
    f=rd(base,"N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826/floquet_spectrum.csv"); f=f[np.isclose(f.fd_eps,1e-7)]
    th=np.linspace(0,2*np.pi,400); ax[i].plot(np.cos(th),np.sin(th),"k--",lw=.8); ax[i].scatter(f.eig_re,f.eig_im,s=40,color="#c0392b" if i==0 else "#1f4e79"); ax[i].set_aspect("equal"); ax[i].set_title(f"{lab}：回転系 1 step Jacobian の固有値（max|μ|={f.modulus.max():.6f}）"); ax[i].set_xlim(-1.2,1.2); ax[i].set_ylim(-1.2,1.2)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"cmp4_floquet_spectrum.png"),dpi=160); plt.close(fig)
# 5) amplitude equalization: 修正版の |z_k|² 時間発展（N5 physical phase test の phase_by_edge に振幅列あり）と原本
fig,ax=plt.subplots(1,2,figsize=(12,4.6))
for i,(base,lab,col) in enumerate(((O,"原本","#c0392b"),(F,"修正版","#1f4e79"))):
    p=rd(base,"N5_complex_simplex_complete_analysis_20260826/N5_phase_by_edge_5000steps.csv"); 
    for e in range(10): d=p[p.edge_index==e]; ax[i].plot(d.step,d.amplitude**2,lw=.8,color=col,alpha=.7)
    amp2=p.pivot(index="step",columns="edge_index",values="amplitude")**2; pr=amp2.div(amp2.sum(axis=1),axis=0); S=-(pr*np.log(pr.where(pr>0,1))).sum(axis=1)/np.log(10)
    ax[i].set_title(f"{lab}：10 本の |z_k|²（最終 S/lnM={S.iloc[-1]:.4f}）"); ax[i].set_xlabel("step")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"cmp5_amplitudes_N5.png"),dpi=160); plt.close(fig)
# 6) N=3..16 final |z|² spread and H_perp/H_total
rows=[]
pk={3:"N3_N4",4:"N3_N4",6:"N6_N7",7:"N6_N7",8:"N8_N9",9:"N8_N9",10:"N10_N11",11:"N10_N11",12:"N12_N13",13:"N12_N13",14:"N14_N15",15:"N14_N15"}
for N,pf in pk.items():
    for base,lab in ((O,"原本"),(F,"修正版")):
        g=rd(base,f"{pf}_complex_simplex_complete_analysis_20260826/N{N}_global_summary.csv"); r=g.iloc[-1]; rows.append((N,lab,(r.r2_max-r.r2_min)/((r.r2_max+r.r2_min)/2),r.H_perp/r.H_total,r.simplex_rank))
df=pd.DataFrame(rows,columns=["N","run","r2_rel_spread","Hperp_frac","rank"])
fig,ax=plt.subplots(1,2,figsize=(12,4.4))
for lab,col in (("原本","#c0392b"),("修正版","#1f4e79")):
    d=df[df.run==lab]; ax[0].semilogy(d.N,np.maximum(d.r2_rel_spread,1e-16),"o-",color=col,label=lab); ax[1].plot(d.N,d.Hperp_frac,"o-",color=col,label=lab)
ax[0].set_title("step 5000 の |z|² 相対ばらつき（0=等分配）"); ax[0].set_xlabel("N"); ax[0].legend(); ax[1].set_title("step 5000 の H⊥/H_total"); ax[1].set_xlabel("N"); ax[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(OUT,"cmp6_N_series_final_state.png"),dpi=160); plt.close(fig); df.to_csv(os.path.join(OUT,"cmp6_N_series_final_state.csv"),index=False)
print("figures:",sorted(os.listdir(OUT)))
