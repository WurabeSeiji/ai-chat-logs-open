# -*- coding: utf-8 -*-
"""取得済データのみで、N=5,8,10,16,20 の H⊥/H_total を 3 色で重ね描き（実験は行わない）。
 図1: L=1000・5000 step（…linear1000_steps5000_equimodular_selfconsistent_directHperp_20260828）
 図2: L=124・40000 step（…linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828）
入力列: data/treatment_linear124_amplitude_aware_timeseries.csv の H_perp, H_total（H⊥ は直交成分の直接計算）。
出力: figures/compare_N_L1000_5000.png, figures/compare_N_L124_40000.png, figures/compare_N_both.png, results.json"""
import os, json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Hiragino Sans","font.size":11})  # 日本語題名の文字化け防止
H=os.path.dirname(os.path.abspath(__file__)); P=os.path.dirname(H)
SETS={"L1000_5000":{N:f"N{N}_linear1000_steps5000_equimodular_selfconsistent_directHperp_20260828" for N in (5,8,10,16,20)},
      "L124_40000":{N:f"N{N}_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828" for N in (5,8,10,16,20)}}
COL={5:"#d7263d",8:"#e67e22",10:"#2e8b57",16:"#1f5fd8",20:"#6c3483"}; out={}
def load(pkg):
    d=pd.read_csv(os.path.join(P,pkg,"data","treatment_linear124_amplitude_aware_timeseries.csv")); return d.step.to_numpy(), (d.H_perp/d.H_total).to_numpy(), d.H_perp.to_numpy()
fig_all,axs=plt.subplots(1,2,figsize=(17,6.5))
for ai,(key,pk) in enumerate(SETS.items()):
    L=1000 if "1000" in key else 124; fig,ax=plt.subplots(figsize=(11,6.5)); out[key]={}
    for N,pkg in pk.items():
        st,f,hp=load(pkg)
        for a in (ax,axs[ai]): a.semilogy(st,np.maximum(f,1e-40),color=COL[N],lw=1.6,label=f"N={N}  (H_total={hp[0]/f[0] if f[0]>0 else np.nan:.3g})" if False else f"N={N}")
        # 数値要約：指数域の傾き（1e-10<f<1e-3 があれば）と、後半（最後の 1/4）の対数傾き
        m=(f>1e-10)&(f<1e-3); q=st>st.max()*0.75
        sl=float(np.polyfit(st[m],np.log(f[m]),1)[0]) if m.sum()>10 else None; sl_late=float(np.polyfit(st[q],np.log(f[q]),1)[0])
        out[key][f"N{N}"]={"f_step1":float(f[1]),"f_final":float(f[-1]),"f_max":float(f.max()),"slope_1e-10..1e-3_per_step":sl,"slope_last_quarter_per_step":sl_late,"slope_last_quarter_per_rad":sl_late*L/(2*np.pi)}
    for a in (ax,axs[ai]):
        a.set_xlabel("step"); a.set_ylabel("H⊥ / H_total  (direct: |Z − P Z|² / |Z|²)"); a.set_title(f"等モジュラー自己無撞着親からの seedless 走行  Δ=2π/{L}, {int(st.max())} step"); a.grid(True,which="both",alpha=.25); a.legend()
    fig.tight_layout(); fig.savefig(os.path.join(H,"figures",f"compare_N_{key}.png"),dpi=160); plt.close(fig)
fig_all.tight_layout(); fig_all.savefig(os.path.join(H,"figures","compare_N_both.png"),dpi=140)
json.dump(out,open(os.path.join(H,"results.json"),"w"),indent=1); print(json.dumps(out,indent=1))
# 図3：横軸を累積刻み角 τ = step·2π/L にして、L=124（実線）と L=1000（破線）を同じ軸に重ねる（同じ N は同色）。刻み依存性と「同じカーブか」の直接比較
fig,ax=plt.subplots(figsize=(11,6.5))
for key,pk in SETS.items():
    L=1000 if "1000" in key else 124; ls="--" if L==1000 else "-"
    for N,pkg in pk.items():
        st,f,_=load(pkg); tau=st*2*np.pi/L; ax.semilogy(tau,np.maximum(f,1e-40),color=COL[N],lw=1.6,ls=ls,label=f"N={N}, Δ=2π/{L}")
ax.set_xlim(0,60); ax.set_xlabel("累積刻み角 τ = step·Δ [rad]"); ax.set_ylabel("H⊥ / H_total"); ax.set_title("同じ τ 軸での比較（実線 Δ=2π/124、破線 Δ=2π/1000）：0〜60 rad"); ax.grid(True,which="both",alpha=.25); ax.legend(ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(H,"figures","compare_N_tau_axis_0-60rad.png"),dpi=160); plt.close(fig)
fig,ax=plt.subplots(figsize=(11,6.5))
for N,pkg in SETS["L124_40000"].items():
    st,f,_=load(pkg); ax.semilogy(st*2*np.pi/124,np.maximum(f,1e-40),color=COL[N],lw=1.6,label=f"N={N}, Δ=2π/124")
ax.set_xlabel("累積刻み角 τ [rad]"); ax.set_ylabel("H⊥ / H_total"); ax.set_title("τ 軸での全区間（Δ=2π/124、40000 step = 2027 rad）"); ax.grid(True,which="both",alpha=.25); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(H,"figures","compare_N_tau_axis_full.png"),dpi=160); plt.close(fig)
# 数値：同じ τ での比較（τ=31.4 rad ＝ L=1000 の 5000 step ＝ L=124 の 620 step）
cmp={}
for N in (5,8,10,16,20):
    s1,f1,_=load(SETS["L1000_5000"][N]); s2,f2,_=load(SETS["L124_40000"][N]); cmp[f"N{N}"]={"f(τ=31.4) L=1000":float(f1[5000]),"f(τ=31.4) L=124(step620)":float(f2[620]),"ratio L124/L1000":float(f2[620]/f1[5000]),"(1000/124)^2":float((1000/124)**2)}
out["same_tau_comparison"]=cmp; json.dump(out,open(os.path.join(H,"results.json"),"w"),indent=1); print(json.dumps(cmp,indent=1))
