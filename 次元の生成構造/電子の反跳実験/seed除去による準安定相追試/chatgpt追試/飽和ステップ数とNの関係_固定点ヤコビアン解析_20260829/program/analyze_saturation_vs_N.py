# -*- coding: utf-8 -*-
"""飽和 step 数と N の関係——相対平衡の共回転モノドロミーによる定式化。

等モジュラー自己無撞着親 v は 1 step 写像 F(z)=exp((2π/124)K_amp(z))·z の相対平衡: F(v)=e^{iφ}v, φ=ANGLE·σ_max。
線形安定性は DF(v) ではなく共回転モノドロミー G = R(−φ)·DF(v)（R は全体位相回転の実 2M 表現）で決まる。
振幅の成長率 λ_G = ln|μ_G|max、H⊥/H_total（振幅の 2 乗）の成長率 λ_f = 2·λ_G。
飽和 step: t50 ≈ (ln 0.5 − ln f_seed)/λ_f、f_seed は指数窓を step 0 に外挿した切片（丸め誤差床 ~1e-20）。

(1) 実測 7 走行（N=5,6,7,8,10,16,20）の指数窓成長率 λ_emp と、その走行の親 Z[0] で計算した 2λ_G の突合
(2) N=4..20 × 親 5 実現（seed）の λ_G 分布
(3) t50 予測と実測の比較、図
出力: data/empirical_growth.csv, data/lambda_G_at_run_parent.csv, data/lambda_G_vs_N_seeds.csv, data/comparison.csv, figures/lambda_and_t50_vs_N.png"""
import os, sys, json, math, importlib.util, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"]=["Hiragino Sans","Arial Unicode MS","DejaVu Sans"]
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); UP=os.path.dirname(ROOT)
spec=importlib.util.spec_from_file_location("eng",os.path.join(HERE,"original_engine.py")); eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress=lambda m: None
ANGLE=eng.ANGLE
RUNS={5:"N5_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828",6:"N6_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260829",7:"N7_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260829",
      8:"N8_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828",10:"N10_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828",16:"N16_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828",20:"N20_linear124_equimodular_selfconsistent_directHperp_treatment_only_20260828"}
def monodromy_lambda(N,v,h=1e-6):
    s=eng.LowRankSystem(N); A=eng._adjacency(s); M=len(v)
    F=lambda z: eng._exp_step(eng._K_amplitude_aware(A,z),z,ANGLE)
    Fv=F(v); phi=float(np.angle(np.vdot(v,Fv))); defect=float(np.linalg.norm(Fv-np.exp(1j*phi)*v)/np.linalg.norm(v))
    to_r=lambda z: np.concatenate([z.real,z.imag]); J=np.zeros((2*M,2*M))
    for k in range(2*M):
        e=np.zeros(2*M); e[k]=h; d=e[:M]+1j*e[M:]; J[:,k]=(to_r(F(v+d))-to_r(F(v-d)))/(2*h)
    c,sn=np.cos(-phi),np.sin(-phi); R=np.block([[c*np.eye(M),-sn*np.eye(M)],[sn*np.eye(M),c*np.eye(M)]]); G=R@J
    mu=np.linalg.eigvals(G); a=np.sort(np.abs(mu))[::-1]
    return dict(phi=phi,relative_equilibrium_defect=defect,mu_G_max=float(a[0]),lambda_G=float(np.log(a[0])),lambda_f_pred=float(2*np.log(a[0])),n_unstable=int((np.abs(mu)>1+1e-9).sum()),mu_top4=", ".join(f"{x:.6f}" for x in a[:4]),lambda_DF_max=float(np.log(np.abs(np.linalg.eigvals(J)).max())))
# ---------- (1) 実測と、その親での 2λ_G
emp=[]
for N,d in RUNS.items():
    t=pd.read_csv(os.path.join(UP,d,"data","treatment_linear124_amplitude_aware_timeseries.csv")); s=t.step.to_numpy(); f=(t.H_perp/t.H_total).to_numpy(); j=json.load(open(os.path.join(UP,d,"data","summary.json"))); r=j.get("parent_residual",float("nan"))
    on=lambda th: int(s[f>=th][0]) if (f>=th).any() else None
    tau=s*ANGLE; ball=(r*tau)**2; m=(f>100*np.maximum(ball,1e-300))&(f<1e-3)&(s>0)
    c=np.polyfit(s[m],np.log(f[m]),1) if m.sum()>=20 else (float("nan"),float("nan"))
    late=s>=s[-1]-10000; lam_late=float(np.polyfit(s[late],np.log(np.maximum(f[late],1e-300)),1)[0])
    Z0=np.load(os.path.join(UP,d,"data","states_treatment.npz"))["Z"][0]
    mono=monodromy_lambda(N,Z0)
    emp.append(dict(N=N,M=N*(N-1)//2,parent_residual=r,steps=int(s[-1]),f_final=float(f[-1]),f_max=float(f.max()),onset_1e8=on(1e-8),step_50pct=on(0.5),saturated=bool(f.max()>=0.5),lambda_emp=float(c[0]),ln_f_seed=float(c[1]),n_window=int(m.sum()),lambda_last10000=lam_late,**{k:mono[k] for k in ("lambda_f_pred","lambda_G","n_unstable","relative_equilibrium_defect","lambda_DF_max","mu_top4")}))
    print(f"N={N}: λ_emp={c[0]:.4e}  2λ_G(run parent)={mono['lambda_f_pred']:.4e}  末尾傾き={lam_late:.2e}  ln f_seed={c[1]:.1f}  t50={on(0.5)}",flush=True)
E=pd.DataFrame(emp); E.to_csv(os.path.join(ROOT,"data","empirical_growth.csv"),index=False)
# ---------- (2) N=4..20 × 親 5 実現
rows=[]
for N in range(4,21):
    for seed in range(5):
        s=eng.LowRankSystem(N); rng=np.random.default_rng(40260721+1000*N+seed)
        try: v,res,sig=eng.make_parent(s,rng)
        except RuntimeError as ex: rows.append(dict(N=N,seed=seed,M=N*(N-1)//2,parent_residual=float("nan"),lambda_G=float("nan"))); print(f"N={N} seed={seed}: parent FAIL"); continue
        mono=monodromy_lambda(N,v); rows.append(dict(N=N,seed=seed,M=N*(N-1)//2,parent_residual=float(res),sigma_max=float(sig[0]),**mono))
    sub=[r["lambda_G"] for r in rows if r["N"]==N and not np.isnan(r["lambda_G"])]; print(f"N={N}: 2λ_G over seeds = {np.round(2*np.array(sub),5)}",flush=True)
S=pd.DataFrame(rows); S.to_csv(os.path.join(ROOT,"data","lambda_G_vs_N_seeds.csv"),index=False)
# ---------- (3) 予測
LNF=float(E.ln_f_seed.mean()); agg=S.groupby("N").agg(lambda_f_min=("lambda_G",lambda x: 2*x.min()),lambda_f_med=("lambda_G",lambda x: 2*x.median()),lambda_f_max=("lambda_G",lambda x: 2*x.max())).reset_index()
agg["t50_pred_med"]=(np.log(0.5)-LNF)/agg.lambda_f_med; agg["t50_pred_min"]=(np.log(0.5)-LNF)/agg.lambda_f_max; agg["t50_pred_max"]=(np.log(0.5)-LNF)/agg.lambda_f_min
C=agg.merge(E[["N","lambda_emp","lambda_f_pred","step_50pct","saturated","lambda_last10000"]],on="N",how="left"); C["t50_pred_run_parent"]=(np.log(0.5)-LNF)/C.lambda_f_pred
C.to_csv(os.path.join(ROOT,"data","comparison.csv"),index=False); print("mean ln f_seed =",LNF); print(C.to_string())
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for N,g in S.groupby("N"): ax[0].semilogy([N]*len(g),2*g.lambda_G,"o",c="#1f77b4",ms=4,alpha=.6)
ax[0].semilogy(agg.N,agg.lambda_f_med,"-",c="#1f77b4",label="2λ_G 中央値（親 5 実現）"); ax[0].semilogy(E.N,E.lambda_f_pred,"D",c="#2ca02c",ms=8,mfc="none",label="2λ_G（走行の親）"); ok=E.lambda_emp.notna(); ax[0].semilogy(E.N[ok],E.lambda_emp[ok],"s",c="#d62728",ms=9,mfc="none",label="λ_emp（実測 指数窓）"); ax[0].semilogy(E.N,E.lambda_last10000,"^",c="gray",ms=6,mfc="none",label="末尾 10000 step の傾き")
ax[0].set_xlabel("N"); ax[0].set_ylabel("H⊥/H_total の成長率 /step（2π/124）"); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8); ax[0].set_title("線形成長率：共回転モノドロミー vs 実測")
ax[1].fill_between(agg.N,agg.t50_pred_min,agg.t50_pred_max,alpha=.2,label="予測 t50 範囲（親 5 実現）"); ax[1].semilogy(agg.N,agg.t50_pred_med,"-",label="予測 t50 中央値"); ax[1].semilogy(C.N,C.t50_pred_run_parent,"D",c="#2ca02c",ms=8,mfc="none",label="予測 t50（走行の親）"); sat=C.saturated==True; ax[1].semilogy(C.N[sat],C.step_50pct[sat],"s",c="#d62728",ms=9,mfc="none",label="実測 50% 到達 step"); ax[1].axhline(40000,ls="--",c="gray",lw=.8); ax[1].text(4.2,44000,"40000",fontsize=8,color="gray")
ax[1].set_xlabel("N"); ax[1].set_ylabel("step"); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8); ax[1].set_title(f"飽和 step 予測 t50=(ln0.5−ln f_seed)/λ_f, ln f_seed={LNF:.1f}")
fig.tight_layout(); fig.savefig(os.path.join(ROOT,"figures","lambda_and_t50_vs_N.png"),dpi=150); print("figure saved")
