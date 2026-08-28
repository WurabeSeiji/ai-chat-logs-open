# -*- coding: utf-8 -*-
"""物理量レベルの突合：軌道そのもの（丸め誤差で発散し位相は巻き付く）でなく、論文が主張する不変量・構造量を original / rerun で並べる。出力: results/physics_comparison.md"""
import os, json, numpy as np, pandas as pd
H=os.path.dirname(os.path.abspath(__file__)); O=os.path.join(H,"original"); R=os.path.join(H,"rerun"); L=[]
def P(s=""): L.append(s)
def both(pkg,f): return pd.read_csv(os.path.join(O,pkg,f)), pd.read_csv(os.path.join(R,pkg,f))
def bothj(pkg,f): return json.load(open(os.path.join(O,pkg,f))), json.load(open(os.path.join(R,pkg,f)))
def flat(o,p=""):
    if isinstance(o,dict):
        for k,v in o.items(): yield from flat(v,p+"/"+str(k))
    elif isinstance(o,list) and o and isinstance(o[0],dict):
        for i,v in enumerate(o): yield from flat(v,p+f"[{i}]")
    else: yield p,o
P("# 物理量レベル突合（original vs rerun）"); P()
# ---- N=3..16 complete analysis
pk={3:"N3_N4",4:"N3_N4",6:"N6_N7",7:"N6_N7",8:"N8_N9",9:"N8_N9",10:"N10_N11",11:"N10_N11",12:"N12_N13",13:"N12_N13",14:"N14_N15",15:"N14_N15"}
P("## 1. N=3〜15 complete analysis：summary.json の全キー（original → rerun）"); P()
for N,pf in pk.items():
    pkg=f"{pf}_complex_simplex_complete_analysis_20260826"; a,b=bothj(pkg,f"N{N}_summary.json"); A=dict(flat(a)); B=dict(flat(b))
    rows=[]
    for k,v in A.items():
        w=B.get(k); 
        if isinstance(v,(int,float)) and not isinstance(v,bool): rows.append(f"{k}: {v:.6g} → {w:.6g}" + ("" if abs(v-w)<1e-9*max(1,abs(v)) else f"  (Δ={w-v:.3g})"))
        elif v!=w: rows.append(f"{k}: {v} → {w}")
    P(f"### N={N}"); P("- "+"\n- ".join(rows)); P()
P("## 2. N=3〜16 global_summary：保存量と最終状態"); P(); P("| N | H_total drift (o/r) | final |ZᵀZ| (o/r) | final r2 rel spread (o/r) | final simplex_rank (o/r) | H_perp onset step ≥1e-8 (o/r) |"); P("|---|---|---|---|---|---|")
for N in list(pk)+[16]:
    pkg=f"{pk[N]}_complex_simplex_complete_analysis_20260826" if N!=16 else "N16_complex_simplex_complete_analysis_20260826"; a,b=both(pkg,f"N{N}_global_summary.csv")
    def st(d):
        drift=float(d.H_total.max()-d.H_total.min()); z=float(d.abs_ZtZ.iloc[-1]); sp=float((d.r2_max.iloc[-1]-d.r2_min.iloc[-1])/((d.r2_max.iloc[-1]+d.r2_min.iloc[-1])/2)); rk=int(d.simplex_rank.iloc[-1])
        on=int(d.step[d.H_perp/d.H_total>=1e-8].iloc[0]) if "H_perp" in d and (d.H_perp/d.H_total>=1e-8).any() else None; return drift,z,sp,rk,on
    x=st(a); y=st(b); P(f"| {N} | {x[0]:.1e} / {y[0]:.1e} | {x[1]:.1e} / {y[1]:.1e} | {x[2]:.1e} / {y[2]:.1e} | {x[3]} / {y[3]} | {x[4]} / {y[4]} |")
P(); P("## 3. time_milestones（秩序化の到達 step、original → rerun）"); P()
for N,pf in pk.items():
    pkg=f"{pf}_complex_simplex_complete_analysis_20260826"; a,b=both(pkg,f"N{N}_time_milestones.csv"); m=a.merge(b,on="metric",suffixes=("_o","_r"))
    P(f"- N={N}: "+", ".join(f"{r.metric}={r.step_o:g}→{r.step_r:g}" for r in m.itertuples()))
P(); P("## 4. 最終クラス数・等分配（final_classes / step5000_final_edges）"); P(); P("| N | classes tol1e-8 (o/r) | |z|² min..max original | |z|² min..max rerun |"); P("|---|---|---|---|")
for N,pf in pk.items():
    pkg=f"{pf}_complex_simplex_complete_analysis_20260826"; fn=f"N{N}_final_classes_tol1e-8.csv" if N<5 else f"N{N}_final_classes_tol1e-08.csv"; a,b=both(pkg,fn); e,f=both(pkg,f"N{N}_step5000_final_edges.csv")
    c=[c for c in e.columns if c in("r2","abs_z2","z2_abs")]; 
    if c: P(f"| {N} | {len(a)} / {len(b)} | {e[c[0]].min():.10f}..{e[c[0]].max():.10f} | {f[c[0]].min():.10f}..{f[c[0]].max():.10f} |")
    else: P(f"| {N} | {len(a)} / {len(b)} | cols={list(e.columns)[:8]} | |")
# ---- followup
P(); P("## 5. N5 followup（記事図 2〜6・8 の源）"); P()
pkg="N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826"
a,b=both(pkg,"tol_sweep_summary.csv"); P("### tol sweep（記事図3）"); P("original:\n```\n"+a.to_string()+"\n```\nrerun:\n```\n"+b.to_string()+"\n```")
a,b=bothj(pkg,"tol_sweep_fit.json"); P(f"fit original: {a}"); P(f"fit rerun: {b}"); P()
a,b=both(pkg,"floquet_spectrum.csv"); P("### Floquet（記事図4）上位乗数 |μ|（fd_eps ごと）")
for eps in sorted(a.fd_eps.unique()):
    x=np.sort(a[np.isclose(a.fd_eps,eps)].modulus.to_numpy())[::-1][:4]; y=np.sort(b[np.isclose(b.fd_eps,eps)].modulus.to_numpy())[::-1][:4]
    P(f"- eps={eps:g}: original {np.array2string(x,precision=9)}  rerun {np.array2string(y,precision=9)}")
a,b=bothj(pkg,"floquet_summary.json"); P(f"floquet_summary original: {a}"); P(f"floquet_summary rerun: {b}"); P()
a,b=both(pkg,"N5_moduli_seed_sweep.csv"); P("### moduli seed sweep（記事図8）"); cols=[c for c in a.columns if c in("seed","parent_residual","counts","family1_mod","family2_mod","relative_phase_mod_pi_rad","cluster_sse")]
P("original:\n```\n"+a[cols].to_string()+"\n```\nrerun (run_moduli_sweep_fast):\n```\n"+b[cols].to_string()+"\n```")
c=pd.read_csv(os.path.join(R,pkg,"N5_moduli_seed_sweep_from_run_followup_20seeds.csv")); cc=[x for x in c.columns if x in("seed","parent_residual","counts","cluster_sse","family_A_modulus","family_B_modulus","relative_phase_mod_pi_rad")]; P("rerun (run_followup_experiments, 20 seeds):\n```\n"+c[cc].to_string()+"\n```"); P()
a,b=both(pkg,"N5_spectral_entropy_timeseries.csv"); P(f"### spectral entropy（記事図6）: final original {a.iloc[-1].to_dict()} / rerun {b.iloc[-1].to_dict()}")
a,b=bothj(pkg,"pump_depletion_summary.json"); P(f"### pump depletion（記事図2）: original {a} / rerun {b}"); P()
# ---- K_sigma
P("## 6. K/σ 正規化比較（N4, N5）：H_perp 成長率と onset"); P()
pkg="K_sigma_normalization_artifact_test_N4_N5_20260826"
for N in (4,5):
    for br in ("normalized","raw"):
        a,b=both(pkg,f"N{N}_{br}_K_raw_observables.csv")
        def g(d):
            f=(d.H_perp/d.H_total).to_numpy(); s=d.step.to_numpy(); on=int(s[f>=1e-8][0]) if (f>=1e-8).any() else None
            m=(f>1e-6)&(f<1e-2); sl=np.polyfit(s[m],np.log(f[m]),1)[0] if m.sum()>3 else float("nan"); return on,sl,float(f[-1])
        x=g(a); y=g(b); P(f"- N={N} {br}: onset {x[0]}→{y[0]}, growth rate {x[1]:.6f}→{y[1]:.6f} /step, final H_perp/H {x[2]:.4f}→{y[2]:.4f}")
a=json.load(open(os.path.join(O,pkg,"summary.json"))); P(f"- original summary.json: {a}"); P()
# ---- decompactification
P("## 7. decompactification N5/N16（記事図1の源）summary.json"); P()
pkg="complex_simplex_decompactification_N5_N16_20260826"
for N in (5,16):
    a,b=bothj(pkg,f"results/N{N}_summary.json"); A=dict(flat(a)); B=dict(flat(b))
    P(f"### N={N}"); P("- "+"\n- ".join(f"{k}: {v} → {B.get(k)}" for k,v in A.items() if not (isinstance(v,(int,float)) and isinstance(B.get(k),(int,float)) and abs(v-B[k])<1e-9*max(1,abs(v))))); P()
    g,h=both(pkg,f"results/N{N}_geometry_summary.csv"); P(f"- H_perp: start {g.H_perp.iloc[0]:.3e}→{h.H_perp.iloc[0]:.3e}, min {g.H_perp.min():.3e}→{h.H_perp.min():.3e}, final {g.H_perp.iloc[-1]:.4f}→{h.H_perp.iloc[-1]:.4f}, max {g.H_perp.max():.4f}→{h.H_perp.max():.4f}; onset(≥1e-8) {int(g.step[g.H_perp>=1e-8].iloc[0])}→{int(h.step[h.H_perp>=1e-8].iloc[0])}"); P()
# ---- gamma
P("## 8. γ 掃引（論文§12）"); pkg="N5_gamma_continuum_test_bundle_20260825"
for f in ("N5_gamma_continuum_rates.csv","N5_gamma_continuum_stats.csv"): a,b=both(pkg,f); P(f"### {f}\noriginal:\n```\n{a.to_string()}\n```\nrerun:\n```\n{b.to_string()}\n```")
# ---- N16 physics + N5 physical phase
P(); P("## 9. N16 physics 標準出力 JSON（log）"); P("original COMBINED_ANALYSIS の値は md 参照。rerun log 末尾:"); P("```"); P(open(os.path.join(H,"results/logs/N16_complex_simplex_complete_analysis_20260826.log")).read()[-900:]); P("```")
P(); P("## 10. N5 physical phase step test：metrics 最終行"); pkg="N5_complex_simplex_complete_analysis_20260826"; a,b=both(pkg,"N5_physical_phase_metrics_5000steps.csv"); P(f"original: {a.iloc[-1].to_dict()}"); P(f"rerun: {b.iloc[-1].to_dict()}")
# ---- closure search
P(); P("## 11. N14–16 閉包探索（k=2..4 exact, k=5,6 MITM）"); pkg="N14_N16_complete_nontrivial_zero_closure_search_20260826"
try:
    a=pd.read_csv(os.path.join(O,pkg,"N14_N16_all_cardinalities_best_results.csv"))
    for N in (14,15,16):
        for f in (f"rerun_exact_k234_N{N}.csv",f"rerun_mitm56_N{N}.csv",f"rerun_search_subsets_N{N}.csv"):
            p=os.path.join(R,pkg,f)
            if os.path.exists(p):
                b=pd.read_csv(p)
                for r in b.itertuples():
                    o=a[(a.N==N)&(a.k==r.k)]; ov=float(o.best_residual.iloc[0]) if len(o) else float("nan"); P(f"- N={N} k={r.k} [{f.split('_')[1]}]: original best {ov:.6e}, rerun {r.best_residual:.6e}, edges rerun {r.edges}")
except Exception as e: P(f"(closure search: {e})")
open(os.path.join(H,"results","physics_comparison.md"),"w",encoding="utf-8").write("\n".join(L)+"\n"); print("\n".join(L))
