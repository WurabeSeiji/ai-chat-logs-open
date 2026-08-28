# -*- coding: utf-8 -*-
"""原本（旧エンジン：Cayley・位相のみ K・正規化あり）／修正版 treatment（線形回転・振幅込み K・正規化なし）／修正版 baseline（線形回転・位相のみ K・正規化なし）の
／等モジュラー（修正版と同じ力学だが親を等モジュラー自己無撞着親にしたもの）の物理量四者比較。出力: results/four_way_comparison.md, results/four_way_comparison.json"""
import os, json, numpy as np, pandas as pd
H=os.path.dirname(os.path.abspath(__file__)); O=os.path.join(os.path.dirname(H),"論文v1_全再現テスト_20260828","original"); F=os.path.join(H,"fixed"); B=os.path.join(H,"fixed_baseline"); E=os.path.join(H,"fixed_equimodular")
ROOTS={"原本":O,"修正版":F,"baseline":B,"等モジュラー":E}; L=[]; J={}
def P(s=""): L.append(s)
def rd(root,pkg,f):
    p=os.path.join(root,pkg,f); return pd.read_csv(p) if os.path.exists(p) else None
def js(root,pkg,f):
    p=os.path.join(root,pkg,f); return json.load(open(p)) if os.path.exists(p) else None
def fmt(x):
    if x is None: return "–"
    if isinstance(x,float): return f"{x:.6g}"
    return str(x)
def row(name,vals): P(f"| {name} | "+" | ".join(fmt(v) for v in vals)+" |")
P("# 四者比較：原本 / 修正版（振幅込み K） / baseline（位相のみ K・線形回転） / 等モジュラー（振幅込み K＋等モジュラー自己無撞着親）"); P(); P("列の意味：原本＝公開エンジン（Cayley γ=tan(π/144)、位相のみ K、親正規化あり）、修正版＝判断 1〜6 適用（exp((2π/144)K)、K=Im(z̄ᵢzⱼ)、正規化なし）、baseline＝同じく判断適用だが K は位相のみ（回転の修正だけの効果を分離する対照）、等モジュラー＝修正版と同じ力学で make_parent だけを 3 段階等モジュラー自己無撞着親（段階2 は 2π/124 の位相のみ線形回転、段階3 は振幅込み K の −σ_max モード磨き、残差<1e-10、段階2 の step 上限 40000）に置換したもの。N=3 は段階2 に約 29000 step を要した（他は 500〜7500 step）。"); P()
# ---- N=3..16
pk={3:"N3_N4",4:"N3_N4",6:"N6_N7",7:"N6_N7",8:"N8_N9",9:"N8_N9",10:"N10_N11",11:"N10_N11",12:"N12_N13",13:"N12_N13",14:"N14_N15",15:"N14_N15",16:"N16"}
P("## 1. N=3〜16：急拡大・保存量・最終状態（global_summary / summary.json）"); P()
P("| N | 量 | 原本 | 修正版 | baseline | 等モジュラー |"); P("|---|---|---|---|---|---|")
for N,pf in pk.items():
    pkg=f"{pf}_complex_simplex_complete_analysis_20260826"; J[f"N{N}"]={}
    def stats(root):
        d=rd(root,pkg,f"N{N}_global_summary.csv")
        if d is None: return None
        Hp=d.H_perp if "H_perp" in d else None; f=(Hp/d.H_total) if Hp is not None else None
        on=int(d.step[f>=1e-8].iloc[0]) if f is not None and (f>=1e-8).any() else None
        on5=int(d.step[f>=0.05].iloc[0]) if f is not None and (f>=0.05).any() else None
        # 成長率：f∈(1e-10,1e-3) の対数傾き
        gr=None
        if f is not None:
            m=(f>1e-10)&(f<1e-3); 
            if m.sum()>=10: gr=float(np.polyfit(d.step[m],np.log(f[m]),1)[0])
        return dict(H_total_0=float(d.H_total.iloc[0]),H_total_drift=float(d.H_total.max()-d.H_total.min()),ZtZ_final=float(d.abs_ZtZ.iloc[-1]),
            Hperp_frac_step0=float(f.iloc[0]) if f is not None else None,Hperp_frac_min=float(f.min()) if f is not None else None,Hperp_frac_max=float(f.max()) if f is not None else None,Hperp_frac_final=float(f.iloc[-1]) if f is not None else None,
            onset_1e8=on,onset_5pct=on5,growth_rate=gr,r2_rel_spread_final=float((d.r2_max.iloc[-1]-d.r2_min.iloc[-1])/((d.r2_max.iloc[-1]+d.r2_min.iloc[-1])/2)),r2_min_final=float(d.r2_min.iloc[-1]),r2_max_final=float(d.r2_max.iloc[-1]),rank_final=int(d.simplex_rank.iloc[-1]))
    S={k:stats(r) for k,r in ROOTS.items()}; J[f"N{N}"]["global"]=S
    for key,lab in [("H_total_0","H_total（=|Z₀|²、保存）"),("H_total_drift","H_total ドリフト"),("ZtZ_final","final |ZᵀZ|"),("Hperp_frac_step0","H⊥/H step0"),("Hperp_frac_max","H⊥/H 最大"),("Hperp_frac_final","H⊥/H final"),("onset_1e8","onset(H⊥/H≥1e-8)"),("onset_5pct","onset(≥5%)"),("growth_rate","指数成長率 /step（1e-10<f<1e-3）"),("r2_rel_spread_final","final |z|² 相対幅"),("r2_min_final","final |z|² min"),("r2_max_final","final |z|² max"),("rank_final","final rank")]:
        row(f"{N} | {lab}" if key=="H_total_0" else f"  | {lab}",[S[k][key] if S[k] else None for k in ROOTS])
P(); P("## 2. N=3〜15：クラス数と秩序化到達 step（final_classes / time_milestones）"); P(); P("| N | 量 | 原本 | 修正版 | baseline | 等モジュラー |"); P("|---|---|---|---|---|---|")
for N,pf in pk.items():
    if N==16: continue
    pkg=f"{pf}_complex_simplex_complete_analysis_20260826"; fn=f"N{N}_final_classes_tol1e-8.csv" if N<5 else f"N{N}_final_classes_tol1e-08.csv"
    cls=[len(rd(r,pkg,fn)) if rd(r,pkg,fn) is not None else None for r in ROOTS.values()]; row(f"{N} | クラス数 tol1e-8（M={N*(N-1)//2}）",cls)
    ms={k:rd(r,pkg,f"N{N}_time_milestones.csv") for k,r in ROOTS.items()}
    if all(m is not None for m in ms.values()):
        for met in ms["原本"].metric:
            vals=[float(m.step[m.metric==met].iloc[0]) if (m.metric==met).any() else None for m in ms.values()]; vals=[None if (v is None or np.isnan(v)) else int(v) for v in vals]; row(f"  | {met}",vals)
    J[f"N{N}"]["classes"]=cls
# ---- decompactification
P(); P("## 3. decompactification N5 / N16（記事図 1 の源）"); P(); P("| N | 量 | 原本 | 修正版 | baseline | 等モジュラー |"); P("|---|---|---|---|---|---|")
pkg="complex_simplex_decompactification_N5_N16_20260826"
for N in (5,16):
    G={k:rd(r,pkg,f"results/N{N}_geometry_summary.csv") for k,r in ROOTS.items()}; S={k:js(r,pkg,f"results/N{N}_summary.json") for k,r in ROOTS.items()}
    def g(d,key):
        if d is None: return None
        c=d[key]; return dict(start=float(c.iloc[0]),min=float(c.min()),max=float(c.max()),final=float(c.iloc[-1]),onset=int(d.step[c>=1e-8].iloc[0]) if (c>=1e-8).any() else None)
    HP={k:g(d,"H_perp") for k,d in G.items()}
    row(f"{N} | H⊥ start",[HP[k]["start"] if HP[k] else None for k in ROOTS]); row("  | H⊥ max",[HP[k]["max"] if HP[k] else None for k in ROOTS]); row("  | H⊥ final",[HP[k]["final"] if HP[k] else None for k in ROOTS]); row("  | onset(H⊥≥1e-8)",[HP[k]["onset"] if HP[k] else None for k in ROOTS])
    for key in ("R_perp_log_growth_rate_per_step","R_perp_fit_window","A_perp_final","R_perp_takagi_final","parent_residual"): row(f"  | {key}",[S[k].get(key) if S[k] else None for k in ROOTS])
    J[f"decompact_N{N}"]={"Hperp":HP,"summary":{k:{kk:S[k].get(kk) for kk in ("R_perp_log_growth_rate_per_step","A_perp_final","R_perp_takagi_final")} if S[k] else None for k in ROOTS}}
# ---- followup
P(); P("## 4. N5 followup（記事図 2〜6・8）"); P()
pkg="N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826"
for k,r in ROOTS.items():
    t=rd(r,pkg,"tol_sweep_summary.csv"); P(f"### tol 掃引 — {k}"); P("```"); P(t.to_string() if t is not None else "–"); P("```")
    fit=js(r,pkg,"tol_sweep_fit.json"); P(f"fit: {fit}"); P()
P("| 量 | 原本 | 修正版 | baseline | 等モジュラー |"); P("|---|---|---|---|---|")
FL={k:rd(r,pkg,"floquet_spectrum.csv") for k,r in ROOTS.items()}
def topm(d):
    if d is None: return None
    x=d[np.isclose(d.fd_eps,1e-7)].modulus.to_numpy(); x=np.sort(x)[::-1]; return ", ".join(f"{v:.6f}" for v in x[:4])
row("Floquet |μ| 上位4（eps=1e-7）",[topm(d) for d in FL.values()])
FS={k:js(r,pkg,"floquet_summary.json") for k,r in ROOTS.items()}; row("不安定実次元",[s["unstable_real_dimensions"] if s else None for s in FS.values()]); row("最大乗数",[s["largest_multiplier"] if s else None for s in FS.values()])
FM={k:js(r,pkg,"floquet_meta.json") for k,r in ROOTS.items()}; row("固定点欠損 |R(v)−v|",[m["relative_fixed_point_defect"] if m else None for m in FM.values()]); row("親残差",[m["parent_residual"] if m else None for m in FM.values()])
SE={k:rd(r,pkg,"N5_spectral_entropy_timeseries.csv") for k,r in ROOTS.items()}; row("S/lnM final",[float(d.entropy_over_lnM.iloc[-1]) if d is not None else None for d in SE.values()]); row("S/lnM min",[float(d.entropy_over_lnM.min()) if d is not None else None for d in SE.values()]); row("p_min final",[float(d.p_min.iloc[-1]) if d is not None else None for d in SE.values()]); row("p_max final",[float(d.p_max.iloc[-1]) if d is not None else None for d in SE.values()])
PD={k:js(r,pkg,"pump_depletion_summary.json") for k,r in ROOTS.items()}
for key in ("Hparallel_initial","Hperp_initial","Hparallel_final","Hperp_final","Htotal_min","Htotal_max","max_abs_Hsum_minus_Htotal"): row(f"pump {key}",[p.get(key) if p else None for p in PD.values()])
MD={k:rd(r,pkg,"N5_moduli_seed_sweep.csv") for k,r in ROOTS.items()}
row("moduli: counts（8 seed）",["; ".join(d.counts.astype(str)) if d is not None else None for d in MD.values()]); row("moduli: family1_mod 範囲",[f"{d.family1_mod.min():.4f}..{d.family1_mod.max():.4f}" if d is not None else None for d in MD.values()]); row("moduli: family2_mod 範囲",[f"{d.family2_mod.min():.4f}..{d.family2_mod.max():.4f}" if d is not None else None for d in MD.values()])
J["followup"]={"floquet":{k:topm(d) for k,d in FL.items()},"entropy_final":{k:(float(d.entropy_over_lnM.iloc[-1]) if d is not None else None) for k,d in SE.items()},"pump":PD}
# ---- K_sigma raw
P(); P("## 5. K_sigma（raw 枝）：N4/N5 の H⊥ 成長"); P(); P("| N | 量 | 原本 | 修正版 | baseline | 等モジュラー |"); P("|---|---|---|---|---|---|")
pkg="K_sigma_normalization_artifact_test_N4_N5_20260826"
for N in (4,5):
    D={k:rd(r,pkg,f"N{N}_raw_K_raw_observables.csv") for k,r in ROOTS.items()}
    def gg(d):
        if d is None: return None
        f=(d.H_perp/d.H_total).to_numpy(); s=d.step.to_numpy(); on=int(s[f>=1e-8][0]) if (f>=1e-8).any() else None; m=(f>1e-6)&(f<1e-2); sl=float(np.polyfit(s[m],np.log(f[m]),1)[0]) if m.sum()>3 else None
        return dict(f0=float(f[0]),onset=on,growth=sl,f_final=float(f[-1]),f_max=float(f.max()),sigma1_0=float(d.sigma_exact.iloc[0]),sigma1_final=float(d.sigma_exact.iloc[-1]))
    S={k:gg(d) for k,d in D.items()}
    for key,lab in [("f0","H⊥/H step0"),("onset","onset(≥1e-8)"),("growth","成長率 /step（1e-6<f<1e-2）"),("f_max","H⊥/H max"),("f_final","H⊥/H final"),("sigma1_0","σ₁(step0)"),("sigma1_final","σ₁(final)")]: row(f"{N} | {lab}" if key=="f0" else f"  | {lab}",[S[k][key] if S[k] else None for k in ROOTS])
    J[f"Ksigma_N{N}"]=S
# ---- gamma
P(); P("## 6. γ 掃引（刻み角 2π/n_den の連続極限）"); P()
pkg="N5_gamma_continuum_test_bundle_20260825"
for k,r in ROOTS.items():
    d=rd(r,pkg,"N5_gamma_continuum_rates.csv"); s=rd(r,pkg,"N5_gamma_continuum_stats.csv"); P(f"### {k}"); P("```"); P((d.to_string() if d is not None else "–")+"\n"+(s.to_string() if s is not None else "")); P("```")
# ---- N5 physical
P(); P("## 7. N5 physical phase step test：最終 metrics"); P(); pkg="N5_complex_simplex_complete_analysis_20260826"
for k,r in ROOTS.items():
    d=rd(r,pkg,"N5_physical_phase_metrics_5000steps.csv"); P(f"- {k}: {d.iloc[-1].to_dict() if d is not None else '–'}")
# ---- closure search / partial / nontrivial
P(); P("## 8. 閉包探索（k=2..4 exact、k=5,6 MITM）と部分閉包（最終状態が変わるので再計算）"); P(); P("| N | k | 原本 best | 修正版 best | baseline best | 等モジュラー best |"); P("|---|---|---|---|---|---|")
pkg="N14_N16_complete_nontrivial_zero_closure_search_20260826"; orig=rd(O,pkg,"N14_N16_all_cardinalities_best_results.csv")
for N in (14,15,16):
    for f in (f"rerun_exact_k234_N{N}.csv",f"rerun_mitm56_N{N}.csv"):
        a=rd(F,pkg,f); b=rd(B,pkg,f); c=rd(E,pkg,f)
        if a is None: continue
        for i,r in a.iterrows():
            k=int(r.k); o=orig[(orig.N==N)&(orig.k==k)]; ov=float(o.best_residual_found.iloc[0]) if len(o) else None
            bv=float(b[b.k==k].best_residual.iloc[0]) if b is not None and (b.k==k).any() else None
            cv=float(c[c.k==k].best_residual.iloc[0]) if c is not None and (c.k==k).any() else None
            row(f"{N} | {k}",[ov,float(r.best_residual),bv,cv])
pkg="N3_N16_partial_zero_closure_analysis_20260826"
P(); P("### 部分閉包 summary（total / best2 / best3 残差、N ごと）"); P(); P("| N | 原本 total/best2/best3 | 修正版 | baseline | 等モジュラー |"); P("|---|---|---|---|---|")
PS={k:rd(r,pkg,"N3_N16_partial_closure_summary.csv") if k=="原本" else rd(r,pkg,"out/N3_N16_partial_closure_summary.csv") for k,r in ROOTS.items()}
if all(v is not None for v in PS.values()):
    for N in range(3,17):
        vals=[]
        for d in PS.values():
            q=d[d.N==N]; vals.append(f"{q.total_closure_residual.iloc[0]:.2e} / {q.best_pair_residual.iloc[0]:.2e} / {q.best_triple_residual.iloc[0]:.2e}" if len(q) else None)
        row(str(N),vals)
open(os.path.join(H,"results","four_way_comparison.md"),"w",encoding="utf-8").write("\n".join(L)+"\n"); json.dump(J,open(os.path.join(H,"results","four_way_comparison.json"),"w"),indent=1,ensure_ascii=False,default=str); print("written", len(L), "lines")
