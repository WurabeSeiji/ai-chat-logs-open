#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Post-processing only. Never imported by the generator.
Analytic reference files are read only after all generator runs are complete.
"""
from pathlib import Path
import json, math, hashlib
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path('/mnt/data/charged8_strict_5case_20260926')
CASES=BASE/'cases'
ANA=Path('/mnt/data/charged_binary_integer_multiple_analytic_sweep')
IDS=['n1_attractive','n1_repulsive','n3_attractive','n3_repulsive','n4_attractive']
HSTEP=2*math.pi/4000.0
COL={'step':0,'P':1,'E':2,'U_re':3,'U_im':4,'H_re':5,'H_im':6,'Q':7,'N':8,'C':9,'D':10,'t':11,'x':12,'y':13,'q':14}

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def parts(cdir):
    m=json.loads((cdir/'part_manifest.json').read_text())
    return [cdir/x['file'] for x in m['parts']],m

def blocks(paths,n=250000):
    for p in paths:
        with h5py.File(p,'r') as f:
            d=f['raw_macro_trajectory']
            for i in range(0,d.shape[0],n): yield d[i:min(i+n,d.shape[0])]

def sampled(paths,total,limit=10000):
    stride=max(1,total//limit); out=[]
    for b in blocks(paths):
        st=b[:,0].astype(np.int64); m=(st%stride)==0
        if np.any(m): out.append(b[m])
    with h5py.File(paths[-1],'r') as f: last=f['raw_macro_trajectory'][-1:]
    if not out or int(out[-1][-1,0])!=int(last[0,0]): out.append(last)
    return np.concatenate(out)

def compare(cid,paths,ref,meta,gen):
    rr=ref.r.to_numpy(); tr=ref.t.to_numpy(); pr=ref.phi.to_numpy()
    srp=sxp=srt=sxt=0.0; mrp=mxp=mrt=mxt=0.0; np_=nt=0
    c0=gen['initial_state']['C']; d0=gen['initial_state']['D']; dn=dc=dd=0.0
    for b in blocks(paths):
        st=b[:,0]; r=b[:,1]; x=b[:,12]; y=b[:,13]; tt=b[:,11]; ph=st*HSTEP
        mp=ph<=pr[-1]+1e-12
        if np.any(mp):
            ri=np.interp(ph[mp],pr,rr); xr=ri*np.cos(ph[mp]); yr=ri*np.sin(ph[mp]); er=np.abs(r[mp]-ri); ex=np.hypot(x[mp]-xr,y[mp]-yr)
            mrp=max(mrp,float(er.max())); mxp=max(mxp,float(ex.max())); srp+=float(er@er); sxp+=float(ex@ex); np_+=er.size
        mt=np.isfinite(tt)&(tt<=tr[-1])
        if np.any(mt):
            ri=np.interp(tt[mt],tr,rr); pi=np.interp(tt[mt],tr,pr); xr=ri*np.cos(pi); yr=ri*np.sin(pi); er=np.abs(r[mt]-ri); ex=np.hypot(x[mt]-xr,y[mt]-yr)
            mrt=max(mrt,float(er.max())); mxt=max(mxt,float(ex.max())); srt+=float(er@er); sxt+=float(ex@ex); nt+=er.size
        dn=max(dn,float(np.max(np.abs(b[:,8]-0.25)))); dc=max(dc,float(np.max(np.abs(b[:,9]-c0)))); dd=max(dd,float(np.max(np.abs(b[:,10]-d0))))
    rs=meta['results_summary']; tg=gen['t_cross_r20']; pg=gen['phi_cross_r20']; rt=rs['t_final']; rp=rs['phi_final']
    return {
      'case_id':cid,
      'generator_crossing':{'t':tg,'phi':pg},'reference_crossing':{'t':rt,'phi':rp},
      'abs_crossing_errors':{'t':None if tg is None else abs(tg-rt),'phi':abs(pg-rp)},
      'same_phase_residuals_full_run':{'sample_count':np_,'max_abs_r':mrp,'rms_abs_r':math.sqrt(srp/np_),'max_xy_norm':mxp,'rms_xy_norm':math.sqrt(sxp/np_)},
      'same_time_residuals_finite_readout_only':{'sample_count':nt,'max_abs_r':mrt,'rms_abs_r':math.sqrt(srt/nt) if nt else None,'max_xy_norm':mxt,'rms_xy_norm':math.sqrt(sxt/nt) if nt else None},
      'identity_drifts':{'N_max':dn,'C_max':dc,'D_max':dd},
      'time_readout':{'first_nonfinite_step':gen.get('first_nonfinite_time_step'),'first_nonfinite_P':gen.get('first_nonfinite_time_P'),'crossing_time_finite':tg is not None},
      'comparison_direction':'generator -> saved analytic reference only; no feedback'
    }

def make_figures(cid,cdir,paths,ref,gen,cmp):
    fd=cdir/'figures'; fd.mkdir(exist_ok=True)
    g=sampled(paths,gen['macro_steps']+1); ir=np.linspace(0,len(ref)-1,min(10000,len(ref)),dtype=int)
    st=g[:,0]; r=g[:,1]; t=g[:,11]; x=g[:,12]; y=g[:,13]; cyc=st/4000.0
    rr=ref.r.to_numpy(); tr=ref.t.to_numpy(); pr=ref.phi.to_numpy(); xr=ref.x_rel.to_numpy(); yr=ref.y_rel.to_numpy()
    def save(fig,name):
        fig.tight_layout(); fig.savefig(fd/(name+'.svg')); fig.savefig(fd/(name+'.png'),dpi=180); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,7)); ax.plot(xr[ir],yr[ir],label='analytic reference',lw=1.1); ax.plot(x,y,'--',label='strict state generator',lw=.9); ax.set_aspect('equal',adjustable='box'); ax.set_xlabel('x / M'); ax.set_ylabel('y / M'); ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+' orbit'); save(fig,'figure01_orbit_overlay')
    fig,ax=plt.subplots(figsize=(8,5)); ax.plot(tr[ir],rr[ir],label='analytic reference'); m=np.isfinite(t); ax.plot(t[m],r[m],'--',label='strict generator (finite t)'); ax.set_xlabel('t / M'); ax.set_ylabel('r / M'); ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+' radius vs time'); save(fig,'figure02_radius_vs_time')
    fig,ax=plt.subplots(figsize=(8,5)); m=np.isfinite(t)&(t<=tr[-1]);
    if np.any(m):
        ri=np.interp(t[m],tr,rr); pi=np.interp(t[m],tr,pr); exr=np.abs(r[m]-ri); exy=np.hypot(x[m]-ri*np.cos(pi),y[m]-ri*np.sin(pi)); ax.semilogy(t[m],np.maximum(exr,1e-18),label='|Δr|'); ax.semilogy(t[m],np.maximum(exy,1e-18),label='||Δ(x,y)||')
    ax.set_xlabel('t / M'); ax.set_ylabel('absolute residual'); ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+' finite-time residuals'); save(fig,'figure03_time_residuals')
    fig,ax=plt.subplots(figsize=(8,5)); ax.plot(st,g[:,8],label='N'); ax.plot(st,g[:,9],label='C'); ax.plot(st,g[:,10],label='D'); ax.set_xlabel('macro step'); ax.set_ylabel('state value'); ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+' identity-propagated states'); save(fig,'figure04_identity_states')
    mdf=pd.read_csv(cdir/'raw_first_macro_microsteps.csv'); fig,ax=plt.subplots(figsize=(8,5));
    for nm,l in [('k1p','k1_P'),('k2p','k2_P'),('k3p','k3_P'),('k4p','k4_P'),('dP','dP')]: ax.plot(mdf.micro,mdf[nm],marker='o',label=l)
    ax.set_xlabel('microstep'); ax.set_ylabel('work-state value'); ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+' first macrostep work states'); save(fig,'figure05_microphase_work_states')
    ph=st*HSTEP; m=ph<=pr[-1]+1e-12; ri=np.interp(ph[m],pr,rr); exr=np.abs(r[m]-ri); exy=np.hypot(x[m]-ri*np.cos(ph[m]),y[m]-ri*np.sin(ph[m])); fig,ax=plt.subplots(figsize=(8,5)); ax.semilogy(cyc[m],np.maximum(exr,1e-18),label='|Δr|'); ax.semilogy(cyc[m],np.maximum(exy,1e-18),label='||Δ(x,y)||'); ax.set_xlabel('orbital cycles'); ax.set_ylabel('absolute residual'); ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+' full phase-domain residuals'); save(fig,'figure06_phase_residuals')
    (cdir/'FIGURE_INDEX_ja.md').write_text('# 図一覧\n\n1. figure01_orbit_overlay — 解析基準と生成器の軌道\n2. figure02_radius_vs_time — 半径対時間（有限時刻読み出し範囲）\n3. figure03_time_residuals — 同一時刻残差\n4. figure04_identity_states — N,C,D の保持\n5. figure05_microphase_work_states — 最初の macrostep の work state\n6. figure06_phase_residuals — 全走行の同一位相残差\n\nPNG と SVG を生成。Google Drive には指示どおり PNG を保存しない。\n',encoding='utf-8')

def analysis_md(cid,cdir,meta,gen,cmp):
    rs=meta['results_summary']; te=cmp['abs_crossing_errors']['t']; ov=gen.get('first_nonfinite_time_step')
    tline='評価不能（Q/t 読み出しが途中で非有限化）' if te is None else f'{te:.6e}'
    extra='' if ov is None else f'''\n## 時刻読み出しの数値的限界\n\n`Q` が有限精度範囲を超え、step **{ov:,}**, $P\\approx{gen['first_nonfinite_time_P']:.9g}$ で `t = TAU0 log(Q)` が非有限化した。ルールを変更せず、P/H/N/C/D の生成は同じ `transition(z)` のまま r=20 まで継続した。したがって r=20 での時刻比較だけは不能だが、位相・軌道比較は全走行で可能である。\n'''
    p=cmp['same_phase_residuals_full_run']; q=cmp['same_time_residuals_finite_readout_only']
    txt=f'''# {cid} — 厳格8状態数値実験 分析\n\n## 実験条件\n\n- $P_0=50$, target $P=20$\n- $N={gen['initial_state']['N']}$, $C={gen['initial_state']['C']}$, $D={gen['initial_state']['D']}$\n- $\\lambda_A={gen['lambda_A']}$, $\\lambda_B={gen['lambda_B']}$\n- 永続物理状態: $(U,P,E,H,Q,N,C,D)$\n- 生成写像: `transition(z)` のみ\n\n5ケースで相互作用写像・11相機械・4000 macrostep/orbit は同一であり、変更したケース依存物理量は初期状態 $C,D$ だけである。解析基準は生成完了後の比較にのみ使用した。\n\n## 実行結果\n\n- macro steps: **{gen['macro_steps']:,}**\n- microsteps: **{gen['microsteps']:,}**\n- generator $\\phi(r=20)$: `{gen['phi_cross_r20']:.15g}`\n- analytic $\\phi(r=20)$: `{rs['phi_final']:.15g}`\n- $|\\Delta\\phi|$: `{cmp['abs_crossing_errors']['phi']:.6e}`\n- $|\\Delta t|$ at r=20: **{tline}**\n- N,C,D 最大ドリフト: `{cmp['identity_drifts']}`\n\n## 全走行・同一位相比較\n\n- max $|\\Delta r|$: `{p['max_abs_r']:.6e}`\n- RMS $|\\Delta r|$: `{p['rms_abs_r']:.6e}`\n- max $||\\Delta(x,y)||$: `{p['max_xy_norm']:.6e}`\n- RMS $||\\Delta(x,y)||$: `{p['rms_xy_norm']:.6e}`\n\n## 有限時刻読み出し範囲の同一時刻比較\n\n- samples: `{q['sample_count']:,}`\n- max $|\\Delta r|$: `{q['max_abs_r']:.6e}`\n- RMS $|\\Delta r|$: `{q['rms_abs_r']:.6e}`\n- max $||\\Delta(x,y)||$: `{q['max_xy_norm']:.6e}`\n- RMS $||\\Delta(x,y)||$: `{q['rms_xy_norm']:.6e}`\n{extra}\n## 解釈上の制限\n\n解析基準と生成器は同じ leading-order 物理方程式を異なる表現で評価している。この一致は新しい物理法則の独立導出ではなく、ケース値を永続状態 N,C,D に内部化した同一の一回写像表現が保存済み解析基準を再現するかの検証である。\n'''
    (cdir/'analysis_ja.md').write_text(txt,encoding='utf-8')

def main():
    combined=[]
    for cid in IDS:
        cdir=CASES/cid; gen=json.loads((cdir/'generator_summary.json').read_text()); paths,mani=parts(cdir)
        ref=pd.read_csv(ANA/cid/'charged_binary_analytic_reference_v1_raw.csv',usecols=['r','t','phi','x_rel','y_rel','cycles']); meta=json.loads((ANA/cid/'charged_binary_analytic_reference_v1_metadata.json').read_text())
        cmp=compare(cid,paths,ref,meta,gen); (cdir/'reference_comparison.json').write_text(json.dumps(cmp,ensure_ascii=False,indent=2)+'\n'); make_figures(cid,cdir,paths,ref,gen,cmp); analysis_md(cid,cdir,meta,gen,cmp)
        combined.append({'case_id':cid,'generator':gen,'comparison':cmp,'reference_summary':meta['results_summary']})
    (BASE/'combined_results.json').write_text(json.dumps(combined,ensure_ascii=False,indent=2)+'\n')
    lines=['# 整数倍クーロン5ケース — 厳格8状態数値実験まとめ','','5ケースすべてで experiment 11 と同一の `transition(z)` / `macrostep(z)` / 8永続状態 / 11相機械を使用した。変更したケース依存物理量は初期状態 $C,D$ のみ。解析基準は生成完了後の比較にのみ使用した。','','| case | C | D | macro steps | |Δphi| at r=20 | |Δt| at r=20 | time readout |','|---|---:|---:|---:|---:|---:|---|']
    for a in combined:
        g=a['generator']; c=a['comparison']; te=c['abs_crossing_errors']['t']; ts='N/A' if te is None else f'{te:.3e}'; st='finite' if g.get('first_nonfinite_time_step') is None else f"nonfinite from step {g['first_nonfinite_time_step']:,}"; lines.append(f"| {g['case_id']} | {g['initial_state']['C']} | {g['initial_state']['D']} | {g['macro_steps']:,} | {c['abs_crossing_errors']['phi']:.3e} | {ts} | {st} |")
    lines += ['','## 主要所見','','- `n1_attractive` は既存 experiment 11 と同一ケースで交差読み出し値を再現した。','- 同符号の `n1_repulsive`, `n3_repulsive` では $D=(lambda_A-lambda_B)^2=0$ のため、現行 LO 表現の EM dipole 放射項はゼロになる。Coulomb 保存力がゼロという意味ではない。','- `n3_repulsive` は $C=0.19$ で非常に長い走行となり、Q/t 読み出しだけが途中で有限精度範囲を超えた。P/H/N/C/D の生成は同一写像で r=20 まで継続した。','- `n4_attractive` は C=2.44, D=5.76 で最も少ない macrostep で r=20 に到達した。','','## raw データ','','全 macrostep 行を HDF5 part ファイルに保存した。HDF5 化と part 分割はシリアライズだけであり、状態更新則は変更していない。`n3_repulsive` の part 境界では完全41成分状態だけを checkpoint し、次 part でそのまま再開した。']
    (BASE/'combined_analysis_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    for cid in IDS:
        cdir=CASES/cid; entries=[]
        for p in sorted(cdir.rglob('*')):
            if p.is_file() and p.suffix.lower()!='.png' and p.name!='SHA256SUMS.json': entries.append({'path':str(p.relative_to(cdir)),'sha256':sha256(p),'size':p.stat().st_size})
        (cdir/'SHA256SUMS.json').write_text(json.dumps(entries,ensure_ascii=False,indent=2)+'\n')
if __name__=='__main__': main()
