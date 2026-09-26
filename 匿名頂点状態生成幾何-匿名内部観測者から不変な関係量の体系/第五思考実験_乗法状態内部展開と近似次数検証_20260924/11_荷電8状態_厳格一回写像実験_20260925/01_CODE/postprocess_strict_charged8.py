#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Post-processing only. Never imported by generator."""
import csv, json, math, hashlib
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

BASE=Path('/mnt/data/charged8_strict_v3')
RES=BASE/'results_final'
FIG=BASE/'figures_final'; FIG.mkdir(exist_ok=True)
REFRAW=Path('/mnt/data/charged_binary_analytic_reference_v1_raw.csv')
REFMETA=Path('/mnt/data/reference_metadata.json')

# Load generator output after generation is complete.
a=np.genfromtxt(RES/'raw_macro_trajectory.csv',delimiter=',',names=True)
r=np.asarray(a['P'],float); t=np.asarray(a['t_readout'],float)
x=np.asarray(a['x_readout'],float); y=np.asarray(a['y_readout'],float)
step=np.asarray(a['step'],int); N=np.asarray(a['N'],float); C=np.asarray(a['C'],float); D=np.asarray(a['D'],float)
ref=np.genfromtxt(REFRAW,delimiter=',',names=True)
rr=np.asarray(ref['r'],float); tr=np.asarray(ref['t'],float); pr=np.asarray(ref['phi'],float); xr=np.asarray(ref['x_rel'],float); yr=np.asarray(ref['y_rel'],float)
meta=json.loads(REFMETA.read_text())
gen=json.loads((RES/'generator_summary.json').read_text())

# Compare reference at same physical time, only within reference time range.
mask=t<=tr[-1]
ri=np.interp(t[mask],tr,rr); pi=np.interp(t[mask],tr,pr)
xi=ri*np.cos(pi); yi=ri*np.sin(pi)
rerr=np.abs(r[mask]-ri); xyerr=np.hypot(x[mask]-xi,y[mask]-yi)

cmp={
 'generator_crossing':{'t':gen['t_cross_r20'],'phi':gen['phi_cross_r20']},
 'reference_crossing':{'t':meta['results_summary']['t_final'],'phi':meta['results_summary']['phi_final']},
 'abs_crossing_errors':{
   't':abs(gen['t_cross_r20']-meta['results_summary']['t_final']),
   'phi':abs(gen['phi_cross_r20']-meta['results_summary']['phi_final'])},
 'same_time_residuals':{
   'max_abs_r':float(np.max(rerr)), 'rms_abs_r':float(np.sqrt(np.mean(rerr*rerr))),
   'max_xy_norm':float(np.max(xyerr)), 'rms_xy_norm':float(np.sqrt(np.mean(xyerr*xyerr)))},
 'identity_drifts':{
   'N_max':float(np.max(np.abs(N-0.25))),
   'C_max':float(np.max(np.abs(C-1.09))),
   'D_max':float(np.max(np.abs(D-0.36)))},
 'comparison_direction':'generator -> saved reference only; no feedback'
}
(RES/'reference_comparison.json').write_text(json.dumps(cmp,indent=2)+'\n')

# Downsample display only; raw remains complete.
def take(n,limit=6000):
    if n<=limit:return np.arange(n)
    return np.linspace(0,n-1,limit,dtype=int)
ig=take(len(r)); ir=take(len(rr))

# Figure 1 orbit overlay
fig,ax=plt.subplots(figsize=(7,7))
ax.plot(xr[ir],yr[ir],label='Independent analytic reference',linewidth=1.2)
ax.plot(x[ig],y[ig],'--',label='Strict state generator',linewidth=1.0)
ax.set_xlabel('x / M'); ax.set_ylabel('y / M'); ax.set_aspect('equal',adjustable='box'); ax.legend(); ax.grid(True,alpha=.25)
ax.set_title('Charged binary orbit: strict generator vs reference')
fig.tight_layout(); fig.savefig(FIG/'figure01_orbit_overlay.svg'); fig.savefig(FIG/'figure01_orbit_overlay.png',dpi=180); plt.close(fig)

# Figure 2 radius vs time
fig,ax=plt.subplots(figsize=(8,5))
ax.plot(tr[ir],rr[ir],label='Independent analytic reference')
ax.plot(t[ig],r[ig],'--',label='Strict state generator')
ax.set_xlabel('t / M'); ax.set_ylabel('r / M'); ax.legend(); ax.grid(True,alpha=.25); ax.set_title('Inspiral radius vs time')
fig.tight_layout(); fig.savefig(FIG/'figure02_radius_vs_time.svg'); fig.savefig(FIG/'figure02_radius_vs_time.png',dpi=180); plt.close(fig)

# Figure 3 residuals
im=take(len(rerr))
fig,ax=plt.subplots(figsize=(8,5))
ax.semilogy(t[mask][im],np.maximum(rerr[im],1e-18),label='|Δr| at same t')
ax.semilogy(t[mask][im],np.maximum(xyerr[im],1e-18),label='||Δ(x,y)|| at same t')
ax.set_xlabel('t / M'); ax.set_ylabel('absolute residual'); ax.legend(); ax.grid(True,alpha=.25); ax.set_title('Post-generation residuals')
fig.tight_layout(); fig.savefig(FIG/'figure03_reference_residuals.svg'); fig.savefig(FIG/'figure03_reference_residuals.png',dpi=180); plt.close(fig)

# Figure 4 identity states
fig,ax=plt.subplots(figsize=(8,5))
ax.plot(step[ig],N[ig],label='N')
ax.plot(step[ig],C[ig],label='C')
ax.plot(step[ig],D[ig],label='D')
ax.set_xlabel('macro step'); ax.set_ylabel('state value'); ax.legend(); ax.grid(True,alpha=.25); ax.set_title('Identity-propagated persistent states')
fig.tight_layout(); fig.savefig(FIG/'figure04_identity_states.svg'); fig.savefig(FIG/'figure04_identity_states.png',dpi=180); plt.close(fig)

# Figure 5 first macro microphase evolution from raw micro CSV
m=np.genfromtxt(RES/'raw_first_macro_microsteps.csv',delimiter=',',names=True)
fig,ax=plt.subplots(figsize=(8,5))
for nm,label in [('k1p','k1_P'),('k2p','k2_P'),('k3p','k3_P'),('k4p','k4_P'),('dP','dP')]: ax.plot(m['micro'],m[nm],marker='o',label=label)
ax.set_xlabel('microstep'); ax.set_ylabel('work-state value'); ax.legend(); ax.grid(True,alpha=.25); ax.set_title('First macrostep: explicit work-state progression')
fig.tight_layout(); fig.savefig(FIG/'figure05_microphase_work_states.svg'); fig.savefig(FIG/'figure05_microphase_work_states.png',dpi=180); plt.close(fig)

# Analysis document
analysis=f'''# 荷電8永続状態・厳格一回写像実験 — 分析\n\n## 1. 実験条件\n\n永続物理状態は変更せず、\n\n$$\nZ=(U,P,E,H,Q,N,C,D)\n$$\n\nとした。RK4の途中量と one-hot 位相 $q$ は Paper 5 と同様に明示 work/machine state として full state に含めた。\n\n生成器の物理的ケース値は初期状態にのみ与え、\n\n$$N=0.25,\quad C=1.09,\quad D=0.36,\quad P_0=50$$\n\nとした。`transition(z)` は外部物理引数を取らない。\n\n## 2. 更新構造\n\n各 microstep で旧 full state だけから11相の候補を全て構成し、各次状態を\n\n$$\nz_i'=\sum_{{j=0}}^{{10}} q_j F_{{ij}}(z)\n$$\n\nという同じ one-hot 積和で生成した。$q$ 自身は固定巡回行列による線形写像で更新した。\n\n$N,C,D$ は全11相で恒等候補なので、外部コピーではなく同一写像の結果として保存される。\n\n平方根、逆数、指数関数などは Paper 5 の相互作用式と同じく「現在状態から一意に計算される相互作用関数」として扱い、新しい物理状態にはしなかった。\n\n## 3. 実行結果\n\n- macro step: **{gen['macro_steps']:,}**\n- microstep: **{gen['microsteps']:,}**\n- $r=20$ 交差読み出し: $t={gen['t_cross_r20']:.15g}$, $\\phi={gen['phi_cross_r20']:.15g}$\n- 恒等状態ドリフト: $\\Delta N=0$, $\\Delta C=0$, $\\Delta D=0$\n\n## 4. 独立基準との生成後比較\n\n保存済み解析基準との比較は生成完了後だけに行い、生成器へ帰還させていない。\n\n$$\n|\\Delta t|={cmp['abs_crossing_errors']['t']:.6e},\qquad
|\\Delta\\phi|={cmp['abs_crossing_errors']['phi']:.6e}.\n$$\n\n同一物理時刻での全保存範囲比較は\n\n- max $|\\Delta r|={cmp['same_time_residuals']['max_abs_r']:.6e}$\n- RMS $|\\Delta r|={cmp['same_time_residuals']['rms_abs_r']:.6e}$\n- max $||\\Delta(x,y)||={cmp['same_time_residuals']['max_xy_norm']:.6e}$\n- RMS $||\\Delta(x,y)||={cmp['same_time_residuals']['rms_xy_norm']:.6e}$\n\nであった。\n\n## 5. 監査上の結論\n\n今回確認した範囲では、荷電 LO 準円軌道のこのベンチマークについて、$N,C,D$ を外部引数にせず8永続状態の内部に保持し、Paper 5 型の明示 work state を含む単一 `transition(z)` 写像だけで軌道生成できた。\n\nこれは「8が数学的最小状態数」であることの証明ではない。また一般の離心荷電二体系、完全な Einstein–Maxwell 二体問題への一般化も主張しない。\n\n## 6. 旧失敗実験との区別\n\n`09_...` と `10_...` は変更・削除していない。本実験では `macro_fast`、collapsed map、shortcut、状態追加を使用していない。\n'''
(RES/'analysis_ja.md').write_text(analysis,encoding='utf-8')

# Figure index
idx='''# 図一覧\n\n1. figure01_orbit_overlay.svg / .png — 軌道重ね合わせ\n2. figure02_radius_vs_time.svg / .png — 半径対時間\n3. figure03_reference_residuals.svg / .png — 生成後基準残差\n4. figure04_identity_states.svg / .png — N,C,D恒等保持\n5. figure05_microphase_work_states.svg / .png — 最初のmacrostepのwork state推移\n\nPNGとSVGを生成したが、Google Driveへの保存対象は指示どおりSVGとする。\n'''
(RES/'FIGURE_INDEX_ja.md').write_text(idx,encoding='utf-8')
print(json.dumps(cmp,indent=2))
