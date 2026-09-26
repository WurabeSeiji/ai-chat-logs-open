from pathlib import Path
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT=Path('/mnt/data/paper5_experiment_archive/第五思考実験_乗法状態内部展開と近似次数検証_20260924')
FIG=ROOT/'figures'; DATA=ROOT/'plot_data'
FIG.mkdir(exist_ok=True,parents=True); DATA.mkdir(exist_ok=True,parents=True)
SRC=Path('/mnt/data')

def save(fig, stem):
    fig.tight_layout()
    fig.savefig(FIG/f'{stem}.svg')
    fig.savefig(FIG/f'{stem}.png', dpi=220)
    plt.close(fig)

# ------------------------------------------------------------
# Figure 00: hidden-state audit structure
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13,7))
ax.set_xlim(0,13); ax.set_ylim(0,7); ax.axis('off')

def box(x,y,w,h,text,fc='#eeeeee',fs=10,lw=1.2):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.03,rounding_size=0.08',fc=fc,ec='black',lw=lw)
    ax.add_patch(p); ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs)
    return p

def arrow(x1,y1,x2,y2,txt=None):
    a=FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='->',mutation_scale=14,lw=1.3)
    ax.add_patch(a)
    if txt: ax.text((x1+x2)/2,(y1+y2)/2+0.16,txt,ha='center',va='center',fontsize=9)

ax.text(0.5,6.65,'Apparent Paper 4 description',fontsize=14,fontweight='bold')
box(0.6,5.0,2.2,0.9,'2-component state\nX = (a,b)',fc='#d9edf7',fs=11)
box(4.1,5.0,2.0,0.9,'state-dependent\ninteraction S_n',fc='#dff0d8',fs=11)
box(7.4,5.0,2.2,0.9,"next state\nX' = S_n X",fc='#d9edf7',fs=11)
box(10.5,5.0,1.8,0.9,'readout\nPhi(X)',fc='#f5e6cc',fs=11)
arrow(2.8,5.45,4.1,5.45); arrow(6.1,5.45,7.4,5.45); arrow(9.6,5.45,10.5,5.45)

ax.text(0.5,3.85,'Audit of actual recurrence/program state',fontsize=14,fontweight='bold')
# persistent explicit states
labels=[('chi','processing phase'),('p','semilatus rectum'),('e','eccentricity'),('phi','orbital phase'),('t','clock readout state')]
xs=[0.5,2.6,4.7,6.8,8.9]
for x,(s,desc) in zip(xs,labels):
    box(x,2.55,1.7,0.95,f'{s}\n{desc}',fc='#fde0dd',fs=9)
box(11.0,2.55,1.5,0.95,'derived\nX, r, theta',fc='#f5e6cc',fs=9)
for x in xs:
    arrow(x+0.85,2.55,6.5,1.7)
box(4.5,0.65,4.0,1.0,'closed synchronous state map\nPsi_n -> Psi_{n+1}\n(current implementation: 5 persistent states)',fc='#ccebc5',fs=11)
arrow(6.5,1.7,6.5,1.65)
arrow(8.9+1.7,3.0,11.0,3.0)

ax.text(0.55,0.15,'Audit conclusion: the two-component X is a derived representation, not the full persistent state.  Five states close the audited implementation; minimality is not proved.',fontsize=9)
ax.set_title('Hidden-state audit of Paper 4',fontsize=16)
save(fig,'figure00_hidden_state_audit_structure')

# ------------------------------------------------------------
# Figure 00b: interaction constraints and B implementation
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13,6.5))
ax.set_xlim(0,13); ax.set_ylim(0,6.5); ax.axis('off')
ax.set_title('Interaction constraints used in Paper 5 audit',fontsize=16)
box(0.5,4.9,2.5,1.0,'read all current states\nsimultaneously',fc='#d9edf7',fs=11)
box(3.7,4.9,2.5,1.0,'products and\nsums of products',fc='#dff0d8',fs=11)
box(6.9,4.9,2.5,1.0,'write all next states\nsimultaneously',fc='#d9edf7',fs=11)
box(10.1,4.9,2.4,1.0,'no updated-state\nread in same step',fc='#f2dede',fs=11)
arrow(3.0,5.4,3.7,5.4); arrow(6.2,5.4,6.9,5.4); arrow(9.4,5.4,10.1,5.4)
box(1.0,2.8,2.6,1.0,'no external step counter\nfor selecting S_n',fc='#f2dede',fs=10)
box(5.2,2.8,2.6,1.0,'internal phase state\nq -> cyclic permutation',fc='#fff2cc',fs=10)
box(9.4,2.8,2.6,1.0,'select S_n from\ninternal state only',fc='#dff0d8',fs=10)
arrow(3.6,3.3,5.2,3.3); arrow(7.8,3.3,9.4,3.3)
box(2.0,0.7,3.6,1.0,'Method A: exp/log wrapper\nworks numerically but violates\nstate-update rule',fc='#f2dede',fs=10)
box(7.4,0.7,3.6,1.0,'Method B: explicit 11-phase\ninternal state machine\npasses old-state-only audit',fc='#ccebc5',fs=10)
arrow(5.6,1.2,7.4,1.2,'adopted direction')
save(fig,'figure00b_interaction_constraints')

# ------------------------------------------------------------
# Table 01: state audit summary
# ------------------------------------------------------------
state_rows=[
    ['X=(a,b)','derived representation','No','constructed from r and phi; not sufficient to close recurrence','Paper 4 / audit'],
    ['chi','processing phase','Yes','required to evaluate cos/sin or U=e^{i chi}','5-state audit'],
    ['p','orbital element','Yes','evolves under local 2.5PN RR recurrence','5-state audit'],
    ['e','orbital element','Yes','evolves under local 2.5PN RR recurrence','5-state audit'],
    ['phi','orbital phase','Yes','accumulates angular recurrence','5-state audit'],
    ['t','clock coordinate','Yes in additive audit','accumulates dt/dchi; later represented by multiplicative Q','5-state audit'],
    ['U=e^{i chi}','multiplicative phase state','Candidate replacement','avoids log/arg inside update','multiplicative audit'],
    ['H=e^{i phi}','multiplicative phase state','Candidate replacement','phase update by multiplication','multiplicative audit'],
    ['Q=e^{t/tau0}','multiplicative clock state','Candidate replacement','t read only by log at observation','multiplicative audit'],
]
pd.DataFrame(state_rows,columns=['quantity','role','persistent_state_status','audit_reason','source_stage']).to_csv(DATA/'table01_state_audit_summary.csv',index=False)

# Table 02: approximation-order summary
order_rows=[
    ['physical dissipative model','leading 2.5PN radiation reaction','higher PN terms not included','physical model order'],
    ['angular recurrence','N=12 baseline; tested N=12,16,20,24','first omitted term O(u^(N+1))','series truncation order'],
    ['Paper 4 time integration','RK4','local O(h^5), global O(h^4)','numerical integration order'],
    ['B 11-phase expansion','exact refactorization of RK4 processing order','no new truncation error','internal state-machine depth'],
    ['Taylor order study','4, 6, 8','empirical convergence approximately 4, 6, 8','numerical approximation depth'],
]
pd.DataFrame(order_rows,columns=['layer','current_or_tested_order','error_or_omitted_term','meaning']).to_csv(DATA/'table02_approximation_orders.csv',index=False)

# Table 03: key numerical audit results
with open(SRC/'five_state_audit_results.json') as f: five=json.load(f)
with open(SRC/'method_B_strict_verification.json') as f: strict=json.load(f)
with open(SRC/'two_path_rk4_fast_audit.json') as f: ab=json.load(f)
ang=pd.read_csv(SRC/'angular_order_convergence.csv')
integ9=pd.read_csv(SRC/'integration_order_9case.csv')
# strongest angular case p=24.3,e=0.55
ang_s=ang[(ang.p0==24.3)&(ang.e0==0.55)].sort_values('N')
# average empirical order, use finest available transition where finite by designed order
emp=[]
for o,g in integ9.groupby('order'):
    vals=g['empirical_order'].dropna()
    # prefer steps >=16, mirrors asymptotic region
    gg=g[g['steps']>=16]['empirical_order'].dropna()
    emp.append((int(o),float(gg.mean() if len(gg) else vals.mean())))
rows=[]
rows.append(['five-state reconstruction','max |Delta t| over audit cases',five['worst']['max_abs_t'],'same Paper4 trajectory reproduced with explicit 5-state map'])
rows.append(['five-state reconstruction','max |Delta xy| over audit cases',five['worst']['max_abs_xy'],'readout-space agreement'])
rows.append(['strict Method B','global max |Delta H| over 12 cases',strict['all12_global_max_abs'][3],'11-phase state machine vs reference'])
rows.append(['strict Method B','global max |Delta Q| over 12 cases',strict['all12_global_max_abs'][4],'11-phase state machine vs reference'])
for k,v in strict['saved_C1_csv_max_abs'].items():
    rows.append(['strict Method B',f'C1 max |Delta {k}| vs saved CSV',v,'saved Paper4 C1 reference'])
for r in ang_s.itertuples():
    rows.append(['angular order',f'p=24.3,e=0.55, N={int(r.N)} phi error',r.phi_abs_error,'3-orbit comparison to exact angular factor'])
for o,val in emp:
    rows.append(['integration order',f'designed order {o}: mean empirical order (steps>=16)',val,'9-case convergence study'])
pd.DataFrame(rows,columns=['experiment','metric','value','interpretation']).to_csv(DATA/'table03_key_numerical_results.csv',index=False)

# Write plan markdown
plan='''# Paper 5 figure/table plan\n\nNo physics simulation was rerun for these figures. All plots and tables are derived from the existing audit JSON/CSV files.\n\n## Recommended main-text figures\n\n1. `figure00_hidden_state_audit_structure` — central motivation: apparent two-state model versus audited five-state persistent state.\n2. `figure00b_interaction_constraints` — explicit computational rules and why Method B was selected.\n3. `figure05_methodB_strict_heatmap` — 12-case strict verification of Method B.\n4. `figure06_angular_order_convergence` — angular truncation convergence N=12,16,20,24.\n5. `figure07_integration_order_convergence` — numerical order 4,6,8 convergence.\n6. `figure08_integration_order_9case_summary` — transfer of order improvement across nine cases.\n\n## Recommended supplementary figures\n\n- `figure02_five_state_reconstruction_errors`\n- `figure03_full_inline_equivalence`\n- `figure04a_methodA_vs_B_internal`\n- `figure04b_methodA_vs_B_observables`\n- `figure05b_methodB_savedC1`\n- `figure01_experiment_flow` (archive/history figure, not essential for main paper)\n\n## Recommended main-text tables\n\n- `table01_state_audit_summary.csv` — what is state, derived quantity, or candidate multiplicative replacement.\n- `table02_approximation_orders.csv` — separate PN order, angular truncation order, integration order, and internal microstep depth.\n- `table03_key_numerical_results.csv` — compact numerical evidence used in the conclusions.\n\n## Storage policy\n\nBoth SVG and PNG are generated locally. PNG files are not intended for Google Drive upload. SVG, plotting programs, source CSV tables, and this index are the archival versions.\n'''
(ROOT/'FIGURE_TABLE_PLAN_ja.md').write_text(plan,encoding='utf-8')
print('done')
