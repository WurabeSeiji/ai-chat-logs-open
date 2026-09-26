from pathlib import Path
import json
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/mnt/data/paper5_experiment_archive/第五思考実験_乗法状態内部展開と近似次数検証_20260924')
FIG = ROOT / 'figures'
DATA = ROOT / 'plot_data'
FIG.mkdir(exist_ok=True, parents=True)
DATA.mkdir(exist_ok=True, parents=True)
SRC = Path('/mnt/data')


def savefig(fig, stem):
    fig.tight_layout()
    fig.savefig(FIG / f'{stem}.svg')
    fig.savefig(FIG / f'{stem}.png', dpi=200)
    plt.close(fig)

# ------------------------------------------------------------
# Figure 01: experiment flow / coverage map
# ------------------------------------------------------------
exp_idx = pd.read_csv(ROOT / 'EXPERIMENT_INDEX.csv')
exp_idx.to_csv(DATA / 'figure01_experiment_index_source.csv', index=False)
alias_map = {
    '00_reference_from_paper4': 'baseline from Paper 4',
    '01_五状態同期写像の再構成': 'five-state synchronous reconstruction',
    '02_乗法状態表現と同期インライン化の経緯': 'multiplicative encoding history',
    '03_A局所exp-log_vs_B内部状態展開': 'method A vs method B',
    '04_B方式_11相同期状態機械_厳密監査': 'strict 11-phase state machine',
    '05_演算ブロック完全展開_Laurent_N12': 'full operator expansion (Laurent N=12)',
    '06_近似次数と状態_S不増加検証': 'order increase without more states',
}
fig, ax = plt.subplots(figsize=(13, 5.5))
ax.axis('off')
ys = np.arange(len(exp_idx))[::-1]
status_color = {
    'reference': '#d9d9d9',
    'numerically_verified': '#9ecae1',
    'mixed_historical': '#fdd0a2',
    'numerically_compared': '#c7e9c0',
    'numerically_verified_algebra': '#c6dbef',
}
for y, (_, row) in zip(ys, exp_idx.iterrows()):
    c = status_color.get(row['status'], '#eeeeee')
    ax.add_patch(plt.Rectangle((0.02, y-0.32), 0.96, 0.62, color=c, ec='black', lw=0.8))
    alias = alias_map.get(row['folder'], row['folder'])
    txt = f"Exp {row['experiment']}: {alias}\nPurpose: {row['purpose']}\nKey output: {row['key_output']}"
    ax.text(0.04, y, txt, va='center', ha='left', fontsize=9)
ax.set_xlim(0, 1)
ax.set_ylim(-0.8, len(exp_idx)-0.2)
ax.set_title('Paper 5 audit flow and experiment coverage', fontsize=14)
fig.text(0.01, 0.01, 'Color indicates status: reference / verified / compared / historical / algebra-verified', fontsize=8)
savefig(fig, 'figure01_experiment_flow')

# ------------------------------------------------------------
# Figure 02: five-state reconstruction errors
# ------------------------------------------------------------
with open(SRC / 'five_state_audit_results.json') as f:
    five = json.load(f)
worst = five['worst']
rows = five['cases']
worst_df = pd.DataFrame([worst])
worst_df.to_csv(DATA / 'figure02_five_state_reconstruction_worst.csv', index=False)
qtys = ['chi','p','e','phi','t','r','X','Y']
vals = [worst[f'max_abs_{q}'] for q in qtys]
fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(qtys, vals)
ax.set_yscale('log')
ax.set_ylabel('max abs error')
ax.set_title('Five-state synchronous reconstruction: worst-case errors')
for i,v in enumerate(vals):
    ax.text(i, v*1.2, f'{v:.1e}', ha='center', va='bottom', fontsize=8, rotation=90)
savefig(fig, 'figure02_five_state_reconstruction_errors')

# ------------------------------------------------------------
# Figure 03: full inline equivalence errors
# ------------------------------------------------------------
with open(SRC / 'five_state_full_inline_audit.json') as f:
    full = json.load(f)
state_order = full['state_order']
seq_inline = full['global_max_sequential_vs_inline']
inline_round = full['global_max_inline_vs_exp_product_roundtrip']
pd.DataFrame({
    'state': state_order,
    'seq_vs_inline': seq_inline,
    'inline_vs_exp_roundtrip': inline_round,
}).to_csv(DATA / 'figure03_full_inline_equivalence.csv', index=False)
fig, ax = plt.subplots(figsize=(10, 5))
x = np.arange(len(state_order))
w = 0.38
ax.bar(x-w/2, seq_inline, width=w, label='sequential vs inline')
ax.bar(x+w/2, inline_round, width=w, label='inline vs exp-product roundtrip')
ax.set_xticks(x)
ax.set_xticklabels(state_order)
ax.set_yscale('log')
ax.set_ylabel('max abs error')
ax.set_title('Inline expansion audit in additive coordinates')
ax.legend()
savefig(fig, 'figure03_full_inline_equivalence')

# ------------------------------------------------------------
# Figure 04: method A vs method B comparison
# ------------------------------------------------------------
with open(SRC / 'two_path_rk4_fast_audit.json') as f:
    ab = json.load(f)
# internal state comparison
internal_labels = ['U','P','E','H','Q']
gA = ab['global_A']
gB = ab['global_B']
pd.DataFrame({'state': internal_labels, 'method_A': gA, 'method_B': gB}).to_csv(DATA / 'figure04a_methodA_vs_B_internal.csv', index=False)
fig, ax = plt.subplots(figsize=(10,5))
x = np.arange(len(internal_labels)); w=0.38
ax.bar(x-w/2, gA, width=w, label='Method A')
ax.bar(x+w/2, gB, width=w, label='Method B')
ax.set_xticks(x)
ax.set_xticklabels(internal_labels)
ax.set_yscale('log')
ax.set_ylabel('global max abs error')
ax.set_title('Method A vs Method B: internal-state agreement')
ax.legend()
savefig(fig, 'figure04a_methodA_vs_B_internal')
# observable comparison
csvd = ab['csv']
obs_labels = ['p','e','t','xy']
A_obs = [csvd['A_p'], csvd['A_e'], csvd['A_t'], csvd['A_xy']]
B_obs = [csvd['B_p'], csvd['B_e'], csvd['B_t'], csvd['B_xy']]
pd.DataFrame({'observable': obs_labels, 'method_A': A_obs, 'method_B': B_obs}).to_csv(DATA / 'figure04b_methodA_vs_B_observables.csv', index=False)
fig, ax = plt.subplots(figsize=(10,5))
x = np.arange(len(obs_labels)); w=0.38
ax.bar(x-w/2, A_obs, width=w, label='Method A')
ax.bar(x+w/2, B_obs, width=w, label='Method B')
ax.set_xticks(x)
ax.set_xticklabels(obs_labels)
ax.set_yscale('log')
ax.set_ylabel('max abs error vs saved CSV')
ax.set_title('Method A vs Method B: observable agreement to saved reference')
ax.legend()
savefig(fig, 'figure04b_methodA_vs_B_observables')

# ------------------------------------------------------------
# Figure 05: strict method B verification over 12 cases
# ------------------------------------------------------------
with open(SRC / 'method_B_strict_verification.json') as f:
    strict = json.load(f)
cases_df = pd.DataFrame(strict['cases'])
# expand max_abs
maxabs = pd.DataFrame(cases_df['max_abs'].tolist(), columns=['U','P','E','H','Q'])
plot_df = pd.concat([cases_df[['p0','e0','q']], maxabs], axis=1)
plot_df.to_csv(DATA / 'figure05_methodB_strict_cases.csv', index=False)
# heatmap in log10 floor units
vals = maxabs.to_numpy(dtype=float)
vals_clipped = np.log10(np.maximum(vals, 1e-18))
fig, ax = plt.subplots(figsize=(9,6))
im = ax.imshow(vals_clipped, aspect='auto')
ax.set_xticks(range(5))
ax.set_xticklabels(['U','P','E','H','Q'])
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels([f"p={r.p0:g}, e={r.e0:g}, q={r.q:g}" for r in plot_df.itertuples()])
ax.set_title('Strict Method B verification across 12 cases (log10 max abs error)')
cb = fig.colorbar(im, ax=ax)
cb.set_label('log10(error)')
savefig(fig, 'figure05_methodB_strict_heatmap')

# strict summary with C1 observable agreement
saved = strict['saved_C1_csv_max_abs']
globalv = strict['all12_global_max_abs']
pd.DataFrame({'state':['U','P','E','H','Q'], 'global_max_abs': globalv}).to_csv(DATA / 'figure05b_methodB_strict_global.csv', index=False)
pd.DataFrame({'observable': list(saved.keys()), 'max_abs': list(saved.values())}).to_csv(DATA / 'figure05c_methodB_strict_savedC1.csv', index=False)
fig, ax = plt.subplots(figsize=(10,5))
labels = ['P','E','t','x','y']
vals = [saved['p'], saved['e'], saved['t'], saved['x'], saved['y']]
ax.bar(labels, vals)
ax.set_yscale('log')
ax.set_ylabel('max abs error')
ax.set_title('Strict Method B: agreement with saved C1 reference CSV')
savefig(fig, 'figure05b_methodB_savedC1')

# ------------------------------------------------------------
# Figure 06: angular truncation convergence
# ------------------------------------------------------------
ang = pd.read_csv(SRC / 'angular_order_convergence.csv')
ang.to_csv(DATA / 'figure06_angular_order_convergence_source.csv', index=False)
fig, ax = plt.subplots(figsize=(10,6))
for (p0, e0), grp in ang.groupby(['p0','e0']):
    grp=grp.sort_values('N')
    ax.plot(grp['N'], grp['phi_abs_error'], marker='o', label=f'p={p0:g}, e={e0:g}')
ax.set_yscale('log')
ax.set_xlabel('angular truncation order N')
ax.set_ylabel('abs error of phi')
ax.set_title('Angular truncation convergence (existing data only)')
ax.legend(fontsize=7, ncol=3)
savefig(fig, 'figure06_angular_order_convergence')

# ------------------------------------------------------------
# Figure 07: integration order convergence (representative case)
# ------------------------------------------------------------
integ = pd.read_csv(SRC / 'integration_order_convergence.csv')
integ.to_csv(DATA / 'figure07_integration_order_convergence_source.csv', index=False)
fig, ax = plt.subplots(figsize=(10,6))
for order, grp in integ.groupby('order'):
    grp = grp.sort_values('steps_per_orbit')
    ax.plot(grp['steps_per_orbit'], grp['combined'], marker='o', label=f'order {order}')
ax.set_xscale('log', base=2)
ax.set_yscale('log')
ax.set_xlabel('steps per orbit')
ax.set_ylabel('combined error metric')
ax.set_title('Integration-order convergence for representative case')
ax.legend()
savefig(fig, 'figure07_integration_order_convergence')

# ------------------------------------------------------------
# Figure 08: 9-case integration summary at finest resolution
# ------------------------------------------------------------
integ9 = pd.read_csv(SRC / 'integration_order_9case.csv')
integ9.to_csv(DATA / 'figure08_integration_order_9case_source.csv', index=False)
finest = integ9.sort_values('steps').groupby(['p0','e0','order']).tail(1).copy()
finest['case'] = finest.apply(lambda r: f"p={r['p0']:g}, e={r['e0']:g}", axis=1)
pivot = finest.pivot(index='case', columns='order', values='combined').sort_index()
pivot.to_csv(DATA / 'figure08_integration_order_9case_finest.csv')
fig, ax = plt.subplots(figsize=(10,6))
pivot.plot(kind='bar', ax=ax)
ax.set_yscale('log')
ax.set_ylabel('combined error at finest steps')
ax.set_xlabel('case')
ax.set_title('Nine-case summary: higher integration order without more states')
ax.legend(title='order')
savefig(fig, 'figure08_integration_order_9case_summary')

# ------------------------------------------------------------
# Figure index markdown
# ------------------------------------------------------------
index_lines = [
    '# Figure index for Paper 5 supplementary archive',
    '',
    '| Figure | File stem | Role | Data source |',
    '|---|---|---|---|',
    '| 1 | figure01_experiment_flow | audit path and experiment coverage | EXPERIMENT_INDEX.csv |',
    '| 2 | figure02_five_state_reconstruction_errors | hidden-state audit reconstructed by 5-state model | five_state_audit_results.json |',
    '| 3 | figure03_full_inline_equivalence | fully inlined additive-coordinate audit | five_state_full_inline_audit.json |',
    '| 4a | figure04a_methodA_vs_B_internal | Method A vs Method B internal-state comparison | two_path_rk4_fast_audit.json |',
    '| 4b | figure04b_methodA_vs_B_observables | Method A vs Method B observable comparison | two_path_rk4_fast_audit.json |',
    '| 5a | figure05_methodB_strict_heatmap | strict 12-case verification heatmap | method_B_strict_verification.json |',
    '| 5b | figure05b_methodB_savedC1 | strict B agreement with saved C1 reference | method_B_strict_verification.json |',
    '| 6 | figure06_angular_order_convergence | angular truncation order study | angular_order_convergence.csv |',
    '| 7 | figure07_integration_order_convergence | integration order study | integration_order_convergence.csv |',
    '| 8 | figure08_integration_order_9case_summary | 9-case finest-resolution summary | integration_order_9case.csv |',
    '',
    'SVG and PNG are both generated from existing stored data only; no physical simulation is rerun in this script.',
]
(ROOT / 'FIGURE_INDEX_ja.md').write_text('\n'.join(index_lines), encoding='utf-8')
print('Generated figures in', FIG)
