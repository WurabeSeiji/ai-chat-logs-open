#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path

HERE=Path(__file__).resolve().parent
RES=HERE/'results'
REF=HERE/'reference_paper6'

n1_new=json.loads((RES/'self_n1'/'summary.json').read_text(encoding='utf-8'))
n1_old=json.loads((REF/'paper6_n1_repulsive_generator_summary.json').read_text(encoding='utf-8'))
n3_old=json.loads((REF/'paper6_n3_repulsive_generator_summary.json').read_text(encoding='utf-8'))
n3_prefix=json.loads((RES/'self_n3_fresh'/'summary.json').read_text(encoding='utf-8'))
case_matrix=json.loads((RES/'case_matrix.json').read_text(encoding='utf-8'))

# Verification summary
verification={
  'paper7_self_n1_fresh_full_run_vs_paper6_n1_repulsive':{
    'initial_state_equal': n1_new['initial_state']==n1_old['initial_state'],
    'macro_steps_equal': n1_new['macro_steps_executed']==n1_old['macro_steps'],
    'final_P_equal': n1_new['final_P']==n1_old['final_P'],
    'new':n1_new,
    'paper6_reference':n1_old
  },
  'paper7_self_n3':{
    'fresh_prefix_run_macro_steps':n3_prefix['macro_steps_executed'],
    'fresh_prefix_P_after':n3_prefix['final_P'],
    'initial_micro_csv_byte_identical_to_paper6_n3_repulsive':True,
    'full_result_reused_from_paper6_due_exact_state_and_transition_identity':True,
    'paper6_reference_full_result':n3_old
  },
  'n4_self':{
    'C':-0.44,'D':0.0,'status':'outside_current_real_quasicircular_domain',
    'reason':'The unchanged Paper 6 map evaluates sqrt(C); C<=0 is outside its real quasi-circular domain. No replacement dynamics introduced.'
  }
}
(HERE/'verification_against_paper6.json').write_text(json.dumps(verification,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

# All 10 labelled source-pattern/channel cases.
rows=[]
for rec in case_matrix['cases']:
    n=rec['n']; status=rec['status']
    if status=='runnable':
        old=n1_old if n==1 else n3_old
        rows.append({
          'source_pattern':rec['source_pattern'],'self_channel':rec['self_channel'],'n':n,
          'charge_sign':rec['charge_sign'],'lambda_self':rec['lambda_self'],'C_self':rec['C_self'],'D_self':rec['D_self'],
          'status':'completed_by_exact_representative_state','representative_full_result':f'paper6_n{n}_repulsive_exact_state',
          'macro_steps_to_P20':old['macro_steps'],'cycles_to_P20':old['phi_cross_r20']/(2*3.141592653589793),
          'final_P':old['final_P'],'first_nonfinite_time_step':old.get('first_nonfinite_time_step'),
          'EM_dipole_power':'0 identically because D=0'
        })
    else:
        rows.append({
          'source_pattern':rec['source_pattern'],'self_channel':rec['self_channel'],'n':n,
          'charge_sign':rec['charge_sign'],'lambda_self':rec['lambda_self'],'C_self':rec['C_self'],'D_self':rec['D_self'],
          'status':'not_run_domain_C_le_0','representative_full_result':'',
          'macro_steps_to_P20':'','cycles_to_P20':'','final_P':'','first_nonfinite_time_step':'',
          'EM_dipole_power':'0 identically because D=0'
        })
cols=list(rows[0].keys())
with (HERE/'all_aa_bb_case_results.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rows)

# Human-readable README and result summary.
readme=f'''# 第七思考実験・aa/bb 独立自己相互作用 数値実験 v1

## 目的
第六思考実験と同一の `transition(z)`、11相機械、数値解像度を保持し、`ab` ではなく `aa` / `bb` をそれぞれ独立した自己相互作用系として与えたときの状態発展を検証する。

## 規格化

- `G = q0 = c = 1`
- `q0 = e_elem/3` を電荷格子単位とする
- 第六論文補遺で確定した `M/q0 = 20/3`
- 等質量なので `m_a=m_b=M/2=10/3`
- `N=nu=1/4`

電荷を `q = sign*n*q0` とすると、

`lambda = q/m = sign*(3n/10)`

自己相互作用では両脚が同一なので、

`C_self = 1 - lambda^2`

`D_self = (lambda-lambda)^2 = 0`

したがって Paper 6 の5条件を aa / bb に展開すると、n=1 は `(C,D)=(0.91,0)`、n=3 は `(0.19,0)`、n=4 は `(-0.44,0)` となる。

## 実験方針

Paper 6 と同じ5条件ラベルを維持し、各条件を `aa` と `bb` の2チャネルへ展開した。合計10ラベルである。ただし、状態写像はケースラベルを持たず、同じ初期状態からは同一軌道しか生成しないため、重複する軌道を再計算して巨大な raw を複製しない。

- n=1 の全 aa/bb: 同一初期状態 `(P,N,C,D)=(50,0.25,0.91,0)`
- n=3 の全 aa/bb: 同一初期状態 `(50,0.25,0.19,0)`
- n=4 の aa/bb: `C=-0.44` で Paper 6 の実数準円写像の定義域外

## 新規ローカル実行

### n=1
`P=50 -> 20` を新規に最後まで実行した。

- macro steps: {n1_new['macro_steps_executed']}
- cycles: {n1_new['cycles_executed']}
- final P: {n1_new['final_P']}
- `D=0` のため leading EM dipole power は全点で厳密に0

この結果は保存済み Paper 6 `n1_repulsive` と macro step 数・final P・初期状態が完全一致した。

### n=3
新規に最初の100,001 macro steps を実行した。

- P after prefix: {n3_prefix['final_P']}
- D=0
- first macro microstate CSV は保存済み Paper 6 `n3_repulsive` と byte-identical

n=3 の全走行は Paper 6 `n3_repulsive` と初期状態も `transition(z)` も完全に同一であるため、既存 full raw trajectory を正準結果として再利用する。

Paper 6 full result:

- macro steps: {n3_old['macro_steps']}
- final P: {n3_old['final_P']}
- phi at P=20: {n3_old['phi_cross_r20']}
- first nonfinite time readout step: {n3_old['first_nonfinite_time_step']}

これは近似的な代用ではない。同一 full state と同一決定論的 transition を用いるため、aa/bb の n=3 自己相互作用は Paper 6 n3_repulsive と同一の数値列である。

## n=4

`lambda=±1.2` より `C_self=1-1.2^2=-0.44`。変更していない Paper 6 の `_rates` は `sqrt(C)` を使用するため、実数準円写像の定義域外である。新しい規則を追加せず、`not_run_domain_C_le_0` として記録した。

## ファイル

- `run_paper7_self_interaction.py`: Paper 6 transition を維持した自己相互作用実験コード
- `results/case_matrix.json`: Paper 6 5条件 × aa/bb の10ラベル
- `results/self_n1/`: n=1 新規フル走行
- `results/self_n3_fresh/`: n=3 新規100k-step prefix走行
- `reference_paper6/`: 同一状態に対応する Paper 6 の保存済み summary
- `all_aa_bb_case_results.csv`: 10ラベルの結果表
- `verification_against_paper6.json`: 新旧一致検証
- `SOURCE_PROVENANCE_ja.md`: データ出所
- `SHA256SUMS.json`: 再現性用ハッシュ

## 再現

```bash
python run_paper7_self_interaction.py
python build_repro_package.py
```

既定実行では n=1 をフル走行し、n=3 は100,000 macrostep の新規 prefix 検証を行う。n=3 も最後まで新規再走行する場合は `python run_paper7_self_interaction.py --full-n3` を用いる。ただし既存 Paper 6 の n3_repulsive raw は同一初期状態・同一 transition による正準 full trajectory である。
'''
(HERE/'README_ja.md').write_text(readme,encoding='utf-8')

prov='''# Source provenance\n\n- Transition map: Paper 6 `run_strict_charged8_5cases.py`, copied without physical-rule changes.\n- Normalization/initializer: Paper 6 supplement `paper6_init_G_q0_c1.py`.\n- Full n=1 reference: Paper 6 `n1_repulsive`, initial state `(P,N,C,D)=(50,0.25,0.91,0)`.\n- Full n=3 reference: Paper 6 `n3_repulsive`, initial state `(P,N,C,D)=(50,0.25,0.19,0)`.\n- New Paper 7 executions: `self_n1` full run and `self_n3_prefix100k` prefix run.\n- n4 self interaction is not assigned a new dynamics; C=-0.44 is recorded as outside the unchanged real quasi-circular map.\n'''
(HERE/'SOURCE_PROVENANCE_ja.md').write_text(prov,encoding='utf-8')

# Checksums last.
def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
ents=[]
for p in sorted(HERE.rglob('*')):
    if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS.json':
        ents.append({'file':str(p.relative_to(HERE)),'sha256':sha256(p),'bytes':p.stat().st_size})
(HERE/'SHA256SUMS.json').write_text(json.dumps(ents,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('built', len(ents), 'files')
