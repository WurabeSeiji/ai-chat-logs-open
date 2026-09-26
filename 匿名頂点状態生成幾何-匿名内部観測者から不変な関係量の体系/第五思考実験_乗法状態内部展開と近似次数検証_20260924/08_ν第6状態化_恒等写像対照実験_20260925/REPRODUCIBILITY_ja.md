# 再現手順 — ν第6状態化・恒等写像対照実験

作成日: 2026-09-25

## 1. 実験本体

```bash
python3 run_nu_state_internalization_control_v1.py
```

このスクリプトには、旧external-ν系、新internal-\(\mathcal N\)系、12ケース対照、長時間対照、保存済みC1照合、負の初期化対照、静的監査の実装を含めている。

主要出力:

- `nu_state_control_raw_all12_4000.csv`
- `nu_state_control_tie_data.csv`
- `nu_state_identity_trace_sampled.csv`
- `nu_state_longrun_three_sampled.csv`
- `saved_C1_tie_data.csv`
- `negative_control_exp_log_q4.csv`
- `nu_state_control_summary.json`
- `structural_audit.json`

実行環境によっては全対照を一括で実行すると数分程度かかる。

## 2. 図の再生成

```bash
python3 plot_nu_state_internalization_control_v1.py
```

`figures/` に PNG と SVG を生成する。

## 3. 旧実装の固定参照

`source_snapshot/` に今回の比較に使用した旧実装と保存済みC1データを固定した。

- `audit_two_paths_rk4.py`
- `audit_method_B_strict.py`
- `experiment02_sn_2p5pn_rr_C1_N12.csv`

これらは編集せず、SHA-256 を `nu_state_control_summary.json` と `SHA256SUMS.txt` に保存する。

## 4. 主要検査

### 4.1 12ケース

旧監査と同じ12ケースを4000 macro step比較する。

合格条件:

```text
max abs difference U,P,E,H,Q = 0
max abs difference x,y,r,t   = 0
N drift                       = 0
bitwise failure count         = 0
```

### 4.2 12000-step 長時間対照

3代表ケースを12000 macro step比較する。

合格条件:

```text
max core difference = 0
max N drift         = 0
bitwise failures    = 0
```

### 4.3 5500 microstep ランダム対照

seed=20260924、500初期状態×11相。

合格条件:

```text
core bitwise failures = 0
N bit failures        = 0
phase state failures  = 0
```

### 4.4 保存済みC1

`source_snapshot/experiment02_sn_2p5pn_rr_C1_N12.csv` と旧・新を同時比較する。

旧対新は0差であることを要求する。保存済みCSVとの既存の丸め差は、旧と新で同一値になることを確認する。

## 5. 数値初期化規則

\(\mathcal N=e^{\log\nu}\) は数学的な定義として用いるが、実装では

```text
Nscale = nu
```

とする。

`exp(log(nu))` の数値往復は、\(q=4\) で1 ULP差を作ることが確認されているため、厳密対照には使用しない。

## 6. 重要な範囲制限

この再現手順は ν の内部状態化だけを検証する。電荷自由度、\(\Gamma\)、EM放射、荷電保存力学はまだこの実験へ入れない。