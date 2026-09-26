# ν第6状態化・恒等写像対照実験 v1

作成日: 2026-09-25

## 目的

第五思考実験で状態外に残っていた対称質量比 \(\nu\) を、第6の乗法状態 \(\mathcal N\) として内部化し、旧5状態系を厳密に再現できるかを対照実験する。

このフォルダでは電荷そのものはまだ導入しない。まず ν だけを内部化し、前実験との完全一致を確認する。

## 主要ファイル

- `run_nu_state_internalization_control_v1.py` — 実験本体
- `plot_nu_state_internalization_control_v1.py` — CSVから図を再生成
- `formulation_and_experiment_design_ja.md` — 定式化・初期条件・合格条件
- `analysis_ja.md` — 結果分析
- `REPRODUCIBILITY_ja.md` — 再現手順
- `nu_state_control_summary.json` — 全結果要約
- `nu_state_control_raw_all12_4000.csv` — 12ケースの生時系列データ
- `nu_state_control_tie_data.csv` — ケース別照合データ
- `nu_state_identity_trace_sampled.csv` — \(\mathcal N\) 恒等保持の時系列
- `nu_state_longrun_three_sampled.csv` — 12000-step長時間照合
- `saved_C1_tie_data.csv` — 保存済みPaper-4 C1との照合
- `negative_control_exp_log_q4.csv` — `exp(log(nu))` 数値往復の負の対照
- `structural_audit.json` — 外部ν引数・恒等係数・microstep監査
- `SHA256SUMS.txt` — 保存物のハッシュ
- `figures/` — PNG/SVG図
- `source_snapshot/` — 比較に使った旧コードと保存済みC1データ

## 最重要結果

12ケース×4000 macro step、3ケース×12000 macro step、ランダム5500 microstepsのすべてで、旧external-ν系と新internal-\(\mathcal N\)系の

```text
(U,P,E,H,Q)
```

は bitwise identical だった。

同時に

```text
Nscale - nu = 0
```

を維持した。

したがって、今回の範囲では ν を第6状態へ移しても前実験の物理出力は変化しない。