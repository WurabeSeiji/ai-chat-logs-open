# 第2予備実験：初期期間 N–分解能掃引（完全無seed）

（§16 の制限に従い、実験名・実行方法・構成・設定場所・出力一覧・再実行・依存のみを記す。結果の要約・意味・解釈・結論は書かない。）

## 1. 実験名

initial_period_N_resolution_scan_no_seed（初期期間 N–分解能掃引、完全無seed）。

## 2. 実行方法

```
cd code
python3 audit_resolution_scan_dependencies_v1.py      # §3 依存監査（read-only）
python3 run_initial_resolution_scan_v1.py ALL          # 全156 execution
python3 compute_prefixed_regressions_v1.py             # §11 固定band回帰
python3 aggregate_resolution_scan_v1.py                # §12 集計 + §14 再現性差分
python3 make_initial_resolution_scan_figures_v1.py     # §13 必須図
```

単一runは `python3 run_initial_resolution_scan_v1.py <N> <p> <Delta_ref> <exec_idx>`。
基準run（分解能OFF）は `python3 run_initial_resolution_scan_v1.py <N> OFF OFF <exec_idx>`。

## 3. フォルダ構成

```
preliminary_02_resolution_scan/
  instructions/   実験仕様書
  code/           実行・監査・回帰・集計・図の各スクリプト
  config/         experiment_manifest.json, source_file_hashes.json
  raw/<run_id>/   timeseries.csv, local_growth.csv, regression_by_fixed_band.csv,
                  run_config.json, run_diagnostics.json, state_vectors.npz, stdout.log, stderr.log
  summary/        all_runs_manifest.csv, all_runs_final_values.csv,
                  all_fixed_band_regressions.csv, all_stop_reasons.csv, all_diagnostics.csv
  figures/        per_run_*.png, cmp_Np_*.png, cmp_pD_*.png（auto/fixedaxis）
  diagnostics/    reproducibility_diff_<run_id>.csv, reproducibility_overview.csv
  logs/           full_sweep_stdout.log, full_sweep_stderr.log
  reports/        resolution_scan_dependency_audit.md, execution_report.md
```

## 4. 設定ファイルの場所

- `config/experiment_manifest.json`（全78 config を事前固定、Δ_actual を含む）
- `config/source_file_hashes.json`（再利用した既存コードの SHA-256）

## 5. 出力ファイル一覧

各 run: `raw/<run_id>/timeseries.csv`（§10 全列）, `local_growth.csv`, `regression_by_fixed_band.csv`,
`run_config.json`, `run_diagnostics.json`, `state_vectors.npz`, `stdout.log`, `stderr.log`。
集計: `summary/*.csv`。図: `figures/*.png`。再現性: `diagnostics/reproducibility_*.csv`。報告: `reports/execution_report.md`。

## 6. 再実行方法

`code/` で §2 の 5 コマンドを順に再実行する。乱数不使用のため各 run は決定論的に再現する。

## 7. 依存パッケージ

Python 3、NumPy、SciPy、Matplotlib。既存エンジン（`run_n_scaling_lowrank_v1.py` ほか）と
第7論文 `run_transverse_stability_v1.py`（polar retraction）を read-only import で再利用する（不変更）。
