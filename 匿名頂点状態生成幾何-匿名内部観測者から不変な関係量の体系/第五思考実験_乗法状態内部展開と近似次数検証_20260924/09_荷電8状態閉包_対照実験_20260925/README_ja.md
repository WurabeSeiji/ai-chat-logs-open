# 荷電8状態閉包・対照実験 v1

## 目的

重力-only で内部化した第6状態 `N=eta` に加え、荷電準円軌道に必要な

- `C = Z = 1 - lambda_A lambda_B`
- `D = (lambda_A - lambda_B)^2`

を第7・第8状態として内部化し、外部物理引数なしの状態閉包を検査する。

## 持続状態

`(U, P, E, H, Q, N, C, D)`

`N,C,D` は11相 one-hot 相互作用の恒等行で保存する。外部コピーは使わない。

## 主要結果

- 外部パラメータ対照 vs 内部8状態: 4ケース×4000 macro step で bitwise mismatch 0。
- ランダム500状態×11 microstep: mismatch 0、N/C/D identity failure 0。
- 明示11-microstep vs 折り畳みmacro: 12000 macro step で8状態 bitwise mismatch 0。
- 6状態および1個だけ荷電状態を追加した7状態は反例により非閉包。
- 基準軌道 `r=50 -> 20` の交点で `|Delta t|=1.46e-8`, `|Delta phi|=8.44e-11 rad`。
- 全軌道で N/C/D drift = 0。

## 主要ファイル

- `run_charged_eight_state_closure_v1.py` 実験本体
- `plot_charged_eight_state_closure_v1.py` 図生成
- `charged_eight_state_summary.json` 結果要約
- `charged_eight_state_benchmark_sampled.csv` 基準軌道のサンプル時系列
- `closure_counterexamples.csv` 6/7状態非閉包反例
- `convergence.csv` 刻み依存確認
- `analysis_ja.md` 詳細分析
- `REPRODUCIBILITY_ja.md` 再現手順
- `reference_raw.csv`, `reference_metadata.json` 独立基準の固定スナップショット
- `figures/` PNG/SVG

## 注意

本実験は leading-order adiabatic quasi-circular charged benchmark に限定する。full Einstein-Maxwell closure および R4 完全乗法閉包を主張しない。
