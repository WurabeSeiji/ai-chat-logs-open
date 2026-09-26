# 再現性記録

## 正規生成器
`run_strict_charged8_final.py`

- 永続物理状態: `(U,P,E,H,Q,N,C,D)`
- full microstate: 上記8状態 + RK work registers + one-hot `q`
- 生成経路: `transition(z)` のみ
- 1 macrostep = `transition(z)` を正確に11回適用
- external physical runtime arguments: なし
- fast/collapsed/shortcut generator: なし
- reference data import: generator側にはなし

## 実験条件
- `STEPS_PER_ORBIT = 4000`
- `P0 = 50`
- `N0 = 0.25`
- `C0 = 1.09`
- `D0 = 0.36`
- 終了条件: `P <= 20`

## 出力
- `raw_macro_trajectory.csv`: 全macrostep raw data
- `raw_first_macro_microsteps.csv`: 最初のmacrostepの全microstate
- `generator_summary.json`: 生成結果
- `preexec_audit.json`: A1-A8実行前監査
- `reference_comparison.json`: 生成後の独立基準比較
- `analysis_ja.md`: 分析
- `figures/*.svg`: 図
- `SHA256SUMS.txt`: ハッシュ

## 実行順序
1. `python run_strict_charged8_final.py`
2. 生成完了後のみ `python postprocess_strict_charged8.py`

独立基準は2の後処理でのみ読み、生成器へのフィードバックは行わない。
