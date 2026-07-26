# 論文7形式・無seed自然軌道 図3/新図4

論文7の原本エンジンを変更せず、状態への明示的介入をすべて止めた
単一自然軌道から、図3用個別占有列と新図4用自然観測量を同時記録する。

## 介入条件

- 初期kernel seed: OFF
- 初期状態: `Z0 = v.copy()`
- kernel seed `g` は生成せず、乱数も消費しない
- 準安定域の横摂動seed: OFF
- Benettin再射影・再規格化・再投入: OFF
- 観測量から状態更新へのフィードバック: OFF

親状態 `v` の生成規則と、力学計算用の冪反復warm-startは既存条件Aと
同じに保つ。したがって、既存の `condition_A_no_seed.csv` と同じ自然軌道になる。

## 追加記録列

- `direction_1_occupation`
- `direction_2_occupation`
- `direction_3_occupation`
- `direction_4_occupation`
- `other_rotating_occupation`
- `occupation_sum`
- `plane_1_occupation`
- `plane_2_occupation`

既存条件Aの `f_outside_parent`, `q1`〜`q4`, `rank_Q` など27列も
同じCSVに保持する。

## 実行

```bash
python3 run_seedless_natural_figures3_4_v1.py verify
python3 run_seedless_natural_figures3_4_v1.py run 5
python3 run_seedless_natural_figures3_4_v1.py compare-existing 5
python3 run_seedless_natural_figures3_4_v1.py run 40
python3 run_seedless_natural_figures3_4_v1.py compare-existing 40
python3 run_seedless_natural_figures3_4_v1.py run 300
python3 run_seedless_natural_figures3_4_v1.py compare-existing 300
python3 run_seedless_natural_figures3_4_v1.py figures
```

`compare-existing` は追加前の27列を既存条件A CSVと文字列単位で照合し、
追加した個別占有成分の和が1に閉じることも確認する。

### t=110000 長時間測定

```bash
python3 run_seedless_natural_long_horizon_v1.py run 40
python3 run_seedless_natural_long_horizon_v1.py compare-prefix 40
python3 run_seedless_natural_long_horizon_v1.py run 300
python3 run_seedless_natural_long_horizon_v1.py compare-prefix 300
python3 run_seedless_natural_long_horizon_v1.py analyze
python3 make_initial_geometric_growth_residual_v1.py
```

長時間ラッパは同じ無seed初期状態と発展則を再利用し、観測を
`f`, `q1`〜`q4` と閉鎖診断に限定する。状態 `Z` と冪反復warm-start
`wp` を `t=55000,110000` で保存する。

## 出力

- 時系列: `outputs/raw/N*/paper7_long_timeseries.csv`
- 要約: `outputs/summary/N*_5color_meta.json`
- 論文7形式の図3: `outputs/figures/figure3_*`
- 無介入自然図4: `outputs/figures/figure4_seedless_natural_f_q3_q4_compare.*`
- 準安定期間min-max図:
  `outputs/figures/figure4_seedless_natural_minmax_metastable_compare.*`
- 終端50000〜55000 min-max図:
  `outputs/figures/figure4_seedless_natural_minmax_late_50000_55000_compare.*`
- min-max前の実振幅と増減反転数:
  `outputs/summary/figure4_seedless_natural_minmax_ranges.json`
- N=40,300の漸近曲線中心・残差20倍図:
  `outputs/figures/figure4_seedless_asymptotic_centered_x20_N40_N300.*`
- 指数漸近曲線の係数・適合度・残差振幅:
  `outputs/summary/figure4_seedless_asymptotic_centered_x20_fit.json`
- 比較結果: `comparison/compare_existing_A_N*.json`
- 長時間CSV:
  `outputs/long_horizon_110000/raw/N*_seedless_f_q3_q4_t110000.csv`
- 長時間状態:
  `outputs/long_horizon_110000/checkpoints/N*_state_t*.npz`
- 単一指数残差20倍図:
  `outputs/long_horizon_110000/figures/figure_long_horizon_one_exp_residual_x20.*`
- 単一指数・二重指数の基準線感度図:
  `outputs/long_horizon_110000/figures/figure_long_horizon_baseline_sensitivity_x20.*`
- フィットなし終端接近図:
  `outputs/long_horizon_110000/figures/figure_long_horizon_terminal_raw_relative.*`
- 現在標本による初期成長図:
  `outputs/long_horizon_110000/figures/figure_initial_growth_existing_samples.*`
- 長時間解析値:
  `outputs/long_horizon_110000/summary/long_horizon_asymptotic_analysis.json`
- 初期成長の実測極値:
  `outputs/long_horizon_110000/summary/initial_growth_existing_sample_assessment.json`
- `t=0`から幾何級数的成長終端までの生値:
  `outputs/initial_geometric_growth/figures/figure_initial_raw_to_geometric_growth.*`
- 幾何級数的成長曲線中心・対数残差20倍図:
  `outputs/initial_geometric_growth/figures/figure_initial_geometric_growth_residual_x20.*`
- 線形・二次対数成長の基準線感度図:
  `outputs/initial_geometric_growth/figures/figure_initial_growth_baseline_sensitivity_x20.*`
- 初期成長解析値:
  `outputs/initial_geometric_growth/summary/initial_geometric_growth_analysis.json`

min-max図は波形の形を拡大する表示であり、単調緩和も0〜1に拡大される。
振動判定にはJSONに記録した実振幅と増減反転数を併用する。

漸近曲線中心図は `t=20000..55000` に
`c + a exp(-(t-20000)/tau)` を適合し、
`20 * (data-fit)` を表示する。破線のゼロ線が適合した漸近曲線に相当し、
灰線で実残差（1倍）も併記する。したがって、20倍表示だけを物理振幅と
誤認せずに微細構造を読める。

長時間測定では、単一指数残差に見えたN=300の超長期の谷は二重指数基準線で
消失した。さらに `t=55000..110000` の生データはN=40,300とも反転せず
単調に終端値へ接近した。この範囲では、持続する超長期振動よりも、複数の
減衰成分と初期の減衰振動の重なりが支持される。

初期成長の既存標本では、`f>=1e-6` 以後にN=40で `f:4`, `q3:8`,
`q4:6` 個の実測極値を確認した。N=300は100-step標本で全て0個であり、
高周波振動の不存在ではなく現在の時間分解能で未分解という判定に留める。

幾何級数的成長域を `1e-10<=f<=1e-2` と限定すると、N=40の
`t=1125..1725` では `q3` に極大 `t=1650`、極小 `t=1700` が残る。
同区間の `f,q4` は単調である。N=300の `t=2000..3100` は12標本しかなく、
生の `f,q3,q4` は全て単調である。線形対数成長からの大きな残差は、
二次対数成長を許すと大幅に減るため、現在標本だけでは振動と判定しない。

## 検証状況

- N=5: 実行・既存条件Aとの共通27列一致・追加列閉鎖・描画を確認済み
- N=40: 実行・既存条件Aとの共通27列一致・追加列閉鎖を確認済み
- N=300: 図3追加列を含む39列版は未実行
- N=40,300長時間版: `t=110000` まで実行し、既存 `t<=55000` の
  `f,q1..q4` が全行で文字列一致

詳細は `validation_report.md` を参照。
