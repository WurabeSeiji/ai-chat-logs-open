# ABCベースライン定常波C 準安定傾斜遅延カーネル注入予備実験 v1

## 目的

遅延カーネル分解器が本当に lag 1 以降の残渣を検出できるかを確認するため、既知の遅延重みを `R*a` 候補へ人工的に注入する。

現行モデルでは lag 0 が支配的であった。そこで本実験では、lag 1, lag 2, 指数尾、純遅延尾、符号反転尾を陽性対照として与え、回帰で既知係数が回収されるかを調べる。

## 統合判定

- kernel_count: `6`
- row_count: `10152`
- single_gauge_only_used: `False`
- all_dominant_lags_recovered: `True`
- all_delayed_fractions_recovered: `True`
- all_coefficients_recovered: `True`
- max_coefficient_error_abs: `9.506478398674303e-19`
- min_lag_kernel_r2: `1.0`
- delayed_positive_controls_detected: `True`
- metastable_delay_kernel_injection_preliminary_valid: `True`

## カーネル別サマリー

| kernel | kind | expected lag | recovered lag | expected delayed | recovered delayed | max coef err | R2 | valid |
|---|---|---:|---:|---:|---:|---:|---:|---|
| instant_control | control | 0 | 0 | 0.0000000000000000e+00 | 4.5820506724848277e-15 | 9.5064783986743026e-19 | 1.0000000000000000e+00 | `True` |
| lag1_dominant | one_step_delay | 1 | 1 | 7.5000000000000000e-01 | 7.5000000000000155e-01 | 4.8999530465866874e-19 | 1.0000000000000000e+00 | `True` |
| lag2_dominant | two_step_delay | 2 | 2 | 8.0000000000000004e-01 | 8.0000000000000071e-01 | 6.5052130349130266e-19 | 1.0000000000000000e+00 | `True` |
| exponential_tail | tail | 1 | 1 | 7.5000000000000000e-01 | 7.5000000000000044e-01 | 2.9730417700753874e-19 | 1.0000000000000000e+00 | `True` |
| pure_delayed_tail | tail_no_current | 1 | 1 | 1.0000000000000000e+00 | 9.9999999999999944e-01 | 4.2678941409788804e-19 | 1.0000000000000000e+00 | `True` |
| alternating_delay | signed_tail | 1 | 1 | 7.5609756097560976e-01 | 7.5609756097561009e-01 | 7.4260160566241350e-19 | 1.0000000000000000e+00 | `True` |

## 解釈

- 既知の遅延尾を入れた場合に lag 1 以降が回収されるなら、前段の `lag 0` 支配は検出器の鈍さではなく現行モデル側の性質と読める。
- 純遅延尾や符号反転尾も回収できるなら、単に正の指数尾だけを見ているのではなく、符号付きカーネルを識別できる。
- 本実験は標準重力の導出ではなく、遅延カーネル分解法の陽性対照である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_preliminary_result_v1.json` |
| rows CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_rows_v1.csv` |
| summaries CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_summaries_v1.csv` |
| coefficients CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_coefficients_v1.csv` |
