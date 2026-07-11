# ABCベースライン定常波C 準安定窓ブラインド抽出予備実験 v1

## 目的

`G1`: `R*aχ ≈ -∂χE_read` 候補の評価窓を、G1誤差や加速度を見ずに先に決める。

窓抽出に使うのは、同一スナップショットから読まれる `C_memory`, `ΔC_memory`, `R2_memory`, `R3_memory`, `Q_raw` だけである。

本実験では、C効果増分が反力記憶に対して十分大きい期間を blind window とし、その後で初めて G1 residual を評価する。

## 統合判定

- case_count: `17`
- single_gauge_only_used: `False`
- blind_signal_ratio_floor: `1e-05`
- all_blind_windows_nonempty: `True`
- all_blind_windows_inside_true_windows: `True`
- all_blind_residual_windows_valid: `True`
- all_C_only_blind_windows_valid: `True`
- all_R2_R3_controls_rejected: `True`
- shifted_window_applicable_count: `15`
- shifted_windows_rejected_or_worse_all_applicable: `True`
- late_window_applicable_count: `15`
- late_windows_rejected_or_worse_all_applicable: `True`
- min_blind_window_length: `5`
- max_blind_window_length: `46`
- max_blind_residual_error: `6.69638462219788e-11`
- blind_metastable_window_preliminary_valid: `True`

## ケース別サマリー

| case | blind window | true window | inside | blind valid | shifted control | late control | R2/R3 rejected | max error |
|---|---:|---:|---|---|---|---|---|---:|
| base | 19 | 27 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 1.3023989807560276e-11 |
| C_decay_0_78 | 7 | 11 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 4.6176043169160496e-12 |
| C_decay_0_86 | 11 | 16 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 7.7647453854253427e-12 |
| C_decay_0_97 | 46 | 46 | `True` | `True` | `not_applicable` | `not_applicable` | `True` | 1.9655368855409219e-11 |
| R_decay_0_78 | 27 | 43 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 4.2019267513608912e-12 |
| R_decay_0_96 | 15 | 20 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 1.8781756211222797e-11 |
| C_source_half | 12 | 21 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 8.5204621366729347e-12 |
| C_source_double | 26 | 34 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 1.6072064287883415e-11 |
| C_return_half | 12 | 21 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 8.5204621366729347e-12 |
| C_return_double | 26 | 34 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 1.6072064287883415e-11 |
| reaction_quarter | 34 | 46 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 1.1329247345318230e-12 |
| reaction_half | 26 | 42 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 4.0178910444296286e-12 |
| reaction_double | 12 | 15 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 3.4082195218814314e-11 |
| slow_C_weak_reaction | 46 | 46 | `True` | `True` | `not_applicable` | `not_applicable` | `True` | 1.2197241151235587e-12 |
| slow_C_strong_reaction | 28 | 30 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 6.6963846221978794e-11 |
| fast_C_weak_reaction | 11 | 20 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 4.8440303336684058e-13 |
| fast_C_strong_reaction | 5 | 7 | `True` | `True` | `rejected_or_worse` | `rejected_or_worse` | `True` | 1.2284171388773090e-11 |

## 解釈

- blind window は、G1 の当たり外れを見ずに、C残渣遅延と反力記憶の比から決めた窓である。
- blind window 内で G1 residual が成立し、shifted/late window が悪化するなら、窓選択は後付けではなく、C残渣遅延から先に指定できる。
- R2/R3 単独が同じ blind window で棄却されるなら、反力候補を G1 と誤認していない。
- 本実験の主張は、G1 が全期間成立することではない。G1 が、観測前に指定可能な準安定窓で読めることである。
- これは標準重力の導出ではなく、加速度読出しへ進む前の窓選択バイアス検査である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_blind_metastable_window_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_blind_metastable_window_cases_v1.csv` |
| signal CSV | `abc_baseline_stationary_wave_blind_metastable_window_signal_rows_v1.csv` |
| eval CSV | `abc_baseline_stationary_wave_blind_metastable_window_eval_rows_v1.csv` |
