# ABCベースライン定常波C G1 時間位相勾配予備実験 v1

## 目的

`G1`: `τ`方向 `E_read` 勾配による `χ`方向補償候補を、C残渣の時間変化から読む。

ここでは `E_read` 側を位置加速度から定義しない。C残渣の時間変化 `ΔC_memory` から `∂χE_read` proxy を構成し、`R*aχ ≈ -∂χE_read` と整合するかを見る。

## 統合判定

- case_count: `8`
- all_cases_valid: `False`
- single_gauge_only_used: `False`
- g1_time_gradient_preliminary_valid: `False`
- max_match_error: `1.967715446588832e-11`
- min_match_sign_ratio: `0.5869565217391305`
- reaction_controls_rejected: `True`
- early_window_residual_valid: `True`
- full_window_residual_valid: `False`

## ケース別判定

| label | mode | expected | max R*a | max -∂χE | max error | sign ratio | valid |
|---|---|---|---:|---:|---:|---:|---|
| C_mediated_only_persistent | single_C_memory_time_gradient | match | 1.5087216506023456e-10 | 1.5087216506023461e-10 | 4.3658818912747546e-23 | 1.000000 | `True` |
| C_mediated_only_reset | single_C_memory_time_gradient | match | 9.5388813878390045e-20 | 9.5384496399095256e-20 | 4.3820397328424214e-24 | 1.000000 | `True` |
| R2_only_control | reaction_control_without_C_memory, break=1 | reject | 1.8943201400821984e-07 | 0.0000000000000000e+00 | 1.8943201400821984e-07 | 0.000000 | `True` |
| R3_only_control | reaction_control_without_C_memory, break=1 | reject | 1.1839490100581431e-07 | 0.0000000000000000e+00 | 1.1839490100581431e-07 | 0.000000 | `True` |
| combined_minus_R2_R3_residual | combined_minus_R2_R3_residual, break=28 | match | 1.5029864898061050e-10 | 1.5087273853441775e-10 | 1.9677154423536674e-11 | 0.586957 | `False` |
| combined_minus_R2_R3_residual_early_window | combined_minus_R2_R3_residual_steps_1_to_16 | match | 1.5029864898061050e-10 | 1.5087273853441775e-10 | 1.1270077814154828e-11 | 1.000000 | `True` |
| combined_minus_R2_R3_residual_mirrored | combined_minus_R2_R3_residual, break=28 | match | 1.5029864898110185e-10 | 1.5087273853441775e-10 | 1.9677154465888321e-11 | 0.586957 | `False` |
| combined_minus_R2_R3_residual_mirrored_early_window | combined_minus_R2_R3_residual_steps_1_to_16 | match | 1.5029864898110185e-10 | 1.5087273853441775e-10 | 1.1270077819237025e-11 | 1.000000 | `True` |

## 解釈

- C媒介のみのケースで `R*aχ ≈ -∂χE_read` が成立する場合、C残渣の時間変化は位置位相加速度候補と整合する。
- R2/R3 単独で成立しないことは陰性対照である。反力候補を G1 と誤認しないために必要である。
- 合成状態から R2/R3 を差し引いた residual は、full window と early window に分けて判定する。
- full window で符号が崩れ、early window で整合する場合、G1 は定常的な同一視ではなく、準安定遷移窓に限って残る候補として扱う。
- これは標準重力の導出ではなく、`G1` 候補の読出し系列整合性の予備検査である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_g1_time_gradient_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_g1_time_gradient_cases_v1.csv` |
| rows CSV | `abc_baseline_stationary_wave_g1_time_gradient_rows_v1.csv` |
