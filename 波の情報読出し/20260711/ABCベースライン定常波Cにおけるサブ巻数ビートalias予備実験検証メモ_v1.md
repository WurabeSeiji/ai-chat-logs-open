# ABCベースライン定常波Cにおけるサブ巻数ビート alias 予備実験 v1

## 目的

整数巻数に近い進行位相を、主値位相と unwrap 位相で読むと、巻数以下の残差が巨視的な逆向きビートとして見えるかを確認する。

これは G3 候補、すなわちサブ巻数ビート/alias による見かけの接近・反転を、C媒介の引力様応答から分離するための予備実験である。

## 統合判定

- case_count: `4`
- all_cases_valid: `True`
- single_gauge_only_used: `False`
- subwinding_alias_preliminary_valid: `True`
- max_beat_slope_error: `4.8031023602845835e-14`
- max_gauge_slope_std: `6.941913964043779e-15`

## ケース別判定

| case | winding | residual | true step phase | beat slope | expected reverse | observed reverse | valid |
|---|---:|---:|---:|---:|---|---|---|
| above_integer_winding | 5 | 2.9999999999999999e-02 | 3.1445926535897932e+01 | 3.0000000000000675e-02 | `False` | `False` | `True` |
| below_integer_winding | 5 | -2.9999999999999999e-02 | 3.1385926535897930e+01 | -3.0000000000001810e-02 | `True` | `True` | `True` |
| above_high_winding | 41 | 1.2500000000000001e-02 | 2.5762309759436300e+02 | 1.2499999999951970e-02 | `False` | `False` | `True` |
| below_high_winding | 41 | -1.2500000000000001e-02 | 2.5759809759436303e+02 | -1.2500000000025832e-02 | `True` | `True` | `True` |

## 解釈

- 真の進行位相は全ケースで正方向である。
- しかし整数巻数よりわずかに小さい場合、主値位相で見た beat は負方向へ進む。
- したがって、後続の加速度実験で距離位相縮小が見えた場合、unwrap 系列と beat 系列を分離しないと、引力様ドリフトと alias 反転を混同する。
- 本予備実験は重力的効果を主張しない。G3 を独立した観測方法依存候補として保持するための基準実験である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_subwinding_alias_preliminary_result_v1.json` |
| case CSV | `abc_baseline_stationary_wave_subwinding_alias_cases_v1.csv` |
| gauge CSV | `abc_baseline_stationary_wave_subwinding_alias_gauge_rows_v1.csv` |
