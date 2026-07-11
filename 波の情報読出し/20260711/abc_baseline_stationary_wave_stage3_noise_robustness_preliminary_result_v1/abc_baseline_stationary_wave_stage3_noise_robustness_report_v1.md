# ABCベースライン定常波C Stage III 別読出しノイズ頑健性予備実験 v1

## 目的

Stage III 別読出し照合が、無雑音の同一式を見ているだけではないかを検査する。

同一状態スナップショットから得る位置位相二階差分、`p_read` 差分、G1傾斜読出しに対して、ゼロ平均ゲージ揺らぎと全ゲージ共通バイアスを加える。

ゼロ平均揺らぎは多ゲージ平均で相殺され、共通バイアスは相殺されず検出される必要がある。

## 統合判定

- case_count: `8`
- zero_mean_level_count: `5`
- common_bias_level_count: `3`
- gauge_count: `97`
- single_gauge_only_used: `False`
- zero_mean_multigauge_valid_all: `True`
- common_bias_detection_floor: `1e-13`
- common_bias_detected_all: `True`
- zero_mean_single_gauge_failures_exist: `True`
- zero_mean_max_mean_error_abs: `8.271806125530277e-25`
- common_bias_min_detected_cross_error_abs: `3.475609756200959e-13`
- stage3_noise_robustness_preliminary_valid: `True`

## ノイズケース別サマリー

| mode | level | gauge | single sign failures | max mean err | max cross err | max std | valid/detected |
|---|---:|---:|---:|---:|---:|---:|---|
| zero_mean | 0.0e+00 | 97 | 0 | 8.2718061255302767e-25 | 1.0339757656912846e-23 | 8.2718061255302767e-25 | `True` |
| zero_mean | 5.0e-02 | 97 | 3336 | 2.0679515313825692e-25 | 1.0344088096653177e-23 | 7.7658053790662592e-11 | `True` |
| zero_mean | 2.0e-01 | 97 | 8312 | 4.1359030627651384e-25 | 1.0357079415874169e-23 | 3.1063221516265032e-10 | `True` |
| zero_mean | 5.0e-01 | 97 | 10078 | 2.0679515313825692e-25 | 1.0355746972877145e-23 | 7.7658053790662584e-10 | `True` |
| zero_mean | 1.0e+00 | 97 | 10982 | 3.1019272970738538e-25 | 1.0401583011974800e-23 | 1.5531610758132517e-09 | `True` |
| common_bias | 1.0e-04 | 97 | 0 | 2.3170731707367557e-13 | 3.4756097562009591e-13 | 8.2718061255302767e-25 | `True` |
| common_bias | 1.0e-03 | 97 | 0 | 2.3170731707322062e-12 | 3.4756097561079015e-12 | 8.2718061255302767e-25 | `True` |
| common_bias | 1.0e-02 | 97 | 0 | 2.3170731707317513e-11 | 3.4756097560985950e-11 | 8.2718061255302767e-25 | `True` |

## 解釈

- ゼロ平均ゲージ揺らぎでは、個々のゲージ値が大きく揺れても、多ゲージ平均は元の Stage III 照合へ戻る。
- 共通バイアスは多ゲージ平均では消えないため、読出し器由来の系統偏差として検出される。
- これは標準重力の導出ではなく、準安定傾斜候補の Stage III 読出しが単一ゲージ依存ではないことを確認する予備検査である。

## 出力

| 種類 | ファイル |
|---|---|
| JSON | `abc_baseline_stationary_wave_stage3_noise_robustness_preliminary_result_v1.json` |
| cases CSV | `abc_baseline_stationary_wave_stage3_noise_robustness_cases_v1.csv` |
| eval CSV | `abc_baseline_stationary_wave_stage3_noise_robustness_eval_rows_v1.csv` |
| noisy gauge CSV | `abc_baseline_stationary_wave_stage3_noise_robustness_noisy_gauge_rows_v1.csv` |
