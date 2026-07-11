# ABC Multigauge Readout Integration Summary v1

## Purpose

This summary collects the executed ABC multigauge readout experiments for the 2026-07-11 series.
It is an index of executed numerical evidence, not an additional physical assumption.

## Aggregate Verdict

- experiment_count: `9`
- all_experiments_valid: `True`
- single_gauge_only_used_any: `False`
- integration_summary_valid: `True`

## Experiment Table

| experiment | purpose | count | valid | single gauge only |
|---|---|---:|---|---|
| single_collision_multigauge_readout | 単回ABC衝突で p/E/R を多ゲージ干渉読出しする |  | `True` | `False` |
| multi_collision_multigauge_readout | 対称ABC衝突の反復で p/E/R 読出しを維持する | 8 | `True` | `False` |
| readout_robustness_sweep | 複数の読出し器構成で p/E/R 再構成が安定する | 5 | `True` | `False` |
| asymmetric_amplitude_diagnostic | 非対称Rで単純反転が保存を破ることを検出する | 8 | `True` | `False` |
| generalized_elastic_collision_readout | 非対称Rで R*p と R*p^2 を保存する一般化写像を読む | 8 | `True` | `False` |
| generalized_velocity_sweep | 非単位・非対称位相勾配でも一般化写像が成立する | 9 | `True` | `False` |
| generalized_multi_collision | 一般化写像を複数回AB衝突へ反復適用する | 4 | `True` | `False` |
| generalized_noise_robustness | ゼロ平均読出しノイズの相殺と共通バイアス検出を確認する | 4 | `True` | `False` |
| generalized_extreme_R_sweep | 極端なR比でも一般化写像と読出しが維持されるか調べる | 12 | `True` | `False` |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_readout_integration_summary_result_v1.json` |
| CSV | `abc_multigauge_readout_integration_summary_v1.csv` |
