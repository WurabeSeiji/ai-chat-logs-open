# 白猫・黒猫・灰色猫 C読出し窓 予備実験結果 v1

## 1. 実験条件

```text
C = on
D = off
steps = 4096
g_C_values = (0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
backaction_scale_values = (0.0, 1e-05, 5e-05, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01)
readout_kappa = 0.02
```

## 2. 判定

AB二体準安定候補にCを加え、A/B配分を読めるが一方選択を起こさない結合窓を確認した。

C_window_count = 247 / 1080
C_informative_window_count = 144 / 1080
C_nonzero_backaction_window_count = 114 / 1080

## 3. C後の相分類

| phase_after_C | count |
|---|---:|
| gray_eigen | 90 |
| gray_metastable | 786 |
| large_oscillation | 141 |
| natural_selection | 2 |
| unstable_or_drifting | 61 |

## 4. 上位C読出し窓

| kind | epsilon | phi/pi | s0 | base_gain | g_C | c_gain | C_rel_err | C_bias_delta | phase_after_C |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| gray_metastable | 0.001 | 0 | 0.01 | 0 | 1 | 0 | 0.0196078 | 0 | gray_metastable |
| gray_metastable | 0.001 | 1 | 0.01 | 0 | 1 | 0 | 0.0196078 | 0 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0 | -0.002 | 1 | 5e-05 | 0.0196078 | 0.00122196 | gray_metastable |
| gray_metastable | 0.003 | 0 | 0.01 | 0 | 1 | 0 | 0.0196078 | 0 | gray_metastable |
| gray_metastable | 0.003 | 1 | 0.01 | 0 | 1 | 0 | 0.0196078 | 0 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0 | -0.002 | 1 | 0 | 0.0196078 | 0 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.01 | -0.002 | 1 | 0 | 0.0196078 | 0 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.000220996 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.001 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.000222243 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.01 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.000234393 | gray_metastable |
| gray_metastable | 0.003 | 0 | 0.01 | 0 | 1 | 1e-05 | 0.0196078 | 0.00034922 | gray_metastable |
| gray_metastable | 0.001 | 1 | 0.01 | 0 | 1 | 1e-05 | 0.0196078 | 0.000367551 | gray_metastable |
| gray_metastable | 0.01 | 1 | 0.01 | 0 | 1 | 1e-05 | 0.0196078 | 0.000398998 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.001 | -0.002 | 1 | 5e-05 | 0.0196078 | 0.00113414 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.0001 | -0.002 | 1 | 5e-05 | 0.0196078 | 0.00120902 | gray_metastable |

## 5. 非ゼロCバックアクションを持つ読出し窓

| kind | epsilon | phi/pi | s0 | base_gain | g_C | c_gain | C_rel_err | C_bias_delta | phase_after_C |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| gray_metastable | 0.003 | 0.0833333 | 0 | -0.002 | 1 | 5e-05 | 0.0196078 | 0.00122196 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.000220996 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.001 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.000222243 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.01 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.000234393 | gray_metastable |
| gray_metastable | 0.003 | 0 | 0.01 | 0 | 1 | 1e-05 | 0.0196078 | 0.00034922 | gray_metastable |
| gray_metastable | 0.001 | 1 | 0.01 | 0 | 1 | 1e-05 | 0.0196078 | 0.000367551 | gray_metastable |
| gray_metastable | 0.01 | 1 | 0.01 | 0 | 1 | 1e-05 | 0.0196078 | 0.000398998 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.001 | -0.002 | 1 | 5e-05 | 0.0196078 | 0.00113414 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.0001 | -0.002 | 1 | 5e-05 | 0.0196078 | 0.00120902 | gray_metastable |
| gray_metastable | 0.003 | 1 | 0.01 | 0 | 1 | 5e-05 | 0.0196078 | 0.00180792 | gray_metastable |
| gray_metastable | 0.001 | 0 | 0.01 | 0 | 1 | 5e-05 | 0.0196078 | 0.00190118 | gray_metastable |
| gray_metastable | 0.01 | 0 | 0.01 | 0 | 1 | 5e-05 | 0.0196078 | 0.00207535 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.001 | -0.002 | 1 | 0.0001 | 0.0196078 | 0.00261115 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.0001 | -0.002 | 1 | 0.0001 | 0.0196078 | 0.00272785 | gray_metastable |
| gray_metastable | 0.003 | 0.0833333 | 0.0001 | -0.002 | 1 | 1e-05 | 0.0196078 | 0.00022113 | gray_metastable |

## 6. 次段階

C読出し窓を固定し、Dを加えた観測実験で、灰色猫固有相・灰色猫準安定相・大振幅分離領域の応答を比較する。
