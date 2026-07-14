# 白猫・黒猫・灰色猫 D選択境界 予備実験結果 v1

## 1. 実験条件

```text
target = gray_metastable candidates only
d_steps = 2048
pre_steps_values_count = 73
D_gain_values = (0.0, 0.005, 0.01, 0.015, 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.0725, 0.075, 0.0775, 0.08, 0.0825, 0.085, 0.0875, 0.09, 0.095, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0)
C_modes = ('record_only', 'weak_C_window')
```

各観測開始ステップごとにDなし対照を並走させ、Dなしで選択しない条件に限ってD起因選択と判定した。

## 2. 全体集計

total_rows = 15768
boundary_points = 438
selection_possible_boundary_points = 438
no_selection_boundary_points = 0
min_D_gain_overall = 0.0225
max_min_D_gain_overall = 0.065

## 3. 候補別境界

| case_id | boundary_points | selection_possible | no_selection | min_gain | max_min_gain | sign_counts |
|---|---:|---:|---:|---:|---:|---|
| gray_metastable_0_eps0.01_phi0_s0.01_g0 | 146 | 146 | 0 | 0.065 | 0.065 | {'A': 90, 'B': 56} |
| gray_metastable_1_eps0.01_phi1_s0.01_g0 | 146 | 146 | 0 | 0.065 | 0.065 | {'A': 90, 'B': 56} |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 146 | 146 | 0 | 0.0225 | 0.0225 | {'B': 114, 'A': 32} |

## 4. 低しきい値の代表点

| case_id | pre | C_mode | S_start | C_sign | min_D_gain | sign | S_after | baseline |
|---|---:|---|---:|---|---:|---|---:|---|
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 0 | record_only | 0 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 0 | weak_C_window | 0 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 1 | record_only | -0.0015498 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 1 | weak_C_window | -0.0015498 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 1600 | record_only | 0.00222348 | gray | 0.0225 | A | 0.952642 | gray_kept_eigen |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 1600 | weak_C_window | 0.002311 | gray | 0.0225 | A | 0.952642 | gray_kept_eigen |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 2 | record_only | -0.00309644 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 2 | weak_C_window | -0.00309646 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 3 | record_only | -0.00463988 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 3 | weak_C_window | -0.00463993 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 4 | record_only | -0.00618007 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 4 | weak_C_window | -0.00618016 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 5 | record_only | -0.00771695 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 5 | weak_C_window | -0.0077171 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 6 | record_only | -0.00925049 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 6 | weak_C_window | -0.00925071 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 7 | record_only | -0.0107806 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 7 | weak_C_window | -0.0107809 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 8 | record_only | -0.0123073 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 8 | weak_C_window | -0.0123077 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 9 | record_only | -0.0138305 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 9 | weak_C_window | -0.013831 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 10 | record_only | -0.0153502 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 10 | weak_C_window | -0.0153508 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 11 | record_only | -0.0168662 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 11 | weak_C_window | -0.016867 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 550 | record_only | 0.0173109 | gray | 0.0225 | A | 0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 550 | weak_C_window | 0.0174219 | gray | 0.0225 | A | 0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 12 | record_only | -0.0183787 | gray | 0.0225 | B | -0.952642 | unresolved |
| gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002 | 12 | weak_C_window | -0.0183796 | gray | 0.0225 | B | -0.952642 | unresolved |

## 5. 選択しなかった代表点

| case_id | pre | C_mode | S_start | C_sign | baseline |
|---|---:|---|---:|---|---|

## 6. 判定

灰色猫準安定候補に対して、観測開始ステップとD利得の境界表を得た。
今回の掃引範囲では、全ての観測開始ステップでD起因選択が可能だった。
候補ごとの最小D利得は二段に分かれ、弱いしきい値候補では `0.0225`、強いしきい値候補では `0.065` だった。
