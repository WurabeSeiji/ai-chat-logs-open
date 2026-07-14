# 白猫・黒猫・灰色猫 D観測応答 予備実験結果 v1

## 1. 実験条件

```text
Stage 1: AB準安定界面探索済み
Stage 2: C読出し窓確認済み
Stage 3: D強観測応答
pre_steps_values = (0, 1, 2, 5, 10, 20, 50, 100, 250, 500, 1000, 2000)
d_steps = 2048
c_modes = ('record_only', 'weak_C_window')
g_D_values = (0.0, 0.05, 0.1, 0.2, 0.5, 1.0)
d_backaction_scale_values = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
```

Dあり条件は、同じD開始状態からDなし対照を並走させて判定した。
Dなし対照でも同じ白猫または黒猫選択が起こる場合は、D起因選択とは数えない。

## 2. 全体集計

total_cases = 7056
D_induced_selection_count = 2016

### D結果分類

| D_outcome | count |
|---|---:|
| black_selected | 772 |
| gray_kept_eigen | 1062 |
| unresolved | 3978 |
| white_selected | 1244 |

### C読出し符号との対応

| D_vs_C_agreement | count |
|---|---:|
| D_selected_from_small_S | 3482 |
| gray_kept | 2566 |
| opposite_sign | 348 |
| same_sign | 660 |

## 3. D起因選択の代表例

| kind | pre | C_mode | eps | phi/pi | s0 | gain | S_start | g_D | D_gain | outcome | baseline | S_mean_after_D | agreement |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---|
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 1 | black_selected | unresolved | -1 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 1 | black_selected | unresolved | -1 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 0.5 | black_selected | unresolved | -1 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 0.5 | black_selected | unresolved | -1 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.5 | 0.5 | black_selected | unresolved | -0.999999 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.5 | 0.5 | black_selected | unresolved | -0.999999 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.5 | 0.25 | black_selected | unresolved | -0.999878 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.5 | 0.25 | black_selected | unresolved | -0.999878 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 0.2 | black_selected | unresolved | -0.99977 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 0.2 | black_selected | unresolved | -0.99977 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.2 | 0.2 | black_selected | unresolved | -0.999714 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.2 | 0.2 | black_selected | unresolved | -0.999714 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 0.1 | black_selected | unresolved | -0.998586 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 1 | 0.1 | black_selected | unresolved | -0.998586 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.5 | 0.1 | black_selected | unresolved | -0.998518 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.5 | 0.1 | black_selected | unresolved | -0.998518 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.2 | 0.1 | black_selected | unresolved | -0.998305 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.2 | 0.1 | black_selected | unresolved | -0.998305 | D_selected_from_small_S |
| gray_metastable | 0 | record_only | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.1 | 0.1 | black_selected | unresolved | -0.997917 | D_selected_from_small_S |
| gray_metastable | 0 | weak_C_window | 0.003 | 0.0833333 | 0 | -0.002 | 0 | 0.1 | 0.1 | black_selected | unresolved | -0.997917 | D_selected_from_small_S |

## 4. 灰色猫固有相の強D応答

| pre | C_mode | S_start | outcome | S_mean_after_D | S_amp_after_D | Q_err |
|---:|---|---:|---|---:|---:|---:|
| 0 | record_only | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 0 | weak_C_window | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 1 | record_only | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 1 | weak_C_window | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 2 | record_only | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 2 | weak_C_window | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 5 | record_only | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 5 | weak_C_window | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 10 | record_only | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 10 | weak_C_window | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 20 | record_only | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |
| 20 | weak_C_window | 0 | gray_kept_eigen | 0 | 0 | 2.22e-16 |

## 5. 大振幅分離領域の強D応答

| pre | C_mode | S_start | C_sign | outcome | S_mean_after_D | agreement |
|---:|---|---:|---|---|---:|---|
| 0 | record_only | 0.06 | A | white_selected | 1 | same_sign |
| 0 | weak_C_window | 0.06 | A | white_selected | 1 | same_sign |
| 1 | record_only | 0.0556528 | A | white_selected | 1 | same_sign |
| 1 | weak_C_window | 0.0556534 | A | white_selected | 1 | same_sign |
| 2 | record_only | 0.0513121 | A | white_selected | 1 | same_sign |
| 2 | weak_C_window | 0.0513133 | A | white_selected | 1 | same_sign |
| 5 | record_only | 0.0383313 | gray | white_selected | 1 | D_selected_from_small_S |
| 5 | weak_C_window | 0.0383338 | gray | white_selected | 1 | D_selected_from_small_S |
| 10 | record_only | 0.0168429 | gray | white_selected | 1 | D_selected_from_small_S |
| 10 | weak_C_window | 0.0168469 | gray | white_selected | 1 | D_selected_from_small_S |
| 20 | record_only | -0.0255206 | gray | black_selected | -1 | D_selected_from_small_S |
| 20 | weak_C_window | -0.025517 | gray | black_selected | -1 | D_selected_from_small_S |

## 6. 判定

灰色猫固有相、灰色猫準安定相、大振幅分離領域を同じD写像で比較できるデータが得られた。
次段階では、D起因選択が現れた準安定候補を中心に、D結合強度と観測開始位相の境界を細かく調べる。
