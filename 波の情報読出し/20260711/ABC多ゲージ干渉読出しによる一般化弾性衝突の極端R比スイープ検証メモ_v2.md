# ABC多ゲージ干渉読出しによる一般化弾性衝突の極端R比スイープ検証メモ v1

**副題:** `R_B/R_A` の大きな非対称性に対する `R*p`・`R*p^2` 保存写像の境界検証  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、一般化弾性衝突写像を、`R_read` が大きく非対称な条件へ適用した。

前段階では、非対称振幅および非対称位相勾配に対して、

```text
R_A p_A + R_B p_B
R_A p_A^2 + R_B p_B^2
```

が保存されることを確認した。

本検証では、`R_B/R_A` を広く掃引し、`R` 比が大きく偏った場合でも、多ゲージ干渉読出しと一般化衝突写像が安定に成立するかを調べる。

ここで `R` は標準物理の質量そのものではない。位相系内で安定軸として読まれる振幅二乗量である。

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v2.py
```

出力先は次である。

```text
abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v2.py
```

---

## 3. 実験条件

初期位相勾配は全ケースで、

```text
p_A=1.0
p_B=-0.5
```

とした。

振幅比を変えることで、

```text
R_A=A_A^2
R_B=A_B^2
```

を変えた。

掃引範囲は次である。

| 量 | 値 |
|---|---:|
| 最小 `R_B/R_A` | `0.015625` |
| 最大 `R_B/R_A` | `64.0` |
| 最大動的範囲 | `64.0` |
| ケース数 | `12` |

---

## 4. 判定量

判定量は次である。

### 4.1 多ゲージ個別読出し

```text
p_read
E_read
R_read
```

が複数ゲージから再構成できるかを確認する。

### 4.2 `R*p` 保存

```text
P_R = R_A p_A + R_B p_B
```

### 4.3 `R*p^2` 保存

```text
K_R = R_A p_A^2 + R_B p_B^2
```

### 4.4 相対位相勾配反転

```text
(p_A' - p_B') = -(p_A - p_B)
```

### 4.5 `R_total` 保存

```text
R_A+R_B
```

が衝突前後で保存されるかを確認する。

---

## 5. 実行結果

全体判定は次である。

```text
case_count: 12
all_cases_valid: true
collision_reached_all_cases: true
individual_readout_valid_all_cases: true
P_R_preserved_all_cases: true
K_R_phase_preserved_all_cases: true
relative_gradient_flipped_all_cases: true
E_tau_R_preserved_all_cases: true
R_total_preserved_all_cases: true
single_gauge_only_used: false
extreme_R_sweep_valid: true
```

主要誤差は次である。

| 量 | 値 |
|---|---:|
| `max_R_dynamic_range` | `64.0` |
| `min_R_ratio_B_over_A` | `0.015625` |
| `max_R_ratio_B_over_A` | `64.0` |
| `max_p_abs_error` | `1.000310945187266e-13` |
| `max_E_abs_error` | `2.2315482794965646e-14` |
| `max_R_abs_error` | `2.842170943040401e-14` |
| `max_R_gauge_std` | `1.7404671430534633e-14` |
| `max_within_particle_separation_ratio_time` | `7.143466159887607e-30` |
| `max_P_R_conservation_error` | `6.465938895416912e-13` |
| `max_K_R_phase_conservation_error` | `1.2789769243681803e-12` |
| `max_relative_flip_error` | `2.6867397195928788e-14` |
| `max_E_tau_R_conservation_error` | `8.881784197001252e-16` |
| `max_R_total_conservation_error` | `3.552713678800501e-15` |

---

## 6. ケース別結果

| case | `R_B/R_A` | dynamic range | `p_A,p_B` after | `R*p` err | `R*p^2` err | valid |
|---|---:|---:|---|---:|---:|---|
| `c01_A1.000_B0.125_u1.00_v-0.50` | `1.5625000000000000e-02` | `6.4000000000000000e+01` | `0.95384615 / 2.4538462` | `1.0103029524088925e-14` | `1.9984014443252818e-14` | `true` |
| `c02_A1.000_B0.250_u1.00_v-0.50` | `6.2500000000000000e-02` | `1.6000000000000000e+01` | `0.82352941 / 2.3235294` | `1.0436096431476471e-14` | `1.4432899320127035e-14` | `true` |
| `c03_A1.000_B0.500_u1.00_v-0.50` | `2.5000000000000000e-01` | `4.0000000000000000e+00` | `0.4 / 1.9` | `5.4400928206632670e-15` | `8.4376949871511897e-15` | `true` |
| `c04_A1.000_B1.000_u1.00_v-0.50` | `1.0000000000000000e+00` | `1.0000000000000000e+00` | `-0.5 / 1` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `true` |
| `c05_A1.000_B2.000_u1.00_v-0.50` | `4.0000000000000000e+00` | `4.0000000000000000e+00` | `-1.4 / 0.1` | `4.5519144009631418e-15` | `2.9531932455029164e-14` | `true` |
| `c06_A1.000_B4.000_u1.00_v-0.50` | `1.6000000000000000e+01` | `1.6000000000000000e+01` | `-1.8235294 / -0.32352941` | `2.3714363805993344e-13` | `2.0161650127192843e-13` | `true` |
| `c07_A1.000_B8.000_u1.00_v-0.50` | `6.4000000000000000e+01` | `6.4000000000000000e+01` | `-1.9538462 / -0.45384615` | `1.3145040611561853e-13` | `1.9895196601282805e-13` | `true` |
| `c08_A0.250_B1.000_u1.00_v-0.50` | `1.6000000000000000e+01` | `1.6000000000000000e+01` | `-1.8235294 / -0.32352941` | `1.4821477378745840e-14` | `1.2601031329495527e-14` | `true` |
| `c09_A0.500_B1.000_u1.00_v-0.50` | `4.0000000000000000e+00` | `4.0000000000000000e+00` | `-1.4 / 0.1` | `1.1379786002407855e-15` | `7.3829831137572910e-15` | `true` |
| `c10_A2.000_B1.000_u1.00_v-0.50` | `2.5000000000000000e-01` | `4.0000000000000000e+00` | `0.4 / 1.9` | `2.1760371282653068e-14` | `3.3750779948604759e-14` | `true` |
| `c11_A4.000_B1.000_u1.00_v-0.50` | `6.2500000000000000e-02` | `1.6000000000000000e+01` | `0.82352941 / 2.3235294` | `1.6697754290362354e-13` | `2.3092638912203256e-13` | `true` |
| `c12_A8.000_B1.000_u1.00_v-0.50` | `1.5625000000000000e-02` | `6.4000000000000000e+01` | `0.95384615 / 2.4538462` | `6.4659388954169117e-13` | `1.2789769243681803e-12` | `true` |

---

## 7. 解釈

本検証では、`R_B/R_A` が `0.015625` から `64.0` まで変化しても、全ケースで相互作用セルに到達し、一般化弾性衝突写像が成立した。

特に、最大動的範囲 `64.0` の条件でも、

```text
R*p
R*p^2
relative p flip
R_total
```

が保存された。

この結果は、`R` が等しい場合だけでなく、非対称な `R` 読出しを持つ局所波同士でも、`R` 重み付き保存量として衝突読出しが構成できることを示す。

---

## 8. 本検証で主張しないこと

本検証は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準質量比の定量予言 | `R` は位相系内の読出し量である |
| 任意に大きい `R` 比での成立 | 本検証範囲は `0.015625 <= R_B/R_A <= 64.0` である |
| 実在粒子衝突の定量再現 | 本検証は内部読出し写像の数値構成である |
| 単一ゲージ測定の成立 | 多ゲージ干渉読出しを用いる |

---

## 9. 結論

本検証では、一般化弾性衝突写像を極端な `R` 非対称条件へ適用した。

`R_B/R_A` は `0.015625` から `64.0` まで掃引され、12ケースすべてで判定が `true` となった。

最大 `R*p` 保存誤差は `6.465938895416912e-13`、最大 `R*p^2` 保存誤差は `1.2789769243681803e-12`、最大相対位相勾配反転誤差は `2.6867397195928788e-14` であった。

したがって、本検証の範囲では、`R` 比が大きく非対称な場合でも、多ゲージ干渉読出しによる `R` 重み付き一般化弾性衝突写像は保存関係を維持した。

---

# 付録A. 出力ファイル

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_report_v2.md](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_report_v2.md) |
| 結果 JSON | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2.json](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2.json) |
| ケース CSV | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_cases_v2.csv](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_cases_v2.csv) |
| stage quantity CSV | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_stage_quantities_v2.csv](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_stage_quantities_v2.csv) |
| gauge CSV | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_gauge_rows_v2.csv](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_gauge_rows_v2.csv) |
| 誤差図 | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_errors_v2.png](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_errors_v2.png) |
| 出力勾配図 | [abc_multigauge_generalized_elastic_collision_extreme_R_sweep_outputs_v2.png](abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_outputs_v2.png) |
