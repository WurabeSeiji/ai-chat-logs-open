# ABC多ゲージ干渉読出しによる一般化弾性衝突の非対称速度スイープ検証メモ v1

**副題:** 非単位・非対称な初期位相勾配に対する `R*p`・`R*p^2` 保存写像の数値検証  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、前段階で構成した `R_read` 重み付き一般化弾性衝突写像が、初期位相勾配

```text
p_A, p_B
```

を `+1,-1` に固定しない場合にも成立するかを検査する。

ここで `p` は標準物理の速度または運動量そのものではない。

本実験における `p` は、空間位相方向の多ゲージ干渉読出し勾配である。

前段階の実験では、非対称振幅条件でも、

```text
R_A p_A + R_B p_B
R_A p_A^2 + R_B p_B^2
```

を保存する写像が構成できた。

ただし、初期条件は主に

```text
p_A=+1
p_B=-1
```

であった。

本実験では、初期位相勾配を複数に変え、対向条件だけでなく、同方向の追いつき衝突も含めて検証する。

---

## 2. 本実験で主張しないこと

本実験は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準力学の速度そのものを導出したこと | `p` は空間位相勾配読出しである |
| 標準運動量・標準エネルギーの完全導出 | 対応写像は別途必要である |
| 任意次元衝突への拡張 | 本実験は1次元位相勾配読出しである |
| 任意の初期条件で衝突すること | 本実験では `p_A > p_B` となるケースを選ぶ |

本実験の主張は、非単位・非対称な初期位相勾配でも、ABC 多ゲージ干渉読出し上で `R` 重み付き一般化衝突写像が保存量を保つかを検査することである。

---

## 3. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1.py
```

出力先は次である。

```text
abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1.py
```

---

## 4. 実験条件

初期位相勾配ケースは次である。

| case | `A_A` | `A_B` | `p_A` before | `p_B` before | 意味 |
|---|---:|---:|---:|---:|---|
| c01 | `1.0` | `1.0` | `1.0` | `-1.0` | 等振幅・基準対向 |
| c02 | `1.0` | `1.0` | `1.4` | `-0.6` | 等振幅・非対称対向 |
| c03 | `1.0` | `1.0` | `0.5` | `-1.7` | 等振幅・B側が速い対向 |
| c04 | `1.0` | `2.0` | `1.0` | `-1.0` | 非対称R・基準対向 |
| c05 | `1.0` | `2.0` | `1.4` | `-0.6` | 非対称R・非対称対向 |
| c06 | `1.0` | `2.0` | `1.2` | `0.2` | 非対称R・同方向追いつき |
| c07 | `2.0` | `1.0` | `0.8` | `-1.5` | A側が大きい非対称R |
| c08 | `1.5` | `1.0` | `1.8` | `-0.2` | A側が大きい非対称R |
| c09 | `1.0` | `3.0` | `0.8` | `-0.4` | B側が非常に大きい非対称R |

すべてのケースで、

```text
p_A > p_B
```

を満たすようにした。

これは、初期位置で A が左、B が右にあるため、相互作用セルへ到達する条件である。

---

## 5. 判定量

本実験の判定量は次である。

### 5.1 `R*p` 保存

```text
P_R = R_A p_A + R_B p_B
```

### 5.2 `R*p^2` 保存

```text
K_R = R_A p_A^2 + R_B p_B^2
```

### 5.3 相対位相勾配の符号反転

```text
(p_A' - p_B') = -(p_A - p_B)
```

### 5.4 `R*E_tau` 保存

```text
E_{\tau,R} = R_A E_A + R_B E_B
```

### 5.5 `R_total` 保存

```text
R_total = R_A + R_B
```

---

## 6. 実行結果

全体判定は次である。

```text
case_count: 9
collision_reached_all_cases: true
individual_readout_valid_all_cases: true
P_R_preserved_all_cases: true
K_R_phase_preserved_all_cases: true
relative_gradient_flipped_all_cases: true
E_tau_R_preserved_all_cases: true
R_total_preserved_all_cases: true
single_gauge_only_used: false
velocity_sweep_generalized_collision_valid: true
```

主要誤差は次である。

| 量 | 値 |
|---|---:|
| `max_p_abs_error` | `4.1033842990145786e-13` |
| `max_E_abs_error` | `2.6423307986078726e-14` |
| `max_R_abs_error` | `1.865174681370263e-14` |
| `max_R_gauge_std` | `8.078732199757252e-15` |
| `max_within_particle_separation_ratio_time` | `2.365357211965411e-28` |
| `max_P_R_conservation_error` | `2.8910207561239076e-13` |
| `max_K_R_phase_conservation_error` | `1.5258905250448151e-12` |
| `max_relative_flip_error` | `1.2434497875801753e-14` |
| `max_E_tau_R_conservation_error` | `4.884981308350689e-15` |
| `max_R_total_conservation_error` | `4.884981308350689e-15` |

---

## 7. ケース別結果

| case | `R_B/R_A` | `p_A,p_B` before | `p_A,p_B` after | `R*p` err | `R*p^2` err | relative flip err |
|---|---:|---|---|---:|---:|---:|
| `c01_A1.00_B1.00_u1.00_v-1.00` | `1.0000000000000000e+00` | `1 / -1` | `-1 / 1` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` |
| `c02_A1.00_B1.00_u1.40_v-0.60` | `1.0000000000000000e+00` | `1.4 / -0.6` | `-0.6 / 1.4` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` |
| `c03_A1.00_B1.00_u0.50_v-1.70` | `1.0000000000000000e+00` | `0.5 / -1.7` | `-1.7 / 0.5` | `1.8052226380405045e-13` | `6.2705396430828841e-13` | `0.0000000000000000e+00` |
| `c04_A1.00_B2.00_u1.00_v-1.00` | `4.0000000000000000e+00` | `1 / -1` | `-2.2 / -0.2` | `1.2878587085651816e-14` | `2.3980817331903381e-14` | `1.2434497875801753e-14` |
| `c05_A1.00_B2.00_u1.40_v-0.60` | `4.0000000000000000e+00` | `1.4 / -0.6` | `-1.8 / 0.2` | `1.8651746813702630e-14` | `1.4210854715202004e-14` | `2.2204460492503131e-15` |
| `c06_A1.00_B2.00_u1.20_v0.20` | `4.0000000000000000e+00` | `1.2 / 0.2` | `-0.4 / 0.6` | `2.4424906541753444e-14` | `1.5543122344752192e-15` | `1.1324274851176597e-14` |
| `c07_A2.00_B1.00_u0.80_v-1.50` | `2.5000000000000000e-01` | `0.8 / -1.5` | `-0.12 / 2.18` | `5.2846615972157451e-14` | `8.8817841970012523e-16` | `1.0214051826551440e-14` |
| `c08_A1.50_B1.00_u1.80_v-0.20` | `4.4444444444444442e-01` | `1.8 / -0.2` | `0.569231 / 2.56923` | `2.8910207561239076e-13` | `1.5258905250448151e-12` | `5.3290705182007514e-15` |
| `c09_A1.00_B3.00_u0.80_v-0.40` | `9.0000000000000000e+00` | `0.8 / -0.4` | `-1.36 / -0.16` | `7.9936057773011271e-14` | `3.8191672047105385e-14` | `4.4408920985006262e-16` |

---

## 8. 解釈

本実験では、初期位相勾配を `+1,-1` に固定せず、複数の非対称条件へ広げた。

結果として、全9ケースで相互作用セルに到達し、多ゲージ干渉読出しが成立した。

さらに、

```text
R*p
R*p^2
R*E_tau
R_total
```

が全ケースで保存された。

また、相対位相勾配

```text
p_A - p_B
```

は、衝突後に符号反転した。

特に、ケース `c06` では、

```text
p_A = 1.2
p_B = 0.2
```

であり、両者が同方向へ進む追いつき衝突である。

この場合でも、

```text
p_A' = -0.4
p_B' = 0.6
```

となり、`R*p` と `R*p^2` は保存された。

したがって、一般化写像は、単純な等速度対向衝突だけでなく、非単位・非対称な位相勾配条件にも対応した。

---

## 9. 本実験の位置づけ

本実験により、ABC 多ゲージ干渉読出しで得られる `R_read` を質量的重みとして使うと、1次元位相勾配読出し上の一般化弾性衝突写像が構成できることが確認された。

ただし、これは標準力学の導出そのものではない。

本実験は、閉じた位相系の内部読出しとして、

```text
質量的量 R_read
運動量的量 p_read
位相勾配二乗型保存量 R*p^2
```

を同時に扱えることを示す数値構成である。

---

## 10. 主要出力

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_generalized_elastic_collision_velocity_sweep_report_v1.md](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_report_v1.md) |
| 結果 JSON | [abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1.json](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1.json) |
| ケース CSV | [abc_multigauge_generalized_elastic_collision_velocity_sweep_cases_v1.csv](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_cases_v1.csv) |
| ステージ量 CSV | [abc_multigauge_generalized_elastic_collision_velocity_sweep_stage_quantities_v1.csv](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_stage_quantities_v1.csv) |
| ゲージ読出し CSV | [abc_multigauge_generalized_elastic_collision_velocity_sweep_gauge_rows_v1.csv](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_gauge_rows_v1.csv) |
| 誤差図 | [abc_multigauge_generalized_elastic_collision_velocity_sweep_errors_v1.png](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_errors_v1.png) |
| 位相勾配図 | [abc_multigauge_generalized_elastic_collision_velocity_sweep_q_outputs_v1.png](abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_q_outputs_v1.png) |

