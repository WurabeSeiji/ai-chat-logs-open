# ABC多ゲージ干渉読出しによる一般化弾性衝突の複数回衝突検証メモ v1

**副題:** `R*p`・`R*p^2` 保存写像の反復衝突安定性と壁戻し条件下での多ゲージ読出し検証  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、前段階で構成した `R_read` 重み付き一般化弾性衝突写像を、複数回の AB 衝突へ反復適用した。

単回検証では、非対称振幅および非対称初期位相勾配に対して、

```text
R_A p_A + R_B p_B
R_A p_A^2 + R_B p_B^2
```

が保存され、相対位相勾配

```text
p_A - p_B
```

が符号反転することを確認した。

本検証では、同じ一般化写像を複数回の AB 衝突へ適用し、壁反射による再遭遇を含む反復過程でも、多ゲージ干渉読出しが安定に成立するかを調べた。

ここで `p` は標準物理の速度または運動量そのものではない。空間位相方向の多ゲージ干渉読出し勾配である。また `R` は固定された背景座標軸ではなく、読出し振幅安定軸として再構成される量である。

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_multi_collision_v2.py
```

出力先は次である。

```text
abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_generalized_elastic_collision_multi_collision_v2.py
```

---

## 3. 実験条件

反復検証ケースは次である。

| case | `A_A` | `A_B` | `p_A` before | `p_B` before | 備考 |
|---|---:|---:|---:|---:|---|
| c01 | `1.0` | `1.0` | `1.0` | `-1.0` | 等振幅・基準対向 |
| c02 | `1.0` | `2.0` | `1.4` | `-0.6` | 非対称R・非対称対向 |
| c03 | `1.0` | `2.0` | `1.2` | `0.2` | 非対称R・同方向追いつき |
| c04 | `1.5` | `1.0` | `1.8` | `-0.2` | A側が大きい非対称R |

各ケースで AB 衝突目標回数を

```text
6
```

とした。

壁反射は粒子を再遭遇させるための補助条件であり、保存判定の対象ではない。保存判定は、各 AB 衝突の直前と直後だけに限定した。

---

## 4. 判定量

本検証の判定量は次である。

### 4.1 個別読出し

複数ゲージから各粒子の

```text
p_read
E_read
R_read
```

を再構成する。

単一ゲージ値だけでは読出し成立とは判定しない。

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

### 4.5 `R*E_tau` 保存

```text
E_{\tau,R} = R_A E_A + R_B E_B
```

### 4.6 `R` 安定性

各 AB 衝突の前後で、

```text
R_A
R_B
```

が保存されるかを確認する。

---

## 5. 実行結果

全体判定は次である。

```text
case_count: 4
all_cases_valid: true
completed_target_collisions_all_cases: true
individual_readout_valid_all_cases: true
closure_preserved_all_cases: true
P_R_preserved_each_collision_all_cases: true
K_R_preserved_each_collision_all_cases: true
relative_flip_each_collision_all_cases: true
E_tau_R_preserved_each_collision_all_cases: true
R_preserved_each_collision_all_cases: true
single_gauge_only_used: false
generalized_multi_collision_valid: true
```

主要誤差は次である。

| 量 | 値 |
|---|---:|
| `max_ab_collision_count` | `6` |
| `max_wall_reflection_count` | `8` |
| `max_p_abs_error` | `1.0336176359260207e-13` |
| `max_E_abs_error` | `3.341771304121721e-13` |
| `max_R_abs_error` | `2.255973186038318e-13` |
| `max_R_gauge_std` | `9.768953270904304e-14` |
| `max_within_particle_separation_ratio_time` | `1.4212404920023786e-27` |
| `max_P_R_error` | `3.552713678800501e-14` |
| `max_K_R_error` | `1.056932319443149e-13` |
| `max_relative_flip_error` | `2.220446049250313e-14` |
| `max_E_tau_R_error` | `0.0` |
| `max_R_A_error` | `4.440892098500626e-16` |
| `max_R_B_error` | `8.881784197001252e-16` |

---

## 6. ケース別結果

| case | AB collisions | wall reflections | max `R*p` err | max `R*p^2` err | max relative err | valid |
|---|---:|---:|---:|---:|---:|---|
| `c01_A1.00_B1.00_u1.00_v-1.00` | `6` | `5` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `true` |
| `c02_A1.00_B2.00_u1.40_v-0.60` | `6` | `6` | `3.5527136788005009e-14` | `7.9936057773011271e-14` | `1.6875389974302379e-14` | `true` |
| `c03_A1.00_B2.00_u1.20_v0.20` | `6` | `8` | `2.4646951146678475e-14` | `3.5527136788005009e-14` | `1.7097434579227411e-14` | `true` |
| `c04_A1.50_B1.00_u1.80_v-0.20` | `6` | `6` | `3.3750779948604759e-14` | `1.0569323194431490e-13` | `2.2204460492503131e-14` | `true` |

---

## 7. 解釈

本検証では、単回衝突で構成した一般化弾性衝突写像を、4種類の初期条件で6回ずつ反復した。

等振幅の基準対向だけでなく、非対称R、非対称位相勾配、同方向追いつき条件を含めても、各 AB 衝突の直前直後で

```text
R*p
R*p^2
R*E_tau
R_A
R_B
```

が保存された。

また、相対位相勾配

```text
p_A - p_B
```

は各 AB 衝突で符号反転した。

これは、一般化弾性衝突写像が単発の代入規則としてだけでなく、反復過程でも多ゲージ干渉読出し量を保存することを示す。

---

## 8. 本検証で主張しないこと

本検証は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準運動量・標準エネルギーを直接導出したこと | 本稿の `p,E,R` は位相系内の読出し量である |
| 壁反射を実在境界条件として一般化すること | 壁反射は複数回 AB 衝突を生成する補助条件である |
| 単一ゲージ測定の成立 | 判定には複数ゲージ干渉読出しを用いる |
| 任意多体系への拡張 | 本検証は ABC セル内の二局所波 AB 衝突に限定する |

---

## 9. 結論

本検証では、`R_read` 重み付き一般化弾性衝突写像を、複数回の AB 衝突へ反復適用した。

4ケースすべてで目標6回の AB 衝突が完了し、全ケースで

```text
P_R_preserved_each_collision
K_R_preserved_each_collision
relative_flip_each_collision
E_tau_R_preserved_each_collision
R_preserved_each_collision
```

が `true` となった。

最大 `R*p` 保存誤差は `3.552713678800501e-14`、最大 `R*p^2` 保存誤差は `1.056932319443149e-13`、最大相対位相勾配反転誤差は `2.220446049250313e-14` であった。

したがって、本検証の範囲では、非対称Rおよび非対称位相勾配を含む一般化弾性衝突写像は、複数回の AB 衝突反復においても、多ゲージ干渉読出し上の保存関係を維持した。

---

# 付録A. 出力ファイル

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_generalized_elastic_collision_multi_collision_report_v2.md](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_report_v2.md) |
| 結果 JSON | [abc_multigauge_generalized_elastic_collision_multi_collision_result_v2.json](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_result_v2.json) |
| ケース CSV | [abc_multigauge_generalized_elastic_collision_multi_collision_cases_v2.csv](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_cases_v2.csv) |
| 衝突読出し CSV | [abc_multigauge_generalized_elastic_collision_multi_collision_readouts_v2.csv](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_readouts_v2.csv) |
| ゲージ CSV | [abc_multigauge_generalized_elastic_collision_multi_collision_gauge_rows_v2.csv](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_gauge_rows_v2.csv) |
| `p` 誤差図 | [abc_multigauge_generalized_elastic_collision_multi_collision_p_errors_v2.png](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_p_errors_v2.png) |
| 要約図 | [abc_multigauge_generalized_elastic_collision_multi_collision_summary_v2.png](abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_summary_v2.png) |
