# ABC多ゲージ干渉読出しによる一般化弾性衝突の読出しノイズ頑健性検証メモ v1

**副題:** ゼロ平均ゲージ揺らぎの相殺と共通読出しバイアスの検出可能性  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、一般化弾性衝突写像の多ゲージ干渉読出しに対して、読出し器側のノイズを加えたときの頑健性を検証した。

ここでノイズは、物理状態そのものへ入れるのではない。状態をシミュレーションした後、各ゲージで得た

```text
p_read
E_read
R_read
```

の読出し行へ、決定論的な擬似ノイズとして加えた。

本検証の問いは次である。

> 単一ゲージ値が揺らいでも、ゼロ平均のゲージ間揺らぎであれば、多ゲージ干渉平均によって `p,E,R` と `R*p`, `R*p^2` の保存読出しは回復するか。また、全ゲージに共通するバイアスは検出できるか。

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_noise_robustness_v2.py
```

出力先は次である。

```text
abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_generalized_elastic_collision_noise_robustness_v2.py
```

---

## 3. 実験条件

検証ケースは次である。

| case | `A_A` | `A_B` | `p_A` before | `p_B` before | 備考 |
|---|---:|---:|---:|---:|---|
| c01 | `1.0` | `1.0` | `1.0` | `-1.0` | 等振幅・基準対向 |
| c02 | `1.0` | `2.0` | `1.4` | `-0.6` | 非対称R・非対称対向 |
| c03 | `1.0` | `2.0` | `1.2` | `0.2` | 非対称R・同方向追いつき |
| c04 | `1.5` | `1.0` | `1.8` | `-0.2` | A側が大きい非対称R |

各ケースで、116個の読出しゲージを用いた。

ノイズレベルは次である。

```text
0
1e-12
1e-10
1e-8
1e-6
1e-4
```

ノイズモードは次の2種類である。

| モード | 意味 | 期待 |
|---|---|---|
| `zero_mean_gauge_noise` | 各 stage/particle ごとにゲージ間平均がゼロになる揺らぎ | 多ゲージ平均で相殺される |
| `common_bias_control` | 全ゲージに共通して入るバイアス | 多ゲージ平均では消えず検出される |

---

## 4. 判定量

### 4.1 ゼロ平均ノイズの相殺

ゼロ平均ゲージ揺らぎに対して、

```text
p_mean_abs_error
E_mean_abs_error
R_mean_abs_error
```

が小さく保たれるかを確認する。

さらに、一般化弾性衝突読出しの保存量として、

```text
R*p
R*p^2
relative p flip
```

が維持されるかを測定する。

### 4.2 共通バイアスの検出

全ゲージに共通するバイアスは、多ゲージ平均では相殺されない。

本検証では検出床を

```text
1e-10
```

とし、それ以上の共通バイアスが全ケースで検出されるかを確認する。

---

## 5. 実行結果

全体判定は次である。

```text
case_count: 4
noise_mode_count: 2
noise_level_count: 6
total_summary_rows: 48
max_gauge_count: 116
zero_mean_multigauge_valid_all: true
common_bias_detection_floor: 1e-10
common_bias_detected_all_above_floor: true
single_gauge_only_used: false
noise_robustness_valid: true
```

ゼロ平均ノイズに対する主要誤差は次である。

| 量 | 値 |
|---|---:|
| `zero_mean_max_p_mean_abs_error` | `1.4477308241112041e-13` |
| `zero_mean_max_E_mean_abs_error` | `6.439293542825908e-15` |
| `zero_mean_max_R_mean_abs_error` | `3.552713678800501e-14` |
| `zero_mean_max_P_R_error` | `1.674216321134736e-13` |
| `zero_mean_max_K_R_error` | `8.322231792590173e-13` |
| `zero_mean_max_relative_flip_error` | `1.687538997430238e-14` |

共通バイアス制御の最大平均誤差は次である。

| 量 | 値 |
|---|---:|
| `biased_control_max_p_mean_abs_error` | `2.1624903162420495e-04` |
| `biased_control_max_E_mean_abs_error` | `9.971756487914263e-05` |
| `biased_control_max_R_mean_abs_error` | `3.9675791940574356e-04` |

---

## 6. ノイズレベル別結果

ゼロ平均ゲージ揺らぎでは、最大 `1e-4` の読出し側ノイズを入れても、多ゲージ平均後の主要誤差は機械精度近傍に保たれた。

| noise level | max `p` mean err | max `R` mean err | max `R*p^2` err | valid |
|---:|---:|---:|---:|---|
| `0.0e+00` | `1.4432899320127035e-13` | `3.5527136788005009e-14` | `8.2867046558021684e-13` | `true` |
| `1.0e-12` | `1.4477308241112041e-13` | `3.5527136788005009e-14` | `8.2867046558021684e-13` | `true` |
| `1.0e-10` | `1.4477308241112041e-13` | `3.5527136788005009e-14` | `8.3044682241961709e-13` | `true` |
| `1.0e-08` | `1.4477308241112041e-13` | `3.5527136788005009e-14` | `8.3133500083931722e-13` | `true` |
| `1.0e-06` | `1.4477308241112041e-13` | `3.5527136788005009e-14` | `8.3222317925901734e-13` | `true` |
| `1.0e-04` | `1.4432899320127035e-13` | `3.5527136788005009e-14` | `8.2778228716051672e-13` | `true` |

共通バイアス制御では、`1e-10` 以上のバイアスが全ケースで検出された。

| noise level | max `p` mean err | max `R` mean err | detected |
|---:|---:|---:|---|
| `1.0e-12` | `2.1973534103381098e-12` | `3.7871927816013340e-12` | `false` |
| `1.0e-10` | `1.8466961293484019e-10` | `3.8134029267666847e-10` | `true` |
| `1.0e-08` | `2.3795240000623608e-08` | `3.9659822093085495e-08` | `true` |
| `1.0e-06` | `1.7035717845281795e-06` | `3.9197052124073650e-06` | `true` |
| `1.0e-04` | `2.1624903162420495e-04` | `3.9675791940574356e-04` | `true` |

`1e-12` は本検証で設定した検出床 `1e-10` 未満であるため、共通バイアス検出の主張対象から外す。

---

## 7. 解釈

本検証では、単一ゲージ値が揺らぐこと自体は許容した。

重要なのは、その揺らぎがゲージ間でゼロ平均であれば、多ゲージ干渉平均によって

```text
p_read
E_read
R_read
R*p
R*p^2
```

の読出しが回復することである。

一方、全ゲージに共通して入るバイアスは、平均では消えない。したがって、共通バイアスは保存読出しの破れとして検出される。

この結果は、多ゲージ干渉読出しが単に多数の値を平均しているだけではなく、ゲージ間で相殺可能な揺らぎと、相殺できない共通バイアスを分離する読出し器として働くことを示す。

---

## 8. 本検証で主張しないこと

本検証は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 任意ノイズ下で測定が成立すること | 本検証はゼロ平均ゲージ揺らぎと共通バイアス制御に限る |
| 単一ゲージ測定の成立 | 判定は多ゲージ干渉平均に基づく |
| 実験装置ノイズの定量モデル | ノイズは読出し側の決定論的擬似ノイズである |
| 標準物理量との直接対応 | `p,E,R` は位相系内の読出し量である |

---

## 9. 結論

本検証では、一般化弾性衝突写像の多ゲージ干渉読出しに、読出し器側ノイズを加えた。

ゼロ平均ゲージ揺らぎでは、最大 `1e-4` のノイズを加えても、116ゲージの多ゲージ平均によって主要読出し誤差は機械精度近傍に戻った。

最大 `p` 平均誤差は `1.4477308241112041e-13`、最大 `R` 平均誤差は `3.552713678800501e-14`、最大 `R*p^2` 保存誤差は `8.322231792590173e-13` であった。

共通バイアス制御では、検出床 `1e-10` 以上のバイアスが全ケースで検出された。

したがって、本検証の範囲では、一般化弾性衝突における `p,E,R` 読出しは、ゼロ平均のゲージ間揺らぎに対して頑健であり、同時に、全ゲージ共通の読出しバイアスを検出可能であった。

---

# 付録A. 出力ファイル

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_generalized_elastic_collision_noise_robustness_report_v2.md](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_report_v2.md) |
| 結果 JSON | [abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2.json](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2.json) |
| summary CSV | [abc_multigauge_generalized_elastic_collision_noise_robustness_summary_v2.csv](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_summary_v2.csv) |
| stage quantity CSV | [abc_multigauge_generalized_elastic_collision_noise_robustness_stage_quantities_v2.csv](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_stage_quantities_v2.csv) |
| gauge CSV | [abc_multigauge_generalized_elastic_collision_noise_robustness_gauge_rows_v2.csv](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_gauge_rows_v2.csv) |
| zero mean 図 | [abc_multigauge_generalized_elastic_collision_noise_robustness_zero_mean_v2.png](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_zero_mean_v2.png) |
| bias control 図 | [abc_multigauge_generalized_elastic_collision_noise_robustness_bias_control_v2.png](abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_bias_control_v2.png) |
