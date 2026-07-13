# ABC多ゲージ干渉読出しによる一般化弾性衝突写像検証メモ v1

**副題:** 質量的読出し `R_read` を重みとする `R*p` および `R*p^2` 保存写像の数値構成  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、ABC 完全弾性衝突モデルにおいて、非対称振幅条件でも保存する一般化衝突写像を構成できるかを検査する。

前段階の非対称振幅診断では、次が分かった。

```text
個体別 p/E/R 読出しは非対称振幅でも成立する。
R_total と R*E は保存される。
しかし単純な q -> -q 反転写像では、R*p 全運動量的読出しが保存しない。
```

したがって、本実験では、単純な `q` 符号反転を用いない。

多ゲージ干渉読出しで得た質量的量

```text
R_A, R_B
```

を重みとして用い、衝突前後で

```text
R_A p_A + R_B p_B
```

および

```text
R_A p_A^2 + R_B p_B^2
```

が保存されるように、衝突後の位相勾配 `p_A', p_B'` を構成する。

ここで `p` は標準物理の運動量そのものではなく、空間位相方向の干渉読出し勾配である。

---

## 2. 本実験で主張しないこと

本実験は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準力学の完全導出 | 対応写像は別途必要である |
| `R_read` が標準質量そのものであること | 本稿では質量的読出し候補として扱う |
| `R*p` が標準運動量そのものであること | 本稿では質量的重み付き位相勾配読出しである |
| `R*p^2` が標準エネルギーそのものであること | 本稿では位相勾配二乗の保存読出しである |
| 外部質量パラメータの導入 | `R` は多ゲージ干渉読出しから得る |

本実験の主張は、内部読出しで得た `R` を重みとして、非対称振幅でも保存する衝突写像を数値構成できるかを調べることである。

---

## 3. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_readout_v2.py
```

出力先は次である。

```text
abc_multigauge_generalized_elastic_collision_readout_result_v2/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_generalized_elastic_collision_readout_v2.py
```

---

## 4. 一般化衝突写像

質量的読出しを、

```text
R_A = A_A^2
R_B = A_B^2
```

として多ゲージ干渉読出しから再構成する。

衝突前の位相勾配を、

```text
u_A = p_A
u_B = p_B
```

とする。

一般化写像では、衝突後の位相勾配を次で与える。

```text
v_A =
  ((R_A - R_B) / (R_A + R_B)) u_A
  + (2 R_B / (R_A + R_B)) u_B
```

```text
v_B =
  (2 R_A / (R_A + R_B)) u_A
  + ((R_B - R_A) / (R_A + R_B)) u_B
```

この写像は、等振幅条件

```text
R_A = R_B
```

では、

```text
v_A = u_B
v_B = u_A
```

となる。

初期条件が

```text
u_A = +1
u_B = -1
```

であれば、これは前段階の単純反転

```text
q_A -> -q_A
q_B -> -q_B
```

と一致する。

---

## 5. 実験条件

振幅ケースは次である。

| `A_A` | `A_B` |
|---:|---:|
| `1.0` | `1.0` |
| `1.0` | `1.10` |
| `1.0` | `1.25` |
| `1.0` | `1.50` |
| `1.0` | `2.00` |
| `1.0` | `3.00` |
| `1.50` | `1.0` |
| `2.00` | `1.0` |

各ケースで、多ゲージ読出しにより、

```text
p_read
E_read
R_read
```

を再構成する。

そのうえで、衝突前後の保存量を評価する。

---

## 6. 判定量

判定量は次である。

### 6.1 質量的重み付き位相勾配

```text
P_R = R_A p_A + R_B p_B
```

### 6.2 質量的重み付き位相勾配二乗

```text
K_R = R_A p_A^2 + R_B p_B^2
```

### 6.3 時間位相読出しの重み付き量

```text
E_tau_R = R_A E_A + R_B E_B
```

### 6.4 質量的総量

```text
R_total = R_A + R_B
```

---

## 7. 実行結果

全体判定は次である。

```text
case_count: 8
individual_readout_valid_all_cases: true
generalized_P_R_preserved_all_cases: true
generalized_K_R_phase_preserved_all_cases: true
E_tau_R_preserved_all_cases: true
R_total_preserved_all_cases: true
single_gauge_only_used: false
generalized_elastic_collision_readout_valid: true
```

主要誤差は次である。

| 量 | 値 |
|---|---:|
| `max_p_abs_error` | `2.8688162956314045e-13` |
| `max_E_abs_error` | `2.6423307986078726e-14` |
| `max_R_abs_error` | `1.865174681370263e-14` |
| `max_R_gauge_std` | `8.078732199757252e-15` |
| `max_within_particle_separation_ratio_time` | `2.365357211965411e-28` |
| `max_P_R_conservation_error` | `2.3803181647963356e-13` |
| `max_K_R_phase_conservation_error` | `1.4086509736443986e-12` |
| `max_E_tau_R_conservation_error` | `3.552713678800501e-15` |
| `max_R_total_conservation_error` | `5.329070518200751e-15` |

---

## 8. ケース別結果

| case | `R_B/R_A` | `q_A` after | `q_B` after | `R*p` err | `R*p^2` err | simple flip `R*p` err |
|---|---:|---:|---:|---:|---:|---:|
| `A_1.00_B_1.00` | `1.0000000000000000e+00` | `-1.0000000000000000e+00` | `1.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` |
| `A_1.00_B_1.10` | `1.2100000000000002e+00` | `-1.1900452488687785e+00` | `8.0995475113122162e-01` | `1.2212453270876722e-14` | `1.4654943925052066e-14` | `4.2000000000000037e-01` |
| `A_1.00_B_1.25` | `1.5625000000000000e+00` | `-1.4390243902439024e+00` | `5.6097560975609762e-01` | `6.7723604502134549e-15` | `1.8651746813702630e-14` | `1.1250000000000000e+00` |
| `A_1.00_B_1.50` | `2.2500000000000000e+00` | `-1.7692307692307692e+00` | `2.3076923076923078e-01` | `5.5511151231257827e-15` | `1.3322676295501878e-15` | `2.5000000000000000e+00` |
| `A_1.00_B_2.00` | `4.0000000000000000e+00` | `-2.2000000000000002e+00` | `-1.9999999999999996e-01` | `1.2878587085651816e-14` | `2.3980817331903381e-14` | `6.0000000000000000e+00` |
| `A_1.00_B_3.00` | `9.0000000000000000e+00` | `-2.6000000000000001e+00` | `-6.0000000000000009e-01` | `2.3803181647963356e-13` | `1.4086509736443986e-12` | `1.6000000000000000e+01` |
| `A_1.50_B_1.00` | `4.4444444444444442e-01` | `-2.3076923076923078e-01` | `1.7692307692307692e+00` | `1.0436096431476471e-14` | `1.4654943925052066e-14` | `2.5000000000000000e+00` |
| `A_2.00_B_1.00` | `2.5000000000000000e-01` | `1.9999999999999996e-01` | `2.2000000000000002e+00` | `9.7699626167013776e-15` | `6.2172489379008766e-15` | `6.0000000000000000e+00` |

---

## 9. 解釈

本実験では、非対称振幅条件でも、個体別の多ゲージ干渉読出しは成立した。

さらに、質量的読出し `R_read` を重みとして用いた一般化衝突写像により、

```text
R*p
R*p^2
R*E_tau
R_total
```

が全ケースで保存された。

一方、単純な `q -> -q` 反転写像は、等振幅ケースでは一般化写像と一致するが、非対称振幅ケースでは `R*p` 保存と整合しない。

したがって、前段階で観測された非対称振幅での保存破れは、読出し不能によるものではなく、等振幅専用写像を非対称条件へ適用したことによる写像不足であった。

本実験により、次の整理が得られた。

| 対象 | 結果 |
|---|---|
| 個体別 `p/E/R` 読出し | 非対称条件でも成立 |
| 単純 `q` 反転写像 | 等振幅の特殊例 |
| 非対称条件での保存写像 | `R` 重み付き一般化写像で成立 |
| `R*p` 保存 | 全ケースで成立 |
| `R*p^2` 保存 | 全ケースで成立 |

---

## 10. 本実験の位置づけ

本実験は、標準力学の質量、運動量、エネルギーをそのまま導出したものではない。

しかし、閉じた位相系内で、

```text
質量的量としての R_read
運動量的量としての p_read
位相勾配二乗型の保存量
```

を多ゲージ干渉読出しから構成し、それらを用いた非対称衝突保存写像を実行できることを示した。

これは、前段階の ABC 完全弾性衝突モデルを、等振幅・等質量的条件から非対称質量的条件へ拡張するための基礎実験である。

---

## 11. 主要出力

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_generalized_elastic_collision_readout_report_v2.md](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_readout_report_v2.md) |
| 結果 JSON | [abc_multigauge_generalized_elastic_collision_readout_result_v2.json](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_readout_result_v2.json) |
| ケース CSV | [abc_multigauge_generalized_elastic_collision_cases_v2.csv](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_cases_v2.csv) |
| ステージ量 CSV | [abc_multigauge_generalized_elastic_collision_stage_quantities_v2.csv](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_stage_quantities_v2.csv) |
| ゲージ読出し CSV | [abc_multigauge_generalized_elastic_collision_gauge_rows_v2.csv](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_gauge_rows_v2.csv) |
| 保存診断図 | [abc_multigauge_generalized_elastic_collision_conservation_v2.png](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_conservation_v2.png) |
| 衝突後位相勾配図 | [abc_multigauge_generalized_elastic_collision_q_outputs_v2.png](abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_q_outputs_v2.png) |

