# ABC多ゲージ干渉読出しの非対称振幅診断スイープ検証メモ v1

**副題:** 等振幅完全反射写像の成立範囲と、質量的重み付き運動量保存の破れを分離する数値診断  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、ABC 完全弾性衝突モデルで構成した多ゲージ干渉読出しを、非対称振幅条件へ拡張して診断する。

前段階では、等振幅条件

```text
A_A = A_B = 1
```

において、複数ゲージから

```text
p_read
E_read
R_read
```

を再構成し、完全反射写像と整合することを確認した。

本実験では、`A_A` と `A_B` を意図的に変える。

目的は、非対称条件でも全てが成立すると主張することではない。

目的は、次を分離して測ることである。

| 対象 | 判定 |
|---|---|
| 個体別の `p/E/R` 読出し | 非対称振幅でも成立するか |
| `R_read` の保存 | 各個体および全体で保存されるか |
| `R_read * p_read` の全体系保存 | 単純な `q` 反転写像と整合するか |

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1.py
```

出力先は次である。

```text
abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1.py
```

---

## 3. 実験条件

振幅ケースは次である。

| `A_A` | `A_B` | 意味 |
|---:|---:|---|
| `1.0` | `1.0` | 等振幅基準 |
| `1.0` | `1.10` | B がわずかに大きい |
| `1.0` | `1.25` | B が中程度に大きい |
| `1.0` | `1.50` | B が大きい |
| `1.0` | `2.00` | B がかなり大きい |
| `1.0` | `3.00` | B が非常に大きい |
| `1.50` | `1.0` | A が大きい |
| `2.00` | `1.0` | A がかなり大きい |

質量的読出し候補は、

```text
R_read ~= A^2
```

として再構成される。

したがって、非対称振幅では、

```text
R_A != R_B
```

となる。

---

## 4. 判定上の注意

非対称振幅では、A と B の `R_read` が本当に異なる。

したがって、等振幅実験で使ったような「全粒子をまとめた `R` 分散」を、そのまま不安定性判定に使ってはいけない。

全体の `R` 分散が大きくなることは、非対称振幅では失敗ではなく、むしろ質量的差を読めていることの診断信号である。

そのため、本実験では `t/R` 分離を、次のように個体ごとに評価する。

```text
within_particle_separation_ratio_time
  = max over P in {A,B} Var_s(R_read(P,s)) / Var_s(t_read(P,s))
```

ここで見るのは、A と B の間の `R` 差ではなく、同じ個体の `R_read` が時間発展とゲージ変更を通じて安定しているかである。

---

## 5. 実行結果

全体判定は次である。

```text
case_count: 8
asymmetric_case_count: 7
individual_multigauge_valid_all_cases: true
weighted_energy_preserved_all_cases: true
R_total_preserved_all_cases: true
equal_case_weighted_momentum_preserved: true
asymmetric_cases_detect_weighted_momentum_failure: true
asymmetric_amplitude_diagnostic_valid: true
```

主要値は次である。

| 量 | 値 |
|---|---:|
| `max_p_abs_error` | `3.6193270602780103e-14` |
| `max_E_abs_error` | `2.6423307986078726e-14` |
| `max_R_abs_error` | `3.552713678800501e-15` |
| `max_R_gauge_std` | `2.082963028648268e-15` |
| `max_global_R_contrast_ratio_time` | `1195.4814536929925` |
| `max_within_particle_separation_ratio_time` | `3.683878850562613e-30` |
| `max_weighted_p_collision_error` | `16.000000000000036` |

---

## 6. ケース別結果

| case | `R_B/R_A` | 個体別読出し | 個体内 `Var(R)/Var(t)` | `R*p` 全運動量誤差 | `R*E` 誤差 | `R_total` 誤差 | 単純反転と整合 |
|---|---:|---|---:|---:|---:|---:|---|
| `A_1.00_B_1.00` | `1.0000000000000000e+00` | `true` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `true` |
| `A_1.00_B_1.10` | `1.2100000000000002e+00` | `true` | `0.0000000000000000e+00` | `4.2000000000000137e-01` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |
| `A_1.00_B_1.25` | `1.5625000000000000e+00` | `true` | `3.6838788505626130e-30` | `1.1249999999999964e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |
| `A_1.00_B_1.50` | `2.2500000000000000e+00` | `true` | `0.0000000000000000e+00` | `2.5000000000000124e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |
| `A_1.00_B_2.00` | `4.0000000000000000e+00` | `true` | `0.0000000000000000e+00` | `5.9999999999999858e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |
| `A_1.00_B_3.00` | `9.0000000000000000e+00` | `true` | `0.0000000000000000e+00` | `1.6000000000000036e+01` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |
| `A_1.50_B_1.00` | `4.4444444444444442e-01` | `true` | `0.0000000000000000e+00` | `2.5000000000000124e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |
| `A_2.00_B_1.00` | `2.5000000000000000e-01` | `true` | `0.0000000000000000e+00` | `5.9999999999999858e+00` | `0.0000000000000000e+00` | `0.0000000000000000e+00` | `false` |

---

## 7. 解釈

本実験では、非対称振幅でも、個体別の多ゲージ干渉読出しは成立した。

すなわち、

```text
p_read
E_read
R_read
```

は全ケースで再構成された。

また、`R_total` と `R*E` は全ケースで保存された。

一方、`R*p` による全運動量的読出しは、等振幅ケースでのみ保存され、非対称振幅ケースでは保存されなかった。

これは失敗ではなく、重要な診断である。

単純な `q -> -q` 反転写像は、等振幅・等質量的条件では完全弾性反射として読める。しかし、質量的読出し `R_read` が非対称になると、同じ写像は `R*p` の全体系保存とは整合しない。

したがって、本実験は次を分離した。

| 分離対象 | 結果 |
|---|---|
| 多ゲージ干渉読出しの成立 | 非対称振幅でも成立 |
| `R_read` の個体別保存 | 成立 |
| `R_total` の保存 | 成立 |
| `R*E` の保存 | 成立 |
| 単純 `q` 反転写像と `R*p` 全運動量保存の整合 | 等振幅でのみ成立 |

この結果は、非対称質量的条件では、読出しの問題ではなく、衝突写像そのものを拡張する必要があることを示す。

---

## 8. 本実験で主張しないこと

本実験は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 非対称質量衝突の完成写像 | 今回は単純 `q` 反転写像の成立範囲を診断した |
| 標準力学の衝突公式の導出 | 対応写像は別途必要である |
| `R_read` が標準質量そのものであること | 本稿では質量的読出し候補である |
| 非対称条件で完全弾性衝突が不可能であること | 必要なのは写像の拡張であり、不可能性ではない |

本実験の主張は、ABC 多ゲージ干渉読出しにおいて、非対称振幅でも個体別 `p/E/R` は読めるが、等振幅用の単純反転写像は `R*p` 全運動量保存と整合しない、という診断結果である。

---

## 9. 主要出力

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_interference_readout_asymmetric_amplitude_report_v1.md](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_report_v1.md) |
| 結果 JSON | [abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1.json](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1.json) |
| ケース CSV | [abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.csv](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.csv) |
| ステージ量 CSV | [abc_multigauge_interference_readout_asymmetric_amplitude_stage_quantities_v1.csv](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_stage_quantities_v1.csv) |
| ゲージ読出し CSV | [abc_multigauge_interference_readout_asymmetric_amplitude_gauge_rows_v1.csv](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_gauge_rows_v1.csv) |
| 保存診断図 | [abc_multigauge_interference_readout_asymmetric_amplitude_conservation_v1.png](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_conservation_v1.png) |
| ケース図 | [abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.png](abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.png) |

