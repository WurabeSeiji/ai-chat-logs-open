# ABC完全弾性衝突における多ゲージ干渉読出し数値検証メモ v1

**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、「ABC完全弾性衝突における多ゲージ干渉読出し仕様書 v1」に基づき、既存の ABC 完全弾性衝突モデルから、質量的量、運動量的量、エネルギー的量に見える保存読出しを、多ゲージ干渉読出しとして再構成できるかを段階的に検査した。

ここでいう読出し量は、標準物理の質量、運動量、エネルギーそのものとの同一視ではない。

本検証で確認するのは、閉じた位相系の内部観測から、次の三種類の保存読出しが構成できるかである。

| 読出し | 検査内容 |
|---|---|
| `p_read` | 空間位相方向の複素干渉相比から読まれ、完全反射で符号反転するか |
| `E_read` | 時間位相方向の複素干渉相比から読まれ、衝突前後で保存されるか |
| `R_read` | 複数ゲージで安定に残る校正済み強度残差として読まれ、衝突前後で保存されるか |

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_v2.py
```

出力先は次である。

```text
abc_multigauge_interference_readout_result_v2/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_interference_readout_v2.py
```

---

## 3. 段階構成

本検証は、次の三段階に分けた。

### 3.1 段階1: ABC基準衝突の再現

既存の ABC 完全弾性衝突条件を用い、A と B が相互作用セルへ到達し、進行方向読出し量が反転するかを確認した。

```text
q_A: +1 -> -1
q_B: -1 -> +1
```

同時に、識別振動 `m_A=1`, `m_B=2` が保存され、補償付き二乗閉鎖残差が閾値内に保たれるかを確認した。

### 3.2 段階2: 多ゲージ干渉読出し

観測機 `C` の参照波を複数ゲージに変え、各ゲージで複素相関を読む。

運動量的読出しは、

```text
p_read = arg(O(chi+h) / O(chi-h)) / (2h)
```

として読む。

エネルギー的読出しは、

```text
E_read = arg(O(tau+h) / O(tau-h)) / (2h)
```

として読む。

単一ゲージの値は測定成立とはしない。複数ゲージで同じ保存関係が再構成されることを要求した。

### 3.3 段階3: R軸安定残差と t/R 分離

質量的読出し候補 `R_read` は、単一の内部値ではなく、ゲージ変更後も安定に残る校正済み強度残差として読む。

また、`t` と `R` は第一原理的な固有軸ではない。

本検証では、変動の大きい連続読出しを `t` 的成分、ゲージを変えても安定して残る微小分散成分を `R` 的成分として分類し、次を測定した。

```text
separation_ratio_time = Var(R_read) / Var(t_read)
```

---

## 4. 実行結果

主要判定は次である。

| 項目 | 結果 |
|---|---:|
| `baseline_collision_valid` | `true` |
| `label_modes_preserved` | `true` |
| `closure_preserved` | `true` |
| `p_reconstructed_all_gauges` | `true` |
| `E_reconstructed_all_gauges` | `true` |
| `R_reconstructed_all_gauges` | `true` |
| `p_reflection_valid` | `true` |
| `E_preserved` | `true` |
| `R_preserved` | `true` |
| `R_gauge_stable` | `true` |
| `t_R_separation_valid` | `true` |
| `single_gauge_only_used` | `false` |
| `multigauge_measurement_valid` | `true` |

主要数値は次である。

| 量 | 値 |
|---|---:|
| `p_max_abs_error` | `2.5202062658991053e-14` |
| `E_max_abs_error` | `2.2315482794965646e-14` |
| `R_max_abs_error` | `4.4408920985006262e-16` |
| `R_max_gauge_std` | `2.2204460492503131e-16` |
| `closure_residual_abs` | `0.0000000000000000e+00` |
| `separation_ratio_time` | `3.6838616474030686e-30` |

---

## 5. 読出し結果の概要

A の `p_read` は、衝突前に `+1` と読まれ、衝突写像後に `-1` と読まれた。

B の `p_read` は、衝突前に `-1` と読まれ、衝突写像後に `+1` と読まれた。

一方、`E_read` は A/B ともに衝突前後で `1` と読まれ、保存された。

`R_read` は A/B ともに複数ゲージで `1` と読まれ、ゲージ分散は `2.2204460492503131e-16` に収まった。

したがって、本 v1 実験では、既存の ABC 完全弾性衝突モデルに対し、運動量的読出しの反転、エネルギー的読出しの保存、質量的読出し候補である R 軸安定残差の保存が、多ゲージ干渉読出しとして同時に再構成された。

---

## 6. R_gain 検査

R軸残差は微小であるため、読出し検証用に `R_gain` を変えた。

| `R_gain` | `R_mean` | `R_std` |
|---:|---:|---:|
| `1` | `1.0000000000000000e+00` | `2.5639502485114184e-16` |
| `10` | `1.0000000000000000e+01` | `2.0511601988091347e-15` |
| `100` | `1.0000000000000000e+02` | `2.4613922385709617e-14` |

このゲインは読出し検証用であり、物理的質量の増大とは主張しない。

---

## 7. 本検証で主張しないこと

本検証は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準物理の質量・運動量・エネルギーの完全導出 | 対応写像は別途必要である |
| `R_read` が実在時空の質量そのものであること | 本検証では質量的読出し候補である |
| 単一ゲージ値で測定が成立すること | 本検証では複数ゲージ再構成を要求した |
| 外部計量符号の導入 | 本検証は全正符号ゼロ閉鎖の読出し表示として扱う |

---

## 8. 結論

本 v1 実験では、既存の ABC 完全弾性衝突モデルに対し、複数ゲージの干渉相関読出しを追加した。

結果として、次が同時に成立した。

1. 基準 ABC 完全弾性衝突が再現された。
2. 識別振動 `m_A,m_B` が保存された。
3. 補償付き二乗閉鎖残差が `0` に保たれた。
4. `p_read` が複数ゲージで再構成され、完全反射で符号反転した。
5. `E_read` が複数ゲージで再構成され、衝突前後で保存された。
6. `R_read` が複数ゲージで安定に再構成され、衝突前後で保存された。
7. `separation_ratio_time=3.6838616474030686e-30` により、t/R 分離が確認された。

したがって、本検証の範囲では、質量的量、運動量的量、エネルギー的量に見える保存読出しは、単一内部値ではなく、多ゲージ干渉相関から再構成できた。

---

# 付録A. 出力ファイル

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_interference_readout_report_v2.md](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_report_v2.md) |
| 結果 JSON | [abc_multigauge_interference_readout_result_v2.json](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_result_v2.json) |
| ステージ CSV | [abc_multigauge_interference_readout_timeline_v2.csv](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_timeline_v2.csv) |
| ゲージ CSV | [abc_multigauge_interference_readout_gauge_sweep_v2.csv](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_gauge_sweep_v2.csv) |
| ステージ要約 CSV | [abc_multigauge_interference_readout_stage_summary_v2.csv](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_stage_summary_v2.csv) |
| イベント CSV | [abc_multigauge_interference_readout_events_v2.csv](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_events_v2.csv) |
| R_gain CSV | [abc_multigauge_interference_readout_r_gain_sweep_v2.csv](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_r_gain_sweep_v2.csv) |
| p/E/R 図 | [abc_multigauge_interference_readout_invariants_v2.png](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_invariants_v2.png) |
| t/R 分離図 | [abc_multigauge_interference_readout_tr_separation_v2.png](abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_tr_separation_v2.png) |

