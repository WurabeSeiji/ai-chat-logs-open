# ABC多ゲージ干渉読出しの読出し器頑健性スイープ検証メモ v1

**副題:** 読出し中心・参照位相・読出し幅・観測ゲインを変えたときの質量的量・運動量的量・エネルギー的量の再構成安定性  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、ABC 完全弾性衝突における多ゲージ干渉読出しが、特定の読出し器設定にだけ依存していないかを検査する。

前段階では、単発衝突および複数回衝突において、複数ゲージから

```text
p_read
E_read
R_read
```

を再構成できることを確認した。

本実験では、衝突軌道そのものは固定し、読出し器側のゲージ群だけを変える。

検査する変更は次である。

| 変更対象 | 内容 |
|---|---|
| 読出し中心 | `delta_chi`, `delta_tau` を変える |
| 参照位相 | `delta_phi` を変える |
| 差分幅 | `h_chi`, `h_tau` を変える |
| 観測機幅 | `nh_chi_c`, `nh_tau_c` を変える |
| 観測ゲイン | `c_gain` を変える |

単一ゲージ値だけを測定値とはしない。各ケースは、複数ゲージ全体で同じ保存読出しを再構成できる場合にのみ成立とする。

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_robustness_sweep_v1.py
```

出力先は次である。

```text
abc_multigauge_interference_readout_robustness_sweep_result_v1/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_interference_readout_robustness_sweep_v1.py
```

---

## 3. ゲージケース

今回の頑健性スイープでは、次の5ケースを用いた。

| ケース | ゲージ数 | 目的 |
|---|---:|---|
| `baseline_default` | `8` | 既定ゲージ群の再確認 |
| `phase_center_grid` | `45` | 参照位相と読出し中心の同時変更 |
| `width_gain_grid` | `36` | 読出し幅、観測機幅、観測ゲインの変更 |
| `near_lobe_offset` | `25` | 局在波カーネルの中心からやや離れた読出し |
| `mixed_readout_grid` | `16` | 位相、中心、幅、ゲインの混合変更 |

総ゲージ数は `130` である。

---

## 4. 判定条件

各ケースで、次を要求した。

1. ABC 衝突セルで完全弾性反射写像が成立する。
2. 識別振動 `m_A,m_B` が保存される。
3. 補償付き二乗閉鎖が保存される。
4. 全ゲージで `p_read` が再構成される。
5. 全ゲージで `E_read` が再構成される。
6. 全ゲージで `R_read` が再構成される。
7. `p_read` が衝突写像で符号反転する。
8. `E_read` が衝突前後で保存される。
9. `R_read` が衝突前後で保存される。
10. `R_read` のゲージ分散が閾値内に収まる。
11. `Var(R)/Var(t)` が十分小さく、t/R 分離が成立する。
12. 単一ゲージだけを判定に使わない。

---

## 5. 実行結果

全体判定は次である。

```text
case_count: 5
total_gauge_count: 130
all_cases_valid: true
single_gauge_only_used: false
robustness_sweep_valid: true
```

主要な最大誤差は次である。

| 量 | 最大値 |
|---|---:|
| `max_p_abs_error_all_cases` | `3.0331293032759277e-13` |
| `max_E_abs_error_all_cases` | `3.0331293032759277e-13` |
| `max_R_abs_error_all_cases` | `1.5765166949677223e-14` |
| `max_R_gauge_std_all_cases` | `5.288392122597181e-15` |
| `max_separation_ratio_time_all_cases` | `1.3338999651354898e-27` |

ケース別結果は次である。

| case | gauges | p max err | E max err | R max err | R std | Var(R)/Var(t) | valid |
|---|---:|---:|---:|---:|---:|---:|---|
| `baseline_default` | `8` | `2.5202062658991053e-14` | `2.2315482794965646e-14` | `4.440892098500626e-16` | `2.220446049250313e-16` | `3.6838616474030686e-30` | `true` |
| `phase_center_grid` | `45` | `2.3137047833188262e-13` | `2.1471713296250527e-13` | `5.773159728050814e-15` | `1.77758999250913e-15` | `2.8276279738980642e-28` | `true` |
| `width_gain_grid` | `36` | `1.0469403122215226e-13` | `6.550315845288424e-14` | `4.440892098500626e-16` | `2.5639502485114184e-16` | `4.911838467416816e-30` | `true` |
| `near_lobe_offset` | `25` | `3.0331293032759277e-13` | `3.0331293032759277e-13` | `1.5765166949677223e-14` | `5.288392122597181e-15` | `1.3338999651354898e-27` | `true` |
| `mixed_readout_grid` | `16` | `1.603162047558726e-13` | `1.3566925360919413e-13` | `1.1102230246251565e-14` | `2.5727487310015434e-15` | `1.0176441594089155e-27` | `true` |

---

## 6. 解釈

本実験では、衝突軌道を変えず、読出し器の設定だけを変えた。

その結果、すべてのゲージケースで、

```text
p_read
E_read
R_read
```

が閾値内で再構成された。

さらに、`R_read` のゲージ分散は最大でも `5.288392122597181e-15` であり、`Var(R)/Var(t)` は最大でも `1.3338999651354898e-27` であった。

したがって、この実装範囲では、質量的量に見える `R_read` は、読出し中心、参照位相、読出し幅、観測機幅、観測ゲインの変更に対して安定である。

これは、`R_read` が単一ゲージの内部値ではなく、多ゲージ干渉読出しで安定に残る成分として読めることを支持する。

---

## 7. 本実験で主張しないこと

本実験は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準物理の質量、運動量、エネルギーを完全導出したこと | 対応写像は別途必要である |
| R軸が先験的な実体軸であること | R は安定残差として読まれる局所表示である |
| 任意の読出し器で成立すること | 本実験は有限個のゲージ族に対する数値検証である |
| 非対称質量衝突まで含むこと | 今回は等振幅・等速度の ABC 完全弾性反射を対象とした |

本実験の主張は、ABC 完全弾性衝突モデルにおいて、読出し器設定を広げても、多ゲージ干渉読出しによる `p/E/R` の再構成と `t/R` 分離が保存された、という数値構成結果である。

---

## 8. 主要出力

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_interference_readout_robustness_sweep_report_v1.md](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_report_v1.md) |
| 結果 JSON | [abc_multigauge_interference_readout_robustness_sweep_result_v1.json](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_result_v1.json) |
| ケース CSV | [abc_multigauge_interference_readout_robustness_sweep_cases_v1.csv](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_cases_v1.csv) |
| ゲージ読出し CSV | [abc_multigauge_interference_readout_robustness_sweep_gauge_rows_v1.csv](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_gauge_rows_v1.csv) |
| ステージ要約 CSV | [abc_multigauge_interference_readout_robustness_sweep_stage_summary_v1.csv](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_stage_summary_v1.csv) |
| 誤差図 | [abc_multigauge_interference_readout_robustness_sweep_errors_v1.png](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_errors_v1.png) |
| 安定性図 | [abc_multigauge_interference_readout_robustness_sweep_stability_v1.png](abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_stability_v1.png) |

