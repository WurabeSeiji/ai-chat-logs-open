# ABC多ゲージ干渉読出しの複数回衝突検証メモ v1

**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、単回衝突で成立した多ゲージ干渉読出しを、複数回の ABC 完全弾性衝突へ拡張して検査した。

前段階では、単回の AB 衝突において、複数ゲージから

```text
p_read
E_read
R_read
```

を再構成し、`p_read` は反射で符号反転し、`E_read` と `R_read` は保存されることを確認した。

本検証では、同じ読出し器を 8 回の AB 衝突へ適用し、衝突反復と壁反射を含んでも、同じ読出し関係が保たれるかを調べた。

---

## 2. 実行スクリプト

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_multi_collision_v2.py
```

出力先は次である。

```text
abc_multigauge_interference_readout_multi_collision_result_v2/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_interference_readout_multi_collision_v2.py
```

---

## 3. 実験構成

実験条件は、既存の ABC 複数回衝突モデルに合わせた。

| 項目 | 値 |
|---|---:|
| AB 衝突目標回数 | `8` |
| 実際の AB 衝突回数 | `8` |
| 壁反射回数 | `7` |
| `A_A,A_B` | `1,1` |
| `q_A,q_B` | `+1,-1` |
| `omega_A,omega_B` | `1,1` |
| `m_A,m_B` | `1,2` |

各 AB 衝突について、衝突直前と衝突直後の状態を保存し、複数ゲージで `p_read`, `E_read`, `R_read` を再構成した。

---

## 4. 判定

主要判定は次である。

| 項目 | 結果 |
|---|---:|
| `completed_target_collisions` | `true` |
| `q_reversed_each_collision` | `true` |
| `label_preserved_all` | `true` |
| `closure_preserved_all` | `true` |
| `p_reconstructed_all_gauges` | `true` |
| `E_reconstructed_all_gauges` | `true` |
| `R_reconstructed_all_gauges` | `true` |
| `p_reflection_each_collision` | `true` |
| `E_preserved_each_collision` | `true` |
| `R_preserved_each_collision` | `true` |
| `R_gauge_stable` | `true` |
| `t_R_separation_valid` | `true` |
| `single_gauge_only_used` | `false` |
| `multi_collision_multigauge_valid` | `true` |

主要数値は次である。

| 量 | 値 |
|---|---:|
| `p_max_abs_error` | `2.5202062658991053e-14` |
| `E_max_abs_error` | `3.3417713041217212e-13` |
| `R_max_abs_error` | `5.6399329650957952e-14` |
| `max_p_reflection_error` | `4.4408920985006262e-16` |
| `max_E_preservation_error` | `0.0000000000000000e+00` |
| `max_R_preservation_error` | `0.0000000000000000e+00` |
| `R_max_gauge_std` | `2.4422383177260761e-14` |
| `closure_max_residual_abs` | `0.0000000000000000e+00` |
| `separation_ratio_time` | `2.7083289874897587e-28` |

---

## 5. 読出しの挙動

各 AB 衝突で、A の `p_read` は `+1` から `-1` へ、B の `p_read` は `-1` から `+1` へ反転した。

最大反転誤差は、

```text
4.4408920985006262e-16
```

であった。

一方、`E_read` は各衝突の前後で保存され、最大保存誤差は

```text
0.0
```

であった。

`R_read` も各衝突の前後で保存され、最大保存誤差は

```text
0.0
```

であった。

複数ゲージにおける `R_read` の最大標準偏差は、

```text
2.4422383177260761e-14
```

であった。

---

## 6. 本検証で主張しないこと

本検証は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準物理量との完全対応 | 対応写像は別途必要である |
| 壁反射を実在境界条件として一般化すること | 本検証では複数回衝突を生成するための数値条件である |
| 単一ゲージ値の測定成立 | 本検証でも複数ゲージ再構成を要求する |
| 非対称質量衝突の成立 | 本検証は対称 ABC 条件の反復検査である |

---

## 7. 結論

本検証では、多ゲージ干渉読出しを 8 回の AB 完全弾性衝突へ適用した。

結果として、衝突反復および壁反射を含んでも、次が保たれた。

1. AB 衝突は 8 回完了した。
2. 各 AB 衝突で `q_A,q_B` は反転した。
3. 識別振動 `m_A,m_B` は保存された。
4. 補償付き二乗閉鎖残差は `0` に保たれた。
5. `p_read` は各衝突で符号反転した。
6. `E_read` は各衝突で保存された。
7. `R_read` は各衝突で保存された。
8. `R_read` は複数ゲージで安定に読まれた。
9. `separation_ratio_time=2.7083289874897587e-28` により、t/R 分離は保たれた。

したがって、本検証の範囲では、単回衝突で成立した多ゲージ干渉読出しは、複数回の ABC 完全弾性衝突に対しても保存的に維持された。

---

# 付録A. 出力ファイル

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_interference_readout_multi_collision_report_v2.md](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_report_v2.md) |
| 結果 JSON | [abc_multigauge_interference_readout_multi_collision_result_v2.json](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_result_v2.json) |
| タイムライン CSV | [abc_multigauge_interference_readout_multi_collision_timeline_v2.csv](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_timeline_v2.csv) |
| イベント CSV | [abc_multigauge_interference_readout_multi_collision_events_v2.csv](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_events_v2.csv) |
| ゲージ CSV | [abc_multigauge_interference_readout_multi_collision_gauge_sweep_v2.csv](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_gauge_sweep_v2.csv) |
| ステージ要約 CSV | [abc_multigauge_interference_readout_multi_collision_stage_summary_v2.csv](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_stage_summary_v2.csv) |
| 衝突読出し CSV | [abc_multigauge_interference_readout_multi_collision_readouts_v2.csv](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_readouts_v2.csv) |
| p/E/R 図 | [abc_multigauge_interference_readout_multi_collision_invariants_v2.png](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_invariants_v2.png) |
| t/R 分離図 | [abc_multigauge_interference_readout_multi_collision_tr_separation_v2.png](abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_tr_separation_v2.png) |

