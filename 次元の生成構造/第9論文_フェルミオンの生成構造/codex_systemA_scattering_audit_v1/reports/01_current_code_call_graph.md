# 01 現行 System A 呼出グラフ

## 1. モジュール読込

**[コード上の事実]**

```text
System A CLI
  run_system_A_localization_exchange_R_sweep_preliminary_v1.py
    └─ importlib により同配下 20260713 散乱源をロード
       run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
```

行番号:

```text
System A:18-45
20260713散乱源:14-32
```

System A を通常 import すると、次の副作用が発生する。

1. System A の既定出力ディレクトリと `.matplotlib` を作る。`System A:20-27`
2. 20260713 散乱源を動的 import する。`System A:35-45`
3. 散乱源側の既定出力ディレクトリと `.matplotlib` も作る。`散乱源:14-32`

このため Stage A では原本を import せず、静的読取りだけを行った。

## 2. System A 主実行経路

**[コード上の事実]**

```text
__main__                                      System A:1379-1466
  ├─ parse_args                              System A:540-579
  ├─ selected_cases_from_args                System A:417-459
  │   ├─ odd_kernel_case                     System A:221-239
  │   ├─ explicit_packet_case                System A:242-280
  │   └─ built_in_packet_cases               System A:287-292
  ├─ selected_r_values_from_args             System A:460-477
  └─ run                                     System A:1305-1376
      ├─ build_source_params                 System A:295-301
      ├─ MetricContext                       System A:304-345
      ├─ run_case                            System A:718-737
      │   ├─ delta_from_reflection_rate      散乱源:140-143
      │   ├─ scattering_coefficients         散乱源:134-137
      │   ├─ make_case_state                 System A:631-659
      │   │   ├─ make_state                  散乱源:88-101
      │   │   └─ make_explicit_packet_state  System A:611-628
      │   ├─ row_for_state                   System A:662-715
      │   └─ 衝突ループ                      System A:729-736
      ├─ summarize_case                      System A:747-828
      ├─ best_rows_for_pair                  System A:862-932
      ├─ terrain_rows                        System A:1025-1103
      └─ CSV/JSON/Markdown/PNG 保存          System A:1340-1375
```

## 3. 入力状態の生成

### 3.1 odd-kernel 系

**[コード上の事実]** `make_case_state` は `state_family == "odd_kernel"` のとき、散乱源の `make_state` を呼ぶ。

```text
System A:631-659
散乱源:88-101
```

`make_state` は、

```text
chi_part = odd_harmonic_kernel(chi - chi_center, n_chi)
phase_chi = exp(i q p0 (chi - chi_center))
eta_phase = exp(i m eta)
psi = chi_part * phase_chi * eta_phase
```

を作り、平坦化した後に単独ノルムを1へ正規化する。

### 3.2 explicit-packet 系

**[コード上の事実]** `make_explicit_packet_state` は指定倍音から、実コード上は正確に

\[
K(u)=
\frac{\sum_j w_j\cos\!\left(\frac{n_j}{\lambda_j}(u-s)+\phi_j\right)}
{\sqrt{\sum_j w_j^2}}
\]

として作り、同じ搬送位相と \(\eta\) 位相を掛けて単独正規化する。

```text
System A:593-628
```

## 4. 散乱ループ

**[コード上の事実]**

```text
R_input
  → delta_f
  → (t, r, T, R)
  → a, b を生成
  → collision = 0 ... max_collision
       ├─ a, b の評価行を作る
       └─ 同じ更新前 a,b から
          a_next = normalize(r a + t b)
          b_next = normalize(t a + r b)
          (a,b) = (a_next,b_next)
```

根拠:

```text
System A:718-737
```

衝突開始条件は、選択された各 `(case, R)` について `collision=0` の状態を最初に記録し、その後更新することである。衝突回数上限は既定256で、CLI `--max-collision` により上書きされる。

## 5. 評価経路

**[コード上の事実]**

各衝突・各チャネルで `row_for_state` が次を計測する。

- `L`: \(\sum_i p_i^2\) 型の局在指標
- `N_eff`, `N_eff_2`: FFT絶対モード分布の一次・二次モーメント
- `p_chi`: FFT周波数の期待値
- `origin_A`, `origin_B`: 初期状態への射影重み
- `sim_to_A0`, `sim_to_B0`: 倍音パワー分布同士のコサイン類似度

```text
System A:662-715
MetricContext:304-345
散乱源 localization:264-269
散乱源 projection_weight:387-394
```

`B_to_A_transfer` は振幅経路ノルムではない。Aチャネルの `sim_to_B0`、すなわち初期Bの**倍音パワー分布との類似度**である。

```text
System A:683-714
System A:761-777
```

反射率・透過率として記録される `R`,`T` は、出力状態から測った量ではなく、`r,t` から返された \(|r|^2,|t|^2\) である。生の反射経路、透過経路、交換経路の配列やノルムは保存されない。

## 6. System Aの保存物

**[コード上の事実]** 基底版は次を保存する。

| 形式 | 内容 | 箇所 |
|---|---|---|
| CSV | 全衝突行、case要約、R別best、collision terrain | System A:1339-1371 |
| JSON | パラメータ、case、R列、best、summary、出力名 | System A:1354-1373 |
| Markdown | best表を中心とする自動報告 | System A:1227-1285,1374-1375 |
| PNG | scoreおよびgap-depth分布（有効時） | System A:1131-1219,1346-1350 |

SVGは保存しない。基底版は複素状態配列も保存しない。

## 7. 最終判定

**[コード上の事実]**

各 `(case,R)` で、

- `L_gap` 最小
- `N_eff_gap` 最小
- `B_to_A_transfer` 最大
- 正規化した上記3量から作る `joint_score` 最小

を探索する。

```text
System A:747-828
```

各caseのR選択列は主に次である。

```text
R_star_L
R_star_N
R_star_transfer
R_star_joint
joint_score_min
R_band_width_5
R_band_width_10
```

```text
System A:862-932
```

## 8. System Bとの共有

**[コード上の事実]** 第9論文配下のSystem B比較コードは、同配下のSystem Aを動的ロードし、その `src.scattering_coefficients` を使用する。

```text
run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py:18-46
同:90-109
```

System Bの状態は空間配列ではなく、複素スカラー対 `(a,b)` である。したがって現在のSystem Bには、波形の偶奇倍音構造を入力する状態自由度自体がない。

System Aには `gray_error`,`gray_depth` は存在しない。これらは比較対象System Bだけで次のように定義される。

\[
\mathrm{gray\_error}
=|S_{\mathrm{mean}}|
+|S_{\mathrm{amp}}-S_{\mathrm{amp,target}}|
+S_{\mathrm{drift}}
+\mathrm{phase\_penalty},
\]

\[
\mathrm{gray\_depth}
=-\log_{10}\max(\mathrm{gray\_error},10^{-300}).
\]

```text
run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py:232-279
```

## 9. 呼出経路外

**[コード上の事実]** 配下外の `run_minimal_system_B_gray_direct_check_v5.py` は同じ係数式を自己完結で再定義するが、現行System Aから呼ばれない。`phase5_eigenphase_resonance_v2.py` と `run_two_physical_roots_multiprecision_v1.py` も現行System Aのimport graph外である。

これらは本監査で実行・importしていない。
