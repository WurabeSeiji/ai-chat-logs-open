# ABC閉鎖位相系における多ゲージ干渉読出し保存量の構成実験

**副題:** `p_read`, `E_read`, `R_read` と `R*p`, `R*p^2` 保存写像による質量的・運動量的・エネルギー的量の数値構成  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・追加論文  
**Version DOI:** pending  
**Concept DOI:** pending  

---

## 要旨

本稿では、背景座標を先験的に仮定しない ABC 閉鎖位相系において、質量的量、運動量的量、エネルギー的量に相当する保存読出しが、多ゲージ干渉から構成できるかを数値的に検証した。

本シリーズでは、空間・時間に相当する `χ,τ` を外部時空座標ではなく、閉じた位相系からの読出し量として扱ってきた。本稿では同じ方針を、質量的量、運動量的量、エネルギー的量へ拡張する。

局所波 `A,B` と観測機 `C` からなる ABC 系に対して、単一ゲージ値を測定値とはみなさず、複数の読出しゲージから干渉相関を取得する。空間位相方向の相関勾配を

```text
p_read
```

時間位相方向の相関勾配を

```text
E_read
```

複数ゲージで安定して残る振幅二乗残差を

```text
R_read
```

として読む。

本稿でいう `p,E,R` は、標準物理の運動量、エネルギー、質量そのものではない。閉鎖位相系内の多ゲージ干渉読出し量である。本稿の目的は、それらが標準量と同一であることを主張することではなく、保存量に見える読出し構造が ABC 閉鎖位相系内に構成できるかを調べることである。

実行結果では、単回 ABC 衝突において `p,E,R` が多ゲージから再構成され、最大誤差はそれぞれ `2.5202062658991053e-14`, `2.2315482794965646e-14`, `4.440892098500626e-16` であった。8回の反復衝突でも `p` 反転、`E,R` 保存、識別振動保存、補償付き閉鎖が維持された。

非対称 `R` 条件では、単純な `q -> -q` 反転が `R*p` 保存を破ることが検出された。そこで、`R_A p_A + R_B p_B` および `R_A p_A^2 + R_B p_B^2` を保存する一般化衝突写像を構成した。8種類の非対称振幅ケースでは、最大 `R*p` 保存誤差 `2.3803181647963356e-13`、最大 `R*p^2` 保存誤差 `1.4086509736443986e-12` で保存が確認された。

さらに、9種類の非単位・非対称初期位相勾配、4種類の複数回衝突ケース、読出しノイズ頑健性、`R_B/R_A=0.015625` から `64.0` までの極端 `R` 比スイープを実行した。統合サマリーでは9本の実験すべてが `valid` となり、単一ゲージのみの判定は一つも用いなかった。

以上により、本稿の数値構成範囲では、ABC 閉鎖位相系から、質量的・運動量的・エネルギー的量に見える保存読出しを、多ゲージ干渉により構成できることが確認された。

**キーワード:** 多ゲージ干渉読出し、ABC閉鎖位相系、全正符号ゼロ閉鎖、質量的読出し、運動量的読出し、エネルギー的読出し、一般化弾性衝突、`R*p` 保存、`R*p^2` 保存

---

## 1. 序論

### 1.1 背景

本シリーズでは、無名性、全正符号ゼロ閉鎖、非自明存在を基本公理として採用している。中心となる閉鎖条件は、

```math
\sum_n x_n^2=0
```

である。

この閉鎖条件は、共役ノルム

```math
\sum_n |x_n|^2
```

ではない。各成分をそのまま二乗し、全正符号で総和する閉鎖条件である。

前段階の「全正符号ゼロ閉鎖の読出し多重性に関する定義補足」では、この公理1を変更せず、同じ閉鎖条件が複数の読出し表示を持つことを整理した。半径的表示を

```math
a^2+b^2=\rho^2
```

として読む場合も、第一原理層では外部負符号付きの

```math
a^2+b^2-\rho^2=0
```

を置かない。全正符号ゼロ閉鎖の表示として、

```math
a^2+b^2+(i\rho)^2=0
```

と読む。

したがって、半径、時間、質量、エネルギー、運動量に見える量も、第一原理的に固有名を持つ軸として先に置くのではなく、読出し表示として扱う。

### 1.2 本稿の問い

本稿の問いは次である。

> 背景座標を先験的に仮定しない ABC 閉鎖位相系において、質量的量、運動量的量、エネルギー的量に見える保存読出しを、多ゲージ干渉から構成できるか。

ここで重要なのは、単一ゲージ値を測定値とみなさないことである。

`χ,τ` に相当する読出しが干渉相関から構成されたのと同様に、`p,E,R` に相当する読出しも、複数の参照波、複数の読出し窓、複数のゲージから再構成されなければならない。

### 1.3 本稿の最小主張

本稿の主張は、次の範囲に限定する。

1. ABC 閉鎖位相系内で、`p_read`, `E_read`, `R_read` を多ゲージ干渉から再構成する。
2. 対称衝突では、`p_read` が反転し、`E_read`, `R_read` が保存されることを確認する。
3. 非対称 `R` 条件では、単純反転が `R*p` 保存を破ることを検出する。
4. `R*p` と `R*p^2` を保存する一般化衝突写像を構成する。
5. 反復衝突、読出し器変更、読出しノイズ、極端 `R` 比に対して、読出し保存が維持されるかを検証する。

---

## 2. 本稿で主張しないこと

本稿は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準運動量・標準エネルギー・標準質量の導出 | 本稿の `p,E,R` は閉鎖位相系内の読出し量である |
| 標準力学の完全な再導出 | 対応写像は別途構成が必要である |
| 実在粒子衝突の定量予言 | 本稿は内部位相モデル上の数値構成実験である |
| 単一ゲージ測定の成立 | 本稿では多ゲージ干渉読出しを要求する |
| `R` を先験的な質量軸とすること | `R` は複数ゲージで安定して残る読出し名である |
| 背景時空座標の仮定 | `χ,τ` は位相系内の読出し変数である |
| 任意条件での一般化衝突成立 | 検証範囲は本稿の数値条件に限る |

本稿の主張は、ABC 閉鎖位相系内で、保存量に見える読出し構造を多ゲージ干渉から構成したことである。

---

## 3. 基本構造

### 3.1 全正符号ゼロ閉鎖

基本公理系 v2 の公理1は、

```math
Q(x)=\sum_n x_n^2=0
```

である。

これは、外部の負符号計量を先に置く条件ではない。

最小補償対

```math
A,\qquad iA
```

では、

```math
A^2+(iA)^2=0
```

となる。

負符号は外部係数ではなく、

```math
i^2=-1
```

という内部位相の二乗から現れる。

### 3.2 読出し多重性

同じ閉鎖条件

```math
\sum_n x_n^2=0
```

は、局所読出しにおいて複数の表示を持つ。

例えば、局所表示で

```math
a^2+b^2+(i\rho)^2=0
```

と読める場合、`a,b,\rho` は第一原理的に異なる種類の成分ではない。ある読出し窓、ある参照波、ある射影において便宜的に付けたラベルである。

本稿の `p,E,R` も同様である。

`p` は運動量そのものではなく、空間位相方向の相関勾配読出しである。

`E` はエネルギーそのものではなく、時間位相方向の相関勾配読出しである。

`R` は質量そのものではなく、複数ゲージで安定して残る振幅二乗読出しである。

---

## 4. ABC閉鎖位相系

### 4.1 局所波 A, B

局所波 `A,B` は、空間位相 `χ`、時間位相 `τ`、内部識別位相 `η`、代表振幅、および内部識別モードを持つ。

本稿では、`A,B` の位相中心を `χ_A,χ_B`、時間位相中心を `τ_A,τ_B` とする。

内部識別モードは、

```text
m_A=1,
m_B=2
```

である。

### 4.2 観測機 C

観測機 `C` は、外部観測者ではない。

同じ閉鎖位相系内にある参照波であり、局所波 `A,B` との干渉相関を読むために用いる。

各ゲージは、読出し中心、読出し幅、位相シフト、参照波利得などを変える。

### 4.3 単一ゲージ値を測定値としない

本稿では、単一ゲージで得られた数値を測定成立とはみなさない。

測定成立には、複数ゲージに対して同じ読出し量が再構成されることを要求する。

したがって、測定対象は、個々のゲージ値ではなく、ゲージ群を通じて安定に残る量である。

---

## 5. 多ゲージ干渉読出し

### 5.1 空間位相勾配読出し

空間位相方向の干渉相関を、微小変位 `h_χ` に対して

```math
p_{\mathrm{read}}
=
\frac{\arg\left(O(\chi+h_\chi)/O(\chi-h_\chi)\right)}
{2h_\chi}
```

として読む。

これは標準運動量そのものではない。空間位相方向の相関勾配である。

### 5.2 時間位相勾配読出し

時間位相方向の干渉相関を、微小変位 `h_τ` に対して

```math
E_{\mathrm{read}}
=
\frac{\arg\left(O(\tau+h_\tau)/O(\tau-h_\tau)\right)}
{2h_\tau}
```

として読む。

これは標準エネルギーそのものではない。時間位相方向の相関勾配である。

### 5.3 R読出し

読出し振幅をゲージごとに較正し、振幅二乗として

```math
R_{\mathrm{read}}
=
\gamma_g A_{\mathrm{read}}^2
```

を読む。

ここで `γ_g` はゲージ側の読出し利得である。

`R_read` は、複数ゲージを変えても安定して残る量であることを要求する。

### 5.4 t/R分離

`t` と `R` は先験的に名前付きの別軸ではない。

変動が大きく連続的に読まれる成分を `t` と読み、変動が小さく安定に残る成分を `R` と読む。

本稿では、作業的に

```math
\frac{\operatorname{Var}(R_{\mathrm{read}})}
{\operatorname{Var}(t_{\mathrm{read}})}
```

を用いて、`t/R` 分離を評価する。

この比が十分小さいとき、`R` は時間的変動から分離された安定読出しとして扱える。

---

## 6. 対称ABC衝突における基礎読出し

### 6.1 単回衝突

まず、等振幅条件

```text
A_A=A_B=1
```

で、単回 ABC 衝突を実行した。

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_v1.py
```

結果は次であった。

| 量 | 値 |
|---|---:|
| `p_max_abs_error` | `2.5202062658991053e-14` |
| `E_max_abs_error` | `2.2315482794965646e-14` |
| `R_max_abs_error` | `4.440892098500626e-16` |
| `p_reflection_error_A` | `3.3306690738754696e-16` |
| `p_reflection_error_B` | `3.3306690738754696e-16` |
| `E_preservation_error_A` | `0.0` |
| `E_preservation_error_B` | `0.0` |
| `R_preservation_error_A` | `0.0` |
| `R_preservation_error_B` | `0.0` |
| `separation_ratio_time` | `3.6838616474030686e-30` |

判定は次であった。

```text
multigauge_measurement_valid: true
single_gauge_only_used: false
```

### 6.2 複数回衝突

次に、同じ読出し器を8回の AB 衝突へ適用した。

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_multi_collision_v1.py
```

結果は次であった。

| 量 | 値 |
|---|---:|
| `ab_collision_count` | `8` |
| `wall_reflection_count` | `7` |
| `p_max_abs_error` | `2.5202062658991053e-14` |
| `E_max_abs_error` | `3.341771304121721e-13` |
| `R_max_abs_error` | `5.639932965095795e-14` |
| `max_p_reflection_error` | `4.440892098500626e-16` |
| `max_E_preservation_error` | `0.0` |
| `max_R_preservation_error` | `0.0` |
| `separation_ratio_time` | `2.7083289874897587e-28` |

判定は次であった。

```text
multi_collision_multigauge_valid: true
single_gauge_only_used: false
```

### 6.3 読出し器構成の頑健性

読出しゲージの中心、幅、位相、利得を変えた5種類のゲージ群を用い、合計130ゲージで読出し器頑健性を検査した。

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_robustness_sweep_v1.py
```

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `5` |
| `total_gauge_count` | `130` |
| `max_p_abs_error_all_cases` | `3.0331293032759277e-13` |
| `max_E_abs_error_all_cases` | `3.0331293032759277e-13` |
| `max_R_abs_error_all_cases` | `1.5765166949677223e-14` |
| `max_R_gauge_std_all_cases` | `5.288392122597181e-15` |
| `max_separation_ratio_time_all_cases` | `1.3338999651354898e-27` |

判定は次であった。

```text
robustness_sweep_valid: true
single_gauge_only_used: false
```

---

## 7. 非対称Rにおける単純反転の破綻診断

### 7.1 問題設定

等振幅条件では、`p_A=+1`, `p_B=-1` の単純反転

```text
p_A -> -p_A
p_B -> -p_B
```

は、反射読出しとして成立する。

しかし、`R_A` と `R_B` が異なる場合、同じ単純反転が

```math
R_A p_A+R_B p_B
```

を保存するとは限らない。

したがって、非対称 `R` 条件では、単純反転が保存読出しを破るかを診断する必要がある。

### 7.2 実行結果

実行スクリプトは次である。

```text
run_abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1.py
```

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `8` |
| `asymmetric_case_count` | `7` |
| `individual_multigauge_valid_all_cases` | `true` |
| `weighted_energy_preserved_all_cases` | `true` |
| `R_total_preserved_all_cases` | `true` |
| `equal_case_weighted_momentum_preserved` | `true` |
| `asymmetric_cases_detect_weighted_momentum_failure` | `true` |
| `max_weighted_p_collision_error` | `16.000000000000036` |

この結果により、非対称 `R` 条件では、単純反転が `R*p` 保存を破ることが検出された。

これは破綻ではなく、次の一般化写像が必要であることを示す診断である。

---

## 8. R重み付き一般化衝突写像

### 8.1 保存条件

非対称 `R` 条件では、保存読出しとして次を要求する。

```math
P_R
=
R_Ap_A+R_Bp_B
```

```math
K_R
=
R_Ap_A^2+R_Bp_B^2
```

ここで `P_R` は運動量的読出し、`K_R` は位相勾配二乗の保存読出しである。

標準物理の運動量・エネルギーと同一視しない。

### 8.2 写像

`P_R` と `K_R` を保存し、相対位相勾配を反転する写像を、

```math
p_A'
=
\frac{R_A-R_B}{R_A+R_B}p_A
+
\frac{2R_B}{R_A+R_B}p_B
```

```math
p_B'
=
\frac{2R_A}{R_A+R_B}p_A
+
\frac{R_B-R_A}{R_A+R_B}p_B
```

とする。

この写像は、結果値として `R,T` を代入するものではない。多ゲージ干渉から読まれた `R_read` と `p_read` に対して、保存関係が成立するかを確認するための局所写像である。

### 8.3 一般化写像の検証

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_readout_v1.py
```

8種類の振幅条件で検証した。

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `8` |
| `individual_readout_valid_all_cases` | `true` |
| `generalized_P_R_preserved_all_cases` | `true` |
| `generalized_K_R_phase_preserved_all_cases` | `true` |
| `E_tau_R_preserved_all_cases` | `true` |
| `R_total_preserved_all_cases` | `true` |
| `max_P_R_conservation_error` | `2.3803181647963356e-13` |
| `max_K_R_phase_conservation_error` | `1.4086509736443986e-12` |

判定は次であった。

```text
generalized_elastic_collision_readout_valid: true
single_gauge_only_used: false
```

---

## 9. 非単位・非対称位相勾配スイープ

### 9.1 目的

前節では、主に `p_A=+1`, `p_B=-1` を基準とした。

本節では、初期位相勾配を非単位・非対称にし、対向衝突だけでなく同方向追いつき条件も含めて検証する。

### 9.2 実行結果

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1.py
```

9種類の初期条件を検証した。

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `9` |
| `collision_reached_all_cases` | `true` |
| `P_R_preserved_all_cases` | `true` |
| `K_R_phase_preserved_all_cases` | `true` |
| `relative_gradient_flipped_all_cases` | `true` |
| `E_tau_R_preserved_all_cases` | `true` |
| `R_total_preserved_all_cases` | `true` |
| `max_P_R_conservation_error` | `2.8910207561239076e-13` |
| `max_K_R_phase_conservation_error` | `1.5258905250448151e-12` |
| `max_relative_flip_error` | `1.2434497875801753e-14` |

判定は次であった。

```text
velocity_sweep_generalized_collision_valid: true
single_gauge_only_used: false
```

---

## 10. 一般化写像の複数回衝突

### 10.1 目的

単回で成立した一般化衝突写像が、反復過程でも維持されるかを確認した。

壁反射は粒子を再遭遇させるための補助条件であり、保存判定は各 AB 衝突の直前直後に限定した。

### 10.2 実行結果

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_multi_collision_v1.py
```

4種類の条件で、各6回の AB 衝突を実行した。

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `4` |
| `max_ab_collision_count` | `6` |
| `max_wall_reflection_count` | `8` |
| `P_R_preserved_each_collision_all_cases` | `true` |
| `K_R_preserved_each_collision_all_cases` | `true` |
| `relative_flip_each_collision_all_cases` | `true` |
| `E_tau_R_preserved_each_collision_all_cases` | `true` |
| `R_preserved_each_collision_all_cases` | `true` |
| `max_P_R_error` | `3.552713678800501e-14` |
| `max_K_R_error` | `1.056932319443149e-13` |
| `max_relative_flip_error` | `2.220446049250313e-14` |

判定は次であった。

```text
generalized_multi_collision_valid: true
single_gauge_only_used: false
```

---

## 11. 読出しノイズ頑健性

### 11.1 目的

測定成立には、単一ゲージ値ではなく多ゲージ干渉が必要である。

したがって、読出し器側の揺らぎがある場合、ゼロ平均ゲージ揺らぎは相殺され、全ゲージ共通のバイアスは検出される必要がある。

本節では、状態シミュレーション後の読出し行に決定論的擬似ノイズを加えた。

### 11.2 実行結果

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_noise_robustness_v1.py
```

4ケース、2種類のノイズモード、6段階のノイズレベルで検証した。

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `4` |
| `noise_mode_count` | `2` |
| `noise_level_count` | `6` |
| `max_gauge_count` | `116` |
| `zero_mean_multigauge_valid_all` | `true` |
| `common_bias_detection_floor` | `1e-10` |
| `common_bias_detected_all_above_floor` | `true` |
| `zero_mean_max_p_mean_abs_error` | `1.4477308241112041e-13` |
| `zero_mean_max_R_mean_abs_error` | `3.552713678800501e-14` |
| `zero_mean_max_K_R_error` | `8.322231792590173e-13` |

判定は次であった。

```text
noise_robustness_valid: true
single_gauge_only_used: false
```

---

## 12. 極端R比スイープ

### 12.1 目的

`R_read` が質量的量に見える場合、`R` 比が大きく非対称な条件でも保存読出しが維持されるかを調べる必要がある。

本節では、

```text
R_B/R_A = 0.015625
```

から

```text
R_B/R_A = 64.0
```

まで掃引した。

### 12.2 実行結果

実行スクリプトは次である。

```text
run_abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v1.py
```

結果は次であった。

| 量 | 値 |
|---|---:|
| `case_count` | `12` |
| `max_R_dynamic_range` | `64.0` |
| `min_R_ratio_B_over_A` | `0.015625` |
| `max_R_ratio_B_over_A` | `64.0` |
| `P_R_preserved_all_cases` | `true` |
| `K_R_phase_preserved_all_cases` | `true` |
| `relative_gradient_flipped_all_cases` | `true` |
| `max_P_R_conservation_error` | `6.465938895416912e-13` |
| `max_K_R_phase_conservation_error` | `1.2789769243681803e-12` |
| `max_relative_flip_error` | `2.6867397195928788e-14` |

判定は次であった。

```text
extreme_R_sweep_valid: true
single_gauge_only_used: false
```

---

## 13. 統合サマリー

9本の実験を統合した。

実行スクリプトは次である。

```text
run_abc_multigauge_readout_integration_summary_v1.py
```

統合判定は次であった。

```text
experiment_count: 9
all_experiments_valid: true
single_gauge_only_used_any: false
integration_summary_valid: true
```

実験一覧は次である。

| 実験 | 目的 | valid |
|---|---|---|
| `single_collision_multigauge_readout` | 単回ABC衝突で `p/E/R` を多ゲージ干渉読出しする | `true` |
| `multi_collision_multigauge_readout` | 対称ABC衝突の反復で `p/E/R` 読出しを維持する | `true` |
| `readout_robustness_sweep` | 複数の読出し器構成で `p/E/R` 再構成が安定する | `true` |
| `asymmetric_amplitude_diagnostic` | 非対称Rで単純反転が保存を破ることを検出する | `true` |
| `generalized_elastic_collision_readout` | 非対称Rで `R*p` と `R*p^2` を保存する一般化写像を読む | `true` |
| `generalized_velocity_sweep` | 非単位・非対称位相勾配でも一般化写像が成立する | `true` |
| `generalized_multi_collision` | 一般化写像を複数回AB衝突へ反復適用する | `true` |
| `generalized_noise_robustness` | ゼロ平均読出しノイズの相殺と共通バイアス検出を確認する | `true` |
| `generalized_extreme_R_sweep` | 極端なR比でも一般化写像と読出しが維持されるか調べる | `true` |

---

## 14. 評価

本稿の結果を、探究型物理学者ロールの分類で整理する。

| 対象 | 分類 | 判定 |
|---|---|---|
| `p_read` が多ゲージ干渉から再構成される | 数値構成済み帰結 | 保持 |
| `E_read` が多ゲージ干渉から再構成される | 数値構成済み帰結 | 保持 |
| `R_read` が複数ゲージで安定残差として読まれる | 数値構成済み帰結 | 保持 |
| 対称衝突で `p` 反転、`E,R` 保存が成り立つ | 数値構成済み帰結 | 保持 |
| 非対称Rで単純反転が `R*p` 保存を破る | 数値構成済み帰結 | 保持 |
| `R*p`, `R*p^2` 保存写像が非対称Rで成立する | 数値構成済み帰結 | 保持 |
| 読出しノイズのゼロ平均成分が多ゲージ平均で相殺される | 数値構成済み帰結 | 保持 |
| 共通読出しバイアスが検出される | 数値構成済み帰結 | 保持 |
| 標準質量・標準運動量・標準エネルギーとの同一視 | 対応写像の課題 | 未主張 |

---

## 15. 考察

### 15.1 質量・運動量・エネルギーを先に置かない

本稿では、質量、運動量、エネルギーを先に置かない。

先に置くのは、閉鎖位相系、局所波、観測機、ゲージ群、干渉相関である。

そこから、

```text
p_read
E_read
R_read
```

が読まれる。

したがって、本稿の順序は、

```text
質量・運動量・エネルギーを仮定する
```

ではなく、

```text
干渉読出しから、質量的・運動量的・エネルギー的に振る舞う保存量が出るかを調べる
```

である。

### 15.2 Rは測りにくい

`p_read` は空間位相勾配であり、符号反転や相対勾配反転として比較的読みやすい。

`E_read` は時間位相方向の勾配であり、時間窓を変えた相関から読める。

一方、`R_read` は安定残差である。

これは、変動が小さいからこそ安定している。しかし、変動が小さいために測りにくい。

本稿の `t/R` 分離は、この点を数値的に確認するための作業指標である。

### 15.3 単一ゲージでは足りない

本稿の全実験で、

```text
single_gauge_only_used: false
```

であった。

これは重要である。

`p,E,R` は、単一の局所値として直接見える量ではない。複数ゲージを通じて安定に再構成される量である。

したがって、質量的量、運動量的量、エネルギー的量の読出しは、単一ゲージ測定ではなく、多ゲージ干渉再構成として扱う必要がある。

### 15.4 単純反転から一般化保存写像へ

等振幅条件では、単純な方向反転が保存的に見える。

しかし、非対称 `R` 条件では、単純反転は `R*p` 保存を破る。

この破れは、モデルの破綻ではない。むしろ、`R` が質量的に読まれるなら、保存写像を `R` 重み付きに一般化する必要があることを示す診断である。

本稿では、その一般化写像が `R*p` と `R*p^2` を保存することを確認した。

### 15.5 標準理論との接続は次の課題である

本稿の構成は、標準物理の運動量・エネルギー・質量を直接導出するものではない。

しかし、次の構造は得られた。

```text
空間位相勾配 -> 運動量的読出し
時間位相勾配 -> エネルギー的読出し
安定振幅二乗残差 -> 質量的読出し
R*p 保存 -> 運動量的保存に見える関係
R*p^2 保存 -> 二乗量保存に見える関係
```

したがって、標準理論との対応写像を構成するための読出し側の土台は得られた。

---

## 16. 結論

本稿では、ABC 閉鎖位相系において、質量的量、運動量的量、エネルギー的量に相当する保存読出しが、多ゲージ干渉から構成できるかを数値的に検証した。

単回 ABC 衝突では、`p_read`, `E_read`, `R_read` が複数ゲージから再構成され、`p` は反転し、`E,R` は保存された。最大読出し誤差は、`p` で `2.5202062658991053e-14`、`E` で `2.2315482794965646e-14`、`R` で `4.440892098500626e-16` であった。

8回の対称反復衝突でも、`p` 反転、`E,R` 保存、識別振動保存、補償付き閉鎖が保たれた。

非対称 `R` 条件では、単純な `q -> -q` 反転が `R*p` 保存を破ることが検出された。これを受けて、`R*p` と `R*p^2` を保存する一般化衝突写像を構成した。8種類の非対称振幅ケースでは、最大 `R*p` 保存誤差 `2.3803181647963356e-13`、最大 `R*p^2` 保存誤差 `1.4086509736443986e-12` で保存が確認された。

さらに、非単位・非対称初期位相勾配、同方向追いつき衝突、複数回衝突、読出しノイズ、極端 `R` 比スイープを実行した。統合サマリーでは、9本の実験すべてが `valid` となり、単一ゲージのみの判定は一つも用いなかった。

以上により、本稿の数値構成範囲では、ABC 閉鎖位相系から、

```text
p_read
E_read
R_read
R*p
R*p^2
```

に相当する保存読出しを、多ゲージ干渉により一貫して構成できることが確認された。

これは標準物理量との同一視ではない。しかし、質量的量、運動量的量、エネルギー的量を、先験的な実体量ではなく、閉鎖位相系からの保存読出しとして扱うための数値構成を与える。

---

# 付録A. 実行済みプログラムと出力

## A.1 単回ABC多ゲージ読出し

```text
python3 run_abc_multigauge_interference_readout_v1.py
```

出力:

```text
abc_multigauge_interference_readout_result_v1/
```

主要ファイル:

| 種類 | ファイル |
|---|---|
| レポート | [abc_multigauge_interference_readout_report_v1.md](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_report_v1.md) |
| JSON | [abc_multigauge_interference_readout_result_v1.json](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_result_v1.json) |
| gauge CSV | [abc_multigauge_interference_readout_gauge_sweep_v1.csv](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_gauge_sweep_v1.csv) |
| p/E/R 図 | [abc_multigauge_interference_readout_invariants_v1.png](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_invariants_v1.png) |
| t/R 分離図 | [abc_multigauge_interference_readout_tr_separation_v1.png](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_tr_separation_v1.png) |

## A.2 対称複数回衝突

```text
python3 run_abc_multigauge_interference_readout_multi_collision_v1.py
```

出力:

```text
abc_multigauge_interference_readout_multi_collision_result_v1/
```

主要ファイル:

| 種類 | ファイル |
|---|---|
| レポート | [abc_multigauge_interference_readout_multi_collision_report_v1.md](abc_multigauge_interference_readout_multi_collision_result_v1/abc_multigauge_interference_readout_multi_collision_report_v1.md) |
| JSON | [abc_multigauge_interference_readout_multi_collision_result_v1.json](abc_multigauge_interference_readout_multi_collision_result_v1/abc_multigauge_interference_readout_multi_collision_result_v1.json) |
| 衝突読出し CSV | [abc_multigauge_interference_readout_multi_collision_readouts_v1.csv](abc_multigauge_interference_readout_multi_collision_result_v1/abc_multigauge_interference_readout_multi_collision_readouts_v1.csv) |

## A.3 読出し器頑健性

```text
python3 run_abc_multigauge_interference_readout_robustness_sweep_v1.py
```

出力:

```text
abc_multigauge_interference_readout_robustness_sweep_result_v1/
```

## A.4 非対称振幅診断

```text
python3 run_abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1.py
```

出力:

```text
abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/
```

## A.5 一般化弾性衝突写像

```text
python3 run_abc_multigauge_generalized_elastic_collision_readout_v1.py
```

出力:

```text
abc_multigauge_generalized_elastic_collision_readout_result_v1/
```

## A.6 非対称速度スイープ

```text
python3 run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1.py
```

出力:

```text
abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/
```

## A.7 一般化複数回衝突

```text
python3 run_abc_multigauge_generalized_elastic_collision_multi_collision_v1.py
```

出力:

```text
abc_multigauge_generalized_elastic_collision_multi_collision_result_v1/
```

## A.8 読出しノイズ頑健性

```text
python3 run_abc_multigauge_generalized_elastic_collision_noise_robustness_v1.py
```

出力:

```text
abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1/
```

## A.9 極端R比スイープ

```text
python3 run_abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v1.py
```

出力:

```text
abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v1/
```

## A.10 統合サマリー

```text
python3 run_abc_multigauge_readout_integration_summary_v1.py
```

出力:

```text
abc_multigauge_readout_integration_summary_result_v1/
```

---

# 付録B. 実行済み検証メモ

| 種類 | ファイル |
|---|---|
| 定義補足 | [全正符号ゼロ閉鎖の読出し多重性に関する定義補足.md](全正符号ゼロ閉鎖の読出し多重性に関する定義補足.md) |
| 仕様方針 | [現在チャットメモ_多ゲージ干渉読出し仕様方針.md](現在チャットメモ_多ゲージ干渉読出し仕様方針.md) |
| 単回検証 | [ABC完全弾性衝突における多ゲージ干渉読出し数値検証メモ_v1.md](ABC完全弾性衝突における多ゲージ干渉読出し数値検証メモ_v1.md) |
| 対称複数回衝突 | [ABC多ゲージ干渉読出しの複数回衝突検証メモ_v1.md](ABC多ゲージ干渉読出しの複数回衝突検証メモ_v1.md) |
| 読出し器頑健性 | [ABC多ゲージ干渉読出しの読出し器頑健性スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しの読出し器頑健性スイープ検証メモ_v1.md) |
| 非対称振幅診断 | [ABC多ゲージ干渉読出しの非対称振幅診断スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しの非対称振幅診断スイープ検証メモ_v1.md) |
| 一般化写像 | [ABC多ゲージ干渉読出しによる一般化弾性衝突写像検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突写像検証メモ_v1.md) |
| 非対称速度 | [ABC多ゲージ干渉読出しによる一般化弾性衝突の非対称速度スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の非対称速度スイープ検証メモ_v1.md) |
| 一般化複数回 | [ABC多ゲージ干渉読出しによる一般化弾性衝突の複数回衝突検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の複数回衝突検証メモ_v1.md) |
| ノイズ頑健性 | [ABC多ゲージ干渉読出しによる一般化弾性衝突の読出しノイズ頑健性検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の読出しノイズ頑健性検証メモ_v1.md) |
| 極端R比 | [ABC多ゲージ干渉読出しによる一般化弾性衝突の極端R比スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の極端R比スイープ検証メモ_v1.md) |
| 統合サマリー | [ABC多ゲージ干渉読出し実験群の統合サマリー_v1.md](ABC多ゲージ干渉読出し実験群の統合サマリー_v1.md) |

---

# 参考文献

## 自己引用

1. 木原範昭「無名等振幅複合波モデル 基本公理系 v2」2026-07-10.
2. 木原範昭「全正符号ゼロ閉鎖の読出し多重性に関する定義補足」2026-07-11.
3. 木原範昭「背景空間を仮定しない閉じた位相系におけるフェルミオン的二局所波の完全弾性反射の構成実験」2026.
4. 木原範昭「フェルミオン的逆相核による完全反射写像の干渉構成」Version DOI: `10.5281/zenodo.21295480`, Concept DOI: `10.5281/zenodo.21295479`, 2026.
5. 木原範昭「曲率付き閉鎖定常波による曲率繰り込みと完全反射安定性」Version DOI: `10.5281/zenodo.21304040`, Concept DOI: `10.5281/zenodo.21304039`, 2026.

## 外部参考文献

6. S. Pancharatnam, “Generalized theory of interference, and its applications,” *Proceedings of the Indian Academy of Sciences A*, 44, 247–262, 1956.
7. M. V. Berry, “Quantal phase factors accompanying adiabatic changes,” *Proceedings of the Royal Society of London A*, 392, 45–57, 1984. DOI: `10.1098/rspa.1984.0023`.
