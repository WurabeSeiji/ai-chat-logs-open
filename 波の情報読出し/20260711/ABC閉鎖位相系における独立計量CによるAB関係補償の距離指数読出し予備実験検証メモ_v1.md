# ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 目的

本メモでは、仕様書

[ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し実験仕様書 v1.md](ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し実験仕様書%20v1.md)

に基づく最初の予備実験を記録する。

AB 二体閉鎖位相系では、相対距離の変化は読めるが、その変化が

```text
L_AB^0
1/L_AB
1/L_AB^2
```

のどれに従うかを独立に計量するゲージがなかった。

そこで、本実験では第三波 `C` を導入し、`C` が AB の位置変化量と時間位相変化量を読むための独立計量ゲージになれるかを検査した。

ただし、本予備実験では、逆比例項や逆二乗項を実装しない。

目的は、まず

```text
C ゲージに有効窓があるか
```

を調べることである。

---

## 2. 実行コマンド

```text
python3 run_abc_c_gauge_ab_distance_exponent_preliminary_v1.py
```

出力:

```text
abc_c_gauge_ab_distance_exponent_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [abc_c_gauge_ab_distance_exponent_preliminary_report_v1.md](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_preliminary_report_v1.md) |
| JSON | [abc_c_gauge_ab_distance_exponent_preliminary_result_v1.json](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_preliminary_result_v1.json) |
| ゲージ適格性 CSV | [abc_c_gauge_ab_distance_exponent_gauge_validity_v1.csv](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_gauge_validity_v1.csv) |
| ケース CSV | [abc_c_gauge_ab_distance_exponent_cases_v1.csv](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_cases_v1.csv) |
| ゲージ適格性図 | [abc_c_gauge_ab_distance_exponent_gauge_eligibility_map_v1.png](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_gauge_eligibility_map_v1.png) |
| alpha 候補図 | [abc_c_gauge_ab_distance_exponent_alpha_candidates_v1.png](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_alpha_candidates_v1.png) |
| 参照曲線図 | [abc_c_gauge_ab_distance_exponent_reference_curve_v1.png](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_reference_curve_v1.png) |
| フィルタ図 | [abc_c_gauge_ab_distance_exponent_validity_filters_v1.png](abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_validity_filters_v1.png) |

---

## 4. 実験設計

### 4.1 掃引パラメータ

本実験では次を掃引した。

```text
R_C / R_A
C coupling
C offset
initial AB deviation
```

`R_C / R_A` は、`C` の計量セル密度にも対応させた。

すなわち、本予備モデルでは、

```text
R_C が小さい
= C の振動数が低い
= C の波長が長い
= C のセル間隔が広い
```

として扱った。

このため、`R_C` が小さすぎると、`C` は AB の位置変化量を複数セルで読めない。

### 4.2 禁止したこと

本実験では次を実装していない。

```text
1/L_AB
1/L_AB^2
1/A_chi_tau
F = G m_A m_B / L_AB^2
```

したがって、逆比例型または逆二乗型が出る場合は、実装で作ったものではなく、C ゲージ読出しの結果でなければならない。

### 4.3 判定条件

C ゲージを距離指数判定に使うには、次をすべて満たす必要がある。

| 判定 | 内容 |
|---|---|
| `resolution_valid` | C が AB の位置変化量を複数セルで読める |
| `clock_valid` | C の時間位相ゲージが粗すぎない |
| `disturbance_valid` | C 由来の `f_AC`, `f_BC`, `f_ABC` が AB 主読出しを支配しない |
| `closure_valid` | 閉鎖残差が許容範囲に戻る |

---

## 5. 結果サマリー

主要結果は次である。

```text
abc_c_gauge_ab_distance_exponent_preliminary_valid: true
config_count: 160
gauge_valid_count: 41
resolution_valid_count: 60
clock_valid_count: 100
disturbance_valid_count: 101
inverse_like_alpha1_count: 0
inverse_square_like_alpha2_count: 0
proportional_like_alpha_minus1_count: 41
```

関係性時間ごとの alpha 分類は次であった。

```text
tau_ABC: proportional_like_alpha_minus1 = 41
tau_AB:  other = 41
tau_AC:  proportional_like_alpha_minus1 = 41
tau_BC:  proportional_like_alpha_minus1 = 41
```

主時間 `tau_ABC` からの alpha 差の最大値は次であった。

```text
tau_AB: 0.6588265529372437
tau_AC: 2.220446049250313e-16
tau_BC: 4.440892098500626e-16
```

ゲージ有効ケースの `alpha` 範囲は次であった。

```text
alpha_min: -1.0000000000000002
alpha_max: -0.9569628567496087
```

ゲージ有効ケースは、`R_C / R_A` ごとに次のように分布した。

```text
R_C / R_A = 1.0 : 19
R_C / R_A = 2.0 : 14
R_C / R_A = 4.0 : 8
```

一方、`R_C / R_A <= 0.5` では、分解能条件を満たさないケースが支配的であった。

```text
R_C が小さい場合、C セルが粗すぎて AB 位置変化量を読めない。
```

---

## 6. 読み

### 6.1 C ゲージには有効窓がある

本実験により、`C` は任意に軽くすればよいわけではないことが確認された。

`R_C` が小さすぎると、C の波長が長くなり、セル間隔が広がる。

この場合、AB の位置変化量が C のセル内に埋もれ、巻数・距離ゲージとして使えない。

一方、`R_C` が大きくなると分解能は上がるが、`C` 由来の関係性

```text
f_AC
f_BC
f_ABC
```

が強くなり、AB 主読出しを汚染する可能性が増える。

したがって、`C` には

```text
粗すぎず、強すぎない
```

有効窓がある。

### 6.2 最小 AB 調和モデルでは比例型が回収された

ゲージ有効ケースでは、距離指数候補はすべて

```text
alpha ≈ -1
```

であった。

ここで、本実験では

```text
I_AB(L) ∝ L^{-alpha}
```

としているため、`alpha=-1` は

```text
I_AB(L) ∝ L
```

を意味する。

つまり、C ゲージが有効な範囲では、既存 AB 調和モデルの比例型が回収された。

逆比例型 `alpha≈1`、逆二乗型 `alpha≈2` は出ていない。

### 6.3 関係性時間によって指数が変わる

念のため、四つの時間候補

```text
tau_AB
tau_AC
tau_BC
tau_ABC
```

で同じ `L_AB` 読出しを比較した。

主時間 `tau_ABC`、および C との関係時間 `tau_AC`, `tau_BC` では、指数分類は同じ比例型に保たれた。

一方、`tau_AB` で読むと、alpha が `tau_ABC` から最大で約 `0.659` ずれ、分類も `other` になった。

これは、距離指数が時間読出しの選び方に依存しうることを示す。

したがって、今後の ABC 距離指数実験では、

```text
何を時間として読むか
```

を明示しなければならない。

本仕様では、実験者が系の外からログを読む代表時間として、まず `tau_ABC` を主時間とする。

ただし、`tau_AB`, `tau_AC`, `tau_BC` は補助診断として記録し続ける。

### 6.4 これは逆二乗否定ではない

本予備実験は、最小 AB 調和モデルを C ゲージで読んだ場合の回収試験である。

したがって、この結果は

```text
自然界または閉鎖位相系一般に逆二乗型がない
```

ことを意味しない。

本実験で言えるのは、次である。

```text
現在の最小 AB 調和モデルを、C ゲージ適格条件を満たす範囲で読むと、
逆比例型でも逆二乗型でもなく、比例型として読まれる。
```

これは重要な否定対照である。

---

## 7. 判定

本予備実験は `valid` である。

ただし、成功の意味は逆二乗発見ではない。

成功したのは、次の分類である。

```text
C ゲージには有効窓がある。
その有効窓では、最小 AB 調和モデルは比例型として回収される。
逆比例型・逆二乗型は、本モデル範囲では未検出である。
```

---

## 8. 次の課題

次に行うべきことは、二つに分かれる。

### 8.1 C ゲージ適格性の精密化

`R_C` とセル幅の対応を、より公理的・読出し的に定義する必要がある。

今回の実装では、`R_C` が小さいほど C のセルが粗くなるという作業モデルを置いた。

これは自然な仮定だが、今後は、C の倍音、波長、セル幅、`R_C` の関係をより明確に定義する必要がある。

### 8.2 AB 物理側モデルの拡張

今回の AB 側は、前段の調和読出しを最小モデルとして使った。

したがって、比例型が回収されるのは自然である。

逆比例型または逆二乗型を物理側の候補として読むには、追加の閉鎖補償機構、時間位相伝搬、または関係性集合の再定義が必要かを検査する必要がある。

ただし、それを実装で作り込んではならない。

次の実験でも、まず

```text
何がネイティブに読めるのか
```

を優先する。
