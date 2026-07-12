# ABC閉鎖位相系における関係性補償分解距離指数予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 目的

本メモでは、

[ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し実験仕様書 v1.md](ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し実験仕様書%20v1.md)

に基づき、ABC 三体系の関係性補償を

```text
f_AB
f_AC
f_BC
f_ABC
```

へ分解して記録する予備実験を行った。

前回の C 配置対照では、C の配置を変えても、ゲージ有効ケースでは比例型 `alpha≈-1` が保たれた。

しかし、ABC 三体系では、`C` を入れた時点で `f_AC`, `f_BC`, `f_ABC` が生成される。

したがって、次を分けずに距離指数を読むと、C ゲージ汚染と AB 固有の関係性補償が混ざる。

```text
AB 関係性: f_AB
C 由来の二体関係性: f_AC, f_BC
全体系の共通関係性: f_ABC
```

本実験の目的は、逆比例型または逆二乗型を探しに行くことではない。

目的は、まず

```text
f_ABC を代表時間として使うこと
f_ABC を円周方向へ直接足さないこと
f_AB, f_AC, f_BC を円周方向の射影候補として分けること
```

の数値上の意味を確認することである。

---

## 2. 実行コマンド

```text
python3 run_abc_c_gauge_relation_decomposition_preliminary_v1.py
```

出力:

```text
abc_c_gauge_relation_decomposition_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [abc_c_gauge_relation_decomposition_preliminary_report_v1.md](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_preliminary_report_v1.md) |
| JSON | [abc_c_gauge_relation_decomposition_preliminary_result_v1.json](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_preliminary_result_v1.json) |
| 設定 CSV | [abc_c_gauge_relation_decomposition_configs_v1.csv](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_configs_v1.csv) |
| ケース CSV | [abc_c_gauge_relation_decomposition_cases_v1.csv](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_cases_v1.csv) |
| 有効数図 | [abc_c_gauge_relation_decomposition_valid_counts_v1.png](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_valid_counts_v1.png) |
| ペア alpha 図 | [abc_c_gauge_relation_decomposition_pair_alpha_v1.png](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_pair_alpha_v1.png) |
| native/pair 比較図 | [abc_c_gauge_relation_decomposition_native_vs_pair_alpha_v1.png](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_native_vs_pair_alpha_v1.png) |
| C 汚染診断図 | [abc_c_gauge_relation_decomposition_contamination_v1.png](abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_contamination_v1.png) |

---

## 4. 実験設計

### 4.1 関係性補償

本実験では、次の四つの関係性を分けて記録した。

```text
f_AB
f_AC
f_BC
f_ABC
```

このうち、主対象は `f_AB` である。

`f_AC`, `f_BC` は、C 由来の二体関係性である。

`f_ABC` は、ABC 全体系の共通関係性である。

### 4.2 円周方向射影

円周方向に直接射影する候補は、二体関係性

```text
f_AB
f_AC
f_BC
```

である。

各局所波の円周方向候補を、次の作業表示で読む。

```text
a_A^circ = f_AB + sigma_AC f_AC
a_B^circ = f_AB + sigma_BC f_BC
```

`sigma_AC`, `sigma_BC` は射影符号である。

本実験では、次の射影モードを比較した。

```text
c_opposes_ab
c_assists_ab
a_opposes_b_assists
a_assists_b_opposes
```

### 4.3 `f_ABC` の扱い

`f_ABC` は代表時間候補として使う。

```text
tau_ABC
```

ただし、`f_ABC` を AB 円周方向へ直接足さない。

直接足した場合は、主結果ではなく、誤注入対照として記録する。

この区別は重要である。

`f_ABC` は存在しないのではない。

しかし、AB 相対円周方向の読出しでは、共通モードとして見えにくい、または別軸に属する可能性がある。

### 4.4 禁止したこと

本実験でも、次は実装していない。

```text
1/L_AB
1/L_AB^2
1/A_chi_tau
F = G m_A m_B / L_AB^2
```

したがって、逆比例型または逆二乗型が出る場合は、関係性補償の分解読出しから出なければならない。

---

## 5. 結果サマリー

主要結果は次である。

```text
relation_decomposition_preliminary_valid: true
config_count: 240
decomposition_valid_count: 180
AB_dominant_valid_count: 48
non_AB_dominant_count: 132
inverse_or_inverse_square_in_AB_dominant_count: 0
```

分解有効ケースの射影モード分布は次であった。

```text
a_assists_b_opposes: 45
a_opposes_b_assists: 45
c_assists_ab: 45
c_opposes_ab: 45
```

AB 主導ケースでの `tau_ABC` 距離指数分類は次であった。

```text
proportional_like_alpha_minus1: 48
```

全分解有効ケースでの `tau_ABC` 距離指数分類は次であった。

```text
constant_like_alpha0: 3
other: 46
proportional_like_alpha_minus1: 131
```

native な `f_AB` の分類は次であった。

```text
proportional_like_alpha_minus1: 180
```

C バイアス項の分類は次であった。

```text
constant_like_alpha0: 150
unfit: 30
```

C 非対称項の分類は次であった。

```text
constant_like_alpha0: 150
unfit: 30
```

`f_ABC` を円周方向へ直接足した誤注入対照は、次であった。

```text
other: 63
proportional_like_alpha_minus1: 117
```

関係性時間の差は、AB 主導ケースで次であった。

```text
max |tau_AB - tau_ABC| alpha delta: 0.6588265529372439
max |tau_AC - tau_ABC| alpha delta: 6.661338147750939e-16
max |tau_BC - tau_ABC| alpha delta: 6.661338147750939e-16
```

---

## 6. 読み

### 6.1 AB 主導ケースでは比例型が保たれた

AB 主導と判定できる 48 ケースでは、距離指数候補はすべて

```text
alpha ≈ -1
```

であった。

ここで、これまでと同じく

```text
I_AB(L) ∝ L^{-alpha}
```

としているため、`alpha≈-1` は比例型

```text
I_AB(L) ∝ L
```

を意味する。

つまり、`f_AB`, `f_AC`, `f_BC` を分離しても、AB 主導ケースでは逆比例型や逆二乗型は出ていない。

### 6.2 C バイアスと C 非対称は定数型として現れた

C バイアス項と C 非対称項は、分解有効ケースで主に

```text
constant_like_alpha0
```

として分類された。

これは自然である。

本実験では、`f_AC`, `f_BC` を距離依存項として実装していない。

したがって、C 由来の射影バイアスは、AB 距離 `L_AB` に対してほぼ定数として読まれる。

この結果は、C 由来の項が逆比例型や逆二乗型に見えるわけではないことを示す。

### 6.3 `f_ABC` の直足しは分類を濁らせる

`f_ABC` を円周方向へ直接足した誤注入対照では、

```text
other: 63
proportional_like_alpha_minus1: 117
```

となった。

これは、`f_ABC` を AB 円周方向の項として足すと、読出し分類が濁ることを示す。

したがって、現段階では

```text
f_ABC は代表時間 tau_ABC として使う
f_ABC は AB 円周方向へ直接足さない
```

という分離が妥当である。

これは `f_ABC` を無視するという意味ではない。

`f_ABC` は、全体系の共通モードまたは中心方向補償として記録する。

ただし、AB 相対円周方向の原因項としては、まず `f_AB`, `f_AC`, `f_BC` を読む。

### 6.4 `tau_AB` は主時間と一致しない

AB 主導ケースでも、`tau_AB` で読んだ alpha は、`tau_ABC` から最大で

```text
0.6588265529372439
```

ずれた。

一方、`tau_AC`, `tau_BC` は、この予備モデルでは `tau_ABC` とほぼ一致した。

これは、距離指数が

```text
どの時間位相で読むか
```

に依存しうることを再確認している。

本実験の主時間は、外部ログの代表時間としての `tau_ABC` である。

ただし、`tau_AB`, `tau_AC`, `tau_BC` は今後も補助診断として記録する。

---

## 7. 判定

本予備実験は `valid` である。

ただし、成功の意味は、逆二乗読出しの発見ではない。

成功したのは、次の分離である。

```text
f_AB は AB 主導ケースで比例型として読まれる。
f_AC, f_BC の射影バイアスは主に定数型として読まれる。
f_ABC を円周方向へ直足しすると分類が濁る。
```

したがって、次の作業規則が支持された。

```text
代表時間: tau_ABC
円周方向候補: f_AB, f_AC, f_BC
共通モード・中心方向候補: f_ABC
```

この分離を行わないと、距離指数の読出しは C ゲージ汚染と混ざる。

---

## 8. 次の課題

### 8.1 この結果は逆二乗否定ではない

本実験で否定されたのは、次の狭い主張である。

```text
現在の最小関係性分解モデルにおいて、
AB 主導ケースから逆比例型または逆二乗型がネイティブに出る。
```

この主張は支持されなかった。

ただし、これは閉鎖位相系一般に逆二乗型がないことを意味しない。

### 8.2 距離指数を読むには追加の計量構造が必要である

今回の結果では、AB 主導ケースは比例型、C 由来項は定数型に分かれた。

したがって、逆比例型または逆二乗型を読むには、単なる C 配置や関係性分解では不足している可能性が高い。

候補は次である。

```text
時間位相伝搬
面積スイープの計量化
球殻または面積ゲージ
C の読出し窓の再定義
```

ただし、これらを実装で作り込んではならない。

次の実験でも、まず

```text
何がネイティブに読めるのか
何が観測側の作り込みなのか
```

を分けて検査する。

