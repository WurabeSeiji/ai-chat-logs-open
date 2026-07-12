# ABC閉鎖位相系における独立計量Cの配置対照距離指数予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 目的

本メモでは、前回の

[ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し予備実験検証メモ v1.md](ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し予備実験検証メモ_v1.md)

に続き、独立計量波 `C` の配置を変えた対照試験を記録する。

前回の予備実験では、AB 二体だけでは距離指数

```text
L_AB^0
1/L_AB
1/L_AB^2
```

を独立に計量するゲージが足りないため、第三波 `C` を導入した。

ただし、`C` は単なる外部観測器ではない。

`C` 自身も閉鎖位相系の一部であり、`f_AC`, `f_BC`, `f_ABC` という関係性を持つ。

したがって、`C` の配置を変えたときに

```text
距離指数が変わるのか
C 由来の汚染としてゲージ有効窓が狭くなるだけなのか
```

を分けて調べる必要がある。

本実験では、逆比例項・逆二乗項は実装しない。

目的は、C 配置の違いがネイティブな距離指数読出しを作るのか、それともゲージ条件の破綻として現れるのかを確認することである。

---

## 2. 実行コマンド

```text
python3 run_abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_v1.py
```

出力:

```text
abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [abc_c_gauge_c_position_control_preliminary_report_v1.md](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_control_preliminary_report_v1.md) |
| JSON | [abc_c_gauge_c_position_control_preliminary_result_v1.json](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_control_preliminary_result_v1.json) |
| 全設定 CSV | [abc_c_gauge_c_position_control_configs_v1.csv](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_control_configs_v1.csv) |
| ケース CSV | [abc_c_gauge_c_position_control_cases_v1.csv](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_control_cases_v1.csv) |
| A側/B側対称性 CSV | [abc_c_gauge_c_position_control_pair_symmetry_v1.csv](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_control_pair_symmetry_v1.csv) |
| 有効数図 | [abc_c_gauge_c_position_valid_counts_v1.png](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_valid_counts_v1.png) |
| alpha 候補図 | [abc_c_gauge_c_position_alpha_candidates_v1.png](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_alpha_candidates_v1.png) |
| A側/B側対称性図 | [abc_c_gauge_c_position_pair_symmetry_v1.png](abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_pair_symmetry_v1.png) |

---

## 4. 実験設計

### 4.1 C 配置モード

本実験では、次の C 配置を比較した。

```text
symmetric
symmetric_pi_flip
a_side_small
b_side_small
a_side_large
b_side_large
a_side_large_pi_flip
b_side_large_pi_flip
```

`symmetric` は C を AB に対して対称に置く作業モデルである。

`a_side_*`, `b_side_*` は、C を A 側または B 側へ寄せた場合の対照である。

`pi_flip` は、C の位相を反転させた対照である。

### 4.2 禁止したこと

本実験でも、次は実装していない。

```text
1/L_AB
1/L_AB^2
1/A_chi_tau
F = G m_A m_B / L_AB^2
```

したがって、逆比例型または逆二乗型が出る場合は、C 配置そのものから読めなければならない。

### 4.3 判定条件

前回と同じく、距離指数判定には次を満たすケースのみを使う。

```text
resolution_valid
clock_valid
disturbance_valid
closure_valid
```

ここで重要なのは、C 配置が非対称になった場合、その効果を距離指数として読む前に、まずゲージ汚染として除外される可能性があることである。

---

## 5. 結果サマリー

主要結果は次である。

```text
c_position_control_preliminary_valid: true
config_count: 160
gauge_valid_count: 70
gauge_valid_position_pair_count: 23
max_pair_abs_alpha_difference: 0.11776992863447344
```

C 配置ごとのゲージ有効数は次であった。

```text
a_side_large: 6
a_side_large_pi_flip: 6
a_side_small: 11
b_side_large: 6
b_side_large_pi_flip: 6
b_side_small: 11
symmetric: 12
symmetric_pi_flip: 12
```

主時間 `tau_ABC` での alpha 分類は次であった。

```text
proportional_like_alpha_minus1: 70
```

逆比例型 `alpha≈1`、逆二乗型 `alpha≈2` は、ゲージ有効ケースでは出ていない。

---

## 6. 読み

### 6.1 C 配置は逆数・逆二乗を作らなかった

本実験では、C の対称配置、位相反転、A 側配置、B 側配置を比較した。

しかし、ゲージ有効ケースの `tau_ABC` 読出しはすべて

```text
alpha ≈ -1
```

に分類された。

ここで、前回と同じく

```text
I_AB(L) ∝ L^{-alpha}
```

としているため、`alpha≈-1` は比例型

```text
I_AB(L) ∝ L
```

を意味する。

したがって、C の配置を変えても、現在の最小 AB 調和モデルから逆比例型または逆二乗型がネイティブに立ち上がる結果は得られなかった。

### 6.2 大きな非対称 C は指数則ではなく汚染として現れた

大きな A 側配置または B 側配置では、ゲージ有効数が減った。

```text
symmetric: 12
a_side_small / b_side_small: 11
a_side_large / b_side_large: 6
```

これは、C を非対称に置くと新しい安定な距離指数が出るというより、`f_AC`, `f_BC` の偏りが C ゲージの有効条件を狭めることを示している。

つまり、C 配置の非対称性は、今の段階では

```text
距離指数を作る機構
```

ではなく、

```text
距離指数を読むためのゲージ汚染
```

として扱う方が安全である。

### 6.3 A側/B側の対称性は概ね保たれた

A 側配置と B 側配置のペア比較では、ゲージ有効ペアの最大 alpha 差は次であった。

```text
max_pair_abs_alpha_difference: 0.11776992863447344
```

完全一致ではないが、A 側と B 側を入れ替えても、距離指数の主分類は変わらなかった。

このため、本実験範囲では、A 側配置だけが特別な指数則を生む、または B 側配置だけが特別な指数則を生む、という非対称な結果は見られない。

### 6.4 fabc は主時間として使うが、円周方向の原因ではない

本実験では、外部ログの代表時間として `tau_ABC` を主時間にした。

これは、実験者が ABC 系全体を外側から読む場合、全体系の関係性 `f_ABC` が代表時間候補になるためである。

ただし、`f_ABC` をそのまま AB 円周方向の原因として読むわけではない。

円周方向の関係性候補は、まず

```text
f_AB
f_AC
f_BC
```

として分けて扱う必要がある。

したがって、今後の ABC 実験では、

```text
時間ゲージとしての tau_ABC
円周方向補償候補としての f_AB, f_AC, f_BC
```

を混同しない。

---

## 7. 判定

本予備実験は `valid` である。

ただし、成功の意味は、逆二乗読出しの発見ではない。

成功したのは、次の否定対照である。

```text
C 配置を対称・非対称・位相反転で変えても、
ゲージ有効ケースでは比例型 alpha≈-1 が保たれた。

大きな非対称 C 配置は、安定な逆数・逆二乗則ではなく、
主にゲージ汚染として現れた。
```

したがって、距離指数を読むためには、C 配置をいじるだけでは不十分である。

次に必要なのは、C を独立計量ゲージとして使いながら、

```text
AB 関係性 f_AB
AC 関係性 f_AC
BC 関係性 f_BC
全体系時間 tau_ABC
```

を分離して記録する実験である。

---

## 8. 次の課題

### 8.1 C を重くしすぎても軽くしすぎてもいけない

`C` は位置変化量を読むためのセル分解能を持つ必要がある。

しかし、`C` が強すぎると `f_AC`, `f_BC` が AB 関係性を汚染する。

したがって、今後の実験では、C の条件を

```text
十分細かいが、主関係性を支配しない
```

範囲に置く必要がある。

### 8.2 三つの円周方向関係性を明示的に分ける

ABC 三体では、関係性は少なくとも次の四つに分かれる。

```text
f_ABC
f_AB
f_AC
f_BC
```

このうち、全体系の代表時間候補は `f_ABC` である。

一方、AB の円周方向の変位を読むときは、`f_AB`, `f_AC`, `f_BC` の寄与をベクトル的に分ける必要がある。

この分離をしないと、距離指数の読みが C ゲージ汚染と混ざる。

### 8.3 逆二乗の有無はまだ未決である

本実験は、逆二乗を否定したのではない。

否定されたのは、次の狭い主張である。

```text
現在の最小 AB 調和モデルと C 配置操作だけで、
ゲージ有効範囲に安定な逆二乗読出しが現れる。
```

この主張は支持されなかった。

逆二乗型が読めるかどうかは、次の段階で、関係性分離と時間ゲージ分離を入れて検査する。

