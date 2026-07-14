# 白猫・黒猫・灰色猫準安定界面におけるAB-C-D段階実験仕様 v1

**副題:** AB二体準安定界面、C弱読出し、D強観測による白猫・黒猫・灰色猫分岐の閉鎖系数値実験  
**作成日:** 2026-07-14  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ／巨視的猫型準安定状態の段階実験仕様  
**対象実装:** Python 数値実験  
**Version DOI:** 10.5281/zenodo.21353209  
**Concept DOI:** 10.5281/zenodo.21353208

---

## 0. 実験の要約

本実験は、白猫状態 `A` と黒猫状態 `B` が等配分近傍にある状態を、直ちに「選択前の重ね合わせ」と仮定しない。

まず、AB二体系だけで次の三相を判別する。

```text
灰色猫固有相:
  A/B = 0.5/0.5 近傍に固定される。

灰色猫準安定相:
  A/B = 0.5/0.5 近傍で微小振動するが、一方へ自然収束しない。

自然選択相:
  CもDもないのに、AまたはBへ落ちる。
```

その後、Cを置いて準安定状態が読出し可能かを確認し、さらにDを置いて、観測によりABどちらかの状態へ落ちるかを検査する。

本実験の第一関門は、CやDの作用ではない。

```text
AB二体だけで、灰色猫固有相・灰色猫準安定相・自然選択相の界面を探すこと。
```

---

## 1. 目的

本実験の目的は、次の段階を順に確認することである。

1. AB二体だけで、灰色猫固有相・灰色猫準安定相・自然選択相を分類する。
2. Cを加えて、灰色猫固有相または準安定相を選択せずに読めるか確認する。
3. Dを加えて、観測によりA/Bの一方へ落ちる条件を調べる。
4. 灰色猫固有相では、Dを加えても灰色状態が保持されるか確認する。
5. A/Bが明確に分かれている大振幅領域では、Dの結果がC読出しと一致するか確認する。

ここで重要なのは、Dが常に収縮を起こすと仮定しないことである。

灰色猫がすでに固有状態として成立しているなら、Dで白猫または黒猫へ落ちる必要はない。

一方、A/Bが明確に分かれる領域では、Dによる結果は、Dが作った選択ではなく、すでに系内に現れているA/B分岐を読んだものとして扱う。

---

## 2. 状態変数

最低限、次の変数を記録する。

| 変数 | 意味 |
|---|---|
| `a` | 白猫Aの複素振幅 |
| `b` | 黒猫Bの複素振幅 |
| `p_A = |a|^2` | 白猫配分 |
| `p_B = |b|^2` | 黒猫配分 |
| `Q = p_A + p_B` | A/B総量 |
| `S = p_A - p_B` | A/B選択秩序変数 |
| `S_amp` | `S` の振動振幅 |
| `S_drift` | `S` の長時間ドリフト |
| `C_A, C_B` | CによるA/B読出し |
| `D_A, D_B` | DによるA/B読出し |
| `L_A, L_B` | A/Bの局在性 |
| `N_eff_A, N_eff_B` | A/Bの有効倍音次数 |
| `E_total` | 全保存量または保存レジストリ量 |

判定の中心は、

```text
Q が保存されるか。
S が固定されるか、微小振動するか、自然成長するか。
CとDがSの相分類を変えるか。
```

である。

---

## 3. AB二体準安定界面実験

## 3.1 目的

AB二体だけで、灰色猫がどのような状態として成立するかを調べる。

```text
C = off
D = off
```

この段階では、観測による選択を入れない。

## 3.2 初期状態

基準初期値は、

```math
a_0 = \frac{1}{\sqrt{2}}
```

```math
b_0 = \frac{e^{i\phi_0}}{\sqrt{2}}
```

である。

このとき、

```math
Q_0 = |a_0|^2 + |b_0|^2 = 1
```

```math
S_0 = |a_0|^2 - |b_0|^2 = 0
```

である。

微小な初期偏りを入れる場合は、

```math
p_{A,0} = \frac12 + s_0
```

```math
p_{B,0} = \frac12 - s_0
```

とする。

## 3.3 AB交換写像

最小のAB交換写像は、

```math
\begin{pmatrix}
a_{k+1}\\
b_{k+1}
\end{pmatrix}
=
U_\epsilon
\begin{pmatrix}
a_k\\
b_k
\end{pmatrix}
```

```math
U_\epsilon
=
\begin{pmatrix}
\cos\epsilon & i\sin\epsilon\\
i\sin\epsilon & \cos\epsilon
\end{pmatrix}
```

とする。

この写像だけでは、単純なユニタリ交換振動になる。

必要に応じて、閉鎖系内の安定性を調べるため、弱い復元項または弱い非線形項を追加する。

ただし、Cなし・Dなしで一方へ自然収束する項は、猫型選択実験には使わない。

## 3.4 掃引パラメータ

最低限、次を掃引する。

| パラメータ | 意味 |
|---|---|
| `epsilon` | AB交換振動の強さ |
| `phi_0` | 初期相対位相 |
| `s_0` | 初期A/B微小偏り |
| `noise_amp` | 数値揺らぎまたは微小外乱 |
| `K` | 発展ステップ数 |

調和振動として灰色猫準安定状態を扱う場合、`S` の振幅は小さく保つ。

目標値は、

```text
S_amp < S_gray_limit
```

である。

初期値としては、

```text
S_gray_limit = 0.05
S_amp target = 0.01 から 0.03 程度
```

を使う。

## 3.5 相分類

AB二体実験の結果は、次の三相に分類する。

### 灰色猫固有相

```text
S_mean ≈ 0
S_amp -> 0
S_drift ≈ 0
Q ≈ 1
```

この場合、A/B等配分状態は、新しい灰色猫固有状態として扱う。

### 灰色猫準安定相

```text
S_mean ≈ 0
0 < S_amp < S_gray_limit
S_drift ≈ 0
Q ≈ 1
```

この相が、Dによる選択実験の本命領域である。

### 自然選択相

```text
abs(S) grows without C or D
S -> +1 or S -> -1
```

この相では、観測による選択を主張できない。

---

## 4. C読出し実験

## 4.1 目的

AB二体で得られた灰色猫固有相または灰色猫準安定相に、Cを加える。

```text
C = on
D = off
```

Cの目的は、A/B配分を読むことである。

この段階では、CがA/Bの一方へ落としてはならない。

## 4.2 C読出し量

Cは次を読む。

```math
C_A \approx p_A
```

```math
C_B \approx p_B
```

誤差は、

```math
E_C = |C_A-p_A| + |C_B-p_B|
```

で記録する。

## 4.3 C合格条件

灰色猫準安定相で、

```text
E_C < tol_C
abs(S_afterC - S_beforeC) < epsilon_C
no selection
```

を満たすこと。

灰色猫固有相では、

```text
C_A ≈ 0.5
C_B ≈ 0.5
S remains near 0
```

を満たすこと。

大振幅領域では、Cはその時点のA/B偏りを読む。

```text
C_A > C_B if S > 0
C_B > C_A if S < 0
```

を確認する。

---

## 5. D観測実験

## 5.1 目的

次にDを加える。

```text
C = optional
D = on
```

Dは、A/Bの一方へ落ちる観測が起こるかを調べるための強い観測写像である。

ただし、Dによる結果は状態相ごとに分類する。

## 5.2 灰色猫固有相でのD

灰色猫固有相では、Dを加えても灰色状態が保持されるかを調べる。

期待される判定は、

```text
D_A ≈ 0.5
D_B ≈ 0.5
S remains near 0
```

である。

この場合、Dで白猫または黒猫へ落ちないことが結果である。

## 5.3 灰色猫準安定相でのD

灰色猫準安定相では、DがA/Bの一方へ落とすかを調べる。

判定は、

```text
D_A ≈ 1, D_B ≈ 0
```

または、

```text
D_A ≈ 0, D_B ≈ 1
```

である。

ただし、Dなしで同じ落下が起こる場合は、Dによる選択とは判定しない。

## 5.4 大振幅分離領域でのD

ABの振幅が大きく、A/Bが明確に分かれている領域では、Dはその時点のA/B偏りを読むだけでよい。

この場合、Dの読出しはCの読出しと一致するはずである。

```text
sign(D_A - D_B) = sign(C_A - C_B)
```

この一致が得られる場合、Dが選択を作ったのではなく、すでに系内に現れているA/B偏りを読んだと判定する。

---

## 6. 四段階本実験

本実験は、次の四段階で実施する。

## Stage 1: AB準安定界面探索

```text
C = off
D = off
```

出力：

```text
gray eigen phase
gray metastable phase
natural selection phase
```

## Stage 2: C読出し確認

```text
C = on
D = off
```

出力：

```text
C_read_error
C_induced_bias
C_selection_triggered
```

## Stage 3: D観測確認

```text
C = off or recorded
D = on
```

出力：

```text
D_result
D_induced_selection
D_vs_C_agreement
```

## Stage 4: C-D比較

同じAB初期条件に対して、C読出しとD観測を比較する。

確認事項：

```text
灰色猫固有相:
  CもDも 0.5/0.5 を読む。

灰色猫準安定相:
  Cは読むが選ばない。
  Dが一方選択を起こすかを判定する。

大振幅分離領域:
  CとDは同じA/B偏りを読む。
  観測が作った選択とは扱わない。
```

---

## 7. 合格条件

本実験の最低合格条件は次である。

1. AB二体だけで、灰色猫固有相・灰色猫準安定相・自然選択相を分類できる。
2. 灰色猫準安定相では、CがA/B配分を読めるが、一方選択を起こさない。
3. 灰色猫固有相では、Dを加えても白猫または黒猫へ落ちない。
4. 灰色猫準安定相では、Dによる選択が起こるか、起こらないかをDなし対照と比較して判定できる。
5. 大振幅分離領域では、D結果がC読出しと一致する。
6. `Q` または全保存レジストリが許容範囲内で保存される。

---

## 8. 失敗条件

次の場合は、後段のC/D実験へ進まない。

1. AB二体だけで全条件が自然選択相へ落ちる。
2. `Q` が保存されない。
3. Cが常に選択を起こす。
4. CがA/B配分を読めない。
5. Dなし対照とDあり条件の差が分類できない。

---

## 9. 実験結果欄

## 9.1 Stage 1: AB準安定界面探索

出力:

```text
gray_cat_ab_metastable_interface_preliminary_result_v1/
```

結果:

```text
total_cases = 5600
gray_eigen = 1248
gray_metastable = 733
large_oscillation = 470
natural_selection = 1070
unstable_or_drifting = 2079
```

AB二体だけで、灰色猫固有相、灰色猫準安定相、大振幅振動相、自然選択相を分離できた。

## 9.2 Stage 2: C読出し窓確認

出力:

```text
gray_cat_c_readout_window_preliminary_result_v1/
```

結果:

```text
total_cases = 1080
C_window_count = 247
C_informative_window_count = 144
C_nonzero_backaction_window_count = 114
```

灰色猫準安定候補に対して、A/B配分を読めるが一方選択を起こさないC読出し窓を確認した。

## 9.3 Stage 3: D観測応答

出力:

```text
gray_cat_d_observation_response_preliminary_result_v1/
```

結果:

```text
total_cases = 7056
D_induced_selection_count = 2016
gray_kept_eigen = 1062
white_selected = 1244
black_selected = 772
```

同じD開始状態からDなし対照を並走させ、Dなしでは白猫または黒猫へ落ちない条件で、Dありの場合に選択が起こる領域を確認した。

灰色猫固有相では、強D条件でも灰色状態が保持された。

大振幅分離領域では、D結果がC読出し符号と一致する条件を確認した。

## 9.4 Stage 4: D選択境界

出力:

```text
gray_cat_d_selection_boundary_preliminary_result_v1/
```

結果:

```text
boundary_points = 438
selection_possible_boundary_points = 438
no_selection_boundary_points = 0
min_D_gain_overall = 0.0225
max_min_D_gain_overall = 0.065
```

灰色猫準安定候補に対して、観測開始ステップとD利得の境界表を作成した。

今回の掃引範囲では、対象とした全ての観測開始ステップでD起因選択が可能だった。

候補ごとの最小D利得は二段に分かれた。

```text
weak threshold candidate:
  min_D_gain = 0.0225

strong threshold candidates:
  min_D_gain = 0.065
```
