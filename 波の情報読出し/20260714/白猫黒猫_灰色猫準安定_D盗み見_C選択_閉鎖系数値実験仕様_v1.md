# 白猫・黒猫準安定共存状態に対する盗み見観測D・選択観測C・再確認実験仕様 v1

**副題:** 灰色猫状態の準安定保持、弱観測による等配分確認、フェルミオン的選択写像による一状態選択の閉鎖系数値実験  
**作成日:** 2026-07-14  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ／巨視的シュレーディンガー猫型思考実験・数値実験仕様  
**対象実装:** Claude Code  
**Version DOI:** pending  
**Concept DOI:** pending

---

## 0. 実験の要約

本実験は、白猫状態 `A` と黒猫状態 `B` が等配分で準安定に共存する状態を「灰色猫」と定義し、次の三段階を閉じた数値系で実証する。

```text
1. 観測Dにより、選択前の状態が A=0.5, B=0.5 であることを弱く確認する。
2. 観測Cとの相互作用により、AまたはBの一方だけが残る選択写像を発生させる。
3. 再度Dで読み、A=1, B=0 または A=0, B=1 になったことを確認する。
```

本実験では、観測と相互作用を区別しない。  
観測Cおよび観測Dは、いずれも系内に存在する波との干渉写像である。

本実験の核心は次である。

```text
灰色猫は「第三の中間状態」ではない。
白猫写像Aと黒猫写像Bが等配分で準安定に共存する複合状態である。
```

---

# 1. 目的

本実験の目的は、以下を同時に満たす閉鎖系写像を構成することである。

1. `A` と `B` が等配分で準安定に共存する。
2. 観測Dは、A/Bを選択せずに `A≈0.5, B≈0.5` と読み出せる。
3. 観測Cは、A/Bを識別できる名前写像を読み取る。
4. 観測Cとの相互作用後、A/Bの一方だけが残る。
5. 選択後、観測Cの内部状態に白猫／黒猫の名前が記録される。
6. 選択後に観測Dで再確認すると、`(1,0)` または `(0,1)` と読める。
7. 観測Dだけでは選択が起きない。
8. ボゾン的線形写像では、同じ条件でも一方選択が完成しない。
9. フェルミオン的非線形写像でのみ、微小な偏りが一状態選択へ増幅される。
10. 全系の保存量は追跡され、非選択状態の量が単純消滅したように実装しない。

---

# 2. 用語定義

## 2.1 白猫状態 A

白猫状態を表す波形・チャネル・名前写像。

```text
A = white cat state
```

Aは次の構造を持つ。

```text
A_wave
A_localization
A_name = white
A_weight
```

---

## 2.2 黒猫状態 B

黒猫状態を表す波形・チャネル・名前写像。

```text
B = black cat state
```

Bは次の構造を持つ。

```text
B_wave
B_localization
B_name = black
B_weight
```

---

## 2.3 灰色猫状態

灰色猫は第三の独立状態ではない。

```math
\Psi_{\mathrm{gray}}
=
a\Psi_A+b\Psi_B
```

初期条件は、

```math
|a|^2 \approx |b|^2 \approx \frac12
```

とする。

ただし、A/Bは静止させず、選択前に微小な交換振動を許す。

```math
|a_k|^2
=
\frac12+\epsilon\cos(\Omega k)
```

```math
|b_k|^2
=
\frac12-\epsilon\cos(\Omega k)
```

ここで、

```text
epsilon << 0.5
```

とする。

灰色猫とは、

```text
白猫写像Aと黒猫写像Bが、
等配分近傍で微小振動しながら共存する準安定複合状態
```

である。

---

## 2.4 観測D：盗み見写像

Dは、選択前のA/B等配分を読み出す弱観測写像である。

Dは次を読み出す。

```text
D_A ≈ 0.5
D_B ≈ 0.5
```

ただし、DはA/Bの一方を選択してはならない。

Dは完全な非破壊測定ではない。  
D自身も相互作用であるため、ABをわずかに変化させてもよい。

許容条件は、

```math
|\Delta_D(A-B)| < \varepsilon_D
```

である。

---

## 2.5 観測C：選択・識別写像

Cは、A/Bとの相互作用により、AまたはBの一方を選択し、同時にその名前を内部記録する写像である。

Cは次の二つの記録チャネルを持つ。

```text
C_white
C_black
```

選択後は、

```text
white selected:
C_white ≈ 1
C_black ≈ 0
```

または、

```text
black selected:
C_white ≈ 0
C_black ≈ 1
```

となる。

---

# 3. 実験全体の状態空間

最低限、状態は次を持つ。

```text
state = {
    A_wave,
    B_wave,
    A_weight,
    B_weight,
    A_name_vector,
    B_name_vector,
    C_state,
    D_state,
    residual_state,
    conserved_registry
}
```

推奨する主要変数は次である。

| 変数 | 意味 |
|---|---|
| `a` | 白猫Aの複素振幅係数 |
| `b` | 黒猫Bの複素振幅係数 |
| `p_A = |a|^2` | 白猫配分 |
| `p_B = |b|^2` | 黒猫配分 |
| `S = p_A - p_B` | 選択秩序変数 |
| `Q = p_A + p_B` | A/B総量 |
| `C_white` | Cに記録された白猫成分 |
| `C_black` | Cに記録された黒猫成分 |
| `D_white` | Dが読んだ白猫成分 |
| `D_black` | Dが読んだ黒猫成分 |
| `L_A, L_B` | 局在性指標 |
| `N_eff_A, N_eff_B` | 有効倍音次数 |
| `R_residual` | 非選択側から移った残余量 |
| `E_total` | 全保存量 |

---

# 4. 名前写像

白猫と黒猫は、局在性だけではなく、Cが識別可能な名前写像を持つ。

名前ベクトルは、相互に識別可能かつ規格化可能な有限複素配列とする。

例：

```math
h_{\mathrm{white}}
=
\frac{1}{2}
(1,1,1,1)
```

```math
h_{\mathrm{black}}
=
\frac{1}{2}
(1,-1,1,-1)
```

または複素位相を使って、

```math
h_{\mathrm{white}}
=
\frac{1}{2}
(1,i,-1,-i)
```

```math
h_{\mathrm{black}}
=
\frac{1}{2}
(1,-i,-1,i)
```

とする。

必要条件：

```math
\langle h_{\mathrm{white}},h_{\mathrm{black}}\rangle
\approx 0
```

ただし完全直交でなくてもよい。  
識別誤差を掃引できるようにする。

Cの名前読出しは、

```math
M_{C,\mathrm{white}}
=
|\langle C,h_{\mathrm{white}}\rangle|^2
```

```math
M_{C,\mathrm{black}}
=
|\langle C,h_{\mathrm{black}}\rangle|^2
```

で計算する。

---

# 5. 初期状態

## 5.1 A/B準安定初期値

初期値は、

```math
a_0
=
\frac{1}{\sqrt2}
```

```math
b_0
=
\frac{e^{i\phi_0}}{\sqrt2}
```

とする。

初期の選択量は、

```math
S_0
=
|a_0|^2-|b_0|^2
=
0
```

である。

初期の総量は、

```math
Q_0
=
|a_0|^2+|b_0|^2
=
1
```

とする。

---

## 5.2 微小交換振動

CおよびDがない場合、A/Bは次の弱い交換写像で準安定に保つ。

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

ここで、

```text
epsilon = small
```

とする。

期待される挙動：

```text
p_A ≈ 0.5 ± small oscillation
p_B ≈ 0.5 ∓ small oscillation
S ≈ 0
```

CなしでA/Bの一方へ収束してはならない。

---

# 6. 観測D：盗み見写像

## 6.1 Dの目的

Dは、C作用前に、

```text
A ≈ 0.5
B ≈ 0.5
```

を読み出す。

またC作用後に、

```text
A ≈ 1, B ≈ 0
```

または、

```text
A ≈ 0, B ≈ 1
```

を再確認する。

---

## 6.2 Dの読み出し量

Dは、A/Bの総量と差分を読む。

```math
Q_D
=
D_A+D_B
```

```math
S_D
=
D_A-D_B
```

復元値は、

```math
P_A^{(D)}
=
\frac{Q_D+S_D}{2}
```

```math
P_B^{(D)}
=
\frac{Q_D-S_D}{2}
```

とする。

選択前の理想値：

```text
Q_D ≈ 1
S_D ≈ 0
P_A^(D) ≈ 0.5
P_B^(D) ≈ 0.5
```

---

## 6.3 Dの弱結合

Dの結合強度を `g_D` とする。

```text
g_D_min < g_D < g_D_select
```

である必要がある。

- `g_D < g_D_min`  
  読み出し不能
- `g_D > g_D_select`  
  D自身が選択を起こす
- 中間域  
  等配分を読めるが選択を起こさない

Dによる状態更新は、例として、

```math
a'
=
a\left(1+\eta_D \xi_D\right)
```

```math
b'
=
b\left(1-\eta_D \xi_D\right)
```

とする。

ただし、

```math
|\eta_D \xi_D| \ll 1
```

とする。

Dの作用前後で、

```math
|S_{\mathrm{afterD}}-S_{\mathrm{beforeD}}|
<
\varepsilon_D
```

を満たすこと。

---

## 6.4 Dの合格条件

Dは以下を満たしたとき合格とする。

```text
before C:
abs(D_A - 0.5) < tol_D
abs(D_B - 0.5) < tol_D
abs(S_afterD - S_beforeD) < epsilon_D
```

さらにD単独で十分長時間再帰しても、

```text
abs(S) < S_stable_limit
```

を保つこと。

---

# 7. 観測C：選択・名前読出し写像

## 7.1 Cの揺らぎ

Cは、わずかな `cos^2` 揺らぎを持つ。

```math
w_A(\theta_C)
=
\cos^2\theta_C
```

```math
w_B(\theta_C)
=
\sin^2\theta_C
```

したがって、

```math
w_A+w_B=1
```

である。

偏りは、

```math
\delta_C
=
\eta_C
\left(
\cos^2\theta_C-\frac12
\right)
```

とする。

ここで `eta_C` は小さい。

---

## 7.2 Cの初期作用

CはABへ微小な非対称を与える。

```math
a'
=
a(1+\delta_C)
```

```math
b'
=
b(1-\delta_C)
```

ただし、単純正規化のみで選択を人工生成しないこと。

Cを含む全保存量を追跡する。

---

## 7.3 フェルミオン的増幅写像

微小偏りを一状態選択へ増幅する。

秩序変数

```math
S_k
=
p_{A,k}-p_{B,k}
```

に対し、最小モデルとして、

```math
S_{k+1}
=
S_k
+
\mu_C
S_k
(1-S_k^2)
+
\delta_C
```

を候補とする。

ただし、この式だけではCなしでも数値誤差で分岐する可能性がある。

より安全な形として、

```math
S_{k+1}
=
S_k
+
\left(
\mu I_C-\gamma
\right)
S_k
(1-S_k^2)
+
\delta_C
```

を使う。

条件：

```text
Cなし:
mu * I_C - gamma < 0
```

```text
Cあり:
mu * I_C - gamma > 0
```

Cは法則を切り替えるのではなく、干渉強度 `I_C` により実効安定性を変える。

---

## 7.4 波形レベルの更新

重みだけでなく、波形も追跡する。

```math
\Psi_A'
=
F_A(
\Psi_A,
\Psi_B,
C,
r,
t,
\theta_C
)
```

```math
\Psi_B'
=
F_B(
\Psi_A,
\Psi_B,
C,
r,
t,
\theta_C
)
```

交換干渉散乱行列を使用する場合、

```math
A_{k+1}
=
\operatorname{norm}_{\mathrm{registry}}
\left(
rA_k+tB_k+\epsilon_C C_{\mathrm{white}}
\right)
```

```math
B_{k+1}
=
\operatorname{norm}_{\mathrm{registry}}
\left(
tA_k+rB_k+\epsilon_C C_{\mathrm{black}}
\right)
```

とする。

ただし、`norm_registry` は単純な各チャネル独立正規化ではなく、全保存レジストリを維持する処理とする。

---

## 7.5 Cの名前記録

選択後、Cの内部記録は、

```math
C_{\mathrm{white}}
\propto
|\langle C,h_{\mathrm{white}}\rangle|^2
```

```math
C_{\mathrm{black}}
\propto
|\langle C,h_{\mathrm{black}}\rangle|^2
```

で読む。

白猫選択時：

```text
C_white > C_select_threshold
C_black < C_reject_threshold
```

黒猫選択時：

```text
C_black > C_select_threshold
C_white < C_reject_threshold
```

---

# 8. ボゾン的対照写像

同一初期条件・同一C・同一Dで、フェルミオン的非線形増幅を外し、線形写像のみを使う。

```math
\begin{pmatrix}
a_{k+1}\\
b_{k+1}
\end{pmatrix}
=
U_C
\begin{pmatrix}
a_k\\
b_k
\end{pmatrix}
```

ここで `U_C` は線形かつ保存的とする。

期待：

```text
A/Bは偏っても再び戻る
一方が消えない
Sが±1へ固定されない
Cに一意の名前記録が残らない
```

---

# 9. 実験シナリオ

## Scenario 0: 自由発展対照

```text
C = off
D = off
mapping = baseline
```

確認事項：

```text
A≈0.5
B≈0.5
small oscillation only
no selection
```

---

## Scenario 1: Dのみ

```text
C = off
D = on
mapping = baseline
```

確認事項：

```text
D reads A≈0.5, B≈0.5
D does not trigger selection
```

---

## Scenario 2: Cのみ・フェルミオン的写像

```text
C = on
D = off
mapping = fermionic
```

確認事項：

```text
A/Bの一方が選択される
Cにwhiteまたはblackが記録される
```

---

## Scenario 3: D → C → D・フェルミオン的写像

完全猫型実験。

```text
Phase 1: stabilize gray cat
Phase 2: weak peek D
Phase 3: selection observation C
Phase 4: recursive fermionic evolution
Phase 5: result readout C
Phase 6: confirmation D
```

理想出力：

```text
D_before = (0.5, 0.5)
C_result = white or black
D_after  = (1,0) or (0,1)
```

---

## Scenario 4: D → C → D・ボゾン的写像

同一条件で写像だけ線形にする。

期待：

```text
D_before = (0.5,0.5)
C_result = unresolved or oscillatory
D_after  != stable (1,0) or (0,1)
```

---

## Scenario 5: D強度掃引

`g_D` を掃引する。

目的：

```text
readable but nonselective window
```

を探す。

記録：

```text
g_D
D_read_error
D_induced_bias
selection_triggered
```

---

## Scenario 6: C強度掃引

`eta_C`, `I_C`, `mu`, `gamma` を掃引する。

目的：

```text
minimum C strength that triggers stable selection
```

を探す。

---

## Scenario 7: C位相掃引

```text
theta_C in [0, 2pi)
```

を掃引する。

期待される選択境界：

```text
cos^2(theta_C) > 0.5 -> white tendency
cos^2(theta_C) < 0.5 -> black tendency
```

ただし境界付近では準安定・長時間選択・再振動が出る可能性を記録する。

---

# 10. 保存則

閉鎖系なので、非選択状態を単純削除しない。

最低限、次を追跡する。

```math
E_{\mathrm{total}}
=
E_A+E_B+E_C+E_D+E_{\mathrm{residual}}
```

期待：

```math
|E_{\mathrm{total}}(k)-E_{\mathrm{total}}(0)|
<
\varepsilon_E
```

A/Bの片方が見えなくなる場合でも、

```text
name coherence
phase order
harmonic order
localization
```

のどれがどこへ移ったかを記録する。

候補：

```text
selected channel
C memory channel
D memory channel
residual phase-dispersed modes
```

---

# 11. 判定指標

## 11.1 選択量

```math
S
=
p_A-p_B
```

判定：

```text
gray:
abs(S) < S_gray_limit
```

```text
white:
S > S_select_limit
```

```text
black:
S < -S_select_limit
```

推奨初期値：

```text
S_gray_limit = 0.05
S_select_limit = 0.95
```

---

## 11.2 D読出し誤差

```math
E_D
=
|D_A-p_A|+|D_B-p_B|
```

---

## 11.3 C名前識別誤差

```math
E_C
=
1-\max(C_{\mathrm{white}},C_{\mathrm{black}})
```

---

## 11.4 局在性

既存定義を利用する。

```math
L_A
=
\int |\Psi_A|^4 d\chi
```

```math
L_B
=
\int |\Psi_B|^4 d\chi
```

または既存の離散版を使う。

---

## 11.5 有効倍音次数

既存定義 `N_eff` を使用する。

```text
N_eff_A
N_eff_B
```

---

## 11.6 干渉可能性

D前後でA/Bの交差項を記録する。

```math
I_{AB}
=
2\operatorname{Re}
\left(
a^*b
\langle\Psi_A|\Psi_B\rangle
\right)
```

選択前：

```text
I_AB != 0
```

選択後：

```text
I_AB -> 0 or below observation limit
```

ただしDがこの量を直接選択しない範囲で読むこと。

---

# 12. 合格条件

本実験が成功したと判定する最低条件は次である。

## 12.1 灰色猫準安定

```text
C off, D off:
p_A and p_B remain within 0.5 ± delta_gray
for at least K_stable steps
```

---

## 12.2 盗み見成功

```text
D before C:
D_A ≈ 0.5
D_B ≈ 0.5
no selection triggered
```

---

## 12.3 C選択成功

```text
fermionic mapping:
C triggers stable white or black selection
```

---

## 12.4 C名前記録成功

```text
C result matches selected A/B name
```

---

## 12.5 D事後確認成功

```text
D after C:
reads (1,0) or (0,1)
within tolerance
```

---

## 12.6 ボゾン対照成立

```text
bosonic mapping:
same initial state and same C do not produce stable one-branch selection
```

---

## 12.7 保存則成立

```text
total registry conserved within tolerance
```

---

# 13. 失敗条件

次の場合は猫型実験として不成立とする。

1. CなしでA/Bの一方へ自然収束する。
2. Dだけで選択が起きる。
3. Dが0.5/0.5を読めない。
4. Cが名前を識別できない。
5. Cが単に外力として直接A/Bを書き込んでいる。
6. フェルミオン写像とボゾン写像で結果差がない。
7. 非選択側の量が保存先なしに消える。
8. 単純な各チャネル独立正規化が選択を人工生成している。
9. 数値丸め誤差だけで選択が決まる。
10. D前にすでに一方へ偏っている。

---

# 14. 推奨パラメータ掃引

```text
epsilon_exchange:
1e-4, 1e-3, 1e-2

g_D:
log sweep from 1e-6 to 1e-1

eta_C:
log sweep from 1e-6 to 1e-1

theta_C:
0 to 2pi, at least 256 divisions

mu:
0.01 to 2.0

gamma:
0.001 to 1.0

R:
0.0 to 1.0

T:
1.0 - R

steps:
128, 256, 512, 1024

random_seed:
multiple fixed seeds only for numerical-noise robustness tests
```

原則として、選択原因に外部乱数を使わない。

---

# 15. 出力仕様

以下を必ず生成する。

## 15.1 実行スクリプト

```text
run_gray_cat_white_black_peek_select_confirm_v1.py
```

---

## 15.2 JSON

```text
gray_cat_white_black_peek_select_confirm_result_v1/
gray_cat_white_black_peek_select_confirm_result_v1.json
```

---

## 15.3 CSV

最低限、各ステップで次を記録する。

```text
step
phase
mapping_type
p_A
p_B
S
Q
D_A
D_B
C_white
C_black
L_A
L_B
N_eff_A
N_eff_B
I_AB
E_A
E_B
E_C
E_D
E_residual
E_total
```

---

## 15.4 図

必須図：

1. `p_A`, `p_B` の時間発展
2. `S = p_A-p_B` の時間発展
3. D_before / C / D_after のイベント図
4. C_white / C_black の記録
5. フェルミオン型とボゾン型の比較
6. D強度掃引
7. C強度掃引
8. C位相 `theta_C` と選択結果の対応
9. `L_A`, `L_B`
10. `N_eff_A`, `N_eff_B`
11. 全保存レジストリ誤差
12. A/B波形スナップショット

推奨ファイル名：

```text
gray_cat_probabilities_v1.png
gray_cat_selection_order_parameter_v1.png
gray_cat_D_C_D_timeline_v1.png
gray_cat_C_name_memory_v1.png
gray_cat_fermionic_vs_bosonic_v1.png
gray_cat_D_strength_sweep_v1.png
gray_cat_C_strength_sweep_v1.png
gray_cat_C_phase_selection_map_v1.png
gray_cat_localization_exchange_v1.png
gray_cat_harmonic_exchange_v1.png
gray_cat_conservation_error_v1.png
gray_cat_waveform_snapshots_v1.png
```

---

# 16. 実験結果の要約テンプレート

```markdown
## 結果要約

- Cなし・Dなしでは、A/Bは0.5近傍で準安定に維持された。
- Dは選択前にA=..., B=...と読み、選択を誘発しなかった。
- C作用後、フェルミオン型写像では...ステップで白／黒が選択された。
- Cの名前記録は white=..., black=... であった。
- 選択後Dは A=..., B=... と読んだ。
- ボゾン型対照では、一状態選択は成立しなかった。
- 全保存レジストリ誤差は最大 ... であった。
```

---

# 17. Claude Codeへの実装指示

1. 既存の交換干渉散乱行列コードを可能な限り再利用する。
2. 既存の `L`, `N_eff`, `rho_chi`, 名前の毛、観測機Cの実装を再利用する。
3. まず重みレベルの最小モデルで成立条件を探索する。
4. 次に実波形レベルへ拡張する。
5. 単純な各チャネル独立正規化で選択を作らない。
6. 保存量の移動先を必ず追跡する。
7. Dは弱観測窓を探索し、選択を起こさない最大読出し強度を求める。
8. Cは名前読出しと選択を同時に行う。
9. フェルミオン型とボゾン型を同一初期条件で比較する。
10. 外部乱数を選択原因として使わない。
11. C位相 `theta_C` と選択先の対応を必ず掃引する。
12. 数値誤差依存性を、刻み幅・配列長・精度・seed変更で検証する。
13. 失敗結果も削除せずJSON/CSVへ残す。
14. 最後に、成功条件・不成功条件を自動判定する。
15. 実験メモをMarkdownで自動生成する。

---

# 18. 最終的に確認したい写像

本実験で確認したい最終写像は次である。

```text
準安定灰色猫
A=0.5, B=0.5
        |
        | weak peek D
        v
D reads A≈0.5, B≈0.5
        |
        | observation/interaction C
        | weak cos^2 bias
        | fermionic recursive amplification
        v
white selected
A≈1, B≈0
C records white
```

または、

```text
準安定灰色猫
A=0.5, B=0.5
        |
        | weak peek D
        v
D reads A≈0.5, B≈0.5
        |
        | observation/interaction C
        | weak cos^2 bias
        | fermionic recursive amplification
        v
black selected
A≈0, B≈1
C records black
```

最後に、

```text
D confirms selected state
```

とする。

---

# 19. 本実験の主張範囲

本実験は、現実の猫、生命、生死、意識を扱わない。

本実験が扱うのは、

```text
二つの巨視的に識別可能な名前付き波形状態が、
等配分準安定状態として共存し、
弱観測で等配分を確認され、
別の観測相互作用によって一方へ選択され、
結果が観測機内部に記録される写像
```

である。

本実験の成功時にも、標準量子力学全体の証明・否定を主張しない。

ただし、

```text
観測前の等配分準安定状態
観測相互作用による一状態選択
観測機への名前記録
選択前後の内部確認
```

を、閉じた数値系で同時に実装できるかを検証する。

---

# 20. 期待される最重要対照

```text
fermionic:
gray -> white or black

bosonic:
gray -> gray-like oscillation remains
```

この差が再現されるかを最重要判定とする。
