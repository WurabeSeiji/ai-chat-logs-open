# Codex向け作業指示
## Stage G：純パリティ区間で退化しない関係依存散乱の実装と、既存System A反復実験への限定接続

---

# 0. 目的

Stage Fでは、反転Candidate 1

\[
\theta_{\mathrm{eff}}
=
\theta_0
-
\kappa\rho(\theta_0)
\frac{c_A+c_B}{2}
\]

を既存System A実験コピーへ接続した。

その結果、純奇数倍音条件では

\[
c_A=c_B=-1
\]

が保存されるため、

\[
\theta_{\mathrm{eff}}
=
\theta_0+\kappa\rho(\theta_0)
\]

となり、毎衝突の状態を読んでいても、実質的には一定の反射率

\[
R_0\mapsto R_{\mathrm{eff}}
\]

へ置き換えただけになることが確認された。

本Stageでは追加観察を続けない。

次の一つの動的候補を実装する。

\[
\boxed{
\text{自己パリティの符号}
\times
\text{二波間の関係強度}
}
\]

を散乱角補正へ使う。

これにより、

- 純フェルミオン型同士では反射を強める方向
- 純ボゾン型同士では透過を強める方向
- 異型対では基準散乱
- 同じ純パリティ区間内でも、二波の関係が変化すれば散乱率が動く

という最小の関係依存散乱を実装する。

本Stageでは候補探索をしない。  
以下で定義する **relational_C1** 一候補だけを実装する。

---

# 1. 作業範囲

作業先を次に固定する。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/
codex_systemA_scattering_audit_v1/
stage_G_relational_scattering/
```

参照元はStage Fの独立統合コピーとする。

既存System A / System B原本は変更しない。

## 禁止事項

- 既存原本の変更
- Stage F成果物の上書き
- 配下外への出力
- Git操作
- Candidate 2・3の再実装
- 複数の新候補比較
- 結果を見た式変更
- 結果を見たパラメータ探索
- N体系への組込み
- 論文本文の変更
- Stage G後の自動継続

---

# 2. 最初に数式として確定すること

実装前に、次を短い数学ノートとして証明する。

## 2.1 反転Candidate 1の退化

ある反復区間で、

\[
c_A(n)=c_A^\ast,
\qquad
c_B(n)=c_B^\ast
\]

が一定なら、

\[
\theta_{\mathrm{eff}}(n)
=
\theta_0
-
\kappa\rho
\frac{c_A^\ast+c_B^\ast}{2}
\]

も一定である。

したがって、その区間の反転Candidate 1は、C0を一定の別反射率で走らせることと等価である。

純奇数対では、

\[
c_A^\ast=c_B^\ast=-1
\]

なので、

\[
\theta_{\mathrm{eff}}
=
\theta_0+\kappa\rho
\]

となる。

純偶数対では、

\[
c_A^\ast=c_B^\ast=+1
\]

なので、

\[
\theta_{\mathrm{eff}}
=
\theta_0-\kappa\rho
\]

となる。

この証明を、

```text
reports/00_reversed_C1_degeneracy_proof.md
```

へ保存する。

この証明後は、反転Candidate 1の追加R掃引を行わない。

---

# 3. 新しい関係依存量

復調した一次元カーネル波を、

\[
a(u),\qquad b(u)
\]

とする。

## 3.1 複素正規化重なり

\[
z_{AB}
=
\frac{\langle a,b\rangle}
{\|a\|\|b\|}
\]

とする。

## 3.2 関係強度

位相符号に依存しない関係強度を、

\[
\Gamma_{AB}
=
|z_{AB}|^2
=
\frac{|\langle a,b\rangle|^2}
{\|a\|^2\|b\|^2}
\]

と定義する。

Cauchy–Schwarz不等式より、

\[
0\le\Gamma_{AB}\le1
\]

である。

この量は、

- 同一波形なら1
- 直交波形なら0
- スペクトル交換・局在性交換・相対位相変化により中間値

を取る。

## 3.3 パリティ符号因子

\[
\bar c
=
\frac{c_A+c_B}{2}
\]

とする。

## 3.4 relational response

\[
F_{\mathrm{rel}}
=
-\bar c\,\Gamma_{AB}
\]

と定義する。

したがって、

### 純フェルミオン型同士

\[
c_A=c_B=-1
\]

より、

\[
F_{\mathrm{rel}}
=
+\Gamma_{AB}
\]

となり、関係強度に応じて反射を強める。

### 純ボゾン型同士

\[
c_A=c_B=+1
\]

より、

\[
F_{\mathrm{rel}}
=
-\Gamma_{AB}
\]

となり、関係強度に応じて反射を弱める。

### 純異型対

\[
c_A=-c_B
\]

より、

\[
\bar c=0
\]

したがって、

\[
F_{\mathrm{rel}}=0
\]

となり、基準散乱へ戻る。

## 3.5 散乱角

\[
\Delta\theta_{\mathrm{rel}}
=
\kappa\rho(\theta_0)F_{\mathrm{rel}}
\]

\[
\theta_{\mathrm{eff}}
=
\theta_0+\Delta\theta_{\mathrm{rel}}
\]

とする。

すなわち、

\[
\boxed{
\theta_{\mathrm{eff}}
=
\theta_0
-
\kappa\rho(\theta_0)
\frac{c_A+c_B}{2}
\Gamma_{AB}
}
\]

である。

散乱係数は既存と同じユニタリ形式を使う。

\[
t_{\mathrm{eff}}
=
e^{i\theta_{\mathrm{eff}}}
\cos\theta_{\mathrm{eff}}
\]

\[
r_{\mathrm{eff}}
=
-i e^{i\theta_{\mathrm{eff}}}
\sin\theta_{\mathrm{eff}}
\]

---

# 4. この候補を採用する理由

本候補は、反転Candidate 1の符号構造を保ちながら、一定パリティ区間で補正が一定になる問題を除く最小拡張である。

反転Candidate 1では、

\[
F=-\bar c
\]

であった。

relational_C1では、

\[
F=-\bar c\,\Gamma_{AB}
\]

とし、二波間の関係強度だけを追加する。

新しい自由な窓関数、接触作用素、任意のモード交換作用素は導入しない。

このStageでは、この式を採用して実装する。  
別の関係量を追加提案・探索してはならない。

---

# 5. 実装モード

同じSystem A独立コピーで、次の3モードを切替可能にする。

```python
scattering_mode = "C0"
scattering_mode = "reversed_C1"
scattering_mode = "relational_C1"
```

## C0

\[
\theta_{\mathrm{eff}}=\theta_0
\]

## reversed_C1

Stage Fの式をそのまま保持する。

\[
\theta_{\mathrm{eff}}
=
\theta_0-\kappa\rho\bar c
\]

## relational_C1

\[
\theta_{\mathrm{eff}}
=
\theta_0-\kappa\rho\bar c\,\Gamma_{AB}
\]

既存のC0とreversed_C1を書き換えず、比較対照として保持する。

---

# 6. 二層表現

Stage D〜Fで確立した二層表現をそのまま使用する。

## 6.1 全状態散乱

散乱・経路振幅・干渉・rawノルムは、

\[
512\times16
\]

全状態で計算する。

## 6.2 関係量読出し

\[
c_A,\quad c_B,\quad z_{AB},\quad\Gamma_{AB}
\]

は、A/B由来をそれぞれ正しい \(\eta\) モードへ射影し、対応搬送波を除去した一次元カーネル波から計算する。

全状態を単一搬送波で復調してはならない。

## 6.3 復調ゼロノルム

復調波のノルムが数値閾値未満の場合、黙って

\[
\Gamma_{AB}=0
\]

としてはならない。

その実行を数値不成立として記録し、条件・衝突番号・ノルムを報告する。

---

# 7. 事前登録パラメータ

\[
\kappa\in\{0.01,\ 0.1,\ 1\}
\]

とする。

結果に応じて追加・削除・微調整してはならない。

基準反射率は、まず既存代表条件、

\[
R_0=0.55
\]

を使う。

次に限定R集合、

\[
R_0\in\{0.55,\ 0.70\}
\]

だけを使用する。

広範囲R掃引は行わない。

目的は谷探しではなく、関係依存性が実際に動くかの検証である。

---

# 8. 実行順序

---

## Stage G-A：単体検証

以下の人工入力で、式の理論値を確認する。

### A1 同一純奇数波

\[
a=b=F_K
\]

理論値：

\[
c_A=c_B=-1,\qquad
\Gamma_{AB}=1
\]

\[
F_{\mathrm{rel}}=+1
\]

### A2 直交純奇数波

異なる直交奇数倍音集合を使う。

理論値：

\[
c_A=c_B=-1,\qquad
\Gamma_{AB}=0
\]

\[
F_{\mathrm{rel}}=0
\]

### A3 同一純偶数波

\[
a=b=B_K
\]

理論値：

\[
c_A=c_B=+1,\qquad
\Gamma_{AB}=1
\]

\[
F_{\mathrm{rel}}=-1
\]

### A4 純奇数×純偶数

理論値：

\[
\bar c=0
\]

\[
F_{\mathrm{rel}}=0
\]

### A5 位相差だけが異なる同一波形

\[
b=e^{i\phi}a
\]

\[
\phi\in
\left\{
0,\frac{\pi}{2},\pi
\right\}
\]

理論値：

\[
\Gamma_{AB}=1
\]

位相だけで関係強度が変化しないことを確認する。

単体検証が失敗した場合は、反復実験へ進まない。

---

## Stage G-B：既存代表条件

既存System A代表条件を使う。

\[
N_A=1,\qquad N_B=63,\qquad R_0=0.55
\]

既存と同じ初期波形、反復回数、評価関数を使う。

比較：

```text
C0
reversed_C1
relational_C1
```

各

\[
\kappa\in\{0.01,0.1,1\}
\]

で実行する。

更新方式は、Stage Fで差が機械精度程度であることを確認済みなら、主系列をraw_updateに固定してよい。

ただしC0の既存正規化系列との再現確認を1ケースだけ残す。

### 最重要判定

純奇数セクターで、

\[
c_A=c_B=-1
\]

が一定でも、

\[
\Gamma_{AB}(n)
\]

が変動し、

\[
R_{\mathrm{eff}}(n)
\]

が時間変動するかを確認する。

以下を機械判定する。

```text
gamma_constant
gamma_dynamic
R_eff_constant
R_eff_dynamic
```

閾値は実行前に固定する。

推奨：

\[
\max\Gamma-\min\Gamma > 10^{-10}
\]

なら `gamma_dynamic`。

\[
\max R_{\mathrm{eff}}-\min R_{\mathrm{eff}} > 10^{-10}
\]

なら `R_eff_dynamic`。

---

## Stage G-C：31系列の限定再実行

Stage Fと同じ31系列条件を使う。

新しい系列探索はしない。

比較：

```text
C0
reversed_C1
relational_C1
```

評価点：

\[
31,62,93,124,155,186,217,247,248,279
\]

最低限保存：

```text
iteration
return_error
exchange_measure
L_difference
N_eff_difference
c_A
c_B
c_mean
Gamma_AB
theta_eff
R_eff
```

目的は、relational_C1が248帰還を回復するかどうかに合わせることではない。

関係強度が時間変化した結果として、既存帰還構造が

- 保存
- 位相移動
- 崩壊
- 別周期化

のどれになるかを記録する。

---

# 9. 毎衝突保存項目

最低限、以下を保存する。

```text
collision_index
scattering_mode
kappa
R0
theta0
rho

c_A
c_B
c_mean

overlap_complex_real
overlap_complex_imag
overlap_abs
Gamma_AB

candidate_response
delta_theta
theta_eff
R_eff
T_eff

L_A
L_B
N_eff_A
N_eff_B

spectral_similarity_A_to_initial_A
spectral_similarity_A_to_initial_B
spectral_similarity_B_to_initial_A
spectral_similarity_B_to_initial_B

path_A_to_A_norm
path_B_to_A_norm
path_B_to_B_norm
path_A_to_B_norm

interference_A
interference_B

raw_norm_A
raw_norm_B
next_state_norm_A
next_state_norm_B

boson_weight_A
fermion_weight_A
boson_weight_B
fermion_weight_B

unitarity_residual
orthogonality_residual
path_sum_residual_A
path_sum_residual_B
total_norm_residual
demodulation_reconstruction_residual_A
demodulation_reconstruction_residual_B
```

旧 `B_to_A_transfer` 名称を使う場合は、

```text
spectral similarity; not path flux
```

と必ず明記する。

---

# 10. 比較指標

以下をC0、reversed_C1、relational_C1で比較する。

## 10.1 動的性

\[
\Delta\Gamma
=
\max_n\Gamma_{AB}(n)
-
\min_n\Gamma_{AB}(n)
\]

\[
\Delta R_{\mathrm{eff}}
=
\max_nR_{\mathrm{eff}}(n)
-
\min_nR_{\mathrm{eff}}(n)
\]

## 10.2 周期・帰還

- 最小帰還誤差
- 最小帰還誤差の反復回
- 自己相関ピーク
- 既存32回交換との偏差
- 248回帰還誤差

## 10.3 局在性交換

- \(\min|L_A-L_B|\)
- その反復回
- \(\min|N_{\mathrm{eff},A}-N_{\mathrm{eff},B}|\)
- その反復回
- A/Bスペクトル類似度の交差回数

## 10.4 関係量との相関

次を計算する。

\[
\operatorname{corr}(\Gamma_{AB},R_{\mathrm{eff}})
\]

\[
\operatorname{corr}(\Gamma_{AB},|L_A-L_B|)
\]

\[
\operatorname{corr}(\Gamma_{AB},|N_{\mathrm{eff},A}-N_{\mathrm{eff},B}|)
\]

定数系列で相関係数が未定義の場合は、0と置かず `not_defined_constant_series` と記録する。

---

# 11. 判定分類

結果を次の事前登録分類で整理する。

```text
relational_term_inactive
constant_relation_reparameterization
dynamic_relation_dynamic_scattering
dynamic_relation_no_observable_change
dynamic_relation_period_shift
dynamic_relation_amplitude_shift
dynamic_relation_period_and_amplitude_shift
new_quasistable_cycle
fixed_point_convergence
divergent
numerically_unstable
not_reproducible
```

## 中心判定

本Stageの主判定は一つである。

\[
\boxed{
\Gamma_{AB}(n)
\text{ が変動し、それに伴って }
R_{\mathrm{eff}}(n)
\text{ が変動するか}
}
\]

変動しない場合は、relational_C1もこの実験条件では一定Rへの再パラメータ化に退化したと結論する。

その場合、新候補を自動生成して続行してはならない。

---

# 12. 数値制約

全ケースで、

\[
0\le\Gamma_{AB}\le1
\]

を数値許容誤差内で確認する。

\[
0\le\theta_{\mathrm{eff}}\le\frac{\pi}{2}
\]

を破る場合は、自動クリップせず、その条件を範囲違反として停止・報告する。

保存則：

\[
|r|^2+|t|^2=1
\]

\[
r^\ast t+t^\ast r=0
\]

経路和、全ノルム、復調再構成残差を保存する。

NaN、Inf、ゼロノルム、範囲違反を黙って除外してはならない。

---

# 13. 出力構成

```text
stage_G_relational_scattering/
├── code/
│   ├── relational_scattering.py
│   ├── system_A_stage_G_copy.py
│   ├── run_stage_G_unit_tests.py
│   ├── run_stage_G_repeated_systemA.py
│   ├── run_stage_G_31_series.py
│   ├── analyze_stage_G.py
│   ├── test_stage_G.py
│   └── build_manifest.py
├── source_copy/
│   ├── README.md
│   └── <Stage F integration snapshot>
├── data/
│   ├── reference_hashes_before.json
│   ├── reference_hashes_after.json
│   ├── reference_hash_comparison.json
│   ├── stage_G_unit_test_results.csv
│   ├── stage_G_collision_results.csv
│   ├── stage_G_run_summary.csv
│   ├── stage_G_31_series_results.csv
│   ├── stage_G_correlation_results.csv
│   ├── stage_G_numerical_residuals.csv
│   ├── stage_G_summary.json
│   └── stage_G_final_states.npz
├── figures/
│   ├── Gamma_and_R_eff_by_collision.png
│   ├── C0_C1_relational_L_exchange.png
│   ├── C0_C1_relational_N_eff_exchange.png
│   ├── relation_vs_localization_difference.png
│   ├── cycle_and_return_error_comparison.png
│   └── return_error_31_series.png
├── reports/
│   ├── 00_reversed_C1_degeneracy_proof.md
│   ├── 01_relational_candidate_definition.md
│   ├── 02_unit_test_report.md
│   ├── 03_existing_systemA_comparison.md
│   ├── 04_dynamic_relation_analysis.md
│   ├── 05_31_series_comparison.md
│   ├── 06_numerical_invariants.md
│   └── Stage_G_report.md
├── logs/
└── manifest.json
```

---

# 14. 報告書で分離する項目

## コード上の事実

実際の実装と既存System A構造。

## 数学的帰結

反転Candidate 1の退化、\(\Gamma\)の範囲、純状態での符号。

## モデル定義

relational_C1の式。

## 作業仮説

パリティ符号と関係強度が散乱角を制御するという物理的解釈。

## 数値観測

\(\Gamma\)、\(R_{\mathrm{eff}}\)、周期、局在性交換、帰還誤差。

## 未導出

\(\kappa\)、包絡 \(\rho\)、重なり二乗を相互作用強度と読む根拠。

## 棄却・保留

関係量が一定、応答が実質一定R、範囲違反、数値不安定など。

---

# 15. 解釈制限

次の主張を自動的に行わない。

- フェルミオンを導出した
- パウリ排他原理を証明した
- ボゾン・フェルミオン散乱を再現した
- 31周期を説明した
- 関係量が自然界の相互作用そのものである

本Stageで言えるのは、

\[
\Gamma_{AB}
\]

を関係強度として導入したとき、純パリティ区間内で散乱率が動的に変化するか、その変化が既存System Aの局在性交換へどう作用したかまでである。

---

# 16. 原本保全

実行前後に、参照した全既存ファイルについて、

```text
path
size
mtime
sha256
```

を保存する。

一つでも差異があればStage G失敗として停止する。

Git操作は行わない。

---

# 17. 実行順序

1. 参照元特定
2. 事前SHA-256
3. 数学的退化証明作成
4. relational_C1実装
5. 単体テスト
6. C0再現確認
7. Stage G-B代表条件
8. Stage G-C 31系列
9. 動的性・周期・相関解析
10. 数値残差集計
11. 事後SHA-256比較
12. 報告書・manifest作成
13. 停止

単体テストまたはC0再現に失敗した場合は、反復実験へ進まない。

---

# 18. 最終停止条件

以下を完了した時点で必ず停止する。

- 反転Candidate 1の一定パリティ区間での退化証明
- relational_C1一候補の実装
- 単体検証
- 既存System A代表条件での3モード比較
- 31系列の限定比較
- \(\Gamma_{AB}\) と \(R_{\mathrm{eff}}\) の動的性判定
- 局在性交換・周期・帰還への影響評価
- 数値残差確認
- 原本SHA-256前後一致
- `Stage_G_report.md`
- `manifest.json`

停止時に明示する。

```text
Stage G 完了。
既存System A / System B原本は変更していない。
新規実装はrelational_C1一候補のみ。
Candidate 2・3の追加実装は行っていない。
N体系へは組み込んでいない。
論文本文は変更していない。
人間の承認待ち。
```
