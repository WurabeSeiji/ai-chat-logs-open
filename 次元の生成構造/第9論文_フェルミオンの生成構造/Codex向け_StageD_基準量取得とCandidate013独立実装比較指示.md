# Codex向け作業指示  
## Stage D：現行基準量の正式取得と Candidate 0・1・3 の独立実装・比較

### 0. 目的

本作業の目的は、偶数倍音ボゾン型・奇数倍音フェルミオン型という作業仮説を既存 System A / System B 本体へ直ちに組み込むことではない。

まず、Stage B で確認された現行散乱核の基準量を正式な共通APIとして取得・保存し、その上で、以下の三候補を**独立テスト環境**に実装して比較する。

- Candidate 0：現行のパリティ盲基準
- Candidate 1：単波パリティ指標が散乱角を制御する候補
- Candidate 3：二波の積の \(\mathbb Z_2\) 合成則が散乱角を制御する候補

既存 System A / System B 本体、既存データ、既存結果、既存論文原本は変更しない。

---

# 1. 作業範囲

作業ディレクトリは、以下の専用ディレクトリ配下に限定する。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/
codex_systemA_scattering_audit_v1/stage_D_candidate_implementation/
```

このディレクトリが存在しない場合は作成してよい。

配下外のファイルは、読み取り・参照・SHA-256比較のみ許可する。

禁止事項：

- 配下外ファイルの変更
- 既存 System A / System B 本体への候補実装
- 既存CSV・JSON・図・報告書の上書き
- 既存出力ディレクトリへの書込み
- 外部コピーの import 実行
- 仮説に都合のよいパラメータ探索
- Candidate 2 の実装
- N体系への組込み
- Stage D 完了後の自動継続

---

# 2. Stage Bで確定した基準事実

以下を再解釈せず、確定済みの基準事実として扱う。

## 2.1 現行散乱核

現行散乱核は、

\[
\widetilde a = r a + t b,
\qquad
\widetilde b = t a + r b
\]

である。

散乱係数は、

\[
\theta_0 = \arcsin\sqrt{R}
\]

\[
t_0 = e^{i\theta_0}\cos\theta_0,
\qquad
r_0 = -i e^{i\theta_0}\sin\theta_0
\]

である。

したがって、

\[
|r_0|^2 = R,
\qquad
|t_0|^2 = 1-R
\]

\[
|r_0|^2 + |t_0|^2 = 1
\]

\[
r_0^*t_0+t_0^*r_0=0
\]

を満たす。

## 2.2 現行核の性質

現行核は入力波の偶数倍音・奇数倍音構成を参照しない。

したがって、

\[
\|r_0 a\|^2 = R\|a\|^2
\]

\[
\|t_0 b\|^2 = (1-R)\|b\|^2
\]

であり、同じノルムの純奇数波・純偶数波に同じ経路応答を返す。

現行核は型を読まず、型依存散乱差を生成しない。

## 2.3 `B_to_A_transfer`

既存の `B_to_A_transfer` は、経路量

\[
\|t b\|^2
\]

ではない。

A出力の倍音パワー分布と初期B倍音パワー分布のスペクトル類似度である。

今後も既存名称を参照する場合は、必ず

```text
spectral_similarity_to_initial_B
```

と併記し、経路交換量と混同しないこと。

## 2.4 現行設定の干渉項

既存 System A の基準設定では、A/Bの \(\eta\) モードが直交するため、

\[
2\operatorname{Re}\langle r a, t b\rangle
\]

はほぼゼロである。

これは一般定理ではなく、既存入力構成に依存する数値事実である。

---

# 3. 今回正式に取得・保存する基準量

新候補実装前後を同一物差しで比較するため、全候補で同じ `ScatteringResult` を返すこと。

最低限、以下を保存する。

## 3.1 入力状態

```python
input_a
input_b
input_norm_a
input_norm_b
```

## 3.2 パリティ関連量

半周期移動を、

\[
(P\psi)(u)=\psi(u+\pi)
\]

とする。

射影は、

\[
\Pi_B=\frac{I+P}{2},
\qquad
\Pi_F=\frac{I-P}{2}
\]

とする。

非正規化半周期相関：

\[
C_\pi^{\mathrm{raw}}[\psi]
=
\langle\psi,P\psi\rangle
\]

正規化実数指標：

\[
c_\pi[\psi]
=
\frac{\operatorname{Re}\langle\psi,P\psi\rangle}
{\|\psi\|^2}
\]

ボゾン型・フェルミオン型比率：

\[
p_B[\psi]
=
\frac{\|\Pi_B\psi\|^2}{\|\psi\|^2}
\]

\[
p_F[\psi]
=
\frac{\|\Pi_F\psi\|^2}{\|\psi\|^2}
\]

返却フィールド：

```python
parity_correlation_raw_a: complex
parity_correlation_raw_b: complex
parity_indicator_a: float
parity_indicator_b: float
boson_weight_a: float
boson_weight_b: float
fermion_weight_a: float
fermion_weight_b: float
```

パリティ判定は、必ず搬送波を逆シフトした後のカーネル座標 \(u\) で行うこと。

密度 \(|\psi|^2\) だけで判定してはならない。

## 3.3 状態依存散乱角

```python
theta_0
rho
candidate_response
delta_theta
theta_eff
```

## 3.4 散乱係数

```python
r_eff
t_eff
reflection_probability
transmission_probability
unitarity_residual
orthogonality_residual
```

定義：

\[
R_{\mathrm{eff}}=|r_{\mathrm{eff}}|^2
\]

\[
T_{\mathrm{eff}}=|t_{\mathrm{eff}}|^2
\]

\[
\varepsilon_U
=
\left|
|r_{\mathrm{eff}}|^2+|t_{\mathrm{eff}}|^2-1
\right|
\]

\[
\varepsilon_{rt}
=
\left|
r_{\mathrm{eff}}^*t_{\mathrm{eff}}
+
t_{\mathrm{eff}}^*r_{\mathrm{eff}}
\right|
\]

## 3.5 経路別振幅

```python
path_a_to_a_amplitude = r_eff * a
path_b_to_a_amplitude = t_eff * b
path_b_to_b_amplitude = r_eff * b
path_a_to_b_amplitude = t_eff * a
```

## 3.6 経路別ノルム

```python
path_a_to_a_norm
path_b_to_a_norm
path_b_to_b_norm
path_a_to_b_norm
```

定義：

\[
N_{A\to A}=\|r_{\mathrm{eff}}a\|^2
\]

\[
N_{B\to A}=\|t_{\mathrm{eff}}b\|^2
\]

\[
N_{B\to B}=\|r_{\mathrm{eff}}b\|^2
\]

\[
N_{A\to B}=\|t_{\mathrm{eff}}a\|^2
\]

## 3.7 干渉項

```python
interference_in_a
interference_in_b
```

\[
I_A
=
2\operatorname{Re}
\langle r_{\mathrm{eff}}a,t_{\mathrm{eff}}b\rangle
\]

\[
I_B
=
2\operatorname{Re}
\langle r_{\mathrm{eff}}b,t_{\mathrm{eff}}a\rangle
\]

## 3.8 正規化前出力

```python
raw_output_a
raw_output_b
raw_output_norm_a
raw_output_norm_b
```

\[
a_{\mathrm{raw}}
=
r_{\mathrm{eff}}a+t_{\mathrm{eff}}b
\]

\[
b_{\mathrm{raw}}
=
t_{\mathrm{eff}}a+r_{\mathrm{eff}}b
\]

整合残差：

\[
\varepsilon_A
=
\left|
\|a_{\mathrm{raw}}\|^2
-
\left(
N_{A\to A}+N_{B\to A}+I_A
\right)
\right|
\]

\[
\varepsilon_B
=
\left|
\|b_{\mathrm{raw}}\|^2
-
\left(
N_{B\to B}+N_{A\to B}+I_B
\right)
\right|
\]

```python
path_sum_residual_a
path_sum_residual_b
```

## 3.9 正規化後出力

```python
normalized_output_a
normalized_output_b
normalized_output_norm_a
normalized_output_norm_b
```

正規化後状態は反復計算用の補助状態とする。

物理的な出力比較の正本は、正規化前の `raw_output_*` とする。

## 3.10 出力パリティ

正規化前出力について、

```python
raw_output_parity_indicator_a
raw_output_parity_indicator_b
raw_output_boson_weight_a
raw_output_boson_weight_b
raw_output_fermion_weight_a
raw_output_fermion_weight_b
```

を保存する。

## 3.11 非線形半周期同変性残差

候補1・3では散乱係数が入力依存となるため、通常の線形交換子ではなく、状態依存同変性残差を使う。

\[
\varepsilon_P(a,b)
=
\frac{
\left\|
\mathcal S(Pa,Pb)
-
(P\oplus P)\mathcal S(a,b)
\right\|
}{
\|\mathcal S(a,b)\|
}
\]

返却名：

```python
half_shift_equivariance_residual
```

`commutator` という名称は使用しない。

---

# 4. 共通散乱API

独立モジュールとして、最低限以下の構造を作る。

```python
from dataclasses import dataclass
from typing import Literal
import numpy as np


CandidateName = Literal["C0", "C1", "C3"]


@dataclass(frozen=True)
class ScatteringResult:
    candidate: str
    kappa: float

    theta_0: float
    rho: float
    candidate_response: float
    delta_theta: float
    theta_eff: float

    r_eff: complex
    t_eff: complex
    reflection_probability: float
    transmission_probability: float
    unitarity_residual: float
    orthogonality_residual: float

    parity_correlation_raw_a: complex
    parity_correlation_raw_b: complex
    parity_indicator_a: float
    parity_indicator_b: float
    boson_weight_a: float
    boson_weight_b: float
    fermion_weight_a: float
    fermion_weight_b: float

    path_a_to_a_amplitude: np.ndarray
    path_b_to_a_amplitude: np.ndarray
    path_b_to_b_amplitude: np.ndarray
    path_a_to_b_amplitude: np.ndarray

    path_a_to_a_norm: float
    path_b_to_a_norm: float
    path_b_to_b_norm: float
    path_a_to_b_norm: float

    interference_in_a: float
    interference_in_b: float

    raw_output_a: np.ndarray
    raw_output_b: np.ndarray
    raw_output_norm_a: float
    raw_output_norm_b: float

    path_sum_residual_a: float
    path_sum_residual_b: float

    normalized_output_a: np.ndarray
    normalized_output_b: np.ndarray
    normalized_output_norm_a: float
    normalized_output_norm_b: float

    raw_output_parity_indicator_a: float
    raw_output_parity_indicator_b: float
    raw_output_boson_weight_a: float
    raw_output_boson_weight_b: float
    raw_output_fermion_weight_a: float
    raw_output_fermion_weight_b: float

    half_shift_equivariance_residual: float
```

必要な補助メタデータは追加してよいが、上記を削除してはならない。

---

# 5. 共通の状態依存散乱形式

基準角：

\[
\theta_0=\arcsin\sqrt R
\]

共通包絡：

\[
\rho(\theta_0)
=
\frac{2}{\pi}
\theta_0
\left(
\frac{\pi}{2}-\theta_0
\right)
\]

候補応答を \(F_j\) として、

\[
\Delta\theta_j
=
\kappa\,\rho(\theta_0)\,F_j
\]

\[
\theta_{\mathrm{eff}}
=
\theta_0+\Delta\theta_j
\]

とする。

ただし、

\[
0\le\theta_{\mathrm{eff}}\le\frac{\pi}{2}
\]

を満たす必要がある。

単純なクリップを使う場合は、クリップ前後の値を別々に保存すること。

より望ましいのは、今回指定する \(\kappa\) と入力で範囲外にならないことを事前検証することである。

散乱係数：

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

この形式によりユニタリ性を保つ。

---

# 6. Candidate 0

現行基準。

\[
F_0=0
\]

\[
\Delta\theta_0=0
\]

\[
\theta_{\mathrm{eff}}=\theta_0
\]

Candidate 0 は Stage B の現行核を再現しなければならない。

再現誤差を必ず報告する。

---

# 7. Candidate 1

Candidate 1は、単波の半周期パリティ指標が散乱角を直接制御する最小作業仮説とする。

\[
F_1
=
\frac{c_\pi[a]+c_\pi[b]}{2}
\]

ここで、

\[
-1\le F_1\le1
\]

である。

純偶数型同士では、

\[
F_1=+1
\]

純奇数型同士では、

\[
F_1=-1
\]

純偶数型と純奇数型では、

\[
F_1=0
\]

となる。

Candidate 1には、交差半周期相関を追加しない。

今回の最初の実装では、

\[
F_1=\frac{c_A+c_B}{2}
\]

のみを使用する。

これは条件分岐ではないが、パリティ固有値を散乱角へ直接結び付ける仮説であることを報告書へ明記する。

---

# 8. Candidate 3

Candidate 3は、二波の点ごとの積が持つ \(\mathbb Z_2\) パリティ合成則を散乱角へ入力する候補とする。

まず、

\[
\zeta(u)=a(u)b(u)
\]

を作る。

積状態のパリティ指標を、

\[
F_3
=
c_\pi[\zeta]
=
\frac{
\operatorname{Re}
\langle
a b,
P(a b)
\rangle
}{
\|ab\|^2
}
\]

と定義する。

\[
-1\le F_3\le1
\]

である。

純状態では、

\[
F\times F \Rightarrow F_3=+1
\]

\[
B\times B \Rightarrow F_3=+1
\]

\[
F\times B \Rightarrow F_3=-1
\]

\[
B\times F \Rightarrow F_3=-1
\]

となる。

Candidate 3は、フェルミオン型同士とボゾン型同士を区別する候補ではない。

同型対と異型対を区別する候補である。

この点を報告書で明確にする。

積 \(ab\) のノルムが数値的にゼロまたは極小の場合は、例外処理で黙って0を返してはならない。

閾値と処理規則を明記し、そのケースを別途報告すること。

---

# 9. 入力波形

全入力は同一格子、同一領域、同一ノルムで作る。

## 9.1 純奇数局在波

\[
F_K(u)
=
\frac{1}{\sqrt K}
\sum_{j=0}^{K-1}
\cos((2j+1)u)
\]

## 9.2 純偶数局在波

\[
B_K(u)
=
\frac{1}{\sqrt K}
\sum_{j=1}^{K}
\cos(2ju)
\]

ゼロモードは含めない。

## 9.3 混合波

\[
M_{50,\varphi}(u)
=
\frac{
B_K(u)+e^{i\varphi}F_K(u)
}{
\left\|
B_K+e^{i\varphi}F_K
\right\|
}
\]

位相は、

\[
\varphi
\in
\left\{
0,\frac{\pi}{2},\pi
\right\}
\]

を使用する。

## 9.4 成分数

最低限、

\[
K\in\{1,4,8,16\}
\]

を実行する。

計算量が軽い場合は \(K=32\) を追加してよい。

ただし結果に都合のよい \(K\) のみを選んではならない。

---

# 10. 入力組合せ

最低限、以下を実行する。

```text
F x F
B x B
F x B
B x F
M(phi=0) x M(phi=0)
M(phi=pi/2) x M(phi=pi/2)
M(phi=pi) x M(phi=pi)
F x M(phi=0)
F x M(phi=pi/2)
F x M(phi=pi)
B x M(phi=0)
B x M(phi=pi/2)
B x M(phi=pi)
```

A/B交換対称性を確認するため、

```text
input_a, input_b
input_b, input_a
```

を両方実行すること。

---

# 11. 基準反射率と感度係数

基準反射率は、最低限、

\[
R\in\{0,\ 0.1,\ 0.5,\ 0.9,\ 1\}
\]

を使う。

既存研究で重要な特定の \(R\) を追加する場合は、基準集合を削除せず追加する。

感度係数は、

\[
\kappa\in\{0,\ 0.01,\ 0.1,\ 1\}
\]

とする。

\(\kappa=0\) は全候補で Candidate 0 と一致しなければならない。

結果を見て \(\kappa\) を追加・変更してはならない。

範囲外の \(\theta_{\mathrm{eff}}\) が発生した場合は、件数・入力条件・クリップ前値を報告する。

---

# 12. 必須検証

## 12.1 Candidate 0再現

Stage B基準値との差を計算する。

```text
max_abs_error
max_relative_error
```

を報告する。

## 12.2 ユニタリ性

全ケースで、

\[
\varepsilon_U
\]

\[
\varepsilon_{rt}
\]

を保存し、最大値を報告する。

## 12.3 経路和整合性

\[
\varepsilon_A
\]

\[
\varepsilon_B
\]

の最大値を報告する。

## 12.4 A/B交換対称性

入力を交換したとき、対応する出力・経路量が適切に交換されることを確認する。

非対称が出た場合は、候補式・実装・数値誤差のどれに由来するかを切り分ける。

## 12.5 半周期同変性

\[
\varepsilon_P
\]

を全ケースで計測する。

## 12.6 端点

\[
R=0
\]

\[
R=1
\]

で共通包絡 \(\rho=0\) となり、Candidate 1・3が基準核へ戻ることを確認する。

## 12.7 \(\kappa=0\)

すべての候補で Candidate 0 と一致すること。

## 12.8 純状態の理論値

Candidate 1：

```text
B x B -> F1 = +1
F x F -> F1 = -1
B x F -> F1 = 0
F x B -> F1 = 0
```

Candidate 3：

```text
B x B -> F3 = +1
F x F -> F3 = +1
B x F -> F3 = -1
F x B -> F3 = -1
```

数値誤差を報告する。

## 12.9 混合状態

混合状態について、\(F_1\) と \(F_3\) がどのように異なるかを表と図で示す。

特に相対位相 \(\varphi\) 依存性を確認する。

---

# 13. 出力ファイル

最低限、以下を作成する。

```text
stage_D_candidate_implementation/
├── code/
│   ├── scattering_api.py
│   ├── parity_metrics.py
│   ├── candidate_responses.py
│   ├── wave_generators.py
│   ├── run_stage_D_candidate_comparison.py
│   └── test_stage_D_candidates.py
├── data/
│   ├── stage_D_full_results.csv
│   ├── stage_D_summary.json
│   ├── stage_D_candidate0_baseline.csv
│   ├── stage_D_candidate1_results.csv
│   └── stage_D_candidate3_results.csv
├── figures/
│   ├── candidate_response_comparison.png
│   ├── reflection_probability_comparison.png
│   ├── path_norm_comparison.png
│   ├── parity_output_comparison.png
│   └── mixed_phase_dependence.png
├── reports/
│   ├── 00_execution_scope_and_hashes.md
│   ├── 01_common_api_and_measurement_definition.md
│   ├── 02_candidate0_baseline_reproduction.md
│   ├── 03_candidate1_results.md
│   ├── 04_candidate3_results.md
│   ├── 05_candidate1_vs_candidate3.md
│   ├── 06_numerical_invariants_and_residuals.md
│   └── Stage_D_report.md
└── manifest.json
```

SVGも作成可能ならPNGと同時に作成してよい。

既存ファイルを上書きしない。

---

# 14. 報告書の分類規則

報告書では、必ず以下を分離する。

## コード上の事実

実際のコードから確認したこと。

## 数学的帰結

定義と式から必ず従うこと。

## モデル定義

今回採用した定義。

## 作業仮説

物理的同定や散乱角との対応。

## 数値観測

実行結果として得た値。

## 未導出

現時点で説明できていない接続。

## 棄却候補

数値的不整合、恒等的ゼロ、対称性破れ、範囲破れなどにより採用困難な候補。

Candidate 1または3を「正しい」と結論してはならない。

このStageでは、候補の数理的・数値的性質を比較するだけである。

---

# 15. 判定基準

## Candidate 0

Stage Bを再現できなければ実装失敗。

## Candidate 1

以下を満たすこと。

- 純B/Bと純F/Fで反対方向の角度補正
- B/Fでは自己項がゼロ
- A/B交換対称
- ユニタリ性維持
- 半周期同変性が数値精度内

## Candidate 3

以下を満たすこと。

- 同型対で \(F_3=+1\)
- 異型対で \(F_3=-1\)
- 混合状態で連続的中間値を返す
- A/B交換対称
- ユニタリ性維持
- 半周期同変性が数値精度内

## 比較上の重要点

Candidate 1は各単波のパリティを読む。

Candidate 3は二波の積パリティを読む。

この違いを最終報告で明確に比較する。

---

# 16. 原本保全

実行前後で、参照した既存コード・既存文書・既存データについてSHA-256を取得する。

変更が検出された場合は直ちに停止する。

新規作成ファイルのみを `manifest.json` に記録する。

Gitへの追加、コミット、ステージングは行わない。

---

# 17. 最終停止条件

以下が完了した時点で必ず停止する。

- Candidate 0の基準再現
- Candidate 1の独立実装と検証
- Candidate 3の独立実装と検証
- 共通APIによる全量保存
- 数値不変量と残差の報告
- Candidate 1対Candidate 3の比較
- Stage_D_report.md の作成
- manifest.json の作成
- 原本SHA-256前後一致確認

停止時に、明示的に次を報告する。

```text
Stage D 完了。
既存 System A / System B 本体への組込みは行っていない。
Candidate 2は実装していない。
N体系への組込みは行っていない。
人間の承認待ち。
```
