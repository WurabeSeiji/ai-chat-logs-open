# Codex向け指示書
## System A 散乱コード監査と、偶数倍音ボゾン型・奇数倍音フェルミオン型仮説に基づく新散乱関数設計

作成日: 2026-07-28  
対象: 現行 System A の交換散乱コード群を主対象とし、System B は散乱核共有の比較参照対象とする  
目的: 既存実装の事実を完全に監査した上で、入力波形の偶数倍音・奇数倍音構造を実際の散乱応答へ反映できる新しい散乱関数を設計する

---

# 0. 最重要方針

この作業では、既存仮説をコードへ直ちに埋め込んではならない。

まず、現在の散乱コードが、

- どの入力を受け取り、
- どこで散乱係数を作り、
- どこで反射・透過・交換を計算し、
- どこで正規化し、
- どの配列を出力し、
- どの評価量を最終判定に使っているか

を、実コードに基づいて完全に追跡すること。

推測、既存文書の要約、過去の説明の再利用だけで済ませてはならない。

既存コードの監査が終わる前に、偶数倍音・奇数倍音の分岐、ボゾン型・フェルミオン型の条件分岐、反射抑制、交換抑制などを実装してはならない。

---

# 1. 背景仮説

本検討では、閉じた周期系上の波について、半周期移動演算子を

\[
(P\psi)(u)=\psi(u+\pi)
\]

と定義する。

偶数倍音のみからなる波をボゾン型、奇数倍音のみからなる波をフェルミオン型と仮定する。

\[
P\psi_B=+\psi_B
\]

\[
P\psi_F=-\psi_F
\]

対応する射影は、

\[
\Pi_B=\frac{I+P}{2}
\]

\[
\Pi_F=\frac{I-P}{2}
\]

である。

混合波は、

\[
\psi=\psi_B+\psi_F
\]

と分解する。

ただし、この分類は現時点では作業仮説である。

重要なのは、この仮説をコードへ手作業で条件分岐として埋め込むことではない。

必要なのは、同一の散乱法則の内部から、入力波形の偶奇構造による応答差が自然に生じる散乱関数を設計できるかを検討することである。

---

# 2. 現時点で判明している既存散乱核

既存文書では、System A の中核更新として、概ね次の形が報告されている。

```python
a_next = normalize(r * a + t * b)
b_next = normalize(t * a + r * b)
a, b = a_next, b_next
```

また、散乱係数は反射率パラメータ `R` から、概ね次の形で作られている。

```python
def scattering_coefficients(reflection_rate: float):
    if reflection_rate < 0.0 or reflection_rate > 1.0:
        raise ValueError("reflection rate must be in [0, 1]")

    delta_f = 2.0 * math.asin(math.sqrt(reflection_rate))
    half_delta = 0.5 * delta_f
    phase = complex(math.cos(half_delta), math.sin(half_delta))

    t = phase * math.cos(half_delta)
    r = -1j * phase * math.sin(half_delta)
    return complex(t), complex(r), abs(t) ** 2, abs(r) ** 2
```

すなわち、現行実装では、

\[
R=|r|^2,
\qquad
T=|t|^2=1-R
\]

が外部パラメータとして先に与えられ、入力波形の偶数倍音・奇数倍音構造から散乱係数が導かれているわけではない可能性が高い。

この理解を、必ず実コードで検証すること。

---

# 3. 対象コードの特定

最初に、以下のコード、およびこれらが import する依存コードを特定すること。

主対象候補:

```text
run_system_A_localization_exchange_R_sweep_preliminary_v1.py
run_minimal_system_B_gray_direct_check_v5.py
phase5_eigenphase_resonance_v2.py
run_two_physical_roots_multiprecision_v1.py
```

加えて、以下の定義元を必ず追跡すること。

```text
normalize
normalize_pair
scattering_coefficients
row_for_state
kernel generator
harmonic basis generator
carrier shift / inverse carrier shift
B_to_A evaluator
gray_error / gray_depth evaluator
```

同名関数が複数存在する場合は、実際に実行対象から参照される定義を import graph と呼出経路で確定すること。

---

# 4. Stage A: 現行コード監査

## 4.1 変更禁止

Stage A では、既存コードを一切変更してはならない。

許可されるのは、

- 読み取り
- 呼出関係の整理
- コピー先での診断用ラッパー作成
- 非侵襲的ログ出力
- 独立した監査スクリプト作成

のみである。

原本コード、原本データ、既存結果は上書きしないこと。

## 4.2 実行経路の完全追跡

以下を行番号付きで報告すること。

1. 入力波 `a`, `b` の生成箇所
2. 基本波・倍音集合の定義箇所
3. 搬送波シフトの有無
4. 搬送波逆シフトの有無
5. 散乱係数 `r`, `t`, `R`, `T` の生成箇所
6. 散乱更新式の本体
7. 衝突開始条件
8. 反復散乱のループ
9. `normalize` または `normalize_pair` の定義
10. 正規化前後で保存される量、失われる量
11. 出力チャネル配列の生成箇所
12. 反射量、透過量、交換量の計測箇所
13. `B_to_A` の定義
14. `gray_error`, `gray_depth` の定義
15. CSV、JSON、PNG、SVG 等の保存箇所
16. 最終判定が参照する列

## 4.3 データ構造の確定

各状態の実型を明記すること。

例:

```text
scalar complex
1D complex ndarray
2D complex ndarray
harmonic coefficient vector
sampled field vector
channel pair
```

特に、`a` と `b` が、

- 空間サンプル列なのか
- 倍音係数列なのか
- 局在波形そのものなのか
- それらを混在させた表現なのか

を確定すること。

## 4.4 現行散乱核の数学化

実コードから、現行更新を可能な限り厳密に数式へ転記すること。

例えば、

\[
\begin{pmatrix}
a'\\
b'
\end{pmatrix}
=
\mathcal N
\left[
\begin{pmatrix}
r&t\\
t&r
\end{pmatrix}
\begin{pmatrix}
a\\
b
\end{pmatrix}
\right]
\]

の形なら、正規化 \(\mathcal N\) の具体的定義まで含めること。

## 4.5 現行コードが偶奇を読んでいるか

次を実コード上で判定すること。

- FFT または倍音番号を参照しているか
- 偶数モードと奇数モードを区別しているか
- 半周期移動演算子 \(P\) を参照しているか
- 入力波の半周期相関を参照しているか
- 交差項の偶奇を参照しているか
- `r`, `t` が入力波形依存か
- `r`, `t` が外部 `R` のみで決まるか

現行コードが全倍音へ同一係数を作用させているなら、その事実をコード引用と数式で明示すること。

## 4.6 正規化監査

`normalize` が、偶数型と奇数型の散乱差を消している可能性を調べること。

以下を測定する。

```text
norm before scattering
norm after raw scattering
norm after normalization
channel-wise norm
combined norm
phase relation before/after normalization
harmonic parity weight before/after normalization
```

正規化前には存在した差が、正規化後に消えていないかを確認すること。

## 4.7 Stage A の成果物

以下を作成すること。

```text
01_current_code_call_graph.md
02_current_scattering_kernel_exact_form.md
03_state_and_array_schema.md
04_normalization_audit.md
05_current_parity_blindness_report.md
06_candidate_insertion_points.md
```

最後に、

```text
Stage A 完了。コード変更は未実施。
```

と明記して停止すること。

---

# 5. Stage B: 既存コードの再現試験

Stage A の報告を提出した時点で必ず停止する。人間の明示承認後に限り、Stage B の再現試験へ進む。Stage A と Stage B を同一実行内で連続実施してはならない。

## 5.1 最小再現ケース

次の入力を用意する。

### F1: 基本奇数波

\[
\psi_{F1}(u)=\cos u
\]

### FK: 等振幅奇数倍音局在波

\[
\psi_{F,K}(u)
=
\frac1K
\sum_{m=0}^{K-1}
\cos((2m+1)u)
\]

### BK: 等振幅偶数倍音局在波

\[
\psi_{B,K}(u)
=
\frac1K
\sum_{m=1}^{K}
\cos(2mu)
\]

### MIX: 奇偶混合波

\[
\psi_M(u)
=
\sqrt p\,\psi_{F,K}(u)
+
 e^{i\phi}\sqrt{1-p}\,\psi_{B,K}(u)
\]

## 5.2 交絡統制

奇数型・偶数型の比較では、可能な限り以下を揃えること。

```text
sample count
spatial domain
component count
combined norm
peak amplitude
RMS amplitude
localization width
spectral RMS wavenumber
highest wavenumber
carrier treatment
phase origin
```

完全一致できない項目は、差分を数値で報告すること。

## 5.3 現行散乱核の結果

同一 `R`、同一衝突回数、同一正規化条件で、次を比較する。

```text
F x F
B x B
F x B
mixed x mixed
```

出力として最低限、

```text
raw reflected norm
raw transmitted/exchanged norm
normalized reflected norm
normalized transmitted/exchanged norm
p_B before/after
p_F before/after
phase shift
localization width before/after
```

を保存すること。

## 5.4 期待される監査結果

現行散乱係数が入力波形に依存しない場合、偶数型と奇数型に同じ `r,t` が適用されるはずである。

その場合、

```text
現行コードは、型を保存し得るが、型ごとの散乱応答差を生成しない
```

と結論づけること。

ここで、

```text
型が保存されること
```

と、

```text
型によって反射・透過・交換応答が異なること
```

を混同してはならない。

---

# 6. Stage C: 新散乱関数の数学設計

Stage C では、まだ既存本体へ実装しない。

新しい散乱関数の候補を独立モジュールとして設計する。

## 6.1 必要条件

新散乱関数は、少なくとも次を満たすこと。

1. 入力波に人為的な `if boson` / `if fermion` 分岐を置かない
2. 同一数式を全入力へ適用する
3. 入力波の内部構造から応答差が生じる
4. 総保存量を明示する
5. 反射・透過・交換配列を返す
6. 純偶数、純奇数、混合を同一APIで処理する
7. 搬送波を除いたカーネル座標 \(u\) でパリティを評価する
8. 密度 \(|\psi|^2\) ではなく符号付き場または複素場を使う
9. 線形応答と非線形応答を分離する
10. 既存 `R` 掃引との互換性を残すか、非互換なら明示する

## 6.2 パリティ分解

入力波に対して、

\[
\psi_B=\Pi_B\psi
\]

\[
\psi_F=\Pi_F\psi
\]

を計算する。

混合比は、

\[
p_B
=
\frac{\|\psi_B\|^2}
{\|\psi_B\|^2+\|\psi_F\|^2}
\]

\[
p_F
=
\frac{\|\psi_F\|^2}
{\|\psi_B\|^2+\|\psi_F\|^2}
\]

とする。

ただし、これらの値をそのまま条件分岐へ使うだけでは不十分である。

## 6.3 相互作用相関行列

二入力 \(\psi_A,\psi_B\) に対して、候補として次の相互作用相関を計算する。

\[
C_{XY}
=
\left\langle
\Pi_X\psi_A,
\Pi_Y\psi_B
\right\rangle,
\qquad
X,Y\in\{B,F\}
\]

これにより、

\[
\mathbf C
=
\begin{pmatrix}
C_{BB}&C_{BF}\\
C_{FB}&C_{FF}
\end{pmatrix}
\]

を得る。

さらに、位置依存相関として、

\[
C_\pi^{\mathrm{raw}}[\psi]
=
\int
\psi(u)^*\psi(u+\pi)\,du
\]

を計算する。

純偶数型では、

\[
C_\pi^{\mathrm{raw}}[\psi_B]=+\|\psi_B\|^2
\]

純奇数型では、

\[
C_\pi^{\mathrm{raw}}[\psi_F]=-\|\psi_F\|^2
\]

となる。ここで \(C_\pi^{\mathrm{raw}}\) はノルム二乗と同じ次元を持つ非正規化量である。後述の無次元正規化指標 \(c_\pi\in[-1,1]\) と同一記号で扱ってはならない。

ただし、これだけで散乱係数を決めてよいとは限らない。

候補関数を複数設計し、比較すること。

## 6.4 散乱関数候補

最低でも、次の四系列を設計すること。

### Candidate 0: 現行基準

\[
\mathbf U_0(R)
=
\begin{pmatrix}
r(R)&t(R)\\
t(R)&r(R)
\end{pmatrix}
\]

入力非依存の基準モデル。

### Candidate 1: 半周期相関依存

\[
\mathbf U_1
=
\mathbf U
\left(
R,
c_\pi[\psi_A],
c_\pi[\psi_B],
c_\pi^{\mathrm{sym}}[\psi_A,\psi_B]
\right)
\]

ここで単波の半周期相関は、

\[
c_\pi[\psi]
=
\frac{\operatorname{Re}\langle\psi,P\psi\rangle}{\|\psi\|^2}
\]

とし、二波の交差半周期相関は、

\[
c_\pi^{\mathrm{sym}}[\psi_A,\psi_B]
=
\frac{\operatorname{Re}
\left(
\langle\psi_A,P\psi_B\rangle
+\langle\psi_B,P\psi_A\rangle
\right)}
{2\|\psi_A\|\|\psi_B\|}
\]

と定義する。ゼロノルム入力は明示的にエラーまたは未定義扱いとする。同一式だが、内部半周期相関により応答が変わる候補。

### Candidate 2: 双線形パリティ相関依存

\[
\mathbf U_2
=
\mathbf U
\left(
R,
\mathbf C_K
\right)
\]

とする。標準的な全周期内積では直交射影により \(C_{BF}=C_{FB}=0\) となるため、Candidate 2 では必ず相互作用作用素または局所重みを挿入し、

\[
C^{(K)}_{XY}
=
\left\langle
\Pi_X\psi_A,
K_{\mathrm{int}}\Pi_Y\psi_B
\right\rangle
\]

を用いる。\(K_{\mathrm{int}}=I\) は退化基準として保存し、非自明候補には、局所窓、位置依存結合核、接触領域射影、または畳み込み核を用いる。核を任意調整して結論を作ってはならず、全候補とパラメータを事前列挙すること。

### Candidate 3: モード積・交差項依存

奇数×奇数、偶数×偶数、奇数×偶数の双線形積が、

\[
F\times F\to B
\]

\[
B\times B\to B
\]

\[
B\times F\to F
\]

という \(\mathbb Z_2\) 合成則を持つことを使い、散乱生成項を構成する候補。

## 6.5 禁止事項

以下は禁止する。

```python
if p_F > threshold:
    reflect_more()
else:
    transmit_more()
```

```python
if is_boson:
    reflection = 0
```

```python
if is_fermion:
    apply_pauli_recoil()
```

これは仮説をコードへ直接作り込むだけであり、検証にならない。

## 6.6 保存条件

候補散乱関数は、少なくとも全体ノルム保存を満たすか、破る場合はその理由を明記すること。

\[
\|\psi_A'\|^2+\|\psi_B'\|^2
=
\|\psi_A\|^2+\|\psi_B\|^2
\]

線形ユニタリ形式を維持できない場合、非線形保存則を別途定義すること。

また、正規化処理を散乱関数の外に置くか内に置くかを明示すること。

---

# 7. 新しい散乱配列のAPI設計

新関数は、最低限次の情報を返すこと。

```python
@dataclass
class ScatteringResult:
    a_raw: np.ndarray
    b_raw: np.ndarray
    a_out: np.ndarray
    b_out: np.ndarray

    # 経路別振幅配列。合成後の出力だけから反射・交換を逆算しない。
    path_a_to_a_amplitude: np.ndarray
    path_b_to_a_amplitude: np.ndarray
    path_b_to_b_amplitude: np.ndarray
    path_a_to_b_amplitude: np.ndarray

    # 経路単独ノルム。干渉項を含めない。
    path_a_to_a_norm: float
    path_b_to_a_norm: float
    path_b_to_b_norm: float
    path_a_to_b_norm: float

    # 出力チャネル内の交差干渉項と合成出力ノルム。
    interference_in_a: float
    interference_in_b: float
    output_a_norm: float
    output_b_norm: float

    # 互換的な読み出し名。必ず経路単独ノルムの別名とし、干渉項を配分しない。
    reflection_a: float
    reflection_b: float
    exchange_a_to_b: float
    exchange_b_to_a: float

    p_boson_in_a: float
    p_fermion_in_a: float
    p_boson_in_b: float
    p_fermion_in_b: float

    p_boson_out_a: float
    p_fermion_out_a: float
    p_boson_out_b: float
    p_fermion_out_b: float

    # 非正規化半周期相関 C_pi^raw。複素数のまま保持する。
    parity_correlation_raw_a: complex
    parity_correlation_raw_b: complex

    # 正規化無次元パリティ指標 c_pi = Re<C,P C>/||C||^2。
    # 純偶数型で +1、純奇数型で -1、混合状態では原則として [-1, 1] 内を取る。
    parity_indicator_a: float
    parity_indicator_b: float

    cross_parity_matrix: np.ndarray

    norm_before: float
    norm_after_raw: float
    norm_after_final: float

    diagnostics: dict
```

候補関数のインターフェース例:

```python
def scatter_wave_pair(
    a: np.ndarray,
    b: np.ndarray,
    *,
    reflection_parameter: float,
    coordinate: np.ndarray,
    carrier_mode: int | None,
    model: str,
    normalize_output: bool,
) -> ScatteringResult:
    ...
```

APIは既存コードへ差し込みやすいように設計すること。

---

# 8. 検証実験

## 8.1 単体テスト

以下を実施する。

### パリティ射影

\[
\Pi_B^2=\Pi_B
\]

\[
\Pi_F^2=\Pi_F
\]

\[
\Pi_B\Pi_F=0
\]

\[
\Pi_B+\Pi_F=I
\]

### 対蹠反応

\[
\|\psi_F+P\psi_F\|\approx0
\]

\[
\|\psi_B+P\psi_B\|\approx2\|\psi_B\|
\]

これは物理予言ではなく実装確認である。

### 搬送波逆シフト

搬送波シフト前後で、カーネル座標上の偶奇分類が一致すること。

## 8.2 基準比較

Candidate 0〜3について、同じ入力群で比較する。

```text
F1 x F1
FK x FK
BK x BK
FK x BK
MIX x MIX
```

## 8.3 必須出力

各候補について、

```text
reflection
transmission/exchange
sector weights
sector conversion
phase shifts
localization width
norm conservation
commutator-like parity response
```

をCSVとMarkdown表で保存すること。

## 8.4 判定基準

良い候補は、

1. 共通数式である
2. 入力構造により応答差が出る
3. 保存則が明確である
4. 数値安定である
5. 既存結果を不必要に破壊しない
6. 偶数型・奇数型の差が単なるラベルではない
7. 交絡統制後も差が残る

こと。

---

# 9. 既存コードへの実装禁止

Stage C まで終了しても、既存 System A / System B 本体へ新散乱関数を組み込んではならない。

最終成果物を提出し、次の承認を待つこと。

```text
承認待ち:
- どの候補散乱関数を採用するか
- 保存則をどれにするか
- 正規化位置をどこにするか
- 既存R掃引を維持するか
- System AのみかSystem Bにも適用するか
```

---


# 9.1 実行前の確定事項

本節は、指示解釈上の曖昧さを除くための確定仕様である。

## 9.1.1 Stage A 後の停止

最初の実作業では Stage A のみを実施する。

Stage A の監査報告を提出した時点で必ず停止し、人間の明示承認を待つこと。Stage B、Stage C を自動開始してはならない。

## 9.1.2 監査範囲

Stage A の主対象は System A とする。

System B は、System A と散乱核、係数生成、正規化、または評価関数を共有している箇所を比較参照する。System B 全体を System A と同じ深さで監査する必要はない。ただし、共有核の実装差、同名関数の別定義、または System A の挙動解釈に影響する差が見つかった場合は、その範囲に限って追跡する。

## 9.1.3 反射量、交換量、干渉項の分離

出力ノルムを反射量と交換量へ一意に分解してはならない。

現行形

\[
a_{\mathrm{out}}=ra+tb
\]

に対しては、最低限、経路別振幅と干渉項を分離して返す。

\[
A_{a\to a}=ra,
\qquad
A_{b\to a}=tb
\]

\[
N_{a\to a}=\|ra\|^2,
\qquad
N_{b\to a}=\|tb\|^2
\]

\[
I_a=2\operatorname{Re}\langle ra,tb\rangle
\]

したがって、

\[
\|a_{\mathrm{out}}\|^2
=N_{a\to a}+N_{b\to a}+I_a
\]

である。B 側も同様に返す。

API は、少なくとも次を別項目として保持すること。

```text
path_a_to_a_amplitude
path_b_to_a_amplitude
path_a_to_a_norm
path_b_to_a_norm
interference_in_a
output_a_norm
path_b_to_b_amplitude
path_a_to_b_amplitude
path_b_to_b_norm
path_a_to_b_norm
interference_in_b
output_b_norm
```

`reflection_a` や `exchange_b_to_a` という単一スカラーを返す場合、それは経路単独ノルムであることを明記し、干渉項を含めない。干渉項を任意配分してはならない。

## 9.1.4 commutator-like parity response の定義

散乱写像を \(S\)、二チャネル双方へ作用する半周期演算子を

\[
\mathcal P=\operatorname{diag}(P,P)
\]

とする。

線形または線形化可能な候補では、状態依存指標を

\[
\eta_{[S,P]}(\Psi)
=
\frac{\|S(\mathcal P\Psi)-\mathcal P S(\Psi)\|}
{\|\Psi\|}
\]

と定義する。

非線形写像でも同じ有限差分定義を用いる。これは抽象作用素ノルムではなく、指定入力に対する状態依存非可換応答である。必要に応じ、複数入力集合上の最大値・平均値も併記する。

さらに、散乱前後のパリティ重量変化

\[
\Delta p_B=p_B'-p_B,
\qquad
\Delta p_F=p_F'-p_F
\]

を別指標として返し、非可換応答と混同しないこと。

## 9.1.5 成果物配置先

既存成果物との混在を避けるため、成果物の親ディレクトリを次に固定する。別の場所へ自動変更してはならない。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/codex_systemA_scattering_audit_v1/
```

指定パスが存在しない場合は、親ディレクトリの実在位置を検索して候補を報告し、書込みを開始せず停止する。

その下を次のように分ける。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/codex_systemA_scattering_audit_v1/
  stage_A_audit/
  stage_B_reproduction/
  stage_C_candidates/
  reports/
  data/
  figures/
  logs/
  manifests/
```

最初の実行では `stage_A_audit/`, `reports/`, `logs/`, `manifests/` のみ使用する。原本コード・原本データ・既存結果へ書き込んではならない。

## 9.1.6 変更記録

Stage A では原本変更ゼロを証明するため、対象原本のパス、サイズ、更新時刻、SHA-256 を実行前後で保存する。

```text
manifests/source_manifest_before.json
manifests/source_manifest_after.json
```

差があれば Stage A を失敗として停止する。


# 10. 最終成果物

以下を作成すること。

```text
00_audit_scope_and_source_hashes.md
01_current_code_call_graph.md
02_current_scattering_kernel_exact_form.md
03_state_and_array_schema.md
04_normalization_audit.md
05_current_parity_blindness_report.md
06_existing_behavior_reproduction.md
07_new_scattering_function_requirements.md
08_candidate_scattering_models.md
09_candidate_api_design.md
10_candidate_unit_tests.md
11_candidate_comparison_results.md
12_recommended_next_step.md
```

コード成果物:

```text
audit_current_scattering.py
parity_projection.py
candidate_scattering_models.py
test_parity_projection.py
test_candidate_scattering_models.py
run_candidate_comparison.py
```

データ成果物:

```text
current_behavior_baseline.csv
candidate_comparison.csv
candidate_diagnostics.json
```

図:

```text
reflection_comparison.png
exchange_comparison.png
sector_weight_comparison.png
norm_error_comparison.png
localization_width_comparison.png
```

---

# 11. 報告時の必須区分

報告書では、各記述を必ず以下へ分類すること。

```text
[コード上の事実]
[数学的帰結]
[モデル定義]
[作業仮説]
[数値観測]
[未導出]
[実装候補]
[棄却候補]
```

仮説を事実として書いてはならない。

既存コードの挙動を、仮説に都合よく再解釈してはならない。

---

# 12. 最後の停止条件

次の一文で終了すること。

```text
現行散乱コードの監査と新散乱関数候補の独立設計まで完了した。
既存System A / System B本体への組込みは未実施であり、人間の承認を待つ。
```

