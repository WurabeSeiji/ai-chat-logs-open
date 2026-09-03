# Claude Code 追加実験指図：初期シード精度と時間発展精度の切り分け（64bit vs 100桁）

作成日: 2026-09-02  
対象実験: `N7_N8_DeltaTau_denominator_wide_sweep_20260902`  
目的: 現在観測されている seedless インフレーションについて、**初期状態に含まれる float64 級の微小ずれ**と、**時間発展演算中に毎 step 新たに入る float64 丸め誤差**を、条件を一度に一つだけ変えて切り分ける。

---

## 0. 最重要原則

この追加実験では、**初期状態の精度と時間発展の精度を同時に変更して比較してはならない。**

比較は必ず次の順で行う。

```text
A: IC64  + Dynamics64    現行基準
B: IC64  + Dynamics100   初期値は A と完全同一、時間発展だけ100桁
C: IC100 + Dynamics100   時間発展は B と完全同一、初期値だけ100桁生成
```

比較の意味は次の通り。

```text
A vs B : 時間発展演算中の有限精度誤差の寄与
B vs C : 初期状態に含まれる有限精度シードの寄与
```

**A→B では初期値を1 bitたりとも変えない。**  
**B→C では時間発展アルゴリズム、N、D、振幅規約、観測量、停止条件を変えない。**

---

## 1. 実験の順序

広域 D sweep は行わない。まず精度問題だけを切り分ける。

### Primary

```text
N = 8
D = 8
DeltaTau = 2*pi/8
```

### Replication

Primary が完了した後にのみ、同じ A/B/C を

```text
N = 7
D = 7
DeltaTau = 2*pi/7
```

で再現確認する。

**N=7 と N=8 を同時に混ぜて結論を出さない。まず N=8 単独で A/B/C の因果切り分けを完了し、その後 N=7 を独立追試として扱う。**

D=124 や他の D はこの実験では走らせない。精度問題が確定してから別実験とする。

---

## 2. A: IC64 + Dynamics64（同環境基準再走）

既存の親 `parent_v.npz` を使用する。

- N=8: 既存 sweep が参照したものと同一 SHA256
- N=7: 同上

状態型:

```text
numpy.complex128
numpy.float64
```

力学は既存正本と完全同一。

\[
H_{ef}(z)=A_{ef}\,\overline{z_e}z_f
\]

\[
z' = \exp\left[-i\left(\frac{2\pi}{D}\right)H(z)\right]z
\]

各 step:

```python
H = A * (conj(z)[:, None] * z[None, :])
w, V = np.linalg.eigh(H)
z = V @ (exp(-1j*(2*pi/D)*w) * (V.conj().T @ z))
```

- seed 追加なし
- clipping なし
- renormalization なし
- 状態の丸め直しなし

A は既存結果を参照するだけでなく、**今回の同じ実行環境で再走**し、既存 sweep の N=8,D=8 / N=7,D=7 と初期状態 SHA256、主要時系列、onset が再現するか確認する。

A が既存正本と再現しない場合は B/C を続行せず、環境差として停止・報告する。

---

## 3. B: IC64 + Dynamics100

### 3.1 最重要条件：IC64 は「同じ数値」ではなく「同じ binary64 値」を保存する

A で読み込んだ `complex128` の各実部・虚部を、その IEEE-754 binary64 値と**完全に同じ有理数**として100桁環境へ持ち上げる。

禁止:

- 100桁の理論式から再生成する
- `str(float)` の短い10進表示を「元値」として使う
- 位相を100桁で再計算する
- norm を100桁で再調整する

推奨方法:

各 binary64 実数 `x` について Python の

```python
num, den = float(x).as_integer_ratio()
```

を使い、100桁側で

```python
mp.mpf(num) / mp.mpf(den)
```

として**exact binary64 lift**する。

複素数は実部・虚部を別々に exact lift する。

これにより B の初期状態は数学的に A の binary64 配列と完全同一であり、違うのは以後の演算精度だけになる。

### 3.2 100桁力学

任意精度ライブラリを使用し、**少なくとも decimal precision = 100 digits** を確保する。

推奨:

```text
mpmath mp.dps = 100
```

ただし、Hermitian 固有値分解が100桁のまま実行されることを開始時 self-test で確認すること。

- NumPy / SciPy の float64/complex128 へ途中で落とさない
- `numpy.linalg.eigh` を呼ばない
- 固有値・固有ベクトル・位相指数・行列積のどこにも binary64 を混入させない
- plotting / CSV 出力のための変換は観測後のコピーにだけ許可し、次 step の状態には戻さない

100桁 Hermitian eigensolver が利用できない場合、**別アルゴリズムへ勝手に変更して実験を継続しない。** 使用予定ライブラリと代替案を報告して停止する。

### 3.3 π も100桁

B の初期状態は IC64 の exact lift だが、時間発展側の

\[
\Delta\tau=2\pi/D
\]

に使う π は100桁で評価する。

A と B の差には「発展演算の高精度化」がすべて含まれるため、これは意図した条件差である。

---

## 4. C: IC100 + Dynamics100

B と同じ100桁力学をそのまま使う。

変更してよいのは**初期状態の生成精度だけ**。

### 4.1 初期 hm を100桁で再生成

N=8 / N=7 の hm 位相構成を、既存の `hm_mp_free_N3_N40_20260901` の設計規則から読み取り、100桁で再生成する。

- N=8: even N の既存 1-factor phase construction
- N=7: odd N の既存 cyclic-distance phase construction

重要:

- 新しい物理仮定を入れない
- 位相規則を変更しない
- 頂点・辺順序を変更しない
- A/B と同じ edge ordering を使う
- 振幅スケールは既存親と同じ数学的規約を再現する
- A/B の float64 parent を「近似目標」として fit しない

C は既存の解析的設計から100桁で直接構築する。

### 4.2 振幅規約の注意

N=7,8 は既存 historical-compatibility 領域に属するため、振幅スケールの扱いを勝手に共通 `r^2=1/15` へ変更しない。

既存親が持つ振幅規約をコードから特定し、その数学的規約を100桁で再現する。

もし既存親の振幅が frozen legacy norm の数値値しか定義されておらず、解析式が存在しない場合は、

1. その事実を明記する
2. A/B と同じ振幅スケールを exact lift した版を C の振幅として使う
3. 位相のみ100桁解析生成する

とし、振幅まで新規推定しない。

---

## 5. 実行長

### N=8,D=8

A/B/C ともまず

```text
STEPS = 2000
```

まで走らせる。

### N=7,D=7

同じく

```text
STEPS = 2000
```

### 自動延長

C が step 2000 まで primary onset

\[
H_\perp/H > 0.05
\]

に到達しない場合、状態を checkpoint し、

```text
4000 -> 8000
```

まで段階的に延長してよい。

ただし B と C の比較では同一 step / 同一 τ の観測値を必ず残す。

100桁計算が極端に遅い場合は、力学を変更して高速化せず、checkpoint を保存して停止・報告する。

---

## 6. 観測量

A/B/C すべてで同一定義を使う。

初期 z0 から

```python
p = Re(z0) / ||Re(z0)||
q = Im(z0) - <Im(z0),p> p
q = q / ||q||
```

を各 run 自身の精度で一度だけ作る。

主観測量:

\[
f_\perp = \frac{H_\perp}{H}
\]

その他:

```text
step
tau = step * 2*pi/D
Hperp_frac
Hparallel_frac
H_total
H_total_rel_drift
global_closure = abs(z@z)/H
PR
PR_over_M
amp_min
amp_max
amp_std
```

高精度値を CSV に保存するとき、float64 へ丸めない。

最低でも科学表記で **110 significant digits 程度**を文字列保存する。

---

## 7. 特に記録すべき「seed floor」

A/B/C について以下を記録する。

```text
initial_Hperp_frac
step1_Hperp_frac
minimum_Hperp_frac_before_growth
minimum_Hperp_step
first_step_above_1e-180
first_step_above_1e-150
first_step_above_1e-120
first_step_above_1e-90
first_step_above_1e-60
first_step_above_1e-30
first_step_above_1e-12
first_step_above_1e-6
first_step_above_1e-3
first_step_above_0.05
```

A/B では 1e-180 等は当然 NA になり得る。

C では100桁計算により Hperp/H が非常に小さくなるため、**log10(Hperp/H)** を高精度で直接計算し、underflow を起こす float 変換をしない。

---

## 8. 成長率 fit

### 8.1 共通 fit 窓

A/B/C 間で比較可能な範囲として、可能なら

\[
10^{-20} \le f_\perp \le 10^{-6}
\]

で

\[
\ln f_\perp = a + \gamma_{\rm step}\,step
\]

\[
\ln f_\perp = a_\tau + \gamma_\tau\,\tau
\]

を fit する。

### 8.2 C 専用の深部 fit

100桁 seed の初期増幅を見るため、C では追加で

\[
10^{-160} \le f_\perp \le 10^{-40}
\]

など、実データに応じた深部指数域を fit する。

窓は結果をよく見せるために恣意的に選ばず、使用した範囲・点数・R²を全て保存する。

### 8.3 最重要比較

```text
gamma_tau_A
gamma_tau_B
gamma_tau_C
```

を比較する。

特に、B と C で onset が大きく違っても

\[
\gamma_{\tau,B} \approx \gamma_{\tau,C}
\]

なら、初期 seed の大きさだけが onset を移動させ、増幅率は同じ力学が支配している可能性が高い。

---

## 9. 判定ロジック

### Case 1: A ≈ B、C の onset だけ大幅に遅れる

最も単純な解釈:

```text
初期 IC64 seed が主要な発火源
時間発展中の新規 float64 丸め注入は二次的
100桁 IC では同じ不安定力学がより小さい seed から成長
```

ただし `gamma_B ≈ gamma_C` を確認するまで断定しない。

### Case 2: A と B が大きく違う

```text
時間発展演算中の float64 丸め誤差が軌道 / onset に重要
```

この場合 B/C 比較だけで「初期 seed が原因」と結論しない。

### Case 3: B ≈ C

```text
初期 IC64 seed の寄与は小さく、発展演算または別の決定論的ずれが支配
```

ただし C の初期残差、closure、self-consistency residual を精査する。

### Case 4: C が100桁でも最終的に発火し、growth rate が B と一致

これは

```text
有限精度は seed の大きさを決めるが、指数不安定性そのものは力学側に存在
```

という解釈を強く支持する。

### Case 5: C が十分長時間でも指数増幅しない

```text
float64 で観測した巨視化が有限精度依存現象である可能性
```

を検討する。ただし「十分長時間」の根拠を growth upper bound とともに示す。

---

## 10. 100桁初期状態の資格審査

C を走らせる前に、IC100 について以下を高精度で記録する。

```text
N
M_edge
norm
mean_amp2
global_closure
local_closure_max_abs
H_eigen_residual = ||H(v)v - mu v|| / ||v||
mu
H_hermiticity_error
```

可能なら theoretical mu と比較する。

N>=4 の hm について既存解析どおり

\[
H(v)v = -2r^2 v
\]

が成立する設計なら、その residual を100桁で確認する。

**この資格審査に失敗した IC100 を「高精度 seed」として使用しない。**

---

## 11. 数値実装監査

100桁 run の開始時に、実際の backend が100桁であることを証明する簡単な self-test を保存する。

例:

```text
precision_digits_requested
precision_digits_effective
pi_string_110digits
1 + 10^-80 != 1  が True
1 + 10^-110 == 1 になるか否か（実効精度確認）
Hermitian eigensolver backend
```

さらに

```text
Python version
mpmath / arbitrary-precision library version
platform
CPU architecture
```

を保存する。

NumPy 配列への変換は、図化専用のコピー以外は禁止。

---

## 12. 必須図

### Fig P1: N=8 A/B/C 全時系列

- x: τ
- y: Hperp/H
- y log scale
- A/B/C overlay

### Fig P2: N=8 初期増幅域

- x: τ
- y: log10(Hperp/H)
- A/B/C
- onset までの直線域を見やすくする

### Fig P3: N=8 平行移動比較

B と C の指数域について、onset または fit intercept 差だけ横方向に平行移動したとき slope が一致するかを可視化する。

**平行移動後の図は補助図であり、原データ図を置き換えない。**

### Fig P4-P6: N=7 replication

N=8 と同じ3図。

### Fig P7: onset vs initial seed floor

A/B/C および N=7/8 の点について、

```text
x = -log10(initial effective seed floor)
y = onset_tau_0p05
```

を散布図化する。

点数が少ないため回帰を「法則」と呼ばない。

---

## 13. 出力先

既存実験フォルダ内に新規サブフォルダを作る。

```text
precision_seed_dynamics_isolation_100digit_20260902/
```

構成:

```text
precision_seed_dynamics_isolation_100digit_20260902/
├── README.md
├── ANALYSIS.md
├── RUN_METADATA.json
├── SHA256SUMS.txt
├── program/
│   ├── run_float64_control.py
│   ├── run_mp100_same_ic64.py
│   ├── build_ic100.py
│   ├── run_mp100_ic100.py
│   ├── analyze_precision.py
│   └── plot_precision.py
├── data/
│   ├── N8_D8/
│   │   ├── A_IC64_DYN64/
│   │   ├── B_IC64_DYN100/
│   │   └── C_IC100_DYN100/
│   └── N7_D7/
│       ├── A_IC64_DYN64/
│       ├── B_IC64_DYN100/
│       └── C_IC100_DYN100/
├── results/
│   ├── precision_summary.csv
│   ├── growth_fits.csv
│   └── qualification_ic100.csv
└── figures/
    ├── figP1_N8_ABC_tau.png
    ├── figP2_N8_ABC_log_growth.png
    ├── figP3_N8_shifted_growth.png
    ├── figP4_N7_ABC_tau.png
    ├── figP5_N7_ABC_log_growth.png
    ├── figP6_N7_shifted_growth.png
    └── figP7_onset_vs_seed_floor.png
```

既存 sweep のファイルは読み取り専用。

---

## 14. ANALYSIS.md で必ず答える問い

1. A は既存 float64 sweep を再現したか。
2. A と B の初期状態は exact binary64 lift により数学的に同一か。その検証方法は何か。
3. A と B の onset、growth rate、飽和状態はどの程度一致／不一致か。
4. B と C の初期 Hperp/H、closure、self-consistency residual は何桁違うか。
5. B と C の onset 差は何 step / 何 τ か。
6. B と C の `gamma_tau` は一致するか。
7. C の初期 seed floor から巨視的 onset まで、単一指数則で説明できる範囲はどこか。
8. インフレーションの「seed」と「不安定性そのもの」を数値的に分離できたと言えるか。
9. N=7 replication は N=8 の結論を再現するか。
10. 100桁にして初めて見えた新しい数値的不安定性、solver 依存性、残差問題はあるか。

---

## 15. 完了条件

以下を全て満たすまで完了としない。

- N=8,D=8 の A/B/C が完了
- A の既存結果再現確認済み
- B の IC64 exact lift を `as_integer_ratio()` 等で検証済み
- B の全時間発展が100桁から一度も float64 に落ちていない
- C の IC100 資格審査済み
- B/C の growth rate と onset 差を定量化
- 必須図 P1-P3 完成
- N=7,D=7 replication 完了
- 必須図 P4-P7 完成
- `ANALYSIS.md` の10問回答済み
- SHA256 / environment / precision self-test 保存済み
- 既存正本を一切変更していない

---

## 16. Claude Code への最終指示

最初に実装予定の arbitrary-precision backend と Hermitian eigensolver が**本当に complex Hermitian を100桁で処理するか**を小行列で self-test すること。

その後、N=8,D=8 の

```text
A -> B -> C
```

をこの順に完了し、A/B/C の初期状態比較・onset・growth rate を中間報告として `ANALYSIS.md` に記録する。

**N=8 の切り分けが成立する前に N=7 を走らせない。**

N=8 完了後、同じコード・同じ判定条件で N=7,D=7 を replication として実行する。

最重要の判定は次である。

> **現在の float64 seedless インフレーションにおいて、巨視化を始動する微小 seed は「初期 IC64 に既に含まれる誤差」なのか、「時間発展中に毎 step 注入される有限精度誤差」なのか。そして seed の大きさを100桁級まで下げても、同じ growth rate を持つ指数不安定性が残るのか。**
