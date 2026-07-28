# 06 現行System A散乱核の独立再現 — Stage B報告

## 1. 実施境界

**[モデル定義]** Stage Bでは、監査済みの式を専用コードへ独立転記した。既存System A／System B原本はimport・実行・変更せず、出力先にも使用していない。

実行した独立コード:

```text
stage_B_reproduction/reproduce_current_system_A_stage_B.py
stage_B_reproduction/verify_stage_B_outputs.py
stage_B_reproduction/build_stage_B_source_manifest.py
```

**[数値観測]** 事前・事後マニフェストでSystem A本体、散乱源、System B比較コードの3原本を比較した。

```text
source_count_before = 3
source_count_after  = 3
unchanged           = true
changed_path_count  = 0
```

原本に差はない。

## 2. 固定再現条件

**[モデル定義]**

| 項目 | 値 |
|---|---:|
| \(\chi\) サンプル数 | 512 |
| \(\eta\) サンプル数 | 16 |
| 空間領域 | \([-\pi,\pi)\) |
| \(R_{\mathrm{input}}\) | 0.6971778791282474 |
| 更新回数 | 32 |
| A搬送波・\(\eta\)モード | \(q_A=+1,\ m_A=1\) |
| B搬送波・\(\eta\)モード | \(q_B=-1,\ m_B=2\) |
| FK/BK成分数 | \(K=4\) |
| MIX | \(p=0.5,\ \phi=\pi/2\) |

入力は次の5ケースとした。

```text
F1_x_F1
FK_x_FK
BK_x_BK
FK_x_BK
MIX_x_MIX
```

\[
\psi_{F1}(u)=\cos u,
\]

\[
\psi_{F,K}(u)=\frac14\sum_{m=0}^{3}\cos((2m+1)u),
\]

\[
\psi_{B,K}(u)=\frac14\sum_{m=1}^{4}\cos(2mu),
\]

\[
\psi_M(u)=\frac{1}{\sqrt2}\psi_{F,K}(u)
+\frac{i}{\sqrt2}\psi_{B,K}(u).
\]

## 3. 交絡統制

**[数値観測]** 主比較FK×BKでは次を一致させた。

| 項目 | FK | BK | 差 |
|---|---:|---:|---:|
| サンプル数 | 512 | 512 | 0 |
| 成分数 | 4 | 4 | 0 |
| 正規化カーネルノルム二乗 | 1 | 1 | \(2.22\times10^{-16}\)以下 |
| 正規化ピーク振幅 | 0.125 | 0.125 | 0 |
| RMS振幅 | 0.0441941738242 | 0.0441941738242 | 0 |
| 位相原点 | 0 | 0 | 0 |

初期A/Bの各状態ノルム二乗は1、結合ノルム二乗は2である。搬送波処理はA/Bの向きを除いて同一規則とした。

完全一致しない量は次のとおり。

| 項目 | FK | BK | 差 |
|---|---:|---:|---:|
| IPR局在幅 | 1.16896470831 | 1.35852655290 | 0.189561844591 |
| スペクトルRMS波数 | 4.58257569496 | 5.47722557505 | 0.894649880096 |
| 最高波数 | 7 | 8 | 1 |

**[未導出]** この3差は入力波形の構成差であり、散乱応答差ではない。したがって局在幅やスペクトル類似度の出力差だけから、型依存散乱を主張することはできない。

## 4. 現行核と位相

**[コード上の事実]** 独立コードへ転記した現行核は、

\[
\delta=2\arcsin\sqrt{R},
\quad
r=-ie^{i\delta/2}\sin\frac{\delta}{2},
\quad
t=e^{i\delta/2}\cos\frac{\delta}{2}
\]

および、

\[
\widetilde a=ra+tb,
\qquad
\widetilde b=rb+ta
\]

である。最後に各出力チャネルを独立に単位ノルムへ正規化した。

**[数値観測]**

```text
delta        = 1.9761630154079322
r            = 0.6971778791282475 - 0.45947892659238054 i
t            = 0.30282212087175253 + 0.45947892659238054 i
arg(r)       = -0.5827148190909305
arg(t)       =  0.9880815077039661
arg(t/r)     =  pi/2
R = |r|^2    = 0.6971778791282475
T = |t|^2    = 0.30282212087175253
```

この係数位相は全入力・全衝突で同一であり、入力パリティによる位相分岐はない。各状態の初期A/Bへの射影位相は衝突別CSVに別列で保存した。

## 5. 経路振幅・ノルム・干渉

**[コード上の事実]** 各衝突について、次の完全な複素配列をケース別NPZへ保存した。

```text
path_a_to_a = r*a
path_b_to_a = t*b
path_b_to_b = r*b
path_a_to_b = t*a
```

各配列の形状は、

```text
(32 collisions, 512 chi samples, 16 eta samples)
```

である。また、

```text
interference_in_a_density
interference_in_b_density
```

を同じ形状で保存した。スカラー干渉項は密度配列の全和である。

**[数値観測]** 5ケース×32衝突のすべてで、

```text
path_a_to_a_norm_raw ≈ R
path_b_to_a_norm_raw ≈ T
path_b_to_b_norm_raw ≈ R
path_a_to_b_norm_raw ≈ T
```

となった。4経路ノルムのケース間最大幅は \(3.89\times10^{-15}\) 以下である。

```text
最大 |干渉項|                    = 1.50e-15
rawチャネルノルム二乗の1からの誤差 = 2.22e-15
正規化後ノルム二乗の1からの誤差    = 4.44e-16
経路分解恒等式の誤差               = 2.11e-15
```

**[数学的帰結]** 現行設定では、A/Bを識別する \(m_A=1,m_B=2\) の\(\eta\)モードが直交するため、二経路の干渉項は数値丸め範囲でゼロである。rawチャネルノルムが既に1なので、チャネル別正規化後の経路ノルムもraw値と一致する。

## 6. 逆搬送波とパリティ量

**[モデル定義]** 反復後の各チャネルは、A由来の \(q=+1,m=1\) 成分とB由来の \(q=-1,m=2\) 成分を含み得る。出力チャネルへ単一の逆搬送波を掛けることはできないため、直交する\(\eta\)モードで二由来成分を射影し、それぞれの \(q\) を除去してから再結合した。

この由来別逆搬送波の再構成残差は最大 \(4.87\times10^{-16}\) だった。

各入力、raw出力、正規化後出力について、

\[
C_\pi^{\mathrm{raw}}=\langle\psi,P\psi\rangle,
\qquad
c_\pi=\frac{\operatorname{Re}C_\pi^{\mathrm{raw}}}{\|\psi\|^2},
\]

\[
p_B=\frac{\|\Pi_B\psi\|^2}{\|\psi\|^2},
\qquad
p_F=\frac{\|\Pi_F\psi\|^2}{\|\psi\|^2}
\]

を衝突別CSVへ保存した。`C_pi_raw_real` と `C_pi_raw_imag` は分離している。

**[数値観測]**

| ケース | 保存結果 |
|---|---|
| F1×F1 | 最大偶数漏れ \(4.45\times10^{-32}\) |
| FK×FK | 最大偶数漏れ \(4.49\times10^{-31}\) |
| BK×BK | 最大奇数漏れ \(5.34\times10^{-31}\) |
| MIX×MIX | \(p_B=p_F=0.5\) からの最大誤差 \(1.55\times10^{-15}\) |
| FK×BK | 二チャネル合計の各セクター保存誤差 \(5.11\times10^{-15}\) |

\[
c_\pi=p_B-p_F
\]

の数値誤差は最大 \(6.66\times10^{-15}\) だった。

半周期移動との有限差分可換応答は、

```text
max eta_[S,P] = 3.07e-16
```

だった。

## 7. 型保存と型依存応答差の分離

**[数学的帰結]** 型保存について確定したのは次である。

- 純奇数×純奇数は、各出力チャネルでも純奇数を保つ。
- 純偶数×純偶数は、各出力チャネルでも純偶数を保つ。
- 50:50混合×50:50混合は、各出力チャネルでも50:50を保つ。
- 純奇数×純偶数では、偶奇セクターは相互変換されないが、セクターを担う振幅がA/Bチャネル間を移動するため、各出力チャネルは混合になる。二チャネル全体の偶奇重量は保存される。

純奇数×純偶数の例:

| 衝突 | A出力 \(p_B\) | A出力 \(p_F\) | B出力 \(p_B\) | B出力 \(p_F\) |
|---:|---:|---:|---:|---:|
| 1 | 0.3028221209 | 0.6971778791 | 0.6971778791 | 0.3028221209 |
| 2 | 0.8444835359 | 0.1555164641 | 0.1555164641 | 0.8444835359 |
| 8 | 0.9974347043 | 0.0025652957 | 0.0025652957 | 0.9974347043 |
| 32 | 0.0405204292 | 0.9594795708 | 0.9594795708 | 0.0405204292 |

**[数値観測]** 型依存応答差については、全ケースの各衝突で同じ \(r,t\) が作用し、4経路ノルムのケース間差は最大 \(3.89\times10^{-15}\) 以下だった。

したがって、

```text
現行コードは、型を保存し得るが、型ごとの散乱応答差を生成しない。
```

型保存は \([S,P]=0\) によるセクター非変換であり、型依存散乱応答差は \(r,t\) または経路ノルムが入力型によって変わることである。本結果では前者だけが成立し、後者は成立しない。

## 8. `B_to_A_transfer` の意味境界

**[コード上の事実]** 全CSV行、全NPZ、診断JSON、検証JSON、マニフェストには次の意味注記を保存した。

```text
B_to_A_transfer is spectral cosine similarity of the A-channel state
to the initial B spectrum; it is NOT a path-exchange norm.
```

**[数値観測]** 真のB→A経路単独ノルム `path_b_to_a_norm_raw` は、全ケース・全衝突で約0.3028221209のまま変わらない。一方、`B_to_A_transfer` は次のように振る舞った。

| ケース | `B_to_A_transfer` 範囲 |
|---|---:|
| F1×F1 | 1.0000000000 |
| FK×FK | 1.0000000000 |
| BK×BK | 1.0000000000 |
| MIX×MIX | 1.0000000000 |
| FK×BK | 0.0006422305〜0.9999966927 |

FK×BKでは、

| 衝突 | 経路 `path_b_to_a_norm_raw` | `B_to_A_transfer` |
|---:|---:|---:|
| 1 | 0.3028221209 | 0.3983955931 |
| 2 | 0.3028221209 | 0.9834627923 |
| 8 | 0.3028221209 | 0.9999966927 |
| 16 | 0.3028221209 | 0.0103401424 |

となった。両者は同じ量ではない。

**[数学的帰結]** `B_to_A_transfer` の増減は、A出力の正規化された絶対倍音パワー分布が初期B分布へ近づいたかを示す。B→A経路の確率、交換率、経路ノルムとして使用してはならない。

## 9. 保存成果物

**[コード上の事実]**

```text
data/stage_B/current_behavior_baseline.csv
data/stage_B/state_parity_metrics.csv
data/stage_B/input_control_metrics.csv
data/stage_B/current_behavior_diagnostics.json
data/stage_B/path_arrays/F1_x_F1_collision_paths.npz
data/stage_B/path_arrays/FK_x_FK_collision_paths.npz
data/stage_B/path_arrays/BK_x_BK_collision_paths.npz
data/stage_B/path_arrays/FK_x_BK_collision_paths.npz
data/stage_B/path_arrays/MIX_x_MIX_collision_paths.npz
logs/stage_B_output_verification.json
manifests/stage_B_source_manifest_before.json
manifests/stage_B_source_manifest_after.json
manifests/stage_B_source_manifest_comparison.json
```

`current_behavior_baseline.csv` は160衝突行、`state_parity_metrics.csv` は初期状態を含む330状態行、`input_control_metrics.csv` は5ケース行を持つ。

永続化後のNPZとCSVを独立に再読込みした検証では、6条件がすべて成立した。

```text
persisted path norm error             <= 1.33e-15
persisted interference sum error      = 0
persisted raw-output reconstruction   <= 1.33e-15
```

## 10. Stage B停止

**[未導出]** Candidate 1〜3の数学設計、API設計、実装、比較は開始していない。現行核の再現結果だけを確定した。

Stage B 完了。既存System A／System B原本への変更は未実施。Candidate 1〜3には進まず停止する。
