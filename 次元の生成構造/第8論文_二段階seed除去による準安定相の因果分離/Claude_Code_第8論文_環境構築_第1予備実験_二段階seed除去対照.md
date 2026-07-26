# Claude Code 実行指示書
# 第8論文：二段階seed除去対照実験

## 0. 最重要命令

この作業では、以下に明記した範囲だけを実行すること。

- 独自の物理解釈を追加しない。
- 指示されていない追加解析を行わない。
- 指示されていない図を作らない。
- 指示されていないパラメータ掃引を行わない。
- seed振幅掃引、seed方向掃引、Lyapunov解析、自己相関解析、FFT解析を勝手に追加しない。
- 既存の第7論文フォルダ、プログラム、CSV、図、Markdownを変更しない。
- 実験結果を先に解釈して、コードや判定条件を変更しない。
- 期待した結果に合わせて閾値、時間範囲、縦軸、横軸、表示範囲を変更しない。
- 不足、実装上の曖昧さ、再現不能箇所があれば、推測で補わず報告して停止する。

本指示の目的は、二段階のseed投入について、seedを投入しない対照軌道を作り、現象のベースラインを確定することだけである。

---

# 1. 新規論文フォルダ

以下のルートフォルダを新規作成する。

```text
次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/
```

本指示書では、このフォルダを `PAPER8_ROOT` と呼ぶ。

第7論文のフォルダ名は、実際のリポジトリ内で確認し、以後 `PAPER7_ROOT` と呼ぶ。

第7論文の成果物は参照専用とし、書き換えない。

---

# 2. 第8論文の中心課題

第7論文の実験系には、役割の異なる二段階のseed投入がある。

## 2.1 第1seed：初期方向生成seed

初期二方向状態に対して、

\[
Z_0
=
\frac{v+\delta g}{\|v+\delta g\|},
\qquad
\delta=10^{-15}
\]

として投入される微小seedである。

過去の対照実験では、この初期seedを投入しない場合、幾何級数的発展が開始しなかった。

## 2.2 第2seed：準安定域への横摂動seed

第1seedによる幾何級数的発展が終了し、三方向準安定閉包が形成された後、活性部分空間の外側へ、

\[
\widetilde Z(t_1)
=
\frac{Z(t_1)+\epsilon\eta_\perp}
{\|Z(t_1)+\epsilon\eta_\perp\|}
\]

として投入される横摂動である。

第8論文の最初の課題は、二つのseedを別々に除去し、

\[
\boxed{
\begin{aligned}
S_0 &: \text{幾何級数的方向生成を起動するのか},\\
S_1 &: \text{準安定振動を生成するのか、既存運動を励起するだけか}
\end{aligned}
}
\]

を分離することである。

本指示では、そのための環境構築と第1予備実験だけを行う。

---

# 3. 作成するフォルダ構成

`PAPER8_ROOT` の直下に、以下を作成する。

```text
第8論文_二段階seed除去による準安定相の因果分離/
├── README.md
├── instructions/
│   └── 第1予備実験_二段階seed除去対照.md
├── code/
│   ├── run_preliminary_seed_ablation_v1.py
│   ├── make_preliminary_seed_ablation_figures_v1.py
│   └── audit_paper7_dependencies_v1.py
├── config/
│   ├── experiment_manifest.json
│   └── source_file_hashes.json
├── raw/
│   ├── N00005/
│   ├── N00040/
│   └── N00300/
├── summary/
├── figures/
├── diagnostics/
├── logs/
└── reports/
```

`instructions/第1予備実験_二段階seed除去対照.md` には、本指示書をそのまま保存する。

---

# 4. 第7論文から利用する実験環境

## 4.1 原則

第7論文で使用した以下を、そのまま再利用する。

- 親状態の生成法
- 初期生成子の構築法
- 状態依存生成子の更新則
- Cayley更新
- 正規化方法
- 零二乗閉鎖条件
- Nごとの既存親状態
- Nごとの既存乱数seed生成規則
- crossing定義
- 準安定開始時刻の定義
- 時刻刻み
- 記録時刻
- N=5、40、300の既存実験パラメータ
- 第7論文で用いた自然軌道の観測量

seedの有無以外を変更してはならない。

## 4.2 コピー方針

第7論文コードを直接編集しない。

必要なモジュールは、次のいずれかで利用する。

1. 第7論文コードをread-only importする。
2. importが困難な場合のみ、必要最小限のファイルを第8論文 `code/vendor_paper7/` へコピーする。

コピーした場合は、元ファイルとコピー後ファイルのSHA-256を `config/source_file_hashes.json` に保存し、コピー時点で同一であることを確認する。

第7論文の力学コードを変更して第8論文用に作り替えてはならない。seedのON/OFFは、第8論文側のラッパーから明示的に切り替える。

---

# 5. 依存関係監査

最初に `code/audit_paper7_dependencies_v1.py` を作成して実行する。

監査項目は以下に限定する。

1. 第7論文の実行コードの場所
2. N=5、40、300のパラメータ設定
3. 初期seedを加えている正確なコード位置
4. 準安定域で横摂動seedを加えている正確なコード位置
5. 初期状態、crossing、準安定開始時刻、最終時刻の定義
6. 既存の自然軌道データの場所
7. 既存の二段階seedあり実験データの場所
8. 再現に必要なPythonパッケージとバージョン
9. 乱数生成器とseed値
10. 第7論文コードおよび設定ファイルのSHA-256

監査結果を、

```text
reports/paper7_dependency_audit.md
```

へ出力する。

不明点または不足がある場合は、実験を開始せず、監査報告だけを出力して停止する。

---

# 6. 第1予備実験の全体設計

## 6.1 実験対象

必須対象は以下である。

```text
N = 5
N = 40
N = 300
```

## 6.2 比較する条件

第1予備実験では、次の三条件だけを扱う。

### 条件A：完全無seed対照

```text
initial_seed = OFF
metastable_seed = OFF
```

初期状態は厳密に、

\[
Z_0=v
\]

とする。

明示的な初期seedを加えない。

準安定域への横摂動も加えない。

過去の結果どおり幾何級数的発展が起きなければ、準安定開始判定は成立しない。その場合も、既存実験と同じ絶対最終stepまで自然発展を継続する。

### 条件B：初期seedあり・準安定seedなし

```text
initial_seed = ON
metastable_seed = OFF
```

初期状態は第7論文と同一に、

\[
Z_0
=
\frac{v+\delta g}{\|v+\delta g\|},
\qquad
\delta=10^{-15}
\]

とする。

幾何級数的発展、第三方向生成、準安定域への移行までは、第7論文と同じ自然軌道を生成する。

準安定開始後も、

\[
\widetilde Z(t_1)=Z(t_1)
\]

として、横摂動seedを一切加えず、そのまま自然発展を継続する。

### 条件D：既存の二段階seedあり基準

```text
initial_seed = ON
metastable_seed = ON
```

これは第7論文で既に実施した条件である。

原則として再計算しない。既存CSVを比較基準として参照する。

既存データが不足している、時刻軸が一致しない、または必要列が存在しない場合のみ、同一コード・同一設定で再生成する。その場合は理由を報告書へ明記する。

## 6.3 実施しない条件

以下は本予備実験では実施しない。

```text
initial_seed = OFF
metastable_seed = ON
```

初期seedなしで準安定状態が形成されない場合、準安定seed投入時刻を自然に定義できないためである。

この条件を人工的な時刻指定で作ってはならない。

---

# 7. seedをOFFにする厳密な実装

## 7.1 初期seed OFF

初期seed OFFでは、以下を禁止する。

- `delta = 0` とした後もseedベクトルを状態へ加えること
- seedベクトル生成処理によって乱数列を消費し、その後の乱数状態を変えること
- 正規化誤差を意図的seedとして加えること
- 微小ノイズを代替seedとして加えること
- 機械epsilonを明示的に加えること

初期状態は、第7論文の親状態 `v` をそのまま用いる。

\[
Z_0=v
\]

`v` が既に規格化・零二乗閉鎖されていることを診断出力する。

## 7.2 準安定seed OFF

準安定seed OFFでは、準安定開始時刻に状態を分岐・変更しない。

禁止事項：

- `epsilon = 0` とした横摂動ルーチンを通すこと
- 横方向ベクトル生成によって乱数列を消費すること
- 摂動軌道を別途作ること
- Benettin再正規化を実行すること
- warm-start状態を複製・変更すること
- 活性部分空間外への射影を状態更新へ使用すること

条件Bは、最初から最後まで一つの自然軌道だけを進める。

---

# 8. 共通時間範囲

N=5、40、300で、条件A・B・Dを同じ絶対step範囲で比較する。

第7論文で採用した共通最終stepを監査で取得し、その値を、

```text
COMMON_FINAL_STEP
```

として `config/experiment_manifest.json` に保存する。

各条件でcrossingや準安定開始が起こらなくても、`COMMON_FINAL_STEP` まで走行する。

crossingを0へ平行移動した図だけを作ってはならない。主比較は絶対stepで行う。

---

# 9. 保存する時系列

各N、各条件について、最低限、以下を同一列構成で保存する。

```text
step
time
N
condition
initial_seed_enabled
metastable_seed_enabled
initial_seed_amplitude
metastable_seed_amplitude
parent_plane_occupation
f_outside_parent
q1
q2
q3
q4
rank_Q
dominant_plane_occupation
non_dominant_occupation
kernel_occupation
residual_occupation
norm_Z
dagger_norm_error
zero_square_real
zero_square_imag
zero_square_abs
projection_closure_error
crossing_detected
metastable_start_detected
```

出力先：

```text
raw/N00005/condition_A_no_seed.csv
raw/N00005/condition_B_initial_only.csv
raw/N00005/condition_D_existing_two_seed.csv

raw/N00040/condition_A_no_seed.csv
raw/N00040/condition_B_initial_only.csv
raw/N00040/condition_D_existing_two_seed.csv

raw/N00300/condition_A_no_seed.csv
raw/N00300/condition_B_initial_only.csv
raw/N00300/condition_D_existing_two_seed.csv
```

既存条件DのCSV列が不足する場合は、存在する列だけを勝手に補間しない。再生成するか、不足として報告する。

---

# 10. 本予備実験で作成する図

図は次の4種類だけを作成する。

## 図1：親平面外占有

各Nについて、条件A・B・Dの

\[
f(t)=1-E_{P_1}(t)
\]

を同一図へ重ねる。

```text
figures/fig01_f_compare_N00005.png
figures/fig01_f_compare_N00040.png
figures/fig01_f_compare_N00300.png
```

## 図2：q3、q4

各Nについて、条件A・B・Dの `q3`, `q4` を比較する。

```text
figures/fig02_q3q4_compare_N00005.png
figures/fig02_q3q4_compare_N00040.png
figures/fig02_q3q4_compare_N00300.png
```

## 図3：rank Q

各Nについて、条件A・B・Dの `rank_Q` を比較する。

```text
figures/fig03_rankQ_compare_N00005.png
figures/fig03_rankQ_compare_N00040.png
figures/fig03_rankQ_compare_N00300.png
```

## 図4：準安定域の自然軌道と追加seed軌道

条件Bと条件Dについて、第7論文で用いた準安定振動の代表観測量を同一図へ重ねる。

使用する観測量は、第7論文の図で準安定振動を示した既存量をそのまま使う。新しい合成指標を作らない。

```text
figures/fig04_metastable_B_vs_D_N00005.png
figures/fig04_metastable_B_vs_D_N00040.png
figures/fig04_metastable_B_vs_D_N00300.png
```

図4は、準安定開始後の区間を示してよいが、絶対stepを併記する。

---

# 11. 作成する集計表

以下を作成する。

```text
summary/preliminary_seed_ablation_summary.csv
```

列：

```text
N
condition
initial_seed_enabled
metastable_seed_enabled
crossing_detected
crossing_step
metastable_start_detected
metastable_start_step
max_f
final_f
max_q3
max_q4
final_q3
final_q4
max_rank_Q
final_rank_Q
mean_metastable_observable
std_metastable_observable
max_norm_error
max_zero_square_error
max_projection_closure_error
```

`mean_metastable_observable` と `std_metastable_observable` は、第7論文で準安定振動の代表として既に使用した観測量に対してのみ計算する。

準安定開始が検出されない条件では、該当列を `NaN` とし、別の区間を勝手に準安定域と定義しない。

---

# 12. 判定は観測事実だけを記載する

最終報告では、以下の問いに数値で答える。

## 12.1 初期seed除去

条件Aについて、

1. 幾何級数的発展が起きたか
2. crossingが検出されたか
3. `rank_Q` が2から増えたか
4. `q3`, `q4` が数値床を超えたか
5. 準安定開始が検出されたか

を報告する。

## 12.2 準安定seed除去

条件Bと条件Dを比較し、

1. 条件Bでも準安定振動が継続したか
2. 条件BとDで振動振幅が異なるか
3. 条件BとDで振動中心値が異なるか
4. 条件BとDで方向数が異なるか
5. 条件Dだけに現れる応答があるか

を報告する。

ただし、本予備実験では以下を実施しない。

- FFTによる周波数同定
- 自己相関時間の算出
- 位相同期解析
- Lyapunov指数の再計算
- seed振幅依存則の推定
- seed方向依存性の推定
- スケール不変性の判定
- 初期潜伏相と後期準安定相の自己相似性判定

これらは次段階の実験である。

---

# 13. 数値診断

全条件で、最低限、以下を確認する。

\[
\|Z\|^2
\]

\[
Z^T Z
\]

\[
\|K+K^T\|
\]

\[
\left|E_{\mathrm{dom}}+E_{\mathrm{nondom}}+E_{\mathrm{ker}}+E_{\mathrm{res}}-1\right|
\]

N=40では、第7論文で検証済みの厳密法と低ランク法の一致条件を壊していないことを確認する。ただし、本予備実験のために新しい精度比較実験を追加してはならない。

診断出力：

```text
diagnostics/N00005_condition_A.json
diagnostics/N00005_condition_B.json
diagnostics/N00005_condition_D.json
...
```

---

# 14. 実行順序

以下の順序を厳守する。

## Phase 1：第7論文環境監査

`audit_paper7_dependencies_v1.py` を実行し、監査報告を作成する。

不足があれば停止する。

## Phase 2：N=5 条件A

初期seedなし、準安定seedなしで実行する。

結果と数値診断を確認する。

## Phase 3：N=5 条件B

初期seedあり、準安定seedなしで実行する。

第7論文の自然軌道と、準安定開始まで一致することを確認する。

## Phase 4：N=5 条件Dデータ接続

既存データを読み込み、A・B・Dの比較CSVと図を作る。

この段階で実装上の問題があれば、N=40へ進まない。

## Phase 5：N=40

N=5で確定した同一コードを用い、条件A・Bを実行する。

条件Dは既存データを使用する。

## Phase 6：N=300

N=5、40と同じ条件・同じ判定規則で、条件A・Bを実行する。

条件Dは既存データを使用する。

## Phase 7：比較図・集計表・報告書

全NのCSVが揃ってからのみ作成する。

---

# 15. 完了条件

以下をすべて満たした場合のみ完了とする。

1. 新規第8論文フォルダが指定名で作成されている
2. 第7論文フォルダは変更されていない
3. 第7論文依存関係監査報告がある
4. 参照・コピーしたコードのSHA-256が保存されている
5. N=5、40、300の条件Aがある
6. N=5、40、300の条件Bがある
7. N=5、40、300の条件D既存データが接続されている
8. A・B・Dが同一絶対step範囲で比較されている
9. 全条件で同一列構成のCSVがある
10. 指定した図1〜4だけが作成されている
11. 集計CSVがある
12. 数値診断JSONがある
13. 初期seed OFF時にseedベクトル生成で乱数を消費していない
14. 準安定seed OFF時に摂動軌道とBenettin処理を作っていない
15. 指示外の解析を追加していない
16. 観測結果を解釈で改変していない
17. 最終報告書がある

---

# 16. 最終報告書

以下へ保存する。

```text
reports/preliminary_seed_ablation_report.md
```

構成は次に限定する。

```text
1. 実行環境
2. 第7論文から再利用したコードと設定
3. 条件A・B・Dの定義
4. N=5の観測結果
5. N=40の観測結果
6. N=300の観測結果
7. 二つのseed除去による差分
8. 数値診断
9. 未実施項目
10. 生成ファイル一覧
```

結論文は、観測された内容だけを記載する。

許可される書き方の例：

> 初期seedを除去した条件Aでは、観測終端までcrossingは検出されず、rank Qは2を維持した。

> 初期seedを保持し、準安定域への追加seedを除去した条件Bでも、準安定観測量の時間変動は継続した。

禁止される書き方の例：

> seedは宇宙創生を引き起こす。

> 準安定振動は量子揺らぎである。

> この結果はスケール不変宇宙を証明する。

本予備実験の目的は、二段階seedの有無による事実上の差分を確定することであり、物理的名称付けや一般化ではない。

---

# 17. 実験マニフェスト

`config/experiment_manifest.json` に、最低限以下を保存する。

```json
{
  "paper": 8,
  "experiment": "preliminary_two_stage_seed_ablation_v1",
  "paper8_root": "次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離",
  "N_values": [5, 40, 300],
  "conditions": {
    "A": {
      "initial_seed": false,
      "metastable_seed": false
    },
    "B": {
      "initial_seed": true,
      "metastable_seed": false,
      "initial_delta": 1e-15
    },
    "D": {
      "initial_seed": true,
      "metastable_seed": true,
      "source": "existing_paper7_data"
    }
  },
  "forbidden_extra_analysis": [
    "seed_amplitude_sweep",
    "seed_direction_sweep",
    "fft",
    "autocorrelation",
    "lyapunov_reanalysis",
    "self_similarity_analysis",
    "physical_interpretation"
  ]
}
```

`COMMON_FINAL_STEP`、実際のseed値、コードパス、SHA-256、Python環境は監査後に追記する。

---

# 18. 実行開始時の確認表示

Claude Codeは実行前に、以下だけを表示する。

```text
PAPER7_ROOT = ...
PAPER8_ROOT = 次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離
N = 5, 40, 300
Condition A = initial seed OFF / metastable seed OFF
Condition B = initial seed ON / metastable seed OFF
Condition D = existing paper 7 two-seed baseline
No extra analysis will be performed.
```

その後、Phase 1の監査から開始する。
