# Claude Code 実行指示書
# 第8論文・第2予備実験
# 初期フラット期間における N–分解能掃引（完全無seed）

## 0. 最重要命令

この作業では、本指示書に明記した実験だけを実行すること。

### 禁止事項

- 実験結果の物理解釈を行わない。
- 結果の意味、理論的意義、期待との一致・不一致を記述しない。
- 「成功」「失敗」「支持」「否定」「インフレーション」「零点揺らぎ」「真空揺らぎ」等の評価語・解釈語を使用しない。
- 期待される結論を推測しない。
- 結果に合わせて実験条件、分解能、停止条件、回帰区間、閾値、図の軸範囲を変更しない。
- 指示されていない追加実験、追加解析、FFT、自己相関、Lyapunov解析、パワースペクトル解析を行わない。
- 最も見栄えのよい区間を事後選択しない。
- 外部seed、乱数摂動、横摂動、初期方向摂動を加えない。
- 第7論文および第8論文第1予備実験の既存コード・データ・図・Markdownを変更しない。
- 不明点を推測で補わない。不足または再現不能箇所がある場合は、監査報告だけを出力して停止する。

### Claude Codeの役割

Claude Codeが行うのは、以下だけである。

1. 既存コードと既存条件の監査
2. 指定された完全無seed条件の実装
3. 指定された N と分解能の組合せの実行
4. 指定された観測量の逐次記録
5. 事前固定された方法による数値回帰値の計算
6. CSV、JSON、PNG、実行ログ、監査報告の出力

結果の分析・解釈・文章化は行わない。

---

# 1. 実験目的

初期状態に明示的なseedを一切加えず、初期フラット期間だけを観測する。

変更する制御変数は、次の二つだけとする。

```text
N
有限分解能 Δ
```

N と Δ の組合せごとに、初期平面外振幅の時間系列を記録する。

本実験では、準安定相の長時間観測を目的としない。初期フラット期間から初期成長期間までを高密度で保存する。

---

# 2. 新規フォルダ

第8論文ルートを `PAPER8_ROOT` とする。

以下を新規作成する。

```text
PAPER8_ROOT/
└── preliminary_02_resolution_scan/
    ├── README.md
    ├── instructions/
    │   └── 第2予備実験_初期期間_N分解能掃引_完全無seed.md
    ├── code/
    │   ├── audit_resolution_scan_dependencies_v1.py
    │   ├── run_initial_resolution_scan_v1.py
    │   ├── compute_prefixed_regressions_v1.py
    │   └── make_initial_resolution_scan_figures_v1.py
    ├── config/
    │   ├── experiment_manifest.json
    │   └── source_file_hashes.json
    ├── raw/
    ├── summary/
    ├── figures/
    ├── diagnostics/
    ├── logs/
    └── reports/
```

本指示書を、次へそのまま保存する。

```text
instructions/第2予備実験_初期期間_N分解能掃引_完全無seed.md
```

---

# 3. 既存環境の監査

最初に、次を監査する。

1. 第8論文第1予備実験で使用した実行コード
2. N=5、40、300の親状態生成法
3. 初期状態 `Z0 = v` の構築位置
4. 初期seedをOFFにした正確なコード位置
5. 準安定seedをOFFにした正確なコード位置
6. 1 stepの更新順序
7. Cayley更新、閉鎖条件への射影、正規化の順序
8. `f_outside`、`q3`、`q4`、`rank_Q` の算出法
9. crossing判定式
10. 時刻刻みと既存の記録間隔
11. 使用する数値型、線形代数ライブラリ、丸め方式
12. 参照する全コード・設定ファイルのSHA-256

監査結果を以下へ保存する。

```text
reports/resolution_scan_dependency_audit.md
config/source_file_hashes.json
```

既存コードを直接編集してはならない。

---

# 4. 固定する力学条件

第1予備実験の条件Aを基準とする。

```text
initial_seed    = OFF
metastable_seed = OFF
random_kick     = OFF
external_noise  = OFF
```

初期状態は厳密に、

\[
Z_0=v
\]

とする。

以下は第1予備実験から変更しない。

- 親状態
- 初期生成子
- 状態依存生成子
- Cayley更新
- 時刻刻み
- 零二乗閉鎖条件
- 正規化方法
- crossing定義
- 観測量の定義
- Nごとの既存パラメータ

本実験で新たに変更してよいのは、有限分解能写像の設定だけである。

---

# 5. 分解能写像の定義

## 5.1 分解能パラメータ

絶対分解能を、

\[
\Delta>0
\]

とする。

複素状態の実部・虚部に対して、決定論的な最近接偶数丸めを適用する。

\[
Q_\Delta(x)
=
\Delta\,\operatorname{round}_{\mathrm{half\ to\ even}}
\!\left(\frac{x}{\Delta}\right)
\]

複素数について、

\[
Q_\Delta(a+ib)=Q_\Delta(a)+iQ_\Delta(b)
\]

とする。

乱数丸めは使用しない。

## 5.2 1 stepの固定順序

各stepは、以下の順序で一度だけ実行する。

```text
1. 既存の状態依存生成子を計算
2. 既存のCayley更新を実行
3. 既存手順で零二乗閉鎖条件へ射影
4. 既存手順でノルムを正規化
5. Q_Δ を状態の実部・虚部へ一度適用
6. 既存手順で零二乗閉鎖条件へ再射影
7. 既存手順でノルムを再正規化
8. 観測量と診断量を保存
9. この状態を次stepの入力とする
```

この順序を途中で変更しない。

## 5.3 分解能写像による残差の記録

各stepで、少なくとも次を保存する。

```text
quantization_l2
closure_residual_before_quantization
closure_residual_after_quantization
closure_residual_after_reprojection
norm_before_quantization
norm_after_quantization
norm_after_reprojection
```

ここで、

\[
\text{quantization\_l2}
=
\|Q_\Delta(Z)-Z\|_2
\]

とする。

---

# 6. N と分解能の組合せ

## 6.1 N

必須のNは以下とする。

```text
N = 5
N = 40
N = 300
```

関係成分数を、

\[
M(N)=\frac{N(N-1)}{2}
\]

として記録する。

## 6.2 N依存則を仮定しないためのパラメータ化

基準Nを、

\[
N_{\mathrm{ref}}=40
\]

とする。

基準絶対分解能を \(\Delta_{\mathrm{ref}}\)、スケーリング指数を \(p\) とし、各Nの絶対分解能を、

\[
\boxed{
\Delta(N;\Delta_{\mathrm{ref}},p)
=
\Delta_{\mathrm{ref}}
\left(
\frac{M(N)}{M(N_{\mathrm{ref}})}
\right)^{-p/2}
}
\]

で機械的に生成する。

この式は実験条件を配置するためだけに使用し、正しい関係であるとは記述しない。

スケーリング指数は、以下を同格に扱う。

```text
p = 0.0
p = 0.5
p = 1.0
p = 1.5
p = 2.0
```

意味の解釈は行わない。

## 6.3 基準絶対分解能

第1段階では、以下を固定する。

```text
Δ_ref = 1e-4
Δ_ref = 1e-6
Δ_ref = 1e-8
Δ_ref = 1e-10
Δ_ref = 1e-12
```

したがって必須実験数は、

```text
3 N values × 5 p values × 5 Δ_ref values = 75 runs
```

である。

各runの実際の \(\Delta\) は、実行前にmanifestへ保存する。

## 6.4 無分解能写像基準

各Nについて、分解能写像を適用しない既存条件Aを1本だけ再実行する。

```text
resolution_operator = OFF
```

合計3本とする。

この基準runにも外部seedを一切加えない。

---

# 7. 実行順序

結果依存の順序変更を避けるため、実行順序を事前に固定する。

```text
1. Nの昇順: 5, 40, 300
2. pの昇順: 0.0, 0.5, 1.0, 1.5, 2.0
3. Δ_refの降順: 1e-4, 1e-6, 1e-8, 1e-10, 1e-12
4. 最後に各Nのresolution_operator=OFF基準run
```

並列実行する場合も、run IDはこの順序で確定する。

---

# 8. 観測期間と停止条件

## 8.1 初期期間の上限

Nごとの最大stepは、第1予備実験・条件Aの既知crossing時刻を参照し、以下に固定する。

```text
N = 5   : max_step = 2500
N = 40  : max_step = 4500
N = 300 : max_step = 10000
```

## 8.2 早期停止

以下のいずれかを満たした最初のstepで停止する。

```text
A. f_outside >= 1e-2
B. 数値例外、NaN、Infが発生
C. closure_residual_after_reprojection > 1e-8
D. norm_after_reprojection が [1-1e-10, 1+1e-10] を外れる
E. max_stepに到達
```

停止理由を機械可読な文字列で保存する。

```text
f_outside_limit
numerical_exception
closure_residual_limit
norm_limit
max_step
```

停止条件をrun途中で変更しない。

---

# 9. 保存間隔

初期期間を観測するため、以下を固定する。

```text
step 0〜1000    : 毎step保存
step 1001以降   : 5 stepごとに保存
停止step        : 必ず保存
```

全状態ベクトル `Z` は、以下だけ保存する。

```text
step = 0, 1, 2, 5, 10, 20, 50, 100,
200, 500, 1000,
および停止step
```

ただし、既存状態サイズにより保存不能な場合は、監査報告に必要容量を記載して停止し、勝手に間引き規則を変更しない。

---

# 10. 必須観測量

各保存stepで、少なくとも以下をCSVへ出力する。

```text
run_id
N
M
p
Delta_ref
Delta_actual
resolution_operator
step
time
f_outside
a_outside
log_a_outside
q3
q4
rank_Q
E_dom
closure_residual_before_quantization
closure_residual_after_quantization
closure_residual_after_reprojection
norm_before_quantization
norm_after_quantization
norm_after_reprojection
quantization_l2
nonzero_real_count
nonzero_imag_count
min_nonzero_abs_component
max_abs_component
stop_reason
```

ここで、

\[
a_{\mathrm{outside}}=\sqrt{f_{\mathrm{outside}}}
\]

とする。

`a_outside = 0` の場合、`log_a_outside` は空欄またはNaNとし、任意の微小値を加えてはならない。

---

# 11. 事前固定回帰

Claude Codeは回帰値を計算してよいが、解釈・選別・評価をしてはならない。

## 11.1 回帰対象

各runについて、以下の固定された振幅帯ごとに、

\[
\log a_{\mathrm{outside}}(t)
=
\alpha+\gamma t
\]

の通常最小二乗回帰を行う。

```text
Band B1: 1e-14 <= a_outside < 1e-12
Band B2: 1e-12 <= a_outside < 1e-10
Band B3: 1e-10 <= a_outside < 1e-8
Band B4: 1e-8  <= a_outside < 1e-6
Band B5: 1e-6  <= a_outside < 1e-4
Band B6: 1e-4  <= a_outside < 1e-2
```

各bandで有効点数が20未満の場合、回帰を行わず `insufficient_points` と記録する。

## 11.2 出力する数値

各bandについて、以下だけを保存する。

```text
n_points
first_step
last_step
alpha
gamma
exp_gamma
r_squared
rmse
min_a
max_a
status
```

`exp_gamma = exp(gamma)` とする。

回帰区間を追加・結合・移動しない。

## 11.3 隣接成長率

保存間隔を考慮して、各隣接保存点について、

\[
\gamma_{\mathrm{local}}
=
\frac{
\log a(t_2)-\log a(t_1)
}{t_2-t_1}
\]

を計算する。

両端の `a_outside` が正の場合だけ保存する。

---

# 12. 出力CSV

各runについて、以下を作成する。

```text
raw/<run_id>/timeseries.csv
raw/<run_id>/local_growth.csv
raw/<run_id>/regression_by_fixed_band.csv
raw/<run_id>/run_config.json
raw/<run_id>/run_diagnostics.json
raw/<run_id>/stdout.log
raw/<run_id>/stderr.log
```

全runの結合表を作成する。

```text
summary/all_runs_manifest.csv
summary/all_runs_final_values.csv
summary/all_fixed_band_regressions.csv
summary/all_stop_reasons.csv
summary/all_diagnostics.csv
```

CSVの行を結果に応じて削除してはならない。

---

# 13. 必須図

図は指定されたものだけを機械的に作成する。

## 13.1 各runの初期時間系列

各runについて1枚作成する。

```text
横軸: step
縦軸: a_outside（対数軸）
```

- 全保存点を描く。
- 縦軸範囲を手動変更しない。
- 0値は描画せず、欠測として扱う。
- 回帰線は描かない。

## 13.2 N別・p別比較

同じN、同じpについて、5つの `Delta_ref` の `a_outside(step)` を1枚に重ねる。

合計、

```text
3 N × 5 p = 15 figures
```

とする。

## 13.3 p別・N比較

同じp、同じ `Delta_ref` について、N=5、40、300を1枚に重ねる。

合計、

```text
5 p × 5 Delta_ref = 25 figures
```

とする。

## 13.4 固定軸版

比較図には、以下の固定軸版も作る。

```text
横軸: step 0〜10000
縦軸: a_outside 1e-16〜1e-2（対数軸）
```

データが範囲外でも軸を変更しない。

図中に「成功」「指数成長」「相転移」等の文字を入れない。

---

# 14. 再現性確認

各runは同一条件で2回実行する。

乱数は使用しない。

二つの実行について、以下を比較する。

```text
停止step
停止理由
全保存stepのf_outside
全保存stepのa_outside
全保存stepのquantization_l2
最終状態ZのSHA-256（保存可能な場合）
```

完全一致しない場合、差分を保存する。

```text
diagnostics/reproducibility_diff_<run_id>.csv
```

不一致を解釈しない。

したがって実行本数は、基準runを含め、

```text
(75 + 3) × 2 = 156 executions
```

である。

---

# 15. 検算用メタデータ

各runの `run_config.json` に以下を必ず保存する。

```text
run_id
execution_index
N
M
N_ref
M_ref
p
Delta_ref
Delta_actual
initial_seed=false
metastable_seed=false
random_kick=false
external_noise=false
resolution_operator
rounding_mode=half_to_even
max_step
save_schedule
stop_conditions
python_version
numpy_version
scipy_version
platform
source_hashes
start_timestamp
end_timestamp
```

浮動小数点環境について、取得可能なら以下も保存する。

```text
float type
mantissa bits
machine epsilon
smallest normal
smallest subnormal
BLAS/LAPACK information
```

---

# 16. READMEの記載制限

READMEには以下だけを書く。

1. 実験名
2. 実行方法
3. フォルダ構成
4. 設定ファイルの場所
5. 出力ファイル一覧
6. 再実行方法
7. 依存パッケージ

結果の要約、意味、解釈、結論は書かない。

---

# 17. 最終実行報告

以下を作成する。

```text
reports/execution_report.md
```

記載項目は、以下に限定する。

```text
- 実行開始・終了日時
- 実行済みrun数
- 未実行run数
- 正常終了run数
- 各停止理由の件数
- 例外発生run一覧
- 欠損ファイル一覧
- 再現性一致run数
- 再現性不一致run一覧
- 生成ファイル一覧
- SHA-256一覧
```

数値結果の意味を文章で説明しない。

---

# 18. 完了条件

以下をすべて満たした場合だけ完了とする。

- 依存関係監査が完了している。
- 既存コードが変更されていない。
- 全runの設定が事前にmanifestへ固定されている。
- 全156 executionが実行済みである。
- 各runに必要なCSV、JSON、ログが存在する。
- 固定band回帰が機械的に計算されている。
- 指定図がすべて生成されている。
- 再現性比較が完了している。
- 実行報告に解釈が含まれていない。

未完了項目がある場合は、未完了のまま明示して終了する。条件を緩和して完了扱いにしてはならない。
