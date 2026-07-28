# Codex向け作業指示
## Stage F：反転Candidate 1を既存System A実験コピーへ限定統合し、従来実験を同条件で再実行する

---

# 0. このStageの目的

Stage A〜Eにより、以下は確認済みである。

1. 現行System A散乱核はパリティ盲である。
2. Stage Dの二層表現により、
   - 物理散乱は \(512\times16\) 全状態
   - パリティ読出しは復調した一次元カーネル
   として矛盾なく分離できる。
3. 反転Candidate 1
   \[
   \theta_{\mathrm{eff}}
   =
   \theta_0
   -
   \kappa\rho(\theta_0)
   \frac{c_A+c_B}{2}
   \]
   は、System A型の反復全状態散乱へ数値的に安定して実装できる。
4. 純状態では
   \[
   F\times F:\ R_{\mathrm{eff}}>R_0,
   \qquad
   B\times B:\ R_{\mathrm{eff}}<R_0,
   \qquad
   F\times B:\ R_{\mathrm{eff}}=R_0
   \]
   が成立する。
5. ユニタリ性、経路和、ノルム保存、復調再構成は倍精度誤差内で成立する。

したがって本Stageでは、追加候補の探索や単体波の追加観察は行わない。

**既存System Aで過去に実行していた本来の局在性交換・再帰散乱・R掃引実験を、反転Candidate 1を組み込んだ独立実験コピーで再実行し、C0との差を直接比較する。**

このStageは、仮説の正式採用ではなく、既存研究系列へ接続したときに何が変わるかを確認する限定統合実験である。

---

# 1. 作業範囲

作業先を次に固定する。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/
codex_systemA_scattering_audit_v1/
stage_F_original_systemA_integration/
```

この配下にのみ、新規コード、データ、図、報告書、ログを作成する。

参照してよい主対象は、Stage Eで保存した原本同一SHA-256のSystem Aスナップショットである。

例：

```text
stage_E_reversed_candidate1/
source_copy/
run_system_A_localization_exchange_R_sweep_preliminary_v1_ORIGINAL_SNAPSHOT.py
```

実際の既存配置を確認し、正確な参照元パスを報告書に記録すること。

## 禁止事項

- 既存System A原本の変更
- 既存System B原本の変更
- Stage E成果物の上書き
- 既存CSV、図、報告書の上書き
- 配下外への出力
- Git add、commit、push、stage
- Candidate 2またはCandidate 3の追加
- 新しい散乱候補の探索
- \(\kappa\) の結果合わせ
- N体系への組込み
- 論文本文の自動修正
- Stage F終了後の自動継続

---

# 2. 統合対象

既存System A実験コピーへ、次の二つの散乱モードを切替可能にして統合する。

```python
scattering_mode = "C0"
scattering_mode = "reversed_C1"
```

## 2.1 C0

既存散乱核をそのまま再現する。

\[
\theta_{\mathrm{eff}}=\theta_0
\]

\[
t_0=e^{i\theta_0}\cos\theta_0,
\qquad
r_0=-ie^{i\theta_0}\sin\theta_0
\]

## 2.2 reversed_C1

毎衝突直前の入力状態からパリティ指標を再計算する。

\[
\bar c_n
=
\frac{c_{A,n}+c_{B,n}}{2}
\]

\[
\theta_{\mathrm{eff},n}
=
\theta_0
-
\kappa
\rho(\theta_0)
\bar c_n
\]

\[
\rho(\theta_0)
=
\frac{2}{\pi}
\theta_0
\left(
\frac{\pi}{2}-\theta_0
\right)
\]

\[
t_{\mathrm{eff},n}
=
e^{i\theta_{\mathrm{eff},n}}
\cos\theta_{\mathrm{eff},n}
\]

\[
r_{\mathrm{eff},n}
=
-i e^{i\theta_{\mathrm{eff},n}}
\sin\theta_{\mathrm{eff},n}
\]

\[
R_{\mathrm{eff},n}
=
|r_{\mathrm{eff},n}|^2
\]

\[
T_{\mathrm{eff},n}
=
|t_{\mathrm{eff},n}|^2
\]

条件分岐で波を「ボゾン」「フェルミオン」と分類して係数を選んではならない。

同一式へ復調パリティ指標を入力すること。

---

# 3. 二層表現

Stage Eで検証済みの二層表現をそのまま使用する。

## 3.1 全状態物理散乱層

既存System Aと同じ、

- \(512\times16\) 全状態
- 搬送波
- A/Bの直交 \(\eta\) モード
- 複素配列
- 既存のチャネル構造

を使用する。

散乱は、

\[
a_{n+1}^{\mathrm{raw}}
=
r_{\mathrm{eff},n}a_n
+
t_{\mathrm{eff},n}b_n
\]

\[
b_{n+1}^{\mathrm{raw}}
=
t_{\mathrm{eff},n}a_n
+
r_{\mathrm{eff},n}b_n
\]

で行う。

## 3.2 パリティ読出し層

各衝突直前のA/B入力状態を、A由来・B由来の \(\eta\) モードごとに射影し、対応する搬送波を除去して一次元カーネルへ復調する。

\[
(P\psi)(u)=\psi(u+\pi)
\]

\[
c_\pi[\psi]
=
\frac{
\operatorname{Re}\langle\psi,P\psi\rangle
}{
\|\psi\|^2
}
\]

を使う。

A/Bチャネル全体を単一搬送波で復調してはならない。

由来別成分を保持し、対応する復調を行った後にパリティ量を集計すること。

---

# 4. raw更新と既存正規化の扱い

Stage Eではraw出力をそのまま次状態に使用し、数値的に安定した。

一方、既存System Aには衝突後の個別正規化処理が存在する可能性がある。

本Stageでは、以下の二系列を明確に分離する。

## 系列F0：既存挙動再現系列

既存System Aの更新手順を可能な限り保持する。

既存コードにチャネル別正規化がある場合は、そのまま使用する。

目的は、過去のSystem A結果をC0で再現することである。

## 系列F1：raw物理更新系列

正規化前の

\[
a_{n+1}^{\mathrm{raw}},
\qquad
b_{n+1}^{\mathrm{raw}}
\]

を次状態として使用する。

目的は、状態依存散乱が作る強度差を個別正規化で消さずに読むことである。

## 重要

C0/reversed_C1と、既存正規化/raw更新を混同しない。

最低限、次の4条件を比較可能にする。

```text
C0 + existing_normalization
reversed_C1 + existing_normalization
C0 + raw_update
reversed_C1 + raw_update
```

ただし、既存System Aの基準条件で既存正規化とraw更新が数値的に完全一致する場合も、両者を別系列として記録し、一致を報告する。

結果を見て片方を削除してはならない。

---

# 5. 最初に行う再現ゲート

新しい本実験を開始する前に、C0で既存System A結果を再現する。

最低限、既存報告にある次の基準ケースを特定する。

- \(N_A=1\)
- \(N_B=63\)
- 基準反射率 \(R=0.55\)
- 既存の反復回数
- 既存の初期条件
- 既存の \(L\)
- 既存の \(N_{\mathrm{eff}}\)
- 既存のスペクトル類似度
- 既存のA/B交換挙動

既存記録として、少なくとも以下が報告されている。

- 1回散乱後：
  - minus側 \(N_{\mathrm{eff}}\approx14.95\)
  - plus側 \(N_{\mathrm{eff}}\approx18.05\)
- 反復中：
  - 約16回で両者 \(N_{\mathrm{eff}}\approx17\)
  - 約32回で入替傾向
- 収束ではなく周期的交換

実際の既存コード・既存データを確認し、丸めた記憶値ではなく、正本の数値を再現基準として使用すること。

## 再現判定

C0再現について、以下を比較する。

```text
L_A
L_B
N_eff_A
N_eff_B
spectral_similarity_to_initial_A
spectral_similarity_to_initial_B
channel_norm_A
channel_norm_B
```

各衝突系列の最大絶対誤差と最大相対誤差を報告する。

再現できない場合は、新散乱実験へ進まず停止する。

---

# 6. 実験群

再現ゲート通過後、以下を実行する。

---

## 6.1 実験F-A：既存代表条件での反復比較

固定条件：

- \(N_A=1\)
- \(N_B=63\)
- \(R_0=0.55\)
- 既存と同じ初期波形
- 既存と同じ搬送波・\(\eta\)モード
- 既存と同じ反復回数
- C0とreversed_C1
- existing_normalizationとraw_update

\(\kappa\) は、

\[
\kappa\in\{0.01,0.1,1\}
\]

とする。

ただし、\(\kappa=1\) で

\[
0\le\theta_{\mathrm{eff}}\le\frac{\pi}{2}
\]

を破る場合は、クリップして結果を続行せず、該当条件を明示して停止または除外する。

自動クリップで物理結果を作ってはならない。

### 取得値

各衝突で最低限、次を保存する。

```text
collision_index
scattering_mode
normalization_mode
kappa
R0
theta0
rho
c_A
c_B
c_mean
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
path_sum_residual_A
path_sum_residual_B
demodulation_reconstruction_residual
```

---

## 6.2 実験F-B：既存R掃引の再実行

既存System Aで用いたR掃引範囲と刻みを正本コードから読み取る。

勝手に新しい範囲へ変更しない。

C0で既存R掃引を再現した後、同じR点をreversed_C1で実行する。

最低限、

\[
R\in\{0,\ 0.55,\ 0.70,\ 1\}
\]

を含むことを確認する。

既存点に含まれていなければ、正本掃引を保持したまま追加点として明示する。

各Rについて、

- 最小 \(|L_A-L_B|\)
- その衝突回数
- 最小 \(|N_{\mathrm{eff},A}-N_{\mathrm{eff},B}|\)
- その衝突回数
- 最大・最小 \(R_{\mathrm{eff}}\)
- \(c_A,c_B,\bar c\) の範囲
- 周期性指標
- 交換回数
- 準安定性
- 発散・収束・周期交換の分類

を保存する。

過去の代表観測として、

\[
R=0.70
\]

付近で

\[
|L_A-L_B|
\]

および

\[
|N_{\mathrm{eff},A}-N_{\mathrm{eff},B}|
\]

が小さくなる記録があるため、C0再現値とreversed_C1との差を明示する。

---

## 6.3 実験F-C：31系列・帰還周期の限定確認

既存研究で観測済みの31系列、

\[
31,62,93,124,155,186,217,248,279
\]

および、247/248近傍の帰還・誤差最小構造を、既存正本コード・データから確認する。

このStageで新しい周期を探索してはならない。

既存の31系列評価方法をそのまま使用して、

- C0
- reversed_C1

を同じ条件で比較する。

最低限、以下を出す。

```text
iteration
C0_return_error
reversed_C1_return_error
C0_exchange_measure
reversed_C1_exchange_measure
C0_L_difference
reversed_C1_L_difference
C0_N_eff_difference
reversed_C1_N_eff_difference
c_mean
R_eff
```

31系列が既存System Aの別条件・別コードに依存し、本統合コピーへ直接適用できない場合は、無理に再現しない。

その場合は、

- 何が不足しているか
- どのコード・データが必要か
- なぜ今回の範囲では接続不能か

を報告し、F-Cのみ未実施としてよい。

配下外のコードを修正して接続してはならない。

---

# 7. 主要な比較問い

本Stageでは、次の問いに答える。

## 問い1

反転Candidate 1は、既存の局在性交換周期を変えるか。

## 問い2

既存の \(L\) 交換と \(N_{\mathrm{eff}}\) 交換は、

\[
R_0
\]

ではなく、毎衝突変動する

\[
R_{\mathrm{eff},n}
\]

によって説明しやすくなるか。

## 問い3

奇数倍音主体の既存 \(N_A=1,N_B=63\) 条件では、

\[
c_A,\ c_B,\ \bar c
\]

はどのように動くか。

両者とも奇数型から始まる場合、反転Candidate 1が単に一定の高反射率へ固定されるのか、それとも局在性交換に伴いパリティ純度が変化して \(R_{\mathrm{eff}}\) が時間変動するのか。

## 問い4

既存の準安定振動は、

- C0でも同じ
- reversed_C1で周期だけ変わる
- 振幅が変わる
- 固定点へ近づく
- 発散する
- 新しい準安定相へ移る

のどれか。

## 問い5

既存正規化を残した場合とraw更新の場合で、状態依存散乱差が消えるか、残るか。

---

# 8. 事前登録された判定分類

結果を見て分類名を変更しない。

各条件を次で分類する。

```text
baseline_equivalent
period_shift_only
amplitude_shift_only
period_and_amplitude_shift
new_quasistable_cycle
fixed_point_convergence
divergent
numerically_unstable
not_reproducible
not_applicable
```

## 周期判定

既存コードに周期判定がある場合はそれを使用する。

ない場合は、状態距離または主要観測ベクトル

\[
X_n=
(L_A,L_B,N_{\mathrm{eff},A},N_{\mathrm{eff},B},
c_A,c_B,R_{\mathrm{eff}})
\]

について、自己相関と帰還距離を併用する。

新しい閾値は実行前に固定し、報告書に記載する。

結果を見て閾値を調整してはならない。

---

# 9. 物理量と診断量の区別

以下を物理的正本として保存する。

- raw全状態
- rawチャネルノルム
- 経路別振幅・ノルム
- \(R_{\mathrm{eff}}\)
- \(T_{\mathrm{eff}}\)
- \(c_A,c_B,\bar c\)
- \(p_B,p_F\)
- \(L\)
- \(N_{\mathrm{eff}}\)

以下は診断量であり、経路フラックスと混同しない。

- スペクトル類似度
- `B_to_A_transfer` という旧名称
- 最終状態距離
- 帰還誤差
- 自己相関

旧名称を出力する場合は必ず、

```text
B_to_A_transfer
(spectral_similarity_to_initial_B; not path flux)
```

と注記する。

---

# 10. 必須数値検査

全実行で以下を保存する。

\[
\varepsilon_U
=
\left|
|r|^2+|t|^2-1
\right|
\]

\[
\varepsilon_{rt}
=
\left|
r^*t+t^*r
\right|
\]

\[
\varepsilon_A
=
\left|
\|a_{\mathrm{raw}}\|^2
-
\left(
\|ra\|^2+\|tb\|^2+I_A
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
\|rb\|^2+\|ta\|^2+I_B
\right)
\right|
\]

さらに、

- 全ノルム保存残差
- 復調・再構成残差
- パリティ射影和残差
- NaN/Inf件数
- \(\theta_{\mathrm{eff}}\) 範囲違反
- rawノルムの最大・最小
- 既存正規化倍率

を保存する。

数値不安定が発生した条件を黙って除外してはならない。

---

# 11. 出力構成

最低限、次を作成する。

```text
stage_F_original_systemA_integration/
├── code/
│   ├── system_A_stage_F_copy.py
│   ├── state_dependent_scattering.py
│   ├── parity_demodulation.py
│   ├── run_stage_F_reproduction_gate.py
│   ├── run_stage_F_repeated_comparison.py
│   ├── run_stage_F_R_sweep.py
│   ├── run_stage_F_31_series_check.py
│   ├── analyze_stage_F_cycles.py
│   ├── build_stage_F_manifest.py
│   └── test_stage_F.py
├── source_copy/
│   ├── README.md
│   └── <original System A snapshot>
├── data/
│   ├── reference_hashes_before.json
│   ├── reference_hashes_after.json
│   ├── reference_hash_comparison.json
│   ├── stage_F_reproduction_gate.csv
│   ├── stage_F_repeated_collision_results.csv
│   ├── stage_F_run_summary.csv
│   ├── stage_F_R_sweep_results.csv
│   ├── stage_F_cycle_metrics.csv
│   ├── stage_F_31_series_results.csv
│   ├── stage_F_numerical_residuals.csv
│   ├── stage_F_summary.json
│   └── stage_F_final_states.npz
├── figures/
│   ├── C0_vs_reversed_C1_L_exchange.png
│   ├── C0_vs_reversed_C1_N_eff_exchange.png
│   ├── parity_and_R_eff_by_collision.png
│   ├── R_sweep_minimum_differences.png
│   ├── cycle_period_comparison.png
│   ├── normalization_vs_raw_update.png
│   └── return_error_31_series.png
├── reports/
│   ├── 00_scope_and_source_hashes.md
│   ├── 01_C0_reproduction_gate.md
│   ├── 02_integration_definition.md
│   ├── 03_repeated_scattering_comparison.md
│   ├── 04_R_sweep_comparison.md
│   ├── 05_cycle_and_quasistability_analysis.md
│   ├── 06_31_series_check.md
│   ├── 07_normalization_effect.md
│   ├── 08_numerical_invariants.md
│   └── Stage_F_report.md
├── logs/
└── manifest.json
```

適用不能な出力は空ファイルで偽装せず、報告書で `not_applicable` とする。

---

# 12. 報告書の分類規則

必ず以下を分離する。

## コード上の事実

実際の既存コードと今回のコピーから確認したこと。

## 再現結果

C0が既存System Aをどの精度で再現したか。

## 数学的帰結

式と保存則から必ず従うこと。

## モデル定義

反転Candidate 1、二層表現、正規化系列など今回採用した定義。

## 数値観測

実行により得た周期、振幅、差、谷、帰還誤差。

## 作業仮説

偶数倍音・奇数倍音とボゾン型・フェルミオン型の物理対応。

## 未導出

\(\kappa\)、包絡 \(\rho\)、パリティ読出し相互作用の起源など。

## 既存結果との一致・不一致

過去の局在性交換、準安定振動、R依存性、31系列がどう変わったか。

## 棄却または保留

数値不安定、再現不能、適用不能となった条件。

---

# 13. 重要な解釈制限

以下の結論を自動的に書いてはならない。

- 「フェルミオンが導出された」
- 「パウリ排他原理を証明した」
- 「自然界のボゾン・フェルミオン散乱を再現した」
- 「31周期を説明した」
- 「観測結果と一致した」

このStageで言えるのは、既存System Aの数値現象へ反転Candidate 1を接続したとき、

- 何が保存されたか
- 何が変化したか
- どの既存現象が維持・変形・消失したか

までである。

---

# 14. 原本保全

実行前後で、参照した全既存ファイルについて以下を記録する。

```text
path
size
mtime
sha256
```

前後差が一つでもあれば、Stage F失敗として停止する。

新規ファイルだけを `manifest.json` に記録する。

既存ファイルのタイムスタンプを変更する操作も禁止する。

---

# 15. 実行順序

必ず次の順で進める。

1. 対象コード・データの特定
2. SHA-256事前記録
3. Stage F専用コピー作成
4. 単体テスト
5. C0再現ゲート
6. 再現成功時のみF-A
7. F-A成功時のみF-B
8. F-B終了後にF-Cの適用可能性判定
9. 適用可能な場合のみF-C
10. 数値残差集計
11. SHA-256事後記録・比較
12. 全報告書・manifest作成
13. 停止

C0再現ゲートが失敗した場合は、F-A以降を実行しない。

---

# 16. 最終停止条件

以下を完了した時点で必ず停止する。

- C0再現ゲート
- 反転Candidate 1の既存System A実験コピーへの限定統合
- 既存代表条件の反復比較
- 既存R掃引との比較
- 31系列の適用可能性判定および可能なら比較
- existing_normalization対raw_update比較
- 数値残差確認
- 原本SHA-256前後一致
- `Stage_F_report.md`
- `manifest.json`

停止時に明示する。

```text
Stage F 完了。
既存System A / System B原本は変更していない。
反転Candidate 1は独立実験コピーにのみ統合した。
Candidate 2・3は追加していない。
N体系へは組み込んでいない。
論文本文は変更していない。
人間の承認待ち。
```
