# Stage A 実行前固定仕様書

## 文書状態

- 作成日: 2026-07-26
- 対象: 第8論文 本実験 Stage A（N=5）のみ
- 状態: **未承認・実装開始禁止**
- この文書の作成時点で実施したもの: 既存ファイルの読み取り、関数監査、SHA-256計算、既存成果物の所在確認
- この文書の作成時点で実施していないもの: Pythonの新規作成、既存Pythonの編集、Python実行、実験、後処理、既存データの再生成

本書の原仕様は次のファイルである。

```text
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/Claude_Code_第8論文_多階層因果検証_段階実行仕様書.md
SHA-256: 08df7fb49c8ecf9ca15e0fa67ba9ce621bd220c63e606f7bebaa47fe738a56e9
```

本書では、既存コードまたは原仕様に機械的定義がない値を推測で補わない。該当箇所は **未定義（実装阻止）** とする。未定義事項が人間承認によって解消されるまで、Stage Aのコード作成・編集・実行を開始しない。

---

## 1. 使用する既存コードの確定

### 1.1 原本力学・親状態生成

#### A. 低ランク原本エンジン

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py
SHA-256:
ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d
```

使用対象となる既存定義:

| 役割 | 既存定義 |
|---|---|
| 辺構築 | `build_edges(n)` |
| 低ランク系 | `LowRankSystem` |
| 位相から生成子情報を構築 | `LowRankSystem.set_theta(theta)` |
| \(Kz\) | `LowRankSystem.kmatvec(z)` |
| \(W^Tz\), \(Wy\) | `LowRankSystem.wt(z)`, `LowRankSystem.w(y)` |
| 厳密 \(\sigma\) スペクトル | `LowRankSystem.sigma_spectrum()` |
| \(\sigma_{\max}\) warm-start推定 | `LowRankSystem.sigma_max_power(wp, iters=3)` |
| Cayley更新 | `LowRankSystem.cayley_step(z, sigma)` |
| float64親状態生成 | `make_parent(sys_lr, rng, iters, beta, tol, restarts)` |
| 親固有モード残差 | `_eigenmode_residual(sys_lr, v)` |
| 第7論文の明示seed生成 | `zero_closure_kernel_seed(sys_lr, rng)` |

Stage Aで維持する原本Cayley更新は

\[
Z_{\tau+1}
=
(I-\gamma K/\sigma)^{-1}(I+\gamma K/\sigma)Z_\tau,
\qquad
\gamma=\tan(\pi/144)
\]

であり、既存関数 `LowRankSystem.cayley_step` を変更しない。Series 1およびSeries 3の状態更新で使用する。Series 2には任意精度版の承認済み既存関数がないため、現時点では使用関数を確定できない。

#### B. 第7論文の分裂量raw生成・増幅率計算

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_spontaneous_splitting_largeN_v1.py
SHA-256:
13baf6f5158c53ee92d1a08a0a5b60832424a83222b0601563a2676a392515ac
```

既存関数:

- `run(n, delta, seed, cap, after, tol)`
- \(f\) の既存計算:
  \[
  f=\frac{\|Z-P_1Z\|^2}{\|Z\|^2}
  \]
- crossingの既存計算: `f > 0.05` を初めて満たすstep
- 増幅率の既存計算: `run` 内のインライン処理。区間
  \[
  \max(10f_0,10^{-300})<f<10^{-3}
  \]
  の点が5点以上ある場合に、`np.polyfit(step, log(f), 1)` の傾きを返す。

このコードには、局所回帰窓、決定係数基準、持続step数、増幅停止区間、停止時刻信頼区間の機械定義はない。したがって、この増幅率計算だけを第8論文のイベント判定へ流用してはならない。

#### C. 準安定系列ラッパー

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_metastable_series_v1.py
SHA-256:
ee99154d54247197d52eb57bb6d9dcc064287b1f1c7a77ae2a657e31dffcc7c0
```

既存関数:

- `verify_core_hashes()`
- `run_one(n)`
- `make_plot()`

既存の準安定要約は `crossing + 100` 以後の \(f\) の中央値・最小・最大・5/95%点を計算するが、データに基づく「準安定開始時刻」を判定する関数ではない。

### 1.2 親平面・方向・q・rank

#### D. N=5の固定親平面分解

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_exact_v1.py
SHA-256:
9cf28ca8c0d2ad8fac2f0f6dae045248695247c5809c21ccb2069ef91a94ab76
```

使用する既存関数:

- `parent_plane_split_exact(sys_lr, v, round_digits=6)`

返り値:

- `p1_sigma`
- `B_p1`: 初期支配平面 \(P_1\) の実正規直交基底
- `B_rot`: その他の回転平面の実正規直交基底
- `spectrum`

N=5では近似版 `parent_plane_split_approx` を使用しない。

#### E. 瞬時支配平面、q1〜q4、rank_Q

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/code/run_n300_dimension_saturation_v2.py
SHA-256:
229938a66631057426f187ed80b17de08cfcb9107cfe509c30f5bbdcca3a03e6
```

使用する既存関数:

| 役割 | 既存関数 |
|---|---|
| \(G=W^TW\) の縮約 | `gram_reduce(sys_lr, Z, tau_G=1e-12)` |
| 瞬時支配回転平面 | `dominant_plane(sys_lr, gr)` |
| \(Q=[B_0\mid B_{\rm dom}]\) の特異値 | `qsv4(B0, Bd)` |
| 親平面と瞬時平面の主角 | `principal_angles(A, B)` |
| raw生成 | `run(n)` |
| 既存q図生成 | `make_figures(n, outdir, figdir, crossing)` |

`qsv4` の機械定義は

\[
Q=[B_0\mid B_{\rm dom}],\qquad
C_4=Q^TQ,\qquad
q_j=\sqrt{\lambda_j(C_4)}
\]

である。

`rank_Q` の独立した関数は存在しない。既存コードでは `run` 内で

\[
\operatorname{rank}_Q
=
\#\{j:q_j>10^{-8}q_1\}
\]

をインライン計算している。

注意: この \(Q\) は4列なので、既存の同一定義から得られるのは `q1`〜`q4`だけである。原仕様が必須列とする `q5`〜`q8` はこの定義からは生成できない。`q5`〜`q8` は **未定義（実装阻止）** である。

#### F. 第7論文の五成分占有・新方向3/4

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py
SHA-256:
fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503
```

使用する既存関数:

- `occ(B, Z)`: 実基底 \(B\) への状態占有
- `build(n)`: 第7論文初期状態・固定基底構築
- `s4_new_dirs(B0, Bdom)`: \(S_4=\operatorname{orth}[B_0|B_{\rm dom}]\) の \(B_0\) 直交補2列
- `align_2d(f_prev, f_new)`: 2次元縮退基底の時間連続整列
- `run(n)`: `paper7_long_timeseries.csv` のraw生成

`direction_3_occupation` と `direction_4_occupation` は、`s4_new_dirs` の2列を固定 `B_rot` へ射影・QR直交化し、`align_2d` で前時刻へ整列した後、`occ` で計算される。

#### G. 第7論文の追加方向萌芽に使われた横摂動解析

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_transverse.py
SHA-256:
ac1073bea329971de3ff4c2fd1588d926029a8502c21e8cc01f406acb86ad60b
```

既存関数:

- `evolve(sys_lr, Z, wp)`
- `s4_basis(sys_lr, B0, Z)`
- `perp(S4, d)`
- `run(n)`

既存コードは `crossing + 3000` を摂動開始 `t0` とし、3方向閉包外の微小摂動に対するBenettin型有限時間横成長率を計算する。`lambda_transverse` の全結果が正なら `further_splitting_lambda_max_positive` と分類する。

これは自然軌道上で追加方向が出現した時刻を返す関数ではない。第7論文の「追加方向萌芽」は、自然軌道上の第5/第6有限占有ではなく、外部から与えた横摂動が増幅可能であるという解釈である。単一の「追加方向萌芽時刻」は既存コードに存在しない。

また、第7論文完成版は、横摂動解析についてwarm-start内部状態の同期不足が定量値へ混入し得ると明記している。この定量値を無条件に第8論文の基準値として扱わない。

### 1.3 量子化・再射影

#### H. 第2予備実験の量子化演算

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/preliminary_02_resolution_scan/code/run_initial_resolution_scan_v1.py
SHA-256:
d82c5a105124a509bdbeeeeaf554aea6e846e6f54cf24b6d2d7073f20a4072b8
```

使用候補となる既存関数:

- `Q_delta(Z, d)`:
  \[
  Q_\Delta(Z)=\Delta\operatorname{round}(\Re Z/\Delta)
  +i\Delta\operatorname{round}(\Im Z/\Delta)
  \]
  `round` はhalf-to-even。
- `observe(...)` は `q3`, `q4`, `rank_Q` を出すが、原仕様の全必須列を満たさないため、そのままStage Aの観測器には使用できない。

`delta_actual` によるNスケーリングは今回使用しない。Stage A原仕様に列挙された \(\Delta\) を直接使用する。

#### I. polar retraction

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_transverse_stability_v1.py
SHA-256:
c868884ba9ab5a10c5b398a57b347afca167e28039425ba3550fc2cda1f65566
```

使用候補となる既存関数:

- `retract(W)`: \([\sqrt2\Re W,\sqrt2\Im W]\) の2列フレームをpolar retractionで直交正規化し、\(Z^TZ=0,\ Z^\dagger Z=1\) の多様体へ戻す。
- `closure_errors(Z)`

Series 3に限り、更新順序を

\[
\text{Cayley}\rightarrow\text{測定}\rightarrow Q_\Delta
\rightarrow\text{測定}\rightarrow\text{retract 1回}
\rightarrow\text{測定}
\]

とする。Series 1/2には量子化とretractionを挿入しない。

### 1.4 図表生成コード

#### J. 第7論文長時間図

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/make_paper7_figures.py
SHA-256:
273fec25c2e4c30c4561a7506ef373829e374d8b0d2eeb813bb7f4e1c06a8000
```

既存関数:

- `load5(n)`
- `crossing(n)`
- `fig1()`
- `fig23()`
- `load_trans(n)`
- `fig_transverse()`
- `fig_lambda_vs_N()`

#### K. q3/q4のN比較図

```text
絶対パス:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/code/make_saturation_comparison.py
SHA-256:
51cee54314e01160156d4fd0a3f8082c42a09a69d05672479aec6b2127c7d3b6
```

Stage Aの図1〜9を生成する既存コードは存在しない。上記コードは第7論文再現図の基準としてのみ使用する。Stage A図表コードを新規作成する場合は、本書承認後でなければならない。

### 1.5 第7論文本文・報告書

```text
第7論文日本語完成版:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/paper7_complete_JP/論文7_N体関係波閉鎖系における三方向空間の創発_日本語完成版.md
SHA-256:
fc721518d0fde4789551966a385d711be4c6213c5d01b14457f0f960bc68fc94

第7論文長時間・横安定性報告:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/reports/paper7_longtime_and_transverse_stability_report.md
SHA-256:
c0ba216b6b1e6346062ab166beb342bbbb2c0377df44d9322a9d120217ea20d9

方向飽和検査報告:
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/reports/observation_report_saturation_N5_40_300.md
SHA-256:
e2b28343f68ca9a81839d58cc4405c143236eadc6e39c477f030eee859281249
```

### 1.6 使用禁止・隔離対象

次の未追跡コードはClaude Codeの中断された試行であり、Stage Aの既存正常コードとして使用しない。

| 絶対パス | SHA-256 | 使用禁止理由 |
|---|---|---|
| `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/main_experiment_multiscale_causality/code/mp_engine_v1.py` | `f2552be19c9121294a626d9978c4ae1f147a9ac068aee12743e7888ccfbb1649` | `cast_parent`、未承認のNewton/LM型 `refine_parent` を含む。承認済み仕様ではない。 |
| `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/main_experiment_multiscale_causality/code/stageA_common_v1.py` | `b182e81623a8ddb31d37d961dfacc68fe8115b4744b660bf48b4c225eaf17b46` | 第7論文にない方向・flux・noise floor定義を独自導入している。承認済み仕様ではない。 |

これらを自動修正、流用、上書きしない。削除・移動も人間の明示承認なしに行わない。

---

## 2. 第7論文N=5の完全再現方法

### 2.1 実行前ゲート

新しい下位階層実験より先に、第7論文N=5の再現を独立工程として行う。既存成果物を上書きしないため、承認後に原本コードと必要入力を隔離した再現用スナップショットへコピーし、そのコピー内だけで実行する。コピー前後に1.1節のSHA-256を検査する。

予定する再現順序:

1. `run_metastable_series_v1.py run 5`
2. `run_n300_dimension_saturation_v2.py 5`
3. `run_paper7_5color_timeseries.py 5`
4. `run_paper7_transverse.py 5`
5. `make_paper7_figures.py`
6. `make_saturation_comparison.py`
7. 既存CSV・JSON・図・報告書と比較
8. 再現報告を作成し、失敗した場合は停止

この手順は本書での予定であり、まだ実行していない。

### 2.2 比較する既存データ

| 対象 | 絶対パス | SHA-256 |
|---|---|---|
| N=5 分裂量 \(f\) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/metastable_series_result_v1/fcurve_N00005_delta1e-15_seed0.csv` | `9220c5f3c1f570c8a52ea24a3cdd95568354cea0943d9bee7d8ed20316d3a9d0` |
| N=5 q1〜q4・rank_Q | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/raw/N00005_dimension_saturation_v2/q_svd_N00005.csv` | `7c16a364c6cc9145293c2625dfe4ebb1f9962655d212679188215e8fad5e5155` |
| N=5 五成分長時間時系列 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/paper7_long_timeseries.csv` | `efeaf9dab753c057ad0c6109b9e4a8919f8d8db1249da186658bfed9fda784e3` |
| N=5 横摂動時系列 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/transverse_stability_timeseries.csv` | `5a3e1e408122b8c7267236f0b708517027778a0abc3de98fff867c82e26d2459` |
| N=5 q/rank診断 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/diagnostics/N00005_saturation.json` | `a4b94ed5f22a2642ee26bc4b2eb744dd16f0ec8c0cc8d50618aae52e55d201b4` |
| N=5 五成分メタ | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/summary/N00005_5color_meta.json` | `5e6d21778778f4f79a4e02d61b59261d0ae3f6fdb31deefccbf7aa30c0edc0ba` |
| N=5 横摂動メタ | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/summary/N00005_transverse_meta.json` | `4bcc284e38b2eb27214afb780e88ab667b37cc5d04c1044a8233e8bd40455c51` |

### 2.3 比較する既存図

| 図 | 絶対パス | SHA-256 |
|---|---|---|
| 分裂量 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/figures/figure1_N00005.png` | `cb1df38404382f3cb36a9dc0347f1bf7736118108da5d5f4f79045d549a351cc` |
| 五成分占有 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/figures/figure2_N00005_5color.png` | `6ae19b7c7ce4930648e4820d3666e95389d99eb8ef3ae11106b60816fdf0154b` |
| 五成分対数図 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/figures/figure3_N00005_5color.png` | `ac7829a4504a29dbcda33c7121a37741f9d7caabdc7b62e08cea421b80f23f83` |
| 横摂動応答 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/figures/transverse_growth_N00005.png` | `5bac0537130bf76e1cb697f4cc1852ab34319c7a48b1a9a5b3aa44bebfd7073d` |

### 2.4 既存結果と、現時点で比較可能な量

| 要求項目 | 既存の機械的基準 | N=5既存値 | 再現可否 |
|---|---|---:|---|
| 幾何級数的拡大開始時刻 | 独立定義なし。既存の「立上がり」は \(f>0.05\) の初回crossing | crossing = 1167 | 要求された意味では未定義 |
| 2方向+1方向の成立時刻 | `rank_Q=4` は定義済みだが、持続・方向連続性を含む成立時刻関数なし | 報告上はcrossing後の記録761点でrank_Q=4 | 時刻は未定義 |
| 幾何級数的拡大停止時刻 | 独立関数・閾値・信頼区間なし | 本文では有限準安定域への移行と記述 | 未定義 |
| 準安定開始時刻 | `run_paper7_transverse.py` が `crossing+3000` をガードとして固定 | 4167 | データ判定ではなく固定ガード |
| 追加方向萌芽 | 横摂動の全条件で有限時間成長率が正かを分類 | \(\lambda_{\perp,\max}=2.9166086965\times10^{-3}\) | 有無は再現可能、時刻は未定義 |
| q3 | `qsv4` | 準安定代表 0.7514116652、終端 0.8139660184 | 比較可能 |
| q4 | `qsv4` | 準安定代表 0.6307332533、終端 0.6893091409 | 比較可能 |
| rank_Q | \(q_j>10^{-8}q_1\) の個数 | 終端4 | 比較可能。ただし底付近で固定閾値揺れあり |

第7論文N=5再現の数値許容誤差、CSV全行一致の基準、図の一致基準は既存仕様に数値で定義されていない。これらは **未定義（実装阻止）** である。

結論: 現時点では、要求された全イベントを同じ定義で再現する方法は確定していない。第7論文を再現できない状態では新実験へ進まない。

---

## 3. 高精度親の構築方法

### 3.1 float64で既存実装されている親アルゴリズム

既存 `make_parent` は、初期位相 \(\theta^{(0)}\) から次を反復する。

1. \(\theta^{(k)}\) で \(G=W^TW\) を構築する。
2. \(JG\) の固有対から、虚部が最も負の固有値 \(-i\sigma_{\max}\) の固有ベクトルを選ぶ。
3. 辺空間へ持ち上げ、
   \[
   v^{(k)}=\frac{Wu^{(k)}}{\|Wu^{(k)}\|}
   \]
   とする。
4. \(\theta_v^{(k)}=\arg v^{(k)}\) を求める。
5. \(\beta=0.5\) のとき、
   \[
   \theta^{(k+1)}
   =
   \arg\left[
   (1-\beta)e^{i\theta^{(k)}}
   +\beta e^{i\theta_v^{(k)}}
   \right]
   \]
   と更新する。

親残差は既存 `_eigenmode_residual` により

\[
\mu=\Re\{v^\dagger(iKv)\},
\qquad
r_{\rm parent}=\|iKv-\mu v\|_2
\]

と定義される。

第7論文N=5では、`make_parent(..., iters=1200, tol=1e-12)`、親PRNG seed `40260722+1000N=40265722` が使われ、既存要約の親残差は \(2.1398983601\times10^{-13}\) である。

### 3.2 80/128/256 bitの現状

| 項目 | 80 bit | 128 bit | 256 bit |
|---|---|---|---|
| 使用する承認済み既存関数 | 存在しない | 存在しない | 存在しない |
| 親生成アルゴリズム | 未定義 | 未定義 | 未定義 |
| float64親の単純キャスト | 禁止 | 禁止 | 禁止 |
| 同じ固定点・親枝の確認法 | 未定義 | 未定義 | 未定義 |
| 別固定点判定量 | 未定義 | 未定義 | 未定義 |
| 位相ゲージ固定 | 既存 `make_parent` には明示固定なし | 同左 | 同左 |
| 縮退固有値の扱い | 既存 `make_parent` には部分空間追跡なし | 同左 | 同左 |
| 親残差 | 上記 \(r_{\rm parent}\) を候補とするが高精度版未実装 | 同左 | 同左 |
| 目標親残差 | 未定義 | 未定義 | 未定義 |
| 閉鎖残差 | \(|v^Tv|\) を候補とするが目標値未定義 | 同左 | 同左 |
| float64親との比較量 | 未定義 | 未定義 | 未定義 |
| 構築失敗時 | 停止・原因報告・代替法へ移行しない | 同左 | 同左 |

同じ親枝を確認するために必要となり得る量として、位相整列後の状態重なり、親平面の主角、\(\sigma_j/\sigma_1\)、固有モード残差、位相差ベクトルが考えられる。しかし、どの量を使い、どの閾値で同一枝と判定するかは原仕様・第7論文に定義されていないため、本書では値を決めない。

### 3.3 新手法の禁止と承認ゲート

- Newton法、Levenberg–Marquardt法、別の固定点反復、連続追跡法、固有部分空間追跡法を自動導入しない。
- 新手法が必要な場合は、実装前に目的関数、未知変数、ゲージ条件、縮退処理、反復式、停止条件、枝判定量を数式で別文書に提示する。
- 人間の明示承認前に実装しない。
- 高精度親が失敗しても、float64親のキャスト、精度水準の削減、別固定点の採用、別ソルバーへの自動切替を行わない。

結論: 高精度親構築法は現時点で **未確定（実装阻止）** である。

---

## 4. 「方向」の定義

### 4.1 第7論文で機械定義された量

| 名称 | 定義 | 既存実装 | 空間方向との区別 |
|---|---|---|---|
| 生成子のrank | \(K=WJW^T\) または \(G=W^TW\) の数値rank。既存コードでは `gram_rank` 等 | `gram_reduce` の `r_G`、`sigma_spectrum` | 回転生成子のrankであり、幾何学的空間方向数ではない |
| 状態占有方向 | 実基底 \(B\) に対する \(\|\!B^T\Re Z\!\|^2+\|\!B^T\Im Z\!\|^2\) | `occ` | 状態ノルムの配分先 |
| 特異方向・瞬時支配平面 | \(iK\) の最大正固有値に対応する複素固有ベクトルの実部・虚部が張る2次元平面 | `gram_reduce`, `dominant_plane` | 生成子が与える回転平面 |
| \(Q\) | \(Q=[B_0|B_{\rm dom}(t)]\) | `qsv4` | 2つの実2次元平面の結合行列 |
| q3, q4 | \(Q\) の第3、第4特異値 | `qsv4` | 状態占有量ではない |
| rank_Q | \(\#\{j:q_j>10^{-8}q_1\}\) | `run_n300_dimension_saturation_v2.run` 内のインライン式 | 生成子rankとも空間座標rankとも異なる |
| 新方向3/4占有 | \(S_4=\operatorname{orth}[B_0|B_{\rm dom}]\) の \(B_0\) 直交補を `B_rot` へ射影し、連続整列した各1列への占有 | `s4_new_dirs`, `align_2d`, `occ` | q3/q4そのものではない |

### 4.2 「2方向+1方向」

第7論文完成版は、数値的には

\[
\operatorname{rank}Q:2\rightarrow4
\]

を測定する。一方、零二乗閉鎖を保つ反対称力学では新方向が共役対として成立するため、その一方を第三空間方向、他方を閉鎖共役方向と読む。このため「初期の実二方向 + 第三空間方向」が「2方向+1方向」に対応する。

これは「rank_Q=4だから幾何学的空間が4次元」という定義ではない。また、回転平面数、生成子rank、状態占有方向数、幾何学的空間方向数を同一視しない。

### 4.3 三方向成立

第7論文本文で機械的に明記される最低条件は

\[
q_3>10^{-8}q_1,\qquad
q_4>10^{-8}q_1,\qquad
\operatorname{rank}Q=4
\]

である。

ただし、第8論文原仕様が要求する方向連続性、持続時間、数値床S/N、複数精度再現性、量子化前後・再射影前後一致まで含む「成立」関数は第7論文コードに存在しない。第7論文の機械的rank条件と第8論文の成立条件を混同しない。

### 4.4 追加方向萌芽

第7論文の追加方向萌芽は、自然軌道上で第5または第6方向の有限占有帯が出現したことではない。自然軌道では55000 stepまで新しい有限占有帯は観測されていない。

第7論文が「萌芽」と呼ぶものは、準安定な \(S_4(t)\) の直交補へ外部微小摂動を与えたとき、複数seed・複数振幅で正の有限時間横成長率が得られたことである。したがって、既存コードには自然軌道上の「追加方向萌芽時刻」は存在しない。

新しい方向定義は本書では作らない。

---

## 5. イベント定義

| イベント | 数式・判定対象列 | 閾値 | 持続step | 状態 |
|---|---|---|---|---|
| 下位幾何級数的拡大開始 | \(g(t)=d\log a_{\rm outside}/dt\)、`local_log_growth`, `smoothed_log_growth`, 回帰 \(R^2\) | \(g\)有意性、\(R^2\)、床比はいずれも未定義 | 未定義 | **実装阻止** |
| 下位幾何級数的拡大停止 | \(g\to0\)、指数回帰離脱、方向配分再編、flux釣合い、準安定振動開始 | 全て未定義 | 未定義。信頼区間法も未定義 | **実装阻止** |
| 第3方向萌芽 | \(q_3>c_{\rm sprout}\eta_{\rm noise}\) かつ短期追跡可能 | \(c_{\rm sprout}\)、\(\eta_{\rm noise}\) の式、連続性閾値が未定義 | 「短期間」が未定義 | **実装阻止** |
| 第3方向成立 | \(q_3/q_1\)、適応床S/N、方向連続性、再現性、精度再現性、操作別一致 | 全て未定義。第7論文の固定rank閾値だけでは不可 | 未定義 | **実装阻止** |
| 上位幾何級数的拡大開始 | 下位と同じ局所回帰定義を上位イベントへ適用する必要がある | 未定義 | 未定義 | **実装阻止** |
| 上位幾何級数的拡大停止 | 下位停止と同じ複合判定を上位イベントへ適用する必要がある | 未定義 | 未定義 | **実装阻止** |
| 準安定開始 | 停止後に有界振動と方向占有の持続を検出する必要がある | 振幅幅、傾き、持続長が未定義 | 未定義 | **実装阻止** |
| 追加方向萌芽 | 第7論文では外部横摂動の `lambda_transverse>0` 分類。第8論文の自然軌道イベントとしては `q4/q1`, `q5/q1` 等を想定するが同一定義ではない | 未定義 | 未定義 | **実装阻止** |

参考として、第7論文の `f>0.05` crossing、`rank_Q=4`、`crossing+3000` は既存値である。しかし、これらを第8論文の拡大開始、拡大停止、方向成立、準安定開始へ自動的に読み替えない。

---

## 6. 実行上限時刻・保存量・見積り

### 6.1 固定状況

| 項目 | Stage A予定 | 固定状況 |
|---|---|---|
| N | 5 | 固定 |
| M | 10 | 固定 |
| 最大step | 未定義 | **実装阻止** |
| 通常停止 | 固定最大step到達 | 最大step未定義 |
| 異常停止 | NaN/Inf、親構築失敗、閉鎖・ノルム基準超過、出力不全 | 数値健全性閾値未定義 |
| 保存間隔 | 初期から終端まで毎step | 固定 |
| 毎step保存区間 | \(0\le t\le\text{max_step}\) | max_step未定義 |
| 測定点 | 8節の6点 | 固定。ただしSeries 1/2の量子化・再射影点はN/A |
| 状態ベクトル | 少なくともstep開始、Cayley後、量子化後、再射影後の状態変更点 | 保存形式と精度保持方式が未定義 |
| 追加方向萌芽3回条件 | 原仕様に存在 | 萌芽イベント定義が未定義のため停止条件に使用不能 |

### 6.2 容量・時間の既存参考値

既存N=5 float64走行には次の実績がある。

- 力学中心の既存走行: 21167 step、実行時間 1.6686201096秒。
- 第2予備実験の一例: 2500 step、保存1301行、runディレクトリ約548 KiB。
- 既存q観測の先頭10点: 1観測あたり概ね \(1.7\times10^{-4}\)〜\(3.1\times10^{-4}\) 秒。

これらはStage Aの6測定点、約60必須列、全step状態保存、高精度演算を含まないため、Stage Aの見積りへ直接外挿しない。

| 精度系列 | 予想ファイル容量 | 予想実行時間 |
|---|---|---|
| float64 | max_step・状態保存形式未定義のため算出不能 | 観測器完成前のため算出不能 |
| 80 bit | 高精度実装・max_step未定義のため算出不能 | ベンチマーク未実施のため算出不能 |
| 128 bit | 同上 | 同上 |
| 256 bit | 同上 | 同上 |
| 各Delta | max_step・状態保存形式未定義のため算出不能 | 観測器完成前のため算出不能 |

結果を見て上限値を変更しないためには、max_step、高精度backend、状態保存形式、容量上限、時間上限を実装前に人間承認で固定する必要がある。

---

## 7. Stage A全run一覧

共通出力root:

```text
/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/main_experiment_multiscale_causality/stage_A_N5/raw/
```

全runで、明示的物理seed、metastable seed、random kick、external noiseは使用しない。`seed` 列の `40265722` はfloat64既存親生成のPRNG seedであり、状態へ加える明示seedではない。exec 1/2は同じ条件を独立構築・独立実行する。

| run_id | series | precision | Delta | exec_id | max_step | seed | 初期状態生成法 | 量子化 | 再射影 | 出力先 |
|---|---|---:|---:|---:|---|---|---|---|---|---|
| A_S1_f64_e1 | Series 1 | float64 | OFF | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | なし | なし | `.../raw/A_S1_f64_e1/` |
| A_S1_f64_e2 | Series 1 | float64 | OFF | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | なし | なし | `.../raw/A_S1_f64_e2/` |
| A_S2_p080_e1 | Series 2 | 80 bit | OFF | 1 | **未定義** | 未定義 / explicit seedなし | **高精度親法未定義** | なし | なし | `.../raw/A_S2_p080_e1/` |
| A_S2_p080_e2 | Series 2 | 80 bit | OFF | 2 | **未定義** | 未定義 / explicit seedなし | **高精度親法未定義** | なし | なし | `.../raw/A_S2_p080_e2/` |
| A_S2_p128_e1 | Series 2 | 128 bit | OFF | 1 | **未定義** | 未定義 / explicit seedなし | **高精度親法未定義** | なし | なし | `.../raw/A_S2_p128_e1/` |
| A_S2_p128_e2 | Series 2 | 128 bit | OFF | 2 | **未定義** | 未定義 / explicit seedなし | **高精度親法未定義** | なし | なし | `.../raw/A_S2_p128_e2/` |
| A_S2_p256_e1 | Series 2 | 256 bit | OFF | 1 | **未定義** | 未定義 / explicit seedなし | **高精度親法未定義** | なし | なし | `.../raw/A_S2_p256_e1/` |
| A_S2_p256_e2 | Series 2 | 256 bit | OFF | 2 | **未定義** | 未定義 / explicit seedなし | **高精度親法未定義** | なし | なし | `.../raw/A_S2_p256_e2/` |
| A_S3_d1e-04_e1 | Series 3 | float64 | \(10^{-4}\) | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-04_e1/` |
| A_S3_d1e-04_e2 | Series 3 | float64 | \(10^{-4}\) | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-04_e2/` |
| A_S3_d1e-06_e1 | Series 3 | float64 | \(10^{-6}\) | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-06_e1/` |
| A_S3_d1e-06_e2 | Series 3 | float64 | \(10^{-6}\) | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-06_e2/` |
| A_S3_d1e-08_e1 | Series 3 | float64 | \(10^{-8}\) | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-08_e1/` |
| A_S3_d1e-08_e2 | Series 3 | float64 | \(10^{-8}\) | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-08_e2/` |
| A_S3_d1e-10_e1 | Series 3 | float64 | \(10^{-10}\) | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-10_e1/` |
| A_S3_d1e-10_e2 | Series 3 | float64 | \(10^{-10}\) | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-10_e2/` |
| A_S3_d1e-12_e1 | Series 3 | float64 | \(10^{-12}\) | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-12_e1/` |
| A_S3_d1e-12_e2 | Series 3 | float64 | \(10^{-12}\) | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-12_e2/` |
| A_S3_d1e-14_e1 | Series 3 | float64 | \(10^{-14}\) | 1 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-14_e1/` |
| A_S3_d1e-14_e2 | Series 3 | float64 | \(10^{-14}\) | 2 | **未定義** | parent RNG 40265722 / explicit seedなし | `make_parent`, \(Z_0=v\) | `Q_delta` | `retract` 1回 | `.../raw/A_S3_d1e-14_e2/` |

run数確認:

- Series 1: 2 run
- Series 2: 80/128/256 bit × 2 = 6 run
- Series 3: 6 Delta × 2 = 12 run
- 合計: **20 run**

全精度・全Delta・各2回は一覧上欠落していない。ただし、max_stepと高精度親が未定義なので実行条件はまだ固定完了していない。

---

## 8. 各stepの測定点

| 測定点 | 状態 | 測定対象 | 状態変更 | 備考 |
|---|---|---|---|---|
| step開始時 | \(Z_t\) | \(Z^TZ\), \(Z^\dagger Z\), \(a_{\rm outside}\), q, rank候補, 方向基底・連続性、成分統計 | なし | 次の生成子構築前 |
| 生成子構築後 | \(Z_t\) | \(G\), \(\sigma\), 主要特異方向、生成子低ランク表現、固有対残差 | なし | `set_theta` は力学上必要な生成子構築。観測結果は更新へ戻さない |
| Cayley更新直後 | \(Z_c\) | 全必須状態量、更新差 \(\|Z_c-Z_t\|\)、閉鎖・ノルム | Cayleyのみ | Series 1/2/3共通 |
| 量子化直後 | \(Z_q=Q_\Delta(Z_c)\) | 全必須状態量、`quantization_l2` | Series 3のみ量子化 | Series 1/2はN/Aとして記録 |
| 再射影直後 | \(Z_r=\operatorname{retract}(Z_q)\) | 全必須状態量、`retraction_correction_l2` | Series 3のみretraction | Series 1/2はN/Aとして記録 |
| 次step直前 | 実際に次stepへ渡す状態 | 全必須状態量、イベントフラグ | なし | Series 1/2は \(Z_c\)、Series 3は \(Z_r\) |

非干渉規則:

1. SVD、rank、方向追跡、event判定、plot用量を `cayley_step`、量子化、retractionへ入力しない。
2. 力学更新に使うのは、現在状態から原本どおり構築した生成子と \(\sigma_{\max}\) だけである。
3. 観測は状態のコピーまたは読み取り参照に対して行う。
4. 観測失敗時に状態を修正して継続しない。停止して報告する。
5. 既存 `align_2d` は表示・占有の連続整列に限り、状態更新へ使わない。

未解決:

- `direction1_continuity`〜`direction4_continuity` の既存同一定義がない。
- `flux_parent_to_d3` 等の既存第7論文定義がない。
- `eta_noise(t)` の既存式がない。
- `q5`〜`q8` の第7論文同一定義がない。

---

## 9. 出力図表の固定状況

各図はPNGとSVGを作る。入力の基本ファイル名は各runの `timeseries.csv` とする。イベント集約は `tables/event_times.csv` を想定するが、イベント定義が未確定なので生成可能状態ではない。

### 9.1 図1〜9

| 図 | 入力CSV | 使用列 | 横軸 | 縦軸 | スケール | イベント表示 | ファイル名 | 固定状況 |
|---|---|---|---|---|---|---|---|---|
| 図1 全時間域概要 | 全run `timeseries.csv` | `step`, `a_outside`, `q2`〜`q5`, event列 | step | \(a_{\rm outside},q_2,q_3,q_4,q_5\) | y対数 | 全イベント縦線 | `fig01_full_time_overview.{png,svg}` | q5・event未定義 |
| 図2 初期下位領域拡大 | `timeseries.csv` | `step`, `log_a_outside`, `local_log_growth`, `smoothed_log_growth`, `eta_noise` | step | \(\log a\)、増幅率、回帰線、床 | \(\log a\)は線形表示、振幅対応は対数 | lower growth start/end | `fig02_lower_growth.{png,svg}` | 回帰窓・床・event未定義 |
| 図3 方向成立 | `timeseries.csv` | `q2_over_q1`〜`q5_over_q1`, `rank_noise_adaptive`, `rank_persistent`, `eta_noise` | step | q比、床、成立度 | y対数 | sprout/established、持続区間 | `fig03_direction_establishment.{png,svg}` | 閾値・持続長未定義 |
| 図4 方向連続性 | `timeseries.csv` | `direction3_continuity`, `direction4_continuity`, `direction3_rotation_angle`, `direction4_rotation_angle` | step | 内積・回転角 | 線形 | sprout/established | `fig04_direction_continuity.{png,svg}` | 方向追跡定義未定義 |
| 図5 拡大停止同時性 | `timeseries.csv`, `event_times.csv` | `smoothed_log_growth`, rank/成立度, flux列 | step | 増幅率、成立度、flux | 線形 | 停止信頼区間、成立時刻 | `fig05_growth_stop_coincidence.{png,svg}` | 停止・flux未定義 |
| 図6 追加方向萌芽 | `timeseries.csv` | `q4_over_q1`, `q5_over_q1`, event列、振動位相 | step | q比・振動位相 | q比対数、位相線形 | metastable start、sprout | `fig06_additional_direction_sprout.{png,svg}` | q5・振動位相列・event未定義 |
| 図7 更新操作別比較 | `timeseries.csv` | `measurement_point`, `q3`, `q4`, `quantization_l2`, `retraction_correction_l2` | step | q3/q4・補正量 | 原則y対数 | 操作位置を色・線種で区別 | `fig07_update_point_comparison.{png,svg}` | 列は固定可能、q観測精度版未定義 |
| 図8 精度比較 | 全run `timeseries.csv` | `precision_bits`, `Delta`, `step`, `a_outside`, q列、event列 | stepまたは承認済み平行移動時間 | 振幅・q | y対数 | 各精度event | `fig08_precision_comparison.{png,svg}` | 時間平行移動規則未定義 |
| 図9 第7論文とのスケール比較 | Stage A `event_times.csv`, 第7論文基準CSV | event列、特徴時間、正規化振幅 | 正規化時間 | 正規化振幅・event列 | 線形/対数は量ごと | 上位・下位event | `fig09_paper7_scale_comparison.{png,svg}` | 特徴時間・第7論文eventが未定義 |

### 9.2 表1〜7

| 表 | 入力 | 列 | ファイル名 | 固定状況 |
|---|---|---|---|---|
| 表1 run一覧 | 全 `run_config.json`, `run_diagnostics.json` | run_id, N, precision, Delta, exec_id, runtime, status, stop_reason | `table01_runs.csv`, `table01_runs.md` | max_step等未定義 |
| 表2 イベント時刻 | `event_times.csv` | 8イベント、下限・代表・上限 | `table02_event_times.csv`, `table02_event_times.md` | event定義未定義 |
| 表3 時間差 | 表2 | \(t_{3,\rm est}-t_{\rm lower,end}\), \(t_{4,\rm sprout}-t_{\rm meta}\), 正規化値 | `table03_time_differences.csv`, `table03_time_differences.md` | event・特徴時間未定義 |
| 表4 増幅率 | `timeseries.csv`, 回帰区間表 | 下位/上位傾き、CI、\(R^2\)、精度、Delta | `table04_growth_rates.csv`, `table04_growth_rates.md` | 回帰仕様未定義 |
| 表5 方向統計 | `timeseries.csv` | q比最大、持続時間、連続性、再現性、床比 | `table05_direction_statistics.csv`, `table05_direction_statistics.md` | 閾値・追跡定義未定義 |
| 表6 数値健全性 | `timeseries.csv`, diagnostics | closure最大、norm最大、量子化補正、再射影補正、exec差 | `table06_numerical_health.csv`, `table06_numerical_health.md` | 合否閾値未定義 |
| 表7 仮説判定 | 表2〜6 | run_id, H1/H2/H0/未判定, 根拠 | `table07_hypothesis.csv`, `table07_hypothesis.md` | event定義確定後のみ作成可能 |

---

## 10. 失敗時の停止規則

次の場合は自動修正せず停止する。

| 失敗条件 | 停止時に報告するもの | 禁止する自動対応 |
|---|---|---|
| 高精度親が同じ固定点・親枝へ収束しない | 精度、反復履歴、残差、枝判定量、固有値スペクトル | 別固定点採用、float64親キャスト、別ソルバー切替 |
| 親残差が目標へ達しない | 到達残差、目標残差、反復数、停滞状況 | 目標緩和、反復法変更 |
| 第7論文N=5再現に失敗 | 不一致ファイル、列、step、最大差、環境 | 新実験開始、既存データ修正 |
| 必要な既存関数が見つからない | 必要な量、探索済みファイル、欠落関数 | 代替関数の即時作成 |
| 原仕様と既存コードが矛盾 | 仕様箇所、コード箇所、影響run/列 | 独自解釈による片方の採用 |
| 出力列が不足 | 欠落列、影響図表・判定 | 空欄を0や推定値で補完 |
| 数値健全性基準を超過 | run_id、最初の超過step、量、値 | retraction追加、閾値緩和、状態修正 |
| 実行時間・容量が見積りを大幅超過 | 実績、見積り、比率、残容量 | 保存間引き、精度削減、run削減 |
| NaN/Infまたは線形代数例外 | run_id、step、measurement point、例外 | 条件数変更、別アルゴリズム切替 |
| 測定が状態へ影響した可能性 | 影響経路、比較結果 | そのまま継続 |

「大幅超過」の倍率、数値健全性閾値、目標親残差は未定義であり、実装前に固定する。

停止時は、原因と不足情報だけを報告する。代替実装、パラメータ変更、次run、次Stageを開始しない。

---

## 11. 未解決事項一覧

以下はすべて実装前に人間判断または明示承認が必要である。

1. 第7論文N=5再現の数値許容誤差、CSV一致基準、図一致基準。
2. 第7論文の「幾何級数的拡大開始」と \(f>0.05\) crossingの関係。
3. 第7論文の「幾何級数的拡大停止」の機械的定義と時刻。
4. 第7論文の「2方向+1方向成立」の持続・連続性を含む時刻定義。
5. 第7論文の「準安定開始」を `crossing+3000` とみなすか、データ判定を新設するか。
6. 第7論文の外部横摂動による萌芽と、第8論文の自然軌道上の追加方向萌芽を同一イベントと比較できるか。
7. 80/128/256 bitの数値backend。
8. 80/128/256 bitの親生成アルゴリズム。
9. 高精度親の位相ゲージ固定法。
10. 高精度親の縮退固有値・縮退部分空間の追跡法。
11. 高精度親がfloat64親と同一固定点・同一枝であることの判定量と閾値。
12. 各精度の親残差目標。
13. 各精度の親閉鎖残差・ノルム残差目標。
14. 高精度 \(\sigma_{\max}\)、Cayley解、Gram縮約、q計算の承認済み実装。
15. `eta_noise(t)` の数式と係数。
16. `c_sprout`、萌芽持続step、方向連続性閾値。
17. 成立用のq比閾値、持続step、exec間相対時刻許容差、精度間再現許容差。
18. 局所回帰の窓幅、回帰次数、\(g>0\)有意性、\(R^2\)閾値、最小持続step。
19. 拡大停止の複合基準、各閾値、信頼区間構成法。
20. 準安定開始の有界性・振動・持続基準。
21. 第7論文とStage Aを正規化する「特徴時間」の定義。
22. `q5`〜`q8` の第7論文と同一な定義。既存 `Q=[B0|Bdom]` は最大4特異値しか持たない。
23. `direction1_continuity`〜`direction4_continuity` の同一定義。
24. 符号反転、順序交換、縮退2次元内回転を含む方向追跡規則。
25. `flux_parent_to_d3`, `flux_d3_to_parent`, `flux_d3_to_d4`, `net_flux_d3`, `net_flux_d4` の既存根拠ある定義。
26. 図6が要求する振動位相の列と計算方法。
27. 図8の時間平行移動の可否と整列規則。
28. 各runの固定 `max_step`。
29. 状態ベクトルの保存形式、各測定点で保存する状態の範囲、高精度値を丸めず保存する形式。
30. 予想容量、容量上限、残容量ゲート。
31. 予想実行時間、実行時間上限、「大幅超過」の倍率。
32. closure・norm・固有対残差・直交性・保存量の合否閾値。
33. exec 1/2で同一PRNG seedを用いる場合のbitwise一致要件と、異なる独立乱数を要求するか否か。
34. 高精度親の初期位相をfloat64 PRNG出力から与えることが許容されるか。
35. 中断試行の `mp_engine_v1.py`, `stageA_common_v1.py` を今後どの場所へ隔離するか。

未解決事項を「妥当と思われる値」で埋めない。

---

## 12. 最終確認

- [ ] 第7論文N=5再現法が確定
- [ ] 方向定義が第7論文と同一
- [ ] 高精度親構築法が確定
- [ ] 別固定点判定法が確定
- [ ] 全run条件が固定
- [ ] 全イベント定義が固定
- [ ] 実行上限が固定
- [ ] 図表仕様が固定
- [ ] 未定義事項がゼロ、または人間承認済み

### 実行許可状態

**不許可。**

本書提示後に停止する。人間が未解決事項を解消し、本書を明示承認するまで、Pythonコード作成・編集、実験、後処理、Stage A実行、Stage B/C準備を開始しない。
