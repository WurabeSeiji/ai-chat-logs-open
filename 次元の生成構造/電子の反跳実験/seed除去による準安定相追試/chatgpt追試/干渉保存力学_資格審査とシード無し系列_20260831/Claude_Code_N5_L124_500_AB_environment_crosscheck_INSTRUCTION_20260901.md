# Claude Code 指図書 --- N=5 / L=124 / 500 step / 64-bit 二系統クロスチェック

作成日: 2026-09-01

## 目的

ChatGPT 環境で、次の二つの独立な実装を **N=5, L=124, 500 step, 64-bit**
で実行したところ、0〜500 step の結果は実質完全一致した。

1.  過去の正本 `pass2_run.py` を忠実にコピーし、`STEPS=40000` だけを
    `STEPS=500` に変更した正本コピー系。
2.  ChatGPT が新しく作成した denominator-control 系。

Claude Code
の実行環境で、この二系統をそれぞれ実行し、両者が一致するかを検証すること。同時に
Claude Code 自身の Python / NumPy / BLAS / LAPACK / CPU / OS /
スレッド環境を記録し、ChatGPT 環境との差を報告すること。

**物理解釈は不要。まず数値再現性の監査だけを行うこと。**

------------------------------------------------------------------------

## 比較対象フォルダ A --- 正本コピー系

Google Drive 上のパス:

``` text
/Google Drive/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_N5_L124_500_exact_canonical_crosscheck_20260901/
```

主要ファイル:

``` text
program/pass2_run_original.py      # 過去の正本そのもの
program/pass2_run_N5_500.py        # 正本から STEPS=40000 -> 500 のみ変更
program/original_engine.py         # 正本が使用する engine
analysis/pass2_only_change.diff    # 上記変更が STEPS のみであることを示す diff
data/hm_N5/parent_v.npz            # 正本の保存済み初期値
data/hm_N5/states_treatment.npz    # ChatGPT環境で正本コピーを500 step実行した状態列
data/hm_N5/treatment_linear124_amplitude_aware_timeseries.csv
analysis/environment.json          # ChatGPT実行環境
analysis/sha256.json
ANALYSIS.md
SELF_INSTRUCTION.md
```

注意: `pass2_run_N5_500.py` は力学部分を変更していない。正本は
`key_steps` に 750 以上を固定列挙しているため、500 step
の状態・CSV・NPZ・summary 保存後、付随する `key_steps.csv` 作成時に
IndexError
になる。このエラーを避けるために力学コードを改変してはならない。比較対象となる
step 0〜500 の501状態は、そのエラーより前に生成済みである。

------------------------------------------------------------------------

## 比較対象フォルダ B --- ChatGPT 新規 denominator-control 系

Google Drive 上のパス:

``` text
/Google Drive/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_denominator_controls_64bit_with124_20260901/
```

主要ファイル:

``` text
run_and_plot.py                     # ChatGPTが新規作成した実行・図化コード
timeseries_64bit_with124.csv        # N=3..16, denominator=N-2..N+2,124 の全時系列
hm_N5_den_124_states_500.npz        # 今回比較する N=5, L=124 の501状態
summary_64bit_with124.csv
fig_Hperp_denominator_controls_with_124_64bit.png
RUN_METADATA.json
ANALYSIS.md
SHA256SUMS.txt
```

今回のクロスチェックでは、このフォルダのうち **N=5, denominator=124,
step=0..500 だけ**を使用すること。他の N や denominator
は実験対象にしない。

------------------------------------------------------------------------

## ChatGPT 側で確認済みの計算条件

両系統とも対象条件は:

``` text
N = 5
L = 124
STEPS = 500
state dtype = numpy.complex128
real dtype = numpy.float64
64-bit 浮動小数点系
固有値分解 = numpy.linalg.eigh
初期値 = 保存済み parent/state の complex128 配列
正規化追加なし
新規ノイズ追加なし
```

発展写像は正本では概略:

``` text
H = A * (conj(Z)[:,None] * Z[None,:])
w, V = np.linalg.eigh(H)
Z_next = V @ (exp(-1j*(2*pi/L)*w) * (V.conj().T @ Z))
```

比較量:

``` text
H_total
H_perp / H
Z(step) 全 complex128 状態
```

------------------------------------------------------------------------

## ChatGPT 実行環境 --- 必ず Claude Code 環境と比較すること

ChatGPT 側で記録された環境:

``` text
Python: 3.13.5
NumPy: 2.3.5
OS/platform: Linux 6.18.35 x86_64, glibc 2.41
state dtype: complex128
real dtype: float64
BLAS/LAPACK: scipy-openblas / OpenBLAS 0.3.30
OpenBLAS: USE64BITINT, DYNAMIC_ARCH, NO_AFFINITY
OpenBLAS target reported: Haswell
MAX_THREADS=64
compiler: GCC 14.2.x
CPU architecture: x86_64
SIMD available includes AVX2 and AVX-512 families
```

この情報の正本はフォルダ A の:

``` text
analysis/environment.json
```

に保存されている。

------------------------------------------------------------------------

# Claude Code が実施する作業

## 1. 実行前ファイル監査

最初に、上記二フォルダが存在し、主要ファイルを読めることを確認する。

特に A について:

``` bash
diff -u program/pass2_run_original.py program/pass2_run_N5_500.py
```

相当の比較を行い、**意図した差が `STEPS=40000 -> 500`
だけであることを確認**すること。

また SHA256 を取得し、報告書に記録すること。

## 2. 初期値監査

A の:

``` text
data/hm_N5/parent_v.npz
```

から正本初期値を読み、B が N=5/L=124 に使用する初期値と `np.array_equal`
で比較すること。

報告する値:

``` text
shape
dtype
np.array_equal
max(abs(delta Z0))
SHA256 of source files
```

初期値を再生成してはならない。保存済み binary64/complex128
配列をそのまま使用すること。

## 3. Claude Code 環境で系 A を実行

A の正本コピー系を N=5 / L=124 / 500 step で実行する。

力学式、`np.linalg.eigh`、初期値、dtype を変更しない。

step 0〜500 の全状態を保存する。

## 4. Claude Code 環境で系 B を実行

B の `run_and_plot.py` の同じ N=5 / denominator=124 / 500 step
条件を実行する。

他の N や denominator を走らせる必要はない。必要なら
**実験式を変えずに実行対象だけ N=5/L=124
に限定したラッパー**を作成し、そのラッパーも保存すること。

## 5. A と B を step 単位・状態配列単位で比較

step=0..500 全点について最低限:

``` text
max_abs_Z_difference
H_total_A
H_total_B
abs_delta_H_total
Hperp_frac_A
Hperp_frac_B
abs_delta_Hperp_frac
```

をCSVへ保存する。

次を明記すること:

``` text
最初に Z が bitwise 不一致になる step
最初に H_total が不一致になる step
最初に Hperp/H が不一致になる step
501状態が bitwise 完全一致するか
全501点の max(abs(delta Z))
全501点の max(abs(delta H_total))
全501点の max(abs(delta Hperp/H))
```

## 6. ChatGPT 保存結果とも比較

ChatGPT環境では、A と B は次の値になった:

``` text
step 0 Hperp/H = 4.755163423728637e-69
step 1 Hperp/H = 2.116260046380432e-31
step 500 Hperp/H ~= 2.792977647318284e-27
```

A と B の `H_total` は501点すべて bitwise 一致し、`Hperp/H` の最大差は約
`7.17e-43` だった。

Claude Code 環境で得た A/B の結果が、この ChatGPT
結果と一致するかも比較すること。

なお、過去に保存された元の正本 N=5/L=124 データでは step 1 が:

``` text
4.923807575912660e-31
```

であり、現在の ChatGPT 環境で正本を再実行した値 `2.116260046380432e-31`
と異なる。これが今回の環境クロスチェックの理由である。

## 7. Claude Code 自身の実行環境を完全に記録

最低限、以下を保存・報告すること:

``` text
python --version
Python executable path
platform.platform()
platform.machine()
NumPy version
np.__config__.show()
BLAS implementation and version
LAPACK implementation and version
OpenBLAS/MKL/Accelerate 等の識別
CPU architecture / model
OS version
thread environment variables
OMP_NUM_THREADS
OPENBLAS_NUM_THREADS
MKL_NUM_THREADS
VECLIB_MAXIMUM_THREADS
NUMEXPR_NUM_THREADS
```

macOS なら Apple Accelerate を使用しているかも明記すること。

可能なら `np.show_config()` / `np.__config__.show()`
の全文をテキスト保存する。

## 8. 図化

N=5, L=124, step 0..500 の `Hperp/H` を semilogy で重ねる。

最低3系列:

``` text
Claude Code run A: canonical-copy
Claude Code run B: ChatGPT-new-program
ChatGPT saved reference
```

必要なら過去保存正本を4本目として追加する。

図だけで一致判定せず、必ず数値比較CSVを主判定にする。

------------------------------------------------------------------------

# 保存先

Claude Code
の全成果物は、元の二フォルダを変更せず、新しいフォルダを作って保存すること。

推奨パス:

``` text
/Google Drive/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/ClaudeCode_N5_L124_500_AB_environment_crosscheck_20260901/
```

最低限保存するもの:

``` text
ANALYSIS.md
ENVIRONMENT.txt
ENVIRONMENT.json
NUMPY_CONFIG.txt
SHA256SUMS.txt
comparison_A_vs_B_steps0_500.csv
comparison_vs_chatgpt_reference.csv
A_states_500.npz
B_states_500.npz
A_timeseries_500.csv
B_timeseries_500.csv
figure_A_B_chatgpt_comparison.png
実行に使った全コード
stdout/stderr
```

------------------------------------------------------------------------

# 最終報告で必ず答える質問

1.  Claude Code 環境で A と B は step 0〜500 で一致したか。
2.  一致した場合、bitwise 完全一致か、丸め誤差レベルの一致か。
3.  Claude Code の step 1 `Hperp/H` は `4.923807575912660e-31` と
    `2.116260046380432e-31` のどちらに近いか、または別の値か。
4.  Claude Code の step 500 `Hperp/H` はいくつか。
5.  ChatGPT 環境との差はどの step から生じるか。
6.  Python / NumPy / BLAS / LAPACK / CPU / OS / thread 条件は ChatGPT
    環境とどう違うか。
7.  観測された差がコード差・入力差で説明できるか。説明できない場合のみ、数値ライブラリ/実行環境差を候補として報告すること。

**推測で原因を断定しない。まず A/B
の再実行結果と生データ比較を提示すること。**
