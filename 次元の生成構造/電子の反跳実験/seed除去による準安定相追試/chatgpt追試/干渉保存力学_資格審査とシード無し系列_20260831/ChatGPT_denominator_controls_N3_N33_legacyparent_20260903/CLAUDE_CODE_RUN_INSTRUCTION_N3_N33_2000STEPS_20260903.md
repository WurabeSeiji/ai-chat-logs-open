# Claude Code 実行指示 — N=3..33 legacy parent / 2000 step 延長実験

## 目的

既に正常完了した 500 step 実験

`run_and_plot_N3_N33_legacyparent_20260903.py`

を基準に、**初期値・力学・分母系列・N範囲を一切変更せず、時間ステップだけ 500 → 2000 に延長**して、N が大きい系列でもインフレーションの立上がりから飽和まで観測できるか確認する。

既存の 500 step プログラムと `results/` は正本として保存し、**絶対に編集・上書きしないこと**。

---

## 基準プログラム

ファイル:

`run_and_plot_N3_N33_legacyparent_20260903.py`

期待 SHA256:

`dd2b9f7d740bfe5a6013e575c59234324c09e59708a8ca8e0b4050f2ae9ac342`

最初に SHA256 を確認し、一致しない場合は何も変更・実行せず停止して報告すること。

---

## 作成する実行プログラム

基準プログラムをコピーして、次の名前にする。

`run_and_plot_N3_N33_legacyparent_2000steps_20260903.py`

**基準プログラムそのものは編集禁止。**

コピー後のプログラムに対して、以下の変更だけを行う。

### 許可する変更は以下だけ

1. 出力先を既存 `results` から新規 `results_2000steps` に変更する。

```python
OUT=os.path.join(ROOT,'ChatGPT_denominator_controls_N3_N33_legacyparent_20260903','results_2000steps')
```

2. ステップ数を 500 から 2000 に変更する。

```python
STEPS=2000; OFFSETS=(-2,-1,0,1,2)
```

3. 状態ファイル名を

```text
states_500.npz
```

から

```text
states_2000.npz
```

に変更する。

4. CSV出力名を次に変更する。

```text
timeseries_64bit_with124_N3_N33_2000steps.csv
summary_64bit_with124_N3_N33_2000steps.csv
```

5. 図の x 軸上限を 500 から 2000 に変更する。

```python
ax.set_xlim(0,2000)
```

6. 図ファイル名を次に変更する。

```text
fig_Hperp_denominator_controls_with_124_N3_N33_2000steps.png
```

7. メタデータファイル名を次に変更する。

```text
RUN_METADATA_N3_N33_2000steps.json
```

**これ以外の変更は禁止する。**

---

## 絶対禁止

- N範囲 `range(3,34)` を変更しない。
- N=3..16 の入力を変更しない。
- N=17..33 の legacy parent 入力を変更しない。
- `hm_mp_free_N3_N40_20260901` を使用しない。
- 初期値を再生成しない。
- 正規化・再スケールを追加しない。
- `edges` / `adjacency` / `H_of` / `one_step` / `plane` / `metrics` を変更しない。
- `np.linalg.eigh` を変更しない。
- dtype、演算順序、BLAS/LAPACK設定を変更しない。
- 分母系列 `(N-2,N-1,N,N+1,N+2,124)` を変更しない。
- 並列化、高速化、リファクタリングをしない。
- 警告を消すための修正をしない。
- パッケージのインストール・更新をしない。
- 既存 `results/` のファイルを削除・変更・上書きしない。
- Git commit / push はこの実験指示の範囲外なので行わない。

エラーが発生した場合、自分で修正せず、その場で停止してエラーをそのまま報告すること。

N=8 で前回と同じ matmul RuntimeWarning が出ても、**コード変更は行わず警告を記録するだけ**とする。Python exception または非0 exit code が発生した場合は停止する。

---

## 入力データ

N=3..16:

```text
.../干渉保存力学_資格審査とシード無し系列_20260831/data/hm_N*/states_treatment.npz
```

の `Z[0]`。

N=17..33:

```text
.../干渉保存力学_資格審査とシード無し系列_20260831/data/hm_N*/parent_v.npz
```

の `v`。

これは 500 step 実験と完全に同じ入力である。

---

## 実行手順

### 1. 同じフォルダへ移動

```bash
cd '/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_denominator_controls_N3_N33_legacyparent_20260903'
```

### 2. 基準プログラムのSHA256を確認

```bash
shasum -a 256 run_and_plot_N3_N33_legacyparent_20260903.py
```

上記の期待SHA256と一致した場合だけ続行する。

### 3. 実行環境を変更せず記録

```bash
{
  pwd
  which python3
  python3 -c "import sys,platform,numpy as np; print('python_executable=',sys.executable); print('python=',platform.python_version()); print('numpy=',np.__version__); print('platform=',platform.platform()); np.show_config()"
} | tee EXECUTION_ENVIRONMENT_N3_N33_2000steps_20260903.txt
```

前回と同じ `.venv/bin/python3` をそのまま使う。環境変更は禁止。

### 4. 入力ファイル存在確認

```bash
python3 - <<'PY'
from pathlib import Path
root=Path('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831')
missing=[]
for n in range(3,17):
    p=root/'data'/f'hm_N{n}'/'states_treatment.npz'
    if not p.is_file(): missing.append(str(p))
for n in range(17,34):
    p=root/'data'/f'hm_N{n}'/'parent_v.npz'
    if not p.is_file(): missing.append(str(p))
print('missing=',len(missing))
for p in missing: print(p)
PY
```

`missing=0` の場合だけ続行する。

### 5. 新規出力フォルダを作成

```bash
mkdir -p results_2000steps
```

既存 `results/` は触らない。

### 6. 基準プログラムをコピーし、上記7項目だけ変更

コピー先:

```text
run_and_plot_N3_N33_legacyparent_2000steps_20260903.py
```

変更後、基準プログラムとの差分を表示して確認する。

```bash
diff -u run_and_plot_N3_N33_legacyparent_20260903.py run_and_plot_N3_N33_legacyparent_2000steps_20260903.py
```

**上記「許可する変更」以外の差分が1箇所でもあれば実行せず停止して報告すること。**

### 7. 構文確認

```bash
python3 -m py_compile run_and_plot_N3_N33_legacyparent_2000steps_20260903.py
```

構文エラーがあれば修正せず停止して報告する。

### 8. 実行

```bash
python3 run_and_plot_N3_N33_legacyparent_2000steps_20260903.py
```

途中経過は `done N ...` をそのまま監視する。

---

## 実行後の確認

`ALL DONE` と exit code 0 を確認する。

次のファイルが存在することを確認する。

```text
results_2000steps/timeseries_64bit_with124_N3_N33_2000steps.csv
results_2000steps/summary_64bit_with124_N3_N33_2000steps.csv
results_2000steps/fig_Hperp_denominator_controls_with_124_N3_N33_2000steps.png
results_2000steps/RUN_METADATA_N3_N33_2000steps.json
```

状態ファイル

```text
results_2000steps/hm_N*_den_*_states_2000.npz
```

が **186個（31 N × 6分母）** あることを確認する。

最後に次を表示する。

```bash
shasum -a 256 \
  run_and_plot_N3_N33_legacyparent_20260903.py \
  run_and_plot_N3_N33_legacyparent_2000steps_20260903.py \
  CLAUDE_CODE_RUN_INSTRUCTION_N3_N33_2000STEPS_20260903.md
```

報告するのは、

- `ALL DONE` / exit code
- 完了したN範囲
- 警告・エラーの有無
- 出力4ファイルの存在とサイズ
- 状態NPZの個数
- 上記SHA256

だけとする。

**追加解析・コード修正・再実行・Git操作は行わない。**
