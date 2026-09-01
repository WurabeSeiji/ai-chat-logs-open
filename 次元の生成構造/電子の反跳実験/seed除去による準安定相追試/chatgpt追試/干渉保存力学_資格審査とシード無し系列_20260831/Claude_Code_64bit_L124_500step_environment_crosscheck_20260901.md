# Claude Code 向け：64-bit / L=124 / N=5 / 500 step 環境差・再現性対照試験 指示書

作成日: 2026-09-01

## 0. 目的

同一の初期データ・同一の64-bit数値型・同一の発展式・同一の L=124・同一の N=5 であるにもかかわらず、保存済み正本データと現在のChatGPT実行環境の再計算結果で step 1 以降が一致しない。

既知の値は以下。

- 保存済み正本データ（Claude Code が過去に生成）
  - step 0: H_perp/H = 4.7551634237e-69
  - step 1: H_perp/H = 4.9238075759e-31
- 現在のChatGPT実行環境で、正本コード相当を64-bitで再計算
  - step 0: H_perp/H = 4.7551634237e-69
  - step 1: H_perp/H = 2.1162600464e-31

step 0 は一致し、差は最初の更新 step 0 -> 1 で発生する。

本試験の目的は、Claude Code 自身の現在の実行環境で、

1. 正本プログラムそのもの
2. ChatGPT 側で今回作成した対照実験の L=124 系

を同じ N=5・同じ parent_v.npz・500 step で実行し、

- 保存済み正本データに一致するのか
- ChatGPT 現環境の再計算結果に一致するのか
- それとも第三の軌道になるのか

を厳密に判定する。

物理解釈は行わず、まず数値再現性だけを検証すること。

---

## 1. 絶対に変更してはいけない条件

以下を固定する。

- N = 5
- L = 124
- STEPS = 500
- 初期状態は保存済み正本 `data/hm_N5/parent_v.npz` の `v`
- `Z = v.copy()` から開始
- 数値型は NumPy `complex128` / `float64`
- 固有値分解は `np.linalg.eigh`
- 発展式は正本のまま

  H = A * (conj(z)[:,None] * z[None,:])

  w, V = np.linalg.eigh(H)

  Z_next = V @ (exp(-1j*(2*pi/L)*w) * (V.conj().T @ Z))

- H_perp は正本と同じ直接読出し

  p = Re(v) / ||Re(v)||

  q = Im(v) - (Im(v) @ p) p

  q = q / ||q||

  Z_perp = Z - p*(p@Z) - q*(q@Z)

  H_perp = vdot(Z_perp, Z_perp).real

  H_total = vdot(Z, Z).real

  H_perp/H = H_perp / H_total

- 正規化を追加しない
- seed を追加しない
- 初期値を再生成しない
- `parent_v.npz` を加工しない
- 高精度化しない
- mpmath を使わない
- dt, delta, L, N の意味を変更しない
- BLAS/LAPACKの実装を意図的に変更しない

---

## 2. 正本側の対象

親パス:

`次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831`

正本プログラム:

`program/pass2_run.py`

正本初期データ:

`data/hm_N5/parent_v.npz`

保存済み正本時系列:

`data/hm_N5/treatment_linear124_amplitude_aware_timeseries.csv`

保存済み正本状態:

`data/hm_N5/states_treatment.npz`

### 正本テストA

`pass2_run.py` は 40000 step 固定なので、正本ファイル自体は変更せず、コピーを作る。

コピー名例:

`audit_20260901/pass2_run_N5_L124_500_exactcopy.py`

変更してよいのはただ一箇所だけ:

`STEPS=40000` -> `STEPS=500`

L=124、式、関数、初期値読込、metrics は一切変更しない。

変更前後の unified diff を保存すること。

---

## 3. ChatGPT側で今回作成した対照系の対象

今回ChatGPT側で作成した 64-bit denominator control 実験は、

- N = 3..16
- denominator = N-2, N-1, N, N+1, N+2, 124
- 500 step
- complex128 / float64

で実行した。

そのうち比較対象は N=5, denominator=124 の系列だけ。

ChatGPT側の保存ファイルがGoogle Driveに存在する場合は、その実行コードと時系列を取得して使用する。

候補名:

- `timeseries_64bit_with124.csv`
- `summary_64bit_with124.csv`
- `fig_Hperp_denominator_controls_with_124_64bit.png`

もし実行コードが別名で保存されている場合は、必ず実ファイルを特定してから使用すること。

推測で書き直してはいけない。

### 対照テストB

ChatGPT側の実行コードをそのまま Claude Code 環境で実行し、N=5, denominator=124 のみを抽出する。

正本テストAと同じ `parent_v.npz` を使用していることを bitwise 検証する。

---

## 4. 実行前監査

A/B両方について以下を保存する。

### 4.1 ファイルハッシュ

SHA256:

- 正本 `pass2_run.py`
- 正本コピー `pass2_run_N5_L124_500_exactcopy.py`
- ChatGPT対照コード
- `parent_v.npz`
- 保存済み正本 `states_treatment.npz`
- 保存済み正本 CSV

`audit_hashes.json` に保存。

### 4.2 Python / NumPy / BLAS / LAPACK / CPU 環境

最低限以下を記録する。

- Python version
- NumPy version
- `np.__config__.show()` の全文
- OS
- machine architecture
- processor / CPU
- `OMP_NUM_THREADS`
- `OPENBLAS_NUM_THREADS`
- `MKL_NUM_THREADS`
- `VECLIB_MAXIMUM_THREADS`
- `NUMEXPR_NUM_THREADS`

可能なら以下も記録する。

- BLAS vendor
- LAPACK vendor
- linked shared libraries
- Accelerate / OpenBLAS / MKL の別

`environment_report.txt` に保存。

---

## 5. 初期データの一致確認

以下を確認する。

1. `parent_v.npz['v']`
2. 保存済み正本 `states_treatment.npz['Z'][0]`
3. テストA実行開始時 Z0
4. テストB実行開始時 Z0

全てについて

- dtype
- shape
- `np.array_equal`
- max abs difference
- raw bytes SHA256

を保存する。

期待値:

`np.array_equal == True`

差がある場合は実行を中止し、その時点で報告する。

---

## 6. step 0 -> 1 の中間値監査

ここが最重要。

A/Bそれぞれで step 0 の状態 Z0 から次を保存する。

1. adjacency matrix A
2. H(Z0)
3. `w, V = np.linalg.eigh(H)` の w
4. V
5. `V.conj().T @ Z0`
6. `phase = exp(-1j*(2*pi/124)*w)`
7. `phase * (V.conj().T @ Z0)`
8. Z1
9. H_total(step0)
10. H_perp(step0)
11. H_perp/H(step0)
12. H_total(step1)
13. H_perp(step1)
14. H_perp/H(step1)

各配列について

- dtype
- SHA256(raw bytes)
- min/max
- norm

を記録する。

AとBを同じClaude Code環境内で比較し、最初に一致しなくなる中間量を特定する。

---

## 7. 保存済み正本との比較

保存済み正本データについて、step 0..500 を読み込む。

比較対象:

- test A vs saved original
- test B vs saved original
- test A vs test B

最低限以下を比較する。

- H_parallel
- H_perp
- H_total
- H_perp/H
- PR
- PR_over_M
- abs_ZT_Z
- amp_min
- amp_max
- amp_std

状態NPZがあるものについては Z 自体を比較する。

各 step で

`max_abs_state_difference`

を計算する。

最初の不一致stepを特定する。

---

## 8. 既知の基準値

保存済み正本 N=5, L=124:

- step 0 H_perp/H = 約 `4.7551634237e-69`
- step 1 H_perp/H = 約 `4.9238075759e-31`

ChatGPT現在環境の再計算:

- step 0 H_perp/H = 約 `4.7551634237e-69`
- step 1 H_perp/H = 約 `2.1162600464e-31`

Claude Code環境で、A/Bがどちらに一致するかを判定する。

判定例:

- A = original, B = original
  -> Claude環境では正本値を再現。ChatGPT現環境との差が有力。

- A = current-ChatGPT-value, B = current-ChatGPT-value
  -> 現在のClaude環境もChatGPT現環境と同系統。過去正本生成時の環境差が有力。

- A != B
  -> プログラム差または入力経路差が残っている。step0->1中間監査で特定。

- A = B だが originalともcurrentとも違う
  -> 数値線形代数環境依存がさらに強く疑われる。

---

## 9. 500 step 時系列比較

0..500 の全501点について、H_perp/H を一つのCSVへまとめる。

必須列:

- step
- saved_original_Hperp_frac
- testA_exactcopy_Hperp_frac
- testB_chatgpt_control_Hperp_frac
- A_minus_original
- B_minus_original
- A_minus_B

状態比較可能なら:

- max_abs_Z_A_minus_original
- max_abs_Z_B_minus_original
- max_abs_Z_A_minus_B

も保存。

ファイル名:

`comparison_step0_500.csv`

---

## 10. 図化

同一図に3系列を描く。

1. 保存済み正本
2. Claude Code test A
3. Claude Code test B

縦軸:

`H_perp/H` 対数

横軸:

step 0..500

縦軸範囲は最低でも `1e-70` から `1e-24` を表示し、step 0->1 の38桁ジャンプが見えるようにする。

出力:

`N5_L124_500_saved_vs_testA_vs_testB.png`

さらに step 0..10 拡大図:

`N5_L124_step0_10_detail.png`

---

## 11. 出力先

親パス:

`次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831`

新規フォルダ:

`audit_20260901_claude_vs_chatgpt_L124_N5`

構成:

```text
audit_20260901_claude_vs_chatgpt_L124_N5/
  README.md
  environment_report.txt
  audit_hashes.json
  program/
    pass2_run_N5_L124_500_exactcopy.py
    chatgpt_control_code.py
    diff_pass2_exactcopy.patch
  raw/
    intermediate_step0_to_1_testA.npz
    intermediate_step0_to_1_testB.npz
    states_testA.npz
    states_testB.npz
  data/
    timeseries_testA.csv
    timeseries_testB.csv
    comparison_step0_500.csv
    initial_state_audit.json
    step1_intermediate_comparison.json
  figures/
    N5_L124_500_saved_vs_testA_vs_testB.png
    N5_L124_step0_10_detail.png
  analysis/
    RESULT.md
```

---

## 12. RESULT.md に必ず書くこと

推測ではなく、実測だけを書く。

1. 正本初期値とA/B初期値が bitwise 同一か
2. A/Bのプログラム差
3. A/Bの step 1 が一致するか
4. 保存済み正本の step 1 と一致するのはA/Bのどちらか
5. 最初に差が発生した中間量
6. Python / NumPy / BLAS / LAPACK / CPU環境
7. step 0..500 の最大差
8. 結論を以下のいずれかで分類

- CODE_DIFFERENCE
- INPUT_DIFFERENCE
- LIBRARY_ENVIRONMENT_DIFFERENCE
- PAST_ENVIRONMENT_NOT_REPRODUCIBLE
- UNRESOLVED

原因が特定できない場合は `UNRESOLVED` とする。推測で環境差と断定しない。

---

## 13. 最重要ルール

- 一度に複数条件を変えない。
- 高精度化しない。
- 64-bit以外を混ぜない。
- Nを変えない。
- Lを変えない。
- seedを入れない。
- parentを再生成しない。
- 500 step以外へ変更しない。
- 保存済み正本データを上書きしない。
- 既存正本コードを直接編集しない。
- 結果が予想と違っても条件を変更しない。
- まず生データ差を報告し、その後に解釈する。

この試験では「物理の解釈」ではなく、「同一64-bit計算が異なる環境でどこから分岐するか」の監査を行う。
