# N=5, L=124, 500 step 正本忠実コピー対照実験 分析

作成日: 2026-09-01

## 1. 実験目的
正本 `pass2_run.py` を忠実にコピーし、変更を `STEPS=40000 -> 500` の一箇所だけに限定してN=5を実行した。既存ChatGPTテスト系のN=5,L=124,500-stepデータ、および当時保存された正本40000-stepデータの先頭501行と比較した。

## 2. 正本出所
Google Drive親パス:
`/Google Drive/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831`

使用正本:
- `program/pass2_run.py`
- `program/original_engine.py`
- `data/hm_N5/parent_v.npz`
- `data/hm_N5/treatment_linear124_amplitude_aware_timeseries.csv`（当時保存正本）

## 3. 正本コピーへの変更
`analysis/pass2_only_change.diff` の通り、変更は次の一箇所のみ。

`STEPS=40000` -> `STEPS=500`

N=5, L=124、H生成、`np.linalg.eigh`、更新式、metrics、初期値読込は未変更。

注意: 正本はkey_stepsに750以上を含むため、500 step化したコピーは力学走行と主要出力保存完了後、key_steps書出しでIndexError停止した。`treatment_linear124_amplitude_aware_timeseries.csv` はheader+501行、`states_treatment.npz` は501状態で生成完了している。力学コードを余分に変更しないため、この後処理エラーは修正せず、比較・図化を別処理で行った。

## 4. 比較対象
### A. 今回の正本忠実コピー500-step
ローカル生成元: `exact_reproduction/data/hm_N5/treatment_linear124_amplitude_aware_timeseries.csv`

### B. 既存ChatGPTテスト系 N=5, L=124, 500 step
元データ: `timeseries_den124_500.csv` の `N=5, denominator=124` 501行。
今回の保存パッケージには `data/chatgpt_test_N5_L124_500.csv` として切り出しを保存。

### C. 当時保存された正本
Google Drive: `.../data/hm_N5/treatment_linear124_amplitude_aware_timeseries.csv`
先頭501行を `data/original_saved_N5_L124_first501.csv` として保存。

## 5. 結果
正本忠実コピーAと既存ChatGPTテストBは:
- H_total 最大絶対差 = `0.00000000000000000e+00`
- H_total 全501点 bitwise一致 = `True`
- H_perp/H 最大絶対差 = `7.17464813734306340e-43`

H_perp/HはCSVへの演算・書出し桁の違いにより最大 `7.175e-43` の極微小差があるが、step 0,1,2など主要値は一致し、数値軌道は同一。例:
- step0: A=4.75516342372863681e-69, B=4.75516342372863681e-69
- step1: A=2.11626004638043211e-31, B=2.11626004638043211e-31
- step500: A=2.79297764731828427e-27, B=2.79297764731828462e-27

一方、当時保存正本Cとはstep1から異なる:
- step0 C=4.75516342372863681e-69, A=4.75516342372863681e-69
- step1 C=4.92380757591266043e-31, A=2.11626004638043211e-31
- 0..500での H_perp/H 最大絶対差 = `3.87177296449173966e-27`

## 6. 結論
今回の現在環境では、正本 `pass2_run.py` の忠実コピー（STEPSだけ500に変更）と既存ChatGPTテスト系は、N=5,L=124の500-step軌道を実質同一に再現した。したがって、既存ChatGPTテスト系のL=124線が別の力学式で生成された、という証拠は得られなかった。

当時保存された正本データとの差はstep1から存在する。今回の比較だけから原因を断定しないが、少なくとも「今回正本忠実コピー vs 既存ChatGPTテスト系」の差ではない。

## 7. Google Drive保存先
`/Google Drive/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_N5_L124_500_exact_canonical_crosscheck_20260901/`

保存内容:
- `SELF_INSTRUCTION.md`
- `program/pass2_run_original.py`
- `program/pass2_run_N5_500.py`
- `program/original_engine.py`
- `data/hm_N5/parent_v.npz`
- `data/hm_N5/treatment_linear124_amplitude_aware_timeseries.csv`
- `data/hm_N5/states_treatment.npz`
- `data/chatgpt_test_N5_L124_500.csv`
- `data/original_saved_N5_L124_first501.csv`
- `analysis/pass2_only_change.diff`
- `analysis/comparison_N5_L124_steps0_500.csv`
- `analysis/environment.json`
- `analysis/sha256.json`
- `ANALYSIS.md`
- `figures/N5_L124_500_canonical_vs_chatgpt_vs_original_saved.png`
