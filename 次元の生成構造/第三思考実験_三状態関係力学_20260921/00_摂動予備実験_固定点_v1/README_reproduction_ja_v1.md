# E3-A 摂動予備実験 固定点 v1 — 再現手順

日付: 2026-09-21  
位置づけ: 第三思考実験へ進む前に行った **階層化二体 Kepler 系の弱結合候補の診断実験**。  
重要: 本パッケージは **導出された三状態則ではない**。`(a,b)` と `((ab),c)` という階層を先に与えた候補が何をするかを調べた固定点である。

## 1. 固定する問い

二つの二自由度 Kepler セクターを相互線形結合した候補

```text
Xin[k+1] + Xin[k-1] = tau_in  Xin[k] + epsilon Xout[k]
Xout[k+1] + Xout[k-1] = tau_out Xout[k] + epsilon Xin[k]
```

について、次を確認した。

- 全体面積保存量が維持されるか。
- 生の階層座標では局所保存量がどの次数で交換されるか。
- 正常モードへ変換したとき、各モードが二自由度則へ戻るか。
- 非縮退階層系で外部法則パラメータのシフトが `O(epsilon^2)` になるか。
- 二次摂動近似・楕円型外部モード・内部/外部ラベルの分離が、結合増大に対してどの順番で失われるか。

この候補については、追加の広域走査を行わず、本パッケージを固定点として第三思考実験本体の導出へ戻る。

## 2. 実行環境

固定時の環境:

```text
Python 3.13.5
NumPy 2.3.5
Matplotlib 3.10.8
```

依存パッケージは `requirements_E3A_fixed_point_v1.txt` に固定した。

## 3. 再現順序

同じフォルダで以下を実行する。

```bash
python3 run_E3A0_hierarchical_two_kepler_perturbation_v1.py
python3 analyze_E3A1_perturbation_breakdown_v1.py
python3 plot_E3A_fixed_point_v1.py
python3 verify_E3A_fixed_point_v1.py
```

`plot_E3A_fixed_point_v1.py` は数値実験を再実行せず、保存済み CSV だけを入力として図を作る。

## 4. 固定データ

### E3-A0 数値実験

- `E3A0_hierarchical_two_kepler_perturbation_full_v1.json`  
  位相ごとの全走行を含む完全データ。
- `E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv`  
  図化・解析用の集約表。
- `E3A0_run_stdout_v1.txt`  
  固定時の標準出力。

主系列は `n_in=31`, `n_out=127`、8相対位相、最大 20,000 step。比較として縮退系列 `31,31` も保存している。

### E3-A1 解析

- `E3A1_perturbation_breakdown_scan_v1.csv`
- `E3A1_perturbation_breakdown_analysis_v1.md`
- `E3A1_analysis_stdout_v1.txt`

これは E3-A0 と同じ線形候補の固有値を解析したもので、新しい力や読み出しを追加していない。

## 5. 固定時に確認した主要結果

非縮退系列 `31,127` では、全体面積保存の最大相対ドリフトは数値精度域に留まる。生の内部/外部局所量は `epsilon` にほぼ一次で応答する一方、正常モードの面積保存と円錐曲線残差は数値精度域に留まる。

外部優勢正常モードの法則パラメータは

```text
delta tau_out = epsilon^2 / (tau_out - tau_in) + O(epsilon^4)
```

で始まる。固定パラメータでは二次近似誤差 5% が `epsilon ≈ 0.00881984`、外部正常モードが `lambda_out=2` に達する楕円/放物境界が `epsilon ≈ 0.01000936`、`2 epsilon / delta_tau = 1` の混合クロスオーバーが `epsilon ≈ 0.01924647` である。

この結果から固定する結論は、**この候補は保存量交換を表せるが、正常モードで二つの二自由度則へ厳密分解できるため、三状態既約力学の導出そのものではない**、である。

## 6. 図

すべて保存済み CSV から生成する。

- `fig_E3A_01_raw_hierarchical_response_v1.png`: 生の内部/外部局所面積変調と外部円錐残差。
- `fig_E3A_02_invariants_and_normal_modes_v1.png`: 全体保存、正常モード保存、正常モード円錐残差。
- `fig_E3A_03_outer_law_shift_v1.png`: 厳密な外部法則シフトと二次摂動近似。
- `fig_E3A_04_approximation_error_and_mixing_v1.png`: 二次近似誤差と階層モード混合。
- `fig_E3A_05_outer_mode_type_boundary_v1.png`: 外部正常モード固有値と楕円/放物境界。

## 7. 固定点への移送確認

固定点作成時に、旧作業名 `hierarchical_kepler_perturbation_v2.csv` と新しい E3-A0 CSV が byte-for-byte 一致すること、旧 `hierarchical_kepler_breakdown_v3.csv` と新しい E3-A1 CSV が byte-for-byte 一致することを確認した。E3-A0 JSON は固定点用の `experiment` ラベルだけを変更し、それ以外のデータ内容が一致することを確認した。

## 8. 完全性確認

`SHA256SUMS.txt` が固定ファイルの SHA-256 を保持する。`verify_E3A_fixed_point_v1.py` は一時フォルダで数値実験・解析・図化を再実行し、固定されたデータおよび図の SHA-256 と比較する。
