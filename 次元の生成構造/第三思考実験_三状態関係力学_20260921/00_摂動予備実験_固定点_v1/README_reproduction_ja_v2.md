# E3-A 摂動予備実験 固定点 v2 — 軌道図を含む完全再現手順

日付: 2026-09-21  
位置づけ: 第三思考実験へ進む前に行った **階層化二体 Kepler 系の弱結合候補の診断実験**。  
重要: 本パッケージは **導出された三状態則ではない**。`(a,b)` と `((ab),c)` という階層を先に与えた候補が何をするかを調べた固定点である。

## 1. 固定する候補

```text
Xin[k+1] + Xin[k-1] = tau_in  Xin[k] + epsilon Xout[k]
Xout[k+1] + Xout[k-1] = tau_out Xout[k] + epsilon Xin[k]
```

E3-A0 はこの候補の数値走査、E3-A1 は同じ線形候補の固有値・摂動近似の解析、E3-A2 は同一候補・同一初期条件から軌道座標を決定論的に再生成して固定 CSV に保存したものである。

## 2. 実行環境

固定時:

```text
Python 3.13.5
NumPy 2.3.5
Matplotlib 3.10.8
```

依存は `requirements_E3A_fixed_point_v1.txt` に固定している。

## 3. 完全再現順序

同じフォルダで次を実行する。

```bash
python3 run_E3A0_hierarchical_two_kepler_perturbation_v1.py
python3 analyze_E3A1_perturbation_breakdown_v1.py
python3 plot_E3A_fixed_point_v1.py
python3 export_E3A2_orbit_trajectories_v1.py
python3 plot_E3A_orbits_v1.py
python3 verify_E3A_fixed_point_v2.py
```

`plot_E3A_fixed_point_v1.py` は E3-A0/A1 の保存済み CSV だけを読み、図1〜5を作る。  
`plot_E3A_orbits_v1.py` は保存済み `E3A2_orbit_trajectories_v1.csv` だけを読み、軌道図6〜7を作る。

## 4. 固定データ

### E3-A0 数値実験

- `E3A0_hierarchical_two_kepler_perturbation_full_v1.json`: 152 run の各診断量を含む完全集計。
- `E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv`: 位相集約済み解析表。
- `E3A0_run_stdout_v1.txt`: 固定時標準出力。

主系列: `n_in=31`, `n_out=127`, 8相対位相, 最大20,000 step。比較として縮退系列 `31,31` も保存。

### E3-A1 解析

- `E3A1_perturbation_breakdown_scan_v1.csv`
- `E3A1_perturbation_breakdown_analysis_v1.md`
- `E3A1_analysis_stdout_v1.txt`

### E3-A2 軌道データ

- `E3A2_orbit_trajectories_v1.csv`

軌道図用に、主系列 `n_in=31`, `n_out=127`, `u_in=0.37`, `u_out=1.11`、
`epsilon = 0, 0.001, 0.004, 0.008` を固定し、各 epsilon につき256点を保存する。
CSV には raw outer readout と outer-dominant normal-mode readout の `(x,y)` を同時保存する。
軌道図はこの CSV だけから再描画できる。

## 5. 固定された主要結果

非縮退系列では全体面積保存のドリフトは数値精度域。生の内部/外部局所量は epsilon にほぼ一次で応答する一方、正常モードの面積保存と円錐曲線残差は数値精度域に留まる。

外部優勢正常モードの法則パラメータは

```text
delta tau_out = epsilon^2 / (tau_out - tau_in) + O(epsilon^4)
```

で始まる。固定パラメータでは二次近似誤差5%が `epsilon ≈ 0.00881984`、外部正常モードの楕円/放物境界が `epsilon ≈ 0.01000936`、混合クロスオーバー `2 epsilon / delta_tau = 1` が `epsilon ≈ 0.01924647`。

固定する結論は、**この候補は保存量交換を表せるが、正常モードで二つの二自由度則へ厳密分解できるため、三状態既約力学の導出そのものではない**、である。

## 6. 図化

保存済みデータから次の7図を再生成できる。

- `fig_E3A_01_raw_hierarchical_response_v1.png`: 生の局所面積変調と外部円錐残差。
- `fig_E3A_02_invariants_and_normal_modes_v1.png`: 全体保存、正常モード保存、正常モード円錐残差。
- `fig_E3A_03_outer_law_shift_v1.png`: 外部法則シフトと二次摂動近似。
- `fig_E3A_04_approximation_error_and_mixing_v1.png`: 二次近似誤差と階層混合。
- `fig_E3A_05_outer_mode_type_boundary_v1.png`: 外部正常モードの型境界。
- `fig_E3A_06_orbit_raw_outer_readout_v1.png`: **生の外部読み出し軌道**。
- `fig_E3A_07_orbit_outer_normal_mode_v1.png`: **外部優勢正常モードの軌道**。

図6では epsilon 増大に伴って生の階層座標が単一円錐から大きく外れる様子を直接確認できる。図7では同じ状態を正常モードで読むと各 epsilon で単一円錐軌道が保たれることを確認できる。

## 7. 完全性・再現性

- `SHA256SUMS_v2.txt`: v2 固定点の全主要ファイルの SHA-256。
- `verify_E3A_fixed_point_v2.py`: 一時ディレクトリで E3-A0/A1/A2 と全7図を再生成し、固定版と SHA-256 比較する。
- `verification_v2.txt`: 固定時の照合結果。

追加の広域走査は行わず、この一式を摂動候補の固定点として第三思考実験本体の導出へ戻る。
