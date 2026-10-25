# 対照実験：N=3..40 段1+2+3 スイープの忠実コピー再走行（2026-10-25）

目的：原本パッケージ `電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/`（2026-09-05、Zenodo 22317636）を、プログラム無変更（パス 4 箇所のみ）・同一の Python/numpy で再走行し、原本の生成物と照合する。以後の「分配率 R の読出し」予備実験は、この対照が通ったコピーを起点に最小変更で行う。

## 構成（原本からの差分）

| ファイル | 原本との関係 |
|---|---|
| `run_n_scaling_lowrank_v1.py` | bit 同一（7月正本エンジンのコピー） |
| `make_static_parents_N3_N40_v1.py` | bit 同一（N=40 ゲートの参照 `REF40` は原本の正本静的親を指したまま） |
| `run_N3_N40_stage123_v1.py` | **パスのみ変更**：`ROOT`→本フォルダ、`OUT`→`results/`、`PARENT_DIR`→`parents/`、`NEW_IN`（未使用）は ROOT 経由。差分は `diff_run_N3_N40_stage123_v1_vs_original.txt`。力学 `one_step`（21〜28 行）は無変更 |
| `check_sweep_inputs_v1.py`, `analyze_sweep_summary_v1.py`, `plot_complex_plane_N3_N40_stage123_v1.py` | bit 同一 |
| `compare_with_original_v1.py` | **新規（照合ツール、読み出しのみ）**：親 npz・状態 npz の bit 一致と数値差、CSV/JSON/図の SHA256 |
| `run_all.sh` | 新規：venv の python で 親生成→スイープ→入力ゲート→集計→図化→照合 |
| `original_reference/` | 原本の SHA256SUMS と RUN_METADATA の写し |
| `requirements_frozen.txt` | venv の pip freeze |

## 環境

- venv：`~/.venvs/stage123_control_py39`（`/usr/bin/python3` 3.9.6 で作成、`pip install numpy==2.0.2 matplotlib==3.9.4`）
- numpy 2.0.2 の BLAS/LAPACK：**Accelerate**（原本の記録と同じ）
- macOS：**26.6.2**（原本の記録は 26.3.1）。Darwin 25.6.0、arm64

## 進行記録

### 2026-10-25 親生成（make_static_parents_N3_N40_v1.py）→ GATE FAIL

- 38 親すべて収束（residual 3.5e-14〜9.4e-13、原本と同じ範囲）。
- **N=40 の v/g/Z0 は原本の正本静的親と bit 一致せず（GATE FAIL）**。38 親とも bit 一致 0/38（`results/compare_with_original_v1.json`）。
- 差の大きさ（原本との max|Δ|）：

| N | max\|Δv\| | max\|Δg\| | max\|ΔZ0\| | max\|Δσ\| | 備考 |
|---|---|---|---|---|---|
| 3 | 1.22 | 0.71 | 1.22 | 5.7e-09 | 別の解（同値性は未確認） |
| 4 | 1.00 | 0.98 | 1.00 | 平面数 3 対 2 | 別の解（σ スペクトルの本数が違う） |
| 5 | 0.56 | 4.0e-13 | 0.56 | 2.3e-09 | v だけ大きく違い g は一致。大域位相の違いの可能性（未確認） |
| 6〜40 | 2e-16〜3e-15 | 1e-16〜2e-15 | 2e-16〜3e-15 | 3e-15〜6e-14 | 末尾ビットの差 |

- 読み：プログラム・Python・numpy・BLAS は原本と同じなので、差の出所は macOS 26.3.1→26.6.2 の Accelerate（LAPACK の `eig`）の実装差と考えられる（未確認）。N≥6 の差 1e-16〜1e-15 は**種の振幅 δ=1e-15 と同じ桁**であり、スイープの出発点 H⊥/H(0) と onset はこの差を受ける。N=3,4,5 は反復が別の固定点に落ちた。
- スイープ本体は、この結果を木原に報告したうえで起動の可否を判断する（未起動）。

### 2026-10-25 案B 走行完了（`run_all_subset.sh`、所要 292 s）

木原判断（10/25）：原本との bit 一致は環境差で不可能と報告したうえで、**案 2＝今生成した親から走らせ原本と統計量で比べる**、かつ **5 分以内に収まる縮約（案 B）** で実行。

- 構成：N ∈ {3,4,5,6,8,10,12,15,20,25,30,40} × den ∈ {N, 124} × 500 step ＝ 24 走行。走行本体は `run_subset_stage123_v1.py`（`run_N3_N40_stage123_v1.py` から N リスト・分母系列・経過時間の印字・図題のみ変更、差分 `diff_run_subset_stage123_v1_vs_copy.txt`。`one_step`/`H_of` は diff で無変更を確認）。
- 所要：N=40 が 1 走行 109 s、N=30 が 21 s、N=25 が 8 s、N≤20 は合計 7 s 未満。合計 292 s。
- 入力ゲート PASS（24 状態ファイルの Z[0] が各親の Z0 と bit 一致）。24 走行中 21 が 500 step 内に H⊥/H＞0.05 を交差（未交差は N=3,30,40 の den=124。原本も同じ 3 走行が未交差）。
- 原本との状態 bit 一致：0/24（親が bit 一致しないため当然）。保存量：H_total の相対ドリフト最大 3.6e-13、閉塞 |zᵀz|/H の最大 5.9e-14（原本と同じ桁）。

**原本との曲線照合**（`results/compare_curves_table.md`、`results/compare_curves_with_original_v1.json`、図 `results/fig_compare_curves_control_vs_original.png`）：

| 項目 | 結果 |
|---|---|
| onset（H⊥/H＞0.05 の step） | 24 走行中 20 で原本と同一。差があるのは N=10 den=10（145/144）、N=15 den=15（111/110）、N=30 den=30（288/289）、N=12 den=124（195/189） |
| 指数区間の傾き（桁/step） | 全走行で原本と 3〜4 桁一致（最大の差は N=12 den=12 の 0.2271/0.2283） |
| final（step 500 の H⊥/H） | 概ね 3 桁一致。差が 1% を超えるのは N=5,6,10,12 の den=124 と N=5 den=5（飽和後の揺らぎ領域） |
| max\|Δlog10 f\|（全 step） | 0.01〜1.97。大きい組は N=5（親が別の解）、N=30、N=10、N=4 |
| N=3,4,5（親が原本と別の解） | onset・傾き・final が原本と一致。別の解が同じ力学を与えたことの記録（対称性同値かは未確認） |

**読み（事実の範囲）：** 環境差（macOS 26.3.1→26.6.2 の Accelerate）で親と状態の bit 一致は失われたが、インフレーション的発展の三指標（種スケールの出発、一定傾きの指数増幅、飽和）は 24 走行すべてで原本と同じ (N, den) の曲線に重なる。以後の「分配率 R の読出し」予備実験は、この `run_subset_stage123_v1.py` を起点に最小変更で行う。

### 生成物

- `results/hm_N{N}_den_{den}_states_500.npz` × 24（状態、原本と同名）、`results/timeseries_64bit_with124_N3_N40.csv`、`results/summary_64bit_with124_N3_N40.csv`（原本と同名、中身は 24 走行分）、`results/analysis_subset_summary_v1.json`、`results/RUN_METADATA_N3_N40_stage123.json`
- 図：`results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`（H⊥/H グリッド、未走行の N は空パネル）、`results/fig_compare_curves_control_vs_original.png`（原本との重ね描き）、`fig_complex_plane_{step0,final,final_zoom}_subset_control.png`（複素平面 3 枚、3×4）
- 照合：`results/compare_with_original_v1.json`（bit 一致と数値差）、`results/compare_curves_with_original_v1.json`、`results/compare_curves_table.md`
- ログ：`results/log_*.txt`。台帳：`SHA256SUMS.txt`
