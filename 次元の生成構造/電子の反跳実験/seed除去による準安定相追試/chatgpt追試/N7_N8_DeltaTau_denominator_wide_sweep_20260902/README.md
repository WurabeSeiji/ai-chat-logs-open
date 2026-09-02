# N=7,8 における Δτ = 2π/D の広域分母スイープ

実行日: 2026-09-02
指図: `CLAUDE_CODE_EXPERIMENT_INSTRUCTION_D_over_N_tau_sweep_20260902.md`
実行環境: macOS 26.3.1 arm64 / Python 3.9.6 / NumPy 2.0.2

## 入力（正本は読み取りのみ・本フォルダにコピーして使用）

- `data/N7/parent_v.npz` — 正本 `../干渉保存力学_資格審査とシード無し系列_20260831/hm_mp_free_N3_N40_20260901/data/hm_N7/parent_v.npz` のコピー
  SHA256 `39f2a2419be63b148859f817e226683041228a7d47b1c6d49d098de0c95a840b`
- `data/N8/parent_v.npz` — 同 `hm_N8/parent_v.npz` のコピー
  SHA256 `9412ec6d465325cc51ec969aa9caa12317efea8f9a39d4bddc5a278f8b691696`

コピーと正本の SHA256 一致を実行前に確認済み。既存正本フォルダへの書き込みは一切ない。

## 力学（指図 §2 のまま・変更なし）

H_ef = A_ef conj(z_e) z_f、1 step: z' = exp(-i(2π/D)H)z を `numpy.linalg.eigh` で厳密適用。
complex128/float64。seed・clipping・renormalization・ノイズ・状態分岐なし。

## 走行

- Stage A: N∈{7,8} × D∈{2..256 全整数, 320, 384, 512}（258 値）× 500 step。計 516 run。
- Stage B: N∈{7,8} × 重点 46 D 値 × S_max=⌈500·D/N⌉ step（同一 τ 窓 τ_max=500·2π/N）。計 92 run。
- Stage C: 機械的異常検出（`results/anomaly_followups.csv`）が 0 件のため追試なし。
- 全 608 run が status=ok。ノルム相対 drift 最大 ~2.3e-12（failure 閾値 1e-7 未到達、warning 閾値 1e-10 超は Stage B 長時間走行の一部のみ）。

## フォルダ構成

- `program/run_sweep.py` — 走行（checkpoint 方式: run 単位で再開可能）
- `program/analyze_sweep.py` — 集約（onset 3閾値・成長率 fit・飽和統計・異常検出）
- `program/plot_sweep.py` — 図化（主図は τ / D/N / χ 軸、step 軸は audit_step_axis/ のみ）
- `data/N{7,8}/D{DDDD}/timeseries_stage{A,B}.csv` + `run_info_stage{A,B}.json`
- `results/` — sweep_summary.csv, stageA/B_summary.csv, anomaly_followups.csv, initial_state_audit_stage{A,B}.json
- `figures/fig01..fig08` + `figures/audit_step_axis/`
- `ANALYSIS.md` — 指図 §14 の 10 問への回答と仮説判定

## 再現

```sh
./run_all.sh
```

（Stage A → 集約 → Stage B → 集約 → 図化。既完了 run は run_info の status を見てスキップ）

## 記録済みの数値監査上の注意

- 集約時 `A @ coef`（最小二乗の残差計算、n≥500 程度）で macOS Accelerate BLAS が
  RuntimeWarning（divide by zero / overflow / invalid in matmul）を発することがある。
  該当 4 run（N7 D32/D40/D48, N8 D40 の Stage B fit）について入力・係数・出力の全有限性と
  残差（max|resid| ≤ 0.02）を個別検証済み。擬似 FP フラグであり数値結果に影響しない。
- 時系列 CSV の `herm_rel_max` 列は「その時点までの H−H† 相対ノルム最大値」（監査列）。
- ヒートマップ（fig07/08）は実サンプルの zero-order hold で描画（平滑化・列間補間なし、
  カバー範囲外はマスク）。
