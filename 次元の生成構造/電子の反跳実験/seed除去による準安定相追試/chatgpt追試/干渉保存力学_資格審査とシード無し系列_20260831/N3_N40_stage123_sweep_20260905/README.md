# N=3..40 段1+2+3 スイープ（2026-09-05）

N=40 ラインで確立した段1+2+3 力学（位相のみ・虚部のみ生成子 i·K、実直交回転 exp(Δτ·K)、
固定 Δτ=2π/den）を、make_parent 静的親から N=3..40 の全域でスイープする。

## ファイル

- `run_n_scaling_lowrank_v1.py` — 7月正本エンジンの bit 同一コピー（diff 確認済み）
- `make_static_parents_N3_N40_v1.py` — 静的親生成。7月 N=40 走行と同一の手順・rng 式
  `default_rng(40260722+1000N+0)`・引数（iters=1200, tol=1e-12, δ=1e-15）。既存ファイルは
  上書きせず検証のみ。**N=40 は正本静的親と bit 一致がゲート**
- `run_N3_N40_stage123_v1.py` — スイープ本体。`run_N40_staticparent_imK_v1.py`
  （ChatGPT_denominator_controls_N40_selfcontrol_20260904、3ゲート合格済み）からの差分は
  「ループ range(3,41)・親の per-N 読込・出力先・図名・メタデータ名」のみ。力学は無変更。
  各 N×6分母（N−2..N+2, 124）×500 step、状態 npz 全保存

## 進行記録

### 2026-09-05 静的親生成 完了

`python3 make_static_parents_N3_N40_v1.py`

- 38 親（N=3..40）すべて生成・全収束（residual 3.5e-14〜9.4e-13、全て ok）
- **GATE PASS**: N=40 の v/g/Z0 が正本
  `自発的分裂予備実験_v1_N40対照実験系_20260904/.../parent_static_N40_makeparent_20260904.npz`
  と bit 一致
- 台帳: `parents/parents_summary.csv`（N・M・残差・平面数・status）
- rank_planes は概ね N（N=3:2、N=4:2、N=5:4、N=6:6、N≥7: N）

### 2026-09-05 スイープ本体 完了

`python3 run_N3_N40_stage123_v1.py` → `ALL DONE`（exit 0、実所要 ~45分）

- 出力: `results/` に timeseries CSV（8.7MB）・summary CSV・グリッド図・RUN_METADATA・
  状態 npz **228個**（38N×6分母）
- **INPUT GATE PASS**（check_sweep_inputs_v1.py）: 全228走行の Z[0] が各 N の静的親 Z0 と bit 一致
- **全 N=3..40 で緩和曲線（種スケール→直線指数増幅→飽和 ~10⁻³〜10⁻²）を確認**。
  228走行中 212 が 500 step 内に 0.05 交差
- レートの分母依存に系統構造: 小 N では 2π/124 が最遅、N≳20 では 2π/124 が最速に逆転。
  N+2 系が概ね速く N−1 系が遅い傾向
- 図: `results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`

### 2026-09-05 複素平面読出し図（plot_complex_plane_N3_N40_stage123_v1.py）

Δτ=2π/N の状態 npz から。既存様式（グリッド版＋最大クラスター拡大）を N=3..40 の 8×5 に拡張。

- `fig_complex_plane_step0_N3_N40_stage123.png` — step0（make_parent 静的親）
- `fig_complex_plane_final_N3_N40_stage123.png` — 終了時（step500）
- `fig_complex_plane_final_zoom_N3_N40_stage123.png` — 終了時の最大角クラスター拡大

## 一括再現

```bash
./run_all.sh   # 親生成 → スイープ → 入力ゲート → 図化
```

## 環境

`.venv/bin/python3`（Python 3.9.6、numpy 2.0.2、macOS arm64 Accelerate）
