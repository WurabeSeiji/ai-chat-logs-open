# 実験——最新版の `make_parent` だけを振幅込み K で自己無撞着にした版（N=5、他は最新版と同一）

**作成日:** 2026-08-28　**ベース:** 木原最新版 `N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828`。
**変更:** `program/original_engine.py` の `make_parent`（と補助関数 `_adjacency`, `_K_amplitude_aware`, `_selfconsistency_residual`）のみ。`run_amplitude_only_fix.py` は最新版と同一（diff なし）。差分は `results/diff_engine_vs_latest.patch`（make_parent の区間 1 hunk）。

make_parent：K_amp(v) = Im(v̄ᵢvⱼ) に対して K_amp(v)v = −iσ_max v（旧実装と同じカイラリティ）を、位相・振幅とも未知数として混合反復で解く。正規化なし（スケールは固有値ソルバの ‖u‖=1）。残差はスケール不変、tol=1e-10、iters=2000、restarts=20、全失敗なら abort。

## 結果

| 版 | 親残差 | H_total | baseline onset(5%) / max H⊥/H / final PR/M | treatment onset(5%) / max H⊥/H / final PR/M / 振幅範囲 |
|---|---|---|---|---|
| 最新版 無変更再実行 | 4.87e-01 | 0.8596 | 334 / 0.716 / 1.000 | 76 / 0.943 / 0.296 / 3.9e-3..0.655 |
| (a) make_parent 正規化復元 | 1.82e-09 | 1.0000 | 196 / 0.624 / 1.000 | 65 / 0.991 / 0.307 / 1.9e-3..0.695 |
| **(b) make_parent 自己無撞着（本実験）** | 2.63e-11 | 1.0000 | 7 / 0.708 / 1.000 | **None / 4e-16 / 0.486 / 3.8e-16..0.564** |

- 自己無撞着状態は 1 リスタート・100 反復で収束（残差 2.6e-11）。中身は頂点 4 を切り離した K₄（6 辺）で、4 辺は振幅 0、K_amp(v) のスペクトルは ±0.5 のみ。
- treatment：**40000 step 不動**（H⊥/H ≈ 10⁻¹⁶）。安定な固定点。潜伏も急拡大もない。
- baseline：この v は位相のみ K の固有モードではないので step 7 で離脱し等分配へ（旧力学の性質）。

## ファイル
`program/`、`data/`、`figures/`、`results/`（diff、run.log、run_progress.log）、`run_all.sh`、`SHA256SUMS.txt`
