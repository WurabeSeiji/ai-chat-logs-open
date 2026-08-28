# 対照実験——最新版の `make_parent` だけを公開エンジンの原文（振幅正規化あり）に戻した版（N=5）

**作成日:** 2026-08-28　**ベース:** 木原最新版 `N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828`。
**変更:** `program/original_engine.py` の `make_parent` 内 1 行（`v = v/np.linalg.norm(v)` を復元）のみ。復元後のエンジンは公開エンジン `run_n_scaling_lowrank_v1_no_sigma_norm.py` と **全体が同一**（diff なし）。`run_amplitude_only_fix.py` は最新版と同一（diff なし）。差分は `results/diff_engine_vs_latest.patch`（1 hunk）。

## 結果（最新版との比較）

| 版 | 親残差 | H_total=|v|² | baseline onset(5%) / max H⊥/H / final PR/M | treatment onset(5%) / max H⊥/H(step) / final PR/M / 振幅範囲 |
|---|---|---|---|---|
| 最新版（保存データ、ChatGPT 機） | 4.87e-01 | 0.8596 | 347 / 0.704 / 1.000 | 76 / 0.923 (15526) / 0.302 / 4.0e-3..0.649 |
| 最新版 無変更再実行（この計算機） | 4.87e-01 | 0.8596 | 334 / 0.716 / 1.000 | 76 / 0.943 (29058) / 0.296 / 3.9e-3..0.655 |
| **(a) make_parent 正規化復元（本実験）** | 1.82e-09 | 1.0000 | 196 / 0.624 / 1.000 | 65 / 0.991 (13203) / 0.307 / 1.9e-3..0.695 |

- 親残差が意味のある値（1.8e-9、正規化ありなので残差式が正しく働く）に戻り、|v|=1。
- baseline（位相のみ K）：潜伏→急拡大→等分配（onset 196。最新版の 334/347 と違うのは丸め床からの出発時刻が親のスケール・計算機で変わるため。成長率・等分配は同じ）。
- treatment（振幅込み K）：最新版と同じ挙動（即離脱、H⊥/H → 0.99、PR/M ≈ 0.30、振幅の極端な不均一化）。時間尺度は |v|² が 0.86→1.00 になった分だけ速い（onset 76→65）。
- 結論：**make_parent の正規化を戻しても、treatment の「潜伏なし・即離脱」は変わらない。** 正規化除去（A1）は treatment の挙動の原因ではない。

## ファイル
`program/`（エンジン＝公開エンジン原文、実行スクリプト＝最新版原文）、`data/`、`figures/`、`results/`（diff、run.log、run_progress.log）、`run_all.sh`、`SHA256SUMS.txt`
