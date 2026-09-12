# 全N=6..40 SSBバッチ検証（2026-09-12）

各N（6..40）で floor Z0＋振幅1e-8の異なる向きの微小シード5通り→厳密one_stepで適応ロックまで走行。
「全シードが inflation→等振幅・Δ=−2π/N の同型平衡へ、位相配置(真空の向き)はシード依存で相異」＝SSBを個別確認。

結果: 35/35 で SSB確認。個別分析は per_N/N{N}.txt、集約は ssb_batch_summary.csv。
再現 `bash run_all.sh`（numpy）。読出し＋厳密one_stepの制御実験（初期シードのみ制御・物理式不変）。
