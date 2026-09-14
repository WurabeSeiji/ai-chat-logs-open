# stage1_ssb_terminal_Z_v1

Phase B 第一段（2026-09-14 木原承認）: 正本SSBバッチ（全N6_40_SSBバッチ_20260912）の忠実再走行に
終端Z保存のみを追加し、対照突合（Tier-0=テキスト完全一致で合格）後、5seed有効標本評価を実施。

- run_all_save_terminal_Z.py — 忠実コピー＋terminal_Z保存（物理式・シード・停止条件は不変）
- compare_against_canonical.py — 事前宣言合格条件（Tier-0/Tier-1）による突合
- evaluate_effective_samples.py — n_raw/n_unique・Bob対距離・事前登録予測の判定
- terminal_Z/ — N{6..40}_terminal_Z.npz（175終端状態の正本保存版）
- analysis_stage1.md — 分析（Tier-0合格・予測の部分反証・κ_B終端値）

再現: `sh run_all.sh`（numpy, pandas。約20分）。正本フォルダ・力学正本には一切書き込まない。
