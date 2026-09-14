# main_estimation_5seed_v1

Phase B 主推定（設計書 v1.2 §13 凍結プロトコル・5seed・stage1 正本 terminal_Z 使用）。

- run_main_estimation.py — K飽和×A/B×E1-E4×seed-cluster bootstrap（実装細部は冒頭に実行前宣言）
- make_figures.py — dimension_vs_N / d(K)曲線 / κ_B null 対照図
- main_results_by_N.csv / d_est_N_K_curves.csv / kappaB_null_contrast.csv
- analysis_main_estimation.md — §20 判定：rank-4 非支持、N≳20 平坦域は標本律速、κ_B null 棄却

再現: `sh run_all.sh`（numpy, pandas, matplotlib。数分）。読出し専用・力学正本無変更。
