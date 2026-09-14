# pilot_long10000_N6_v1

Phase B 読出し pilot（設計書 v1.2 §23 手順、次元推定なし・読み取り専用）。
対象: N3_N40_long10000_20260905 の N=6 終状態（step 10000）。

- pilot_readout_N6.py — Bob-star 抽出・P_k・基準次数 r・Bob 間距離の検査
- pilot_bobstar_invariants_N6.csv — Bob ごとの不変量
- pilot_bob_pair_distances_N6.csv — Bob 対距離
- pilot_summary_N6.json — 集約
- analysis_pilot_N6.md — 分析

再現: `sh run_all.sh`（numpy, pandas）。元データ・力学正本には一切触れない。
