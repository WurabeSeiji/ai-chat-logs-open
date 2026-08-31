# 独立辺 k による rank-3 幾何監査 — 2026-08-31

## 目的
Google Drive の hm_N3〜hm_N16 `states_treatment.npz` 生データだけを入力とし、旧 `hm_series_k.csv` や群一様 k の結果を入力に使わず、辺別波長を再測定して 3 次元整数巻き数整合を探索した。

## 方法
- 全 40001 step を読み込む。
- 終端の 8192/4096/2048/1024 step の隣接窓で位相傾斜周波数を比較し、median <=2%, p90 <=5% を満たす最長安定窓を採用。満たさない場合は最小変動窓を fallback 採用。
- 各辺について lambda = max(f)/f を独立に測定。
- `L_ij = k_ij lambda_ij / 2`, `k_ij` は各辺独立整数、1..100、min k=1。
- 3D 座標と k を交互最適化し、最後に整数長 L だけから centered Gram を再構成して rank/PSD を検査。
- PASS は整数長そのものが target rank (N=3 は2、それ以外3)、PSD、座標残差 rms<1e-7, max<1e-6 を同時に満たす場合のみ。

## 実行結果
N=3,4 は PASS。N=5〜16 は今回の有限探索では PASS 未発見。これは不存在証明ではない。特に N>=14 は実行時間制約のため restart 数を段階的に減らしたため、探索強度は低い。

N=5 は rank=3 候補まで近づくが整数長 Gram に負固有値が残り、rms 約3e-3 で不合格。N=6 以降は残差・PSD違反が増え、N>=9 では k=100 上限に達する候補が多い。したがって kmax=100 が探索境界として効いている可能性がある。

## 重要な解釈
1. 旧「全 N で simplex 成立」は rank-3 創発の検証ではない。
2. 今回は同波長群に同じ k を強制していない。
3. N=5〜16 の `PASS=False` は「rank3 が存在しない」ではなく「宣言した有限探索で厳密整合解を発見できなかった」。
4. 次段階では kmax 拡張、複数 seed、探索器の独立実装、波長不確かさを許容した interval/weighted fit、null control が必要。

## 再現性
正本プログラム: `program/independent_edge_k_rank3_search.py`。集約結果: `results/rank3_search_summary_20260831.csv`。個別出力は `results/per_n8/N*/`。
