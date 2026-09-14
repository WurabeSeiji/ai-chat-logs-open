# 低位 N 数値実験 v0.6 再較正 — 分析

匿名頂点状態生成幾何 v0.6 の完全不変量フレームを Phase B 前に低位 N で再較正する。

- N = 2,3,4,5,6
- 各 N 2000 標本
- Re(z_j), Im(z_j) ~ Normal(0,1)
- seed = 20260914

## 結論

G_B 不変性、K=N-1 での Newton 復元、N=3 恒等式群を再検証する。
数値詳細は `lowN_v06_summary.csv` と `N3_identity_summary.json` を正本とする。
