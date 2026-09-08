# 高対称理論床 simplex 次元監査（2026-09-08）

## 目的
高対称理論床の den=N 正本走行について、N=8,16,24,32,40 の step 0〜500 を用い、複素二乗辺値 `D_ij = z_ij^2` から simplex の数値的退化／次元開口を監査する。

## 入力
`inputs/hm_N{N}_den_{N}_states_500.npz`（N=8,16,24,32,40）。各 NPZ は `Z.shape=(501,M)`、辺順は `np.triu_indices(N,k=1)`。
元データは `対称親v2_500step走行_20260906/results/` の正本ファイルを複製したスナップショット。

## 定義
頂点0を基準に、複素二乗距離行列 `D` から

`G_ij = (D_0i + D_0j - D_ij)/2`,  i,j=1..N-1

を構成する。`G` は複素対称Gram型行列であり、最大rankは N-1。rank判定は特異値比 `sigma/sigma_max` に対して 1e-8, 1e-10, 1e-12 の3閾値を併記する。

Cayley–Menger 行列も独立に構成し、退化の有無を確認する。複素・有限精度データなので「rank」は厳密代数rankではなく数値的 effective rank として扱う。

## 実行
```bash
./run_all.sh
```

## 出力
- `simplex_dimension_summary.csv`: step0/500 と開口時刻の要約
- `simplex_dimension_timeseries.csv`: 5系列×501step の全監査値
- `simplex_singular_values_step0_step500.csv`: step0/500 の全特異値
- `input_manifest.csv`: 入力NPZのSHA256
- `高対称理論床_simplex次元保存性監査_20260908.md`: 分析本文
- `SHA256SUMS.txt`: 成果物ハッシュ

## 注意
本監査は `D_ij=z_ij^2` を複素二乗距離として読む本研究系列の定義に基づく。通常の正定値Euclid距離の埋込み次元とは同一視しない。
