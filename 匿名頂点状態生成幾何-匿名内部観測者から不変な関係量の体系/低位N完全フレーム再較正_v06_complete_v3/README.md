# lowN_recalibration_v06 — complete

匿名頂点状態生成幾何 v0.6 の完全不変量フレームを Phase B 前に低位 N で再較正する正本パッケージ。

## 条件

- N = 2,3,4,5,6
- 各 N 2000 標本
- Re(z_j), Im(z_j) ~ Normal(0,1)
- seed = 20260914

## 実行

```sh
./run_all.sh
```

## 同梱

- `run_lowN_v06_recalibration.py` — 数値較正正本
- `verify_kappa_exact.py` — κ 厳密値との比較
- `run_all.sh` — 一括実行
- `make_sha256.py` — SHA256SUMS 再生成
- `lowN_v06_samples.csv`
- `lowN_v06_summary.csv`
- `N3_identity_checks.csv`
- `N3_identity_summary.json`
- `kappa_exact_values.csv`
- `kappa_exact_comparison.csv`
- `analysis.md` — 完全分析
- `kappa_exact_distribution.md` — κ の解析証明
- `SHA256SUMS.txt`

既存の力学正本プログラムには一切触れない。


## 独立再現検証

`verify_beta_law.py` は論文2 §14 の再現性主張をパッケージ内で再現する。

- 同一 seed=20260914 で abs_P1・κ の系列を再生成し保存CSVと突合
- 別実装・別 seed=12345 で不変性・Newton復元・平均κを対照検算
- κ=(λ1−λ2)/(λ1+λ2) 恒等式をランダム構成で機械精度検査
- m=2,3,4,5,6,39 について κ²~Beta(1,(m−1)/2) を KS 検定
- 厳密平均・中央値との突合

既定の KS 標本数は各 m 1,000,000。
