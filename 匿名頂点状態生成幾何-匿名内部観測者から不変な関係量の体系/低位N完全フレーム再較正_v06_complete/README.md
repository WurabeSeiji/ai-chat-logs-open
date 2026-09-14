# lowN_recalibration_v06

匿名頂点状態生成幾何 v0.6 の低位 N 再較正実験。

## 実行
```sh
./run_all.sh
```

## 正本条件
- N = 2..6
- 各 N 2000 標本
- 複素ガウス乱数 Re,Im ~ Normal(0,1)
- seed = 20260914

## 同梱
- run_lowN_v06_recalibration.py
- run_all.sh
- make_sha256.py
- lowN_v06_samples.csv
- lowN_v06_summary.csv
- N3_identity_checks.csv
- N3_identity_summary.json
- analysis.md
- README.md
- SHA256SUMS.txt

既存の力学正本プログラムには一切触れない。
