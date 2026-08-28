# N=5 長時間追試（5000 → 40000 step）

## 変更点
直前の seedless + make_parent 内振幅正規化除去版から、`STEPS` だけを 5000 から 40000 へ変更した。その他のコード・初期条件・相互作用則は変更していない。

## 条件
- N=5, M=10
- L=124
- SEED=0, DELTA=0
- Z0=v（外部 perturbation seed なし）
- 初期 Z 正規化なし
- make_parent 内 v 正規化なし
- 線形指数回転
- treatment: K_ij = Im(conj(Z_i) Z_j)

## 実測結果
make_parent residual = 0.4868999851906601。

### phase-only baseline
- H_perp/H_total > 0.05: step 347
- H_perp 最大: 0.6052085420734356 (step 435)
- 最大 H_perp/H_total: 0.7040180999628793
- step 40000: H_perp = 0.3707776260740153
- step 40000: PR/M = 1.0
- step 40000: |Z^T Z| = 2.8531058002805865e-13
- H_total 最大ドリフト = 3.716804641840099e-12

### amplitude-aware treatment
- H_perp/H_total > 0.05: step 76
- H_perp 最大: 0.7932642902113481 (step 15526)
- 最大 H_perp/H_total: 0.9227768273882487
- step 40000: H_perp = 0.5921311844072588
- step 40000: PR/M = 0.3020287734000109
- step 40000: |Z^T Z| = 1.2049262565475316e-13
- step 40000: 振幅範囲 = 0.004004484091925404 .. 0.6487521752952705
- H_total 最大ドリフト = 7.494005416219807e-13

## 観察
5000 step では見切れなかった長時間変動が続いている。amplitude-aware 系は固定点・単純飽和には入らず、step 15526 で H_perp/H_total ≈ 0.923 まで再上昇する。その後も step 40000 まで大きな状態再配分が継続し、PR/M は 0.302 まで低下している。

この実験だけから U^n=I の次数は決定しない。全複素状態の recurrence 距離を別途評価する必要がある。
