# 初期AB状態の探索環境 v1

## 設計境界

このフォルダは散乱本体を変更しない。
`../run_ab_invariant_theta_toy_v1.py` の `theta_from_ab(A,B)` と
`rotate_ab(A,B,theta)` をそのまま利用する。

目標 $R$ は独立した初期状態探索で候補を比較するためだけに使う。
前進散乱へ渡すのは探索済みの初期A/B配列だけであり、目標 $R$ や
目標 $\theta$ は渡さない。

## 探索する初期値

- A: 単位振幅の広域基本波
- B: 等振幅奇数倍音 $1,3,\ldots,63$ の形は固定
- 探索変数: Bの初期振幅だけ

探索器は候補振幅で初期B配列を作り、既存の `theta_from_ab(A,B)` が
返した $R$ と目標値の差だけを見て二分探索する。読出し式を探索器へ
複製した解析的な逆算ではない。

## 実行

既定では $R=0.50$ と

$$
R_\alpha
=1-\sqrt{\frac{4\pi}{137.035999084}}
$$

の近傍になる初期状態を探索し、その後、目標値を受け取らない前進処理で
衝突回数 $0,1,2,3,5,10,20,42$ を図化する。

```bash
python3 \
  ab_invariant_theta_toy_v1/inverse_initial_conditions_v1/search_initial_conditions_and_plot_v1.py
```

任意の比較目標は `--target-r` を繰り返して指定できる。
