# make_parent 段階1: 120度三値構造の導出（第4検証）

## 結果
第3検証で未証明だった

`exp(3 i (phi_e-theta_e)) = const`

を、収束軌道に現れる二つの「射影的周期条件」に還元した。

N=3..10、20 seed、2000反復のうち第2検証で収束判定した128ケースすべてについて、連続する二ステップを

- `theta0 -> (phi0) -> theta1`
- `theta1 -> (phi1) -> theta2`

と書くと、次を確認した。

1. 入力位相の射影的2周期

   `exp(2 i (theta2-theta0))` は全辺で共通。

   これは `theta2 = theta0 + global (mod pi edgewise)` を意味する。

2. 出力固有モード位相の射影的不変性

   `exp(2 i (phi1-phi0))` は全辺で共通。

   これは `phi1 = phi0 + global (mod pi edgewise)` を意味する。

3. beta=1/2 の中点恒等式

   `exp(2 i theta1)=exp(i(theta0+phi0))`

   `exp(2 i theta2)=exp(i(theta1+phi1))`。

## 解析的導出
辺に依存しない全体位相は `~` で同一視し、各辺の `mod pi` で計算する。

射影的2周期と出力位相不変性より

`theta2 ~ theta0`, `phi1 ~ phi0`。

中点式を mod pi で書けば

`2 theta1 ~ theta0 + phi0`  ...(1)

`2 theta2 ~ theta1 + phi1`  ...(2)

(2) に周期条件を代入すると

`2 theta0 ~ theta1 + phi0`  ...(3)

(1) から `phi0 ~ 2 theta1-theta0`。これを (3) へ代入すると

`2 theta0 ~ theta1 + 2 theta1-theta0`

したがって

`3(theta1-theta0) ~ 0 (mod pi)`。

ゆえに、ある全体定数 gamma に対して各辺 e で

`theta1_e-theta0_e = gamma + k_e*pi/3 (mod pi)`, `k_e in Z`。

一方 beta=1/2 の中点式から、非対蹠点では

`phi0-theta0 = 2(theta1-theta0) (mod 2pi)`。

よって

`phi0_e-theta0_e = 2 gamma + 2 k_e*pi/3 (mod 2pi)`。

従って

`exp(3 i(phi0_e-theta0_e)) = exp(6 i gamma)`

は全辺で共通となる。

これで第3検証の120度三値構造は、**二つの射影的周期条件 + beta=1/2 の中点恒等式**から解析的に従う。

## 128ケース数値検証
`projective_cycle_proof_converged128.csv` に全ケースを保存した。

検証量:
- `theta_projective_two_cycle_spread`
- `phi_projective_fixed_spread`
- 二つの midpoint identity error
- `sixth_step_spread = spread(exp(6i(theta1-theta0)))`
- `cube_delta_spread = spread(exp(3i(phi0-theta0)))`

すべて倍精度丸め誤差程度で成立した。

## 証明済み / 未証明の境界
今回解析的に閉じたもの:
- 上記二つの射影的周期条件を仮定すれば、120度三値構造は厳密に導かれる。
- 120度三値構造と第3検証の `K(phi)=-2K(theta)` を組み合わせた30/90度量子化、sigma 2倍則も従う。

まだ一般証明が必要なもの:
- なぜ段階1のスペクトル写像が収束時に
  `theta2 ~ theta0 (mod pi)` と `phi1 ~ phi0 (mod pi)`
  を必ず満たすのか。

したがって未解決点は「120度がなぜ出るか」からさらに縮小され、**収束軌道の射影的2周期性そのものの証明**になった。
