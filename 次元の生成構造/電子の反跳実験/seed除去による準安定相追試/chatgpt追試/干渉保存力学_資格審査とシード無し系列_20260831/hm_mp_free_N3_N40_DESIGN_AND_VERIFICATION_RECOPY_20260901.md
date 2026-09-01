# mp を使わない hm 初期値生成器 — 設計方針と検証記録

作成日: 2026-09-01

## 1. 目的

既存 hm 系列の構造を保ったまま、`mp` / `make_parent` に依存せず、任意の N に対して高対称・等振幅・干渉項を保持した自己無撞着な複素初期値を決定論的に生成する。

互換性上の要求は次の通り。

1. N=3..16 は既存 hm 実験の初期値スケールを変えない。
2. N>=17 は mp 由来の調整を持ち込まない。
3. 15 の使い方は元の hm 規約どおり `RBAR2 = 1/15` とする。
4. 位相配置は探索・乱数・Newton 法を使わず、完全グラフの対称構成から直接生成する。
5. 複素数 `v_ij` をそのまま状態とし、振幅・位相・干渉項を捨てない。

## 2. 位相配置

辺状態を

`v_ij = r exp(i theta_ij)`

とする。

### N=3

既存の Z3 特殊配置をそのまま使う。

### 偶数 N

K_N の 1-factor 分解を用いる。N-1 個のクラス c=0..N-2 に

`theta_c = pi c/(N-1)`

を割り当てる。各頂点には各クラスが1本ずつ接続するため、二乗位相は (N-1) 次単位根を一周する。

### 奇数 N >= 5

頂点を Z_N 上に置き、円環距離 d=1..(N-1)/2 で辺を分類する。q=(N-1)/2、クラス c=d-1 として

`theta_c = pi c/q`

を割り当てる。各頂点には各距離クラスが2本ずつ接続するため、二乗位相は q 次単位根を2回ずつ含む。

したがって N>=4 では各頂点 i で

`S_i = sum_{j != i} v_ij^2 = 0`

が解析的に成立する。全体二乗閉塞

`sum_{i<j} v_ij^2 = 0`

も従う。N=3 は局所閉塞ではなく既存 Z3 特殊解だが、全体二乗閉塞を満たす。

## 3. 振幅スケール

### N=3..16: historical compatibility layer

既存 `pass1_parents.py` は、解析的に作った hm を

`v = sp.equimodular(N)`
`v = v * NORM[N] / np.linalg.norm(v)`

で mp 親のノルムに合わせていた。

新生成器では mp を実行せず、既存実験を再計算可能にするため N=3..16 の `NORM[N]` だけを固定互換定数として保持する。これは理論定数ではなく、過去の実験初期値を保存するための compatibility data である。

### N>=17: intrinsic hm rule

調整値を一切使わず、

`r^2 = 1/15`

すなわち

`r = sqrt(1/15)`

を全 N に共通の規約として使う。したがって N>=17 では mp 由来スケールは存在しない。

## 4. 自己無撞着性

生成後、元の干渉 H と同じ

`H_ef = A_ef conj(v_e) v_f`

を直接構成し、

`H v = mu v`

の残差を検査した。

N=17..40 の生成データで得られた最大値:

- normalized global |sum z^2|: 1.43e-16
- normalized local max |S_i|: 7.46e-15
- H self-consistency residual: 8.39e-16

すべて binary64 の丸め誤差水準である。

## 5. N=3..16 互換性検証

新生成器は N=3..16 について、旧 `pass1_parents.py` と同じ

1. 位相生成式
2. `RBAR2=1/15`
3. `v * NORM[N] / np.linalg.norm(v)` という演算順序
4. 保存済み `NORM[N]`

を使う。

旧 `parents_predictions.csv` と新生成値を N=3..16 全点で比較した結果:

- `amp_min` 差: 全 N で 0
- `amp_max` 差: 全 N で 0
- `mean_amp2` 差: 最大 9.72e-17（CSV 表現を介した比較）
- `norm` 差: 最大 4.44e-16（CSV 表現を介した比較）

さらに旧 N=5 `parent_v.npz` と新 N=5 `parent_v.npz` の `v` 配列を直接比較し、

`np.array_equal(old_v, new_v) == True`

`max(abs(old_v-new_v)) == 0.0`

を確認した。

全 N の比較値は `legacy_N3_N16_verification.csv` に保存した。

## 6. 出力

`program/generate_hm_mp_free.py`

- デフォルト N=3..40
- `--n-min`, `--n-max`, `--out` で範囲・出力先を変更可能
- 乱数なし
- mpなし
- scipyなし
- Newtonなし

各 `data/hm_N*/` に以下を保存する。

- `parent_v.npz`: 複素初期値と位相クラス情報
- `parent_v.csv`: 全辺の実部・虚部・振幅・位相
- `parent_checks.json`: 二乗閉塞・局所閉塞・H自己無撞着性の検査値

全 N の一覧は `generation_summary.csv`。

## 7. 重要な扱い

N=3..16 の legacy norm は今後の理論へ外挿しない。N=17 以降は `1/15` の一定規約だけで生成する。これにより、過去の実験を壊さず、将来の N から mp 依存を切断する。
