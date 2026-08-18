# シードあり系列（第7論文 δ=1e-15）—— 30桁の駆け上がりの復元

インフレーション的増幅が第7論文の図から消えた原因の特定と、その復元。
**シードなし系列は `../paper8_a2a_seedless_v1/` にある。混同しないこと。**

## 結論

断絶点は **`run_paper7_5color_timeseries.py:118`** の1行。

```python
f = 1.0 - E_P1 / totZ          # 固定親基底からの引き算。桁落ちする形
```

原本は良い式を**同じ関数の中に持っている**（86–88行 `fval`、射影形）が、
crossing 判定にしか使われず CSV に記録されない。図の `set_ylim(0,1)` /
`clip(1e-6)` は、$10^{-16}$ 床で 0 や負値を含む CSV を対数軸に載せるための
**後追い対処**であって原因ではない。

同一走行・同一ステップ・同一精度での実測（N=5）:

| 式 | 初期値 | 最小値 |
|---|---|---|
| 射影形（記録されない） | 1.06601505920675953e-30 | 同左 |
| 引き算形（CSV記録列） | 2.220446e-16 | **-4.440892e-16**（負値） |

射影形の初期値は原点データ `fcurve_N00005_delta1e-15_seed0.csv` の
f(0)=1.0660150592067595e-30 と一致する。**第7論文の走行は30桁を持っていた。**

## 条件（シードなし系列と共通／相違）

| 項目 | 本系列 | `../paper8_a2a_seedless_v1/` |
|---|---|---|
| N | 5（M=10） | 5（M=10） |
| 親 PRNG seed | 40265722 | 40265722（**同一**） |
| make_parent | iters=1200, tol=1e-12 | iters=1200, tol=1e-12（**同一**） |
| 親残差 | 2.139898e-13 | 2.139898e-13（**同一**） |
| **初期状態** | **Z0 = (v + 1e-15·g)/‖·‖** | **Z0 = v.copy()** |
| f の定義 | ‖Z−p(p·Z)−q(q·Z)‖²/‖Z‖² | 同左（**同一**） |
| f(0) | 1.066015e-30 | 3.274787e-33 |
| 記録範囲 | 0–5000 毎ステップ | 0–5000 毎ステップ（**同一**） |
| 図の ZOOMS | 0-5000 / 0-250 / 0-25 | 同左（**同一**） |

## プログラム（すべて原本を import、変更点は ★ コメント付き）

| ファイル | 役割 | 検証 |
|---|---|---|
| `run_control_paper7_5color_v1.py` | CTRL-2 純再現 | 公開CSV **2202行×16列 バイト一致**（md5 a510ec9e…） |
| `run_paper7_5color_dual_f_v1.py` | 16列＋射影列（se_ev=25, XMAX=55000） | 既存16列が公開CSVと全行一致 |
| `run_paper7_5color_everystep_v1.py` | 毎ステップ 0–5000 | 25の倍数 201行×10列 一致（下記注意） |
| `make_paper7_figures_control_v1.py` | 原本図化の忠実コピー（パス2箇所のみ） | 公開図と **md5 一致**（比較図4枚） |
| `make_paper7_figures_projection_v1.py` | 射影形・原本レイアウト | — |
| `make_paper7_figures_projection_stacked_v1.py` | 射影形・縦2段（間引き25データ） | — |
| `make_paper7_figures_everystep_v1.py` | 毎ステップ・縦2段 | — |
| `make_paper7_figures_everystep_v2.py` | 毎ステップ・**縦3段（A2a と同一仕様）** | — |

### 注意: 記録間隔に依存する列がある

原本 `run()` は方向3・4 の基底を `align_2d(f_prev, f34)` で**前回記録時**の基底へ
整列させる（docstring「縮退平面は連続基底固定」）。f_prev は記録した step でしか
更新されないため、**記録間隔を変えると方向3・4 の割り当てが変わる**。

実測（step 0..3000 の 25 の倍数 121 点、se_ev=25 vs se_ev=1）:

```
splitting_fraction              不一致 0
splitting_fraction_projection   不一致 0
plane_1_occupation              不一致 0
kernel_occupation               不一致 0
direction_3_occupation          ★ 120 / 121
direction_4_occupation          ★ 120 / 121
```

バグではなく原本の設計。毎ステップ版の照合は f_prev 非依存の10列に限定している
（`VERIFY_COLUMNS` / `SAMPLING_DEPENDENT`）。方向3・4 が要るときは間引き25の
`dual_f_timeseries_N00005.csv` を使う。

## 実行順序

```bash
python3 run_control_paper7_5color_v1.py          # 固定点の確保（必ず最初）
python3 run_paper7_5color_dual_f_v1.py 5
python3 run_paper7_5color_everystep_v1.py 5
python3 make_paper7_figures_control_v1.py        # 公開図と md5 一致を確認
python3 make_paper7_figures_projection_v1.py
python3 make_paper7_figures_projection_stacked_v1.py
python3 make_paper7_figures_everystep_v1.py
python3 make_paper7_figures_everystep_v2.py      # 縦3段（A2a と同一仕様）
```

CSV は リポジトリの .gitignore（*.csv）で追跡外。上記で決定論的に再生成できる。

## 主要な図

- `figures_control/` —— 元の図（公開図と md5 一致）
- `figures_everystep_v2/figureZ_stacked_zoom_seeded_N00005.png` —— 30桁の縦3段
