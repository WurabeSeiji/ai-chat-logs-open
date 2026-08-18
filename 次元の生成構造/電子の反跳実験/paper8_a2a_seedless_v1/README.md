# シードなし系列（第8論文 Stage A2a）

明示シードを置かない走行。**シードあり系列は `../paper7_f_projection_v1/` にある。**

## 条件（シードあり系列と共通／相違）

| 項目 | 本系列 | `../paper7_f_projection_v1/` |
|---|---|---|
| N | 5（M=10） | 5（M=10） |
| 親 PRNG seed | 40265722 | 40265722（**同一**） |
| make_parent | iters=1200, tol=1e-12 | iters=1200, tol=1e-12（**同一**） |
| 親残差 | 2.139898e-13 | 2.139898e-13（**同一**） |
| **初期状態** | **Z0 = v.copy()**（シードなし） | **Z0 = (v + 1e-15·g)/‖·‖** |
| f の定義 | ‖Z−p(p·Z)−q(q·Z)‖²/‖Z‖² | 同左（**同一**） |
| f(0) | **3.274787e-33** | 1.066015e-30 |
| 記録範囲 | 0–5000 毎ステップ | 0–5000 毎ステップ（**同一**） |
| 図の ZOOMS | 0-5000 / 0-250 / 0-25 | 同左（**同一**） |

**親状態 v は両系列で完全に同一。違いは初期状態だけ。** f(0) は 2.5桁違う。

## 凍結条件（原本 config_locked.json）

```
n=5, dtype=float64, parent_prng_seed=40265722,
parent_iters=1200, parent_tolerance=1e-12,
initial_state="Z0 = v.copy()", max_step=5000
forbidden: delta=True, zero_closure_kernel_seed=True, high_precision=True
```

`high_precision` が禁止されているのは、増幅が数値誤差の産物でないことが
第6論文の**生成子固定対照**で既に決着しているため（同じ底・同じ精度で、
自己参照の有無だけが30桁の有無を決める）。底の深さを多倍長で掘る実験は立てない。

## プログラム（すべて原本を import、変更点は ★ コメント付き）

| ファイル | 役割 | 検証 |
|---|---|---|
| `run_a2a_seedless_reproduce_v1.py` | データ再現 | 公開 `f_timeseries.csv` の **5001 step 全一致**（.17e 文字列） |
| `make_a2a_figures_control_v1.py` | 原本図化の忠実コピー（パス2箇所のみ） | 公開図 **14枚すべて md5 一致** |
| `make_a2a_figures_stacked_v1.py` | 上記＋縦積み拡大図（★変更3の1箇所） | 既存14枚はそのまま生成 |

原本 `run_seedless.py` は `raw/` が空でないと実行を拒否する設計のため直接は呼べない。
走行ループのうち f 系列を作る部分だけを転写し、力学は import している。

## 実行順序

```bash
python3 run_a2a_seedless_reproduce_v1.py    # 固定点の確保（必ず最初）
python3 make_a2a_figures_control_v1.py      # 公開図14枚と md5 一致を確認
python3 make_a2a_figures_stacked_v1.py
```

## 主要な図

- `figures_control/` —— 原本の図14枚（公開図と md5 一致）
- `figures_stacked/figureZ_stacked_zoom_seedless.png` —— 30桁の縦3段
