# make_parent型初期値：低N構造解析と500 step再現（2026-09-07）

## 1. 目的

以前の「高対称理論床」系列は、重心閉塞や巡回対称性を強く作り込んだため、canonical `make_parent` が自然に選ぶ枝と異なる可能性がある。本検討では、実際に保存されている canonical `make_parent` 親 N=3..7 を直接用い、位相・振幅構造を監査したうえで、同じ stage123 力学で 500 step の低N再現を行う。

重要な方針は、`make_parent` に存在しない拘束を後付けしないことである。特に `sum z = 0` を必須条件にしない。

## 2. 使用した実データ

元データは `N3_N40_stage123_sweep_20260905/parents/parent_static_N{N:05d}_makeparent_20260905.npz` の N=3..7。各ファイルの `Z0` を初期状態として使用した。

親生成の canonical 手順は同梱 `make_static_parents_N3_N40_v1.py` の通りで、

```text
rng = default_rng(40260722 + 1000*N)
v, residual, sig = make_parent(..., iters=1200, tol=1e-12)
g = zero_closure_kernel_seed(...)
Z0 = (v + 1e-15*g) / ||v + 1e-15*g||
```

である。

## 3. 位相構造

実データを直接解析すると、各 N の位相は任意の大域位相 `phi0` を除いてほぼ Z4 四相、

```text
theta_e - phi0 in {0, pi/2, pi, 3pi/2}
```

にロックされる。

ここで重要なのは、四相を手で指定していないことである。`make_parent` の自己無撞着反復がこの枝を選んでいる。

## 4. 振幅構造

振幅は等振幅ではない。実データから

```text
W_ef = A_ef * sin^2(theta_f - theta_e)
```

を作ると、`r_e = |z_e|` は

```text
W r = sigma r
```

を高精度で満たす。したがって canonical `make_parent` 型の本質は、

1. Z4 四相の位相枝
2. その位相枝が作る 90度接続グラフ W
3. W の正固有ベクトルとしての非等振幅

の組である。

単純な「四相＋等振幅」は canonical `make_parent` の再現ではない。

## 5. 二乗閉塞

四相では `z_e^2` の符号が二軸で反転するため、

```text
sum z_e^2 = 0
```

は二軸パワーの等分配

```text
P_even = sum_{q=0,2} r_e^2
P_odd  = sum_{q=1,3} r_e^2
P_even = P_odd
```

として現れる。低N実データではこの関係が機械精度で成立する。

## 6. 重心閉塞は canonical 条件ではない

実 `make_parent` 親では一般に `|sum z| != 0` である。したがって、以前の高対称系列のように `sum z=0` を先に課すと別枝へ移る。この差は、N=3の安定性差を含む枝選択の違いとして今後比較すべきである。

## 7. 500 step 再現条件

同梱 `run_makeparent_lowN_validation_20260907.py` は canonical stage123 の one-step 定義をそのまま用いる。

```text
H = i K,
K_ef = A_ef sin(theta_f-theta_e),
z_{t+1} = exp((2*pi/den) K) z_t
```

N=3..7、分母 `N-2,N-1,N,N+1,N+2,124`、500 step を実行した。

## 8. den=N の結果

| N | Hperp/H > 0.05 初回 step | step500 Hperp/H |
|---:|---:|---:|
| 3 | 45 | 0.166666666666285 |
| 4 | 50 | 0.166666666666227 |
| 5 | 64 | 0.370169547484929 |
| 6 | 70 | 0.235268632883814 |
| 7 | 88 | 0.088341072305993 |

N=3を含む全低Nでインフレーション的な横方向成分の成長が発生した。

これは高対称理論床の N=3 が閉塞面内で強く復元した結果と対照的であり、`N=3だから安定` ではなく `枝依存` であることを示す。

## 9. 複素平面

`lowN_validation/results/fig_complex_plane_step0_step500_denN_N3_N7.png` に den=N の step0/500 を並べた。step0 は canonical `make_parent` 特有の軸偏り・非等振幅を持ち、step500 では位相再配置と半径の再編が生じる。

## 10. 現段階の結論

確認済み：

- canonical `make_parent` の低N実データは Z4 四相枝に落ちる。
- 振幅は等振幅ではなく、90度接続グラフ W の正固有ベクトル構造を持つ。
- 二乗閉塞は成立する。
- 重心閉塞は一般には成立しない。
- N=3..7 の den=N は全て500 step内に Hperp/H > 0.05 を越える。

未証明・今後の課題：

- Z4 四相ロックが一般Nで自己無撞着方程式から解析的に必然か。
- W の正固有ベクトル構造と `make_parent` の元の JG 固有問題の厳密な同値関係。
- 高対称理論床系列との枝分類・線形安定性の統一的説明。

## 11. 同梱ファイル

- `make_static_parents_N3_N40_v1.py` — canonical 静的親生成手順
- `run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py` — canonical 500 step スイープ参照コード
- `run_makeparent_lowN_validation_20260907.py` — 今回の低N再現
- `plot_makeparent_lowN_20260907.py` — 今回の図化
- `analyze_makeparent_parent_structure_20260907.py` — 位相・振幅・W固有ベクトル監査
- `parents_actual_N3_N7/` — 実際に使用した親データ
- `parent_structure_summary_N3_N7.csv` — 構造監査結果
- `lowN_validation/results/` — 500 step 生データ、集計、図
