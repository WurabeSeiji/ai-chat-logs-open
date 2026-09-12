# parity 二部構造による90度整数K床の解析的一本化（検証パッケージ）

日付: 2026-09-12。読出しのみ・物理無変更・新規走行なし。論文3の格上げ根拠。

## 目的

種ラベル $k_e=(i+j)\bmod4$ の parity $p_e=(i+j)\bmod2$ だけで、論文1の $W_{ef}=A_{ef}\sin^2\psi_{ef}$ が
二部グラフ $W_{ef}=A_{ef}\mathbf 1_{p_e\neq p_f}$ になることを使い、床・振幅・全スペクトル・重心の閉形式を
**任意 $N$**（sympy 記号）＋ **$N=3$–$40$**（数値）で確認する。§3 のラベル反復が床の構成に不要であること
（反復床の $W$ が種 parity の $W$ と機械精度一致、種直接構成 $z=D_{\rm seed}r_{\rm Perron}$ が符号ゲージを除き反復床と一致）も確認する。

## 実行

```
bash run_all.sh
```

- `verify_parity_bipartite_analytic_20260912.py`：物理正本 `../../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py`
  （SHA256 `1abf2353…`）から adjacency/one_step を抽出（照合・物理無変更）。既存 38 床
  `../parents_90deg_floor_analytic/*.npz` を入力。numpy, sympy 依存。終了コード 0 = 全 OK。

## 出力

- `results/parity_bipartite_summary.csv`：全 $N$ の検証量（W一致・releq・振幅一致・符号ゲージ・σ²・振幅CF誤差・重心実測/CF・rank/核・スペクトル一致）。
- `results/実行ログ_20260912.log`：数値表＋ sympy 記号検証の全出力。

## 結果（要約）

全 $N=3$–$40$：反復床 $W$ = 種 parity $W$（‖差‖ $\le3\times10^{-27}$）、種直接床は相対平衡（残差 $\le1.2\times10^{-14}$）で
符号ゲージを除き反復床と一致（‖差‖ $\le8\times10^{-15}$）。$\sigma^2$（偶 $N(N-2)$/奇 $N^2-2N-1$）・振幅²（parity クラス割当）・
全スペクトル（偶 $N\ge6$/奇 $N\ge7$）・残差重心（N mod 4）の閉形式が一致。rank$W=2(N-1)$・核 $(N-1)(N-4)/2$。
sympy 記号で縮約 2×2（偶）/3×3（奇）から σ²・振幅² を任意 $N$ で厳密確認。詳細は `分析_parity二部構造_解析的一本化_20260912.md`。
