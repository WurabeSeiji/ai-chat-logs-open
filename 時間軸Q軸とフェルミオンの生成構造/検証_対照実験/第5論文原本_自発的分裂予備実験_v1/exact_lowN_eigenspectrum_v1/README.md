# exact_lowN_eigenspectrum_v1 — N=5, N=40 厳密全固有値・全固有方向 観測（第1段階）

第7論文 第1段階の生データ・図・数値診断。**N=300 には着手していない。近似・反復・上位限定・
閾値削除・縮退平均は用いていない。解釈は書いていない。**

## 対象と条件（変更禁止条件を既存実験と一致）

- N=5（M=10）, N=40（M=780）
- seed = 40260722 + 1000·N、初期微小種 δ=1e-15、零閉鎖層、正規化、生成子構成、σ計算法
  （power反復）、Cayley更新、演算順序、crossing定義（f>0.05）、crossing後 50000 step。
- 原本エンジン `../run_n_scaling_lowrank_v1.py`（不変更 import）。初期親平面 P0 は
  `../run_plane_flow_exact_v1.py` の `parent_plane_split_exact` で確定。
- 軌道は同一条件で再生成し、既存 `../paper6_definitive_control_v1/obs_N*.csv` の f と一致を検証
  （diagnostics の `max_f_ref_dev`）。

## 固有分解（近似なし）

各記録時刻で実反対称生成子 K(t) を明示構成し `numpy.linalg.eig` で**完全**固有分解。
固有値 ±iσ_j、正の σ_j を降順、複素共役固有ベクトル対から各モードの実2次元回転平面 B_j を
正規直交化。ゼロ近傍固有値も削除せず核次元として計上。

## サンプリング時刻（事前固定）

t=0〜crossing+2000 は 5step 毎、crossing+2000〜+10000 は 50step 毎、+10000〜+50000 は 200step 毎、
及び終端。crossing 近傍を高密度化。

## 出力

```
code/
  run_exact_lowN_eigenspectrum_v1.py   # 観測本体
  make_figures_exact_lowN_v1.py        # 図生成
raw/N00005/, raw/N00040/
  eigenvalues.csv          # 全時刻・全順位: σ, σ/σ1, Nσ/σ1, N²σ/σ1, log10, 固有値実虚, solver_residual (%.17e)
  delta.csv                # rank, branch, σ, σ/σ1, overlap_with_parent, delta2, delta
  occupation.csv           # rank, branch, σ, σ/σ1, delta, E_j, occupation_fraction
  residual_svd.csv         # (I-Π0)B_all の全特異値 s_k, N s_k, N² s_k
  q_svd.csv                # Q=[B0|B1] の特異値 q1..q4, Nq3,Nq4,N²q3,N²q4
  branch_tracking.csv      # 隣接時刻 平面重なり, source/target rank・branch, ambiguity_flag
  diagnostics_timeseries.csv  # 全時刻の数値誤差 + f + f_ref_dev
binary/N00005/, binary/N00040/
  timeseries.npz           # 全時刻の σ, ratio, delta, E, s_delta, q（object配列, 丸めなし）
  planes_<label>_step<t>.npz  # 代表6時刻の全平面基底 B_stack, B0, 固有値, Z, s_delta, q
tables/N00005/, tables/N00040/
  fulltable_<label>_step<t>.csv  # 代表時刻の全モード一行表（省略なし）
figures/N00005/, figures/N00040/
  fig1A_ratio_linear, fig1B_ratio_log, fig1C_ratio_heatmap, fig1D_branch_ratio,
  fig2A_N_ratio, fig2B_N2_ratio, fig3_delta_{linear,log},
  fig4_sdelta_{linear,log,N,N2}, fig5_q_svd, fig5_q34_zoom,
  fig6_occ_{linear,log}, fig7_mode_correspondence, fig8_representative_spectra
figures/comparison/
  comparison_N5_N40.png    # N=5,40 の残差最大特異値・非支配最大δ の重ね描き（生データ・無次元）
diagnostics/
  N00005.json, N00040.json # 数値診断の集約（最大誤差・crossing・代表時刻）
observation_report_no_interpretation.md  # 観測値のみ（解釈なし）
```

## 再現とデータの所在

`python3 code/run_exact_lowN_eigenspectrum_v1.py 5`（N=40 は `40`）→ `code/make_figures_exact_lowN_v1.py 5`
（`40`）→ `... comparison 5 40`。生データ CSV から全図を再生成できる。

**raw/ と tables/ の CSV（計約40MB）はリポジトリ規約 `*.csv` により git 追跡外だが、ディスク（Drive）上に
全て存在する。** 全精度の時系列データは binary/`timeseries.npz`（全時刻の σ,ratio,delta,E,s_delta,q）と
`planes_*.npz`（代表6時刻の全平面基底・固有値・状態）に丸めなしで格納され、git で追跡される。
CSV は上記コマンドでいつでも再生成できる。tables/ の代表全モード表は `git add -f` で明示追跡する。

## 数値検証（各時刻保存）

固有対残差 |Kv−λv|/(|K||v|)、実基底直交誤差 |BᵀB−I₂|、射影冪等 |Π²−Π|、固有平面間直交 |ΠⱼΠₖ|、
状態分解閉鎖 ||Z|²−ΣEⱼ−E_ker|、反対称 |K+Kᵀ|。集約は diagnostics/N*.json。
