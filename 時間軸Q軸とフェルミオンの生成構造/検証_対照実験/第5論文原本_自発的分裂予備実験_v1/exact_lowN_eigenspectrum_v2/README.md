# exact_lowN_eigenspectrum_v2 — N=5 厳密固有値／固有方向 修正版（第1段階）

`exact_lowN_eigenspectrum_v1_修正指示.md` に基づく修正実装。**N=5 のみ**。N=40 は人間検収後、
N=300 には着手しない。解釈は書かない（observation_report_no_interpretation.md 参照）。

## v1 からの修正

1. **基本固有分解を `eigh(1j*K)`** に変更（H=iK は Hermitian、U は全体直交規格）。
   v1 の `eig(K)`＋実虚部独立直交化（平面間直交誤差 最大0.9）を主要法から除外。
2. **絶対閾値 σ>1e-9 廃止**。全 M 固有値を保存し、σ/(ε‖K‖) の相対基準（1000, 100 併記）で表示。
3. **近縮退の併合クラスタ化**：平面間重なり(作用素ノルム)>1e-12 の対を union-find 併合。
   非縮退（未併合）平面の相互直交を ≤1e-12 に保証。生固有値は平均せず全保存、方向比較は Π_C。
4. **状態占有は直交射影**（複素直交固有ベクトルの二乗和）、核は明示 Π_ker、閉鎖検証。
5. **新方向は対象別残差**（A 支配／C 占有／D 縮退／all 診断）。R_all 一括は診断のみ。

## 検収（§13, N=5）達成

非縮退平面間直交 ≤9.7×10⁻¹³（≤10⁻¹²）、固有対残差 7.5×10⁻¹⁵、正負対誤差 6.7×10⁻¹⁵、
KB=BJ ≤2.35×10⁻¹⁵、閉鎖 3.3×10⁻¹⁵、全固有値保存（閾値削除なし）、最終スペクトルは v1 と一致。

## 出力

```
code/run_exact_lowN_eigenspectrum_v2.py   # 観測本体（N=5のみ）
code/make_figures_exact_lowN_v2.py        # 図（§10 の13図）
raw/N00005/
  eigenvalues.csv     # 全 M 固有値/時刻: mu, sigma, sigma/sigma1, N/N², abs_sigma_over_eps_normK, pair_id, cluster_id, solver_residual
  clusters.csv        # クラスタ/時刻: sigma, mult, dim, occupation, delta_C, overlap_parent, invariance_residual, kbj, initial_floor_flag
  delta_targets.csv   # 対象別残差(A/C/D/all)の全特異値
  q_svd.csv           # Q=[B0|Bdom] の q1..q4
  cluster_tracking.csv, diagnostics_timeseries.csv
binary/N00005/  timeseries.npz, decomp_<label>_step<t>.npz（代表6時刻: mu,U,cluster_*,Z 等 丸めなし）
tables/N00005/  fulltable_<label>.csv（全固有値）, clusters_<label>.csv（全クラスタ）
figures/N00005/ fig01〜fig13
diagnostics/N00005.json
observation_report_no_interpretation.md
```

raw/ と tables/ の CSV はリポジトリ規約 `*.csv` により git 追跡外だがディスク（Drive）上に全存在。
全精度データは binary/ npz で追跡。`git add -f tables/` で代表表を明示追跡。

## 再現

`python3 code/run_exact_lowN_eigenspectrum_v2.py 5` → `python3 code/make_figures_exact_lowN_v2.py 5`。
