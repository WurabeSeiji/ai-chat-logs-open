# 実験結果——等モジュラー自己無撞着親からの seedless 走行、H⊥ を直交成分で直接読出し、treatment のみ（N=5）

**作成日:** 2026-08-28　**ベース:** `N5_linear124_decisions_applied_make_parent_equimodular_selfconsistent_20260828`（エンジン無変更、diff なし）。
**変更（実行スクリプトのみ、`results/diff_run_vs_equimodular_selfconsistent.patch`、3 hunk）:**
1. H⊥ の読出しを差し引き `htot − hpar`（丸め床 10⁻¹⁵）から、直交成分の直接計算 Z⊥ = Z − p(p·Z) − q(q·Z)、H⊥ = ‖Z⊥‖²（床 10⁻³² 級、非負）に変更。力学には影響しない。
2. baseline（位相のみ K）枝を実行しない・描かない（ノイズにしかならないため）。図の切詰めを 1e-30 → 1e-40 に。
親（make_parent）・相互作用・回転・step 数・seed は無変更。

## 結果

```
N=5: 親残差 3.8e-11、H_total 0.8596、onset(5%) 4708、max H⊥/H 9.888e-01（step 7474）、final PR/M 0.2728
  H⊥/H_total: step 0: 9.6e-32, step 1: 3.7e-24, step 10: 3.7e-22, step 100: 4.0e-20, step 1000: 2.4e-16, step 2000: 1.8e-12, step 5000: 4.2e-01, step 10000: 8.6e-01, step 20000: 8.8e-01, step 30000: 8.5e-01, step 40000: 7.7e-01
  指数域 fit（1e-10<H⊥<1e-3）: slope 0.00890/step, R² 1.000000000, window 2470–4280
  fit 1e-28<f<1e-3 全域 (n=4163): slope 0.008898/step, R² 0.999961
  fit 1e-18<f<1e-14 後期 (n=1052): slope 0.008749/step, R² 0.999843
```

- 直接読出しにより潜伏期の中身が見える：step 0 の 10⁻³²（丸め）から step 1 で 10⁻²⁴ に跳び（固定点からの最初の丸め変位）、以後 **一定率の指数成長 0.0090/step が 10⁻²⁴ から 10⁻¹ まで 23 桁一直線**（全域 fit と 1e-10〜1e-3 の fit が一致、R² ≈ 1）。差し引き読出しで「助走」に見えていた区間は、実際には同じ率で成長していた。
- 約 5000 step で飽和（H⊥/H 0.99）、以後 PR/M ≈ 0.27 の局在状態。
## 図（`figures/`、ファイル名はスクリプト固定で `N5_`・`baseline_vs` を含むが、中身は N=5 の treatment のみ）

- `N5_Hperp_baseline_vs_amplitude_aware.png`：H⊥（対数）vs step
- `N5_PR_baseline_vs_amplitude_aware.png`：PR/M
- `N5_amplitude_std_compare.png`：振幅分散
- `N5_closure_residual_compare.png`：|ZᵀZ|

## ファイル
`program/`、`data/`（treatment 時系列 csv・states npz・summary.json・key_steps.csv）、`figures/`、`results/`（diff・run.log・run_progress.log）、`run_all.sh`、`SHA256SUMS.txt`
