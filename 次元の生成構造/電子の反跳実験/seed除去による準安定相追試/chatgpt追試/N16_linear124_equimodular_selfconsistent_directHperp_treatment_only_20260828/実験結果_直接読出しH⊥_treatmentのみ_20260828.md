# 実験結果——等モジュラー自己無撞着親からの seedless 走行、H⊥ を直交成分で直接読出し、treatment のみ（N=16）

**作成日:** 2026-08-28　**ベース:** `N16_linear124_decisions_applied_make_parent_equimodular_selfconsistent_20260828`（エンジン無変更、diff なし）。
**変更（実行スクリプトのみ、`results/diff_run_vs_equimodular_selfconsistent.patch`、3 hunk）:**
1. H⊥ の読出しを差し引き `htot − hpar`（丸め床 10⁻¹⁵）から、直交成分の直接計算 Z⊥ = Z − p(p·Z) − q(q·Z)、H⊥ = ‖Z⊥‖²（床 10⁻³² 級、非負）に変更。力学には影響しない。
2. baseline（位相のみ K）枝を実行しない・描かない（ノイズにしかならないため）。図の切詰めを 1e-30 → 1e-40 に。
親（make_parent）・相互作用・回転・step 数・seed は無変更。

## 結果

```
N=16: 親残差 9.4e-13、H_total 3.6471、onset(5%) None、max H⊥/H 1.577e-15（step 40000）、final PR/M 1.0000
  H⊥/H_total: step 0: 8.2e-32, step 1: 2.3e-27, step 10: 2.3e-25, step 100: 2.4e-23, step 1000: 2.1e-21, step 2000: 6.8e-21, step 5000: 3.8e-20, step 10000: 2.8e-19, step 20000: 6.2e-18, step 30000: 1.0e-16, step 40000: 1.6e-15
  fit 1e-28<f<1e-3 全域 (n=39900): slope 0.000320/step, R² 0.975247
  fit 1e-18<f<1e-14 後期 (n=26139): slope 0.000279/step, R² 0.999843
```

- **N=16 も動いている**：差し引き読出しでは 40000 step「不動」に見えたが、直接読出しでは 10⁻³² → 10⁻²⁷（step 1）→ 10⁻²¹（step 1000）→ 1.6×10⁻¹⁵（step 40000）と単調に成長している。前半は step の冪的（丸め変位の弾道的伝播）、後半（10⁻¹⁸〜10⁻¹⁴）は指数的で率 ≈ 0.0003/step——N=5 の 0.0090 の約 1/28。
- 40000 step では O(1) に達しない（10⁻¹⁵ から 10⁻¹ までさらに ~10⁵ step 必要）。「N=16 は安定」という前回の記述は**読出しの床による誤読**で、正しくは「不安定だが成長率が小さい」。飽和と局在を見るには step 数を 10 倍以上に延ばす必要がある（未実施）。
## 図（`figures/`、ファイル名はスクリプト固定で `N5_`・`baseline_vs` を含むが、中身は N=16 の treatment のみ）

- `N5_Hperp_baseline_vs_amplitude_aware.png`：H⊥（対数）vs step
- `N5_PR_baseline_vs_amplitude_aware.png`：PR/M
- `N5_amplitude_std_compare.png`：振幅分散
- `N5_closure_residual_compare.png`：|ZᵀZ|

## ファイル
`program/`、`data/`（treatment 時系列 csv・states npz・summary.json・key_steps.csv）、`figures/`、`results/`（diff・run.log・run_progress.log）、`run_all.sh`、`SHA256SUMS.txt`
