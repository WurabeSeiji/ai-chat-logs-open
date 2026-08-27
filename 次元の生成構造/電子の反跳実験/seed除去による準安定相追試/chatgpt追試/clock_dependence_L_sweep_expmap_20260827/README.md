# 時計依存性の検証（L 掃引・指数写像・秩序化・体積収縮）再現パッケージ

作成 2026-08-27〜28。「144 はサンプリング格子」という指摘を受け、γ = tan(π/L) を L = 144, 288, 576, 1152 で変え、自然時間 s = 2γ·step を揃えて、何が時計に依存し何が依存しないかを分離した。エンジンは論文 v6（DOI 10.5281/zenodo.22112009）と同一・無変更。

## 結論（詳細は `検討結果_時計依存性_L掃引と指数写像_20260827.md`）

- 時計に依存しない：onset 時刻、成長率 λ₁、σ_final = N−1、中立 Floquet 角の整数倍性、最終状態。
- 時計に依存する：step 単位の周期、体積収縮 Σln|μ|（自然時間あたり ∝ γ、L 倍で半減）、秩序化（等分配・位相クラス）の速さ（∝ γ）。
- 指数写像 e^{2γK} でも同じ L 依存 → 原因は Cayley の warping ではなく状態依存生成子の 1 step 凍結（一次分割誤差）。
- 論文 v6 §18–19・§26.2 の「秩序化／体積収縮」は数値散逸として読み直す必要がある。

## ファイル

| ファイル | 内容 |
|---|---|
| `sweep_L_and_expmap.py` → `results/sweep_L_and_expmap.{json,log}` | L 掃引（N=5,8）＋指数写像（N=5、L=144,288）：σ_eff²、onset、成長率、N=5 のクラス構造と δ |
| `ordering_vs_L.py` → `results/ordering_vs_L.{json,log}` | s_max=900 で σ² ロック・等分配・（N=5）クラス誤差の到達時刻を自然時間で記録 |
| `floquet_vs_L.py` → `results/floquet_vs_L.{json,log}` | 親まわりの rotating-frame Jacobian を L ごと・指数写像で：λ₁/s、Σln|μ|/s、中立角 |
| `run_n_scaling_lowrank_v1_no_sigma_norm.py` | エンジン原本（`../complex_simplex_decompactification_N5_N16_20260826.zip` 収録と同一） |
| `run_all.sh` | 全再現（約 4 分、python3 + numpy） |
| `SHA256SUMS.txt` | 全ファイルのハッシュ |

注：`ordering_vs_L.py` の N=5 クラス誤差指標（固定クラス割当）は親状態で既に 0 になるため判別力がなく、分析では σ² ロックと等分配を用いた。

## 再現

```bash
cd clock_dependence_L_sweep_expmap_20260827 && bash run_all.sh
```
