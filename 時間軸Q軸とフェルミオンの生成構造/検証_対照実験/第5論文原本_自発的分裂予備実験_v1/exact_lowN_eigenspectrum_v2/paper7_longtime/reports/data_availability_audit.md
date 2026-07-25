# データ監査（フェーズ1）— 論文7 再計算・再図化・横安定性検査

指示書「Claude_Code_論文7_再計算_再図化_横安定性検査_厳密指示書.md」§15フェーズ1。事実のみ列挙。図は捏造しない。

## 1. 論文6 図2・図3（平面流入）の元コード・定義

- 元コード：`run_plane_flow_exact_v1.py`（N=5,40, 密行列 eig）／`run_plane_flow_approx_v1.py`（N=300, 低ランク JG, σ_rel_threshold=1e-6）。
- 3色の定義（固定親基底 B_p1, B_rot を親生成子 K(arg v) から一度だけ確定し、時間発展する Z を射影）：
  1. `frac_P1`：支配平面 P1（親支配平面, 2次元）への占有比
  2. `frac_other_rotation`：その他回転部分空間（0<σ<σ1）への占有比
  3. `frac_kernel`：核（σ=0）への占有比
- 固定親基底次元（`plane_flow_result_v1/*.json` より）：

| N | M | P1 | その他回転 | 核 | crossing |
|---:|---:|---:|---:|---:|---:|
| 5 | 10 | 2 | 4 | 4 | 1167 |
| 40 | 780 | 2 | 78 | 700 | 2011 |
| 300 | 44850 | 2 | 598 | 44250 | 4844 |

## 2. 論文6 図1（黒線 = 分裂量 f）

- `metastable_series_result_v1/fcurve_N000{05,40,300}_delta1e-15_seed0.csv`（列 tau, f）。
- 恒等式 `f = 1 - frac_P1`（論文6 §5.3）。

## 3. 既存の時系列の到達 step

| データ | N=5 | N=40 | N=300 |
|:--|--:|--:|--:|
| 論文6 planeflow CSV（after=20000） | ~crossing+20000 | ~crossing+20000 | ~crossing+20000 |
| v2 saturation q_svd（crossing+50000） | 51167 | 52011 | 54844 |

## 4. 今回の確定共通横軸

- **絶対 step 0 ≤ τ ≤ 55000**（全N共通、目盛り5000刻み、crossing不動）。
- 既存データはいずれも 55000 絶対 step に未達 → **全N（5,40,300）を 0〜55000 で再計算する**。

## 5. 不足項目（再計算・新規作成が必要）

1. **0〜55000 の占有時系列（全N）**：論文6の固定親基底3色（frac_P1, frac_other_rotation, frac_kernel）＋黒線 f を、共通横軸55000まで再取得。→ 再計算。
2. **5色分解（新規）**：確定指示に基づき「その他回転部分空間」を「新方向3・新方向4・残余その他回転」の3色へ分解（P1と核は保持、計5色）。新方向3・4 = 時間依存 S₄(t)=正規直交化[B₀|B_dom(t)] のうち B₀(=P1)直交補の2方向を、固定 other 空間へ射影して構成。縮退時は2次元平面として追跡し表示のみ連続基底固定。固有値順の色割当はしない。→ 新規計算。
3. **横安定性時系列（新規）**：時間依存 S₄(t)=[B₀|B_dom(t)] の S₄(t)^⊥ へ複数seed・複数eps（1e-8/1e-10/1e-12, 可能なら1e-14）の微小摂動、Benettin再正規化で λ⊥,max^(4)。→ 新規計算。
4. **図（新規）**：図1/2/3（各N＋比較3段）、横摂動増幅率、横成長率時系列、λ vs N / 1/N / 1/logN。PNG+SVG。→ 新規作成。

## 6. 変更しない条件（指示書§1.1, §13）

seed（40260722+1000N）、初期微小種 δ=1e-15、零閉鎖、正規化、Cayley更新、演算順序、crossing定義（f>0.05）、論文6の物理量・分類定義（P1/other/kernel）、黒線の意味（f=1-E_P1）。

## 7. 使用する既存コード（不変更 import）

- 原本エンジン `run_n_scaling_lowrank_v1.py`（SHA固定）
- 固定親基底：`run_plane_flow_exact_v1.parent_plane_split_exact`（N=5,40）／`run_plane_flow_approx_v1.parent_plane_split_approx`（N=300）
- 支配平面（B_dom, gram法）：`run_n300_dimension_saturation_v2` の gram_reduce / dominant_plane（N=40で密行列と一致検証済）
