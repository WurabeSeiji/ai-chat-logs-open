# v2補完実験_4生成法_N3toN16_統一プロトコル_20260830

論文 v1（DOI 10.5281/zenodo.22112009）の v2 に必要な不足データを、**統一プロトコル**（振幅込み K_ij=Im(z̄ᵢzⱼ)、厳密線形回転 exp((2π/124)K)、正規化なし、種なし Z₀=v、H⊥ 直交成分の直接読出し、40000 step）で N=3〜16 × 4 生成法（54 親）について取得したもの。予測（共回転モノドロミー ρ）は走行前に固定。

- 生成法：`mp` make_parent 等モジュラー（3 段階、`original_engine.make_parent`、rng 40260721+1000N）／`hm` 手作り等モジュラー（基礎文書構成：偶数 1-因子分解、奇数 距離クラス、N=3 Z3）／`ne` 非等モジュラー（q≥4 クラス重み付き族 ε=0.6、q≤3 多様体上の代表点）／`rb` 乱数均衡親 {S=0, W 均等}（N≥5）。各 N で ‖v‖ を mp 親に揃える（`results/norm_by_N.json`）。
- `program/original_engine.py`：N8 テンプレートの byte コピー（sha fa68d344…）。`run_dynamics.py`：手作り親パッケージ pass2_run.py の 3 行変更コピー。`common.py`／`state_provider.py`：重心閉塞パッケージ／非等モジュラー版基礎文書からのコピー。
- `program/pass1_parents.py` → `data/<tag>/parent_v.{npz,csv}`, `parent_checks.json`；`results/parents_predictions.csv`（走行前予測）、`closure_step0_4methods.csv`（主張 1）
- `program/pass2_embed_random.py` → `results/embed_random{,_summary}.csv`（主張 2：乱数状態 100×14 の複素埋め込み vs 実の半正定値性）
- `program/pass5_analysis.py` → `results/dynamics_summary.csv`, `matrix_N_by_method.md`
- `program/pass6_figures.py` → `figures/fig1〜fig7`
- `data/reference/closure_step0_4systems_20260829.csv`：`../公理見直し_ゼロ閉塞定理と固有時計_20260829/results/closure_step0.csv` のコピー（fig1 の参照点）
- `results/run_dynamics_{A,B}.log`、`pass*.log`
- 走行 CSV は gzip 保存、全 step 状態 `states_treatment.npz` は git 除外（`run_all.sh` で再生成）。

結果の読み方は `実験結果_v2補完_4生成法_N3toN16_20260830.md`。
