# 第1予備実験 二段階seed除去対照 報告（観測値のみ・解釈なし）

指示書「Claude_Code_第8論文_環境構築_第1予備実験_二段階seed除去対照.md」。数値と図のみ。物理解釈・
一般化・名称付けは行わない。FFT・自己相関・Lyapunov・seed掃引・ε掃引・スケール不変性判定は未実施。

## 1. 実行環境

- Python 3、NumPy 2.0.2、SciPy 1.13.1、Matplotlib 3.9.4。
- 決定論: 乱数 `numpy.random.default_rng(40260722 + 1000*N)`。初期微小種 δ=1e-15。Cayley γ=tan(π/144)。
- 共通最終 step COMMON_FINAL_STEP = 55000（絶対 step、crossing 不動、目盛り5000刻み）。
- 記録間隔: N=5,40 は 25 step、N=300 は 100 step。
- 実行時間（概算）: N=5 各条件 ~4s、N=40 各条件 ~19s、N=300 各条件 ~13分（A/B/D 合計 ~40分）。

## 2. 第7論文から再利用したコードと設定（read-only import・不変更）

- 原本エンジン `run_n_scaling_lowrank_v1.py`（SHA `ba0fc19b…`）
- 固定親基底 `run_plane_flow_exact_v1.py`（N=5,40）/ `run_plane_flow_approx_v1.py`（N=300）
- 支配平面 Gram 縮約 `run_n300_dimension_saturation_v2.py`
- 横摂動 η_⊥ 生成・S₄基底・射影 `run_paper7_transverse.py`（seed index 0 の初回生成のみ流用、Benettin 不使用）

SHA-256 は `config/source_file_hashes.json`。第7論文フォルダは書き換えていない。seed の ON/OFF は
第8論文ラッパー `run_preliminary_seed_ablation_v1.py` で明示切替。

## 3. 条件A・B・Dの定義

- 条件A: initial seed OFF / metastable seed OFF。Z₀ = v（kernel seed g を生成せず乱数を消費しない）。
- 条件B: initial seed ON / metastable seed OFF。Z₀ = (v+δg)/‖·‖, δ=1e-15。以後自然発展のみ。
- 条件D: initial seed ON / metastable seed ON。B と t₁=crossing+3000 直前までビット一致。t₁ で
  単一横摂動 Z_D(t₁)=(Z_B(t₁)+ε η_⊥)/‖·‖（ε=1e-8, transverse seed index=0）を一回注入し、
  規格化一回の後、以後 B と同一の自然発展。Benettin・再注入・再正規化なし。既存 transverse CSV は参照専用。

## 4. N=5 の観測結果

| 条件 | crossing | 準安定開始 | max_f | final_f | final_q3 | final_q4 | max/final rank_Q | mean_meta(f) | std_meta(f) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| A | 1166 | 4166 | 0.9655 | 0.8754 | 0.9315 | 0.7202 | 4 / 4 | 0.8744 | 1.29×10⁻² |
| B | 1167 | 4167 | 0.9600 | 0.8053 | 0.8140 | 0.6893 | 4 / 4 | 0.8054 | 2.87×10⁻² |
| D | 1167 | 4167 | 0.9600 | 0.8026 | 0.8110 | 0.6878 | 4 / 4 | 0.8034 | 2.88×10⁻² |

- B/D は t₁=4167 直前まで最大差 0.00×10⁰（ビット一致）。
- 条件A（seed無し）で crossing 検出、rank_Q 2→4、q₃,q₄ が数値床を超えた。

## 5. N=40 の観測結果

| 条件 | crossing | 準安定開始 | max_f | final_f | final_q3 | final_q4 | max/final rank_Q | mean_meta(f) | std_meta(f) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| A | 2011 | 5011 | 0.1938 | 0.1938 | 0.3353 | 0.3031 | 4 / 4 | 0.1934 | 7.38×10⁻⁴ |
| B | 2011 | 5011 | 0.2035 | 0.2035 | 0.3396 | 0.3159 | 4 / 4 | 0.2030 | 1.29×10⁻³ |
| D | 2011 | 5011 | 0.2035 | 0.2035 | 0.3396 | 0.3159 | 4 / 4 | 0.2030 | 1.29×10⁻³ |

- B/D は t₁=5011 直前まで最大差 0.00×10⁰（ビット一致）。final_f, final_q₃, final_q₄, mean/std_meta が
  B と D で一致。
- 条件A（seed無し）で crossing 検出、rank_Q 2→4。

## 6. N=300 の観測結果

| 条件 | crossing | 準安定開始 | max_f | final_f | final_q3 | final_q4 | max/final rank_Q | mean_meta(f) | std_meta(f) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| A | 4849 | 7849 | 0.0898 | 0.0898 | 0.2174 | 0.2114 | 4 / 4 | 0.0882 | 3.53×10⁻³ |
| B | 4844 | 7844 | 0.0861 | 0.0861 | 0.2123 | 0.2074 | 4 / 4 | 0.0848 | 2.94×10⁻³ |
| D | 4844 | 7844 | 0.0861 | 0.0861 | 0.2123 | 0.2074 | 4 / 4 | 0.0848 | 2.94×10⁻³ |

- B/D は t₁=7844 直前まで最大差 0.00×10⁰（ビット一致）。final_f, final_q₃, final_q₄, mean/std_meta が
  B と D で一致（有効桁内）。
- 条件A（seed無し）で crossing 検出、rank_Q 2→4。B の crossing=4844 は第7論文 N=300 と一致。

## 7. 二つのseed除去による差分

### 7.1 初期seed除去（条件A, §12.1）

| 問い | N=5 | N=40 | N=300 |
|:--|:--|:--|:--|
| 幾何級数的発展が起きたか | 起きた | 起きた | 起きた |
| crossing が検出されたか | 検出（1166） | 検出（2011） | 検出（4849） |
| rank_Q が2から増えたか | 2→4 | 2→4 | 2→4 |
| q₃,q₄ が数値床を超えたか | 超えた | 超えた | 超えた |
| 準安定開始が検出されたか | 検出（4166） | 検出（5011） | 検出（7849） |

条件A の crossing（1166/2011/4849）は各N で条件B（1167/2011/4844）とほぼ一致。

### 7.2 準安定seed除去（条件B vs D, §12.2）

| 問い | N=5 | N=40 | N=300 |
|:--|:--|:--|:--|
| B でも準安定振動が継続したか | 継続（std_meta 2.87×10⁻²） | 継続（1.29×10⁻³） | 継続（2.94×10⁻³） |
| B と D で振動振幅（std_meta）が異なるか | 2.87 vs 2.88×10⁻²（差 ~10⁻⁴） | 一致 | 一致 |
| B と D で振動中心（mean_meta）が異なるか | 0.8054 vs 0.8034（差 ~2×10⁻³） | 一致 | 一致 |
| B と D で方向数（rank_Q）が異なるか | 同一（4） | 同一（4） | 同一（4） |
| 条件D だけに現れる応答があるか | final_f/q₃/q₄ に ~10⁻³ の差 | 差なし（B=D） | 差なし（B=D） |

## 8. 数値診断

全N・全条件で（記録時刻の最大値）:

| 量 | 最大値の範囲 |
|:--|:--|
| 規格化誤差 \|\|Z\|²−1\| | ≤ 1.4×10⁻¹⁴ |
| 零二乗閉鎖 \|ZᵀZ\| | ≤ 1.8×10⁻¹⁰（条件D N=5 の注入直後が最大、他は ≤1.3×10⁻¹³） |
| 射影閉鎖誤差 (E_P1+E_other+E_ker)/\|Z\|²−1 | 0.0 |
| 反対称誤差 \|K+Kᵀ\|（N≤40 で密行列確認, 終端） | 倍精度範囲（N=300 は密行列不可のため未算出） |
| v の規格化・零二乗（条件A の Z₀=v） | v_norm_error=0.0、v_zero_square_abs≈1.1×10⁻¹⁵ |

（診断 JSON: `diagnostics/N000{05,40,300}_condition_{A,B,D}.json`）

## 9. 未実施項目（指示書§12.2 により本予備実験では行わない）

FFT による周波数同定、自己相関時間、位相同期解析、Lyapunov 指数、seed 振幅依存則、seed 方向依存性、
スケール不変性判定、初期潜伏相と後期準安定相の自己相似性判定。ε 掃引・seed 掃引も未実施（条件D は
ε=1e-8, seed index=0 の固定一条件）。

## 10. 生成ファイル一覧

- コード: `code/audit_paper7_dependencies_v1.py`, `code/run_preliminary_seed_ablation_v1.py`,
  `code/make_preliminary_seed_ablation_figures_v1.py`
- 設定: `config/experiment_manifest.json`, `config/source_file_hashes.json`
- 生時系列（§9 全27列）: `raw/N000{05,40,300}/condition_{A_no_seed,B_initial_only,D_existing_two_seed}.csv`
- 集計: `summary/preliminary_seed_ablation_summary.csv`
- 図（§10, 各N 4種）: `figures/fig01_f_compare_*`, `fig02_q3q4_compare_*`, `fig03_rankQ_compare_*`,
  `fig04_metastable_B_vs_D_*`
- 診断: `diagnostics/N000{05,40,300}_condition_{A,B,D}.json`
- 監査: `reports/paper7_dependency_audit.md`

---

観測事実の要約（数値のみ）: 条件A（明示的初期seed無し, Z₀=v）は N=5,40,300 のいずれでも crossing が
検出され、rank_Q は 2 から 4 へ増加した。条件B（初期seedのみ）でも準安定観測量の時間変動（std_meta>0）は
最終step まで継続した。条件D（t₁で単一横摂動 ε=1e-8 を一回注入）は、N=40,300 では final_f/q₃/q₄ および
mean/std_meta が条件B と一致し、N=5 では ~10⁻³ の差であった。B と D は t₁ 直前まで全N でビット一致した。
