# 20260717 予備実験環境

## 概要

20260717 CHATGPT思考実験 (Section 18, 27) の検証を目的とした予備実験環境。

20260715 の実験結果から出現した二つのピーク (R_128≈0.688364, R_137≈0.697178) が、
単なる「倍音の共鳴」ではなく、「散乱系の共通安定点」であることを検証する。

---

## 実験構成

### 1. verify_section18_invariants_v1.py

**目的**: Section 18.1 の「一周期散乱行列の不変量比較」

散乱行列 M_c(R) の不変量（Trace, Determinant, 固有値の大きさ）が、
複数の調和条件でも同じ二つの R 値で一致するかを調べる。

**実行例**:
```bash
python3 verify_section18_invariants_v1.py \
  --output-dir results/section18_invariants \
  --num-steps 512 \
  --r-min 0.686 \
  --r-max 0.700 \
  --r-points 141
```

**出力**:
- `section18_invariants_detailed_v1.csv`: 各条件での Tr(M), det(M), |λ_±|
- `section18_common_points_v1.csv`: R 値ごとの不変量の consensus 判定
- `section18_summary_v1.json`: 統計サマリー

**解釈ポイント**:
- `trace_consensus` または `mag_consensus` が true のR値 → 複数条件で不変量が一致
- この一致が R_128 と R_137 で起こるなら、共通安定点仮説が強まる

---

### 2. verify_section18_harmonic_dependence_v1.py

**目的**: Section 18.2-18.4 の検証

- 18.2: 反復回数依存性（構造は初期から存在するか、それとも反復で形成されるか）
- 18.3: N=1（追加倍音なし）vs N=2, 3, ... の比較
- 18.4: 内部状態の同一性検査（二つのピークは別アトラクターか、同じ状態の別観測窓か）

**実行例**:
```bash
python3 verify_section18_harmonic_dependence_v1.py \
  --output-dir results/section18_harmonic
```

**出力**:
- `section18_iteration_dependency_v1.csv`: 反復回数と gray_depth の関係
- `section18_harmonic_comparison_v1.csv`: N=1, N=2, ... での gray_error 比較
- `section18_state_comparison_v1.json`: R_128 と R_137 での終状態の比較

**解釈ポイント**:
- **18.2**: 初期ステップから gray_depth が大きいなら、構造は初期から存在（v1 仮説）
- **18.3**: N=1 では peak が出ず、N≥2 でのみ出るなら、追加倍音が必須
- **18.4**: `Identical_internal_states` なら同じアトラクターを二つの R が支える

---

### 3. verify_section27_paper7_connection_v1.py

**目的**: Section 27 の「Paper7 との接続」を検証

閉塞条件 A²+B²=0 が、Paper7 の R=3 4次元完全包含セル問題と対応するかを確認。

**実行例**:
```bash
python3 verify_section27_paper7_connection_v1.py \
  --output-dir results/section27_paper7 \
  --r-min 0.685 \
  --r-max 0.702 \
  --r-points 171 \
  --steps 512
```

**出力**:
- `section27_closure_residual_sweep_v1.csv`: 各 R での閉塞残差 E_close
- `section27_paper7_alignment_v1.json`: Paper7 ターゲットとの alignment 判定
- `section27_radius_analysis_v1.json`: 4D 状態空間の半径統計

**解釈ポイント**:
- `E_close = |A²+B²| / (|A|²+|B|²)` が R_128, R_137 でのみ小さい
  → 閉塞条件の自動満足が二つの特定R値に限定される
- `closure_min` が Paper7 の 137-セル構造に対応する場合、
  → R=3 の完全包含条件が動力学として実現されている可能性

---

## 予期される結果と解釈

### 成功シナリオ

```
1. verify_section18_invariants_v1.py:
   → R_128, R_137 で trace_consensus = True
   → 複数調和条件でも Tr(M) が一致
   
2. verify_section18_harmonic_dependence_v1.py:
   → Section 18.2: gray_depth が初期から大きい
   → Section 18.3: N=1 では peak なし、N≥2 で peak 出現
   → Section 18.4: interpretation = "Identical_internal_states"
   
3. verify_section27_paper7_connection_v1.py:
   → closure_min が R_128, R_137 で最小化
   → 他の R では E_close が大きい
   → paper7_alignment = "BOTH_PEAKS_ALIGNED"
```

### 理論的意義

全スクリプトでこのシナリオが確認されれば：

> **二つの R 値は、閉塞公理 A²+B²=0 から必然的に現れる、
> 散乱系の共通な安定構造であり、
> Paper7 の 4次元完全包含セル問題の動的実現である。**

という仮説が大きく強化される。

---

## 環境要件

- Python 3.8+
- numpy, scipy

```bash
pip install numpy scipy
```

---

## 実行手順（推奨順序）

### Phase 1: 基本検証（~30分）

```bash
# 1. 不変量検査（最優先）
python3 verify_section18_invariants_v1.py --output-dir results/phase1

# 2. 調和依存性
python3 verify_section18_harmonic_dependence_v1.py --output-dir results/phase1

# 3. 結果確認
ls -la results/phase1/
cat results/phase1/section18_summary_v1.json
```

### Phase 2: Paper7 接続検証（~20分）

```bash
# 閉塞残差スイープ
python3 verify_section27_paper7_connection_v1.py \
  --output-dir results/phase2 \
  --r-points 171
  
# 結果確認
cat results/phase2/section27_summary_v1.json
```

### Phase 3: 詳細分析（オプション）

CSV ファイルを Excel/Python で分析：
- 二つの R ピークの特異性を可視化
- closure_residual プロット
- 不変量の調和依存性を検討

---

## トラブルシューティング

### ImportError: numpy/scipy

```bash
pip install --upgrade numpy scipy
```

### MemoryError（大規模スイープ）

`--r-points` を減らす:
```bash
python3 verify_section27_paper7_connection_v1.py \
  --r-points 85  # 50 に削減
```

### 結果が N/A か異常値

- `--steps` を増やしてみる（デフォルト 512）
- 初期条件 `s0` の値を確認（コードで 0.01 固定）

---

## 次ステップ（予定）

1. **外部仮定なしで R を導出**
   - Section 18.1 の Tr(M), det(M) から R=3 を逆算できるか
   - 散乱行列の安定性分析

2. **流域図（Basin of Attraction）**
   - 初期 (φ, S₀) 平面での流域分離を可視化

3. **有限時間リアプノフ指数**
   - λ_T の計算で吸引性を定量化

---

## 参考文献

- 20260717 CHATGPT思考実験.md (Section 18, 27)
- 20260715 exchange_weight_alpha_correspondence_numerical_experiment_ja.md
- Paper7_Alpha_Identity_ja.md (Paper7 参照用)

---

作成日: 2026-07-17
