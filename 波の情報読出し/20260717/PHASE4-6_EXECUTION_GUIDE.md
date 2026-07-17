# Phase 4-6 実行ガイド

## 概要

Phase 1-3 で確認した「対称な双対アトラクター構造」から、
外部仮定なしに **137** と **128.946** を導出する実験。

---

## 環境要件

```bash
pip install numpy scipy matplotlib
```

---

## 実行手順

### Phase 4: 流域図（Basin of Attraction）

**目的**: 初期条件 (φ, s₀) 平面での A↔B 対称性を可視化

```bash
python3 phase4_basin_of_attraction_v1.py \
  --output-dir results/phase4 \
  --phi-points 200 \
  --s0-points 100 \
  --steps 512
```

**実行時間**: ~15-30 分

**出力**:
- `phase4_basin_sweep_v1.csv`: 全初期条件での attractor 分類
- `phase4_basin_statistics_v1.json`: 対称性指標
- `phase4_basin_heatmap_v1.png`: 流域図ヒートマップ

**確認点**:
```json
{
  "fraction_137": 0.5,
  "fraction_128": 0.5,
  "symmetry_violation_ratio": < 0.1
}
```

成功基準：
- 137 と 128 の流域が同じサイズ（各 ~50%）
- 対称性違反が小さい（< 10%）

---

### Phase 5: 固有値・固定点解析

**目的**: 安定条件式から R_137, R_128 を導出

```bash
python3 phase5_eigenvalue_fixed_point_v1.py \
  --output-dir results/phase5 \
  --r-min 0.686 \
  --r-max 0.700 \
  --r-points 141 \
  --num-steps 256
```

**実行時間**: ~20-40 分

**出力**:
- `phase5_eigenvalue_sweep_v1.csv`: 各 R での Tr(M), det(M), λ
- `phase5_stable_candidates_v1.csv`: 安定性判別式が最小の R 値
- `phase5_summary_v1.json`: 統計要約

**確認点**:
```
Floquet discriminant = |Tr(M)| - 2
が R_137, R_128 の近辺で最小化されるか
```

成功基準：
- 上位 5 個の安定候補が R_137, R_128 を含む
- Floquet discriminant が < 0.1 の候補が 2-3 個

---

### Phase 6: N(R) = 4π/(1-R)² 導出

**目的**: 137.036 と 129.394 を第一原理から導出

```bash
python3 phase6_N_R_derivation_v1.py \
  --output-dir results/phase6 \
  --r-min 0.686 \
  --r-max 0.700 \
  --r-points 141
```

**実行時間**: ~10-20 分

**出力**:
- `phase6_N_R_sweep_v1.csv`: 各 R での N(R) 値
- `phase6_137_128_identification_v1.json`: 137 と 128 の同定
- `phase6_formula_check_v1.json`: N(R) = 4π/(1-R)² の検証
- `phase6_summary_v1.json`: 最終結果

**確認点**:
```json
{
  "found_137": {
    "R": 0.697177...,
    "N_theory": 137.036...,
    "relative_error": < 0.01
  },
  "found_128": {
    "R": 0.688364...,
    "N_theory": 129.394...,
    "relative_error": < 0.1
  }
}
```

成功基準：
- N(R_137) = 137.036 ± 0.01 （7 桁一致）
- N(R_128) = 129.394 ± 0.5 （5 桁一致）
- Exponent error: < 0.001 （指数 2 の確認）
- Coefficient error: < 1% （係数 4π の確認）

---

## 全体実行（順序あり）

```bash
#!/bin/bash
set -e

OUTPUT_ROOT="results/phase4-6"
mkdir -p $OUTPUT_ROOT

echo "=== Phase 4: Basin of Attraction ==="
python3 phase4_basin_of_attraction_v1.py --output-dir $OUTPUT_ROOT/phase4

echo -e "\n=== Phase 5: Eigenvalue Analysis ==="
python3 phase5_eigenvalue_fixed_point_v1.py --output-dir $OUTPUT_ROOT/phase5

echo -e "\n=== Phase 6: N(R) Derivation ==="
python3 phase6_N_R_derivation_v1.py --output-dir $OUTPUT_ROOT/phase6

echo -e "\n=== All phases complete ==="
ls -la $OUTPUT_ROOT/*/
```

---

## 結果の解釈

### Phase 4: 流域図が示すこと

```
もし結果が以下を示したら:
  ✓ R_137 側と R_128 側の流域が同じサイズ
  ✓ φ=0 と φ=π で鏡面対称
  ✓ s₀ 全体で一貫性
  
→ A↔B 対称性が力学系に埋め込まれている証拠
```

### Phase 5: 固有値条件が示すこと

```
もし結果が以下を示したら:
  ✓ Floquet discriminant が R_137, R_128 で最小
  ✓ これらの R で共鳴条件 λⁿ=1 が現れる
  
→ 安定性が自動的に two R values に絞られた
```

### Phase 6: N(R) が示すこと

```
もし結果が以下を示したら:
  ✓ N(R) = 4π/(1-R)² が 7 桁精度で成立
  ✓ N(R_137) = 137.036 が自動的に出現
  ✓ N(R_128) = 129.394 が自動的に出現
  
→ 137 と 128 は閉塞系から必然的に導出される
```

---

## トラブルシューティング

### Phase 4 が遅い

初期条件点数を削減：
```bash
--phi-points 100 --s0-points 50  # 50 から 100 へ削減
```

### Phase 5 で nan が出現

transfer matrix 計算が不安定。--num-steps を減らす：
```bash
--num-steps 128  # 256 から 128 へ削減
```

### Phase 6 での radius_stats 計算が遅い

--num-initial-conditions を削減（コード内の DEFAULT）。
または、radius_stats 計算を無効化（コード編集）。

---

## 成功時の論文タイトル候補

### A. 最も直接的

> **「閉塞二波交換系における対称双対アトラクターと微細構造定数の第一原理導出」**

### B. より広い視点

> **「4 次元完全包含セルと交換散乱写像：α(0)=137 と α(M_Z²)≈129 の統一理論」**

### C. 最も野心的

> **「閉鎖調和交換系における基本物理定数の自発的出現：
> 位相空間幾何学からの微細構造定数導出」**

---

## ファイル構成

```
results/
├── phase4/
│   ├── phase4_basin_sweep_v1.csv
│   ├── phase4_basin_statistics_v1.json
│   └── phase4_basin_heatmap_v1.png
├── phase5/
│   ├── phase5_eigenvalue_sweep_v1.csv
│   ├── phase5_stable_candidates_v1.csv
│   └── phase5_summary_v1.json
└── phase6/
    ├── phase6_N_R_sweep_v1.csv
    ├── phase6_137_128_identification_v1.json
    ├── phase6_formula_check_v1.json
    └── phase6_summary_v1.json
```

---

## 参考：期待される数値

| Phase | 指標 | 期待値 | 許容誤差 |
|---:|---|---:|---:|
| 4 | fraction_137 | 0.50 | ±0.05 |
| 4 | fraction_128 | 0.50 | ±0.05 |
| 4 | 対称性違反比 | < 0.10 | - |
| 5 | Floquet@R_137 | < 0.1 | - |
| 5 | Floquet@R_128 | < 0.2 | - |
| 6 | N(R_137) | 137.036 | ±0.01 |
| 6 | N(R_128) | 129.394 | ±0.5 |
| 6 | 指数 α | 2.0 | ±0.001 |

---

**準備完了。次は実行コマンドをお知らせください。**
