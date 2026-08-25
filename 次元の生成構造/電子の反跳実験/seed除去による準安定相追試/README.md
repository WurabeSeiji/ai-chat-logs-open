# seed除去による準安定相追試（論文8 Stage A2a 再現実験）

## 目的

第8論文 Stage A2a（完全無seed条件、N=5、倍精度浮動小数点）の再現対照実験。同じ結果が厳密に得られるか検証する。

## 実験設計

### 対照対象

- **原本環境**: `../../../次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper8_stage_A2a_seedless_N5`
- **追試環境**: `./original_environment`（原本の完全コピー）
- **新規実行**: `./reproduction_run`（このフォルダで独立実行）

### 検証項目

1. **ビット一致性**: 独立実行の $f(t)$、$q(t)$ がビット一致するか
2. **数値健全性**:
   - 規格化誤差 $\lVert Z\rVert^2$ 偏差
   - 零二乗閉鎖 $|Z^T Z|$ 誤差
   - 射影閉鎖誤差
   - 占有和誤差
3. **物理量の一致**:
   - crossing（$f > 0.05$ 初回到達）
   - 準安定開始（crossing + 3000）
   - 最大分裂フラクション
   - 方向部分空間の一致度

### 設定値（原本から厳密に継承）

- 親: `make_parent(sys, np.random.default_rng(40265722), iters=1200, tol=1e-12)`
- 初期状態: `Z0 = v.copy()`（明示的摂動なし）
- 実行範囲: step 0 ～ 5000
- PRNG seed: 固定（同じseedで複数実行）
- 更新則: Cayley変換（$\gamma = \tan(\pi/144)$）
- q測定: 5 stepごと
- 占有測定: 25 stepごと

## 実行手順

```bash
cd ./reproduction_run
python3 01_run_seedless_independent.py
python3 02_compare_with_original.py
python3 03_numerical_health_check.py
python3 04_generate_verification_report.py
```

各スクリプトの成功記録は `logs/` に保存される。

## 出力構成

```
./reproduction_run/
├── raw/
│   ├── run1_f_series.npy
│   ├── run1_q_series.npy
│   ├── run2_f_series.npy
│   ├── run2_q_series.npy
│   └── ...
├── processed/
│   ├── bit_identity_check.json
│   ├── numerical_health_metrics.json
│   └── crossing_comparison_table.csv
├── reports/
│   └── reproduction_verification_report.md
└── logs/
    ├── run_seedless_log.txt
    ├── comparison_log.txt
    └── verification_log.txt
```

## 合格基準

### 必須条件

1. ✓ 独立実行1・2が f, q でビット一致
2. ✓ 規格化誤差 ≤ 1.4×10⁻¹⁴
3. ✓ 零二乗閉鎖誤差 ≤ 1.3×10⁻¹³
4. ✓ 原本との crossing 一致

### 強い条件

5. ✓ 方向部分空間 $\mathcal{O}(P_{34}, P_{34}^{\text{orig}}) > 0.999$
6. ✓ 最大分裂フラクション一致（±1ULP以内）

## 参考

- **論文**: `../../../次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/論文8_二段階seed除去による準安定相の因果分離_完成論文_v1.md`（§2.1 実験設計、§10 再現性）
- **第7論文原本**: `../../../次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/code/` （SHA-256 固定利用）
- **original_environment**: 第8論文 paper8_stage_A2a_seedless_N5 の完全複製

---

**開始日**: 2026-08-25  
**再現対象 DOI**: https://doi.org/10.5281/zenodo.21614403  
**Concept DOI**: https://doi.org/10.5281/zenodo.21614402
