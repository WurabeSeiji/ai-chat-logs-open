# Phase 4-6 実験計画書

## 目標

**外部仮定なしに、R_137 と R_128 に対応する 137 と 128.946 を、
閉鎖交換系の第一原理から導出する。**

---

## Phase 1-3 からの知見

### 確認された事実

1. **対称な双対アトラクター構造**
   - R_128 と R_137 での内部状態が A↔B 交換で同一
   - 同じ力学系が二つの等価な吸引状態を持つ

2. **正規化条件による制約**
   - |A|² + |B|² = 1 が固定
   - E_close が常に ≈1.0（変動不可）

3. **初期からの構造**
   - 反復回数の少ない段階で既に本質的な分化が見える
   - 反復は「形成」ではなく「増幅」

### 新しい検証が必要な項目

- 流域図（どの初期条件がどのアトラクターへ行くか）
- リアプノフ指数（吸引力の定量化）
- 固有値条件から R を逆算
- 137 の数学的起源

---

## Phase 4: 流域図（Basin of Attraction）

### 目的

初期条件 (φ, s₀) 平面での流域分離を描き、
R_128 と R_137 での A↔B 対称性を可視化。

### 実装戦略

```python
for phi in linspace(0, 2π, 200):
    for s0 in linspace(0.001, 0.05, 100):
        # R_137 で反復 512 ステップ
        final_state_137 = evolve(R_137, phi, s0, 512)
        
        # R_128 で反復 512 ステップ
        final_state_128 = evolve(R_128, phi, s0, 512)
        
        # どちらへ引き寄せられたか判定
        if converges_to_137(final_state_137):
            color = "red"
        elif converges_to_128(final_state_128):
            color = "blue"
        else:
            color = "gray"
        
        plot(phi, s0, color)
```

### 検証項目

1. **流域境界が滑らかか（A↔B 対称か）**
   - φ=0 と φ=π での対称性
   - s₀=0.01 での鏡面対称性

2. **流域面積の比較**
   - R_137 の流域サイズ
   - R_128 の流域サイズ
   - 二つが対称的に同じサイズか

3. **カオス的散乱の証拠**
   - 流域が複雑に入り込んでいるか（フラクタル境界）
   - または滑らかで単純か

### 出力

- `phase4_basin_137_128_heatmap.png`：初期条件平面での流域図
- `phase4_basin_statistics.json`：流域サイズ、対称性指標

---

## Phase 5: 固有値・固定点の第一原理解析

### 目的

散乱行列の固有値条件から、R_137 と R_128 が安定点である理由を
方程式の形で導出する。

### 仮説

反復散乱写像 $T_R(X_n) = S_R M_R X_n$ において、
安定な周期軌道または固定点が存在する条件は：

$$\text{固有値条件} = 0$$

例えば、1-周期固定点なら：
$$\lambda(R) = e^{i\theta(R)}$$

で、$|\lambda(R)| \lesssim 1$ が安定条件。

### 実装戦略

```python
# 1. 伝達行列の固有値を（複数R値で）精密計算
for R in linspace(0.686, 0.700, 1000):
    M = compute_transfer_matrix(R, num_steps=256)
    eigenvalues = eig(M)
    
    # 固有値の絶対値
    mags = abs(eigenvalues)
    
    # 位相
    phases = angle(eigenvalues)
    
    # 固有値の組み合わせ条件
    # 例：λ₁ + λ₂ = ?
    trace = sum(eigenvalues)
    
    # 離散固有値の周期条件
    # λⁿ = 1 となる n を探す
    for n in range(1, 256):
        if abs(eigenvalues[0]**n - 1) < 1e-10:
            print(f"R={R:.10f}: λ^{n}=1 (resonance)")

# 2. 安定条件を定式化
# D(R) = det(M - λI) = 0 の根
# または安定性判別式：|Tr(M)| = 2

stable_R_values = []
for R in linspace(0.686, 0.700, 1000):
    discriminant = abs(trace(M)) - 2
    if abs(discriminant) < threshold:
        stable_R_values.append(R)
```

### 検証項目

1. **固有値の単位円上での配置**
   - R_137 と R_128 で固有値がどのように異なるか
   - 安定性判別式が最小化される R

2. **共鳴条件の発見**
   - λⁿ = 1 となる R と n の対応
   - 31, 124 などの周期との関連

3. **Trace と Determinant の極値**
   - dTr/dR = 0 となる R
   - 最大吸引性の点

### 出力

- `phase5_eigenvalue_spectrum.csv`：各R での固有値
- `phase5_stability_discriminant.csv`：安定性判別式
- `phase5_resonance_map.json`：λⁿ=1 となる (R, n) 対
- `phase5_fixed_point_analysis.md`：導出された安定条件式

---

## Phase 6: N(R) = 4π/(1-R)² の第一原理導出

### 目的

**137 という数値を、閉塞系の力学から自動的に導出する。**

### 仮説的フレーム

#### A. 位相空間の次元性

4次元状態空間 $(A_{\text{Re}}, A_{\text{Im}}, B_{\text{Re}}, B_{\text{Im}})$ で、
クロスセクション（ポアンカレ切断面）をとったとき、
許容状態数（または区別可能セル数）が R に依存する。

$$N(R) = \text{（有効次元数）} \times f(R)$$

#### B. 散乱強度との関係

交換係数 R が小さい（反射が強い）と、
システムが「狭い」領域に限定される。

$$\text{有効容量} \propto \frac{1}{(1-R)^2}$$

因子 4π は、単位球面の表面積。

### 実装戦略

```python
# 1. ポアンカレ断面での状態点の計数
for R in linspace(0.686, 0.700, 141):
    states_poincare = []
    
    for phi in linspace(0, 2π, 1000):
        for s0 in linspace(0.001, 0.05, 500):
            # 反復 256 ステップ
            X_n = evolve(R, phi, s0, 256)
            
            # ポアンカレ切断面での交差を記録
            # （例：Im(A) = 0 の平面）
            if crosses_poincare_section(X_n):
                states_poincare.append(X_n)
    
    # 状態の「広がり」を測定
    # - 占有領域の面積
    # - 区別可能状態の個数（density）
    occupancy_measure = len(states_poincare)
    
    # N(R) 候補を計算
    N_measured = occupancy_measure
    N_theory = 4 * pi / (1 - R)**2
    
    print(f"R={R:.10f}: N_measured={N_measured}, N_theory={N_theory:.1f}")
```

#### 代替法：情報理論的測度

```python
# 2. エントロピー・情報量での次元評価
for R in linspace(0.686, 0.700, 141):
    trajectory = evolve_full_trajectory(R, phi=0, s0=0.01, steps=512)
    
    # 位置空間でヒストグラム作成
    hist_A_real, _ = histogram(real_parts(trajectory['A']), bins=20)
    hist_A_imag, _ = histogram(imag_parts(trajectory['A']), bins=20)
    hist_B_real, _ = histogram(real_parts(trajectory['B']), bins=20)
    hist_B_imag, _ = histogram(imag_parts(trajectory['B']), bins=20)
    
    # Shannonエントロピー
    S_A = entropy(hist_A_real) + entropy(hist_A_imag)
    S_B = entropy(hist_B_real) + entropy(hist_B_imag)
    S_total = S_A + S_B
    
    # エントロピーと N(R) の関連
    N_from_entropy = exp(S_total / log(2))
```

#### 第3法：離散基礎からの計数

```python
# 3. 離散セル（整数格子）での直接計数
for R in linspace(0.686, 0.700, 141):
    # 正規化された状態 (A_Re, A_Im, B_Re, B_Im) に対し
    # 有効半径を定義
    effective_radius = compute_effective_radius(R)
    
    # 単位セル（幅 1）で完全包含される個数
    # Paper7 の「R=3 で 137 個」の類推
    cell_count = count_cells_in_ball(effective_radius, cell_width=1)
    
    # N(R) の計算
    N_computed = cell_count
    N_theory = 4 * pi / (1 - R)**2
    
    match = abs(N_computed - N_theory) / N_theory < 0.01
    print(f"R={R:.10f}: N_computed={N_computed}, match={match}")
```

### 検証項目

1. **N(R) ∝ 1/(1-R)² の確認**
   - スケーリング指数 α が 2 か
   - 係数 4π が出現するか

2. **R_137 での N 値**
   - N(0.697178) = ?
   - 137.036 への一致度

3. **R_128 での N 値**
   - N(0.688364) = ?
   - 129.394 への一致度

4. **中間値での N(R)**
   - 単調性の確認
   - 離散跳躍の有無

### 出力

- `phase6_poincare_section_density.csv`：ポアンカレ断面での状態密度
- `phase6_entropy_analysis.csv`：各R でのエントロピー
- `phase6_cell_counting_result.csv`：離散セル計数
- `phase6_N_R_derivation.md`：N(R) の第一原理導出
- `phase6_137_128_confirmation.json`：137.036 と 129.394 の確認

---

## 実験順序とマイルストーン

| Phase | 目標 | 期待される達成 | 依存関係 |
|---:|---|---|---|
| 4 | 流域図を描く | A↔B 対称性の可視化 | Phase 1-3 基盤 |
| 5 | 固有値解析 | 安定条件式の導出 | Phase 4 並行可能 |
| 6 | N(R) 導出 | 137 の第一原理導出 | Phase 5 の結果を使用 |

### 成功基準

✓ **Phase 4 成功**：流域が A↔B 対称で、同じサイズ

✓ **Phase 5 成功**：固有値条件から R_137, R_128 が導出される

✓ **Phase 6 成功**：
```
N(R_137) = 4π/(1-0.697178)² ≈ 137.036
N(R_128) = 4π/(1-0.688364)² ≈ 129.394
```
が計算から自動的に出現

---

## 全体的な意図

### 層別的な理解

```
外観レベル（Phase 1-3）:
  - 二つのピークが存在する
  - 対称な双対構造
  
力学レベル（Phase 4-5）:
  - 流域の配置
  - 固有値の条件
  
数学的起源（Phase 6）:
  - なぜ 137 なのか
  - 閉塞系から自動的に出現する理由
```

### 論文化への道筋

Phase 4-6 すべて成功すれば、論文タイトル候補：

> **「閉鎖二波交換系における微細構造定数の第一原理導出：
> 対称双対アトラクターと 4 次元完全包含セルの統一」**

---

作成日: 2026-07-17  
次実行予定: 直後
