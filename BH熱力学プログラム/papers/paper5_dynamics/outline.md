# 論文5：動的接続 — $R$ の自己整合的決定と Hawking 蒸発 構成草案

**仮タイトル**: Self-consistent Determination of the Curvature Radius $R$ and Recovery of the Tangherlini Evaporation Scaling

**分量目安**: 10–15 ページ

---

## 目次

| § | タイトル | 分量 | 状態 |
|---|---|---|---|
| §1 | Introduction | 1.5 ページ | 未着手 |
| §2 | Total Energy and Curvature Radius | 2.5 ページ | 草稿要点 |
| §3 | Mass–Radius Relation | 2 ページ | 草稿要点 |
| §4 | Hawking Temperature and Lifetime | 2.5 ページ | 草稿要点 |
| §5 | The Evaporation Rate Discrepancy | 2 ページ | 草稿要点 |
| §6 | Discussion | 1 ページ | 未着手 |

---

## §2 全エネルギーと曲率半径

### §2.1 主観空間内の全エネルギー

論文2 の中心投影由来計量 $g_{\mu\nu}$ における全エネルギーを、Friedmann 類似の構造から導出：

$$E_{\mathrm{total}} = \int \rho \sqrt{-g} \, d^3x \propto R^2$$

詳細：宇宙項 $\Lambda = 3/R^2$ と Friedmann 方程式の構造から、4次元主観空間内の全エネルギーが曲率半径 $R$ の2乗に比例する関係が導かれる。

### §2.2 自己整合方程式

質量 $M$ を持つブラックホールに対して：
$$E_{\mathrm{total}} = M c^2 \quad \Leftrightarrow \quad M \propto R^2$$

逆に解いて：
$$\boxed{\;R \propto M^{1/2}\;}$$

これは 4+1 次元 Tangherlini 解の地平面半径 $r_h \propto M^{1/2}$ と完全一致。

---

## §3 質量–半径関係

### §3.1 Tangherlini 解との一致

4+1 次元 Schwarzschild–Tangherlini：$r_h^2 = 8 G M / (3\pi)$、すなわち $r_h \propto M^{1/2}$。

本論文の主観空間導出：$R \propto M^{1/2}$。

両者が同一の質量–半径スケーリングを与えることは、論文4 で仮定した $R = r_h$ 同定が動的に正当化されることを意味する。

### §3.2 比例定数

比例定数の決定は **後で検討**：5次元 Newton 定数 $G_5$ と Planck 長さの関係を含む正規化に依存。

---

## §4 Hawking 温度と寿命

### §4.1 温度

地平面表面重力 $\kappa = 1/R$（4+1 次元 Tangherlini で $\kappa = 1/r_h$）：
$$T_H = \frac{\kappa}{2\pi} = \frac{1}{2\pi R} \propto R^{-1} \propto M^{-1/2}.$$

これは 4+1 次元 Tangherlini と完全一致。

### §4.2 寿命

Stefan–Boltzmann 様の蒸発レートを仮定：
$$\frac{dM}{dt} = -\sigma T^4 A_3$$

$T \propto M^{-1/2}$, $A_3 \propto R^3 \propto M^{3/2}$ より：
$$\frac{dM}{dt} \propto -M^{-2} \cdot M^{3/2} = -M^{-1/2}$$

→ $\int M^{1/2} dM = -t/\tau$ → $M^{3/2} \propto -t$ → 寿命 $\propto M^2$.

これは 4+1 次元 Tangherlini と一貫。

---

## §5 蒸発レート指数の不一致

### §5.1 単純仮定での結果

§4 の Stefan–Boltzmann 仮定では：
$$\frac{dM}{dt} \propto -M^{-1/2}, \qquad n = -1/2$$
（または別の係数選択で $n = 3/2$, **後で検討**：正確な指数の特定）

### §5.2 標準 Hawking 計算

標準的な半古典的計算：
$$\frac{dM}{dt} \propto -M^{-1}, \qquad n = -1$$

### §5.3 不一致の明示

$n = -1/2$ vs $n = -1$ の差。本論文では完全一致を主張せず、不一致点として明示する。

可能な解決方向（**後で検討**）：
- 離散ドリフト機構の動的な記述
- 反作用効果の取り込み
- 4+1 次元から有効 3+1 次元への射影効果

これらは論文6 または将来論文での課題として位置づける。

---

## §6 議論

- 主要スケーリング（$R \propto M^{1/2}$, $T_H \propto M^{-1/2}$, 寿命 $\propto M^2$）は完全一致
- 蒸発レート指数のみ不一致、これは離散ドリフト機構の精密化が必要
- 部分的成功と残課題の明示

---

**自己引用なし** / **完全一致を主張しない誠実なアプローチ**

**未確定事項（後で検討）**:
- §2 の Friedmann 類似導出の正確な書き方
- §4.1 の Hawking 温度の表面重力からの導出
- §5 の指数 $n$ の正確な値（$n = -1/2$ or $n = 3/2$ など、定義依存）
