# 論文3：離散構造 — 4次元球の単位立方体充填と $\Delta(R)$ の漸近解析 構成草案

**仮タイトル**: Discrete Structure of a Four-Dimensional Ball: Unit-Cube Packing and the Asymptotic Volume Deficit $\Delta(R)$

**分量目安**: 15–25 ページ

---

## 目次

| § | タイトル | 分量 | 状態 |
|---|---|---|---|
| §1 | Introduction | 2 ページ | 未着手 |
| §2 | Problem Setting | 2 ページ | 試案完成（[draft_sections_2_4.md](draft_sections_2_4.md)）|
| §3 | Formulation of the Packing Condition | 2 ページ | 試案完成 |
| §4 | The Sequence $N(k)$: Combinatorial Properties | 3 ページ | 試案完成 |
| §5 | Asymptotic Analysis: $\Delta(R) \sim cR^3$ | 4 ページ | 草稿（math_report.md）|
| §6 | Coefficient $c$ and the Lagrange–Jacobi Connection | 3 ページ | 草稿（math_report.md）|
| §7 | Numerical Verification | 2 ページ | 草稿（math_report.md）|
| §8 | Conclusion | 1 ページ | 未着手 |

---

## §5 漸近解析の概要

包除原理による境界補正：
$$\mathrm{Vol}(\Omega_R) = V_4(R) - \frac{16\pi R^3}{3} + 6\pi R^2 + O(R)$$

格子点誤差 $E(R) = O(R^3)$ で振動成分平均ゼロを仮定（数値で支持）すると：
$$\Delta(R) = V_4(R) - N(k) \sim \frac{16\pi}{3} R^3$$

正規化定数：$c = \frac{16\pi/3}{2\pi^2} = \frac{8}{3\pi}$。

詳細は `computations/results/math_report.md` §3 を参照。

---

## §6 Jacobi 接続の概要

変数変換 $p_i = 2|c_i| + 1$ により：
$$\sum p_i^2 \le (4k+2)^2, \quad p_i \in \{1, 3, 5, \ldots\}$$

Jacobi の四平方定理：
$$r_4(N) = \begin{cases} 8\sigma(N) & 4 \nmid N \\ 24\sigma(N_{\mathrm{odd}}) & 4 \mid N \end{cases}$$

全 $p_i$ が奇数のとき $\sum p_i^2 \equiv 4 \pmod 8$。よって $4 \mid N$ かつ $8 \nmid N$ の表現が関係する。

正確な係数（**後で検討**）は数値検証と整合させる必要がある。基本構造：
$$N(k) \sim (\text{const}) \sum_{N \le (4k+2)^2, N \equiv 4 \pmod 8} \sigma(N_{\mathrm{odd}})$$

---

## §7 数値検証の概要

$k_{\min} = 30$ の3次多項式フィット：
- $c_3$ = 16.7511 vs 理論 $16\pi/3 = 16.7552$（誤差 0.024%）
- $c_2$ = $-22.22$ vs 理論 $-6\pi = -18.85$（オーダー一致、定数差は格子誤差）

詳細データは `computations/results/N_table_k0_60.tsv` および `math_report.md` §4。

---

## §1 序論案

中心的な問い：「半径 $R$ の4次元球に単位立方体は何個詰め込めるか？」

これは以下の点で興味深い：

1. **古典数論との接続**：Lagrange の四平方定理（任意の非負整数は4平方の和）は4次元の特殊性の数論的表現。Jacobi の四平方表現数 $r_4(N)$ は単純な閉形式を持ち、これは4次元固有の現象。

2. **Gauss 円問題の高次元版**：Gauss 円問題は2次元での格子点計数問題。本論文は4次元版で、特に「内側に完全包含される単位立方体」という制約付き問題を扱う。

3. **物理学への動機（背景化）**：4次元球の体積と充填可能な単位立方体数の差 $\Delta(R)$ は、ある種の「離散性に由来する欠損」として解釈可能。本論文ではこの解釈は背景にとどめ、数学的内容に集中する。

中心結果：$\Delta(R) \sim (16\pi/3) R^3$、係数 $c = 8/(3\pi)$ は包除原理から導出され、数値で確認される。

---

## §8 結論案

- 充填条件と $N(k)$ の正確な定義
- 漸近係数 $c = 8/(3\pi) \approx 0.84883$ の理論的決定と数値検証
- Jacobi 接続による $N(k)$ の数論的構造
- 残課題：Jacobi 接続の正確な係数決定、副主要項の数値精度向上

---

**自己引用なし** / **物理的解釈最小限** / **純粋数学の論文として完結**
