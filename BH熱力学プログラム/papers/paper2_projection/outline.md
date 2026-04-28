# 論文2：中心投影と4次元主観計量 構成草案

**仮タイトル**: Gnomonic Projection from $S^4$ and the Schwarzschild-like Form of the Induced Metric on a Four-Dimensional Subjective Chart

**分量目安**: 12–18 ページ

---

## 目次

| § | タイトル | 分量 | 状態 |
|---|---|---|---|
| §1 | Introduction | 1.5 ページ | 未着手 |
| §2 | The Construction: Gnomonic Projection | 2.5 ページ | 草稿要点 |
| §3 | Induced Metric and its Properties | 3 ページ | 草稿要点 |
| §4 | Behaviour at $R$-fixed and the Schwarzschild-like Form | 3 ページ | 草稿要点 |
| §5 | Geodesics and Horizon Structure | 3 ページ | 草稿要点 |
| §6 | Comparison with the Tangherlini Solution | 2 ページ | 草稿要点 |
| §7 | Discussion and Open Problems | 1 ページ | 未着手 |

---

## §2 グノモン正写像の構成

### §2.1 設定

5次元背景空間 $\mathbb{R}^5$ に4次元球面
$$S^4(R) = \{Y \in \mathbb{R}^5 : Y_0^2 + Y_1^2 + \cdots + Y_4^2 = R^2\}$$
を埋め込む。$Y_0$-軸を投影軸として選ぶ。

中心投影（グノモン正写像）$\Phi : \Pi_R \to S^4(R)$ は、4次元接超平面
$$\Pi_R = \{Y \in \mathbb{R}^5 : Y_0 = R\}$$
の点 $x = (x_1, x_2, x_3, x_4) \in \Pi_R$ を、$S^4(R)$ 上の点
$$\Phi(x) = \frac{R}{\ell} (x_1, x_2, x_3, x_4, R/\ell \cdot R) = \cdots$$
（厳密表式は **後で検討**）に写す。

要点：$\ell = \sqrt{R^2 + \|x\|^2}$ として、写像は球面の中心から接超平面への射影。

### §2.2 引き戻し計量

$\Phi$ により $S^4(R)$ 上の標準計量を $\Pi_R$ に引き戻すと
$$g_{\mu\nu}(x) = \frac{R^2}{\ell^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{\ell^2}\right).$$

これは4次元主観計量と呼ぶ。性質：
- $\det g$, 逆計量 $g^{\mu\nu}$, Christoffel 記号, Riemann テンソル, Einstein テンソルを計算可能
- $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ with $\Lambda = 3/R^2$（4次元での宇宙項付き真空 Einstein 方程式）

---

## §3 計量の性質

### §3.1 球対称形

中心からの動径距離 $r = \|x\|$ で書き直すと：
$$ds^2 = -\left(\frac{R^2}{\ell^2}\right) d\tau^2 + \frac{R^2}{\ell^2}\left(dr^2 + r^2 d\Omega_2^2\right) + (\text{cross terms})$$
（**後で検討**：時間化した形式の正確な符号と項の整理）

### §3.2 物理的解釈の最小化

物理的解釈は本論文では扱わず、純粋に幾何学的性質に集中。

---

## §4 Schwarzschild 類似形

$R$-固定の射影部分多様体上で、計量を Schwarzschild 形に整理：
$$ds^2 \to f(r) dt^2 - f(r)^{-1} dr^2 - r^2 d\Omega^2$$
の形を求める。$f(r) = 1 - r^2/R^2$ + 補正項（**後で検討**：補正項の正確な形）。

主要結果：$R$-依存の宇宙項を持つ de Sitter 類似計量が現れる。

---

## §5 測地線と地平面

$f(r) = 0$ となる動径 $r_h$ が地平面候補。Schwarzschild 類似計量の解析。

---

## §6 Tangherlini 解との比較

4+1 次元 Schwarzschild–Tangherlini 計量：
$$ds^2 = -f_T(r) dt^2 + f_T(r)^{-1} dr^2 + r^2 d\Omega_3^2$$
where $f_T(r) = 1 - \mu/r^2$, $\mu \propto M$.

本論文の中心投影由来計量との比較：両者の $r$ 依存性、地平面構造、漸近的挙動。

---

## §7 議論と残課題

- **後で検討**：物理的解釈をどこまで本論文で扱うか
- 4+1 次元 Tangherlini 解への正確な還元手続き
- 質量パラメータ $M$ の中心投影由来の解釈

---

**未確定事項（後で検討）**:
- §2 の写像の正確な表式
- §3, §4 の計量の項の符号と表記の選択
- §5 の地平面構造の正確な記述
- §6 の Tangherlini 解との対応の詳細

執筆方針：論文2 は数学的内容に集中し、物理的解釈は §6 で控えめに触れる程度にとどめる。
