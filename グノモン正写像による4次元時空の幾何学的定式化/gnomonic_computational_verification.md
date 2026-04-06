# グノモン正写像による誘導計量の計算機検証

**著者**: 木原 範昭（Noriaki Kihara）  
**所属**: WF System Co., Ltd. / 大阪大学基礎工学部（卒業）  
**ORCID**: 0009-0004-6753-4020  
**状態**: テクニカルノート（計算検証報告）  
**日付**: 2026-04-06  
**関連論文**:  
- 論文1: Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection (DOI: 10.5281/zenodo.19427780)  
- 論文2: Geometric Origin of Lorentzian Spacetime via Gnomonic Projection (DOI: 10.5281/zenodo.19434932)

---

## 概要

本稿は、前論文（論文1・2）で提示されたグノモン正写像（中心射影）による誘導計量が、Einstein方程式を幾何学的恒等式として満たすという主張を、計算機代数システム **SymPy** により独立に検証した結果を報告する。

検証は以下の項目を含む：

1. 誘導計量 $g_{\mu\nu}$ と逆計量 $g^{\mu\nu}$ の整合性
2. Christoffel 記号の記号計算
3. Ricci テンソル・スカラーの検証
4. Einstein 方程式 $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ の恒等式としての成立
5. 一般 $n$ 次元（$n = 2, 3, 4, 5$）での検証
6. Lorentzian 署名（$+,+,+,-$）での検証

すべての項目において、理論的予測と計算結果が完全に一致することを確認した。

---

## 1. グノモン正写像の定義

### 1.1 写像

$n$ 次元ユークリッド空間 $\mathbb{R}^n$ から半径 $R$ の $n$ 次元超球面 $S^n(R)$ への中心射影（グノモン正写像）を次のように定義する：

$$\Phi: \mathbb{R}^n \to S^n(R), \qquad Y^\mu = \frac{R x^\mu}{l}, \quad Y^{n+1} = \frac{R^2}{l}$$

ここで

$$l = \sqrt{R^2 + \sum_{\mu=1}^{n} (x^\mu)^2}$$

### 1.2 誘導計量

この写像による引き戻し（pullback）計量は：

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\delta_{\mu\nu} - \frac{x^\mu x^\nu}{l^2}\right)$$

逆計量は：

$$g^{\mu\nu} = \frac{l^2}{R^2}\left(\delta^{\mu\nu} + \frac{x^\mu x^\nu}{R^2}\right)$$

### 1.3 Lorentzian 版

4次元の場合に Wick 回転（$\epsilon_4: +1 \to -1$）を適用すると：

$$l^2 = R^2 + x^2 + y^2 + z^2 - t^2$$

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\eta_{\mu\nu} - \frac{\tilde{x}_\mu \tilde{x}_\nu}{l^2}\right)$$

ここで $\eta_{\mu\nu} = \mathrm{diag}(+1, +1, +1, -1)$、$\tilde{x}_\mu = \eta_{\mu\nu}x^\nu$ である。

---

## 2. 理論的予測

$S^n(R)$ は断面曲率 $K = 1/R^2$ の最大対称空間であるから、微分幾何学の一般論より：

- **Ricci テンソル**: $R_{\mu\nu} = \frac{n-1}{R^2}\, g_{\mu\nu}$
- **Ricci スカラー**: $\mathcal{R} = \frac{n(n-1)}{R^2}$
- **Einstein テンソル**: $G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}\mathcal{R}\, g_{\mu\nu} = -\frac{(n-1)(n-2)}{2R^2}\, g_{\mu\nu}$
- **Einstein 方程式**: $G_{\mu\nu} + \Lambda\, g_{\mu\nu} = 0, \qquad \Lambda = \frac{(n-1)(n-2)}{2R^2}$

| 次元 $n$ | $\Lambda$ | Ricci スカラー $\mathcal{R}$ |
|:--------:|:---------:|:---:|
| 2 | 0 | $2/R^2$ |
| 3 | $1/R^2$ | $6/R^2$ |
| **4** | **$3/R^2$** | $12/R^2$ |
| 5 | $6/R^2$ | $20/R^2$ |

**本検証の目的は、グノモン写像の明示的な座標表示から出発し、上記の結果が再現されることを記号計算で確認することである。**

---

## 3. 検証方法

### 3.1 計算機代数システム

Python の記号計算ライブラリ **SymPy** (version 1.14.0) を使用した。すべての計算は記号的（symbolic）に実行され、数値近似は一切含まれない。

### 3.2 検証手順

各次元 $n$ に対し、以下の手順で検証を行った：

**P-1. 計量テンソルの検証**

計量 $g_{\mu\nu}$ と逆計量 $g^{\mu\nu}$ を明示的に構成し、$g_{\mu\alpha} g^{\alpha\nu} = \delta^\nu_\mu$ が成立することを確認。

**P-2. Christoffel 記号の計算**

$$\Gamma^\alpha_{\mu\nu} = \frac{1}{2} g^{\alpha\sigma}\left(\frac{\partial g_{\sigma\mu}}{\partial x^\nu} + \frac{\partial g_{\sigma\nu}}{\partial x^\mu} - \frac{\partial g_{\mu\nu}}{\partial x^\sigma}\right)$$

を記号的に計算。対称性 $\Gamma^\alpha_{\mu\nu} = \Gamma^\alpha_{\nu\mu}$ を利用して計算量を削減。

**P-3. Ricci テンソル・スカラーの計算と検証**

$$R_{\mu\nu} = \partial_\alpha \Gamma^\alpha_{\mu\nu} - \partial_\nu \Gamma^\alpha_{\mu\alpha} + \Gamma^\alpha_{\alpha\beta}\Gamma^\beta_{\mu\nu} - \Gamma^\alpha_{\nu\beta}\Gamma^\beta_{\mu\alpha}$$

を記号的に計算し、$R_{\mu\nu} - \frac{n-1}{R^2}g_{\mu\nu} = 0$ が全成分で成立するか検証。
Ricci スカラー $\mathcal{R} = g^{\mu\nu}R_{\mu\nu}$ を計算し、$n(n-1)/R^2$ と一致するか検証。

**P-4. Einstein 方程式の検証**

Einstein テンソル $G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}\mathcal{R}\, g_{\mu\nu}$ を構成し、

$$G_{\mu\nu} + \frac{(n-1)(n-2)}{2R^2}\, g_{\mu\nu} = 0$$

が全成分で成立するか検証。

**P-5. 一般次元の検証**

上記 P-1〜P-4 を $n = 2, 3, 4, 5$ で実行。

**P-8. Lorentzian 版の検証**

$n = 4$ で署名を $\eta = \mathrm{diag}(+1,+1,+1,-1)$ に変更し、$l^2 = R^2 + x^2 + y^2 + z^2 - t^2$ とした上で、同一の検証を実行。

### 3.3 検証プログラム

検証プログラムは `verification/` フォルダに格納されている：

| ファイル | 検証対象 |
|---------|---------|
| `verification/verify_euclidean.py` | P-1〜P-5（ユークリッド版、$n=2,3,4,5$） |
| `verification/verify_lorentzian.py` | P-8（Lorentzian版、$n=4$） |

実行方法：
```
pip install sympy
python3 verification/verify_euclidean.py
python3 verification/verify_lorentzian.py
```

---

## 4. 検証結果

### 4.1 ユークリッド版（P-1〜P-5）

| 次元 $n$ | P-1 $g \cdot g^{-1} = I$ | P-3 $R_{\mu\nu}$ | P-3 $\mathcal{R}$ | P-4 $G + \Lambda g = 0$ | 計算時間 |
|:--------:|:---:|:---:|:---:|:---:|:---:|
| 2 | ✅ | ✅ $R_{\mu\nu} = \frac{1}{R^2}g_{\mu\nu}$ | ✅ $\frac{2}{R^2}$ | ✅ $\Lambda = 0$ | 0.3秒 |
| 3 | ✅ | ✅ $R_{\mu\nu} = \frac{2}{R^2}g_{\mu\nu}$ | ✅ $\frac{6}{R^2}$ | ✅ $\Lambda = \frac{1}{R^2}$ | 0.8秒 |
| **4** | ✅ | ✅ $R_{\mu\nu} = \frac{3}{R^2}g_{\mu\nu}$ | ✅ $\frac{12}{R^2}$ | ✅ $\Lambda = \frac{3}{R^2}$ | 2.0秒 |
| 5 | ✅ | ✅ $R_{\mu\nu} = \frac{4}{R^2}g_{\mu\nu}$ | ✅ $\frac{20}{R^2}$ | ✅ $\Lambda = \frac{6}{R^2}$ | 4.1秒 |

### 4.2 Lorentzian版（P-8）

4次元 Lorentzian 署名 $(+,+,+,-)$、$l^2 = R^2 + x^2 + y^2 + z^2 - t^2$ での検証結果：

| 項目 | 結果 |
|-----|:----:|
| $g_{\mu\nu} \cdot g^{\mu\nu} = I_4$ | ✅ |
| $R_{\mu\nu} = \frac{3}{R^2}\, g_{\mu\nu}$ | ✅ |
| $\mathcal{R} = \frac{12}{R^2}$ | ✅ |
| $G_{\mu\nu} + \frac{3}{R^2}\, g_{\mu\nu} = 0$ | ✅ |

計算時間：2.2秒

---

## 5. 考察

### 5.1 検証の意義

本検証により、以下の事実が計算機代数により独立に確認された：

1. **グノモン正写像の誘導計量は、$S^n(R)$ の最大対称計量と一致する。** 明示的な座標表示から出発し、Christoffel 記号 → Riemann テンソル → Ricci テンソル → Einstein テンソルの全計算チェーンを通じて、$G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$（$\Lambda = (n-1)(n-2)/(2R^2)$）が任意の座標点で成立する恒等式であることを確認した。

2. **この結果は次元に依存しない。** $n = 2, 3, 4, 5$ のいずれでも、同一の公式 $\Lambda = (n-1)(n-2)/(2R^2)$ が成立する。特に $n = 4$ で $\Lambda = 3/R^2$ となる。

3. **Lorentzian 署名でも成立する。** Wick 回転により署名を $(+,+,+,+)$ から $(+,+,+,-)$ に変更しても、Einstein 方程式 $G_{\mu\nu} + (3/R^2)g_{\mu\nu} = 0$ は恒等式として成立する。

### 5.2 この結果が意味すること

$G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ は、通常の一般相対性理論では**宇宙定数 $\Lambda$ を持つ真空 Einstein 方程式の解**として扱われる。本結果は、同一の式がグノモン正写像の**純粋な幾何学的恒等式**として得られることを確認したものである。

### 5.3 本検証の限界

本検証で確認されたのは、誘導計量が最大対称空間の計量であるという**数学的事実**のみである。これを物理（重力理論）としてどう解釈するかは、本検証の範囲外であり、別途の検討を要する。

特に以下の点は本検証では扱っていない：

- 物質項 $T_{\mu\nu}$ の導入と重力定数 $G$ の関係
- 曲率半径 $R$ と観測量（Hubble定数、宇宙定数）の同定
- 第5次元（電荷次元）の導入とその妥当性

---

## 6. 結論

グノモン正写像の誘導計量に対する Einstein テンソルの計算を、計算機代数システム SymPy による記号計算で独立検証した。4次元を含む $n = 2, 3, 4, 5$ のすべての次元で、また Lorentzian 署名においても、$G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$（$\Lambda = (n-1)(n-2)/(2R^2)$）が恒等式として成立することを確認した。

検証プログラムは公開されており（`verification/` フォルダ）、誰でも再現可能である。

---

## 参考文献

1. Kihara, N. (2026). *Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection*. Zenodo. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)
2. Kihara, N. (2026). *Geometric Origin of Lorentzian Spacetime via Gnomonic Projection*. Zenodo. DOI: [10.5281/zenodo.19434932](https://doi.org/10.5281/zenodo.19434932)
3. do Carmo, M. P. (1992). *Riemannian Geometry*. Birkhäuser.
4. Wald, R. M. (1984). *General Relativity*. University of Chicago Press.

---

*本稿の検証はすべて記号的（exact）であり、数値近似を含まない。*  
*検証環境：Python 3.9, SymPy 1.14.0, macOS*  
*最終更新：2026-04-06*
