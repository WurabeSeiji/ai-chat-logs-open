---
title: "球面投影 (radial projection) の定義と中心投影との関係 ── テクニカルノートの公開"
emoji: "🎯"
type: "idea"
topics: ["数学", "幾何学", "球面投影", "中心投影", "位相幾何学"]
published: true
---

## はじめに

著者の中心投影シリーズ（[基礎論文](https://doi.org/10.5281/zenodo.19427780)、[合成演算論文](https://doi.org/10.5281/zenodo.20060728)）の **基礎写像** を明示的に定義し整理するテクニカルノートを Zenodo に公開しました。

本ノートが扱う写像は、極めて初等的なものです：

$$
\sigma_R: \mathbb{R}^{n+1} \setminus \{0\} \to S^n(R), \quad \sigma_R(x) = \frac{R}{\|x\|} \cdot x
$$

つまり「ベクトルを正規化して半径 $R$ の球面に乗せる」だけの写像です。位相幾何学では **radial projection（動径射影）**、線形代数では **正規化写像 (normalization map)** として古典的に既知の操作で、Hatcher『Algebraic Topology』Chapter 0 で deformation retract の典型例として登場します。

- **Concept DOI**: https://doi.org/10.5281/zenodo.20462569
- **v3.1 DOI**: https://doi.org/10.5281/zenodo.20462570
- **GitHub**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/グノモン正写像による4次元時空の幾何学的定式化

なぜこんなに初等的な写像を改めて定義するのか？　それは、著者の中心投影シリーズで暗黙に使われてきた基礎操作を、**新規定理の主張なしに**、シリーズ内の参照基盤として整理するためです。

## 球面投影 $\sigma_R$ とは

### 定義

$$
\sigma_R: \mathbb{R}^{n+1} \setminus \{0\} \to S^n(R), \quad \sigma_R(x) = \frac{R}{\|x\|} \cdot x
$$

原点以外のあらゆる点を、半径 $R$ の球面に「放射状に」投影します。

![Fig. 1 Radial Projection](https://github.com/WurabeSeiji/ai-chat-logs-open/raw/main/グノモン正写像による4次元時空の幾何学的定式化/fig1_radial_projection.png)

円の内側にある点（A）も、外側にある点（B, C）も、原点を通る同じ半直線上の球面交点に投影されます。

### 基本的な性質

球面投影は次の性質を持ちます：

1. **滑らか**: $C^\infty$ 級
2. **全射**: 像は球面 $S^n(R)$ 全体
3. **冪等**: $\sigma_R \circ \sigma_R = \sigma_R$（だから projection と呼ぶ）
4. **正のスカラーで不変**: $\sigma_R(tx) = \sigma_R(x)$（$t > 0$）
5. **強変形レトラクト**: $\mathbb{R}^{n+1} \setminus \{0\}$ が球面 $S^n(R)$ にホモトピー同値
6. **商空間表示**: $(\mathbb{R}^{n+1} \setminus \{0\}) / \mathbb{R}_{>0} \cong S^n(R)$

つまり、**動径方向の情報を完全に捨てて、角度方向だけを取り出す**写像です。

### 微分構造

各点 $x$ における微分 $D\sigma_R|_x$ は階数 $n$ の沈め込みで、

$$
\ker D\sigma_R|_x = \mathrm{span}\{x\}, \quad \mathrm{Im}(D\sigma_R|_x) = x^\perp = T_{\sigma_R(x)} S^n(R)
$$

動径方向は核に入り、接方向だけが球面に写ります。これが「動径成分を捨てる」ことの微分幾何学的な定式化です。

### 角度成分の保存

球面投影は、二つのベクトル間の角度を保存します：

$$
\angle(\sigma_R(x), \sigma_R(y)) = \angle(x, y)
$$

**注意**: これは **共形性**（曲線間の角度を保存する性質）とは異なる概念です。$\sigma_R$ は次元を落とすので、曲線間の角度の保存は意味を成しません。あくまで「2 つのベクトル間の角度」の話です。

### スケール不変性

異なる半径 $R_1, R_2$ に対しても、原点から $x$ への半直線上に像が乗ります：

![Fig. 3 Radial Projection (double R)](https://github.com/WurabeSeiji/ai-chat-logs-open/raw/main/グノモン正写像による4次元時空の幾何学的定式化/fig3_radial_projection_double_R.png)

半径 $R=3$ の円（黒）と $R=5$ の円（紫）の双方で、各点 $A, B, C$ の像は同一の放射状半直線上に並びます。角度 $\theta_A, \theta_B, \theta_C$ は $R$ に依らず保存されます。

## 中心投影との関係

### 中心投影 $\Phi_R$ の定義

球面 $S^n(R)$ の北極 $N = (0, \ldots, 0, R)$ における接超平面

$$
\Pi_R = \{x \in \mathbb{R}^{n+1} \mid x_{n+1} = R\}
$$

から球面への中心投影 $\Phi_R$ を以下で定めます（[基礎論文](https://doi.org/10.5281/zenodo.19427780)の定義）：

$$
\Phi_R: \Pi_R \to S^n_+(R), \quad \Phi_R(x) = \frac{R}{\|x\|} \cdot x
$$

ここで $S^n_+(R) = \{y \in S^n(R) \mid y_{n+1} > 0\}$ は **開上半球面** です。

### 中心投影は球面投影の制限

これら 2 つの写像、$\sigma_R$ と $\Phi_R$ は、**写像の式が完全に同じ**です：

$$
\sigma_R(x) = \Phi_R(x) = \frac{R}{\|x\|} \cdot x
$$

違いは **定義域** だけです：

$$
\Phi_R = \sigma_R \big|_{\Pi_R}
$$

つまり、**中心投影は球面投影を接超平面 $\Pi_R$ に制限したもの**です。

![Fig. 2 Central Projection](https://github.com/WurabeSeiji/ai-chat-logs-open/raw/main/グノモン正写像による4次元時空の幾何学的定式化/fig2_central_projection.png)

紫線が接超平面 $y = R$ で、その上の点 $A'', B''$ を放射状に球面上の $A', B'$ に投影します。

### 一般化と単射性のトレードオフ

両者を対比すると：

| 性質 | $\sigma_R$ | $\Phi_R$ |
|---|---|---|
| 定義域 | $\mathbb{R}^{n+1} \setminus \{0\}$ | $\Pi_R$ |
| 像 | $S^n(R)$（全球面） | $S^n_+(R)$（**開上半球面のみ**） |
| 単射性 | **非単射**（半直線を 1 点に潰す） | **単射**（実は微分同相） |
| 幾何構造 | 動径方向が微分の核に入る | 豊かな誘導計量を持つ |

球面投影は定義域を全空間に拡張する代償として **単射性を失います**。一方、中心投影は単射性と豊かな幾何構造（誘導計量、曲率テンソル等）を持ちますが、定義域は接超平面、像は開上半球面に限定されます。

これが本ノートの **核心的な対比** です。

## なぜこのテクニカルノートを書いたのか

本ノートは **新規の幾何学的定理を主張しません**。

ベクトルを正規化して球面に乗せる写像は、位相幾何学・微分幾何学では 100 年来の常識です。本ノートが行うのは：

1. この既知の操作に、著者の中心投影フレームワーク内で用いる **記号 $\sigma_R$ を固定する**
2. 中心投影 $\Phi_R$ がこの上位概念の特殊化であることを **明示的に整理する**
3. 既存の中心投影シリーズの **参照基盤として一箇所にまとめる**

これだけです。Zenodo の「プレプリント」よりも **テクニカルノート** という性格付けに近い文書です。

## 4 AI 査読のプロセス

このノートは、初稿（v1）から最終版（v3.1）まで、**Claude.ai / ChatGPT / Gemini / Grok の 4 AI を 2 ラウンド経て精査**しました。

主な修正：

- **v1 → v2**: 「上位概念」という主張過剰の弱化、radial projection との同一性明示、像 $S^n_+(R)$ の厳密化、命題追加（冪等性・商空間・微分の核）、Hatcher 引用追加
- **v2 → v3**: ChatGPT が指摘した命題 2.4 の核証明の論理エラー修正、§3.5 の「well-defined でない」という数学的誤りを「well-defined だが退化する」に修正
- **v3 → v3.1**: Claude.ai 指摘の微分の像 $\mathrm{Im}(D\sigma_R|_x) = x^\perp$ 追記、Gemini 指摘の Fig.2 キャプション座標表記補足

特に v3 で ChatGPT が検出した 2 つの数学的誤りは、Grok と Gemini は見逃していました。**複数 AI 査読の併用が品質確保に有効** という実例にもなっています。

## まとめ

球面投影 $\sigma_R(x) = (R/\|x\|) x$ は、位相幾何学・微分幾何学では古典的に既知の写像です。本ノートはこの写像に著者の中心投影フレームワーク内で用いる記号を固定し、中心投影 $\Phi_R$ が $\sigma_R$ を接超平面 $\Pi_R$ に制限したものに一致することを明示しました。

新規性は記号と位置付けの明確化に限られますが、既存の中心投影シリーズの基礎を固める参照基盤として、シリーズ全体の再現性と引用可能性が向上しました。

論文本体（md / tex / pdf × 日英）+ 図 3 点 = 計 9 ファイルを Zenodo で公開しています：

- **Concept DOI**: https://doi.org/10.5281/zenodo.20462569
- **GitHub**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/グノモン正写像による4次元時空の幾何学的定式化

ライセンスは CC BY 4.0 です。

---

## 関連記事

- [中心投影による次元削減を正しく代数化する ── 1 回の中心投影と球面上の可換切断](https://zenn.dev/noriaki_kihara/articles/central-projection-composition) — 本テクニカルノートが基礎を固める合成演算論文の解説
- [中心投影による4次元空間の幾何学的定式化](https://zenn.dev/noriaki_kihara/articles/gnomonic-projection-spacetime-geometry) — 中心投影シリーズの基礎論文
