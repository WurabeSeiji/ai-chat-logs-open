# 中心投影・球面投影の幾何 ── 多軸モデルと、その物理的舞台への一瞥（観察ノート）

**木原 範昭**（Noriaki Kihara）／ WF System Co., Ltd.／ 大阪大学 基礎工学部 卒
ORCID 0009-0004-6753-4020 ／ 連絡先 kihara.noriaki@gmail.com ／ CC BY 4.0 ／ 2026-06

> **観察ノート。** 新規の幾何学的定理は主張しない。古典的写像（radial projection）を整理し、その上に現れる構造を観察する。物理的導出は行わない。

## 1. 球面投影 $\sigma_R$（radial projection）

$$\sigma_R:\mathbb{R}^{n+1}\setminus\{0\}\to S^n(R),\qquad \sigma_R(x)=\frac{R}{\|x\|}\,x$$

- 位相幾何の **radial projection**（Hatcher *Algebraic Topology* Ch.0 の変形レトラクト）と同一。$C^\infty$ レトラクションで冪等 $\sigma_R\circ\sigma_R=\sigma_R$。
- 商空間 $(\mathbb{R}^{n+1}\setminus\{0\})/\mathbb{R}_{>0}\cong S^n(R)$。微分の核 $=$ 動径方向 $\mathrm{span}\{x\}$、像 $=x^{\perp}=T_{\sigma_R(x)}S^n(R)$（階数 $n$ の沈め込み）。角度保存、半径スケールに対する方向不変。

## 2. 中心投影 $\Phi_R$（gnomonic ＝ 球面投影の接平面への制限）

北極で接する接超平面 $\Pi_R=\{x_{n+1}=R\}\cong\mathbb{R}^n$ への制限：

$$\Phi_R=\sigma_R\big|_{\Pi_R}:\ \Pi_R\ \xrightarrow{\ \sim\ }\ S^n_+(R)\quad(\text{開上半球面・微分同相})$$

- 引き戻し計量 $g_{\mu\nu}=\dfrac{R^2}{\ell^2}\!\left(\delta_{\mu\nu}-\dfrac{x_\mu x_\nu}{\ell^2}\right)$（$\ell=\sqrt{R^2+|x|^2}$）は $G_{\mu\nu}+\Lambda g_{\mu\nu}=0$（$\Lambda=\tfrac{(n-1)(n-2)}{2R^2}$）を満たし、**de Sitter 空間の Beltrami 座標と内在的に一致**。$R\to\infty$ で平坦（Minkowski）に退化。
- **核心**：$\sigma_R$（非単射・全球）と $\Phi_R$（単射・開上半球のみ）の対比。

## 3. 多軸モデル

背景 $\mathbb{R}^{n+1}$ の $n+1$ 本の軸のどれを投影中心に選んでも、計量・曲率テンソルの構造は同一（**軸の対等性**）。$n+1$ 通りの主観空間が同時に立ち、相互変換は $T_{A\to B}=\Phi_B^{-1}\circ\Phi_A$。観測者の立場により「投影中心軸（直接アクセス不可）」と「主観座標軸」の役割が交換される。

## 4. 物理的舞台への一瞥（観察のみ）

外部を持たない真空（内在幾何系）の中心投影として、等質な4次元球面 $S^4(R_{\mathrm U})$（断面曲率 $1/R_{\mathrm U}^2$、$dS_4$ のユークリッド版）を取ると、その上に「広がり位相を持つ局所モード（粒子様状態）」$P_a=(\boldsymbol{\nu}_a,R_a)$、広がり $W_a=2R_a$ を置けることが観察される。ローレンツ計量・因果構造・物理的導出は主張しない。［Paper 14］

## 参考（Concept DOI）

- 球面投影 $\sigma_R$ の定義（基礎写像）：[10.5281/zenodo.20462569](https://doi.org/10.5281/zenodo.20462569)
- 中心投影による4次元空間の幾何学的定式化：[10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)
- 中心投影の合成演算（多軸・合成曲率半径の閉形式）：[10.5281/zenodo.20060728](https://doi.org/10.5281/zenodo.20060728)
- ［Paper 14］真空宇宙の中心投影と、広がり位相を持つ粒子様状態：[10.5281/zenodo.20543044](https://doi.org/10.5281/zenodo.20543044)

GitHub: github.com/WurabeSeiji/ai-chat-logs-open ／ 全文は Zenodo に公開（日英対訳・CC BY 4.0）
