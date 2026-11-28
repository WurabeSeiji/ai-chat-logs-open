---
title: "三次元正方配列の三角形ループ積における電弱対称性破れとパリティ最大破れ（Part IV）"
emoji: "⚡"
type: "idea"
topics: ["physics", "octonion", "standardmodel", "quantum", "mathematics"]
published: true
---

三角形ループ積シリーズ Part I（Concept DOI: `10.5281/zenodo.24622200`）・Part II（Concept DOI: `10.5281/zenodo.24622565`）・Part III（Concept DOI: `10.5281/zenodo.24625962`）の続編として、**電弱対称性破れとパリティ最大破れ（Part IV）**を Zenodo に公開しました。

Part III で特定した外側 $O_8^{\mathrm{ext}}$ の `1+3+4` 分割（余り 4 空間 $\{e_4,e_5,e_6,e_7\}$ 上の純粋 3D $SU(2)$）を起点に、標準模型 1 世代の電弱セクター（$SU(2)_L \times U(1)_Y \to U(1)_{\mathrm{em}}$、パリティ最大破れ、$\mathrm{Tri}(T)$ 動的ヒッグス質量生成、素のワインバーグ角）を乱数なしの決定論的計算（`tri_loop_part4_electroweak_parity_v0.py`）で直接検証しました。

## 論文（Zenodo・日英 PDF・再現コード一式）

- **Concept DOI（常に最新版を指します）**: https://doi.org/10.5281/zenodo.24626592
- **Version DOI（v1.0）**: https://doi.org/10.5281/zenodo.24626593
- **日本語 PDF**: https://zenodo.org/records/24626593/files/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E4%B8%89%E8%A7%92%E5%BD%A2%E3%83%AB%E3%83%BC%E3%83%97%E7%A9%8D%E3%81%AB%E3%81%8A%E3%81%91%E3%82%8B%E9%9B%BB%E5%BC%B1%E5%AF%BE%E7%A7%B0%E6%80%A7%E7%A0%B4%E3%82%8C%E3%81%A8%E3%83%91%E3%83%AA%E3%83%86%E3%82%A3%E6%9C%80%E5%A4%A7%E7%A0%B4%E3%82%8C_Part4.pdf
- **English PDF**: https://zenodo.org/records/24626593/files/triangle_loop_electroweak_symmetry_breaking_parity_part4_en.pdf
- **関連リポジトリ**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E7%A9%8D%E3%81%A8%E8%87%AA%E5%B7%B1%E7%84%A1%E6%92%9E%E7%9D%80_20261127

## 主な検証結果（乱数なし・機械精度）

1. **パリティ最大破れの厳密な代数的起源（$G_2$ カイラル 100%:0% 定理）**
   4 空間 $\{e_4,e_5,e_6,e_7\}$ の $SO(4) \cong SU(2)_{\mathrm{SD}} \times SU(2)_{\mathrm{ASD}}$ に対し八元数 $G_2$（14D）の微分制約を課すと、3 空間不変な内部 $SU(2)_{\mathrm{pure4}}$ は **自己双対（+カイラル）`1.000000000000`・反自己双対（−カイラル）`0.000000000000`** に厳密一致します（$G_2$ 微分残差・リー代数残差とも `0.000e+00`）。逆カイラル側は 100% 外部空間回転 $SO(3)$ にロックされ独立内部ゲージになれません。
2. **標準模型全 12 ゲージ生成子の厳密可換性とパリティ選択則**
   $Y = I_8\otimes Q_{\mathrm{int}} - \frac{1}{2}P_4\otimes I_8$（$Q_{\mathrm{em}} = T_3^L + Y$）により、$[SO(3)_{\mathrm{ext}}, \mathrm{SM}_{12}]$（36 組）・$[SU(3)_C, SU(2)_L]$（24 組）・$[SU(3)_C, U(1)_Y]$（8 組）・$[SU(2)_L, U(1)_Y]$（3 組）のすべてが誤差 `0.000e+00` で可換となります。$W_1$ 回転では左巻きダブレット（$e_4$ 側）が `0.323572` で結合し、右巻きシングレット（$e_0, e_1$ 側）は `0.000e+00`（厳密非結合）です。
3. **$\mathrm{Tri}(T)$ 実動力学によるヒッグス凝縮と $W$ 質量生成 vs 残余無質量**
   真空 $H$ にサイト $A=4$ 局在のヒッグス凝縮体 P4$_4$（$\epsilon_H = 0.05 \to 0.22, \tau=40$）を $\mathrm{Tri}(T)$ で発展させると、$W_1, W_2, W_3$ の対称性破れノルムが `0.012593 → 0.025105 → 0.035323 → 0.035517` へ 3 成分 6 桁完全縮退で立ち上がる一方、$e_4$ 固定の残余 8D 対称群は全振幅で `3.82e-14`（厳密無質量）を維持します。
4. **素のワインバーグ角 $\sin^2\theta_W^{\mathrm{bare}} = 3/8$ と 4×4 電弱質量スペクトル**
   $SU(5)/SO(10)$ を仮定せず 64D トレース比だけで $\mathrm{Tr}(T_3^2)=8, \mathrm{Tr}(Y^2)=40/3, \mathrm{Tr}(Q_{\mathrm{em}}^2)=64/3$ となり、**$\sin^2\theta_W^{\mathrm{bare}} = 3/8 = 0.375000$** が調整パラメータゼロで導出されます。この幾何比 $g'/g = \sqrt{3/5}$ を複素ヒッグス二重項の 4×4 質量行列 $M^2$ に適用すると、固有値 $[m_\gamma^2, m_{W^+}^2, m_{W^-}^2, m_Z^2] = [-2.082\times 10^{-17}, 0.250000, 0.250000, 0.400000]$・$m_W^2/m_Z^2 = 5/8 = 0.625000$ が厳密に得られます。

## 残課題の切り分け

低エネルギー実測値（$\sin^2\theta_W \approx 0.231$）への繰り込み群走行、多殻大格子 $O_8^{\otimes K}$ の連続極限（Part V 予定）、2 粒子動的散乱・CHSH 動的破れ（Part VI 予定）、3 世代構造と精密質量比・結合定数（Part VII 予定）は未検証課題として明確に切り分けています。
