---
title: "三次元正方配列の三角形ループ積における 2 粒子動的散乱、ゲージ媒介相互作用、および空間分離 CHSH ベル不等式の動的破れ（Part VI）"
emoji: "⚛️"
type: "idea"
topics: ["physics", "quantum", "octonion", "bellinequality", "mathematics"]
published: true
---

三角形ループ積シリーズ Part I–V（Concept DOIs: `24622200`, `24622565`, `24625962`, `24626592`, `24627519`）の続編として、**2 粒子動的散乱、ゲージ媒介相互作用、および空間分離 CHSH ベル不等式の動的破れ（Part VI）**を Zenodo に公開しました。

真空 $H = O_8 \otimes M_8$（$N = 64$）上に空間分離した 2 粒子と媒介波束を配置し、非線形写像 $T(\tau+1) = \mathrm{Tri}(T(\tau))/\|\mathrm{Tri}(T(\tau))\|$ のもとでの 2 粒子重力合算・局所ゲージ不変性とゲージリンク電荷分離・遠隔スピン 2 散乱・空間分離 CHSH ベル不等式破れを、追加仮定ゼロ・乱数なしの決定論的計算（`tri_loop_part6_scattering_chsh_v0.py`, `tri_loop_part6_chsh_and_coulomb_v0.py`）で実証しました。

## 論文（Zenodo・日英 PDF・再現コード一式）

- **Concept DOI（常に最新版を指します）**: https://doi.org/10.5281/zenodo.24628297
- **Version DOI（v1.0）**: https://doi.org/10.5281/zenodo.24628298
- **日本語 PDF**: https://zenodo.org/records/24628298/files/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E4%B8%89%E8%A7%92%E5%BD%A2%E3%83%AB%E3%83%BC%E3%83%97%E7%A9%8D%E3%81%AB%E3%81%8A%E3%81%91%E3%82%8B2%E7%B2%92%E5%AD%90%E5%8B%95%E7%9A%84%E6%95%A3%E4%B9%B1%E3%81%A8CHSH%E3%83%99%E3%83%AB%E4%B8%8D%E7%AD%89%E5%BC%8F%E3%81%AE%E5%8B%95%E7%9A%84%E7%A0%B4%E3%82%8C_Part6.pdf
- **English PDF**: https://zenodo.org/records/24628298/files/triangle_loop_two_particle_scattering_chsh_part6_en.pdf
- **関連リポジトリ**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E7%A9%8D%E3%81%A8%E8%87%AA%E5%B7%B1%E7%84%A1%E6%92%9E%E7%9D%80_20261127

## 主な検証結果（乱数なし・機械精度）

1. **2 粒子重力合算崩壊（符号非依存の普遍引力）**
   1 粒子の臨界崩壊振幅 $\epsilon_c^{(1)} = 0.273580$ に対し、2 粒子では $\epsilon_c^{(2)} = 0.269267$（比 `0.9842`、1.58% 低下）。同符号 $(+,+)$ と異符号 $(+,-)$ が 6 桁完全一致します。
2. **局所 $U(1)$ ゲージ冗長性 vs ゲージリンク媒介電荷分離**
   ゲージ場なしの局所 $U(1)$ 位相回転だけでは同種と異種の差が厳密 `0.000e+00` となる一方、共有中間サイト $C=2$ を介した反対称 $U(1)$ ゲージリンク $G_{02}\pm G_{12}$ を導入すると、同種 $q=(+1,+1)$ と異種 $q=(+1,-1)$ で非線形相互作用強度 $|\Delta\|P\|_{\mathrm{int}}|$ が明確に分裂します（$g=0.002$: `0.2311` vs `0.0777`）。
3. **因果遅延を伴う遠隔スピン 2 ($\sigma_x$) 動的散乱と厳密二乗則**
   遠隔サイト $B=1$ のスピン 2 源による標的 $A=0$ の位相シフトは、$\tau=1$ で全振幅厳密 `0.0000 rad`（有限速度因果律）、$\tau=2$ 到着時に $\delta\theta(\tau=2) = -1.4940\,\epsilon_{\mathrm{s2}}^2$（4 桁一定の厳密二乗則＝源エネルギー比例）、$\tau=8$ で最大 `-2.289 rad` に成長します。
4. **空間分離 2 サイト（Alice=1, Bob=2）の CHSH ベル不等式動的破れと自己修復アトラクタ**
   真空 $H$ の任意の空間分離サイト対 $(1,2),(0,1),(1,4)$ は純粋度 `1.0000`・$S_E = 1.0000\text{ bit}$・**$S_{\mathrm{CHSH}} = 2.828358 \approx 2\sqrt{2}$**（Horodecki 最大値と具体的 4 単位ベクトル直接評価が 6 桁完全一致）を自発達成します。さらに初期積状態 $|\uparrow_1\uparrow_2\rangle$ を投入して $S_{\mathrm{CHSH}}(0) = 2.1462$ へ低下させても、$\mathrm{Tri}(T)$ が 3 ステップで `2.8259` へ動的に引き戻します。
5. **事象の地平面崩壊に伴う CHSH 破れの不連続破壊（$2\sqrt{2} \to 0$）**
   中心質量 $\epsilon_{\mathrm{mass}} \le 0.27$（弱結合共存相）では遠隔 2 サイトの $S_{\mathrm{CHSH}} = 2.8284$ が完全維持されますが、臨界点 $\epsilon_c \approx 0.2736$ を超えた $\epsilon_{\mathrm{mass}} \ge 0.28$ では `purity = 0.2500`（完全混合 $\mathbb{I}/4$）・$S_{\mathrm{CHSH}} = 0.0000$ へ不連続崩壊します。

## 残課題の切り分け

多殻 $K \ge 3$ での 2 粒子波束散乱断面積、および 3 世代構造（$J_3(\mathbb{O})$）・精密質量比・結合定数 $\alpha, \alpha_s$ の導出（Part VII 予定）は次段階として切り分けています。
