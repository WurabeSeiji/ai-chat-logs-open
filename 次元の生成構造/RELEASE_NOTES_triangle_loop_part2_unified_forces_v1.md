# リリースノート：三次元正方配列の三角形ループ積における重力・電磁気・QCD・ディラックの統一接続仮説（Part II v1）

- **公開日**：2026-11-28
- **著者**：木原 範昭（株式会社WFシステム、ORCID: `0009-0004-6753-4020`）
- **Zenodo URL**：https://zenodo.org/records/24622566
- **Concept DOI（全版共通）**：`10.5281/zenodo.24622565`（https://doi.org/10.5281/zenodo.24622565）
- **Version DOI（v1）**：`10.5281/zenodo.24622566`（https://doi.org/10.5281/zenodo.24622566）
- **解説記事（note.com）**：
  - 日本語版：https://note.com/kiharanoriaki/n/n09cbc13c6be5
  - 英語版：https://note.com/kiharanoriaki/n/n0b72e078fd3b
- **接続する先行論文**：
  - Part I（統合完成論文：速報 1–7 統合）：`10.5281/zenodo.24622201`（Concept: `10.5281/zenodo.24622200`）
  - 八元数 QCD 閉じ込め・質量ギャップ論文：`10.5281/zenodo.24369465`
  - 八元数ヤン＝ミルズ・ディラック論文：`10.5281/zenodo.24387753`（note 解説：https://note.com/kiharanoriaki/n/nb44bbd7a083a）

## 概要

Part I（仮説ゼロの確立結果）の延長上に、**外部時空と内部ゲージの直積分離仮説**（$O_8^{\mathrm{ext}} \otimes O_8^{\mathrm{int}}$）を明示的に導入し、重力（$SO(3)$ 外部計量歪み）・電磁気（$U(1)_Q$ 可換マクスウェル極限）・強い力 QCD（$SU(3)$ 非可換ヤン＝ミルズ極限）・ディラック物質場（$Spin(8)$ Triality と CAR フェルミオン `1:3:3:1`）を**単一の三角形ループ積 $\mathrm{Tri}(T)$ の枠組みで統一接続する仮説論文（Part II）**である。

数値・代数検証（`tri_loop_part2_unified_forces_verify_v0.py`）により以下を確認した：
1. **外部・内部の完全可換分離**：$[G_{\mathrm{ext}}(SO(3))\otimes I, I\otimes G_{\mathrm{int}}(SU(3)\times U(1)_Q)] = 0$（全 27 組、最大誤差 `0.000e+00`）。
2. **マクスウェルとヤン＝ミルズの同時分離創発**：$e_7$ 可換三重欠損 `0.000e+00`（マクスウェル）に対し、$SU(3)$ 非可換プラケットは $O(\varepsilon^4)$ ヤン＝ミルズ作用へ収束（`0.9955 -> 0.99993`）し、ビアンキ／ヤコビ恒等式 56 組すべてが `0.000e+00` で成立。
3. **ディラック分散と三重読み出し**：ディラック分散誤差 `2.591e-15`、三角形ループ積によるディラック三重読み出し誤差 `1.110e-16`。
4. **標準模型と整合する選択的相互作用**：$\nu$（EM=0, QCD=0）、$\bar{d}$（EM=0.243, QCD=0.504）、$u$（EM=0.482, QCD=0.580）、$e^+$（EM=0.714, QCD=0）の 4 セクター分離を実証。

あわせて、確立済みの代数・数値事実と、今後の連続極限・微細構造定数・アインシュタイン方程式導出など未解決の仮説領域を厳密に切り分けている。

## 公開ファイル（全 9 点、Zenodo 登録済み）

作業ディレクトリ：`次元の生成構造/三次元正方配列の積と自己無撞着_20261127/`

- `三次元正方配列の三角形ループ積における重力・電磁気・QCD・ディラックの統一接続仮説_Part2.pdf`（日本語版 6 ページ、2 図付き）
- `triangle_loop_unified_gravity_em_qcd_dirac_hypothesis_part2_en.pdf`（英語版 5 ページ、2 図付き）
- `三次元正方配列の三角形ループ積における重力・電磁気・QCD・ディラックの統一接続仮説_Part2.tex` / `.md`
- `triangle_loop_unified_gravity_em_qcd_dirac_hypothesis_part2_en.tex` / `.md`
- `build_tex_part2.py`
- `figures_and_programs_part2.zip`（検証スクリプト・作図スクリプト・ログ・日英図一式）
- `SHA256SUMS_part2`
