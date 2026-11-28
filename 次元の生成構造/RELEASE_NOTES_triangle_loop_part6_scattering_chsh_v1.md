# リリースノート：三次元正方配列の三角形ループ積における 2 粒子動的散乱、ゲージ媒介相互作用、および空間分離 CHSH ベル不等式の動的破れ（Part VI v1）

- **公開日**：2026-11-28
- **著者**：木原 範昭（株式会社WFシステム、ORCID: `0009-0004-6753-4020`）
- **Zenodo URL**：https://zenodo.org/records/24628298
- **Concept DOI（全版共通）**：`10.5281/zenodo.24628297`（https://doi.org/10.5281/zenodo.24628297）
- **Version DOI（v1）**：`10.5281/zenodo.24628298`（https://doi.org/10.5281/zenodo.24628298）
- **Zenn 記事**：`articles/triangle-loop-two-particle-scattering-chsh-part6.md`
- **解説記事（note.com）**：Part I〜VII 総括（2026-11-28 公開）
  - 日本語版：https://note.com/kiharanoriaki/n/n6f27afeb5795
  - 英語版：https://note.com/kiharanoriaki/n/nd11ac852843b
- **Facebook 投稿**：Part I〜VII 総括
  - 日本語版：https://www.facebook.com/kihara.noriaki/posts/pfbid0QX3FBTuwMpv2W3FvfC92BsaDpEsa45ERJuTyxk3trJa4wFQ5NAK4iSur2VPFEwz1l
  - 英語版：https://www.facebook.com/kihara.noriaki/posts/pfbid0Q17gHw7gzCko7svW61xx4NW3ernhnqFeKhybhKNSNV5jhtRqcR9LM9fqNMnABu3Tl
- **接続する先行論文**：
  - Part I（統合完成論文）：`10.5281/zenodo.24622200`（v1: `24622201`）
  - Part II（重力・電磁気・QCD・ディラック統一接続仮説）：`10.5281/zenodo.24622565`（v1: `24622566`）
  - Part III（時空・内部空間読出し接続仮説）：`10.5281/zenodo.24625962`（v1: `24625963`）
  - Part IV（電弱対称性破れとパリティ最大破れ）：`10.5281/zenodo.24626592`（v1: `24626593`）
  - Part V（多殻大格子重力・因果伝播・スピン 2 $g_{00}(d)$ 解禁）：`10.5281/zenodo.24627519`（v1: `24627520`）

## 概要

真空 $H = O_8 \otimes M_8$（$N = 64$）上に空間分離した 2 粒子と媒介波束を配置し、非線形反復 $\mathrm{Tri}(T)$ のもとでの 2 粒子重力合算・ゲージ媒介電荷分離・動的散乱・空間分離 CHSH ベル不等式破れを決定論的計算（`tri_loop_part6_scattering_chsh_v0.py`, `tri_loop_part6_chsh_and_coulomb_v0.py`、乱数なし）で実証した論文（Part VI）である。

1. **2 粒子重力合算崩壊（SC-3）**：$\epsilon_c^{(1)} = 0.273580 \to \epsilon_c^{(2)} = 0.269267$（比 `0.9842`、同符号・異符号 6 桁一致）。
2. **局所 $U(1)$ ゲージ冗長性とゲージリンク電荷分離（SC-1, SC-1b）**：ゲージ場なしの局所 $U(1)$ 回転差は厳密 `0.000e+00`。共有サイト $C=2$ 経由の $G_{02}\pm G_{12}$ 媒介で同種 vs 異種が明確分裂（$g=0.002$: `0.2311` vs `0.0777`）。
3. **遠隔スピン 2 ($\sigma_x$) 因果散乱と厳密二乗則（SC-2, SC-2b）**：$\tau=1$ で全振幅厳密 `0.0000 rad`、$\tau=2$ 到着時に $\delta\theta(\tau=2) = -1.4940\,\epsilon_{\mathrm{s2}}^2$（4 桁一定）、$\tau=8$ で最大 `-2.289 rad`。
4. **空間分離 2 サイト CHSH 動的破れと自己修復アトラクタ（SC-4b）**：真空全対で $S_{\mathrm{CHSH}} = 2.828358 \approx 2\sqrt{2}$（Horodecki 最大値と 4 単位ベクトル直接評価が 6 桁一致）。積状態投入（`2.1462`）からも 3 ステップで `2.8259` へ動的回復。
5. **地平面崩壊に伴う CHSH 不連続破壊（SC-5b）**：$\epsilon_{\mathrm{mass}} \le 0.27$ で $S=2.8284$ 維持、$\epsilon_{\mathrm{mass}} \ge 0.28$ で `purity = 0.2500`, $S = 0.0000$ へ不連続崩壊。

## 公開ファイル（全 9 点、Zenodo 登録済み）

作業ディレクトリ：`次元の生成構造/三次元正方配列の積と自己無撞着_20261127/`

- `三次元正方配列の三角形ループ積における2粒子動的散乱とCHSHベル不等式の動的破れ_Part6.pdf`（日本語版、2 図付き）
- `triangle_loop_two_particle_scattering_chsh_part6_en.pdf`（英語版、2 図付き）
- `三次元正方配列の三角形ループ積における2粒子動的散乱とCHSHベル不等式の動的破れ_Part6.tex` / `.md`
- `triangle_loop_two_particle_scattering_chsh_part6_en.tex` / `.md`
- `build_tex_part6.py`
- `figures_and_programs_part6.zip`（検証スクリプト 2 本・作図スクリプト・ログ 2 本・日英図一式）
- `SHA256SUMS_part6`
