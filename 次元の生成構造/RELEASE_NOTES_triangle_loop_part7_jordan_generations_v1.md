# リリースノート：三次元正方配列の三角形ループ積における例外ジョルダン代数 $J_3(\mathbb{O})$ 結晶、3 世代トライアリティ、および素の結合定数幾何統一（Part VII v1）

- **公開日**：2026-11-28
- **著者**：木原 範昭（株式会社WFシステム、ORCID: `0009-0004-6753-4020`）
- **Zenodo URL**：https://zenodo.org/records/24628777
- **Concept DOI（全版共通）**：`10.5281/zenodo.24628776`（https://doi.org/10.5281/zenodo.24628776）
- **Version DOI（v1）**：`10.5281/zenodo.24628777`（https://doi.org/10.5281/zenodo.24628777）
- **Zenn 記事**：`articles/triangle-loop-jordan-three-generations-part7.md`
- **接続する先行論文**：
  - Part I（統合完成論文）：`10.5281/zenodo.24622200`（v1: `24622201`）
  - Part II（重力・電磁気・QCD・ディラック統一接続仮説）：`10.5281/zenodo.24622565`（v1: `24622566`）
  - Part III（時空・内部空間読出し接続仮説）：`10.5281/zenodo.24625962`（v1: `24625963`）
  - Part IV（電弱対称性破れとパリティ最大破れ）：`10.5281/zenodo.24626592`（v1: `24626593`）
  - Part V（多殻大格子重力・因果伝播・スピン 2 $g_{00}(d)$ 解禁）：`10.5281/zenodo.24627519`（v1: `24627520`）
  - Part VI（2 粒子動的散乱・空間分離 CHSH ベル不等式動的破れ）：`10.5281/zenodo.24628297`（v1: `24628298`）

## 概要

27 次元例外ジョルダン代数 $J_3(\mathbb{O})$ と三角形ループ積 $\mathrm{Tri}(T)$ の厳密接続を、乱数なしの決定論的計算（`tri_loop_part7_jordan_generations_v0.py`）により検証した論文（Part VII）である。

1. **非対角 3 世代巡回積 $J_{\mathrm{off}}$（24D）の $\mathrm{Tri}(T)$ 厳密結晶定理（J3-1）**：$\mathrm{Tri}(J_{\mathrm{off}}) = -\frac{1}{64} J_{\mathrm{off}}$（不動点残差 `0.000e+00`、100 ステップ後も `0.000e+00`、非対角摂動に対し `3.536e-05` で自己修復）。
2. **$F_4 \supset SO(8) \supset G_2 \supset SU(3)_C$ 連鎖と 3 世代トライアリティ（J3-2）**：SVD 核空間から $\mathfrak{f}_4$（52D）・$\mathfrak{so}(8)$（28D）・$\mathfrak{so}(7)_1$（21D、固定次元 `[1,0,0]` → $8_v\oplus 8_s\oplus 8_c$）・$\mathfrak{g}_2$（14D）・$SU(3)_C$（8D）を抽出し、全 3 世代で $SU(3)_C$ カシミール固有値 `[0, 0, 8/9×6]` が完全一致することを確認。
3. **3 世代非線形階層増幅と小出比 $Q$ の反証境界（J3-3）**：2 段階局在破れにより `500.4 : 17.2 : 1.0` や `2138.6 : 16.0 : 1.0` の指数階層が自発増幅される一方、小出比 $Q$ は初期比に依存し $2/3$ へ自動ロックしないことを明示。
4. **素の結合定数幾何統一比（J3-4）**：64D 表現上で $\alpha_3^{\mathrm{bare}} : \alpha_2^{\mathrm{bare}} : \alpha_1^{\mathrm{bare}} = 1:1:1$、$\sin^2\theta_W^{\mathrm{bare}} = \alpha_{\mathrm{em}}^{\mathrm{bare}}/\alpha_s^{\mathrm{bare}} = 3/8$（27D 3 世代 $SU(3)_C$ トレースも `1:1:1` 均等）。低エネルギー実測値との差は RG 走行境界として切り分け。

## 公開ファイル（全 9 点、Zenodo 登録済み）

作業ディレクトリ：`次元の生成構造/三次元正方配列の積と自己無撞着_20261127/`

- `三次元正方配列の三角形ループ積における例外ジョルダン代数J3Oと3世代構造_Part7.pdf`（日本語版、2 図付き）
- `triangle_loop_jordan_three_generations_part7_en.pdf`（英語版、2 図付き）
- `三次元正方配列の三角形ループ積における例外ジョルダン代数J3Oと3世代構造_Part7.tex` / `.md`
- `triangle_loop_jordan_three_generations_part7_en.tex` / `.md`
- `build_tex_part7.py`
- `figures_and_programs_part7.zip`（検証スクリプト・作図スクリプト・ログ・日英図一式）
- `SHA256SUMS_part7`
