# リリースノート：三次元正方配列の三角形ループ積における Hurwitz 4 段階局在質量階層と直積準同型指数則（第 6 速報 v1）

- **公開日**：2026-11-28
- **著者**：木原 範昭（株式会社WFシステム、ORCID: `0009-0004-6753-4020`）
- **Zenodo URL**：https://zenodo.org/records/24621454
- **Concept DOI（全版共通）**：`10.5281/zenodo.24621453`（https://doi.org/10.5281/zenodo.24621453）
- **Version DOI（v1）**：`10.5281/zenodo.24621454`（https://doi.org/10.5281/zenodo.24621454）
- **先行速報**：
  - 第 1 速報 v1：`10.5281/zenodo.24616931`
  - 第 2 速報 v1：`10.5281/zenodo.24618400`
  - 第 3 速報 v1：`10.5281/zenodo.24619357`
  - 第 4 速報 v1：`10.5281/zenodo.24620262`
  - 第 5 速報 v1：`10.5281/zenodo.24621139`

## 概要（第 6 速報の確定成果）

第 1 層時間結晶 $H = O_8 \otimes M_8$ の上に外側 $O_8$ の全 4 種の Hurwitz 部分代数（$\mathbb{O}_8, \mathbb{H}_4, \mathbb{C}_2, \mathbb{R}_1$、余次元 $m = \log_2(8/d) \in \{0, 1, 2, 3\}$）を配置して決定論的計測を行い、1 層内の離散質量階層と $K$ 段入れ子による指数関数的階層（$10^5 \sim 10^{16}$）の 4 つの厳密法則を確定した。

1. **衣／芯エネルギー比のメルセンヌ公式**：弱結合極限で $E_{\mathrm{core}}/E_{\mathrm{tot}} \to d/8$、$E_{\mathrm{cloud}}/E_{\mathrm{core}} \to (8-d)/d = 2^m - 1 \in \{0, 1, 3, 7\}$（実測 $\epsilon=0.05$: `0.000, 1.000, 2.988, 6.682`）。
2. **総エネルギー $d/8$ 則およびコアエネルギー $(d/8)^2$ 則**：$\Delta\mathrm{nz} = 64d$、$E_{\mathrm{tot}} \propto d/8$（`1 : 1/2 : 1/4 : 1/8`）、$E_{\mathrm{core}} \propto (d/8)^2$（`1 : 1/4 : 1/16 : 1/64`）。
3. **空間倍音のパスカル二項則**：$\mathbb{Z}_2^3$ 上の空間倍音 $k = 1 \dots m$ の比が厳密に二項係数 $\binom{m}{k}$ に一致（`1`, `2:1`, `3:3:1`）、不確定性超過も `0.0000 → +0.0732 → +0.1759 → +0.3044` と単調増加。
4. **直積準同型定理による桁違いの指数関数的階層 $r^{-K}$**：$\mathrm{Tri}(A \otimes B) = \mathrm{Tri}(A) \otimes \mathrm{Tri}(B)$（相対誤差 `2.855e-16`）により、コアエネルギー比 $1/64$ はわずか $K = 3$ 段で $2.6 \times 10^5$（電子／トップ比級）、$K = 9$ 段で $1.8 \times 10^{16}$（電弱／プランク比級）に到達する。

## 公開ファイル（全 9 点、Zenodo 登録済み）

作業ディレクトリ：`次元の生成構造/三次元正方配列の積と自己無撞着_20261127/`

- `三次元正方配列の三角形ループ積におけるHurwitz 4段階局在質量階層と直積指数則_第6速報.pdf`（日本語版 5 ページ）
- `triangle_loop_hurwitz_mass_hierarchy_exponential_nesting_rapid6_en.pdf`（英語版 4 ページ）
- `三次元正方配列の三角形ループ積におけるHurwitz 4段階局在質量階層と直積指数則_第6速報.tex` / `.md`
- `triangle_loop_hurwitz_mass_hierarchy_exponential_nesting_rapid6_en.tex` / `.md`
- `build_tex_rapid6.py`
- `figures_and_programs_rapid6.zip`（検証・作図スクリプト、数値ログ、日英図 4 枚）
- `SHA256SUMS_rapid6`
