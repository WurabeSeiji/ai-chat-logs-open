# リリースノート：三次元正方配列の三角形ループ積における電弱対称性破れとパリティ最大破れ（Part IV v1）

- **公開日**：2026-11-28
- **著者**：木原 範昭（株式会社WFシステム、ORCID: `0009-0004-6753-4020`）
- **Zenodo URL**：https://zenodo.org/records/24626593
- **Concept DOI（全版共通）**：`10.5281/zenodo.24626592`（https://doi.org/10.5281/zenodo.24626592）
- **Version DOI（v1）**：`10.5281/zenodo.24626593`（https://doi.org/10.5281/zenodo.24626593）
- **Zenn 記事**：`articles/triangle-loop-electroweak-symmetry-breaking-parity-part4.md`
- **解説記事（note.com）**：Part I〜VII 総括（2026-11-28 公開）
  - 日本語版：https://note.com/kiharanoriaki/n/n6f27afeb5795
  - 英語版：https://note.com/kiharanoriaki/n/nd11ac852843b
- **Facebook 投稿**：Part I〜VII 総括
  - 日本語版：https://www.facebook.com/kihara.noriaki/posts/pfbid0QX3FBTuwMpv2W3FvfC92BsaDpEsa45ERJuTyxk3trJa4wFQ5NAK4iSur2VPFEwz1l
  - 英語版：https://www.facebook.com/kihara.noriaki/posts/pfbid0Q17gHw7gzCko7svW61xx4NW3ernhnqFeKhybhKNSNV5jhtRqcR9LM9fqNMnABu3Tl
- **X（Twitter）投稿**：Part I〜VII 総括
  - 日本語版：https://x.com/NoriakiKihara/status/2126693283250532491
  - 英語版：https://x.com/NoriakiKihara/status/2126693505334657414
- **接続する先行論文**：
  - Part I（統合完成論文：速報 1–7 統合）：`10.5281/zenodo.24622200`（v1: `24622201`）
  - Part II（重力・電磁気・QCD・ディラック統一接続仮説）：`10.5281/zenodo.24622565`（v1: `24622566`）
  - Part III（時空・内部空間読出し接続仮説）：`10.5281/zenodo.24625962`（v1: `24625963`）

## 概要

Part III で特定した外側 $O_8^{\mathrm{ext}}$ の `1+3+4` 分割（純粋 4 空間 $SU(2)$）を起点に、標準模型 1 世代の電弱セクター（$SU(2)_L \times U(1)_Y \to U(1)_{\mathrm{em}}$、パリティ最大破れ、$\mathrm{Tri}(T)$ 動的ヒッグス質量生成、素のワインバーグ角）を決定論的検証スクリプト（`tri_loop_part4_electroweak_parity_v0.py`、乱数なし）で直接検証した論文（Part IV）である。

1. **パリティ最大破れの厳密な代数的起源（EW-1）**：$G_2 \supset SO(4)$ において 3 空間不変な内部 $SU(2)_{\mathrm{pure4}}$ は自己双対（+カイラル）`1.000000000000`・反自己双対（−カイラル）`0.000000000000` に厳密一致（$G_2$ 微分残差・リー代数残差とも `0.000e+00`）。逆カイラル側は 100% 外部空間回転 $SO(3)$ にロックされる。
2. **標準模型全 12 生成子の厳密可換性とパリティ選択則（EW-2）**：$Y = I_8\otimes Q_{\mathrm{int}} - \frac{1}{2}P_4\otimes I_8$ により $[SO(3)_{\mathrm{ext}}, \mathrm{SM}_{12}]$（36 組）・$[SU(3)_C, SU(2)_L]$（24 組）・$[SU(3)_C, U(1)_Y]$（8 組）・$[SU(2)_L, U(1)_Y]$（3 組）すべて誤差 `0.000e+00`。左巻きダブレット（$e_4$ 側）は `0.323572` で結合、右巻きシングレット（$e_0, e_1$ 側）は `0.000e+00`。
3. **$\mathrm{Tri}(T)$ 実動力学によるヒッグス凝縮（EW-3a,b）**：サイト $A=4$ 局在 P4$_4$（$\epsilon_H = 0.05 \to 0.22, \tau=40$）で $W_1, W_2, W_3$ の破れノルムが `0.012593 → 0.035517` へ 3 成分 6 桁完全縮退で立ち上がる一方、$e_4$ 固定の残余 8D 対称群は `3.82e-14`（無質量）を維持。
4. **素のワインバーグ角 $\sin^2\theta_W^{\mathrm{bare}} = 3/8$ と 4×4 質量スペクトル（EW-3c, EW-5）**：64D トレース比だけで $\sin^2\theta_W^{\mathrm{bare}} = 8/(64/3) = 3/8 = 0.375000$、複素ヒッグス二重項の 4×4 質量行列固有値 $[m_\gamma^2, m_{W^+}^2, m_{W^-}^2, m_Z^2] = [-2.082\times 10^{-17}, 0.250000, 0.250000, 0.400000]$（$m_W^2/m_Z^2 = 5/8$）。

## 公開ファイル（全 9 点、Zenodo 登録済み）

作業ディレクトリ：`次元の生成構造/三次元正方配列の積と自己無撞着_20261127/`

- `三次元正方配列の三角形ループ積における電弱対称性破れとパリティ最大破れ_Part4.pdf`（日本語版 6 ページ、2 図付き）
- `triangle_loop_electroweak_symmetry_breaking_parity_part4_en.pdf`（英語版 5 ページ、2 図付き）
- `三次元正方配列の三角形ループ積における電弱対称性破れとパリティ最大破れ_Part4.tex` / `.md`
- `triangle_loop_electroweak_symmetry_breaking_parity_part4_en.tex` / `.md`
- `build_tex_part4.py`
- `figures_and_programs_part4.zip`（検証スクリプト・作図スクリプト・ログ・日英図一式）
- `SHA256SUMS_part4`
