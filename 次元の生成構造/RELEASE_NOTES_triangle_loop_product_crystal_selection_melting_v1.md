# Release Notes: 三次元正方配列の三角形ループ積による四元数・八元数結晶の選択、カスケード融解、および視野角スペクトルの二極分裂（速報）

## DOI
- Concept DOI: https://doi.org/10.5281/zenodo.24614040
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.24614041
- Zenodo record: https://zenodo.org/records/24614041

## Zenn
- 記事: `articles/triangle-loop-product-crystal-selection-melting.md`

## note
- 日本語版: https://note.com/kiharanoriaki/n/n6a1c65eb049b
  （元原稿: `次元の生成構造/三次元正方配列の積と自己無撞着_20261127/note_triangle_loop_ja_v1.md`）
- 英語版: https://note.com/kiharanoriaki/n/na8f7739b447d
  （元原稿: `次元の生成構造/三次元正方配列の積と自己無撞着_20261127/note_triangle_loop_en_v1.md`）

## Facebook
- 日本語版: https://www.facebook.com/kihara.noriaki/posts/pfbid02sQT1LuEaH7iJq2gvUrRy1zA467jCyDvyBs8YPLJqhkFG8aAzHghvLxHLJzVvrpQcl
  （元原稿: `次元の生成構造/三次元正方配列の積と自己無撞着_20261127/fb_message_triangle_loop_ja.md`）
- 英語版: https://www.facebook.com/kihara.noriaki/posts/pfbid02m9eqfe1GEYbnEEyehDckaF6Wa3rhCgQJxEwvbcTCZ9X9UvmLggTgswWzTDrc3Za8l
  （元原稿: `次元の生成構造/三次元正方配列の積と自己無撞着_20261127/fb_message_triangle_loop_en.md`）

## 公開ファイル（Zenodo, 9点）
- `三次元正方配列の三角形ループ積による結晶選択と融解_速報.md` / `triangle_loop_product_crystal_selection_melting_rapid_en.md`
- `三次元正方配列の三角形ループ積による結晶選択と融解_速報.tex` / `triangle_loop_product_crystal_selection_melting_rapid_en.tex`
- `三次元正方配列の三角形ループ積による結晶選択と融解_速報.pdf` / `triangle_loop_product_crystal_selection_melting_rapid_en.pdf`
- `build_tex.py`
- `figures_and_programs.zip`（全 18 本の検証スクリプト、`make_figures_en.py`、および `results/` 内の TXT/CSV/NPZ/日英 PNG 全データ）
- `SHA256SUMS`

## 位置づけ

背景時空・外部演算子・乱数・割り算を一切置かない三次元正方配列上の**三角形ループ積** $\mathrm{Tri}(T)_{ijk} = \sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}$ に関する速報論文。2026-11-25 公開の V. Chandrasekaran 論文（arXiv:2611.34924v1）との厳密な対比を含め、複素数数値インフレーション実験の八元数・高次元正方配列への拡張の中間成果を、全再現コード・データ・図表とともに公開した。

## 主な内容

1. **ケーリー＝ディクソン代数の自発的選択（厳密整数計算）**: 四元数（$N=4$）と八元数（$N=8$）のみが厳密周期 2 結晶（$U^2=I$）として成立（$\mathbb{C}$ は step 1 で消滅、$\mathbb{S}$ は欠陥なしで step 36 自壊）。
2. **八元数全 1,024 通りの 1 マス欠陥スペクトル**: ゼロ消滅 0 件（100% 生存）。実数元ゼロ 448 件は 8 マス複素回転行列ペアへ析出、虚数元ゼロ 448 件は厳密 336 : 112（3 : 1）で 32 マス四元数中間相と 8 マス複素数固定点（いずれも $\sum U^2 = 0$、パウリ射影 $I \pm \sigma_y$）へ分岐。
3. **8 マス溶け残りの厳密解析定理と 20,000 ステップ追跡**: $z_j' = (z_0^2 + z_1^2) z_j$、保存量 $z_1/z_0$、$S^1$ 上の 3 倍角回転、外部 504 マス遮断。虚数 32 マス中間相は $8 \to 4 \to 2$ 次元の 2 段階カスケード相転移を経て step 100 で完全飽和し step 20,000 まで不変。
4. **結合大結晶 $N = 64$（$O_8 \otimes O_8$）の 512 マス複合溶け残り**: step 70〜500 安定維持、64 要素の 3 段階離散エネルギー階層。
5. **三角形視野角スペクトルの二極分裂**: 背景方向の幾何級数的急落（$\Delta\ln > 687$、因果分離）と粒子内部の V 字再凝縮（$\rho \sim 1$）が $E_{512} \to 1$ と完全同期。

## 変更履歴

| 日付 | 内容 |
|---|---|
| 2026-11-28 | v1.0 公開。Zenodo 9 ファイル（日英 md/tex/pdf、`build_tex.py`、`figures_and_programs.zip`、`SHA256SUMS`）、DOI 取得、Zenn 記事原稿作成、PAPERS 目次更新 |
| 2026-11-28 | note 日本語版公開（https://note.com/kiharanoriaki/n/n6a1c65eb049b）、ヘッダー画像・note 英語版原稿作成 |
| 2026-11-28 | note 英語版公開（https://note.com/kiharanoriaki/n/na8f7739b447d） |
| 2026-11-28 | Facebook 日本語版公開（https://www.facebook.com/kihara.noriaki/posts/pfbid02sQT1LuEaH7iJq2gvUrRy1zA467jCyDvyBs8YPLJqhkFG8aAzHghvLxHLJzVvrpQcl）、Facebook 英語版原稿作成 |
| 2026-11-28 | Facebook 英語版公開（https://www.facebook.com/kihara.noriaki/posts/pfbid02m9eqfe1GEYbnEEyehDckaF6Wa3rhCgQJxEwvbcTCZ9X9UvmLggTgswWzTDrc3Za8l）、X 日本語版原稿作成 |
