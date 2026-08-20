# Release Notes: 零閉包・有限位数・自己無撞着幾何からの対称性生成 —— 唯一の外部指定パラメータ N と、一般化・動力学導出という残された課題

## DOI
- Concept DOI: https://doi.org/10.5281/zenodo.22028072
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.22028073
- Zenodo record: https://zenodo.org/record/22028073

## Zenn
- 記事: https://zenn.dev/noriaki_kihara/articles/closure-axioms-symmetry-derivation
  （元原稿: `articles/closure-axioms-symmetry-derivation.md`）

## 公開ファイル（Zenodo, 6点）
- `closure_axioms_symmetry_derivation_ja_public_v1.0.md` / `closure_axioms_symmetry_derivation_en_public_v1.0.md`
- `closure_axioms_symmetry_derivation_ja_public_v1.0.tex` / `closure_axioms_symmetry_derivation_en_public_v1.0.tex`
- `closure_axioms_symmetry_derivation_ja_public_v1.0.pdf`（35ページ） / `closure_axioms_symmetry_derivation_en_public_v1.0.pdf`（29ページ）

## 位置づけ

「次元の生成構造」シリーズの公理系整理ノート。4公理（A1: 複素零閉包 Σxₙ²=0、A2: 有限位数 U^N=I、A3: simplex 整合、A4: 自己無撞着 X=F(X)）と唯一の外部離散パラメータ N から、幾何・対称性・統計・読出し構造がどこまで導かれるかを、厳密／条件付き厳密／自己論文で数値確認済み／未導出の4区分で帳簿化した。

自己引用は3本のみ（波の周期表 v2 / 波と場の二層分離 / ゼロ閉塞は4次元だった）。

## 主な内容

- **負符号の統一起源**: 時間・曲率・内部軸の負符号はすべて (iu)²=−u²（観測不能複素軸の実表示）から出る。Minkowski 計量・曲率場の外部導入は不要
- **複素5自由度の厳密導出**: 6レジスタ表現で dim_C C⁶ − 1 = 5
- **標準模型ゲージ群への二経路**:
  - 粗視化経路: Hermitian 3⊕2 分解保存＋全体位相除去 ⟹ S(U(3)×U(2)) ≅ [SU(3)×SU(2)×U(1)]/Z₆（条件付き厳密）
  - 細分化経路: Q²=Q₁²+Q₂²+Q₃² の内部三重項 ⟹ 内部 SU(3)、(R,Q₁,Q₂,Q₃) ⟹ SU(4)⊃SU(3)×U(1)、Euclid 読出し ⟹ Spin(4)≅SU(2)×SU(2)（Pati–Salam 型接続）
- **既知数学による検算経路**: M=6 の射影零集合＝複素二次超曲面 Q⁴⊂CP⁵、実形の選択が SO(3,3)・SO(4,2)・Spin(6)≅SU(4) を与える
- **無名性≠任意性**: 許容 sector は5拘束の交差 S_allowed = S_closure ∩ S_recurrence ∩ S_simplex ∩ S_harmonic ∩ S_self-consistent
- **地平線型境界候補**: t 従属読出し sector では x²+y²+z² ≥ R²+Q² が必要、等号面で t=0
- **det=1 の由来候補**: 全体位相冗長性の除去として導出できる見込み（離散中心の厳密化が残課題)

## 残る未解決問題（論文内 問題1〜6）

1. 三方向自発生成の解析的選択原理と一般化
2. 読出し解像度・stabilizer の自己無撞着選択則（3⊕2 vs 内部三重項 vs SU(4) 型）
3. R(q) からの一般 Riemann 動力学
4. simplex cochain からの局所ゲージ接続
5. 倍音パリティ統計と spin-statistics の一般対応・質量同定
6. chirality / hypercharge / anomaly cancellation
- 最終課題: N 自体の選択則（閉じれば外部自由入力ゼロ）

## 変更履歴

| 日付 | 内容 |
|---|---|
| 2026-08-20 | v1.0 公開。Zenodo 6ファイル（日英 md/tex/pdf）、DOI 取得、Zenn 記事公開 |
