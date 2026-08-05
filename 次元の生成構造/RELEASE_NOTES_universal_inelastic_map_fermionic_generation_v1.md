# Release Notes: フェルミオン的構造の生成は誘導・自己触媒・対相関で起こる——波形のみを入力とする万能非弾性写像の仮定と帰結（第9論文続編・生成の力学編）

Release date: 2026-08-05

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21808091
- Version DOI (v1): https://doi.org/10.5281/zenodo.21808092
- Zenodo record: https://zenodo.org/record/21808092

## Zenn

- 記事: https://zenn.dev/noriaki_kihara/articles/universal-inelastic-map-fermion-gen （2026-08-05、push後に公開）

## note（一般向け）

- 日本語: https://note.com/kiharanoriaki/n/n6718389c48ea （2026-08-05 公開）
  （元原稿: `次元の生成構造/反物質は仕込まなくても出てきた_note完成稿.md`）
- 英語: https://note.com/kiharanoriaki/n/n280577cf8bf8 （2026-08-05 公開）
  （元原稿: `次元の生成構造/antimatter_came_out_unbidden_note_en.md`）

## Facebook

- 日本語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid02JXo3dtdApqvo9V17YjxytyudBCo8JRZsz9oHMUPvuFCfdeDbzZ1YBEHFqQN9hjn8l （2026-08-05、日本語note URLはコメント欄）
- 英語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid02m99mRzYCYdbEsxXRVVNreQotfUFLn69HHPsbfzUyem9AXLFM9cfuxuPtgQ4tBGQtl （2026-08-05、英語note URLはコメント欄）

## X

- 未実施（公開後に追記）

## 公開ファイル（Zenodo, 6点）

- `universal_inelastic_map_fermionic_generation_ja.md` / `_en.md`
- `universal_inelastic_map_fermionic_generation_ja.tex` / `_en.tex`
- `universal_inelastic_map_fermionic_generation_ja.pdf` / `_en.pdf`

（正本: `次元の生成構造/万能非弾性写像_managed_v1/フェルミオン的構造の生成_万能非弾性写像_完成論文_v1.md`、図3枚 fig1–fig3 同フォルダ）

## 主張（三層分離）

- **性格**: 相互作用の正しさを証明する論文ではない。(a) 仮定＝非弾性相互作用のクラス（二波型・点ごと・三次・共通位相不変・自己散乱なし・強度=R=sin²θ流用、IF文なしの万能演算）、(b) 条件付き導出＝閉塞保存が g₂=−g₁ を一意化（クラス内一意性）、(c) 実測——を分離して報告。仮定の価値は整合性と生産性で測る。
- **設計原理（第一柱）**: 写像は波形 (a,b) のみを受け取り、粒子種を判別する分岐を持たない。「ボゾン/フェルミオン/種/相棒」は全て読出し側の分類であり動力学側のラベルではない。
- **閉形式厳密解（§5.1/5.2）**: 頂点は δa=−2R·Im(b̄a)·b, δb=+2R·Im(b̄a)·a と厳密に書き換わる——駆動スカラーは相対位相の虚部ひとつ。s=Im(b̄a) が流れの下で保存されるため各格子点で閉形式解（角度 φ=2R·Im(b̄a) の点ごと二チャネル回転）。閉塞・パワーは点ごとに恒等保存（実測ドリフト 4.1e-14/3000衝突）。生成＝完全保存された回転角の空間非一様性によるスペクトル再配分。弾性部と同一の回転生成方向（統合読みは未検証課題）。
- **パリティ定理**: 純偶ポンプから奇数内容は恒等ゼロ（<1e-25）——生成は必ず誘導過程。偶奇類=Z/2Zの創発的保存量（群論は事後分類）。
- **点火則**: rate=C·f²、C=11.45、f₀の4桁で一定性0.4%。点火時間∝1/f₀²。collapse表示で軌跡全体がf₀²スケール。自己触媒則は先行例未発見（新規性候補・反証歓迎）。
- **運命**: 暴走なし。統計的平衡 f*=0.4690 ≈ マスク位相空間割合0.494（等分配読み、差5%の帰属は未解決課題）。
- **census三判定**: 相棒和則 k*=2k_p−k_s のみに生成（予言外485bin機械ゼロ、比1.4e-27）／毛の帳簿 q*=+4 でコヒーレンス0.832（帳簿は頂点1回作用から機械導出）／厳密ゼロから同時成長。**反粒子は入力ではなく出力**——「相棒=生成頂点が全保存帳簿を閉じるために同時生成する共役な関係状態」という非循環的定義。
- **展望**: 多体接続6項（二因子ゲート・媒介定理のM=2特別場合・C·f²移植p=2.001・50/50くじ＋ラチェット・レジスタ位置＋調和閉鎖運動・census完全移植）。未解決課題4項（f²則の摂動導出・f*の状態数導出・生成子代数の抽出・仮定縮約）。

## 反証と修正の記録（本文§11、6件）

1. v1頂点 g₂=+g₁ の閉塞破れ→一意化定理の発見に転化
2. census v1 熱化後測定→熱化前40衝突窓へ
3. census v2 広帯域分散＋η一様平均→狭帯域単一巻き＋毛分解相関（ゼロは選択則だった）
4. 毛帳簿の手計算誤り（m_B=−1仮定）→頂点1回作用からの機械導出（q*=+4厳密一致）
5. RK2積分バイアス（ドリフト2.9e-3、C=10.4）→RK4較正（2.8e-9、C=11.45）→閉形式解発見で積分器不要が最終確定（4.1e-14）
6. 「閉塞が平均を選ぶ」→対対称性のみが定理、平均選択は動力学の実測事実に降格

## 方法の透明性

- 原本＝第9論文実行環境 ab_invariant_theta_toy_v1（SHA-256 90c7b272…、bit一致確認済み、read-only import）
- 実験正本: `次元の生成構造/万能非弾性写像_managed_v1/`（v1実装＋T1–T5、v3確定版、RK4再走行、閉形式検証 run_ignition_fate_exact_v3.py、census v3、図生成＋census RK4再測定）
- 判定基準は全実験で実行前固定（R1–R3 / E1–E3 / P1–P3）。本文の全数値は結果JSONと一致
- 多体接続の参照: `次元の生成構造/万能相互作用多体接続_v1/`（コミット系列）

## 引用

- 自己引用: 第9論文（Concept 21766706）、開始様式判別論文（21798854）、逆二乗則論文（21441081）、三部作A 相互作用の二文法分解（21763995）
- 外部: Menyuk 1987（結合NLS）、Agrawal NFO（四波混合）、Milonni 1976（半古典論に自発放出なし）、Greene–Kofman 1999/2000＋KLS 1997（preheating）、FPU LA-1940、Baudin et al. PRL 125, 244101 (2020)（RJ凝縮）、Manley–Rowe 1956

## 変更履歴

- 2026-08-05 v1 公開（Version DOI 21808092）。日本語正本＋英訳＋tex/pdf 6点。査読反映済み：閉形式厳密解の発見（§5.1/5.2）、s不変性の途中式明示、Z₂事後分類の位置づけ、共役相棒の非循環的定義、未解決課題4項（§13.1）。
