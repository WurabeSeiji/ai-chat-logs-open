# Release Notes: 幾何級数的急拡大は不安定な自己無撞着閉包に固有である——開始様式の因果判別（次元の生成構造 第8論文続編・線形相完結）

Release date: 2026-08-05

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21798854
- Version DOI (v1): https://doi.org/10.5281/zenodo.21798855
- Zenodo record: https://zenodo.org/record/21798855

## Zenn

- 記事: https://zenn.dev/noriaki_kihara/articles/onset-mode-unstable-equilibrium （2026-08-05、push後に公開）

## note（一般向け）

- 日本語: https://note.com/kiharanoriaki/n/n7b655e3da27a （2026-08-05 公開）
  （元原稿: `次元の生成構造/雑踏からは宇宙は始まらない_note完成稿.md`）
- 英語: https://note.com/kiharanoriaki/n/n1b83f7b50e0e （2026-08-05 公開）
  （元原稿: `次元の生成構造/no_universe_from_noise_note_en.md`）
- 前回記事の続編（「タネがなくても、インフレーションは起きた」 https://note.com/kiharanoriaki/n/n48a02cd70f47 ）

## Facebook

- 日本語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid0vj2owZQ7YFsdawuKdUhKmfBsVk7LKzXiQryyAqZDPzoNQULDmJWAeymKf8MnJsL1l （2026-08-05、日本語note URLはコメント欄）
- 英語投稿: 未公開（公開後にURLを追記）

## 公開ファイル（Zenodo, 6点）

- `nbody_onset_mode_causal_discrimination_ja.md` / `_en.md`
- `nbody_onset_mode_causal_discrimination_ja.tex` / `_en.tex`
- `nbody_onset_mode_causal_discrimination_ja.pdf` / `_en.pdf`

## 主張

第8論文§9.6が事前宣言した系統的対照（閉鎖条件のみを満たす一般状態からの出発）を実行した。開始様式は完全に二分する：自己無撞着閉包＝潜伏バースト型（潜伏288〜695 step→log fが18桁線形、N=5: 0.0494/step, R²=0.9999998）／非平衡零閉塞82本（白色セクター42＋一般状態40）＝全て即時型で増幅過程が不在。例外なし。終着点（三方向準安定・rank 4）は両様式共通で、開始様式を区別する化石は増幅の痕跡だけ。

機構検証3実験（ONS-1/2/3、判定基準・予言は実行前固定）：①一段残差が平衡4状態（≤2.5e-15）と非平衡82状態（≥1.8e-2）を13桁分離・中間なし ②接線写像 λ_max=rate_f/2 が N=5 で厳密一致（0.02468）、不活性だった白色起源親 N=5 は不安定固有値ゼロの安定平衡 ③摂動振幅の対数則（N=5 傾き比0.996）と ε≲3e-13 での飽和→内在床＝親固有モード残差2.14e-13 の力学的測定（第8論文未解決(i)の部分回答）。検査した全86状態で「潜伏指数増幅⟺不安定な相対平衡」が例外なく成立。三分類：非平衡=即時離脱／安定平衡=離脱なし／不安定平衡=潜伏後指数離脱。

先行研究の位置づけ（§10）：独自設計の対照が作業仮説（白色でも増幅するはず）を棄却し、事後の再調査でインフレーション初期条件問題（40年来の schism）と同構造であることが判明。標準理論は開始状態を仮定し（Guth/Linde/Brandenberger/Planck）、「非一様から着火」主張（East–Kleban–Linde–Senatore、Clough–Lim ら）は初期データがプラトー制限でコヒーレンスを密輸——本結果と矛盾しない。宇宙論との対応は構造的類比に限定（§8.6）。

## 反証と修正の記録（本文§7）

- 生成器v2の欠陥（k/N−k交差対相殺による閉塞・自己対セクター非閉塞）→v3で源から射影除去（定理の履行）
- 実験側の閉塞測定の誤り（誤った二次形式）→正しい場の閉塞で白色セクター42本は厳密零閉塞
- いずれも原本SHA一致コピー＋誤結果のbit再現後に修正。中心の二分はv2/v3を跨いで存続

## 方法の透明性

- 力学コード＝第7論文原本 read-only import（SHA-256: 75a10a5b…）、生成器v3（SHA-256: d3217579…）
- 実験正本: `第8論文_二段階seed除去による準安定相の因果分離/paper8_em9r_white_harmonics_inflation_v1/`（E-M9r v1/fix_v2/v3/profile）＋ `同/paper8_onset_mechanism_v1/`（ONS-1/2/3）
- 生成器正本: `make_parent_white_managed_v1/`（v2対照・v3・v4・単体テスト・分類器対照）
- コミット系列: 52c83f13→aeaf8e57→d9346e37→e5232885→94aa76f9→b1cded81→e6681744→2f443467（機構検証）
- 乱数シード全記録（生成器 N=5:2/N=40:1、η=93000+j、直接標本94000+j、混合95000+j、注入補助91000/92000+k/96000+series）

## 引用

- 自己引用: 基本公理系（Concept DOI 21315735）、第7論文（21543070）、第8論文（21614402）、零二乗和解説ノート（21495305）、第9論文（21766706）
- 外部: Guth 1981、Linde 1983、Brandenberger 2017 (arXiv:1601.01918)、Linde 2018 (1710.04278)、East–Kleban–Linde–Senatore 2016 (1511.05143)、Clough ほか 2017 (1608.04408)、Ijjas–Steinhardt–Loeb 2013/2014 (1304.2785/1402.6980)、Planck 2018 X (1807.06211)

## 変更履歴

- 2026-08-05 v1 公開（Version DOI 21798855）。日本語正本 `幾何級数的急拡大は不安定な自己無撞着閉包に固有である_完成論文_v1.md`。査読反映済み：線形無生成は仮説として明示（境界条件7）、表題に「不安定な」を挿入、三分類明文化、§10は模型内結果先行。
