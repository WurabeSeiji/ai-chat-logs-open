# Release Notes: 自己無撞着な関係波閉鎖系におけるインフレーション的急拡大の機構——正規化監査、rank 生成、二乗閉包保存、simplex 対称化および公理系の再構成（次元の生成構造 第8論文続編）

Release date: 2026-08-26

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.22112008
- Version DOI (v1): https://doi.org/10.5281/zenodo.22112009
- Zenodo record: https://zenodo.org/record/22112009

## Zenn

- 記事: https://zenn.dev/noriaki_kihara/articles/self-consistent-inflation-mechanism （2026-08-26 push。Zenn の新規記事レート制限のため公開反映は要確認）

## note（一般向け）

- 日本語: https://note.com/kiharanoriaki/n/n07c3e4c97e3a （2026-08-26 公開）
  （元原稿: `次元の生成構造/インフレーションは倒れかけの波に書き込まれていた_note完成稿.md`、図: `次元の生成構造/note_figs_self_consistent_inflation/`）
- 英語: 原稿 `次元の生成構造/inflation_written_in_the_tilting_wave_note_en.md`（図: 同フォルダ `en/`）、公開後に URL 追記
- 前2記事の続編（「タネがなくても」 n48a02cd70f47 ／「雑踏からは」 n7b655e3da27a）

## Facebook / X

- 未投稿

## 公開ファイル（Zenodo, 21点）

- `nbody_self_consistent_inflation_mechanism_ja.md` / `_en.md`
- `nbody_self_consistent_inflation_mechanism_ja.tex` / `_en.tex`
- `nbody_self_consistent_inflation_mechanism_ja.pdf`（40p, lualatex）/ `_en.pdf`（36p, pdflatex）
- `nbody_self_consistent_inflation_figures_v1.zip`（本文図 PNG 21枚＝日本語図18＋英語版インフォグラフィック3、英語版作図スクリプト `make_N5_infographics_en.py`）
- 再現パッケージ zip 14点（本文「再現性に関する注記」の一覧と同一）：
  `K_sigma_normalization_artifact_test_N4_N5_20260826.zip`, `N5_gamma_continuum_test_bundle_20260825.zip`,
  `N3_N4_…`/`N5_…`/`N6_N7_…`/`N8_N9_…`/`N10_N11_…`/`N12_N13_…`/`N14_N15_…`/`N16_complex_simplex_complete_analysis_20260826.zip`,
  `N3_N16_partial_zero_closure_analysis_20260826.zip`, `N3_N16_nontrivial_zero_closure_analysis_20260826.zip`,
  `N14_N16_complete_nontrivial_zero_closure_search_20260826.zip`, `N5_dynamics_followup_theorems_and_stability_20260826.zip`

## 主張

第8論文が確立した seedなし急拡大（潜伏→指数増大→rank-4→三方向準安定）を、コード監査・N=3〜16 再実験・N=5 線形安定性解析で再構成した。

1. **厳密保存定理**：実反対称生成子の Cayley 更新 $C_t=(I-\gamma K_t)^{-1}(I+\gamma K_t)$ は実直交行列であり、$Z^\dagger Z$ と $Z^TZ$（二乗ゼロ閉包）を厳密保存する。急拡大は $H_\parallel\to H_\perp$ の内部移送（$H_\parallel+H_\perp$ の実装誤差 $4.4\times10^{-16}$）で、$H_\perp\le H_{\rm total}$ が無限増大を運動学的に禁止する。
2. **開始機構＝相対平衡の線形不安定性（三重整合）**：`make_parent` 固定点残差 $3.87\times10^{-7}$〜$2.38\times10^{-13}$ で onset は 72/134/176/238 step、成長率 0.17251/step 不変。実測残差に対する回帰 $t_{\rm onset}=11.6162[-\ln\varepsilon]-99.5631$（$R^2=0.99999201$）。rotating-frame Jacobian の支配不安定乗数は実2重固有値 $\mu_1=1.090086569$ で、$2\ln\mu_1=0.172514$＝実測成長率、$1/\ln\mu_1=11.593$＝回帰傾き 11.616。親 rank-2＋支配不安定2次元＝rank-4。第2対 $\mu_2=1.052603$、$\sum\ln|\mu_i|=-0.0855$（局所体積収縮）。
3. **正規化はアーティファクトでない**：$C(K/\sigma,\gamma)=C(K,\gamma/\sigma)$。有限 step の成長率差 6.8% は刻み収束則 $g(n)=1.1596346-4.10498/n$（残差 $8\times10^{-5}$）と $n_{\rm eff}=144/3.49$ から $1.0600$ と予測され実測 $1.05874$ と 0.13% で一致。144 は 2.5° 刻みの分母。
4. **ヌル錐定理**：全頂点 star 二乗閉包 ⟺ 重心表示で全頂点が複素ヌル錐上（$x_i^2=0$）。等分配と合わせ、準安定状態＝equimodular null complex simplex（Gram rank $N-1$、N=3〜16）。スペクトルエントロピーは step 5000 で $S/\ln M=1.000000$（非単調）。
5. **N=4/N=5 の導出**：N=4 の120°は頂点閉包＋等モジュラスから定理。N=5 の 3+3+2+2 符号対構造から非自明2辺閉包 $3\cdot3+2\cdot2=13$、exact cover $3!2!=12$。8-seed 掃引で 3+3+2+2 と等モジュラスは全 run 再現、二距離族間相対位相は平坦方向。N=14 の6辺 quasi-closure（$10^{-6}$ で停止）は厳密閉包でない。
6. **公理系再構成**：自己無撞着固定点⇒複素回転対⇒$\sum z^2=0$⇒$S^1$ コンパクト位相軌道までは導出。$U^n=I$ は $\Delta\theta/2\pi\in\mathbb Q$ を選ぶ有理ロック機構が別途必要（導出したとは主張しない）。

## 反証と修正の記録

- 要求 tolerance と実測固定点残差の取り違えによる回帰統計の見かけの矛盾（4点整数 onset で R²≈0.9945 に見える）→ 実測残差表を掲載し R²=0.99999201 を再現
- 「支配乗数は複素共役対」という誤った機構説明 → `floquet_spectrum.csv` で虚部 0 を確認し「実2重固有値」へ訂正。対称性起源（Jacobian と複素構造の可換性）は未証明として §27 に登録
- 有限 step 差 6.8% を「追加収束試験に委ねる」とした保留 → §12 の既存系列で 0.13% 予測に決着
- 「大域位相同変性⇒複素線形」の断定を撤回（同変性は $DF$ の共役関係のみを与える）

## 方法の透明性

- 正本 md：`次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/自己無撞着な関係波閉鎖系におけるインフレーション的急拡大の機構_完成論文_v6.md`（公開版 `nbody_self_consistent_inflation_mechanism_ja.md` は図パスのみ変更）
- 図 19 枚と英語版インフォグラフィック3枚（`make_N5_infographics_en.py`、解析 CSV から matplotlib で再構成）は同フォルダ
- 再現パッケージ 14 zip は同フォルダ（`N5_gamma_continuum_test_bundle_20260825.zip` のみ `…/第5論文原本_自発的分裂予備実験_v1/A2a_N5_ab_probe_20260825/`）。各 zip に SHA256SUMS
- 分業：ChatGPT が生成・Drive 上で最終修正、Claude（Fable 5）が査読・検算・英訳・英語版図・組版・公開

## 引用

- 自己引用：基本公理系 v9（K1）、K2〜K7（Concept DOI 21465898/21468959/21486233/21486544）、第8論文（Version 21614403 / Concept 21614402）。Zenodo related_identifiers：第8論文 21614403、開始様式判別論文 21798855、公理系 21315735
- 外部：Aste 2019（arXiv:1905.12894）、Aoki–Hirasawa–Ito–Nishimura–Tsuchiya PTEP 2019 093B03、Hirasawa ほか 2024（arXiv:2407.03491）、H. Kihara–Nitta–Sasaki–Yoo–Zaballa PRD 80, 066004 (2009)

## 生成手順

- tex/PDF は `/tmp/tex_compile/build.py` で生成（日本語＝ltjsarticle+lualatex 2回、英語＝article+pdflatex 2回、LaTeX Error 0・Missing character 0）
- 前処理：YAML フロントマター（title/author/date）、図 alt text と太字キャプション行の統合、長文 `\boxed{}` の `\fbox{\parbox}` 化、数式外の Unicode 記号を数式化、英語版は CJK パスをローマ字化、zip 名は `\path{}`

## 変更履歴

- 2026-08-26 v1 公開（Version DOI 22112009）。
