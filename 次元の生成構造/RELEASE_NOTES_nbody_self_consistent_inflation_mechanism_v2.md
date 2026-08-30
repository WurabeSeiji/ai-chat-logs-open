# Release Notes: 自己無撞着な関係波閉鎖系におけるインフレーション的急拡大の機構——第 2 版：計算条件の訂正と N=3〜16 の再検証（次元の生成構造 第8論文続編、訂正版）

Release date: 2026-08-30

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.22112008（維持。常に最新版へ転送）
- Version DOI (v2): https://doi.org/10.5281/zenodo.22176949
- Zenodo record: https://zenodo.org/records/22176949
- 訂正対象 v1: https://doi.org/10.5281/zenodo.22112009（2026-08-27 公開。レコードとファイルはそのまま保存）

## note（一般向け、v1 記事のバージョンアップ）

- 日本語: https://note.com/kiharanoriaki/n/n0b5cc8dbcee5 （2026-08-30 公開。v1 記事 n07c3e4c97e3a の第 2 版）
- 英語: https://note.com/kiharanoriaki/n/na53b0313ef06 （2026-08-30 公開。v1 英語記事 nf2549985b38a の第 2 版）

## Facebook

- 日本語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid02KLYtqPm2scLcpkEjLi69XwoFttvhp6gV42k6LidQiicYe2USMHyGvpGooFgRq1Dnl （2026-08-30、DOI・note URL はコメント欄）
- 英語投稿: https://www.facebook.com/kihara.noriaki/posts/pfbid04hSezsbUEKJ6PCQm14RcCCqyuQQEjFKo6kZCHRWivxz5499usMAPyVrtt8LdqxPnl （2026-08-30、同上）

## Zenn

- 記事（全面書き換え、v2 の内容）: https://zenn.dev/noriaki_kihara/articles/self-consistent-inflation-mechanism （v1 時点の push は Zenn の投稿数上限で未公開のままだった。v2 push 後に 200 を要確認）

## 公開ファイル（Zenodo v2、36 点）

- `nbody_self_consistent_inflation_mechanism_ja.md` / `_en.md`（v1 と同名で差替え）
- `nbody_self_consistent_inflation_mechanism_ja.tex` / `_en.tex`（lualatex、ja 22p / en 21p、LaTeX Error 0・Missing character 0）
- `nbody_self_consistent_inflation_mechanism_ja.pdf` / `_en.pdf`
- `nbody_self_consistent_inflation_figures_v2.zip`（本文図 11 枚）
- 新規の再現パッケージ zip 14 点（git 管理分。全 step 状態 `states_treatment.npz` は除外、`run_all.sh` で再生成可）：
  `v2補完実験_4生成法_N3toN16_統一プロトコル_20260830.zip`（54 走行、走行前予測、図 2〜8）、
  `公理見直し_ゼロ閉塞定理と固有時計_20260829.zip`、`Nall_linear1000000_steps5_mpmath50_…_20260828.zip`（弾道則、多倍長 50 桁）、
  `飽和ステップ数とNの関係_固定点ヤコビアン解析_20260829.zip`（二層構造）、`手作り自己無撞着親と対称性_倍音と関係数の検討_20260829.zip`、
  `複素シンプレックス基礎_N別全展開_20260830.zip`、`複素シンプレックス基礎_N別全展開_非等モジュラー版_20260830.zip`、`複素シンプレックス_重心閉塞_非等モジュラー族_20260830.zip`、
  `N5_linear124_…directHperp_treatment_only_20260828.zip`、`N16_linear124_…directHperp_treatment_only_20260828.zip`、
  `論文v1_全プログラム修正版_20260828_light.zip`（fixed/fixed_baseline/fixed_equimodular の走行木は除外、results・md・スクリプトのみ）、`論文v1_全再現テスト_20260828_light.zip`
- v1 から継承：図バンドル `nbody_self_consistent_inflation_figures_v1.zip` と再現パッケージ zip 14 点（v1 再現用）

## 主張（3 つ）

1. **ゼロ閉塞は定理**：実反対称生成子の非零固有モード（自己無撞着）から $\sum z^2=0$ が従う（5 行証明）。正規化・振幅分布・$N$・親の生成法に依らない。110 親の実装で数値検証。
2. **複素シンプレックスは制約を与えない**：任意の複素 2 乗距離が Autonne–Takagi 分解で $\mathbb C^{N-1}$ に厳密埋め込み（乱数 1400 状態で検証）。形は符号を忘れた像 $v/(\mathbb Z_2)^M$、符号枝は $K(Sv)=SK(v)S$ で力学的に共役。選択原理ではない。
3. **計算条件変更後もインフレーション的発展（線形不安定性＋非線形飽和）は再現**：成長率は v1 の 10〜30 分の 1、飽和後は局在（PR/M 0.05〜0.31）、発生の有無は $N$ 単独でなく初期状態の構造で決まる。4 生成法 × $N=3$〜16 の 54 走行で走行前予測 53/54 一致、λ 比 0.997〜1.008。高対称系列に奇偶非対称（偶数 6〜16 飽和／奇数 5〜15 床）、ただし乱数均衡親の奇数 $N$ に不安定あり。

## 3 つの変更

1. 隠れた振幅正規化の除去（`make_parent` の `v/‖v‖`、位相のみ相互作用 → 振幅込み $K_{ij}=\mathrm{Im}(\bar z_iz_j)$、外部 seed と初期正規化の除去）
2. Cayley 有理写像 → 凍結生成子の指数写像 $\exp(\Delta K)$（連続流 $dv/d\tau=K(v)v$ の一次積分器）
3. 初期化を厳密な自己無撞着解に（3-1 等モジュラー／3-2 非等モジュラー）

## v1 の主張の維持／修正／撤回（本文 §9 の表）

- 維持：厳密保存、内部移送、ヌル錐定理（用語は複素ヌル錐へ）、$N=4$ の 120°（定理としては）、自己無撞着⇒閉塞⇒$S^1$ 軌道
- 修正：onset–残差則（→弾道則 $(r\tau)^2$＋親依存の $\lambda_f$）、Floquet 三重整合（→1 刻み線形化行列の予測と実測の一致）、記事図 1 の「31 桁」（底は精度と親残差が決める数値的な底で物理の底ではない）、§20 表（rank 維持・等分配撤回、$\sigma=N-1$ は $\mu=-(N-1)r^2$ として維持）
- 撤回：等分配＝等モジュラー・ヌル単体、$N=5$ の 13/12 閉包、$3+3+2+2$、8 seed モジュライ、IIB／動的コンパクト化／preheating との類似
- 訂正：v1 §4「make_parent は監査済み」

## 査読の記録

- ChatGPT による 3 回の査読（2026-08-30）。主要修正：状態→形は 1:1 でなく $v/(\mathbb Z_2)^M$／「非線形→線形」でなく「Cayley 有理近似→凍結生成子の指数写像」／連続流と積分器の分離／$L$ 非本質は連続極限に限定／奇偶非対称は高対称系列の性質で原因は未分離／確認したのは線形不安定性と飽和で開始機構ではない／53/54 が主結果・t50 規則は事後／「選ばない」→「一意に選択できない」／閉塞の起源は 3 条件／Takagi の 2 段追加／§10「正準積分でない」→「$H_{\rm int}$ を厳密保存しない」／弾道則は ≃／ハミルトニアン符号規約／序論に射程限定の一文
- 木原指摘：底の桁数は計算精度の下限であって物理の底ではない（桁が減ったと読ませない）

## 引用

- 自己引用：[K1] 公理系 v9、[K2] 21465898、[K8] 21614402、[K9] 21798854、[v1] 22112009
- 外部：Aste 2019、Takagi 1925（Horn–Johnson Cor. 4.4.4）、Schoenberg 1935。v1 の [E2]〜[E4] は削除
- Zenodo related_identifiers：21614403、21798855、21315735、21465898、21798854

## 生成手順

- 正本 md：`chatgpt追試/自己無撞着な関係波閉鎖系におけるインフレーション的急拡大の機構_第2版_完成論文_v2_draft1.md`（公開版 `nbody_self_consistent_inflation_mechanism_ja.md` は同内容、図は `論文v2_figures/`）
- 英訳：`nbody_self_consistent_inflation_mechanism_en.md`（Claude 訳）
- tex/PDF：pandoc（`-s -f markdown-implicit_figures`、YAML フロントマターで title/author/date、図は `figs/`）→ lualatex 2 回。ja＝ltjsarticle、en＝article＋luatexja（コードスパン内の日本語パッケージ名のため）
- Zenodo：`newversion` で Concept 22112008 を維持して draft 22176949 を作成 → prereserve DOI を md に埋め込み → tex/pdf 生成 → bucket へ PUT → metadata（title/description/version=v2/keywords/related）更新 → publish

## 変更履歴

| 日付 | 内容 |
|---|---|
| 2026-08-30 | v2 公開（Version DOI 22176949）。3 変更・3 主張・維持／修正／撤回表・N=3〜16 結果行列・図 11 枚 |
| 2026-08-27 | v1 公開（Version DOI 22112009） |
