# リリースノート：三次元正方配列の三角形ループ積における時空・内部空間読出し接続仮説（Part III v1）

- **公開日**：2026-11-28
- **著者**：木原 範昭（株式会社WFシステム、ORCID: `0009-0004-6753-4020`）
- **Zenodo URL**：https://zenodo.org/records/24625963
- **Concept DOI（全版共通）**：`10.5281/zenodo.24625962`（https://doi.org/10.5281/zenodo.24625962）
- **Version DOI（v1）**：`10.5281/zenodo.24625963`（https://doi.org/10.5281/zenodo.24625963）
- **Zenn 記事**：`articles/triangle-loop-spacetime-readout-connection-part3.md`
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
  - 四元数の関係読み出しと結晶化インフレーション（2026-10-21）：`10.5281/zenodo.23523237`
  - 複素数を正しく理解していたか（2026-10-23）：`10.5281/zenodo.23585436`
  - 観測者の分解能はゲージである（2026-11-10）：`10.5281/zenodo.24097611`
  - 有界なヒルベルト空間と八元数による時空・素粒子の創発（2026-11-16）：`10.5281/zenodo.24267347`
  - 投影軸の選択と保存則の共変性（2026-11-16）：`10.5281/zenodo.24273723`
  - コインとコップの思考実験（2026-11-26）：`10.5281/zenodo.24562867`

## 概要

Part I・Part II 執筆後に残った 3 つの不満（①状態空間から $(x,y,z,t,R,Q_1,Q_2,Q_3)$ への写像が未接続、②分解能による「見える／見えない」の定義が欠落、③$n$ の過去光円錐扇型座標としての読み直しと `1–7` 八元数の `1+3+4` 役割整理）を動機とし、先行思考実験の読出し幾何学を $\mathrm{Tri}(T)$ の実発展状態へ直接接続した仮説論文（Part III）である。

決定論的検証（`tri_loop_spacetime_log_resolution_bridge_v0.py`、乱数なし）により以下 5 点を確認した：
1. **`(1+3)+(1+3)` ヌル条件と 1,2,4D 投影軸共変性**：真空 $H$ の全 7 種四元数分割で $|a|^2 = |b|^2 = 0.500000$（最大差 `3.331e-16`）、1D/2D/4D 投影軸選択で比率 $d/8$ が全選択間偏差 `5.55e-17 ~ 2.22e-16` で厳密共変。
2. **外側 $O_8^{\mathrm{ext}}$ の `1+3+4` 純粋 $SU(2)$**：$G_2(14) \supset SO(4)(6) \supset SU(2)_{\mathrm{pure4}}(3)$、3 空間 $\{e_1,e_2,e_3\}$ への作用誤差 `8.327e-17`、4 空間 $\{e_4,e_5,e_6,e_7\}$ の $2\pi$ 回転スピン $1/2$ 誤差 `5.778e-15`（$O_8^{\mathrm{ext}}\otimes O_8^{\mathrm{int}}$ 64 自由度で重力＋$SU(2)_L$＋$SU(3)_C\times U(1)_Q$ が重複なく収まる）。
3. **対数変換 $O=e^{\ln A+n\theta}$ による局所曲率半径 $R$ の縮小と地平面**：往復誤差 `1.665e-16`・ボルン／ローレンツ 2 断面誤差 `2.220e-16`。質量振幅 $\epsilon = 0 \to 0.270$ で $R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$ が `1.0000 → 0.8328` へ単調縮小し、$\epsilon = 0.275$ で `0.0000`（事象の地平面切断）へ不連続遷移。
4. **自分中心の過去光円錐と潮汐分裂**：真空ではどの観測者からも全 7 方向完全等方（$0.000000 \pm 10^{-16}$）、質量存在下の遠方観測者（$A_{\mathrm{obs}}=1$）には質量方向（`-0.006094`）vs 接線 6 方向（`-0.034709`）の潮汐分裂 $\Delta D_{\mathrm{tidal}} = +0.028615$ が出現。
5. **分解能ゲージによる「見える／見えない」の定義**：辺長閾値 $\delta_{\mathrm{res}}$ を 18 桁掃引しても見える集合（64 個）と視野角 $\langle\rho\rangle_w = 0.7767719701$ は 10 桁不変。また $\mathbb{Z}_2^3$ 倍音ローパス窓 $K_{\mathrm{cut}} = m$ が Hurwitz 4 段階粒子（$m=0,1,2,3$）の完全可視化境界と厳密一致。

## 公開ファイル（全 9 点、Zenodo 登録済み）

作業ディレクトリ：`次元の生成構造/三次元正方配列の積と自己無撞着_20261127/`

- `三次元正方配列の三角形ループ積における時空・内部空間読出し接続仮説_Part3.pdf`（日本語版 8 ページ、2 図付き）
- `triangle_loop_spacetime_readout_connection_hypothesis_part3_en.pdf`（英語版 7 ページ、2 図付き）
- `三次元正方配列の三角形ループ積における時空・内部空間読出し接続仮説_Part3.tex` / `.md`
- `triangle_loop_spacetime_readout_connection_hypothesis_part3_en.tex` / `.md`
- `build_tex_part3.py`
- `figures_and_programs_part3.zip`（検証スクリプト・作図スクリプト・ログ・日英図一式）
- `SHA256SUMS_part3`
