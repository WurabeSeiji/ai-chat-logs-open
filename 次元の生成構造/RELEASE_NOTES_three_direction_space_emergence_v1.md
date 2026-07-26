# Release Notes: 三方向空間の創発（次元の生成構造 第7論文）

Release date: 2026-07-26

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21578401
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.21578402
- Zenodo record: https://zenodo.org/records/21578402

（公開直後は Version DOI の doi.org 解決に DataCite 伝播の遅延がある。Zenodo レコード API は公開直後から open で解決する。Concept DOI レコードは 302 で最新版へ転送＝正常。）

## Zenn

- https://zenn.dev/noriaki_kihara/articles/three-direction-space-emergence
  （元原稿: `articles/three-direction-space-emergence.md`）

## 公開ファイル（Zenodo, 29点）

**論文（JP/EN）**
- `paper7_ja.md` / `paper7_en.md`
- `paper7_ja.tex` / `paper7_en.tex`
- `paper7_ja.pdf` / `paper7_en.pdf`（各13ページ, LaTeX Error 0）

**図（4点）**
- `figure1_compare_N5_N40_N300.png`（分裂量 f）
- `figure2_compare_N5_N40_N300.png`（五成分占有 stack）
- `figure3_compare_N5_N40_N300.png`（五成分占有 log）
- `transverse_growth_compare_N5_N40_N300.png`（横摂動応答）

**再現プログラム（9点, Python 3）**
- `run_n_scaling_lowrank_v1.py`（原本エンジン, SHA `ba0fc19b…`, 不変更）
- `run_plane_flow_exact_v1.py` / `run_plane_flow_approx_v1.py`（固定親基底）
- `run_n300_dimension_saturation_v2.py`（支配平面 Gram 縮約）
- `run_paper7_5color_timeseries.py`（5色占有時系列）
- `run_paper7_transverse.py` / `run_paper7_transverse_cached.py`（横安定性）
- `run_paper7_exact_vs_approx_N40.py`（§6.2 検証）
- `make_paper7_figures.py`（図生成）

**実験結果（8点）＋報告書＋README**
- `paper7_long_timeseries_N000{05,40,300}.csv`（5色占有）
- `transverse_stability_timeseries_N000{05,40,300}.csv`（横摂動）
- `N_comparison_table.csv` / `transverse_stability_summary.csv`
- `paper7_longtime_and_transverse_stability_report.md`（数値報告）
- `README.md`（実行条件・ファイル一覧）

## 主張

二方向初期状態から独立な第三空間方向が生成され、三方向空間からなる動的準安定閉包が形成される。数値的には初期親平面と瞬時支配平面の結合階数が 2→4 へ増加し、転移後の全記録時刻（761点, N=5,40,300）で rank=4 を維持。準安定域 q₃/q₄ は N=5=0.751/0.631, N=40=0.338/0.311, N=300=0.200/0.195。調査時間（0〜55000）では三方向を超える有限占有方向は自然発展で観測されず、三方向閉包の外側には正の横成長率（λ⊥,max/σ₁ = 7.3e-4 / 9.65e-5 / 1.21e-5）を持つ横方向応答が残る。

## 方法の透明性

- 論文6の固定親基底3色（P1/その他回転/核）を維持し、確定指示に基づき「その他回転」を新方向3・新方向4・残余へ分解して5色化。P1と核は不変更。新方向は時間依存 S₄(τ)=orth[B₀|B_dom(τ)] の B₀直交補2方向を固定 other 空間へ射影して構成（縮退平面は2次元として追跡、固有値順の色割当なし）。
- 共通横軸 絶対step 0〜55000（crossing 不動）、目盛り5000刻み。全図（個別・比較）共通。
- 支配平面 B_dom は Gram 縮約 G=WᵀW（≤600次元）。N=40 で密行列法と倍精度一致（最大偏差 1.78e-15）。
- 保存誤差 ≤ 2.0e-15、五成分射影閉鎖 ≤ 2.2e-16。
- §6.4 に横摂動 Benettin の warm-start 内部状態同期の実装課題を明記（λ⊥ の定量値のみに影響、rank Q:2→4・占有・閉包には非影響）。

## 引用

- 自己引用: 生成子ランク線形上界・三方向飽和（21465898）、平面分解読出し（21468959）、自発的分裂の開始と帰結分類（21486233）、波の数は分解能（21486544）、自発的分裂の停止と新しい直交回転平面の創発（第6論文, 21543070）
- 外部: Horn-Johnson（歪対称正準形）、Taghavi-Chabert（純スピノル/ツイスター）、Walker（回転からの時間創発）、Smith / Furey / Todorov / Manogue-Dray-Wilson（D4・八元数・E8 比較対象）

## 生成手順

- 数式は GitHub 記法（`$$` / `$`）。tex/PDF は /tmp/tex_compile で生成（日本語=lualatex、英語=pdflatex、各2回、LaTeX Error 0）。DOI は tex/PDF 生成前に md へ埋め込み（PDF 内にも正しい Version/Concept DOI が載る）。
- 大容量アーティファクト（v2 全体 zip、全固有ベクトル npz）は git 追跡外（Drive 上に保持）。
