# 第七思考実験 再実行：第六と補遺のプログラムを読み込み、補遺の初期化で初期状態を作る

**実行日:** 2026-09-27

## 何をしたか

第六思考実験のプログラムと、第六思考実験 補遺の初期化プログラムを、書き換えずにコピーして読み込んだ。

- 初期状態の 41 成分は、補遺の `init_case(case_id)` が作る。`init_case` は `case_relations` と `init_state_from_relations` を呼ぶ。
- 力学（`transition`、`macrostep`）、行の生成（`generate_chunk`）、保存形式、図の描き方は、第六のプログラムがそのまま行う。
- ラッパーが変えたのは、実行条件、保存先、そして第六の生成器が初期状態を受け取る関数 `init_state` の 1 箇所である。`init_state` は、補遺が作った初期状態を返すものに替えた。
- ラッパーの中に、力学や初期値の計算はない。

`01_独立自己相互作用数値実験_v1/` と `02_論文図_自己相互作用_v1/` は変更していない。

## コピーしたプログラム

`paper6_programs/` に置いた。コピー元とコピー先の SHA-256 は `paper6_programs/COPY_MANIFEST.json` に記録した。8 本とも一致している。

| 役割 | ファイル | コピー元 |
|---|---|---|
| 状態生成器 | `run_strict_charged8_5cases.py` | 第五 `13_…/00_CODE_AND_RULES/` |
| ケース別の図 | `postprocess_strict_charged8_5cases.py` | 第五 `13_…/00_CODE_AND_RULES/` |
| 解析基準との比較、ケース別の図、分析 | `postprocess_strict_5case.py` | 第五 `13_…_staging_duplicate/00_CODE/` |
| 独立解析基準の solver | `solve_charged_binary_analytic_reference_v1.py` | 第五 `07_荷電二体系_解析近似基準解_20260925/` |
| 独立解析基準の図 | `plot_charged_binary_analytic_reference_v1.py` | 第五 `07_荷電二体系_解析近似基準解_20260925/` |
| 条件比較の図 | `generate_paper6_orbit_condition_comparison_v1.py` | 第六 `01_軌道条件比較図/` |
| 対数半径の軌道図 | `plot_first_10_orbits_log_radial_scale.py` | 第五 `11_…/01_CODE/` |
| 初期化 G=q0=c=1 | `paper6_init_G_q0_c1.py` | 第六 `04_補講_Gq0c1_初期化同値検証_20260927/` |

## 実行条件

a の電荷は +n·q0、b の電荷は −n·q0。aa は a を両脚に、bb は b を両脚に入れる。規格化は G=q0=c=1、M/q0=20/3。

| ケース | n | (s_A, s_B) | λ_A | λ_B | C | D |
|---|---:|---|---:|---:|---:|---:|
| self_aa_n1 | 1 | (+, +) | 0.3 | 0.3 | 0.91 | 0.0 |
| self_bb_n1 | 1 | (−, −) | −0.3 | −0.3 | 0.91 | 0.0 |
| self_aa_n3 | 3 | (+, +) | 0.9 | 0.9 | 0.19 | 0.0 |
| self_bb_n3 | 3 | (−, −) | −0.9 | −0.9 | 0.19 | 0.0 |

n=4 は C=−0.44 で定義域外なので、実行していない。補遺の `init_case` が返した値（q0、M、m_A、m_B、λ、N、C、D）は `CASE_TABLE.json` にある。第六の `run_case` が `init_state` を呼んだ回数も同じファイルに記録した。4 ケースとも 1 回である。

## ラッパーと、変えた変数

| プログラム | 読み込む第六のプログラム | 変えた変数 |
|---|---|---|
| `step00_copy_programs.py` | — | コピーだけ |
| `wrapper_common.py` | 補遺 | 補遺の `CASES`（条件の表）を今回の 4 条件に替える |
| `wrap10_run_generator.py` | 補遺、状態生成器 | 生成器の `init_state`（補遺の `init_case` が作った初期状態を返す）、`CASES`、`HERE`、`OUTROOT`。ケースごとに生成器の `run_case` を呼ぶ |
| `wrap20_analytic_reference.py` | solver | 引数 `--outdir`、`--lambda-A`、`--lambda-B` |
| `wrap30_postprocess_figures.py` | 後処理 2 本 | `BASE`、`CASES`、`ANA`／`ANABASE`、`IDS`／`CASE_IDS`。実行は `main()` |
| `wrap40_orbit_condition_comparison.py` | 条件比較の図 | `CASES`、`PANEL_TITLES`、`STRICT_FILES`、引数 `--outdir`、`--strict-svg-dir` |
| `wrap50_log_radial_orbits.py` | 対数半径の図 | `RAW`、`SVG_OUT`、`PNG_OUT`。HDF5 の先頭 40,001 行を CSV に書き出して渡す |
| `wrap60_analytic_reference_figures.py` | 解析基準の図 | 引数 `--input`、`--outdir` |
| `compare_with_previous.py` | — | 保存済み出力との比較 |
| `inspect_Q_difference.py` | — | Q と t_readout の違いの詳細 |
| `step90_write_sha256sums.py` | 後処理（`sha256` 関数） | SHA-256 の一覧 |

## 実行結果

| ケース | macro steps | 最終 P | Q が非有限になった step |
|---|---:|---:|---:|
| self_aa_n1 | 1,456,531 | 19.999971790606857 | なし |
| self_bb_n1 | 1,456,531 | 19.999971790606857 | なし |
| self_aa_n3 | 15,266,914 | 19.999998005124105 | 6,316,005 |
| self_bb_n3 | 15,266,914 | 19.999998005124105 | 6,316,005 |

全 macrostep の行を、第六と同じ形式（HDF5、15 列、lzf 圧縮）で `cases/<ケース>/raw_macro_partNNN.h5` に保存した。

## 保存済み出力との比較

数値は `comparison_results/comparison_with_previous_ja.md` と `comparison_results/Q_difference_detail.json` にある。

### 第六思考実験の保存データとの比較（全行）

| 今回 | 第六 | 行数 | 違いがなかった列 | 違いがあった列 |
|---|---|---:|---|---|
| self_aa_n1、self_bb_n1 | n1_repulsive | 1,456,532 | 13 列 | Q、t_readout |
| self_aa_n3、self_bb_n3 | n3_repulsive | 15,266,915 | 13 列 | Q、t_readout |

- 違いがなかった 13 列は step、P、E、U_re、U_im、H_re、H_im、N、C、D、x_readout、y_readout、q_index。
- 最終状態の 41 成分、最初の macrostep の CSV、macro steps、最終 P、t_cross_r20、phi_cross_r20 は同じ。
- Q の違いは、その値の最小刻みで数えて n=1 が −8〜+31 個、n=3 が −80〜+124 個。
- t_readout の違いは、最小刻みで数えて n=1 が −2〜+4 個、n=3 が −2〜+2 個。
- 最初に Q が違う行は、n=1 が step 10,931、n=3 が step 13,616。

Q の更新は `Q × exp(Δt/τ0)` で、P の列は全行で同じなので、Δt も同じである。掛け算と割り算の結果は規格で決まっている。Q の違いは、指数関数 `exp` が返す値が実行環境によって最後の桁で違うために生じている。t_readout は対数関数 `log` を使う。

### 第七思考実験の旧プログラムの出力との比較

| 今回 | 旧 | 旧の保存行数 | 違いがあった列 | 最終状態 41 成分 |
|---|---|---:|---|---|
| self_aa_n1、self_bb_n1 | self_n1 | 1,458 | Q、t_readout | 同じ |
| self_aa_n3、self_bb_n3 | self_n3_fresh | 101 | Q、t_readout | 比べられない（旧は 100,001 step まで） |

旧プログラムだけにある列は P_GW、P_EM。

### 今回の aa と bb

n=1、n=3 とも、全行の全列と最終状態が同じ。

## 図

| 種類 | 場所 | 枚数 |
|---|---|---|
| ケース別 | `cases/<ケース>/figures/` | figure01〜06（後処理 2 本の出力、名前は 8 種類） |
| 対数半径の軌道 | `cases/<ケース>/figures/figure07_first_10_orbits_log_radial_scale.svg` | 1 |
| 条件比較 | `orbit_condition_comparison/figures/` | 4 |
| 独立解析基準（放射量など） | `analytic_reference/<ケース>/figures/` | 8 |

- ケース別の図は、生成器の出力と独立解析基準を重ねて描いている。
- 後処理 2 本は同じ `figures/` に書く。名前が同じ 4 枚は `postprocess_strict_charged8_5cases.py` の出力である。
- 条件比較の図の題名とファイル名にある "Five"、"five_case" は第六のプログラムの固定文字列である。今回は 4 条件なので、3×2 の枠のうち 5 番目が空の枠として描かれる。
- 条件比較の図 figure02 の PNG は作られていない。第六のプログラムが使う cairosvg がこのマシンにないためで、プログラムは警告を出して先へ進む。SVG は作られている。

## 実行環境

| 項目 | 値 |
|---|---|
| Python | 3.9.6 |
| numpy | 2.0.2 |
| numba | 0.60.0 |
| matplotlib | 3.9.4 |
| h5py | 3.14.0（`~/paper6_pydeps_py39`） |
| pandas | 2.3.3（`~/paper6_pydeps_py39`） |
| 機種 | Apple シリコンの Mac |

`postprocess_strict_5case.py` の 54 行目と 58 行目で、`matmul` の RuntimeWarning が 24 件出た。記録は `logs/wrap30_postprocess_figures.log` にある。書き出された値はすべて有限である。

## 再現

```sh
./run_all.sh
```

## 第六の保存フォルダについて分かったこと

| 対象 | 事実 |
|---|---|
| `n3_repulsive/part_manifest.json` | 4 part × 4,000,000 行と書かれているが、フォルダにあるのは 16 本 × 1,000,000 行。合うのは `part_manifest (1).json` |
| `reference_comparison.json` を作るプログラム | 第六の論文 11.3 節のフォルダにない。`13_…_staging_duplicate/00_CODE/postprocess_strict_5case.py` にある |
| 後処理プログラムの保存先 | `/mnt/data/…` の絶対パスが書かれている |
| `(1)` の付いたファイル | 同名ファイルの 2 回目のアップロードが残っている |
