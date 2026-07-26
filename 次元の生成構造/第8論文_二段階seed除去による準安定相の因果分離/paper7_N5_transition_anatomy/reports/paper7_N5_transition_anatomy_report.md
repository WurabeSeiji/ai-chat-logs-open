# 第7論文 N=5 最初の主成長エピソード構造観察報告（Stage A1b）

## 1. 範囲と総合状態

- 入力: Stage A0で完全再現されたN=5 CSV 3件だけ
- 基本観察範囲: absolute step 0〜3000
- 固定拡大範囲: `[('0-500', 0, 500), ('500-1000', 500, 1000), ('800-1400', 800, 1400), ('1000-1800', 1000, 1800), ('1400-2500', 1400, 2500)]`
- 原本力学コード、新しい軌道、第8論文本実験、高精度計算、N=40、N=300: 未実行
- 全時間域候補の直積: 未生成
- 単一の成長開始・終了・方向成立時刻: 未採用
- H1/H2/H0: 未判定
- **総合状態: `TRANSITION_ANATOMY_COMPLETE`**

## 2. 入力SHA-256

| 入力 | 行数 | 絶対パス | SHA-256 | Stage A0照合 |
|---|---:|---|---|---|
| `fcurve` | 21168 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/reproduced/metastable_series_result_v1/fcurve_N00005_delta1e-15_seed0.csv` | `9220c5f3c1f570c8a52ea24a3cdd95568354cea0943d9bee7d8ed20316d3a9d0` | 一致 |
| `q_svd` | 995 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/reproduced/exact_lowN_eigenspectrum_v2/raw/N00005_dimension_saturation_v2/q_svd_N00005.csv` | `7c16a364c6cc9145293c2625dfe4ebb1f9962655d212679188215e8fad5e5155` | 一致 |
| `paper7_long` | 2201 | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_reproduction/reproduced/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/paper7_long_timeseries.csv` | `efeaf9dab753c057ad0c6109b9e4a8919f8d8db1249da186658bfed9fda784e3` | 一致 |

## 3. 時間軸と保存間隔

- 共通時間軸: `absolute_step`
- f: 全範囲 `[0, 21167]`、観察窓内1 step間隔
- q: 全範囲 `[0, 51167]`、観察窓内保存間隔 `[5]` step
- 状態占有: 全範囲 `[0, 55000]`、観察窓内保存間隔 `[25]` step
- 全q保存行で `step-relative_time=[1167]`
- q未保存stepは補間していない。
- 状態占有の解析表・記述統計・初回到達対応は実保存値だけを使用した。
- 状態占有の線形補間は表示専用ファイルへ分離し、図の線にだけ使用した。

## 4. fの各水準初回到達step

| level | source | status | first step | f at step | step-crossing |
|---:|---|---|---:|---:|---:|
| `1e-30` | `decade` | `found` | 0 | 1.06601506e-30 | -1167 |
| `1e-29` | `decade` | `found` | 2 | 2.13645696e-29 | -1165 |
| `1e-28` | `decade` | `found` | 5 | 1.37658972e-28 | -1162 |
| `1e-27` | `decade` | `found` | 14 | 1.07949587e-27 | -1153 |
| `1e-26` | `decade` | `found` | 36 | 1.00706240e-26 | -1131 |
| `1e-25` | `decade` | `found` | 68 | 1.00569367e-25 | -1099 |
| `1e-24` | `decade` | `found` | 107 | 1.01418050e-24 | -1060 |
| `1e-23` | `decade` | `found` | 151 | 1.04311246e-23 | -1016 |
| `1e-22` | `decade` | `found` | 196 | 1.01521643e-22 | -971 |
| `1e-21` | `decade` | `found` | 242 | 1.00116465e-21 | -925 |
| `1e-20` | `decade` | `found` | 289 | 1.02472296e-20 | -878 |
| `1e-19` | `decade` | `found` | 336 | 1.04462827e-19 | -831 |
| `1e-18` | `decade` | `found` | 382 | 1.01234928e-18 | -785 |
| `1e-17` | `decade` | `found` | 429 | 1.03031637e-17 | -738 |
| `1e-16` | `decade` | `found` | 476 | 1.04846946e-16 | -691 |
| `1e-15` | `decade` | `found` | 522 | 1.01551456e-15 | -645 |
| `1e-14` | `decade` | `found` | 569 | 1.03335305e-14 | -598 |
| `1e-13` | `decade` | `found` | 615 | 1.00085718e-13 | -552 |
| `1e-12` | `decade+explicit` | `found` | 662 | 1.01843287e-12 | -505 |
| `1e-11` | `decade` | `found` | 709 | 1.03631679e-11 | -458 |
| `1e-10` | `decade+explicit` | `found` | 755 | 1.00372596e-10 | -412 |
| `1e-9` | `decade` | `found` | 802 | 1.02135144e-09 | -365 |
| `1e-8` | `decade+explicit` | `found` | 849 | 1.03928637e-08 | -318 |
| `1e-7` | `decade` | `found` | 895 | 1.00660166e-07 | -272 |
| `1e-6` | `decade+explicit` | `found` | 942 | 1.02427305e-06 | -225 |
| `1e-5` | `decade` | `found` | 989 | 1.04221190e-05 | -178 |
| `1e-4` | `decade+explicit` | `found` | 1035 | 0.00010089934 | -132 |
| `1e-3` | `decade+explicit` | `found` | 1082 | 0.0010221299 | -85 |
| `1e-2` | `decade+explicit` | `found` | 1130 | 0.01042778 | -37 |
| `0.050000000000000003` | `explicit` | `found` | 1167 | 0.051404687 | 0 |
| `1e-1` | `decade+explicit` | `found` | 1188 | 0.10092563 | 21 |

これらは増幅過程の座標であり、採用イベント閾値ではない。正の最小fを含むdecadeの下端から0.1までを列挙したため、最初の複数水準が同じ保存stepに到達する場合も削除していない。

## 5. decade間の増幅率の安定性

| decade pair | lower step | upper step | Δstep | ln(level ratio)/Δstep |
|---|---:|---:|---:|---:|
| `1.00000000e-30 → 1.00000000e-29` | 0 | 2 | 2 | 1.1512925 |
| `1.00000000e-29 → 1.00000000e-28` | 2 | 5 | 3 | 0.76752836 |
| `1.00000000e-28 → 1.00000000e-27` | 5 | 14 | 9 | 0.25584279 |
| `1.00000000e-27 → 1.00000000e-26` | 14 | 36 | 22 | 0.10466296 |
| `1.00000000e-26 → 1.00000000e-25` | 36 | 68 | 32 | 0.071955784 |
| `1.00000000e-25 → 1.00000000e-24` | 68 | 107 | 39 | 0.059040643 |
| `1.00000000e-24 → 1.00000000e-23` | 107 | 151 | 44 | 0.052331479 |
| `1.00000000e-23 → 1.00000000e-22` | 151 | 196 | 45 | 0.051168558 |
| `1.00000000e-22 → 1.00000000e-21` | 196 | 242 | 46 | 0.050056198 |
| `1.00000000e-21 → 1.00000000e-20` | 242 | 289 | 47 | 0.048991172 |
| `1.00000000e-20 → 1.00000000e-19` | 289 | 336 | 47 | 0.048991172 |
| `1.00000000e-19 → 1.00000000e-18` | 336 | 382 | 46 | 0.050056198 |
| `1.00000000e-18 → 1.00000000e-17` | 382 | 429 | 47 | 0.048991172 |
| `1.00000000e-17 → 1.00000000e-16` | 429 | 476 | 47 | 0.048991172 |
| `1.00000000e-16 → 1.00000000e-15` | 476 | 522 | 46 | 0.050056198 |
| `1.00000000e-15 → 1.00000000e-14` | 522 | 569 | 47 | 0.048991172 |
| `1.00000000e-14 → 1.00000000e-13` | 569 | 615 | 46 | 0.050056198 |
| `1.00000000e-13 → 1.00000000e-12` | 615 | 662 | 47 | 0.048991172 |
| `1.00000000e-12 → 1.00000000e-11` | 662 | 709 | 47 | 0.048991172 |
| `1.00000000e-11 → 1.00000000e-10` | 709 | 755 | 46 | 0.050056198 |
| `1.00000000e-10 → 1.00000000e-09` | 755 | 802 | 47 | 0.048991172 |
| `1.00000000e-09 → 1.00000000e-08` | 802 | 849 | 47 | 0.048991172 |
| `1.00000000e-08 → 1.00000000e-07` | 849 | 895 | 46 | 0.050056198 |
| `1.00000000e-07 → 1.00000000e-06` | 895 | 942 | 47 | 0.048991172 |
| `1.00000000e-06 → 1.00000000e-05` | 942 | 989 | 47 | 0.048991172 |
| `1.00000000e-05 → 0.0001` | 989 | 1035 | 46 | 0.050056198 |
| `0.0001 → 0.001` | 1035 | 1082 | 47 | 0.048991172 |
| `0.001 → 0.01` | 1082 | 1130 | 48 | 0.047970523 |

- 有効なdecade-to-decade平均率数: `28`
- 最小 / Q25 / 中央値 / Q75 / 最大: `0.047970523 / 0.048991172 / 0.050056198 / 0.051459288 / 1.1512925`
- 平均 / 標準偏差: `0.1250029 / 0.24024737`

平均率は水準間で完全一定ではない。分布と各行を提示するだけで、安定区間や主指数を自動選択していない。

## 6. direction 3/4占有が増加する固定観察step帯

次表は人間が指定した固定窓内の最初と最後の実保存値、および差である。増加帯の境界を新たに推定していない。

| 固定窓 | d3 first→last | d3差 | d4 first→last | d4差 |
|---|---|---:|---|---:|
| `0-500` | 1.06166888e-31 → 5.37713245e-25 | 5.37713139e-25 | 7.12594079e-32 → 5.37275173e-25 | 5.37275102e-25 |
| `500-1000` | 5.37713245e-25 → 5.61913219e-11 | 5.61913219e-11 | 5.37275173e-25 → 3.81841855e-11 | 3.81841855e-11 |
| `800-1400` | 1.94342688e-19 → 0.023071065 | 0.023071065 | 1.04208686e-19 → 0.19474398 | 0.19474398 |
| `1000-1800` | 5.61913219e-11 → 0.21300514 | 0.21300514 | 3.81841855e-11 → 0.20984638 | 0.20984638 |
| `1400-2500` | 0.023071065 → 0.31864123 | 0.29557016 | 0.19474398 → 0.045157354 | -0.14958663 |

## 7. q3/q4が増加する固定観察step帯

| 固定窓 | q3/q1 first→last | 差 | q4/q1 first→last | 差 |
|---|---|---:|---|---:|
| `0-500` | 1.58955059e-08 → 1.51219108e-08 | -7.73595166e-10 | 0 → 1.02173735e-08 | 1.02173735e-08 |
| `500-1000` | 1.51219108e-08 → 0.0034834195 | 0.0034834044 | 1.02173735e-08 → 2.39971206e-06 | 2.38949469e-06 |
| `800-1400` | 2.50198306e-05 → 0.60813577 | 0.60811075 | 1.05367121e-08 → 0.24814264 | 0.24814262 |
| `1000-1800` | 0.0034834195 → 0.73324339 | 0.72975997 | 2.39971206e-06 → 0.50660507 | 0.50660267 |
| `1400-2500` | 0.60813577 → 0.75235514 | 0.14421937 | 0.24814264 → 0.54921781 | 0.30107517 |

## 8. A〜Eの重要な区別

| 区別 | データ上の記述 | 同一視しない理由 |
|---|---|---|
| A. rank_q=4 | 最初の実保存rank_q=4はstep `265`。その行のq3/q1=`1.29047841e-08`、q4/q1=`1.05367121e-08` | 既存相対閾値への応答であり、状態占有増加やf増幅の成立を意味しない |
| B. q3/q4の有限・非ゼロ値 | q3は初期保存行から有限値。q4の最初の正の実保存値はstep `15` | 有限値・非ゼロ値と大きさの増加は別である |
| C. direction 3/4占有増加 | 25 step間隔の実保存占有を固定窓ごとに上表へ記載 | qのrankや有限値とは異なる状態占有である |
| D. fの指数増幅 | decade初回到達のstep間隔と平均指数率として記載 | 単一rank行や単一占有行とは異なる区間的挙動である |
| E. f>0.05 | 既存crossing=`1167` | 第7論文既存閾値であり、A〜Dを定義しない |

### rank_qの早期数値床応答

- 最初のrank_q=4保存行: step `265`、f=`3.12731256e-21`。
- その時点を挟む状態占有の実保存step: `250` と `275`。
- beforeのdirection 3/4: `3.91565632e-28` / `4.16967107e-28`。
- afterのdirection 3/4: `7.95559530e-28` / `8.48683563e-28`。
- crossing前のrank_q=4保存行数: `97`。

rank_q=4がfの大域増幅およびdirection 3/4の有限占有より大幅に早く現れるため、既存rank閾値が早期の数値床に反応していることと整合する。ただし、この後処理だけで数値誤差の原因や物理的無効性までは証明しない。

## 9. crossing=1167前後の実保存値

fは毎step、qは5 step、占有は25 step保存なので、同一stepへ補間せず各系列の実レコードで挟む。

| 系列 | before step | before値 | after/crossing step | after/crossing値 |
|---|---:|---|---:|---|
| f | 1166 | f=0.049510489 | 1167 | f=0.051404687 |
| occupation | 1150 | d3=0.00011566435, d4=8.23579784e-05, kernel=0.025606995 | 1175 | d3=0.0008289504, d4=0.00061128974, kernel=0.066833249 |
| q | 1165 | q3/q1=0.17444856, q4/q1=0.0018572455, rank=4 | 1170 | q3/q1=0.19028221, q4/q1=0.0012716289, rank=4 |

crossingの1 step前からcrossingまでにfが既存0.05水準を超える。一方、qと状態占有は異なる保存間隔で既に連続的に変化しており、同一stepで新たに全量が同時成立したとは読まない。

## 10. crossing後に遅れて変化する量

| 固定step | f | d3実保存 | d4実保存 | kernel実保存 | q3/q1実保存 | q4/q1実保存 | rank_q |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1200 | 0.13332104 | 0.0027224709 | 0.0035100285 | 0.12706355 | 0.27944827 | 0.015543983 | 4 |
| 1400 | 0.58549718 | 0.023071065 | 0.19474398 | 0.36693866 | 0.60813577 | 0.24814264 | 4 |
| 1800 | 0.84034596 | 0.21300514 | 0.20984638 | 0.41289221 | 0.73324339 | 0.50660507 | 4 |
| 2500 | 0.80118766 | 0.31864123 | 0.045157354 | 0.43372349 | 0.75235514 | 0.54921781 | 4 |

この固定step表と図9・10は、crossing後もdirection 3/4占有、kernel、q比が同じ速度では変化しないことを示す。遅れの開始・終了時刻は採用していない。

## 11. 増幅から飽和・振動へ移る固定観察帯

- 800〜1400のf: min `9.25337815e-10`、max `0.58549718`、正差分 `600`、負差分 `0`。
- 1400〜2500のf: min `0.50653988`、max `0.9599702`、正差分 `546`、負差分 `554`。
- 1400〜2500のrunning maximum: `0.58549718` → `0.9599702`。

固定窓1400〜2500ではfの正負両方向差分が多数あり、単調なdecade通過だけでは記述できない飽和・振動的挙動が見える。これは移行帯の観察であり、成長終了時刻の決定ではない。

## 12. データから直接言えること

1. fはstep 0〜3000の範囲で多数の10進水準を順次通過し、既存crossing=1167で初めて0.05を超える。
2. decade間のstep差と平均指数率は保存データから直接計算できるが、全decadeで完全一定ではない。
3. 既存rank_q=4はcrossingや大域的状態占有増加より非常に早い保存行で現れる。
4. direction 3/4占有、q3/q4比、kernel、fは、共通absolute step上で異なる推移を示す。
5. crossing前後の各系列は保存間隔が異なるため、同時性は実保存レコードの範囲でしか言えない。
6. 1400〜2500の固定窓ではfに増減が共存し、初期の単調な増幅だけとは異なる。

## 13. データだけでは言えないこと

1. 単一の指数成長開始時刻、終了時刻、方向成立時刻。
2. rank_q=4の早期応答が物理方向成立を意味するかどうか。
3. q未保存stepまたは占有未保存stepでの厳密な同時性。
4. direction 3とdirection 4のどちらを先行方向と解釈すべきか。
5. 増幅から飽和・振動へ移る境界の一意性。
6. H1/H2/H0のいずれが正しいか。

## 14. 人間が次に固定すべき最小定義

1. 「主増幅」をrunning maximumのdecade通過で定義するか、log(f)回帰で定義するか。
2. 増幅の終了をrunning maximumの停滞、局所slope、振動幅のどれで記述するか。
3. q3/q4について数値床を除外する最小振幅・持続・保存レコード規則。
4. direction 3/4状態占有について、増加を認定する振幅と持続の規則。
5. 異なる保存間隔の系列間で、同時・先行・遅延を記述するbracketing規則。
6. 既存crossing=1167を参照座標のまま使うか、他イベントの定義へ組み込むか。

これらを人間が固定するまで、単一イベント時刻を採用しない。

## 15. 出力表

| 表 | CSV | Markdown |
|---|---|---|
| `f_first_passage_levels` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/f_first_passage_levels.csv` (3.55 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/f_first_passage_levels.md` (4.40 KiB) |
| `f_decade_growth_rates` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/f_decade_growth_rates.csv` (4.58 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/f_decade_growth_rates.md` (5.47 KiB) |
| `first_passage_nearest_occupation_records` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/first_passage_nearest_occupation_records.csv` (9.08 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/first_passage_nearest_occupation_records.md` (10.62 KiB) |
| `first_passage_nearest_q_records` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/first_passage_nearest_q_records.csv` (9.21 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/first_passage_nearest_q_records.md` (10.87 KiB) |
| `transition_window_descriptive_statistics` | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/transition_window_descriptive_statistics.csv` (40.74 KiB) | `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/processed/transition_window_descriptive_statistics.md` (47.75 KiB) |

## 16. 出力図

- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure01_f_running_max_0_3000.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure01_f_running_max_0_3000.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure02_log10_f_and_slopes_0_3000.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure02_log10_f_and_slopes_0_3000.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure03_direction_1_to_4_occupation_0_3000.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure03_direction_1_to_4_occupation_0_3000.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure04_other_kernel_splitting_0_3000.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure04_other_kernel_splitting_0_3000.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure05_q1_to_q4_0_3000.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure05_q1_to_q4_0_3000.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure06_q_ratios_and_rank_q_0_3000.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure06_q_ratios_and_rank_q_0_3000.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure07_first_passage_vs_direction_3_4_actual_records.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure07_first_passage_vs_direction_3_4_actual_records.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure08_first_passage_vs_q3_q4_actual_records.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure08_first_passage_vs_q3_q4_actual_records.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure09_all_quantities_800_1400.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure09_all_quantities_800_1400.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure10_all_quantities_1000_1800.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure10_all_quantities_1000_1800.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure11_first_passage_level_step_intervals.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure11_first_passage_level_step_intervals.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure12_adjacent_level_mean_exponential_rates.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure12_adjacent_level_mean_exponential_rates.svg`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/supplement_all_requested_zoom_ranges.png`
- `/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/supplement_all_requested_zoom_ranges.svg`

図の単一イベント参照線は既存crossing=1167だけである。状態占有の線形補間は表示専用、qは実保存点だけを使用した。

## 17. 実行環境と停止

- Python: `3.9.6 (default, Jan  9 2026, 11:03:41) 
[Clang 17.0.0 (clang-1700.6.4.2)]`
- NumPy: `2.0.2`
- OS: `macOS-26.3.1-arm64-arm-64bit`
- 報告生成日時（UTC）: `2026-07-26T05:49:11.297157+00:00`
- 報告生成時間（秒）: `0.029355`
- **最終状態: `TRANSITION_ANATOMY_COMPLETE`**
- Stage A1bはここで停止する。Stage A2へ進まない。
