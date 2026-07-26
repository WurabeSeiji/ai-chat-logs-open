# Stage A1b: 第7論文 N=5 最初の主成長エピソード構造観察

## 目的

Stage A0で完全再現されたN=5 CSV 3件だけを後処理し、最初の大域的なf増幅、既存direction 3/4の状態占有、既存q3/q4、既存crossing=1167、増幅後の飽和・振動への移行を、人間が同じ絶対step軸で読める資料にする。

単一の成長開始・終了・方向成立時刻は採用しない。仮説H1/H2/H0も判定しない。

## 読み取る入力

- Stage A0 `fcurve_N00005_delta1e-15_seed0.csv`
- Stage A0 `q_svd_N00005.csv`
- Stage A0 `paper7_long_timeseries.csv`

`verify_inputs.py` はStage A0報告書の完全再現判定、CSV比較記録、固定SHA-256、実ファイルSHA-256、必須列、絶対step対応を照合する。

原本力学コードはimportも実行もしない。N=40、N=300、Stage A1の候補表も読み取らない。

## 固定観察窓

基本範囲はabsolute step 0〜3000である。追加表示・記述統計窓は次のとおりである。

- 0〜500
- 500〜1000
- 800〜1400
- 1000〜1800
- 1400〜2500

これらは既存crossing=1167を含む最初の遷移を観察する表示窓であり、イベント定義ではない。3000より後を主成長候補として探索しない。

## f記述量

- `f`
- `log10_f`: `f>0` のときだけ計算
- `running_max_f`
- `running_max_log10_f`
- `f_diff_1 = f(step)-f(step-1)`
- `log(f)` の中心回帰slopeとR²: 11, 21, 41, 81, 161 step窓

`f<=0` は変更せず、対数不能値はNaNとする。回帰窓は完全な連続窓だけを使う。成長区間は選ばない。

## 初回到達座標

観察窓内の正の最小fを含む10進decadeの下端から0.1までの各 `10^k` と、指定された追加水準を重複なしで列挙する。各水準の最初の `f>=level` を座標として保存するが、採用閾値とはみなさない。

隣接する列挙水準間について、

```text
step_difference = upper_first_step - lower_first_step
threshold_log_amplitude_difference = ln(upper_level/lower_level)
mean_exponential_rate = threshold_log_amplitude_difference / step_difference
```

を保存する。step差が0または水準未到達なら平均率はNaNとする。

## 実保存値と表示専用補間

状態占有は25 step間隔の実保存レコードを解析表に使う。初回到達stepとの対応は、実保存stepの `before_or_at` と `after_or_at` を列挙する。

図の線を見やすくするためだけに、状態占有の1 step線形補間ファイルを別に作る。列名に `display_interp` を付け、イベント候補・初回到達表・記述統計には使わない。

qは実保存レコードだけを使い、未保存stepを補間しない。図では実保存点をマーカーで示す。

## 区別する既存観測

- A: 既存 `rank_q` が4になった保存レコード
- B: q3/q4が有限または非ゼロの数値として保存されたこと
- C: 既存direction 3/4の状態占有が増加したこと
- D: fのdecade初回到達が進んだこと
- E: 既存定義の `f>0.05` crossing=1167

これらを同一イベントとみなさない。

## 実行順序

1. `verify_inputs.py`
2. `analyze_first_transition.py`
3. `make_transition_figures.py`
4. `make_transition_report.py`

各工程は直前工程の成功を検査し、自身の出力が既にあれば上書きせず停止する。

## 出力

指定5表をCSVとMarkdownで保存する。全時間域候補の直積は生成しない。図1〜12をPNG/SVGで生成し、指定された全拡大範囲を一枚にまとめた補助図もPNG/SVGで保存する。

最終報告書は `reports/paper7_N5_transition_anatomy_report.md` である。完成後はStage A1bで停止し、Stage A2へ進まない。
