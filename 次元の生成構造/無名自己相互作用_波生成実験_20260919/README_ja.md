# 無名・自己相互作用込み複素波生成則 数値実験 2026-09-19

## 目的
初期値として複素波 `z1 = a + i b` を1個だけ与え、粒子名・隣接行列・外部時刻・Δτを使わず、自己相互作用を含む同一演算を全波に適用したときの波生成を愚直に追跡する。

## 固定した生成ルール
現在存在する N 波 `z_i` に対し、順序付き全組 `(i,j)`（`i=j` を含む）を一項ずつ計算する。

`T_ij = (conj(z_i) * z_j) * z_j`

全 `N^2` 項を加算し、

`L_N = sum_i sum_j T_ij`

その符号反転を新しい1波として追加する。

`z_(N+1) = -L_N`

既存波は捨てない。コードの更新部では代数的簡約を使わず、二重ループで全項を実際に加算する。

## 入力
物理的な初期入力は `a`, `b` のみ。`max-generations` と `stop-abs` は計算停止用の数値ガードであり、力学パラメータではない。

## テストした主な初期値
- `1`, `i`
- `1/sqrt(2)`, `i/sqrt(2)`
- `sqrt(2)`, `i sqrt(2)`
- `1.001`, `0.999`
- 位相45度で絶対値 `1`, `1.001`, `0.999`

## 実測上の要点
- `|z1|=1` の実軸・虚軸は、第2波が反対向きとなり、その後の生成値は0になる。
- `|z1|=1/sqrt(2)` は220波まで非零生成が続き、生成振幅は非常に小さく減衰した。
- `|z1|=sqrt(2)` は急速に増大し、数値ガードに到達した。
- `|z1|=1.001` では第3波以降に小振幅の交互符号系列が現れ、generation 3..52 の `|z_new|` は直線近似で R^2=0.994648842、傾き=4.21237596022e-05。ただし長く走ると generation 約200 付近で急増に移るため、永久的な線形則とは未確定。
- `|z1|=0.999` の generation 3..52 は R^2=0.996527460、傾き=-2.54233657998e-05 で緩やかに減衰する。
- 位相0と45度で同じ絶対値を与えたケースは、今回の実測では生成振幅系列が数値精度内で一致した。位相方向自体は保持され、生成点は同一直線上に並んだ。これは観測事実であり、一般証明はこの実験ファイルでは行っていない。

## 主要ファイル
- `anonymous_interaction_wave_generator_20260919.py`: 単一初期値 a,b と組込みテストスイートを実行する本体
- `results/suite_summary.csv`: テストケース一覧と停止状況
- `results/suite_nearcritical_zoom.png`: |z1|≈1 の小振幅枝の拡大図
- `results/suite_amplitude_log.png`: 全テスト比較（対数振幅）
- `results/suite_complex_plane.png`: 複素平面上の生成波
- `linear_fit_generation_3_52.json`: 近臨界ケースの generation 3..52 の直線近似診断
- `anonymous_interaction_wave_generation_full_20260919.zip`: 全CSV・個別図を含む再現パッケージ

## 実行例
`python anonymous_interaction_wave_generator_20260919.py 1.001 0 --outdir run_1p001`

組込みテスト:
`python anonymous_interaction_wave_generator_20260919.py --suite --max-generations 220 --outdir results`
