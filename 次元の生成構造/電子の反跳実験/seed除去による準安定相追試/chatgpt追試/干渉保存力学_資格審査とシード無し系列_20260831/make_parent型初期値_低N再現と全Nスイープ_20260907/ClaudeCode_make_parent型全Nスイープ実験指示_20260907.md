# Claude Code 実験指示：canonical make_parent型親による N=3..40 全N 500 step スイープ

作成日: 2026-09-07

## 0. 目的

`make_parent` が自然生成する Z4 位相＋非等振幅固有ベクトル構造を、以前の高対称理論床系列と同一の stage123 力学・同一の読出しで N=3..40 全域に展開し、インフレーションの発生・成長率・step0/500複素平面を比較する。

これは新しい力学を試す実験ではない。**初期親の枝を canonical make_parent 型に固定した比較実験**である。

## 1. 絶対条件

1. canonical stage123 力学を変更しない。
2. `one_step`、`H_of`、`adjacency`、`plane`、`metrics` の数式を変更しない。
3. `Z0` に重心閉塞 `sum z=0` を追加しない。
4. 振幅を等振幅化しない。
5. 位相を手で四相へ丸めない。四相は `make_parent` が生成した結果としてのみ扱う。
6. 既存の `N3_N40_stage123_sweep_20260905` を上書き・破壊しない。
7. 新しい専用フォルダへ全出力を保存する。
8. float64 / complex128 を維持する。

## 2. 入力親の生成

本フォルダの `make_static_parents_N3_N40_v1.py` を参照し、元プロジェクトの検証済み `run_n_scaling_lowrank_v1.py` から `LowRankSystem`, `make_parent`, `zero_closure_kernel_seed` を import して N=3..40 を生成する。

固定条件:

```text
seed = 0
rng = default_rng(40260722 + 1000*N)
iters = 1200
tol = 1e-12
delta = 1e-15
Z0 = (v + delta*g)/||v+delta*g||
```

生成物は新規 `parents_makeparent/` に保存する。

### 入力ゲート

- N=3..7 は本パッケージ `parents_actual_N3_N7/` と `v,g,Z0` を比較し、一致しない場合は中断して原因を報告。
- 可能なら既存 `N3_N40_stage123_sweep_20260905/parents/` の N=3..40 と全配列比較する。
- 環境差でbit一致しない場合は勝手に丸めず、最大絶対差、残差、位相差を報告して停止判断を求める。

## 3. 親構造監査

全Nについて以下を保存する。

- `abs(sum z^2)`
- `abs(sum z)`
- 位相を大域位相 modulo pi/2 で見た Z4 残差
- 振幅種数（許容誤差を明記）
- 振幅 CV
- `W=A*sin^2(delta theta)`
- `Wr=sigma r` の Rayleigh 固有値と相対残差
- 四相パワー P0,P1,P2,P3
- `P_even=P0+P2`, `P_odd=P1+P3`

CSV: `parent_structure_N3_N40.csv`

この監査は分類のみ。親を変更しない。

## 4. 500 step 本走行

canonical `run_N3_N40_stage123_v1.py` と同じ条件で N=3..40。

```text
STEPS = 500
denominators = N-2, N-1, N, N+1, N+2, 124  (positive only)
```

各状態を保存:

```text
states/hm_N{N}_den_{den}_states_500.npz
```

集計:

- `timeseries_64bit_with124_N3_N40_makeparent.csv`
- `summary_64bit_with124_N3_N40_makeparent.csv`

summaryには最低限、

```text
N, series, denominator,
onset_gt_0.05,
initial, step1, step10, step100, final, max,
closure0, closure500
```

を含める。

## 5. 必須図

### A. インフレーション図

N=3..40 の 8x5 グリッド。

- 横軸 step 0..500
- 縦軸 Hperp/H、log
- 各Nについて N-2,N-1,N,N+1,N+2,124 を重ねる
- 既存図と比較可能な範囲・レイアウトにする

ファイル名:

`fig_Hperp_makeparent_N3_N40_500.png`

### B. step0 複素平面

全N 3..40、den=N の Z[0]。

`fig_complex_plane_step0_makeparent_N3_N40.png`

### C. step500 複素平面

全N 3..40、den=N の Z[500]。

`fig_complex_plane_step500_makeparent_N3_N40.png`

### D. step0/500 対比

少なくとも N=3..10 を左右2列で直接比較。

`fig_complex_plane_step0_step500_makeparent_N3_N10.png`

## 6. 高対称理論床系列との直接比較

既存の `最も対称性の高い初期値_20260906` の全N結果を読取り専用で使用し、同じ N, den=N について比較する。

必須比較量:

- onset_gt_0.05
- step1, step10, step100, step500 Hperp/H
- 初期 `abs(sum z)`
- 初期 `abs(sum z^2)`
- 振幅CV
- 位相四相残差
- 最大初期線形成長率が既に保存されていればそれも

CSV:

`compare_makeparent_vs_theoretical_floor_N3_N40.csv`

図:

`fig_compare_makeparent_vs_theoretical_floor_onset_growth.png`

特に N=3 は、make_parent型で成長し理論床で復元する差を最重要比較点として扱う。

## 7. 判定上の注意

- `make_parent` の Z4 位相は入力仮定ではなく出力構造として記録する。
- `sum z=0` の大小を「合否」にしない。
- `sum z^2=0` と `Wr=sigma r` を別々に監査する。
- N>=5 の枝依存性が既知なので、canonical seed系列の結果を一般N唯一解と呼ばない。
- 以前の人工四相等振幅NG系列を canonical make_parent 型と混同しない。

## 8. 実験完了時に作るレポート

`make_parent型全Nスイープ結果_20260907.md`

内容:

1. 実行環境・SHA256
2. 親生成ゲート結果
3. N=3..40構造監査
4. 500 step全結果
5. 複素平面の観察
6. 高対称理論床との比較
7. N=3の枝差
8. 証明済み / 数値確認 / 仮説 の三分類
9. 次の実験候補

## 9. 既存資産の保護

元フォルダは絶対に変更しない。新規結果は、この指示書のあるフォルダ内に新しい `full_N3_N40_sweep/` を作って保存すること。

実行前に、使う canonical ソースと入力ファイルの SHA256 を `SHA256_INPUTS.txt` に保存する。実行後に全生成物を `SHA256SUMS.txt` にまとめる。
