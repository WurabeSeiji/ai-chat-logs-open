# 第七思考実験・aa/bb 独立自己相互作用 数値実験 v1

## 目的
第六思考実験と同一の `transition(z)`、11相機械、数値解像度を保持し、`ab` ではなく `aa` / `bb` をそれぞれ独立した自己相互作用系として与えたときの状態発展を検証する。

## 規格化

- `G = q0 = c = 1`
- `q0 = e_elem/3` を電荷格子単位とする
- 第六論文補遺で確定した `M/q0 = 20/3`
- 等質量なので `m_a=m_b=M/2=10/3`
- `N=nu=1/4`

電荷を `q = sign*n*q0` とすると、

`lambda = q/m = sign*(3n/10)`

自己相互作用では両脚が同一なので、

`C_self = 1 - lambda^2`

`D_self = (lambda-lambda)^2 = 0`

したがって Paper 6 の5条件を aa / bb に展開すると、n=1 は `(C,D)=(0.91,0)`、n=3 は `(0.19,0)`、n=4 は `(-0.44,0)` となる。

## 実験方針

Paper 6 と同じ5条件ラベルを維持し、各条件を `aa` と `bb` の2チャネルへ展開した。合計10ラベルである。ただし、状態写像はケースラベルを持たず、同じ初期状態からは同一軌道しか生成しないため、重複する軌道を再計算して巨大な raw を複製しない。

- n=1 の全 aa/bb: 同一初期状態 `(P,N,C,D)=(50,0.25,0.91,0)`
- n=3 の全 aa/bb: 同一初期状態 `(50,0.25,0.19,0)`
- n=4 の aa/bb: `C=-0.44` で Paper 6 の実数準円写像の定義域外

## 新規ローカル実行

### n=1
`P=50 -> 20` を新規に最後まで実行した。

- macro steps: 1456531
- cycles: 364.13275
- final P: 19.999971790606857
- `D=0` のため leading EM dipole power は全点で厳密に0

この結果は保存済み Paper 6 `n1_repulsive` と macro step 数・final P・初期状態が完全一致した。

### n=3
新規に最初の100,001 macro steps を実行した。

- P after prefix: 49.88204444925427
- D=0
- first macro microstate CSV は保存済み Paper 6 `n3_repulsive` と byte-identical

n=3 の全走行は Paper 6 `n3_repulsive` と初期状態も `transition(z)` も完全に同一であるため、既存 full raw trajectory を正準結果として再利用する。

Paper 6 full result:

- macro steps: 15266914
- final P: 19.999998005124105
- phi at P=20: 23981.211759437156
- first nonfinite time readout step: 6316005

これは近似的な代用ではない。同一 full state と同一決定論的 transition を用いるため、aa/bb の n=3 自己相互作用は Paper 6 n3_repulsive と同一の数値列である。

## n=4

`lambda=±1.2` より `C_self=1-1.2^2=-0.44`。変更していない Paper 6 の `_rates` は `sqrt(C)` を使用するため、実数準円写像の定義域外である。新しい規則を追加せず、`not_run_domain_C_le_0` として記録した。

## ファイル

- `run_paper7_self_interaction.py`: Paper 6 transition を維持した自己相互作用実験コード
- `results/case_matrix.json`: Paper 6 5条件 × aa/bb の10ラベル
- `results/self_n1/`: n=1 新規フル走行
- `results/self_n3_fresh/`: n=3 新規100k-step prefix走行
- `reference_paper6/`: 同一状態に対応する Paper 6 の保存済み summary
- `all_aa_bb_case_results.csv`: 10ラベルの結果表
- `verification_against_paper6.json`: 新旧一致検証
- `SOURCE_PROVENANCE_ja.md`: データ出所
- `SHA256SUMS.json`: 再現性用ハッシュ

## 再現

```bash
python run_paper7_self_interaction.py
python build_repro_package.py
```

既定実行では n=1 をフル走行し、n=3 は100,000 macrostep の新規 prefix 検証を行う。n=3 も最後まで新規再走行する場合は `python run_paper7_self_interaction.py --full-n3` を用いる。ただし既存 Paper 6 の n3_repulsive raw は同一初期状態・同一 transition による正準 full trajectory である。
