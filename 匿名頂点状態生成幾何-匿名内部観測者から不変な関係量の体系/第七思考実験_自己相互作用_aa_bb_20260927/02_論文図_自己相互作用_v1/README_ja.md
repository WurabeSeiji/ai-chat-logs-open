# 第七思考実験・独立自己相互作用 論文図 v1

## 目的

Paper 6 と同一の leading-order 準円軌道モデルを、独立な自己相互作用 `aa` / `bb` に適用した結果を論文化用に図化する。

規格化:

- `G=q0=c=1`
- `M/q0=20/3`
- `N=1/4`

自己相互作用では

- `n=1`: `|lambda_self|=0.30`, `C=0.91`, `D=0`
- `n=3`: `|lambda_self|=0.90`, `C=0.19`, `D=0`

であり、`aa` と `bb` は同じ初期状態を生成するため、同じ軌道を与える。

## 図

- `figure01_self_orbits_n1_n3_same_scale.svg`
  - n=1, n=3 の自己相互作用軌道を同一空間スケールで表示。
- `figure02_self_radius_vs_cycles.svg`
  - 半径読み出し P と累積周回数。
- `figure03_self_gw_and_em_power.svg`
  - leading GW quadrupole power と leading EM dipole power。
  - `D=0` のため EM dipole power は全域で厳密に 0。
- `figure04_self_generator_reference_checks.svg`
  - fresh strict generator の sample と Paper 6 closed-form reference の比較。
  - n=1 は全走行 sample、n=3 は fresh 100,000 macrostep prefix を使用。

## 再現

```bash
python generate_paper7_self_interaction_figures_v1.py \
  --data-root <experiment-root> \
  --outdir figures
```

このプログラムは各図を SVG と PNG の両方で生成する。Google Drive には論文用 SVG を保存する。
