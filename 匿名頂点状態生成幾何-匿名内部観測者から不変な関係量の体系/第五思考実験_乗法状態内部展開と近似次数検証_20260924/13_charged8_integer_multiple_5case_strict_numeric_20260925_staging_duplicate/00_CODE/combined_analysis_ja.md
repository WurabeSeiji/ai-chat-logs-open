# 整数倍クーロン5ケース — 厳格8状態数値実験まとめ

5ケースすべてで experiment 11 と同一の `transition(z)` / `macrostep(z)` / 8永続状態 / 11相機械を使用した。変更したケース依存物理量は初期状態 $C,D$ のみ。解析基準は生成完了後の比較にのみ使用した。

| case | C | D | macro steps | |Δphi| at r=20 | |Δt| at r=20 | time readout |
|---|---:|---:|---:|---:|---:|---|
| n1_attractive | 1.09 | 0.36 | 488,345 | 8.436e-11 | 1.464e-08 | finite |
| n1_repulsive | 0.91 | 0.0 | 1,456,531 | 1.305e-09 | 3.203e-07 | finite |
| n3_attractive | 1.81 | 3.24 | 66,688 | 4.229e-09 | 8.739e-07 | finite |
| n3_repulsive | 0.19 | 0.0 | 15,266,914 | 3.653e-09 | N/A | nonfinite from step 6,316,005 |
| n4_attractive | 2.44 | 5.76 | 33,396 | 8.323e-09 | 1.545e-06 | finite |

## 主要所見

- `n1_attractive` は既存 experiment 11 と同一ケースで交差読み出し値を再現した。
- 同符号の `n1_repulsive`, `n3_repulsive` では $D=(lambda_A-lambda_B)^2=0$ のため、現行 LO 表現の EM dipole 放射項はゼロになる。Coulomb 保存力がゼロという意味ではない。
- `n3_repulsive` は $C=0.19$ で非常に長い走行となり、Q/t 読み出しだけが途中で有限精度範囲を超えた。P/H/N/C/D の生成は同一写像で r=20 まで継続した。
- `n4_attractive` は C=2.44, D=5.76 で最も少ない macrostep で r=20 に到達した。

## raw データ

全 macrostep 行を HDF5 part ファイルに保存した。HDF5 化と part 分割はシリアライズだけであり、状態更新則は変更していない。`n3_repulsive` の part 境界では完全41成分状態だけを checkpoint し、次 part でそのまま再開した。
