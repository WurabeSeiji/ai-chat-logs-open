# 第八思考実験 8A — 最小共面位相差結合の仮定探索

## 目的

Paper 7 までで独立に扱った `aa`, `ab`, `bb` を、最小の交差結合仮定で結び、GW/EM 放射振幅が同時に相殺される位相差が存在するかを探索する。

本実験は **仮定探索** であり、結合則の物理的正しさを先に仮定しない。最初の最小仮定が相殺条件を持つかどうかだけを検査する。

## 今回だけの明示仮定

1. `aa := a^2`, `ab := a*b`, `bb := b^2` と定義する。
2. 3軌道面の法線は同一方向と仮定する。
3. `phi_a=+delta/2`, `phi_b=-delta/2` と置き、`phi_aa=delta`, `phi_ab=0`, `phi_bb=-delta` とする。
4. 交差結合係数は

   `K(delta)=[[1,cos(delta),cos(2delta)],[cos(delta),1,cos(delta)],[cos(2delta),cos(delta),1]]`

   とする。
5. この最初の screen では K は **放射振幅層** に適用し、Paper 7 の局所 `transition(z)` へはまだ feedback しない。
6. GW は quadrupole phase として `2*phi`、EM dipole は `phi` を使う。

## 初期条件

Paper 6 補遺の `G=q0=c=1`, `M/q0=20/3` を継承し、n=1 異符号 `ab` を用いる。

- aa: C=0.91, D=0
- ab: C=1.09, D=0.36
- bb: C=0.91, D=0

123実数成分（41成分 x 3チャネル）の初期状態も保存した。初期非零成分は30個で、残りは0である。

## 探索

`delta in [0,2*pi]` を 360001 点で走査。評価半径は `r=50,40,30,20`。

結果として全評価半径で

- `delta = 2*pi/3`
- `delta = 4*pi/3`

に同時相殺点が現れた。数値残差は double precision の丸め誤差水準（概ね 1e-16）である。

`delta=2*pi/3` では

`cos(delta)=cos(2delta)=-1/2`

なので K の各行和・各列和が0となる。したがって今回の仮定では、合成放射振幅のゼロは特定の放射強度比に依存せず K の構造自体から生じる。

これは **この仮定に相殺 sector が存在する** ことを示すだけであり、結合則が物理的に正しいことや一意であることは示さない。

## 再現

```bash
python run_paper8_8A_minimal_phase_coupling.py
```

主要出力:

- `results/summary.json`
- `results/delta_scan_compact.npz`
- `results/delta_scan_minima.csv`
- `results/initial_state_nonzero_2pi_over_3.csv`
- `results/initial_state_123_*.npy`
- `results/coupled_state_123_*.npy`
- `results/SHA256SUMS.json`

`paper7_reference_transition.py` は Paper 7 で用いた局所 transition の出典固定用コピーであり、本 8A では変更していない。
