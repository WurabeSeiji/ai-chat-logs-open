# standalone parent/census v1 実装検証記録

実施日: 2026-08-04<br>
身分: 実装スモークテスト。物理結論を与える本実験ではない。

## 1. 白色零閉鎖生成器

条件:

```text
N-body = 5
resolution = 16
seed = 40260804
steps = 64
source = white_null
storage = complex128
```

不変量監査:

```text
initial norm             = 1.0
initial |Z^T Z|          = 2.7887742898292585e-16
maximum norm drift       = 2.220446049250313e-16
maximum |Z^T Z|          = 4.791768344701908e-16
maximum sigma residual   = 9.809054598600952e-13
```

同一条件を別出力ディレクトリへ再実行し、次の二ファイルがSHA-256で一致した。

```text
trajectory.npy
55c05a1878a2ba3bee18955c00173c5441aca2c593802b31d7ab1799fbdf31d0

diagnostics.npz
240fc66616127f2aa11d006e394b20111a1aa87af8832c82cae8590e548bd618
```

従って、seed再現性と零閉鎖・ノルム保存を確認した。

## 2. 単一モード対照による読出し陽性対照

同じ $N=5,\mathcal N=16$, 64ステップを `single_mode` で生成し、生成器とは
別プロセスで一覧表プログラムだけを実行した。

主モード $k=1$ の判定:

```text
address                         = 1/16
state recurrence fidelity       = 1.0
half-cycle overlap              = -1.0 + 1.50e-14 i
single-mode closure residual    = 4.84e-14
status                          = exact_finite_order_closed_mode
spin structure                  = half_integer_internal_layer
physical spin                   = unresolved_parent_cartan_map
```

既知の単一円偏波を、厳密有限位数・零閉鎖・2:1内部被覆として検出しつつ、
物理スピン値へ昇格させていないことを確認した。

## 3. 白色条件による偽陽性対照

同じ短時間の `white_null` 出力では、DFT上の支配ビン $k=1$ は存在したが、
代表空間モードの単独閉鎖残差が $1.58\times10^{-2}$ であった。このため

```text
status = persistent_needs_bundle_closure
```

とされ、厳密有限位数・単独零閉鎖モード数は0になった。DFTビンを置いただけで
粒子を作ったと判定しないことを確認した。

## 4. 分解能144での動作確認

$N=5,\mathcal N=144$, 576ステップ、1 seed の短時間実行では、

```text
maximum norm drift       = 8.881784197001252e-16
maximum |Z^T Z|          = 2.991448473778328e-15
W                        = 10 <= M=10
exact closed modes       = 0
```

となった。これは方式1が分解能144でも不変量を保ち、$W\le M$ を機械的に守る
ことの実装確認に限られる。形成・不形成の物理判定には、固定した本実験時間、
複数seed、$N=5,40,300$ の比較が別途必要である。
