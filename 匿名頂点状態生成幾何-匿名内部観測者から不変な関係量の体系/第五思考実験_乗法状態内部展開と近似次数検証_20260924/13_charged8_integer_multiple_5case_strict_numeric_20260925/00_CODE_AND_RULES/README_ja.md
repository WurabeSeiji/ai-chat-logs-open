# 5ケース厳格数値実験 — 保存構成

この系列は、既存の `11_荷電8状態_厳格一回写像実験_20260925` と同じ
8永続状態 `(U,P,E,H,Q,N,C,D)`、11相機械、`transition(z)` / `macrostep(z)` を用いる。

変更したケース依存物理量は、初期状態に格納する `C,D` のみである。
生成側は解析基準を読まず、解析基準との比較は生成完了後だけに行う。

## raw data

全 macrostep を HDF5 (`raw_macro_partNNN.h5`) に保存する。
各行は:
`step,P,E,U_re,U_im,H_re,H_im,Q,N,C,D,t_readout,x_readout,y_readout,q_index`

HDF5 は保存形式であって状態更新式ではない。
`n3_repulsive` は実行環境の単回実行時間上限に対応するため4 partに分割した。
part間で持ち越した物理状態は完全41成分 full state そのものだけであり、
別の物理量・参照誤差・外部力情報は追加していない。

## figures

各ケースについて PNG と SVG を両方生成した。
Google Drive にはユーザー指示に従い PNG を保存せず、SVGのみ保存する。

## 既知の数値的限界

`n3_repulsive` では、軌道状態 `P,H` は `P=20` まで生成できた一方、
乗法時刻状態 `Q` が有限精度範囲を超え、step 6,316,005 以降の
`t = TAU0 log(Q)` 読み出しが非有限となった。
この挙動は修正せず実験結果として保存する。
