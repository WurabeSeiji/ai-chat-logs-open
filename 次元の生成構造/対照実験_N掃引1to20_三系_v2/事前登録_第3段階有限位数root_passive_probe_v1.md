# 事前登録 — 第3段階・有限位数 root passive probe v1

作成日: 2026-08-10  
状態: コード作成・smoke・非侵襲性監査まで。本走行は未許可・未実施。

## 1. 問い

有限位数共鳴論文の外生パラメータ

$$
\rho\equiv R_{\mathrm{scat}},\qquad
U_\rho^n=I
$$

と、N体母体が各stepで状態から生成する頂点係数を混同せず、次を測る。

1. 非線形頂点が実際に使う直前の関係辺別 $R_e(\Psi_t)$ は、有限位数root付近を訪れるか。
2. seed型はrootの数値住所を動かすのではなく、root付近を訪れる辺数・時刻・滞在を変えるか。
3. 大域診断量 `r_nopump` のroot横断と、実際の頂点係数 $R_e$ のroot訪問は一致するか。

本段階はroot **visibility/visitation** のprobeである。N体力学内の
$U_\rho^n=I$、Jacobian、二チャネル不変部分空間、monodromyは測らない。

## 2. 正本と変更禁止

- 母体runner: `run_nsweep_three_series_v2.py`
- 力学: `../統一万能関数_v1/unified_interaction_v1.py` の F v1
- F/D/G/S、`make_parent`、親 `seed=2`、閾値、時間刻み、シードrecipeは変更禁止。
- probeは新規 `run_finite_order_root_probe_stage3_v1.py` から母体を派生ロードする。
- F v1本体と母体runnerは編集しない。
- 通常の母体NPZへprobe配列を追加しない。別NPZへだけ保存する。

## 3. 三種類の量の定義

### 3.1 有限位数論文の $\rho$

$$
\rho_{n,m}=\cos^{2}\left(\frac{\pi m}{n}\right),\quad
U_{\rho_{n,m}}^{n}=I
$$

固定した同じ二チャネル作用素を反復する場合のパラメータである。

### 3.2 母体の関係辺別 $R_e$

母体 `_readout()` は、関係辺 $e=(i_e,j_e)$ ごとに、その辺自身と両端へ隣接する
関係辺の奇数帯・偶数帯パワーを集約し、

$$
R_e=\mathrm{scale}\,\sin^{2}\left(\mathrm{atan2}
\left(\sqrt{P_{F,e}},\sqrt{P_{B,e}}\right)\right)
$$

を返す。本実験の `strength_scale` は正本の既定値 `1.0` であり、この条件下で
$R_e=\sin^2\theta_e$ となる。返値の長さはnode数 $N$ ではなく、

$$
M=N(N-1)/2
$$

である。従来図中の「node別」は添字の呼称として不正確であり、本probeでは
**relation-edge / 関係辺別**と呼ぶ。

probe時刻は各stepの

```text
共有線形部 → R_e読出し → 非線形頂点 → post-step通常記録
```

の中央、すなわち頂点が消費する同じ配列を親 `_readout()` からコピーした時点である。

### 3.3 `r_nopump`

$$
r_{\mathrm{nopump}}=
\frac{P_{\mathrm{odd}}}
{P_{\mathrm{odd}}+P_{\mathrm{even,np}}}
$$

は全関係辺・全毛を集約し、pump帯を分母から除いた受動診断量である。
力学へ入力されず、$R_e$ または $\rho$ と等置しない。

## 4. 独立に保持する比較値

| 名称 | 値 | 分類 |
|---|---:|---|
| $\rho_{124,23}$ | `cos²(23π/124)=0.6971779275566593` | 有限位数root、n=124 |
| $\rho_{620,117}$ | `cos²(117π/620)=0.6878251911141453` | 有限位数root、n=620 |
| $\rho_{M_Z}$ | `0.687822933884774` | 物理対応比較値。有限位数rootとは置かない |
| $\rho_{122,23}$ | `cos²(23π/122)=0.6883639468175926` | 近傍の有限位数root対照 |

$$
\rho_{620,117}-\rho_{M_Z}=2.2572293713\times10^{-6}.
$$

この二値を判別したと呼べる最小条件は、関係辺の実測値が二値の半間隔
$1.1286146857\times10^{-6}$ より近くへ入ることである。線形補間だけでは実測判別と呼ばない。

数値shamは事前固定する。

- $\rho_{124,23}\pm10^{-3}$
- $\rho_{620,117}\pm10^{-3}$
- $(\rho_{620,117}+\rho_{M_Z})/2$

有限位数rootは稠密なので、`±1e-3` shamを「rootが存在しない点」とは主張しない。
物理的・論文上の候補として選んでいない、対称な数値対照である。

## 5. 記録配列

各matter走行とその内部vacuum対照について次を別NPZへ保存する。

- `step`: 1…T
- `edge_index`: 0…M−1
- `edge_ia`, `edge_ib`: 関係辺の両端node
- `R_prevertex_matter[T,M]`
- `R_prevertex_vacuum[T,M]`
- `r_nopump_prevertex_matter[T]`, `r_nopump_prevertex_vacuum[T]`
- 上記の分子・分母を監査する `odd_power_prevertex_*` と
  `even_power_nopump_prevertex_*`
- 各stepでpre-vertex読出しが1回だけだったことの監査配列

通常の母体NPZは従来と同じkey・dtype・shape・値で保存する。
`r_nopump_prevertex_*` は $R_e$ と同じ状態から同期計算するが、頂点が使う
係数ではなく、定義の違う大域集約診断量である。

## 6. 非侵襲性の合格条件

本走行へ進む前に、短い同条件probeを既存の非instrumented NPZ prefixと比較する。

1. referenceのkeyがcandidateにすべて存在し、全reference配列でdtype一致。
2. shape一致、または宣言された時系列keyに限り、短いsmokeの第0軸だけが
   既存長時間走行のprefix。静的配列のprefixは許さない。
3. NaN配置一致。
4. C連続バイト列が全byte一致。
5. `R_prevertex_*` のshapeが `(T,M)`。
6. 全値が有限かつ $0\le R_e\le1$。
7. 各stepで頂点前読出しを厳密に1回捕捉。
8. 母体runnerとF v1の走行前後SHA-256不変。

一つでも外れれば、本走行を開始しない。
旧版referenceが現行母体の追加keyを持たない場合は、追加key名を明示した
`reference_complete_subset`とし、「同一key集合」とは呼ばない。別に現行の同一source・同一
`T=8`非instrumented対照を作り、`--require-identical-key-set`で全key合格を要求する。

## 7. 最小pilot行列

全条件は既存実験と同じ `N=12, T=42000, Nn=16, Neta=8` とする。

| 順序 | mode | δ | 役割 |
|---:|---|---:|---|
| 1 | mixed | 0.03162277660168379 | `r_nopump`が128側二値を通るが124側へ未到達 |
| 2 | mixed | 0.04357 | `r_nopump`が128側・124側をともに横断 |
| 3 | neutral | 0.03162277660168379 | 保存済みpost-step辺別Rが候補へ最接近した型 |
| 4 | electron | 0.03162277660168379 | 同じ立上り時刻・異なる内部住所の対照 |
| 5 | fermion_family | 0.03162277660168379 | F5型。mixedとの奇数seed power一致対照 |
| 6 | boson_family | 0.03162277660168379 | odd不在・$R_e=0$の負対照 |

逐次実行し、各走行後に報告して停止できる構造にする。mixed–F5同一δ比較は
$P_F$ は一致するが、$P_B$ と $A_{coh}$ が同時に異なるため、ボゾン組成だけの因果対照とは呼ばない。

## 8. 解析量

関係辺ごとに独立して、各比較値への以下を出す。

- edge-time全体の最小絶対距離、そのstep、edge、両端node、$R_e$
- 上向き・下向きの厳密な両側横断数、横断したedge数、最初・最後の線形補間step
  （中心値での接触後に元の側へ戻る接触は横断数から除く）
- 半幅 `1e-2,1e-3,…,1e-8` のedge-step滞在数
- その帯へ入ったedge数、任意edgeが帯内にいる時間step数
- 単一edgeの最長連続滞在、任意edgeの最長連続滞在
- 有限位数rootだけ、まずその名前付き $\rho_{n,m}$ へ最も近い実測sampleを
  $R$ 距離で選び、そのsampleで

$$
\left|e^{in\omega(R_e)}-1\right|,\quad
\omega(R)=\pi+2\arcsin\sqrt{R}
$$

と、選択rootへの一step位相距離を計算する。

最後の量はscalar root proximityであり、N-step作用素積ではない。また
$n$重閉包defectだけを全sample上で最小化すると、同じ $n$ の別の $m$-rootを拾う。
そのため必ず「名前付きrootへ最近のsample」に限定し、閉包defectだけで
$m$ を識別できるとは主張しない。

同期記録した `r_nopump_prevertex` も別系列として、同じ比較値ごとの最近実測
step・厳密横断・滞在をJSON/CSVに出す。これと $R_e[t]$ は同じpre-vertex
状態である。一方、通常母体NPZの `rec_*_r_nopump[t]` はstep後であるため、
同時性の対照には用いず、時刻定義の差を測る別診断として併記する。

## 9. 判定規則と主張分類

### 9.1 現在すでに棄却されている同一視

`r_nopump = 有限位数論文のrho` は定義・分母・作用先が違うため削除する。

### 9.2 pilotで判定できるもの

- `r_nopump`が候補値を横断してもpre-vertex $R_e$ が候補へ入らない
  → 大域比の横断をroot訪問とする写像を反証。
- B3で $R_e=0$ が全辺・全stepで維持
  → odd sectorがvisibilityに必要という負対照を支持。
- 有限位数候補の滞在・減速が両側shamより複数modeで再現して大きい
  → root-specific visibilityの**候補**として保持。物理的α同定ではない。
- 候補と対称shamが同程度
  → 滑らかな通過で説明でき、root specificityを支持しない。
- $\rho_{620,117}$ と $\rho_{M_Z}$ の半間隔内へ実測点がない
  → この時間標本では両者を判別不能。

### 9.3 pilotで判定できないもの

- N体作用素が有限位数である。
- $U_\rho^n=I$ がN体軌道を生成した。
- α=137/128の物理的選択原因。
- seedがrootの住所自体を動かした。

これらには、意味ある固定二チャネル射影、射影外漏れ監査、接線作用素または
monodromyが必要である。射影を構成できない限り、名目だけのJacobian実装は行わない。

## 10. 実行ゲート

1. AST/help/describe-only。
2. `T=8` smokeを一条件だけ実行。
3. 既存NPZ prefixと全reference-key byte一致を監査。
4. 現行同一sourceの非instrumented `T=8` 対照と、同一key集合・全byte一致を監査。
5. probe解析器のJSON/CSV smoke。
6. 結果を報告し、本走行の許可を待つ。

本事前登録の作成だけでは、§7のpilot本走行を許可しない。
