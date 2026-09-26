# 荷電状態内部化・厳格ルール再実装 2026-09-25

## 目的
旧 `09_荷電8状態閉包_対照実験_20260925` は、`macro_fast_eight()` という別経路を本実験に使用し、R2/R3/R11/R14 に違反した。旧フォルダは証拠としてそのまま残し、本フォルダでは別経路を作らず、保存済み `STRICT_EXPERIMENT_RULES_ja` に従ってゼロから再実装した。

## 実行ゲート
A1〜A8をすべてPASSした場合に限りA9数値実験を実行する。実際に監査器の初版は自己誤検出で停止し、その時点ではA9へ進まなかった。監査器を関数定義ベースへ修正後、A1〜A8を再実行し、全PASSを確認してからA9へ進んだ。

## 単一生成経路
動力学更新は `transition(z)` のみ。

$$
z_{n+1}=\sum_j q_j F_j(z_n)
$$

で、全候補は旧状態だけを読む。全状態成分を同じ one-hot 積和で一度に書く。`fast`, `collapsed`, `shortcut`, external-control など別の力学経路は存在しない。

## 状態
厳密な積・積和化のため、8状態に固定せず、平方根と逆数も状態化した。

- U: 位相状態
- X = sqrt(r)
- Y = 1/X
- H: 軌道位相状態
- Q: 時計の乗法状態
- N = eta
- S = sqrt(C)
- T = 1/S
- D
- RK内部状態 k1..k4, stage, delta
- one-hot phase q[0..8]

全体36複素スロット。物理状態数の最小性は主張しない。

## 結果
- macro steps: 488345
- microsteps: 4395105
- t(r=20) error vs saved independent reference: 3.9959559217095375e-08
- phi(r=20) error: 1.3369572116062045e-10
- max |XY-1|: 8.1046280797636427e-14
- max |ST-1|: 0
- N drift: 0
- D drift: 0

## 重要な限定
これは LO adiabatic quasi-circular charged benchmark の再構成であり、一般Einstein-Maxwell二体系や真の最小状態数を証明しない。
