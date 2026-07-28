# Stage G report

## 結論

relational_C1一候補をStage Fの独立System Aコピーへ実装し、単体検証、C0再現、代表条件2点、既存31系列を完了した。

中心判定は否定結果である。代表条件では \(\Gamma=1/32\)、31系列では \(\Gamma=1/2\) が保存され、`Delta Gamma`と`Delta R_eff`はいずれも事前閾値 `1e-10` 未満だった。したがってrelational_C1も今回の対称System A条件では一定Rへの再パラメータ化に退化した。

## 分類件数

- `constant_relation_reparameterization`: 12
- `relational_term_inactive`: 6

## コード上の事実

- 散乱・経路・rawノルムは512×16全状態で計算した。
- 関係波は由来別eta射影・搬送波除去後にコヒーレント合成した。
- 新規散乱候補はrelational_C1だけである。
- C0とreversed_C1はStage F系列を誤差0で保持した。

## 数学的帰結

- 一定パリティ区間のreversed_C1は一定Rへ退化する。
- Cauchy–Schwarz不等式から `0<=Gamma<=1`。
- 等ノルム・実重なりのGram行列とSystem Aの対称散乱行列はともに `I` と `sigma_x` の線形結合で可換なため、Gammaは角度が状態依存でも保存される。

## モデル定義

\[
\theta_{eff}=\theta_0-\kappa\rho(\theta_0)\bar c\,\Gamma_{AB}。
\]

## 作業仮説の判定

「パリティ符号と関係強度が散乱角を制御する」はモデルとして実装済み。「同じ純パリティ区間でも関係強度が変化して散乱率が動く」は、今回のSystem A対称更新では成立しなかった。

## 数値観察

relational_C1はC0ともreversed_C1とも異なる一定Rを与え、局在性交換と帰還位置を移動させた。しかしその変化は動的Gammaではなく一定Gammaによる。

## 未導出

κ、rho、重なり二乗を相互作用強度と読む根拠、自然界のボゾン・フェルミオン対応は未導出である。

## 棄却・保留

本条件における `dynamic_relation_dynamic_scattering` は棄却。別候補、別関係量、非対称散乱、N体系、論文反映は保留し、自動継続しない。


## 完了状態

参照正本12ファイルのpath・size・mtime・SHA-256は事前・事後で全件一致した。

```text
Stage G 完了。
既存System A / System B原本は変更していない。
新規実装はrelational_C1一候補のみ。
Candidate 2・3の追加実装は行っていない。
N体系へは組み込んでいない。
論文本文は変更していない。
人間の承認待ち。
```
