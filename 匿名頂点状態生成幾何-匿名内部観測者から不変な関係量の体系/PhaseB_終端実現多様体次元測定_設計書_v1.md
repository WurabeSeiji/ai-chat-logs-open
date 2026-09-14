# Phase B 設計書 v1
## 匿名 Bob-star 完全不変量フレームによる終端実現多様体次元の直接測定

**系列:** 匿名頂点状態生成幾何  
**著者:** 木原範昭 / Noriaki Kihara  
**ORCID:** 0009-0004-6753-4020  
**版:** v1.0 / 2026-09-14  
**凍結基準コミット:** `5944ad80`  
**対象:** 保存済み自己無撞着関係波データ  
**目的:** `dim(M_terminal(N)/G_B) ?= 4` の直接測定  
**重要規約:** 力学正本 `run_N3_N40_stage123_v1.py` には一切変更を加えない。Phase B は保存済みデータの読出し専用実験とする。

---

# 0. 本設計の位置づけ

generic null 側はクローズ済みである。

確定済み事項：

1. Bob-star generic quotient の実次元

$$
\dim_{\mathbb R}(\mathcal S_B^{\rm gen}/G_B)=2N-4.
$$

2. 完全不変量フレーム

$$
\Phi_B^{(N-1)}
$$

は \(P_1\neq0\) chart 上で Bob-star を \(G_B\) ゲージを除いて完全に分離する。

3. generic 複素ガウス null では

$$
\kappa^2\sim \operatorname{Beta}\left(1,\frac{m-1}{2}\right),\qquad m=N-1.
$$

4. N=2〜40 の generic null 対照実験はクローズ済み。

5. 高 N では Newton/root による明示的 star 復元は数値的に悪条件化するため、Phase B では根復元を使わず、\(\Phi_B\) を直接座標として用いる。

従って Phase B の判別問題は一つに固定する。

$$
\boxed{
d_{\rm terminal}(N)
=
\dim\left(\mathcal M_{\rm terminal}(N)/G_B\right)
\stackrel{?}{=}4
}
$$

---

# 1. 主仮説と反証条件

## H_B1 — rank-4 飽和仮説

既存力学で外部解析から観測された rank-4 が、Bob 内部観測幾何にも対応するなら、

$$
\boxed{
d_{\rm terminal}(N)\approx4
}
$$

が十分大きい N に対して N 非依存に成立するはずである。

## 反証条件

以下のいずれかが安定に観測されれば H_B1 を棄却または修正する。

- \(d_{\rm terminal}(N)\) が N とともに増加する。
- 座標系 A/B、標準化、推定器を変えると推定次元が大きく変わる。
- \(K\to N-1\) で次元が飽和せず増加し続ける。
- 最終 step のみでは4だが、seed/Bob を変えると大きく崩れる。
- \(P_1\) chart 境界が高頻度で現れ、単一 chart の測定が成立しない。

---

# 2. 主標本の定義

主推定は T001/T003 の凍結仕様に従う。

固定した N ごとに、標本集合を

$$
\mathcal D_N
=
\left\{
\Phi_B(z_{N,s,B,t_{\rm final}})
\right\}_{s,B}
$$

とする。

ここで：

- \(N\)：固定
- \(s\)：seed / independent run
- \(B\)：Bob 頂点
- \(t_{\rm final}\)：各走行の最終 step
- 時系列途中は主標本に混ぜない

**異なる N を同一 PCA / intrinsic-dimension 推定に混在させない。**

主推定単位は N ごと。

---

# 3. Bob-star の取り方

各保存済み終端状態から、Bob を1頂点固定し、Bob を端点とする incident relation の無順序多重集合

$$
\mathcal Z_B
=
\{z_{B1},\ldots,z_{B,N-1}\}
$$

を抽出する。

辺の外部 index は読出し時のみ利用し、\(\Phi_B\) 計算では順序に依存しない冪和を用いる。

---

# 4. 不変量フレーム

## 4.1 正規化冪和

$$
P_k
=
\frac{\sum_jz_j^k}
{\left(\sum_j|z_j|^2\right)^{k/2}}.
$$

位相基準として \(P_1\) を使う chart では、

$$
\phi_k
=
\arg\left(P_k\overline{P_1^k}\right),
\qquad k\ge2.
$$

主たる完全フレームは

$$
\Phi_B^{(K)}
=
\{|P_1|\}
\cup
\{|P_k|,\phi_k\}_{k=2}^{K}.
$$

完全性の上限：

$$
K=N-1.
$$

---

# 5. 測定座標 A/B 二系統

## A系 — 直接不変量系

各 \(K\) について、

$$
X_A^{(K)}
=
\left(
|P_1|,
|P_2|,\ldots,|P_K|,
\cos\phi_2,\sin\phi_2,\ldots,
\cos\phi_K,\sin\phi_K
\right).
$$

**位相角 \(\phi_k\) を生の角度として PCA に入れない。**  
\(-\pi/\pi\) 境界を避けるため必ず \((\cos\phi_k,\sin\phi_k)\) に展開する。

## B系 — log 振幅系

$$
X_B^{(K)}
=
\left(
\log(|P_1|+\epsilon),
\log(|P_2|+\epsilon),\ldots,\log(|P_K|+\epsilon),
\cos\phi_2,\sin\phi_2,\ldots,
\cos\phi_K,\sin\phi_K
\right).
$$

終端リング状態では高次 \(|P_k|\) が急減する可能性があるため、B系は小振幅情報を保持しやすい。

---

# 6. 標準化方式

主推定では、各 N・各 K・各座標系ごとに**標本集合内で成分標準化**する。

標準方式 S1：

$$
x_j'
=
\frac{x_j-\mu_j}{\sigma_j}.
$$

ただし、

$$
\sigma_j<\tau_\sigma
$$

のほぼ定数成分は除外し、その除外数を必ず記録する。

感度検査として：

- S0：無標準化
- S1：平均0・標準偏差1
- S2：robust scaling（median / MAD または IQR）

を比較する。

**主結果は S1。S0/S2 は感度検査。**

---

# 7. epsilon の規約

B系では

$$
\log(|P_k|+\epsilon)
$$

を使う。

\(\epsilon\) はデータを見て最適化しない。

主値：

$$
\boxed{\epsilon=10^{-12}}
$$

感度検査：

$$
10^{-10},\ 10^{-12},\ 10^{-14}.
$$

倍精度の数値床より十分上/近傍を跨ぐ範囲とする。

---

# 8. P1 chart 境界

各 Bob-star で

$$
|P_1|
$$

を必ず保存する。

主 chart 条件：

$$
|P_1|>\tau_{P1},
\qquad
\tau_{P1}=10^{-8}.
$$

報告項目：

- N ごとの \(\min|P_1|\)
- \(|P_1|<10^{-3}\) 比率
- \(|P_1|<10^{-6}\) 比率
- \(|P_1|<10^{-8}\) 比率

主解析は \(P_1\) chart 内で行う。

chart 境界標本が無視できない場合のみ、別の非零 \(P_r\) を位相基準とする atlas を Phase B.1 として追加する。最初から複雑化しない。

---

# 9. K 飽和曲線

各 N について、

$$
K=2,3,\ldots,N-1
$$

を走査する。

各 K で intrinsic dimension 推定値

$$
d_{\rm est}(N,K)
$$

を計算する。

必須出力：

- \(d_{\rm est}\) vs K
- estimator ごとの曲線
- A/B 系統比較
- S1 主推定
- 最終値 \(K=N-1\)

rank-4 仮説が成立する場合の期待形：

$$
d_{\rm est}(N,K)
\rightarrow 4
$$

が比較的小さい \(K_*\ll N\) で飽和し、その後ほぼ一定。

重要：飽和 K を理論から先に仮定しない。

---

# 10. 次元推定器

主推定は単一手法に依存させない。

## E1 — PCA / SVD 有効 rank

標準化座標行列 X に対し SVD。

報告：

- 特異値全列
- 累積寄与率
- 95%, 99%, 99.9% 寄与率で必要な成分数
- gap spectrum

ただし PCA 次元は線形包絡次元なので、**単独で intrinsic dimension と同一視しない。**

## E2 — local PCA

各点の kNN 近傍で局所 PCA。

近傍サイズ：

$$
k_{\rm nn}\in\{10,20,30,50\}
$$

を感度検査。

局所 rank の中央値・四分位範囲を報告。

## E3 — TwoNN / nearest-neighbor intrinsic dimension

局所幾何に比較的低仮定な intrinsic dimension estimator として採用。

## E4 — Levina-Bickel MLE

kNN ベース MLE intrinsic dimension。

複数 k を比較。

### 主判定

「4次元」とするには少なくとも、

- local PCA
- TwoNN または MLE

の双方が 4 近傍で一致し、A/B 座標系でも安定することを要求する。

PCA だけで4と判定しない。

---

# 11. 標本数と識別可能性

N 固定での標本数は

$$
n_N
=
(\text{seed 数})
\times
(\text{Bob 数})
$$

となる。

最初に保存データを棚卸しし、N ごとに

- seed 数
- Bob 数
- 有効標本数
- 重複状態の有無

を表にする。

**標本数が少なすぎる N は次元推定対象から除外し、参考値扱いにする。**

最低基準の暫定値：

$$
n_N\ge 10(d+1)
$$

を4次元仮説の最低条件とし、

$$
n_N\ge100
$$

を望ましい目安とする。

これは固定閾値ではなく、bootstrap 安定性と合わせて判定する。

---

# 12. Bob 間・seed 間依存への注意

同一走行内の Bob は完全独立標本ではない可能性がある。

従って bootstrap は二種類行う。

1. point bootstrap：Bob×seed を点として resample
2. seed-cluster bootstrap：seed 単位で resample し、その中の Bob をまとめて保持

主信頼区間は **seed-cluster bootstrap** を優先する。

---

# 13. 主推定

各 N について主推定仕様を固定する。

- step：最終 step のみ
- 座標：A系 / B系の両方
- 標準化：S1
- epsilon：\(10^{-12}\)
- chart：\(|P_1|>10^{-8}\)
- K：2〜N−1 全走査
- 推定器：E1〜E4
- bootstrap：seed-cluster 優先

最終主結果表：

| N | samples | P1境界率 | K* | PCA | local PCA | TwoNN | MLE | consensus dim |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|

---

# 14. 感度検査 — 主推定から分離

凍結済み6軸を主解析と混同しない。

## SENS-1 — 最終 step のみ
主推定そのもの。

## SENS-2 — 終端近傍 step 窓
例：

$$
t_{\rm final}-\Delta t,\ldots,t_{\rm final}
$$

を追加し、横断緩和方向が次元を上振れさせるか確認。

## SENS-3 — 座標 A/B
raw magnitude vs log magnitude。

## SENS-4 — 標準化方式
S0 / S1 / S2。

## SENS-5 — epsilon
\(10^{-10},10^{-12},10^{-14}\)。

## SENS-6 — 推定器
PCA / local PCA / TwoNN / MLE。

**感度検査の結果で主仕様を変更しない。**  
必要なら次版設計として明示的に改訂する。

---

# 15. κ_B の併行読出し

Phase B 主目的は次元測定だが、同じ最終標本について

$$
\kappa_B=|P_2|
$$

を必ず保存する。

各 N で、

- 平均
- 中央値
- 分散
- 最小
- 最大
- generic null exact mean
- generic null CDF に対する percentile
- \(\kappa<\epsilon_\kappa\) 比率

を報告。

generic null：

$$
F_\kappa(q)
=
1-(1-q^2)^{(m-1)/2}.
$$

これにより、終端状態が generic null からどれだけ逸脱しているかを次元推定と独立に測る。

---

# 16. Phase C との境界

Phase B では原則として最終 step を扱う。

時系列

$$
\kappa_B(t),\quad P_k(t)
$$

の追跡は Phase C。

ただし SENS-2 の終端窓だけは Phase B 感度検査として許可する。

インフレーション開始・モード転移などの時系列解釈は Phase B に持ち込まない。

---

# 17. 保存済みデータ確認フェーズ

実装前に、まず読み取り専用で以下を確定する。

1. 正本データフォルダ
2. N=6〜40 の存在範囲
3. seed 数
4. step 範囲
5. final step の定義
6. 辺データの列構造
7. 頂点番号と辺 index の対応
8. 複素値の保存形式
9. N ごとの欠損
10. 同一 run の重複コピー有無

この段階では計算しない。

`DATA_INVENTORY.md` と `data_inventory.csv` を作る。

---

# 18. 読出し専用実装規約

専用フォルダ案：

`PhaseB_terminal_manifold_dimension_v1/`

必須ファイル：

- `DESIGN_PhaseB_terminal_manifold_dimension_v1.md`
- `DATA_INVENTORY.md`
- `data_inventory.csv`
- `extract_bob_stars.py`
- `compute_invariant_frames.py`
- `estimate_intrinsic_dimension.py`
- `run_phaseB.py`
- `run_all.sh`
- `results/`
- `figures/`
- `analysis_PhaseB.md`
- `README.md`
- `SHA256SUMS.txt`

### 禁止

- `run_N3_N40_stage123_v1.py` の変更
- 正本 CSV の上書き
- 正本フォルダへの中間生成物混在
- N 混在 PCA
- Newton/root 再構成を主解析へ使用
- 結果を見て K / epsilon / estimator の主仕様を変更

---

# 19. 出力図

最低限：

1. `dimension_vs_N.png`
2. `dimension_vs_N.svg`
3. `dimension_vs_K_Nxx.png` 主要 N
4. `singular_spectrum_Nxx.png`
5. `coordinate_A_vs_B_dimension.png`
6. `kappa_terminal_vs_null_N.png`
7. `P1_chart_margin_vs_N.png`

必要に応じて全 N の K 飽和曲線 CSV を保存する。

---

# 20. 判定規則

## rank-4 支持

以下を満たす N 群で、

$$
d_{\rm consensus}(N)\in[3.5,4.5]
$$

程度に安定し、

- K 増加で4付近に飽和
- A/B 両系統で一致
- local PCA / NN estimator が一致
- seed-cluster bootstrap CI が4を含む
- N 依存増加がない

場合、rank-4 内部実現多様体仮説を支持する。

## rank-4 非支持

- 5以上へ安定
- N とともに増加
- K とともに増え続ける
- estimator 間で不整合

なら支持しない。

「4に近いから4」と丸めない。

---

# 21. Phase B 成功時の意味

もし

$$
\boxed{
\dim(\mathcal M_{\rm terminal}(N)/G_B)=4
}
$$

が N に依存せず成立すれば、

generic quotient

$$
2N-4
$$

の中から、自己無撞着力学が4実次元の内部可読部分多様体を選択していることになる。

これは既存シリーズの外部 rank-4 観測と、匿名内部観測者 Bob の幾何を初めて接続する。

---

# 22. Phase B 失敗時の意味

4にならなくても実験は失敗ではない。

- 2次元 → Bob-star がさらに強く縮退
- 6,8 → 外部 rank-4 と内部 rank は異なる
- N依存 → rank-4 は外部射影固有の可能性
- estimator 不安定 → 終端状態族の標本設計不足または多様体でない可能性

を判別できる。

従って Phase B は結果に依存せず有効な判別実験である。

---

# 23. 実行順序

1. 本設計書の査読・凍結
2. 保存済みデータ棚卸し
3. DATA_INVENTORY の査読
4. 読出し専用コード作成
5. N=6 のみ dry-run
6. N=6〜10 pilot
7. pilot 検査
8. N=6〜40 本走行
9. 主推定
10. 感度検査
11. 独立対照走行
12. 発見台帳更新
13. コミット
14. 論文3 作成

---

# 24. 現時点では実行しないこと

本設計の査読完了までは：

- 保存データの数値解析を開始しない
- Phase B コードを書かない
- 既存力学正本へ触れない
- rank-4 を前提に結果を解釈しない

---

# 結論

Phase B は、

$$
\boxed{
\text{N固定}
\times
\text{最終step}
\times
\text{seed}
\times
\text{Bob}
}
$$

を主標本とし、

$$
\boxed{
\Phi_B^{(K)}
}
$$

を root 復元せず直接座標として使う。

主判別量は、

$$
\boxed{
d_{\rm est}(N,K)
}
$$

の K 飽和曲線と、

$$
\boxed{
d_{\rm terminal}(N)
\stackrel{?}=4
}
$$

である。

座標 A/B、標準化、epsilon、step 窓、推定器は主推定から分離した感度検査とする。

この仕様を凍結してから、保存済みデータの棚卸しへ進む。
