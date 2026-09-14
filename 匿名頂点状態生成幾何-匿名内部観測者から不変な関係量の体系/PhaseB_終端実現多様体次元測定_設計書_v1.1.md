# Phase B 設計書 v1
## 匿名 Bob-star 完全不変量フレームによる終端実現多様体次元の直接測定

**系列:** 匿名頂点状態生成幾何  
**著者:** 木原範昭 / Noriaki Kihara  
**ORCID:** 0009-0004-6753-4020  
**版:** v1.1 / 2026-09-14  
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

既存力学で外部解析から観測された rank-4 は、単一軌道における幾何学的 span rank である。一方、Phase B が測定するのは seed×Bob を横断する終端 Bob-star アンサンブルの内部商次元であり、両者は定義上同一ではない。

したがって H_B1 は、

$$
\boxed{
\text{外部 single-orbit span rank}
\stackrel{?}{=}
\text{内部 terminal ensemble dimension}
}
$$

という**仮説的リンク**を検定するものである。

このリンクが成立するなら、

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

## 2.1 終端実現多様体の定義

固定した N に対して、seed と Bob を横断した終端 Bob-star の \(G_B\) 商軌道の閉包を

$$
\boxed{
\mathcal M_{\rm terminal}(N)
:=
\overline{
\left\{
[\mathcal Z_B(N,s,B,t_{\rm final})]_{G_B}
\right\}_{s,B}
}
}
$$

と定義する。

既存シリーズの終端では、等振幅化と共通 tail 構造を共有しつつ、seed により相対位相配置が選択される。従って Phase B が測るのは、自己無撞着力学が終端で選択する**真空配置アンサンブルの内部可読連続族の次元**である。

これは単一軌道の外部 span rank とは別概念であり、その一致・不一致自体が判別対象である。

## 2.2 主標本

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
\phi_{k|r}
=
\arg\left(P_k^{\,r}\overline{P_r^{\,k}}\right).
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
\cos\phi_{2|r},\sin\phi_{2|r},\ldots,
\cos\phi_{K|r},\sin\phi_{K|r}
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
\cos\phi_{2|r},\sin\phi_{2|r},\ldots,
\cos\phi_{K|r},\sin\phi_{K|r}
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

# 8. 位相基準 chart と縮退フォールバック規則

終端状態は generic null ではないため、\(P_1\) が十分大きいことを前提にしない。むしろ終端の放射状・近似多回対称配置では、

$$
\sum_j z_j \approx 0
$$

となり、

$$
|P_1|\ll1
$$

が系統的に現れる可能性がある。

さらに近似4回対称性が残る場合、

$$
z\mapsto iz
$$

に近い対称性から、

$$
P_1\approx0,\qquad P_2\approx0
$$

が同時に起こり得る。

従って pilot 後に atlas を追加するのではなく、以下のフォールバック規則を**事前登録**する。

## 8.1 位相基準の選択

各標本について、位相基準次数 \(r\) を

$$
\boxed{
r
=
\min\left\{
k\ge1:\ |P_k|>\tau_{\rm ref}
\right\}
}
$$

と定義する。

主閾値：

$$
\boxed{
\tau_{\rm ref}=10^{-8}
}
$$

とする。

位相不変量は \(P_r\) を基準として構成する。\(P_k\) の U(1) 電荷は \(k\) なので、単純な \(\arg(P_k\overline{P_r^{\,k}})\) は電荷が一致しない。従って一般 chart では、

$$
\boxed{
\phi_{k|r}
=
\arg\left(
P_k^{\,r}\overline{P_r^{\,k}}
\right)
}
$$

を用いる。これは U(1) 電荷 \(kr-kr=0\) であり、完全にゲージ不変である。

\(r=1\) の場合は従来どおり、

$$
\phi_{k|1}
=
\arg(P_k\overline{P_1^k})
$$

に戻る。

## 8.2 記録項目

各 Bob-star について必ず、

- \(|P_1|\)
- \(|P_2|\)
- \(|P_3|\)
- \(|P_4|\)
- 選択された \(r\)
- \(r\) が見つからない標本数

を保存する。

N ごとに \(r\) の分布そのものを報告する。

## 8.3 pilot の事前予測

保存データを見る前に、以下を pilot 予測として登録する。

終端の十字・4回対称構造が内部 star に残っているなら、

$$
\boxed{
|P_1|,\ |P_2|\ \text{は小さく、}\ |P_4|\ \text{は相対的に大きい}
}
$$

可能性がある。

その場合、

$$
\boxed{
r=4
}
$$

が高頻度で選ばれることを予測する。

これは atlas の都合ではなく、終端4回対称構造の**独立な内部観測テスト**である。予測が外れた場合も結果として報告する。

## 8.4 chart 感度

\(\tau_{\rm ref}\) の感度として、

$$
10^{-6},\ 10^{-8},\ 10^{-10}
$$

を別解析で確認する。

主仕様は \(10^{-8}\) から変更しない。

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

また、K が大きくなると埋め込み座標次元は概ね \(O(3K)\) で増加し、固定 N の標本数 \(n_N\) を上回り得る。従って高 K で推定次元が上昇・不安定化した場合、それを直ちに物理次元増加と解釈しない。各 K について

$$
rac{n_N}{D_{m ambient}(K)}
$$

を併記し、標本不足による estimator artifact の可能性を明示する。

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

ただし、同一走行内の Bob-star が対称性により数値的に一致または極近接する可能性がある。距離ゼロまたは数値分解能未満の重複点は TwoNN / MLE 前に除去し、N ごとに

- 原標本数
- 一意標本数
- 除去数
- 除去率

を報告する。

重複判定の主閾値は、標準化後 Euclidean 距離

$$
\boxed{
d<10^{-12}
}
$$

とし、\(10^{-10},10^{-12},10^{-14}\) を感度確認する。

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
- 位相基準：最小 \(r\) with \(|P_r|>10^{-8}\)（§8 の事前登録規則）
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
7. `reference_order_r_distribution_vs_N.png`
8. `P1_P2_P4_terminal_vs_N.png`

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


---

# 改稿履歴 v1.1

- \(\mathcal M_{\rm terminal}(N)\) を seed×Bob 終端 Bob-star の \(G_B\) 商軌道の閉包として明示定義。
- 外部 single-orbit span rank と内部 terminal ensemble dimension の一致は恒等ではなく H_B1 の仮説的リンクであることを明記。
- \(P_1\) 固定位相基準を廃し、\(|P_r|>\tau_{\rm ref}\) を満たす最小次数 \(r\) を用いる事前登録 atlas 規則へ変更。
- 一般参照次数に対する位相不変量を \(\phi_{k|r}=\arg(P_k^r\overline{P_r^k})\) と定義。
- pilot 事前予測として「終端で \(|P_1|,|P_2|\) が縮退し \(|P_4|\) が相対的に大きく、\(r=4\) が浮上する可能性」を登録。
- NN 推定器前の近接重複除去規則と除去率報告を追加。
- 高 K では標本数に対し周囲次元が増大するため、K 飽和曲線の解釈に \(n_N/D_{\rm ambient}\) を併記する規則を追加。
