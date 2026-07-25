# exact_lowN_eigenspectrum_v2 修正指示
## N=5 branch追跡・初期床判定・縮退感度確認

## 0. 作業範囲

`exact_lowN_eigenspectrum_v2` のN=5計算について、固有分解本体は原則として維持し、次の後処理を修正する。

1. branch IDの重複割当
2. `initial_floor_flag` の誤った定義
3. 縮退クラスタ併合閾値への依存性
4. crossing代表時刻の表記
5. branch別図表と報告書の再生成

N=40およびN=300にはまだ進まない。

物理解釈、論文6の評価、新方向創発の結論は書かない。

---

# 1. 現行v2で維持してよい部分

次の実装と結果は維持してよい。

- `H = 1j * K` に対する `numpy.linalg.eigh`
- 全固有値の保存
- 正負固有値対の検証
- 固有対残差の保存
- Hermitian固有ベクトルの直交性検証
- 実回転平面の構成
- 瞬時縮退クラスタ射影
- 瞬時クラスタ占有
- 初期親平面と瞬時支配平面の直接比較
- q特異値
- N=5最終スペクトル

修正版でも、最終スペクトル

\[
1,\quad
0.251458,\quad
0.250026,\quad
0.249974,\quad
0.248533
\]

が再現されることを確認する。

固有値計算自体を、旧v1方式の `numpy.linalg.eig(K)` に戻してはならない。

---

# 2. branch追跡の重複割当を禁止する

## 2.1 現行問題

現在の追跡では、各現在クラスタが独立に最大重なりの過去クラスタを選んでいる。

この方法では、複数の現在クラスタが同じ過去branch IDを取得できる。

実際に最終表で、異なる二つのクラスタに同じbranch IDが割り当てられている。

これは一対一branch追跡として不正である。

## 2.2 必須修正

隣接時刻 \(t_a,t_b\) の過去クラスタ射影を

\[
P_i(t_a)
\]

現在クラスタ射影を

\[
P_j(t_b)
\]

とし、重なり行列を

\[
C_{ij}
=
\frac{
\operatorname{Tr}\!\left[P_i(t_a)P_j(t_b)\right]
}{
\sqrt{
\operatorname{rank}P_i(t_a)\,
\operatorname{rank}P_j(t_b)
}
}
\]

または、既存定義と整合する正規化射影重なりで構成する。

一対一branch対応は、重なり総和を最大化する割当問題として解く。

推奨：

```python
from scipy.optimize import linear_sum_assignment

row_ind, col_ind = linear_sum_assignment(-overlap_matrix)
```

各過去branchは、通常追跡では高々一つの現在クラスタへ割り当てる。

各現在クラスタも、高々一つの過去branchを継承する。

## 2.3 最小重なり条件

一対一割当が得られても、重なりが小さすぎる対応はbranch継承とみなさない。

最低重なり閾値は一つに固定せず、少なくとも次を診断する。

```text
0.50
0.70
0.90
0.99
```

主要branch IDには事前に採用した基準を明記する。

閾値未満の割当は、

```text
tracking_status = unmatched
```

とし、新しいbranch IDを付与する。

## 2.4 branch ID一意性検査

各時刻について、次を自動検査する。

```python
assert len(active_branch_ids) == len(set(active_branch_ids))
```

一対一branch表では、同一時刻に同じbranch IDが二回以上現れてはならない。

違反した場合は処理を失敗終了させる。

---

# 3. 分裂・合流はbranch IDと別にlineageで保存する

一対一branch IDだけでは、縮退クラスタの分裂・合流を表現できない。

そのため、通常branch追跡と、分裂・合流lineageを分ける。

## 3.1 一対一branch

時系列上の代表的な連続対応を表す。

- 一つの過去branchから一つの現在branch
- 同一時刻の重複ID禁止

## 3.2 lineage edge

一対一割当とは別に、重なりが一定値以上の全対応を保存する。

出力例：

```text
source_step
target_step
source_cluster_id
target_cluster_id
source_branch_id
target_branch_id
overlap
relation_type
```

`relation_type` は少なくとも次を持つ。

```text
continuation
split_candidate
merge_candidate
birth
death
ambiguous
```

判定規則は明示する。

例：

- 一つのsourceが複数targetへ高重なり  
  → `split_candidate`
- 複数sourceが一つのtargetへ高重なり  
  → `merge_candidate`
- 対応sourceなし  
  → `birth`
- 対応targetなし  
  → `death`

ただし、これらは幾何的追跡ラベルであり、物理解釈を付けない。

---

# 4. `initial_floor_flag` を初期時刻との直接比較で定義する

## 4.1 現行問題

前時刻との対応が見つからなかったことを、

```text
initial_floor_flag = 1
```

とするのは誤りである。

対応失敗は、枝の出現以外にも、

- 縮退クラスタの分裂
- 縮退クラスタの合流
- 急激な方向変化
- サンプリング間隔
- branch追跡失敗
- 閾値選択

で発生する。

## 4.2 初期床の定義

初期時刻 \(t_0\) の全固有値を保存し、各初期固有値について数値床比を

\[
\rho_j(t_0)
=
\frac{
|\sigma_j(t_0)|
}{
\epsilon_{\mathrm{mach}}\|K(t_0)\|_2
}
\]

として計算する。

初期床候補は、複数の基準で併記する。

```text
rho < 1e2
rho < 1e3
rho < 1e4
rho < 1e6
```

単一閾値だけで最終分類しない。

## 4.3 branchの初期空間重なり

現在クラスタ \(P_C(t)\) と、初期床部分空間 \(P_{\mathrm{floor},0}\) の重なりを

\[
O_{\mathrm{floor},C}(t)
=
\frac{
\operatorname{Tr}
\left[
P_{\mathrm{floor},0}P_C(t)
\right]
}{
\operatorname{rank}P_C(t)
}
\]

として保存する。

同様に、初期非床部分空間との重なりも保存する。

```text
overlap_with_initial_floor_space
overlap_with_initial_nonfloor_space
```

`initial_floor_flag` を単一の0/1へ潰す前に、この連続値を主データとして保存する。

## 4.4 推奨ラベル

二値フラグではなく、次を保存する。

```text
initial_origin_status
```

候補：

```text
initial_nonfloor
initial_floor
mixed
undetermined
```

例：

- 初期床空間重なり ≥ 0.99  
  → `initial_floor`
- 初期非床空間重なり ≥ 0.99  
  → `initial_nonfloor`
- 両方に有限成分  
  → `mixed`
- 追跡または縮退のため判定不能  
  → `undetermined`

閾値感度も併記する。

---

# 5. 縮退クラスタ併合閾値の感度解析

## 5.1 現状

現行実装は、

\[
\mathrm{MERGE\_TOL}=10^{-12}
\]

を使用している。

非縮退平面間重なりの最大値が約

\[
9.65\times10^{-13}
\]

であり、採用閾値に近い。

通常時は \(10^{-15}\) 程度でも、単一時刻のスパイクがあるため、閾値依存性を確認する。

## 5.2 必須スイープ

少なくとも次の併合基準でN=5後処理を再実行する。

```text
1e-10
1e-11
1e-12
1e-13
1e-14
```

可能なら相対基準も追加する。

\[
\tau_{\mathrm{merge}}
=
c\,\epsilon_{\mathrm{mach}}\|K\|_2
\]

## 5.3 各閾値で比較する量

各時刻について、少なくとも次を比較する。

- クラスタ数
- 各クラスタ次元
- 支配クラスタID
- 支配クラスタ固有値帯
- クラスタ射影
- 親平面重なり
- \(\delta_C\)
- クラスタ占有 \(E_C\)
- q特異値
- branch数
- branch対応
- lineage edge数
- 非縮退平面間最大重なり
- 状態分解閉鎖誤差

## 5.4 安定性表

代表時刻ごとに、次の表を作る。

```text
step
merge_tolerance
cluster_count
cluster_dimensions
dominant_cluster_dimension
dominant_delta
dominant_occupation
q1
q2
q3
q4
max_intercluster_overlap
closure_error
```

## 5.5 合格条件

主要観測量、

- 全固有値
- 支配クラスタ
- 支配クラスタ占有
- 支配 \(\delta\)
- q特異値
- 初期床由来候補branch

が、合理的な閾値範囲で定性的・定量的に安定していること。

不安定な場合は、単一閾値の結果を採用せず、縮退部分空間全体で報告する。

---

# 6. crossing時刻の表記修正

実際のcrossingは、

\[
t_{\mathrm{cross}}=1167
\]

である。

保存代表時刻が5step刻みのため、

\[
t=1165
\]

を使用している場合、`crossing` と表記してはならない。

次のいずれかに統一する。

```text
nearest_sample_before_crossing_step1165
```

または、

```text
crossing_neighborhood_step1165
```

報告書本文には、

```text
実際のcrossingはstep=1167。
代表表は保存時刻中の直前点step=1165を使用。
```

と明記する。

可能であれば、step=1167を追加で保存・後処理し、真のcrossing表を出す。

軌道を再計算する場合は、既存軌道との一致を検証する。

---

# 7. 修正後の必須出力

## 7.1 branch一対一追跡

```text
raw/N00005/branch_tracking_bijective.csv
```

必須列：

```text
source_step
target_step
source_cluster_id
target_cluster_id
source_branch_id
target_branch_id
overlap
assignment_cost
accepted
tracking_status
```

## 7.2 lineage

```text
raw/N00005/cluster_lineage.csv
```

必須列：

```text
source_step
target_step
source_cluster_id
target_cluster_id
source_branch_id
target_branch_id
overlap
relation_type
```

## 7.3 初期床空間対応

```text
raw/N00005/initial_space_origin.csv
```

必須列：

```text
step
cluster_id
branch_id
sigma_min
sigma_max
sigma_representative
overlap_with_initial_floor_space
overlap_with_initial_nonfloor_space
initial_origin_status
floor_threshold_label
```

## 7.4 縮退感度解析

```text
raw/N00005/merge_tolerance_sensitivity.csv
```

## 7.5 診断

```text
diagnostics/N00005_branch_revision.json
```

最低限：

```text
duplicate_branch_id_count
unmatched_cluster_count
ambiguous_assignment_count
split_candidate_count
merge_candidate_count
max_intercluster_overlap_by_tolerance
closure_error_by_tolerance
crossing_step
representative_crossing_step
```

---

# 8. 修正後の必須図

## 図A：一対一branch固有値推移

\[
\sigma_{\mathrm{branch}}(t)/\sigma_1(t)
\]

同一時刻にbranch重複がないこと。

## 図B：lineage図

各クラスタの分裂・合流候補を、ノードとエッジで表示する。

物理解釈は書かない。

## 図C：初期床空間重なり

各branchについて、

\[
O_{\mathrm{floor}}(t)
\]

と、

\[
O_{\mathrm{nonfloor}}(t)
\]

を表示する。

## 図D：初期床由来候補branchの固有値

初期床空間との重なりが大きいbranchについて、

\[
\sigma_j(t)/\sigma_1(t)
\]

を表示する。

## 図E：縮退閾値感度

横軸を時刻、または代表時刻ごとの閾値として、

- クラスタ数
- 支配 \(\delta\)
- 支配占有
- q3
- q4

を比較する。

## 図F：branch追跡の最小重なり

各時刻の採用branch対応について、対応重なりの最小値を表示する。

---

# 9. 現行図表の扱い

修正完了まで、次を未確定扱いにする。

- branch IDを使用した全図
- branch別 \(\delta\)
- 初期床由来branchの図
- branchのbirth/death
- branchの分裂・合流
- `initial_floor_flag` を使用した集計
- 最終表のbranch列

次は維持してよい。

- 全固有値スペクトル
- 正負対検証
- 固有対残差
- 瞬時クラスタの生データ
- 瞬時クラスタ占有
- 親平面との瞬時比較
- q特異値
- 全体数値診断

---

# 10. 実装上の禁止事項

- 各現在クラスタが独立に最大重なりbranchを取得する
- 同一時刻にbranch IDを重複させる
- 対応失敗を初期床由来と同一視する
- `initial_floor_flag` を前時刻追跡だけで決める
- 縮退閾値を一値だけ試して安定とみなす
- branch分裂・合流を通常branch IDだけで表す
- branch対応が曖昧なのに強制継承する
- N=40へ先に進む
- N=300へ着手する
- 解釈を書く

---

# 11. N=5再検収条件

以下をすべて満たした後にだけN=40へ進む。

1. 同一時刻のbranch ID重複がゼロ
2. Hungarian法等による一対一対応が実装済み
3. 分裂・合流lineageが別データとして保存済み
4. 初期床判定が初期時刻部分空間との直接比較に基づく
5. `initial_floor_flag` の旧定義を廃止
6. 併合閾値 \(10^{-10}\) から \(10^{-14}\) の感度表がある
7. 主要観測量の閾値安定性が確認されている
8. crossingとcrossing近傍の表記が区別されている
9. branch別全図が修正版生データから再生成可能
10. 固有値スペクトルが現行v2と一致
11. 数値診断が倍精度範囲
12. 観測報告に物理解釈がない

完了後は、N=5の修正版一式だけを提出して停止すること。

N=40は人間による検収後に別途指示する。
