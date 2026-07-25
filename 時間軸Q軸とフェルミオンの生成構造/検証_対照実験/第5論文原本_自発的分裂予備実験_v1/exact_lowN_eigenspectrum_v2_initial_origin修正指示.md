# exact_lowN_eigenspectrum_v2 初期床空間重なり 修正指示
## N=5 最終後処理修正

## 0. 作業範囲

`exact_lowN_eigenspectrum_v2` のN=5計算について、以下だけを修正する。

1. 初期床空間との重なり計算
2. `initial_origin_status` の再判定
3. FigC・FigDおよび関連集計の再生成
4. 旧方式の成果物の隔離または廃止表示
5. 報告書・READMEの更新

固有分解本体、固有値、branch一対一追跡、lineage、crossing処理、縮退感度解析は原則として変更しない。

N=40およびN=300にはまだ進まない。

物理解釈、論文評価、新方向創発の結論は書かない。

---

# 1. 現行問題

現在の `overlap_space` は、クラスタ基底 \(B\) の先頭2列だけを使っている。

現行実装の概念形：

```python
return float(
    np.sum(
        np.abs(Ucols.conj().T @ (B[:, 0] + 0j)) ** 2
        + np.abs(Ucols.conj().T @ (B[:, 1] + 0j)) ** 2
    ) / B.shape[1]
)
```

この実装では、クラスタ次元が4以上の場合でも、分子は2方向分しか計算していない。

一方、分母はクラスタ全次元

```python
B.shape[1]
```

である。

したがって、多次元クラスタでは重なり値が過小評価される。

実データでは、

\[
O_{\mathrm{floor}}
+
O_{\mathrm{nonfloor}}
=
0.5
\]

となる行が存在している。

初期床空間と初期非床空間が全空間を直交分解しているなら、本来は

\[
O_{\mathrm{floor}}
+
O_{\mathrm{nonfloor}}
=
1
\]

でなければならない。

現在の以下の結果は未確定とする。

- `initial_origin_status`
- `initial_floor` 件数
- `mixed` 件数
- `undetermined` 件数
- FigC
- FigD
- 初期床由来branch判定
- 初期床由来branchの集計
- 初期床由来branchに関する報告書記述

---

# 2. 正しい初期部分空間重なり

初期部分空間基底を \(U\)、現在クラスタ基底を \(B\) とする。

両者は列方向に正規直交基底を持つものとする。

重なりは、

\[
O(U,B)
=
\frac{
\|U^\dagger B\|_F^2
}{
\dim B
}
\]

と定義する。

これは、

\[
O(U,B)
=
\frac{
\operatorname{Tr}
\left(
\Pi_U\Pi_B
\right)
}{
\operatorname{rank}\Pi_B
}
\]

と等価である。

ここで、

\[
\Pi_U=UU^\dagger,
\qquad
\Pi_B=BB^\dagger
\]

である。

---

# 3. 必須コード修正

`overlap_space` を、クラスタ基底の全列を使う実装へ変更する。

推奨実装：

```python
def overlap_space(Ucols, B):
    # Ucols: 初期部分空間の正規直交基底, shape=(M, dU)
    # B: 現在クラスタの正規直交基底, shape=(M, dB)
    if B.shape[1] == 0:
        return 0.0

    if Ucols.shape[1] == 0:
        return 0.0

    overlap_matrix = Ucols.conj().T @ B.astype(complex)

    value = (
        np.linalg.norm(overlap_matrix, ord="fro") ** 2
        / B.shape[1]
    )

    return float(np.real_if_close(value))
```

実基底同士なら、次でもよい。

```python
def overlap_space_real(Ucols, B):
    if B.shape[1] == 0:
        return 0.0

    if Ucols.shape[1] == 0:
        return 0.0

    overlap_matrix = Ucols.T @ B

    return float(
        np.linalg.norm(overlap_matrix, ord="fro") ** 2
        / B.shape[1]
    )
```

先頭2列だけを使ってはならない。

---

# 4. 初期床空間と初期非床空間

初期時刻 \(t_0\) において、全空間を次の二つへ分ける。

\[
\mathcal H
=
\mathcal H_{\mathrm{floor},0}
\oplus
\mathcal H_{\mathrm{nonfloor},0}
\]

それぞれの正規直交基底を、

\[
U_{\mathrm{floor},0}
\]

\[
U_{\mathrm{nonfloor},0}
\]

とする。

各現在クラスタ \(B_C(t)\) について、

\[
O_{\mathrm{floor},C}(t)
=
\frac{
\|U_{\mathrm{floor},0}^\dagger B_C(t)\|_F^2
}{
\dim B_C(t)
}
\]

\[
O_{\mathrm{nonfloor},C}(t)
=
\frac{
\|U_{\mathrm{nonfloor},0}^\dagger B_C(t)\|_F^2
}{
\dim B_C(t)
}
\]

を計算する。

全時刻・全クラスタで、

\[
O_{\mathrm{floor},C}
+
O_{\mathrm{nonfloor},C}
\approx 1
\]

を確認する。

---

# 5. 必須診断

各クラスタについて、

```text
overlap_with_initial_floor_space
overlap_with_initial_nonfloor_space
origin_overlap_sum
origin_overlap_closure_error
```

を保存する。

定義：

\[
S_C
=
O_{\mathrm{floor},C}
+
O_{\mathrm{nonfloor},C}
\]

\[
\epsilon_{\mathrm{origin},C}
=
|S_C-1|
\]

全時刻最大値を診断JSONへ保存する。

```text
max_origin_overlap_closure_error
median_origin_overlap_closure_error
count_origin_overlap_error_gt_1e-12
count_origin_overlap_error_gt_1e-10
```

合格条件：

\[
\max_C
\left|
O_{\mathrm{floor},C}
+
O_{\mathrm{nonfloor},C}
-1
\right|
<
10^{-12}
\]

倍精度誤差の影響で厳密に満たせない場合でも、理由と最大値を報告する。

---

# 6. `initial_origin_status` の再判定

重なり再計算後、旧 `initial_origin_status` を全面的に再生成する。

推奨規則：

```text
initial_floor:
    overlap_with_initial_floor_space >= 0.99

initial_nonfloor:
    overlap_with_initial_nonfloor_space >= 0.99

mixed:
    どちらも0.01を超え、
    どちらも0.99未満

undetermined:
    数値診断失敗、
    origin_overlap_closure_errorが許容値超過、
    またはその他の明示的理由
```

境界値は設定ファイルまたは定数として明記する。

少なくとも次の閾値感度も集計する。

```text
0.90
0.95
0.99
0.999
```

単一閾値の件数だけで結果を確定しない。

---

# 7. 必須再生成ファイル

## 7.1 初期空間対応データ

再生成：

```text
raw/N00005/initial_space_origin.csv
```

必須列：

```text
step
time
cluster_id
branch_id
cluster_dimension
sigma_min
sigma_max
sigma_representative
overlap_with_initial_floor_space
overlap_with_initial_nonfloor_space
origin_overlap_sum
origin_overlap_closure_error
initial_origin_status
origin_threshold
floor_threshold_label
```

## 7.2 診断

再生成または追加：

```text
diagnostics/N00005_initial_origin_revision.json
```

必須項目：

```text
max_origin_overlap_closure_error
median_origin_overlap_closure_error
count_origin_overlap_error_gt_1e-12
count_origin_overlap_error_gt_1e-10
initial_floor_count
initial_nonfloor_count
mixed_count
undetermined_count
status_counts_by_origin_threshold
```

## 7.3 集計表

```text
tables/N00005/initial_origin_summary.csv
```

代表時刻ごとに、

```text
step
label
cluster_count
initial_floor_count
initial_nonfloor_count
mixed_count
undetermined_count
max_origin_overlap_closure_error
```

を出力する。

---

# 8. 必須再生成図

## FigC：初期部分空間重なり

各branchまたは各クラスタについて、

\[
O_{\mathrm{floor}}(t)
\]

\[
O_{\mathrm{nonfloor}}(t)
\]

を表示する。

同一対象について、

\[
O_{\mathrm{floor}}+O_{\mathrm{nonfloor}}
\]

も診断表示する。

可能なら別図として、

```text
figC_origin_overlap_closure.png
```

を追加する。

## FigD：初期床由来候補branch

修正後の

\[
O_{\mathrm{floor}}(t)
\]

が指定閾値以上のbranchだけを抽出し、

\[
\sigma_j(t)/\sigma_1(t)
\]

を表示する。

旧 `initial_floor_flag` を使って抽出してはならない。

図中またはキャプションに、使用したorigin閾値を明記する。

---

# 9. 旧方式成果物の整理

同一ZIP内に旧方式と新方式の成果物を混在させない。

少なくとも以下を確認する。

- 旧 `fig09_floor_branches.png`
- 旧 `fig10_delta_branches.png`
- 旧 `cluster_tracking.csv`
- 旧 `initial_floor_flag` を含む表
- 旧初期床由来集計
- 旧報告書の件数

対応方法は次のいずれか。

## 推奨A：削除

再生成可能であり、旧方式が明確に誤っているなら削除する。

## 推奨B：隔離

```text
deprecated/
```

または、

```text
obsolete_v2_before_origin_fix/
```

へ移動する。

## 必須表示

READMEに次を明記する。

```text
旧 initial_floor_flag は、
前時刻branch対応失敗を初期床由来とみなしていたため廃止。

旧 overlap_space は、
多次元クラスタの先頭2列しか評価していなかったため廃止。

旧方式の図表は解析に使用禁止。
```

旧ファイルを残す場合、ファイル名にも、

```text
_obsolete
```

を付ける。

---

# 10. 報告書修正

`observation_report_branch_revision.md` および関連報告書で、次を更新する。

- 初期床空間重なりの定義
- 全列Frobeniusノルムを使用したこと
- 重なり閉鎖誤差
- `initial_origin_status` の新件数
- 閾値感度
- FigC・FigDの差し替え
- 旧方式結果の廃止

以下は書かない。

- 初期床由来branchの物理解釈
- フェルミオン生成の結論
- 新方向創発の肯定・否定
- N=40への外挿
- N=300への外挿
- 論文6の評価

観測値と数値診断だけを報告する。

---

# 11. 再計算範囲

原則として、保存済みの以下を再利用してよい。

- 固有値
- 固有ベクトル
- 瞬時クラスタ基底
- branch ID
- lineage
- crossing情報

再計算が必要なのは主に、

- 初期床空間基底
- 初期非床空間基底
- 各クラスタとの全列重なり
- origin status
- 関連CSV
- FigC
- FigD
- 集計表
- 報告書

である。

固有分解本体を再実行する必要はない。

---

# 12. 最終検収条件

以下をすべて満たした後にだけN=40へ進む。

1. `overlap_space` がクラスタ基底の全列を使用
2. 初期床空間と初期非床空間の和が全空間を構成
3. 全時刻・全クラスタで
   \[
   |O_{\mathrm{floor}}+O_{\mathrm{nonfloor}}-1|<10^{-12}
   \]
4. `initial_origin_status` を全面再生成
5. 旧 `initial_floor_flag` を解析から排除
6. FigCを再生成
7. FigDを再生成
8. origin閾値感度を保存
9. 旧方式成果物を削除または隔離
10. READMEに旧方式使用禁止を明記
11. 報告書の件数・表・図を修正版へ更新
12. 固有値スペクトルとbranch追跡結果は変更されていない
13. N=40には未着手
14. N=300には未着手
15. 物理解釈を書いていない

完了後は、N=5修正版一式のみを提出して停止すること。
