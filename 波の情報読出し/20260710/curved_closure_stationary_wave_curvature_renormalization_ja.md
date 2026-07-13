# 曲率付き閉鎖定常波による曲率繰り込みと完全反射安定性 v2

**副題:** 全正符号ゼロ閉鎖 `Σx_n^2=0` による曲率相対位相漏れの検出・再選別・交換干渉反射回復の数値構成実験  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・追加論文  
**Version DOI:** 10.5281/zenodo.21332874
**Concept DOI:** 10.5281/zenodo.21304039

フェルミオン様反射の計算式を変更したため、同一条件で再計算した V2 とする。

---

## 要旨

本稿では、本論文系の基本公理系で公理1として採用している全正符号ゼロ閉鎖

```math
Q(x)=\sum_n x_n^2=0
```

を満たす奇数倍音複素波が、曲率付き局所セルにおいてどのように安定化されるかを数値構成した。

この閉鎖条件は、本稿の曲率検証のために新たに導入する仮定ではない。無名性、全正符号ゼロ閉鎖、非自明存在からなる既存の基本公理系において、非自明な複素位相波の存在条件として先に置かれていた公理である。本稿は、その公理1を曲率付き局所セルおよび交換干渉反射へ適用した数値構成実験である。

前稿では、フェルミオン的内部逆相核から交換干渉写像を構成し、外部の `q -> -q` 命令を用いずに完全弾性反射に相当する方向反転出力を得た。本稿では、その交換干渉写像が曲率相対位相漏れを受けたとき、完全反射が破れるか、また閉鎖定常波への再選別によって回復するかを調べる。

曲率効果そのものを実在時空の曲率として定量予言するのではなく、局所セルにおける曲率由来の相対位相漏れを

```math
\delta_{K,m}
```

としてモデル化する。閉鎖ペア

```math
x_m,\qquad ix_m
```

に相対位相漏れが入ると、

```math
x_m^2+\left(e^{i\delta_{K,m}}ix_m\right)^2
=
x_m^2\left(1-e^{i2\delta_{K,m}}\right)
```

となり、一般に閉鎖残差が出る。一方、内部位相再選別 `\beta_{K,m}` により、

```math
\delta_{K,m}+\beta_{K,m}=0
```

が満たされると、閉鎖条件は回復する。

最小実験では、曲率相対位相漏れを受けた過渡状態で閉鎖ペア RMS `1.2319416790092972e-02`、通過漏れ `1.1503183254481797e-01` が発生した。内部位相再選別後の定常状態では、閉鎖ペア RMS は `9.4283259783636047e-19`、通過漏れは `0.0` となった。

広範囲検証では、8種類の曲率相対位相モデルと、`none`, `constant`, `linear`, `affine`, `quadratic`, `cubic`, `full` の7種類の補正自由度を比較した。最大曲率相対位相 `1.2` に対して、無補正では最大通過漏れ `1.6202719613622976e-01` が残り、`full` 補正では閉鎖ペア RMS `7.8949412793793227e-19`、通過漏れ `0.0` へ回復した。

さらに、片側入射の局所交換干渉写像へ残留曲率位相を戻した統合検証では、無補正で最大動的通過漏れ `1.6202719613622971e-01` が発生し、`full` 補正では最大動的通過漏れ `1.6608667985580024e-19` まで低下した。動的散乱は二チャンネル期待式と最大誤差 `5.551115123125783e-17` で一致し、ノルム最大誤差は `4.440892098500626e-16` であった。

以上により、本稿の数値構成範囲では、曲率の影響は存在しないのではなく、過渡状態では閉鎖残差および通過漏れとして現れる。しかし、系が `Σx_n^2=0` を満たす閉鎖定常波へ再選別されると、曲率相対位相漏れは内部位相へ吸収され、完全反射読出しは回復する。

**キーワード:** 全正符号ゼロ閉鎖、曲率付き閉鎖定常波、曲率繰り込み、奇数倍音、交換干渉、完全反射、局所定曲率セル、内部位相再選別

---

## 1. 序論

### 1.1 背景

本シリーズでは、無名性、全正符号ゼロ閉鎖、非自明存在を基本公理として採用している。特に、全正符号ゼロ閉鎖は基本公理系 v2 の公理1であり、本稿のために追加された補助仮定ではない。中心となる閉鎖条件は、

```math
\sum_n x_n^2=0
```

である。

これは共役ノルム

```math
\sum_n |x_n|^2
```

ではない。各成分をそのまま二乗し、全正符号で総和する閉鎖条件である。

この条件は、最小二成分では、

```math
A^2+(iA)^2=0
```

として、外部の負符号係数を使わずに符号反転を内部生成する。したがって、複素位相は装飾ではなく、非自明閉鎖のために必要な位相代数として現れる。

前稿「フェルミオン的逆相核による完全反射写像の干渉構成」では、内部逆相核と交換二経路の干渉から局所交換干渉写像を構成し、外部の `q -> -q` 命令を用いずに完全反射に相当する方向反転読出しを生成した。

本稿では、この交換干渉反射が曲率付き局所セルに置かれたとき、曲率相対位相漏れを受けても安定に回復するかを検討する。

### 1.2 本稿の問い

本稿の問いは次である。

> 曲率付き局所セルで曲率相対位相漏れが生じたとき、全正符号ゼロ閉鎖 `Σx_n^2=0` は、その漏れを検出し、閉鎖定常波への内部位相再選別によって完全反射読出しを回復できるか。

ここで重要なのは、曲率の影響を最初から無視しないことである。むしろ本稿では、曲率相対位相漏れが入ると閉鎖残差および通過漏れが出ることを確認し、その後、閉鎖定常波への再選別によってそれらが消えるかを調べる。

---

## 2. 本稿で主張しないこと

本稿は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 実在時空の曲率効果の定量予言 | 標準的な時空計量、場の方程式、実験単位系を用いない |
| 曲率が常に観測不能であること | 過渡状態では閉鎖残差と通過漏れが発生する |
| 曲率相対位相 `δ_K` の唯一の物理式 | 本稿では検証用の局所相対位相モデルとして置く |
| 標準量子論または一般相対論の導出 | 本稿は内部公理系上の数値構成実験である |
| 連続ハミルトニアンの導出 | 本稿は局所写像および閉鎖定常波の存在条件を扱う |

本稿の主張は、`Σx_n^2=0` を満たす閉鎖定常波が、曲率相対位相漏れを検出し、内部位相再選別によって交換干渉反射読出しを回復する最小数値構成である。

---

## 3. 基本公理と閉鎖零錐

### 3.1 全正符号ゼロ閉鎖

本節では、基本公理系 v2 の公理1を再掲する。本稿はこの公理を出発点として用いるが、本稿の曲率繰り込みを成立させるためにここで都合よく導入するものではない。

公理1は、閉鎖条件を

```math
Q(x)=\sum_{n=1}^N x_n^2=0
```

と置く。

閉鎖状態全体を、

```math
\mathcal N
=
\{x\mid Q(x)=0\}
```

と書く。

`\mathcal N` は、共役ノルムがゼロである集合ではない。全正符号二乗和がゼロになる複素閉鎖零錐である。

### 3.2 閉鎖ペア

最小閉鎖ペアを、

```math
x_m,
\qquad
ix_m
```

とする。

このとき、

```math
x_m^2+(ix_m)^2=0
```

である。

奇数倍音複素波では、`m` を倍音ラベルとし、

```math
h_m=2m+1
```

を用いる。

---

## 4. 曲率付き局所セルと相対位相漏れ

### 4.1 曲率作用の二分類

曲率作用を二つに分ける。

第一は、共通因子型である。

```math
x_m\mapsto g_{K,m}x_m,
\qquad
ix_m\mapsto g_{K,m}ix_m
```

この場合、

```math
(g_{K,m}x_m)^2+(g_{K,m}ix_m)^2=0
```

であり、閉鎖ペアは保存される。

第二は、相対位相型である。

```math
x_m\mapsto x_m,
\qquad
ix_m\mapsto e^{i\delta_{K,m}}ix_m
```

この場合、

```math
x_m^2+(e^{i\delta_{K,m}}ix_m)^2
=
x_m^2\left(1-e^{i2\delta_{K,m}}\right)
```

となり、一般に閉鎖残差が発生する。

### 4.2 閉鎖定常波

曲率付き局所セルで安定に残る波を、閉鎖定常波 `x_K` と呼ぶ。

これは、

```math
Q(x_K)=0
```

かつ、

```math
\mathcal U_K x_K=e^{i\alpha_K}x_K
```

を満たす状態である。

ここで `\mathcal U_K` は曲率付き局所セルを通過した後の局所発展写像であり、`e^{i\alpha_K}` は外部読出しで消える共通位相である。

---

## 5. 曲率相対位相の内部再選別

相対位相漏れ `\delta_{K,m}` に対し、内部位相再選別 `\beta_{K,m}` を導入する。

```math
\delta_{K,m}+\beta_{K,m}=0
```

が成立するなら、

```math
e^{i(\delta_{K,m}+\beta_{K,m})}=1
```

となり、閉鎖ペアは回復する。

したがって、曲率相対位相漏れは消えるのではない。閉鎖定常波として存在可能な内部位相配置へ再選別される。

---

## 6. 交換干渉反射との接続

前稿の交換干渉反射では、有効位相を、

```math
\Delta_{\mathrm{eff}}=\Delta_F
```

とし、純粋逆相核 `\Delta_F=\pi` で完全反射が得られた。

本稿では、曲率相対位相漏れを含めて、

```math
\Delta_{\mathrm{eff},m}
=
\Delta_F+\delta_{K,m}+\beta_{K,m}
```

とする。

完全反射条件は、

```math
\Delta_{\mathrm{eff},m}=\pi
```

である。

曲率相対位相漏れが補正されない場合、

```math
\Delta_{\mathrm{eff},m}=\pi+\delta_{K,m}
```

となり、通過漏れが出る。

一方、

```math
\beta_{K,m}=-\delta_{K,m}
```

が成立すれば、

```math
\Delta_{\mathrm{eff},m}=\pi
```

が回復し、完全反射読出しが戻る。

---

## 7. 数値実験

### 7.1 実験1:最小閉鎖定常波検証

実行スクリプトは次である。

```text
run_curved_closure_stationary_wave_v2.py
```

出力先は次である。

```text
curved_closure_stationary_wave_result_v2/
```

曲率相対位相漏れを、

```math
\delta_{K,m}=\kappa h_m
```

として与えた。

最大 `\kappa=0.012`、最大曲率位相 `1.1879999999999999e+00` に対して、次が得られた。

| 状態 | 閉鎖ペア RMS | 通過漏れ |
|---|---:|---:|
| flat | `0.0000000000000000e+00` | `0.0000000000000000e+00` |
| conformal | `9.4283259783636047e-19` | `0.0000000000000000e+00` |
| transient | `1.2319416790092972e-02` | `1.1503183254481797e-01` |
| stationary | `9.4283259783636047e-19` | `0.0000000000000000e+00` |

緩和列では、最終閉鎖ペア RMS `2.0285753500943141e-09`、最終通過漏れ `2.4915322496409796e-15` となった。

判定はすべて `true` であった。

### 7.2 実験2:広範囲曲率位相モデル掃引

実行スクリプトは次である。

```text
run_curved_closure_stationary_wave_broad_sweep_v2.py
```

出力先は次である。

```text
curved_closure_stationary_wave_broad_sweep_result_v2/
```

曲率相対位相モデルを次の8種類とした。

```text
linear
quadratic_area
cubic_high
quartic_edge
alternating_linear
sinusoidal_loop
mixed_smooth
rippled_random_like
```

内部補正自由度を次の7種類とした。

```text
none
constant
linear
affine
quadratic
cubic
full
```

最大曲率相対位相 `1.2` における補正別集計は次である。

| 補正 | 最大閉鎖ペア RMS | 最大通過漏れ | 最大残留位相 |
|---|---:|---:|---:|
| none | `1.3992789770439718e-02` | `1.6202719613622976e-01` | `1.2000000000000000e+00` |
| constant | `1.2323357129304359e-02` | `1.1632365832184562e-01` | `1.1885985027360557e+00` |
| linear | `1.2317721991993871e-02` | `1.1624067781252551e-01` | `1.2122605013795469e+00` |
| affine | `1.2315221532542749e-02` | `1.1619972192010980e-01` | `1.2235371932078578e+00` |
| quadratic | `1.2300337923101595e-02` | `1.1598359742226144e-01` | `1.2740926312051919e+00` |
| cubic | `1.2278927640523116e-02` | `1.1567739316276086e-01` | `1.3311608455010591e+00` |
| full | `7.8949412793793227e-19` | `0.0000000000000000e+00` | `0.0000000000000000e+00` |

共通因子型曲率は全掃引で閉鎖を保存し、無補正では非自明な曲率漏れが検出された。`full` 補正では全位相モデルで閉鎖と完全反射が回復した。限定補正ではモデル依存の残差が残った。

### 7.3 実験3:片側入射散乱への統合検証

実行スクリプトは次である。

```text
run_curved_closure_scattering_integration_v2.py
```

出力先は次である。

```text
curved_closure_scattering_integration_result_v2/
```

実験2で得た残留曲率相対位相を、片側入射の局所交換干渉写像へ戻した。

補正別集計は次である。

| 補正 | 最大動的通過漏れ | 最大閉鎖ペア RMS | 最大期待式誤差 |
|---|---:|---:|---:|
| none | `1.6202719613622971e-01` | `1.3992789770439718e-02` | `5.5511151231257827e-17` |
| constant | `1.1632365832184560e-01` | `1.2323357129304359e-02` | `1.3877787807814457e-17` |
| linear | `1.1624067781252549e-01` | `1.2317721991993871e-02` | `1.3877787807814457e-17` |
| affine | `1.1619972192010981e-01` | `1.2315221532542749e-02` | `1.3877787807814457e-17` |
| quadratic | `1.1598359742226148e-01` | `1.2300337923101595e-02` | `4.1633363423443370e-17` |
| cubic | `1.1567739316276089e-01` | `1.2278927640523116e-02` | `2.7755575615628914e-17` |
| full | `1.6608667985580024e-19` | `7.8949412793793227e-19` | `1.6608667985580024e-19` |

動的散乱は二チャンネル期待式と最大誤差 `5.551115123125783e-17` で一致し、ノルム最大誤差は `4.440892098500626e-16` であった。

---

## 8. 結果の分類

本稿の結果を、探究型物理学者ロールの分類で整理する。

| 対象 | 分類 | 判定 |
|---|---|---|
| `Σx_n^2=0` が非自明複素閉鎖を要求する | 導出済み帰結 | 保持 |
| 曲率相対位相漏れが閉鎖残差を生む | 数値構成済み帰結 | 保持 |
| 共通因子型曲率作用は閉鎖を破らない | 数値構成済み帰結 | 保持 |
| 内部位相再選別が閉鎖を回復する | 数値構成済み帰結 | 保持 |
| 内部位相再選別が完全反射を回復する | 数値構成済み帰結 | 保持 |
| 実在時空の局所平坦性がこの構造で説明される | 既存理論との接続課題 | 未主張 |

---

## 9. 考察

### 9.1 曲率は消えていない

本稿の結果は、「曲率の影響がない」ことを示すものではない。

むしろ、曲率相対位相漏れを入れた過渡状態では、閉鎖残差と通過漏れが明確に発生する。

したがって、曲率効果は検出可能である。

### 9.2 閉鎖定常波が曲率を内部へ吸収する

重要なのは、曲率相対位相漏れを受けた後、系が閉鎖定常波として再選別される場合である。

このとき、内部位相補正 `\beta_{K,m}` により、

```math
\delta_{K,m}+\beta_{K,m}=0
```

が満たされると、閉鎖残差と通過漏れが消える。

これは、曲率が存在しないという意味ではない。曲率相対位相が、安定に存在できる閉鎖定常波の内部位相構成へ繰り込まれるという意味である。

### 9.3 局所平坦性との関係

標準的な局所平坦近似では、十分小さい領域で曲率効果を無視する。

本稿の構成は異なる。

曲率効果を一度入れ、閉鎖条件が破れることを確認し、その後、閉鎖定常波として再選別された状態で曲率残差が読出しから消えることを確認した。

したがって、本稿の内部公理系では、局所平坦的に見える読出しは、曲率を無視した結果ではなく、曲率込みで閉じる定常干渉モードだけが安定に読まれる結果として解釈できる。

これは標準理論との同一視ではない。対応写像を構成するための作業仮説である。

### 9.4 完全弾性反射との接続

前稿で構成した完全反射写像は、内部逆相核によって、

```math
\Delta_F=\pi
```

を与える交換干渉写像であった。

本稿では、曲率相対位相漏れを入れると、

```math
\Delta_{\mathrm{eff},m}
=
\pi+\delta_{K,m}+\beta_{K,m}
```

となる。

無補正では `\delta_{K,m}` が残り、通過漏れが出る。`full` 補正では `\beta_{K,m}=-\delta_{K,m}` となり、

```math
\Delta_{\mathrm{eff},m}=\pi
```

が回復する。

このとき、片側入射散乱でも `R=1,T=0` が回復した。

---

## 10. 結論

本稿では、全正符号ゼロ閉鎖 `Σx_n^2=0` を満たす奇数倍音複素波が、曲率付き局所セルにおいてどのように安定化されるかを数値構成した。

最小実験では、曲率相対位相漏れを受けた過渡状態で閉鎖ペア RMS `1.2319416790092972e-02`、通過漏れ `1.1503183254481797e-01` が発生した。内部位相再選別後の定常状態では、閉鎖ペア RMS `9.4283259783636047e-19`、通過漏れ `0.0` へ回復した。

広範囲検証では、8種類の曲率相対位相モデルと7種類の補正自由度を比較した。無補正では最大通過漏れ `1.6202719613622976e-01` が残り、`full` 補正では閉鎖ペア RMS `7.8949412793793227e-19`、通過漏れ `0.0` へ回復した。

さらに、片側入射の局所交換干渉写像へ曲率残留位相を戻す統合検証では、無補正で最大動的通過漏れ `1.6202719613622971e-01` が発生し、`full` 補正では最大動的通過漏れ `1.6608667985580024e-19` まで低下した。

以上により、本稿の数値構成範囲では、曲率の影響は無視されていない。曲率相対位相漏れは過渡状態で閉鎖残差および通過漏れとして現れる。しかし、系が `Σx_n^2=0` を満たす閉鎖定常波へ再選別されると、曲率相対位相漏れは内部位相構成へ吸収され、交換干渉による完全反射読出しが回復する。

したがって、`Σx_n^2=0` は単なる保存条件ではなく、非自明複素波の存在条件、曲率相対位相漏れの検出条件、閉鎖定常波への再選別条件、および完全反射読出しの安定条件として働く。

---

# 付録A. 実行済みプログラムと出力

## A.1 最小検証

```text
python3 run_curved_closure_stationary_wave_v2.py
```

出力:

```text
curved_closure_stationary_wave_result_v2/
```

主要ファイル:

| 種類 | ファイル |
|---|---|
| レポート | [curved_closure_stationary_wave_report_v2.md](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_report_v2.md) |
| JSON | [curved_closure_stationary_wave_result_v2.json](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_result_v2.json) |
| sweep CSV | [curved_closure_stationary_wave_sweep_v2.csv](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_sweep_v2.csv) |
| relaxation CSV | [curved_closure_stationary_wave_relaxation_v2.csv](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_relaxation_v2.csv) |
| sweep 図 | [curved_closure_stationary_wave_sweep_v2.png](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_sweep_v2.png) |
| relaxation 図 | [curved_closure_stationary_wave_relaxation_v2.png](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_relaxation_v2.png) |

## A.2 広範囲検証

```text
python3 run_curved_closure_stationary_wave_broad_sweep_v2.py
```

出力:

```text
curved_closure_stationary_wave_broad_sweep_result_v2/
```

主要ファイル:

| 種類 | ファイル |
|---|---|
| レポート | [curved_closure_stationary_wave_broad_report_v2.md](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_report_v2.md) |
| JSON | [curved_closure_stationary_wave_broad_sweep_result_v2.json](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_sweep_result_v2.json) |
| sweep CSV | [curved_closure_stationary_wave_broad_sweep_v2.csv](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_sweep_v2.csv) |
| control CSV | [curved_closure_stationary_wave_broad_conformal_control_v2.csv](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_conformal_control_v2.csv) |
| aggregate 図 | [curved_closure_stationary_wave_broad_aggregate_v2.png](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_aggregate_v2.png) |
| closure heatmap | [curved_closure_stationary_wave_broad_closure_heatmap_v2.png](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_closure_heatmap_v2.png) |
| leakage heatmap | [curved_closure_stationary_wave_broad_leakage_heatmap_v2.png](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_leakage_heatmap_v2.png) |

## A.3 片側入射散乱統合検証

```text
python3 run_curved_closure_scattering_integration_v2.py
```

出力:

```text
curved_closure_scattering_integration_result_v2/
```

主要ファイル:

| 種類 | ファイル |
|---|---|
| レポート | [curved_closure_scattering_integration_report_v2.md](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_report_v2.md) |
| JSON | [curved_closure_scattering_integration_result_v2.json](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_result_v2.json) |
| CSV | [curved_closure_scattering_integration_v2.csv](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_v2.csv) |
| 図 | [curved_closure_scattering_integration_v2.png](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_v2.png) |

---

# 参考文献

## 自己引用

1. 木原範昭「無名等振幅複合波モデル 基本公理系 v2」2026-07-10.
2. 木原範昭「論文0:正曲率定曲率空間における測地的単位セルの歪み — 一辺・角・面積・体積の厳密評価」Version DOI: `10.5281/zenodo.21303433`, Concept DOI: `10.5281/zenodo.20680269`, 2026.
3. 木原範昭「フェルミオン的逆相核による完全反射写像の干渉構成」Version DOI: `10.5281/zenodo.21332867`, Concept DOI: `10.5281/zenodo.21295479`, 2026.

## 外部参考文献

4. H. S. M. Coxeter, *Regular Polytopes*, 3rd ed., Dover, 1973.
5. S. Pancharatnam, “Generalized theory of interference, and its applications,” *Proceedings of the Indian Academy of Sciences A*, 44, 247–262, 1956.
6. M. V. Berry, “Quantal phase factors accompanying adiabatic changes,” *Proceedings of the Royal Society of London A*, 392, 45–57, 1984. DOI: `10.1098/rspa.1984.0023`.
