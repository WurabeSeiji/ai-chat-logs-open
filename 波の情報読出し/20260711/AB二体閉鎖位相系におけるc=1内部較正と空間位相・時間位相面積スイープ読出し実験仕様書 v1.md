# AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ読出し実験仕様書 v1

**副題:** 一角度円周位相調和読出しの否定対照を受けた、`s` と `tau_read` の分離、内部 `c=1` 較正、`chi-tau` 面積スイープ、逆冪候補の最小検査  
**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・実験仕様書  
**Version DOI:** pending  
**Concept DOI:** pending  

---

## 1. 目的

本仕様書の目的は、AB 二体閉鎖位相系において、空間位相読出し `chi_read` と時間位相読出し `tau_read` を同時に構成し、内部較正としての

```text
c = 1
```

を用いた場合に、`chi-tau` 面積スイープが読めるかを検査することである。

前段階の一角度円周位相調和読出しでは、`D_AB`, `V_AB`, `f_AB` の調和的読出しと、読出し波による減衰候補が確認対象になった。

しかし、その系では、実験ステップ `s` はあくまで反復番号であり、時間位相読出し `tau_read` そのものではなかった。

したがって、一角度実験で逆比例型または逆二乗型の距離減衰が出ないことは、失敗ではない。

むしろ、

```text
時間方向を読まない位相系では、面積スイープを要求できない
```

という否定対照である。

本仕様書では、その次の最小検査として、

```text
chi_read(s)
tau_read(s)
A_chi_tau(s)
c_internal(s)
```

を同じ AB 閉鎖系から読む。

---

## 2. 前提文書

本仕様書は、次を前提とする。

| 文書 | 本仕様書での役割 |
|---|---|
| [無名等振幅複合波モデル基本公理系v4_純化定義論文.md](../20260710/無名等振幅複合波モデル基本公理系v4_純化定義論文.md) | 公理0、0+、0.5、1、3の正本 |
| [基本公理系 v4.md](../20260710/基本公理系%20v4.md) | 公理系 v4 の現時点での解釈メモ |
| [全正符号ゼロ閉鎖の読出し多重性に関する定義補足.md](全正符号ゼロ閉鎖の読出し多重性に関する定義補足.md) | `Σx_n^2=0` の読出し多重性 |
| [AB二体閉鎖位相系におけるラベルなし二弧相対位相と調和読出しに関する定義補足.md](AB二体閉鎖位相系におけるラベルなし二弧相対位相と調和読出しに関する定義補足.md) | AB 二体系のラベルなし相対位相読出し |
| [AB二体閉鎖位相系における一角度円周位相調和読出し実験仕様書 v1.md](AB二体閉鎖位相系における一角度円周位相調和読出し実験仕様書%20v1.md) | 本仕様書の直接の否定対照 |
| [位置位相自由度の拡張における球殻スイープと逆冪読出しに関する定義補足.md](位置位相自由度の拡張における球殻スイープと逆冪読出しに関する定義補足.md) | 二角度・三角度拡張の後続課題 |

---

## 3. 本仕様書で主張しないこと

本仕様書は、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準相対論の導入 | `c=1` は内部較正条件であり、外部時空定数ではない |
| 標準光速度の導入 | 本仕様書では物理的な光速度を置かない |
| 外部時間座標の仮定 | `tau_read` は閉鎖位相系からの読出しである |
| 実験ステップ `s` を時間と同一視すること | `s` は反復番号であり、`tau_read` ではない |
| 標準重力の導出 | 本仕様書は `chi-tau` 面積スイープの有無を検査する |
| 逆二乗則の要求 | 逆冪候補は、面積スイープが成立した場合のみ検査する |
| 三次元空間の仮定 | ここで加えるのは時間位相読出しであり、外部三次元空間ではない |
| `f_A`, `f_B` の独立導入 | AB 合成補償 `f_AB` の読出し規律を維持する |

---

## 4. 基本方針

### 4.1 `s` と `tau_read` を分ける

実験ステップ `s` は、計算手続きの順序を表す。

一方、時間位相読出しは

```text
tau_read(s)
```

である。

したがって、

```text
tau_read(s) = s
```

を先験的に置かない。

`tau_read(s)` は、AB 閉鎖位相系の内部位相構造から読まれなければならない。

### 4.2 `c=1` は内部較正である

本仕様書でいう `c=1` は、外部時空の光速度ではない。

空間位相読出し `chi_read` と時間位相読出し `tau_read` の交換比を、内部ゲージとして

```text
c_internal = |Delta chi_read| / |Delta tau_read| = 1
```

に較正するという意味である。

したがって、`c=1` は物理名ではなく、読出し単位の内部較正条件である。

### 4.3 面積スイープの有無を先に検査する

逆比例型、逆二乗型、または逆三乗型の読出しを議論する前に、まず

```text
chi-tau 面が本当に形成されているか
```

を検査する。

もし `chi_read` と `tau_read` が独立でなく、同じ一次元曲線の別名にすぎないなら、面積スイープはない。

その場合、逆冪読出しを要求しない。

---

## 5. 実験対象

AB 二体閉鎖対象を

```text
P_AB = A + B
```

とする。

この対象は、定常化後に

```math
Q(P_{AB}) = \sum_n x_n^2 \approx 0
```

を満たすものとして扱う。

ただし、公理系 v4 に従い、外部影響直後から定常化までの準定常・準安定状態の存在は否定しない。

各ステップで、次を分けて記録する。

| 量 | 意味 |
|---|---|
| `Q_raw(s)` | 相互作用直後、再定常化前の閉鎖残差 |
| `Q_closed(s)` | 再定常化後の閉鎖残差 |
| `closure_relaxation(s)` | `Q_raw -> Q_closed` の回復量 |

---

## 6. 読出し量

### 6.1 空間位相読出し

一角度 AB 実験で用いたラベルなし二弧相対位相から、空間位相読出しを

```text
chi_read(s)
```

として読む。

これは外部位置座標ではない。

AB 二体閉鎖系内の円周位相または相対位置位相からの読出しである。

### 6.2 時間位相読出し

AB 閉鎖系内の時間方向候補を

```text
tau_read(s)
```

として読む。

これは実験ステップ `s` ではない。

`tau_read` は、空間位相 `chi_read` と同じ閉鎖対象から、多ゲージ読出しにより再構成される時間位相方向の読出しでなければならない。

### 6.3 内部速度読出し

`chi_read` と `tau_read` から、

```math
v_{\mathrm{read}}(s)
= \frac{\Delta \chi_{\mathrm{read}}(s)}
       {\Delta \tau_{\mathrm{read}}(s)}
```

を読む。

ここで `Delta` は、unwrap された読出し差分を用いる。

### 6.4 c=1 内部較正誤差

内部較正誤差を

```math
\epsilon_c(s) = v_{\mathrm{read}}(s) - 1
```

として記録する。

必要に応じて絶対値読出し

```math
|\epsilon_c(s)|
```

も出力する。

### 6.5 chi-tau 面積スイープ

`chi_read` と `tau_read` が独立な読出し面を形成する場合、離散面積を

```math
A_{\chi\tau}(s)
```

として読む。

実装上は、例えば閉路または時系列窓に対して

```math
A_{\chi\tau}
= \frac{1}{2}\sum_s
\left(
\chi_s \Delta\tau_s
- \tau_s \Delta\chi_s
\right)
```

のような符号付き面積を計算する。

ただし、符号は第一原理名ではない。

必要に応じて、絶対面積または面積包絡を併記する。

### 6.6 有効次元読出し

`chi_read` と `tau_read` の独立性を

```text
d_eff_chi_tau
```

として読む。

作業上は、読出し系列の共分散行列、閉路面積、または独立ゲージ応答の rank により判定する。

---

## 7. 実験モード

本仕様書では、次の段階的モードを置く。

### 7.1 Mode 0: tau disabled control

```text
tau_disabled_control
```

では、時間位相読出しを構成しない。

このモードは、一角度 AB 実験の否定対照である。

期待される結果は、

```text
d_eff_chi_tau = 1
A_chi_tau ≈ 0
```

である。

このモードで逆冪読出しが出ないことは失敗ではない。

### 7.2 Mode 1: tau locked control

```text
tau_locked_control
```

では、`tau_read` を導入するが、`chi_read` と同相または固定比でロックする。

この場合、見かけ上 `tau_read` は存在しても、独立な面は形成されない。

期待される結果は、

```text
d_eff_chi_tau = 1
A_chi_tau ≈ 0
```

である。

このモードは、

```text
tau という名前を付けただけでは時間位相面にならない
```

ことを確認する制御条件である。

### 7.3 Mode 2: tau independent surface

```text
tau_independent_surface
```

では、`chi_read` と `tau_read` を独立な位相方向として読む。

期待される結果は、

```text
d_eff_chi_tau ≈ 2
A_chi_tau != 0
```

である。

このモードで初めて、面積スイープ候補を検査できる。

### 7.4 Mode 3: c1 calibrated surface

```text
c1_calibrated_surface
```

では、`tau_independent_surface` を基礎に、内部較正

```text
|Delta chi_read| / |Delta tau_read| ≈ 1
```

を満たすように読出しゲージを調整する。

期待される結果は、

```text
max |epsilon_c| < tolerance_c
A_chi_tau != 0
```

である。

### 7.5 Mode 4: readout leak sweep

```text
readout_leak_sweep
```

では、一角度 AB 実験と同様に、読出し波の強度を変える。

目的は、`chi-tau` 面積スイープおよび振幅包絡が、読出し波によって減衰または歪むかを検査することである。

| readout_mode | 意味 |
|---|---|
| `readout_off` | 能動読出し波を停止し、最小限の終端読出しのみ行う |
| `readout_weak` | 弱い読出し波で毎ステップ読む |
| `readout_normal` | 標準読出し波で毎ステップ読む |
| `readout_strong` | 強い読出し波で毎ステップ読む |

---

## 8. 初期条件

一角度 AB 実験と同様、まず `pi` 近傍の対称初期偏差を用いる。

| case | 初期二弧相対位相 |
|---|---|
| `near_pi_02deg` | `{178deg, 182deg}` |
| `near_pi_05deg` | `{175deg, 185deg}` |
| `near_pi_10deg` | `{170deg, 190deg}` |
| `near_pi_20deg` | `{160deg, 200deg}` |
| `wide_35deg` | `{145deg, 215deg}` |
| `wide_60deg` | `{120deg, 240deg}` |

時間位相側には、次の作業ゲージを置く。

| tau gauge | 意味 |
|---|---|
| `tau_disabled` | `tau_read` を構成しない |
| `tau_locked_0deg` | `chi_read` と同相にロックする |
| `tau_locked_90deg` | 固定位相差を置くが独立発展させない |
| `tau_independent_slow` | `chi_read` より遅い独立時間位相候補 |
| `tau_independent_c1` | `c_internal≈1` を狙う独立時間位相候補 |
| `tau_independent_fast` | `chi_read` より速い独立時間位相候補 |

ここでいう `slow`, `fast` は外部速度名ではない。

内部読出しゲージの掃引名である。

---

## 9. 実験手順

各ケースで、次を実行する。

1. AB 二体閉鎖初期状態を生成する。
2. `D_AB(0)` を `pi` 近傍の対称二弧位相へ設定する。
3. `Protocol F` と `Protocol B` をそれぞれ選ぶ。
4. `tau_mode` を選ぶ。
5. `readout_mode` を選ぶ。
6. ステップ `s=0...N` について、AB 相互閉鎖の更新を行う。
7. 相互作用直後の `Q_raw(s)` を記録する。
8. 再定常化後の `Q_closed(s)` を記録する。
9. `chi_read(s)` を読む。
10. `tau_read(s)` を読む。ただし `tau_disabled_control` では読まない。
11. `v_read(s)`, `epsilon_c(s)`, `A_chi_tau(s)`, `d_eff_chi_tau(s)` を計算する。
12. `D_AB(s)`, `V_AB(s)`, `rho_AB(s)`, `theta_AB(s)`, `f_AB(s)` を併記する。
13. `Protocol F/B`, `tau_mode`, `readout_mode` 間の差を比較する。

明示的に禁止する更新は次である。

```text
tau_read(s) = s
F = G m_A m_B / r^2
F = k x
x(s+1) = x(s) + v(s) ds
v(s+1) = v(s) + a(s) ds
```

本実験では、これらの標準力学式を実装しない。

---

## 10. 逆冪候補の扱い

本仕様書では、逆冪則を先に要求しない。

ただし、`chi-tau` 面積スイープが成立した場合に限り、候補として次を記録する。

```text
power_candidate_alpha
```

これは、読出し強度候補 `I_read` と有効距離読出し `L_read` の関係を

```math
I_{\mathrm{read}} \propto L_{\mathrm{read}}^{-\alpha}
```

として後処理フィットした値である。

重要なのは、次である。

| 状態 | 解釈 |
|---|---|
| `A_chi_tau≈0` | 逆冪候補を判定しない |
| `A_chi_tau!=0`, `alpha≈0` | 面積はあるが距離減衰なし |
| `A_chi_tau!=0`, `alpha≈1` | 一次元的または線的減衰候補 |
| `A_chi_tau!=0`, `alpha≈2` | 面積スイープ型候補 |
| `A_chi_tau!=0`, `alpha≈3` | 体積スイープ型候補 |

この `alpha` は標準重力指数ではない。

閉鎖位相系内の読出し候補である。

---

## 11. 判定量

最低限、次を JSON と CSV に出力する。

| 量 | 目的 |
|---|---|
| `case_id` | 初期条件識別 |
| `protocol` | `Protocol F` または `Protocol B` |
| `tau_mode` | 時間位相読出し条件 |
| `readout_mode` | 読出し波条件 |
| `step_count` | 実行ステップ数 |
| `chi_read_series` | 空間位相読出し系列 |
| `tau_read_series` | 時間位相読出し系列 |
| `v_read_series` | 内部速度読出し系列 |
| `epsilon_c_series` | `c=1` 内部較正誤差 |
| `A_chi_tau_series` | `chi-tau` 面積スイープ系列 |
| `d_eff_chi_tau_series` | 有効次元読出し系列 |
| `D_AB_series` | ラベルなし二弧相対位相系列 |
| `V_AB_series` | 対称偏差系列 |
| `rho_AB_series` | 合成半径的読出し系列 |
| `theta_AB_series` | 円周位相読出し系列 |
| `f_AB_series` | AB 合成補償読出し系列 |
| `Q_raw_series` | 再定常化前閉鎖残差 |
| `Q_closed_series` | 再定常化後閉鎖残差 |
| `max_Q_closed_abs` | 閉鎖成立精度 |
| `max_epsilon_c_abs` | 内部 `c=1` 較正誤差 |
| `max_A_chi_tau_abs` | 面積スイープ最大値 |
| `rank_chi_tau` | `chi-tau` 読出し rank |
| `power_candidate_alpha` | 逆冪候補指数 |
| `protocol_readout_difference` | `Protocol F/B` の読出し差 |
| `readout_mode_decay_difference` | 読出し波条件による減衰差 |
| `tau_is_step_used` | 必ず `false` |
| `external_c_used` | 必ず `false` |
| `absolute_background_axis_used` | 必ず `false` |
| `f_A_or_f_B_used` | 必ず `false` |

---

## 12. 成功条件

本予備実験は、次を満たすとき `valid` とする。

| 条件 | 判定 |
|---|---|
| `Q_closed(s)` が全ステップで十分小さい | 閉鎖が維持される |
| `tau_disabled_control` で `A_chi_tau≈0` | 一角度否定対照が再現される |
| `tau_locked_control` で `A_chi_tau≈0` | 名前だけの `tau` を排除できる |
| `tau_independent_surface` で `d_eff_chi_tau≈2` | 独立な `chi-tau` 面が形成される |
| `c1_calibrated_surface` で `max_epsilon_c_abs` が小さい | 内部 `c=1` 較正が成立する |
| 逆冪候補を `A_chi_tau` 成立後にだけ判定する | 読出し順序を守る |
| `Protocol F/B` の差を記録する | 表示多重性を隠さない |
| `readout_mode` による減衰差を記録する | 読出し波停止反証テストを継続する |
| `tau_read=s` を使わない | 外部時間を持ち込まない |
| `f_A`, `f_B` を使わない | `f_AB` 作業仮説を守る |

---

## 13. 保留条件

次の場合、本予備実験は `hold` とする。

| 条件 | 理由 |
|---|---|
| `Q_closed(s)` が大きく崩れる | 閉鎖系として読めない |
| `tau_read=s` を使う | 時間位相読出しの検査にならない |
| 外部 `c` を物理定数として入れる | 内部較正ではなくなる |
| `tau_locked_control` で面積が出たと誤判定する | rank 判定が壊れている |
| `A_chi_tau≈0` のまま逆冪指数を主張する | 面積スイープなしに指数を読んでいる |
| `f_A`, `f_B` を独立力として使う | AB 合成補償の規律を破る |
| `readout_off` が実装されない | 読出し波停止反証テストにならない |

---

## 14. 図化方針

本実験では図化が重要である。

最低限、次の SVG と PNG を作成する。

| 図 | 内容 |
|---|---|
| `chi_tau_surface_path` | `chi_read` と `tau_read` の軌跡 |
| `chi_tau_area_sweep` | 面積スイープ `A_chi_tau` の閉路または窓表示 |
| `c1_calibration_error` | `epsilon_c(s)` の時系列 |
| `tau_disabled_vs_independent` | `tau_disabled`, `tau_locked`, `tau_independent` の比較 |
| `power_candidate_slope` | 面積成立後の `alpha` 候補 |
| `readout_leak_chi_tau` | 読出し波強度と面積・減衰の関係 |

図中の軸は読出し補助表示であり、絶対背景座標の存在を示すものではない。

---

## 15. 出力予定

実行スクリプト案:

```text
run_ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_v1.py
```

出力ディレクトリ案:

```text
ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/
```

主要出力案:

| 種類 | ファイル |
|---|---|
| レポート | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_report_v1.md` |
| JSON | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1.json` |
| 時系列 CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv` |
| 条件比較 CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv` |
| `chi-tau` 図 | `ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png` |
| 面積図 | `ab_two_body_c1_internal_calibration_area_sweep_v1.png` |
| `c=1` 較正図 | `ab_two_body_c1_internal_calibration_error_v1.png` |
| 逆冪候補図 | `ab_two_body_c1_internal_calibration_power_candidate_v1.png` |

---

## 16. 次段階

本仕様書の結果が安定した後、次へ進む。

1. 二角度球殻位相 `S^2` への拡張。
2. 三角度超球殻位相 `S^3` への拡張。
3. `chi-tau` 面積スイープと位置位相球殻スイープの比較。
4. 逆冪指数が `rho^{-1}`, `rho^{-2}`, `rho^{-3}` のどれとして読まれるかの検査。
5. ABC 三体系への拡張。

特に、二角度・三角度拡張に進む前に、本仕様書により

```text
時間位相を読んだだけで面積スイープが現れるのか
```

を確認する。

もしここで面積スイープが出るなら、外部三次元空間を仮定せず、`chi-tau` 面として逆冪候補を検査できる。

もしここでも面積スイープが出ないなら、位置位相自由度を増やす `S^2`, `S^3` 拡張が必要になる。

---

## 17. 結論

本仕様書では、AB 二体閉鎖位相系において、`c=1` を外部物理定数ではなく内部較正条件として扱い、空間位相読出し `chi_read` と時間位相読出し `tau_read` から `chi-tau` 面積スイープが構成できるかを検査する条件を定めた。

本実験の要点は、`s` と `tau_read` を同一視しないことである。

実験ステップ `s` は時間ではない。

時間位相は、閉鎖位相系からの読出しとして再構成されなければならない。

また、逆冪候補は、面積スイープが成立した場合にのみ検査する。

これにより、一角度 AB 実験で距離減衰が出なかった事実を失敗としてではなく、`chi-tau` 面を持たない否定対照として位置づける。

本仕様書は、二角度・三角度の位置位相自由度拡張へ進む前に行う、最小の時間位相面積スイープ予備実験である。
