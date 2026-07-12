# AB二体閉鎖位相系におけるchi-tau面積逆数補償診断予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 目的

本メモでは、AB 二体 `c=1` 内部較正 `chi-tau` 面積スイープ予備実験の出力を用い、

```text
A_chi_tau が成立した後、
閉鎖補償側に 1 / A_chi_tau 型の読出しが自然に残るか
```

を診断した結果を記録する。

重要なのは、`1/A_chi_tau` を作ることではない。

`1/A_chi_tau` は、逆二乗型が出ることを確認するための構成済み対照であり、ネイティブな読出し発見としては数えない。

本診断では、次を明確に分ける。

| 種類 | 扱い |
|---|---|
| native readout | 既存実験出力に自然に含まれる読出し量 |
| constructed reciprocal control | `1/A_chi_tau` のように後処理で明示的に作った量 |
| derived ratio | `f_AB/A_chi_tau` など、既存読出しを面積で割った量 |

---

## 2. 実行コマンド

```text
python3 run_ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_v1.py
```

入力:

```text
ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/
```

出力:

```text
ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_report_v1.md](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_report_v1.md) |
| JSON | [ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1.json](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1.json) |
| ケース診断 CSV | [ab_two_body_chi_tau_inverse_area_compensation_diagnostic_cases_v1.csv](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_compensation_diagnostic_cases_v1.csv) |
| フィット CSV | [ab_two_body_chi_tau_inverse_area_compensation_diagnostic_fits_v1.csv](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_compensation_diagnostic_fits_v1.csv) |
| 逆面積対照図 | [ab_two_body_chi_tau_inverse_area_constructed_control_v1.png](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_constructed_control_v1.png) |
| native 候補図 | [ab_two_body_chi_tau_inverse_area_native_candidates_v1.png](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_native_candidates_v1.png) |
| alpha 比較図 | [ab_two_body_chi_tau_inverse_area_alpha_comparison_v1.png](ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_alpha_comparison_v1.png) |

---

## 4. 結果サマリー

主要結果は次である。

```text
inverse_area_compensation_diagnostic_preliminary_valid: true
diagnostic_case_count: 288
area_valid_case_count: 144
fit_count: 576
native_fit_count: 132
derived_fit_count: 84
native_positive2_count: 0
constructed_reciprocal_positive2_count: 2
c1_readout_off_area_min: 0.0038252030346360663
c1_readout_off_area_max: 3.4426827311724595
```

---

## 5. 読み

### 5.1 1/A は alpha=2 になる

`A_chi_tau` は初期位相偏差の二乗でスケールする。

したがって、後処理で

```text
1 / A_chi_tau
```

を作れば、当然

```text
alpha ≈ +2
```

が得られる。

本診断では、この結果を構成済み対照として扱う。

これは逆二乗候補の数学的確認ではあるが、閉鎖位相系が自然にそれを読んだ証拠ではない。

### 5.2 native readout には alpha=2 が出ていない

ネイティブ候補として、次を検査した。

```text
native_max_f_AB
native_max_Q_raw
native_max_closure_relaxation
native_max_envelope_V_AB
native_max_V_AB
native_max_epsilon_c
```

結果は、

```text
native_positive2_count = 0
```

であった。

つまり、今回の出力範囲では、閉鎖位相系の既存読出し量が自然に `1/A_chi_tau` 型の逆二乗読出しを示したとは言えない。

### 5.3 derived ratio は参考値にとどめる

`f_AB/A_chi_tau` などの derived ratio は、後処理で面積除算を行っている。

したがって、これも native readout とは区別する。

今回の診断では、derived ratio を比較対象として記録したが、発見主張には用いない。

---

## 6. 判定

本診断は `valid` である。

ただし、これは逆二乗型が発見されたという意味ではない。

むしろ、厳密な判定は次である。

```text
chi-tau area exists.
1 / chi-tau area has inverse-square scaling by construction.
native inverse-area compensation has not yet been detected.
```

したがって、本診断は、逆二乗候補を無理に作り込まないための重要な否定対照である。

---

## 7. 次の課題

次に進む場合、少なくとも次のどちらかが必要である。

1. `1/A_chi_tau` を native readout として生む閉鎖補償写像を、第一原理層から別途定義し、その後に検証する。
2. 現段階では逆面積読出しを保留し、位置位相自由度を増やす `S^2` 型または `S^3` 型実験へ進む。

現時点では、どちらを選ぶにしても、今回の結果を境界条件にする。

```text
後処理で作った 1/A を、native readout と混同しない。
```

これが本診断の最重要結論である。
