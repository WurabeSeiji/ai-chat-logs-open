# AB二体閉鎖位相系におけるchi-tauネイティブ逆面積拡張スイープ予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 目的

本メモでは、AB 二体 `chi-tau` 面積逆数補償診断を、より広いパラメータ範囲へ拡張した結果を記録する。

前段階では、標準的な `tau_independent_c1` 系列に対して、

```text
native inverse-area compensation has not yet been detected.
```

という結果を得た。

本実験では、初期位相偏差、周波数比、振幅比、位相差、読出しリークを同時に掃引し、

```text
c=1 + rank=2 + A_chi_tau != 0
```

を満たす広い条件のどこかで、ネイティブ読出しが `alpha≈2` を示すかを検査する。

---

## 2. 実行コマンド

```text
python3 run_ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_v1.py
```

出力:

```text
ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_report_v1.md](ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_report_v1.md) |
| JSON | [ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1.json](ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1.json) |
| ケース CSV | [ab_two_body_chi_tau_native_inverse_area_extended_sweep_cases_v1.csv](ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_chi_tau_native_inverse_area_extended_sweep_cases_v1.csv) |
| フィット CSV | [ab_two_body_chi_tau_native_inverse_area_extended_sweep_fits_v1.csv](ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_chi_tau_native_inverse_area_extended_sweep_fits_v1.csv) |
| 参照曲線図 | [ab_two_body_native_inverse_area_extended_reference_curves_v1.png](ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_native_inverse_area_extended_reference_curves_v1.png) |
| alpha scan 図 | [ab_two_body_native_inverse_area_extended_alpha_scan_v1.png](ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_native_inverse_area_extended_alpha_scan_v1.png) |

---

## 4. 結果サマリー

主要結果は次である。

```text
native_inverse_area_extended_sweep_preliminary_valid: true
sweep_case_count: 1323
area_valid_case_count: 1260
c1_surface_like_case_count: 126
fit_count: 4158
native_fit_count: 1056
native_positive2_count: 0
c1_native_positive2_count: 0
constructed_reciprocal_positive2_count: 198
c1_area_min: 0.00024281646266173604
c1_area_max: 3.4426827311724595
```

---

## 5. 読み

### 5.1 拡張範囲でも native alpha=2 は出ない

本スイープでは、次を同時に変えた。

```text
initial deviation
frequency_ratio
amplitude_ratio
phase_shift
readout_leak
```

それでも、

```text
native_positive2_count = 0
c1_native_positive2_count = 0
```

であった。

したがって、現在の AB 二体 `chi-tau` 面積モデルの範囲では、ネイティブ読出しが自然に逆面積型 `alpha≈2` を示す証拠は出ていない。

### 5.2 構成済み対照は alpha=2 を示す

一方、

```text
constructed_reciprocal_positive2_count = 198
```

であった。

これは、`A_chi_tau` が位相偏差の二乗でスケールするため、後処理で `1/A_chi_tau` を作れば逆二乗が出ることを確認している。

しかし、これはネイティブ読出しではない。

### 5.3 作り込みを避けるための境界

今回の結果により、次が明確になった。

```text
chi-tau 面積がある
```

ことと、

```text
閉鎖補償が 1 / chi-tau 面積として自然に読まれる
```

ことは別である。

前者は成立した。

後者は、現時点では未検出である。

---

## 6. 判定

本予備実験は `valid` である。

ただし、これは逆二乗型発見ではない。

むしろ、次の否定対照である。

```text
広い c1/rank/area sweep でも native inverse-area scaling は未検出。
```

これにより、後処理で作った `1/A_chi_tau` を、閉鎖系の自然読出しと混同する危険を避けられる。

---

## 7. 次の課題

次の候補は二つである。

1. 逆面積型を native にするための追加閉鎖補償写像が、第一原理層から必要かを検討する。
2. 逆面積型を保留し、位置位相自由度を増やす `S^2` 型実験へ進む。

現時点では、逆面積型を無理に採用しない。

本メモの結論は次である。

```text
逆二乗は作れる。
しかし、まだ読まれていない。
```
