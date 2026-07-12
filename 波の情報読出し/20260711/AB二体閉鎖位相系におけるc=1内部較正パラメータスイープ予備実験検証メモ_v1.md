# AB二体閉鎖位相系におけるc=1内部較正パラメータスイープ予備実験検証メモ v1

**日付:** 2026-07-12  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・予備実験検証メモ  

---

## 1. 目的

本メモでは、AB 二体 `c=1` 内部較正 `chi-tau` 面積スイープ予備実験に続き、時間位相候補 `tau_read` の

```text
frequency_ratio
amplitude_ratio
phase_shift
readout_leak
```

を掃引した結果を記録する。

目的は、

```text
c=1 較正だけで chi-tau 面が成立するのか
```

を検査することである。

---

## 2. 実行コマンド

```text
python3 run_ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_v1.py
```

出力ディレクトリ:

```text
ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/
```

---

## 3. 出力ファイル

| 種類 | ファイル |
|---|---|
| レポート | [ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_report_v1.md](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_report_v1.md) |
| JSON | [ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1.json](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1.json) |
| ケース要約 CSV | [ab_two_body_c1_internal_calibration_parameter_sweep_case_summary_v1.csv](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_parameter_sweep_case_summary_v1.csv) |
| `c` 誤差 heatmap | [ab_two_body_c1_parameter_sweep_c_error_heatmap_v1.png](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_c_error_heatmap_v1.png) |
| 面積 heatmap | [ab_two_body_c1_parameter_sweep_area_heatmap_v1.png](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_area_heatmap_v1.png) |
| 位相差応答図 | [ab_two_body_c1_parameter_sweep_phase_response_v1.png](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_phase_response_v1.png) |
| 読出しリーク図 | [ab_two_body_c1_parameter_sweep_readout_leak_v1.png](ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_readout_leak_v1.png) |

---

## 4. 結果サマリー

主要結果は次である。

```text
c1_internal_calibration_parameter_sweep_preliminary_valid: true
sweep_case_count: 756
readout_off_case_count: 252
rank2_readout_off_count: 246
c1_surface_like_readout_off_count: 13
c1_locked_like_readout_off_count: 8
min_c_error_readout_off: 1.1102230246251565e-15
max_area_readout_off: 0.19685536479742288
tau_is_step_used_any: false
external_c_used_any: false
f_A_or_f_B_used_any: false
```

---

## 5. 読み

### 5.1 c=1 だけでは足りない

本スイープでは、`c=1` に近い RMS 交換比を持つ設定が複数現れた。

しかし、その中には

```text
locked_like
```

すなわち、`chi-tau` 面積を作らない設定も含まれる。

具体的には、

```text
c1_surface_like_readout_off_count = 13
c1_locked_like_readout_off_count = 8
```

であった。

したがって、内部 `c=1` 較正だけでは、時間位相面が成立したとは言えない。

### 5.2 必要条件

本スイープにより、次の三条件を同時に要求する必要が明確になった。

```text
c=1 calibration
rank_chi_tau = 2
A_chi_tau != 0
```

これは、`tau` という名前だけを付けた一次元表示、または `chi` と固定関係にある `tau` を排除するために必要である。

### 5.3 位相差応答

`frequency_ratio=1`, `amplitude_ratio=1` でも、`phase_shift` により面積は変化した。

`c` 較正誤差はほぼ数値丸め範囲に留まる一方、`A_chi_tau` は大きく変化し、`phase_shift=90deg` 側で面積が消える。

したがって、`c=1` と面積形成は同じ条件ではない。

### 5.4 読出しリーク

`readout_leak` を増やすと、

```text
max |epsilon_c|
|decay_rate_envelope|
```

がともに増加した。

これは、前段階と同様に、読出し波が閉鎖位相系の較正と包絡へ影響する候補である。

---

## 6. 判定

本予備実験は `valid` である。

成立したことは次である。

1. `c=1` 較正候補は複数存在する。
2. `c=1` 較正候補の一部は `chi-tau` 面積を作らない。
3. `c=1`, `rank=2`, `A!=0` の三条件を同時に置く必要がある。
4. 読出しリークは較正誤差と包絡減衰を増やす。

成立していないことは次である。

1. `c=1` だけで時間位相面が成立すること。
2. `chi-tau` 面積から逆冪減衰が自動的に出ること。
3. 標準重力対応。

---

## 7. 次の課題

次の段階では、次のどちらへ進むかを判定する必要がある。

1. `chi-tau` 面に伝搬写像を追加し、面積スイープが距離減衰へ変換されるかを検査する。
2. 位置位相自由度を二つに増やし、`S^2` 型の球殻スイープへ進む。

現時点の最重要結論は、

```text
c=1 は必要条件候補だが、十分条件ではない。
```

である。
