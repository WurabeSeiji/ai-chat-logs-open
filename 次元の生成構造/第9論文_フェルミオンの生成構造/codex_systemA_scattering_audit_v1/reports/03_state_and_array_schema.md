# 03 状態・配列スキーマ

## 1. System A 空間

**[コード上の事実]**

| 量 | 実型・形状 | 定義箇所 |
|---|---|---|
| `chi` | `float64 ndarray`, `(512,)` | 散乱源:82-85 |
| `eta` | `float64 ndarray`, `(16,)` | 散乱源:82-85 |
| `chi_part` / packet kernel | `float64 ndarray`, `(512,)` | 散乱源:59-68、System A:593-608 |
| `phase_chi` | `complex128 ndarray`, `(512,)` | 散乱源:98、System A:625 |
| `eta_phase` | `complex128 ndarray`, `(16,)` | 散乱源:99、System A:626 |
| `psi` | `complex128 ndarray`, `(512,16)` | 散乱源:100、System A:627 |
| `a`, `b` | `complex128 ndarray`, `(8192,)` | 散乱源:101、System A:628 |
| `t`, `r` | Python `complex` scalar | 散乱源:134-137 |
| `T`, `R` | Python `float` scalar | 散乱源:134-137 |
| `h` | `dict[int,float]` | System A:325-338 |
| 衝突評価行 | `dict[str, scalar|string]` | System A:662-715 |

`a`,`b` は倍音係数列ではない。`chi×eta` 格子上の複素場を平坦化した状態ベクトルである。倍音係数は評価時に `chi` 軸FFTから再計算される。

## 2. 入力波の内部構造

### odd-kernel

\[
\psi(\chi,\eta)
=
K_N(\chi-\chi_0)
e^{iqp_0(\chi-\chi_0)}
e^{im\eta},
\]

を平坦化・単独正規化する。

### explicit-packet

\[
\psi(\chi,\eta)
=
K_{\mathrm{packet}}(\chi)
e^{iqp_0(\chi-\chi_0)}
e^{im\eta},
\]

を平坦化・単独正規化する。

**[コード上の事実]** `HarmonicCase` は倍音番号、重み、位相、波長スケール、平行移動を保持するが、これらは状態生成とメタデータにだけ渡され、散乱係数生成には渡されない。

```text
System A:89-160
System A:631-720
```

## 3. 搬送波

**[コード上の事実]** A側は \(q_A=+1\)、B側は \(q_B=-1\) の搬送位相を持つ。

```text
散乱源:41-47
散乱源:88-101
System A:611-628
```

搬送波逆シフトは実装されていない。FFT評価は搬送波を含む全状態に対して行われる。

```text
System A MetricContext:304-345
```

**[数値観測]** 独立診断で、奇数カーネル

\[
\cos u+\cos3u
\]

は \(c_\pi=-1\) だが、\(e^{iu}\) を掛けた全状態は \(c_\pi=+1\) になった。逆搬送波 \(e^{-iu}\) を掛けると \(c_\pi=-1\) に戻った。

```text
kernel_odd_c_pi = -1
full_state_after_q_plus_1_carrier_c_pi ≈ +1
after_inverse_carrier_c_pi = -1
```

## 4. 散乱中の配列

**[コード上の事実]**

```text
a_raw = r*a + t*b       complex128 ndarray (8192,)
b_raw = t*a + r*b       complex128 ndarray (8192,)
a_out = normalize(a_raw)
b_out = normalize(b_raw)
```

現行基底ランナーは `a_raw`, `b_raw`, `a_out`, `b_out` をファイルへ保存しない。各衝突でスカラー評価行だけを作る。

## 5. 評価配列

| 評価 | 中間表現 | 出力 |
|---|---|---|
| `harmonic_distribution` | `chi`軸FFT `(512,16)` | 絶対倍音パワー辞書 |
| `p_chi` | `chi`軸FFTパワー `(512,)` | 実スカラー |
| `localization` | 全8192成分の確率 | 実スカラー |
| `origin_A/B` | 初期状態との複素内積 | 絶対値二乗 |
| `sim_to_A0/B0` | 倍音パワー辞書 | コサイン類似度 |

`B_to_A_transfer` は `origin_B` ではなく `sim_to_B0` を使用する。

## 6. 計装版

**[コード上の事実]** 計装版だけが、各衝突の複素FFT係数を `.npz` に保存する。

```text
coeffs:
  complex128
  (n_records, 2, n_harmonics, eta_grid_n)
collisions:
  int64 (n_records,)
harmonics:
  int64 (n_harmonics,)
coverage:
  float64 (n_records, 2)
```

```text
run_system_A_localization_exchange_R_sweep_instrumented_v1.py:16-23
同:65-166
```

この計装版も基底版と同じ `r,t` と同じチャネル別正規化を使う。

## 7. System B比較

**[コード上の事実]** 同配下System Bの状態は、

```text
a: scalar complex
b: scalar complex
```

であり、`normalize_pair` は

\[
|a|^2+|b|^2=1
\]

へ結合正規化する。

```text
run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py:90-109
```

System BはSystem Aのケース名称とR列を再利用するが、空間波形や倍音配列を散乱状態として使用しない。
