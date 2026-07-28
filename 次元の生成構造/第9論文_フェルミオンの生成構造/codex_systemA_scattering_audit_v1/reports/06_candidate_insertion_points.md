# 06 新散乱候補の挿入候補点

## 0. 状態

**[実装候補]** 本書は将来の挿入位置を示すだけであり、既存System A/System B本体への変更は行っていない。

## 1. カーネル座標への復帰

### 候補位置

```text
System A make_case_state:631-659
System A run_case:718-727
```

### 必要な役割

現在の状態は搬送波 \(e^{iqp_0u}\) を含む。パリティ評価前に、

\[
\psi_{\mathrm{kernel}}(u,\eta)
=e^{-iqp_0u}\psi(u,\eta)
\]

を作り、カーネル座標で \(\Pi_B,\Pi_F,c_\pi\) を計算する。

**[実装候補]** 原本の状態生成は変更せず、将来の独立候補モジュールが `a,b,q_A,q_B,coordinate` を受け取って逆搬送波状態を作る。

## 2. 散乱係数と生の経路分解

### 候補位置

```text
System A run_case:718-736
```

現行の、

```text
delta_from_reflection_rate
scattering_coefficients
r*a + t*b
t*a + r*b
```

が一か所に集まっている。

**[実装候補]** 将来はここを直接編集せず、独立関数

```text
scatter_wave_pair(a, b, reflection_parameter, coordinate, carrier_mode, model, ...)
```

の比較呼出点とする。

この関数は少なくとも次を分離して返す必要がある。

```text
path_a_to_a_amplitude = r*a
path_b_to_a_amplitude = t*b
path_b_to_b_amplitude = r*b
path_a_to_b_amplitude = t*a
interference_in_a
interference_in_b
a_raw
b_raw
a_out
b_out
```

## 3. 正規化境界

### 候補位置

```text
System A:734-735
散乱源:774-775,833-834,978-979
```

**[実装候補]** 新候補では次を別々のモデル設定として扱う。

1. 生のユニタリ出力を使用
2. 二チャネル結合ノルムだけを正規化
3. 現行互換のチャネル別正規化
4. 非線形保存則に基づく正規化

どの設定でも、正規化前の経路配列と干渉項を先に保存する。

## 4. パリティ診断の追加位置

### 候補位置

```text
System A row_for_state:662-715
System A MetricContext:304-345
```

**[実装候補]** 将来の独立候補は、現行評価行とは別に次を返す。

```text
C_pi_raw
c_pi
p_B
p_F
cross_parity_matrix
Delta p_B
Delta p_F
eta_[S,P]
```

現行 `harmonic_distribution` は搬送波を含む絶対倍音パワーなので、カーネルパリティ診断の代用にしない。

## 5. 出力スキーマ

### 候補位置

```text
System A run:1340-1375
System A instrumented runner:130-166
```

基底版はスカラー評価だけを保存し、計装版は複素FFT係数を保存する。

**[実装候補]** 新候補の比較出力は、既存出力に上書きせず、

```text
codex_systemA_scattering_audit_v1/stage_C_candidates/
codex_systemA_scattering_audit_v1/data/
codex_systemA_scattering_audit_v1/figures/
```

へ保存する。

## 6. System B境界

### 候補位置

```text
run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py:90-109
```

**[コード上の事実]** System Bは複素スカラー対だけを状態とし、空間波形を持たない。

**[未導出]** System Aで設計したパリティ依存散乱をSystem Bへ適用するには、

1. System Aで得た内部構造依存係数をSystem Bへ縮約する
2. System Bの状態空間を波形へ拡張する
3. System Bを比較対象のまま変更しない

のいずれかを人間が選ぶ必要がある。Stage Aでは選択しない。

## 7. 変更してはならない場所

**[コード上の事実]** 次は同一コードの配下外コピーまたは独立候補であり、現行System Aの実行経路ではない。

```text
時間軸Q軸とフェルミオンの生成構造/.../finite_order_resonance_v1/src/
波の情報読出し/20260715/
波の情報読出し/20260717/
波の情報読出し/20260718/
次元の生成構造/対照実験_公開α論文_v1/
```

これらへ挿入・同期・実行・出力を行わない。

## 8. Stage C前の推奨境界

**[実装候補]**

```text
既存System A原本
  ↓ 読取り専用入力
独立 parity_projection
  ↓
独立 candidate_scattering_models
  ↓
独立 ScatteringResult
  ↓
専用比較データ
```

この構造なら、既存本体を変更せずCandidate 0〜3を同じ入力上で比較できる。

既存本体への組込み位置は以上のとおり特定したが、組込みは未実施である。
