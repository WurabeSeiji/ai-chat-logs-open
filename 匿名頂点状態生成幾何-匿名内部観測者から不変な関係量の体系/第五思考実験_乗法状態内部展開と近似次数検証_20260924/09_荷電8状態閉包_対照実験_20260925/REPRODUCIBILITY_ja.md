# 荷電8状態閉包・対照実験 v1 — 再現性

作成日: 2026-09-25

## 基準データ

`reference_raw.csv` と `reference_metadata.json` は独立に作成済みの荷電二体系 LO 解析基準の固定スナップショットである。生成器から参照側への一方向比較だけを許し、参照誤差を生成則へ戻さない。

基準 raw CSV SHA-256:

`f4f49789e0555903f77cd76b34a3a617d6bf0970a5a077ef0e8ccf38cc111689`

metadata SHA-256:

`ed9d213735fc9f693728b69c3a7073f8f5cb3127f0bfa5f56ac2438e9800dc82`

## 実験本体

```bash
python run_charged_eight_state_closure_v1.py
```

本体では以下を実行する。

1. 静的監査: `micro_eight(s)` が外部物理引数を受けないこと。
2. 4ケース×4000 macro step の external-vs-internal bitwise 対照。
3. 500ランダム状態×11 microstep の同期・恒等保持監査。
4. 6/7状態非閉包反例の生成。
5. 明示11-microstep と折り畳みmacroの12000-step bitwise 照合。
6. `r=50 -> 20` の長時間基準比較。
7. 500/1000/2000/4000/8000 steps-per-orbit の刻み監査。

## 長時間時系列の高速再構成

`reconstruct_sampled_numba.py` は、長時間サンプルCSVを同じ折り畳みRK4式から高速再生成する補助スクリプトである。これは状態閉包の合否判定には使用しない。

再生成時に得られた主要値は元の実行要約と一致した。

```text
samples      = 4884
cross_step   = 488345
last_r       = 19.999894652364606
max |t-tcl|  = 1.1528754839673638e-08
max |H-Hcl|  = 4.2482729081731716e-11
max |r-rref| = 9.069200501699015e-09
max xy error = 3.754312359683735e-06
```

## 図再生成

```bash
python plot_charged_eight_state_closure_v1.py
```

`figures/` に PNG/SVG を生成する。

## 厳密監査の再確認値

独立再監査で以下を再確認した。

```text
4000 macro external vs internal:
  max core difference = 0
  bitwise mismatch     = 0
  identity failure     = 0

12000 macro explicit micro vs collapsed macro:
  max 8-state difference = 0
  bitwise mismatch       = 0
```
