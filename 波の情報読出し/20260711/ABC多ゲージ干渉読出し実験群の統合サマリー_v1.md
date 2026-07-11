# ABC多ゲージ干渉読出し実験群の統合サマリー v1

**副題:** `p,E,R` 読出しから `R*p`・`R*p^2` 保存写像までの数値検証一覧  
**日付:** 2026-07-11  
**著者:** 木原 範昭  
**位置づけ:** 波の情報読出しシリーズ・数値検証メモ  

---

## 1. 目的

本メモでは、2026-07-11 に実行した ABC 多ゲージ干渉読出し実験群を一覧化する。

本メモは新しい物理仮定を追加するものではない。実行済み数値実験の結果を整理し、次の論文化で参照しやすくするための統合索引である。

---

## 2. 実行スクリプト

統合サマリー生成スクリプトは次である。

```text
run_abc_multigauge_readout_integration_summary_v1.py
```

出力先は次である。

```text
abc_multigauge_readout_integration_summary_result_v1/
```

実行コマンドは次である。

```text
python3 run_abc_multigauge_readout_integration_summary_v1.py
```

---

## 3. 全体判定

全体判定は次である。

```text
experiment_count: 9
all_experiments_valid: true
single_gauge_only_used_any: false
integration_summary_valid: true
```

9本の実験すべてが `valid` であり、単一ゲージのみの読出し判定は用いていない。

---

## 4. 実験一覧

| 実験 | 目的 | count | valid | single gauge only |
|---|---|---:|---|---|
| `single_collision_multigauge_readout` | 単回ABC衝突で `p/E/R` を多ゲージ干渉読出しする |  | `true` | `false` |
| `multi_collision_multigauge_readout` | 対称ABC衝突の反復で `p/E/R` 読出しを維持する | `8` | `true` | `false` |
| `readout_robustness_sweep` | 複数の読出し器構成で `p/E/R` 再構成が安定する | `5` | `true` | `false` |
| `asymmetric_amplitude_diagnostic` | 非対称Rで単純反転が保存を破ることを検出する | `8` | `true` | `false` |
| `generalized_elastic_collision_readout` | 非対称Rで `R*p` と `R*p^2` を保存する一般化写像を読む | `8` | `true` | `false` |
| `generalized_velocity_sweep` | 非単位・非対称位相勾配でも一般化写像が成立する | `9` | `true` | `false` |
| `generalized_multi_collision` | 一般化写像を複数回AB衝突へ反復適用する | `4` | `true` | `false` |
| `generalized_noise_robustness` | ゼロ平均読出しノイズの相殺と共通バイアス検出を確認する | `4` | `true` | `false` |
| `generalized_extreme_R_sweep` | 極端なR比でも一般化写像と読出しが維持されるか調べる | `12` | `true` | `false` |

---

## 5. 実験群の構造

本実験群は、次の順に積み上がっている。

1. ABC 単回衝突で `p,E,R` が多ゲージ干渉読出しできることを確認する。
2. 対称衝突で、反復しても `p` 反転、`E,R` 保存が維持されることを確認する。
3. 読出し器の構成を変えても `p,E,R` が安定に再構成されることを確認する。
4. 非対称 `R` 条件では、単純な `q -> -q` 反転が `R*p` 保存を破ることを診断する。
5. 非対称 `R` 条件に対して、`R*p` と `R*p^2` を保存する一般化衝突写像を構成する。
6. 初期位相勾配を非単位・非対称にしても、一般化写像が成立することを確認する。
7. 一般化写像を複数回衝突へ反復しても、保存読出しが維持されることを確認する。
8. ゼロ平均読出しノイズは多ゲージ平均で相殺され、共通バイアスは検出されることを確認する。
9. `R` 比が大きく非対称な条件でも、一般化写像と多ゲージ読出しが維持されることを確認する。

---

## 6. 本サマリーで主張しないこと

本サマリーは、次を主張しない。

| 主張しないこと | 理由 |
|---|---|
| 標準物理量との完全対応 | 本実験群の `p,E,R` は位相系内の読出し量である |
| 単一ゲージ測定の成立 | すべて多ゲージ干渉読出しを要求している |
| 任意条件での一般化衝突成立 | 検証範囲は各実験の条件内に限る |
| 実在粒子衝突の定量予言 | 対応写像は別途構成する必要がある |

---

## 7. 結論

2026-07-11 の ABC 多ゲージ干渉読出し実験群では、9本の数値検証を実行した。

すべての実験が `valid` となり、単一ゲージのみの判定は一つも用いなかった。

これにより、本検証範囲では、ABC 衝突モデルから

```text
p_read
E_read
R_read
R*p
R*p^2
```

に相当する保存読出しを、多ゲージ干渉により一貫して構成できることが確認された。

---

# 付録A. 出力ファイル

| 種類 | ファイル |
|---|---|
| 実行レポート | [abc_multigauge_readout_integration_summary_report_v1.md](abc_multigauge_readout_integration_summary_result_v1/abc_multigauge_readout_integration_summary_report_v1.md) |
| 結果 JSON | [abc_multigauge_readout_integration_summary_result_v1.json](abc_multigauge_readout_integration_summary_result_v1/abc_multigauge_readout_integration_summary_result_v1.json) |
| CSV | [abc_multigauge_readout_integration_summary_v1.csv](abc_multigauge_readout_integration_summary_result_v1/abc_multigauge_readout_integration_summary_v1.csv) |
