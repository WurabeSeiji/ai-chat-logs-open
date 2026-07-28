# 00 Stage A 監査範囲と原本ハッシュ

## 判定

**[コード上の事実]** Stage A の主対象は、第9論文配下の次の System A コピーである。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/
  対照実験_波束収縮_実行環境_v1/
    20260715/
      run_system_A_localization_exchange_R_sweep_preliminary_v1.py
```

このファイルが直接ロードする散乱源は、同じ第9論文配下の次のファイルだけである。

```text
次元の生成構造/第9論文_フェルミオンの生成構造/
  対照実験_波束収縮_実行環境_v1/
    20260713/
      run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
```

根拠は System A 本体の `SOURCE_PATH` と `load_source_module` である。

```text
run_system_A_localization_exchange_R_sweep_preliminary_v1.py:18-45
```

## 読取り・実行境界

**[モデル定義]** 本監査では次の境界を適用した。

- 第9論文配下の原本: 読取り、静的 AST 解析、SHA-256 比較のみ。原本は変更しない。
- 第9論文配下の独立監査コード: 実行可。出力は専用監査ディレクトリ内だけ。
- 第9論文配下外のプログラム・データ: 読取りと SHA-256 比較だけ。実行、import、修正、出力は禁止。

**[コード上の事実]** 配下外候補は System A の直接・間接 import graph に入っていない。次の3ファイルは独立スクリプトであり、主対象 System A から呼ばれない。

```text
run_minimal_system_B_gray_direct_check_v5.py
phase5_eigenphase_resonance_v2.py
run_two_physical_roots_multiprecision_v1.py
```

配下外の System A 同名コピーは、主対象コピーとバイト単位で同一だった。

```text
SHA-256:
91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec
```

## 実行前ハッシュ

**[数値観測]** 実行前の全11対象について、パス、スコープ、サイズ、更新時刻、行数、SHA-256 を次へ保存した。

```text
manifests/source_manifest_before.json
```

主要4ファイルの実行前 SHA-256 は次のとおり。

| 役割 | SHA-256 |
|---|---|
| System A 本体 | `91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec` |
| 20260713 散乱源 | `f815320f5632ae1b23ccade3a53b01e9110ad770da8407e2b10fa6065ef1695c` |
| System A 計装版 | `8ea54a3a11cbbe98de4dc147ae1e5625561158ab1b07bf6bdaf571137812fb86` |
| 同配下 System B 比較コード | `47cf38edb3a55707137ebf103eabc621443abf4bc2b184517685765c3e510511` |

## 監査方法

**[コード上の事実]** 次の独立監査コードだけを実行した。

```text
stage_A_audit/build_source_manifest.py
stage_A_audit/static_code_audit.py
stage_A_audit/audit_current_scattering.py
stage_A_audit/compare_source_manifests.py
```

これらは監査対象モジュールを import しない。`static_code_audit.py` はソースをテキストとして AST 解析し、`audit_current_scattering.py` は監査で特定した式を独立に転記して数値評価する。

生成ログ:

```text
logs/static_code_inventory.json
logs/current_scattering_diagnostic.json
logs/current_scattering_diagnostic.csv
```

## 最終確認

**[数値観測]** 実行後の全11対象について同じ項目を再取得し、次へ保存した。

```text
manifests/source_manifest_after.json
manifests/source_manifest_comparison.json
```

比較結果は次のとおりである。

```text
source_count_before = 11
source_count_after  = 11
unchanged           = true
changed_path_count  = 0
```

**[コード上の事実]** 第9論文配下の監査対象原本、および配下外の読取り専用対象について、サイズ、更新時刻、行数、SHA-256 の差はなかった。監査対象原本は変更されていない。
