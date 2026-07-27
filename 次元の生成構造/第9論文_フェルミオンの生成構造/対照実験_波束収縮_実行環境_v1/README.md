# 対照実験_波束収縮_実行環境_v1

第9論文「フェルミオンの生成構造」の調査指示（`../波束収縮と交差項の関係性調査.md`、`../フェルミオン的な相互作用の探求指針.md`）に基づく再解析・対照実験用の隔離実行環境。

## 目的

原本のプログラム・実験データを一切変更せずに、波束収縮実験（系統A）・灰色猫準安定実験（系統B）の再実行・対照実験・内部量計測を行う。

## 構成方針：ミラー構造による無修正コピー

原本スクリプトは相対パス（`BASE_DIR.parent / "20260713" / ...`）で相互参照するため、日付フォルダ階層ごと複製した。これにより **コードは1文字も変更していない**（SHA256で原本と同一であることを検証済み）。

```
対照実験_波束収縮_実行環境_v1/
├── 20260711/ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v2/
│   └── ...result_v2.json                    # エンジンが読む加速度V2基底（データ、コピー）
├── 20260713/
│   └── run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py
│                                            # 物理エンジン正本のコピー
└── 20260715/
    ├── run_system_A_localization_exchange_R_sweep_preliminary_v1.py   # 系統Aランナー
    └── run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py    # 系統Bランナー
```

委譲チェーン: 系統B → 系統A → 20260713エンジン → 20260711基底JSON（全て本フォルダ内で閉じる）。

## 原本（正本）の所在と同一性検証

コピー元（リポジトリルート相対）と SHA256（2026-07-28 検証、コピーと完全一致）:

| ファイル | 原本パス | SHA256 |
|---|---|---|
| エンジン | `波の情報読出し/20260713/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py` | `f815320f5632ae1b23ccade3a53b01e9110ad770da8407e2b10fa6065ef1695c` |
| 系統Aランナー | `波の情報読出し/20260715/run_system_A_localization_exchange_R_sweep_preliminary_v1.py` | `91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec` |
| 系統Bランナー | `波の情報読出し/20260715/run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py` | `47cf38edb3a55707137ebf103eabc621443abf4bc2b184517685765c3e510511` |
| 加速度V2基底 | `波の情報読出し/20260711/ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v2/...result_v2.json` | `ab3b2f4d883e4cb5621188b5384128d7b0f8f0480b3c2ab6689c174eb6db097b` |

エンジンの最終git commit: `fccc6561`（Add localization exchange waveform summary）。

## 出力の隔離

ランナーの `DEFAULT_OUT_DIR` は `Path(__file__).parent` 起点のため、出力は自動的に本フォルダ内
`20260715/system_A_localization_exchange_R_sweep_result_v1/`（系統B同様）に書かれる。
原本の結果ディレクトリに書き込まれることはない。`--run-id <名前>` でサブフォルダ分離、`--output-dir` で任意指定も可能。

## 動作検証（parity check、2026-07-28 実施）

コピー環境で原本と同一条件を再実行し、原本CSVと突合した:

- 条件: `--pairs 1:1 --r-values 0.0`（衝突256回デフォルト、odd_kernel|A=1|B=1）
- 突合先: 原本 `odd_kernel_N_1_2_3_5_15_63_collision_terrain_v1.csv` の同一 case_id・R 行
- 結果: **物理量 L_A, L_B, N_eff_A, N_eff_B, B_to_A_transfer, A_to_B_transfer は全257衝突で文字列レベル完全一致**（環境は完全決定論、乱数不使用）
- 唯一の差異: `gap_terrain_score` 系の導出スコア。これはスイープ集合全体の最大値を正規化分母に使うため、掃引点構成が異なると値が変わる（物理量の差ではない）。原本と同じ土俵で比較する場合は `--fixed-l-norm` / `--fixed-n-norm` で分母を固定すること（原本の `*_fixed_global_norm` 実験と同じ扱い）。
- 検証出力: `20260715/system_A_localization_exchange_R_sweep_result_v1/parity_check_R0_v1/`（スモークテスト `smoke_test_v1/` も保存）

## 対照テストスイート（parity_suite_v1、2026-07-28 実施）

原本3実験（odd_kernel本体・femtofocus R137極細・phase10位相対照）のパラメータを
原本JSONから機械復元して再実行し、**原本に存在する全データのセル単位完全一致を確認済み**。
femtofocus はバイト単位で同一。詳細と再現データ一式は `parity_suite_v1/` 参照
（唯一の差異は後年のスクリプト拡張による追加メタデータ列で、データ値の差はゼロ）。

## 計測版ランナー（instrumented、2026-07-28 追加）

`20260715/run_system_A_localization_exchange_R_sweep_instrumented_v1.py`:
基底ランナーを無修正のまま importlib でロードし `run_case` のみ差し替えて、
各衝突時点の**倍音別複素係数**（符号付き n、eta次元込み、位相保存）を
`harmonic_dump_v1/*.npz` に追加保存する。交差相関行列 C_mn = a_m a_n^* はここから
オフライン再構成できる。追加CLI: `--dump-max-n` / `--dump-stride` / `--no-dump`。

検証済み（`20260715/instrument_check_v1/` 参照）:
既存4CSVは基底版とバイト同一、ダンプ係数からCSVの N_eff / p_chi を機械精度（〜1e-15）で再現。

## 使用上の規約

1. **原本（`波の情報読出し/` 配下）は読み取り専用として扱う。** 変更・追記・削除をしない。
2. 本フォルダ内のスクリプトへの計測点追加（倍音別複素係数・交差項のダンプ等）は、
   改変版を別名（例: `*_instrumented_v1.py`）で保存し、無修正コピーはそのまま残す。
3. 改変版を作った場合は、まず無修正コピーとの出力一致（既存列の parity）を確認してから内部量を読む。
4. 数値の独自再実装による差し替えは禁止（シリーズ内完結の再現性規約に従う）。

## 実行例

```bash
cd 20260715
# 系統A: 片側高次倍音 A=1, B=1,2 を R137近傍で
python3 run_system_A_localization_exchange_R_sweep_preliminary_v1.py \
  --packet-a 1 --packet-b 1,2 --r-values 0.697177927 --no-plots --run-id my_experiment_v1
```
