# 第7論文 依存関係監査（第8論文 Phase 1・read-only）

PAPER7 engine root: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1`

## 1. 第7論文 実行コードの場所 と 10. SHA-256

| 役割 | パス | SHA-256(先頭16) |
|:--|:--|:--|
| engine | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py` | `ba0fc19b03caf06d` |
| parent_basis_exact | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_exact_v1.py` | `9cf28ca8c0d2ad8f` |
| parent_basis_approx | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_approx_v1.py` | `a9d247a8070d849f` |
| gram_dominant_plane | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/code/run_n300_dimension_saturation_v2.py` | `229938a666310574` |
| paper7_5color | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py` | `fe5c7cbc33437890` |
| paper7_transverse | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_transverse.py` | `ac1073bea329971d` |
| paper7_transverse_cached | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_transverse_cached.py` | `897a20d0e4d6a2cf` |
| paper7_exact_vs_approx | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_exact_vs_approx_N40.py` | `2127c2f35e33bb60` |
| paper7_figures | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/make_paper7_figures.py` | `273fec25c2e4c30c` |

（完全なSHA-256は `config/source_file_hashes.json`。）

## 2. N=5,40,300 のパラメータ / 9. 乱数生成器と seed値

- 乱数: `numpy.random.default_rng(40260722 + 1000*N)`（build 内, N毎に固定）
- 初期微小種 δ = 1e-15
- 記録間隔 SAMPLE = {5:25, 40:25, 300:100}
- build 関連行:
```
36: DELTA = 1e-15
37: XMAX = 55000                     # 共通横軸（絶対step）
39: SAMPLE = {5: 25, 40: 25, 300: 100}
50:     rng = np.random.default_rng(40260722 + 1000 * n)
```

## 3. 初期seedを加えている正確なコード位置（run_paper7_5color_timeseries.py, build）

```
31: from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
36: DELTA = 1e-15
58:     g = zero_closure_kernel_seed(sys_lr, rng)
59:     Z = v + DELTA * g; Z = Z / np.linalg.norm(Z)
```

## 4. 準安定域で横摂動seedを加えている位置（run_paper7_transverse.py / _cached.py）

run_paper7_transverse.py:
```
6: S4(t)^⊥ で測る。複数seed・複数eps。Benettin型再正規化で λ⊥,max^(4)。§11で機械判定。
9: Benettin: ΔTごとに δ⊥ を S4^⊥ へ再射影し ε へ再正規化、g_k=‖δ⊥‖/ε、λ=Σlog g_k/(KΔT)。
29: DT = 500                    # Benettin 再正規化間隔
75:     S4_t0 = s4_basis(sys_lr, B0, Z0)
82:                 "local_transverse_growth_rate", "renormalization_factor",
89:         eta_r = rng_dir.normal(size=M); eta_i = rng_dir.normal(size=M)
90:         eta_r = eta_r - S4_t0 @ (S4_t0.T @ eta_r); eta_i = eta_i - S4_t0 @ (S4_t0.T @ eta_i)
91:         nrm = np.sqrt(eta_r @ eta_r + eta_i @ eta_i)
92:         eta = (eta_r + 1j * eta_i) / nrm
94:             # 基準・摂動 二軌道を t0 から XMAX まで、DT ごとに Benettin 再正規化
96:             Zt = Z0 + eps * eta; Zt = Zt / np.linalg.norm(Zt)
118:                 # Benettin 再正規化：δ⊥ を S4^⊥ で取り出し ε へ戻す
```
run_paper7_transverse_cached.py:
```
81:             Zt = Z0 + eps * eta; Zt = Zt / np.linalg.norm(Zt); wpt = wp0.copy()
105:                 Zt = Zb_c + eps * (dp / max(gnorm, 1e-300)); Zt = Zt / np.linalg.norm(Zt)
```

## 5. 初期状態 / crossing / 準安定開始 / 最終 の定義

- 初期状態: build() で `Z = (v + δ g)/‖·‖`（δ=1e-15）。無seedなら `Z0 = v`。
- crossing: 分裂量 f = 1 - E_P1 が初めて 0.05 を超える step。
- 準安定開始 t0: crossing + GUARD（transverse で GUARD=3000）。
- 最終: 絶対 step XMAX = 55000。

定義行:
```
[crossing_fval>0.05]
  95:         if fval(Zc) > 0.05:
[GUARD(metastable start offset)]
  30: GUARD = 3000                # crossing 後、準安定域開始までのガード
[XMAX(final step)]
  37: XMAX = 55000                     # 共通横軸（絶対step）
[DT(Benettin interval)]
  29: DT = 500                    # Benettin 再正規化間隔
[SEEDS/EPS]
  32: SEEDS = 3
  33: EPS = [1e-8, 1e-10, 1e-12, 1e-14]
```

## 6. 既存の自然軌道データ の場所 と 列

- N=5: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/paper7_long_timeseries.csv` exists=True
- N=40: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00040/paper7_long_timeseries.csv` exists=True
- N=300: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00300/paper7_long_timeseries.csv` exists=True

N=5 列: ['step', 'time', 'crossing_flag', 'splitting_fraction', 'direction_1_occupation', 'direction_2_occupation', 'direction_3_occupation', 'direction_4_occupation', 'other_rotating_occupation', 'kernel_occupation', 'occupation_sum', 'plane_1_occupation', 'plane_2_occupation', 'norm_error', 'conservation_error', 'projection_closure_error']

## 7. 既存の二段階seedあり実験データ の場所 と 列

- N=5: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/transverse_stability_timeseries.csv` exists=True
- N=40: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00040/transverse_stability_timeseries.csv` exists=True
- N=300: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00300/transverse_stability_timeseries.csv` exists=True

N=5 列: ['step', 'time', 'seed', 'epsilon', 'baseline_norm', 'perturbed_norm', 'total_difference', 'transverse_difference', 'normalized_transverse_amplification', 'local_transverse_growth_rate', 'renormalization_factor', 'active_subspace_dimension', 'projection_closure_error', 'norm_error', 'conservation_error']

## 8. 再現に必要な Python パッケージとバージョン

- numpy: 2.0.2
- scipy: 1.13.1
- matplotlib: 3.9.4

## §9 要求列と既存自然軌道CSVの差分（不足=第8論文側で新規記録が必要）

既存 paper7_long_timeseries.csv に**存在しない** §9 列: ['N', 'condition', 'initial_seed_enabled', 'metastable_seed_enabled', 'initial_seed_amplitude', 'metastable_seed_amplitude', 'parent_plane_occupation', 'f_outside_parent', 'q1', 'q2', 'q3', 'q4', 'rank_Q', 'dominant_plane_occupation', 'non_dominant_occupation', 'residual_occupation', 'norm_Z', 'dagger_norm_error', 'zero_square_real', 'zero_square_imag', 'zero_square_abs', 'crossing_detected', 'metastable_start_detected']

→ 条件A・Bは第8論文ラッパー `run_preliminary_seed_ablation_v1.py` で §9 全列を新規記録する。
→ 条件D（初期ON＋準安定ON）は、既存 transverse CSV が §9 列を持たないため、同一コード・同一設定で
   §9 列を記録して再生成する（指示書§6.2の「必要列が存在しない場合のみ再生成」に該当。理由: 上記列差分）。

## COMMON_FINAL_STEP

- COMMON_FINAL_STEP = 55000（第7論文 XMAX と同一）。

## 監査判定

第7論文コード・自然軌道データ・二段階seedデータの所在を全て確認。欠落なし。
第7論文コードは read-only import で再利用（コピー不要）。seed の ON/OFF は第8論文ラッパーで明示切替。
