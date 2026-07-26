# 第2予備実験 依存監査（read-only）

## 1. 第1予備実験 実行コード と 12. SHA-256

| 役割 | パス | SHA-256(先頭16) |
|:--|:--|:--|
| prelim1_run | `次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/code/run_preliminary_seed_ablation_v1.py` | `75a10a5b951302be` |
| engine | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py` | `ba0fc19b03caf06d` |
| parent_basis_exact | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_exact_v1.py` | `9cf28ca8c0d2ad8f` |
| parent_basis_approx | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_approx_v1.py` | `a9d247a8070d849f` |
| gram_dominant_plane | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/code/run_n300_dimension_saturation_v2.py` | `229938a666310574` |
| retract_source | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_transverse_stability_v1.py` | `c868884ba9ab5a10` |
| paper7_5color | `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py` | `fe5c7cbc33437890` |

（完全なSHA: `config/source_file_hashes.json`）

## 2. N=5,40,300 親状態生成法

- `make_parent(LowRankSystem(N), default_rng(40260722+1000*N), iters=1200, tol=1e-12)` → 親 v。

## 3. 初期状態 Z0=v 構築位置（run_preliminary_seed_ablation_v1）

```
5: 条件A: initial seed OFF / metastable seed OFF（Z0 = v。kernel seed 生成で乱数を消費しない）
84:         Z0 = v + DELTA * g; Z0 = Z0 / np.linalg.norm(Z0)
86:         Z0 = v.copy()                                      # 無seed（乱数を消費しない）
```

## 4,5. 初期seed/準安定seed OFF の位置

```
5: 条件A: initial seed OFF / metastable seed OFF（Z0 = v。kernel seed 生成で乱数を消費しない）
34: from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
51:     "A": {"initial_seed": False, "metastable_seed": False, "file": "condition_A_no_seed"},
52:     "B": {"initial_seed": True, "metastable_seed": False, "file": "condition_B_initial_only"},
53:     "D": {"initial_seed": True, "metastable_seed": True, "file": "condition_D_existing_two_seed"},
69: def build_init(n, initial_seed):
70:     """第7論文 build と同一。ただし initial_seed=False では kernel seed g を生成せず乱数を消費しない。"""
82:     if initial_seed:
83:         g = zero_closure_kernel_seed(sys_lr, rng)         # 乱数消費（B/D）
84:         Z0 = v + DELTA * g; Z0 = Z0 / np.linalg.norm(Z0)
86:         Z0 = v.copy()                                      # 無seed（乱数を消費しない）
98:     sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = build_init(n, c["initial_seed"])
111:     w.writerow(["step", "time", "N", "condition", "initial_seed_enabled", "metastable_seed_enabled",
112:                 "initial_seed_amplitude", "metastable_seed_amplitude", "parent_plane_occupation",
119:     init_amp = DELTA if c["initial_seed"] else 0.0
161:             w.writerow([t, t, n, cond, int(c["initial_seed"]), int(c["metastable_seed"]),
182:             "initial_seed_enabled": c["initial_seed"], "metastable_seed_enabled": c["metastable_seed"],
```

## 6,7. 1 step 更新順序（既存 evolve = Cayley のみ。閉鎖・正規化は Cayley が暗黙保存）

```
91: def evolve(sys_lr, Z, wp):
92:     sys_lr.set_theta(np.angle(Z)); se, wp = sys_lr.sigma_max_power(wp)
93:     return sys_lr.cayley_step(Z, se), wp
```

本実験では量子化後に polar retraction（下記 retract）を1回のみ適用。Cayley 直後は測定のみ。

## 8. f_outside/q3/q4/rank_Q 算出法

```
34: from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
45: Q_REL_TAU = 1e-8             # rank_Q 判定閾値（第7論文と同一）
57: def occ(B, Z):
63: def qsv4(B0, Bd):
113:                 "f_outside_parent", "q1", "q2", "q3", "q4", "rank_Q", "dominant_plane_occupation",
151:             E_P1 = occ(B_p1, Z); E_other = occ(B_rot, Z); E_ker = totZ - E_P1 - E_other
155:             rankQ = int(np.sum(qs > Q_REL_TAU * qs[0]))
163:                         fmt % qs[0], fmt % qs[1], fmt % qs[2], fmt % qs[3], rankQ,
```
- f_outside = 1 - E_P1/|Z|²（E_P1 = 固定親平面占有）。q3,q4 = svals([B0|B_dom]) の3,4。
  rank_Q = #{q_j > 1e-8 q1}。E_dom = 瞬時支配平面占有（Gram）。

## 9. crossing 判定式

```
7: 条件D: initial seed ON  / metastable seed ON （B と t1 直前までビット一致。t1=crossing+3000 で
13: 共通最終 step = 55000。crossing = f>0.05 の最初、t0/t1 = crossing+GUARD(3000)。
43: GUARD = 3000                 # crossing → t1（第7論文と同一）
116:                 "projection_closure_error", "crossing_detected", "metastable_start_detected"])
124:           "max_antisym_error": 0.0, "crossing_step": None, "metastable_start_step": None,
128:     crossing = None
132:         if crossing is None and f > 0.05:
133:             crossing = t; dg["crossing_step"] = t
134:         t1 = (crossing + GUARD) if crossing is not None else None
167:                         fmt % (residual / totZ), int(crossing is not None), met_start])
183:             "v_normalized_zero_square": v_diag, "crossing_step": dg["crossing_step"],
193:     print(f"[N={n} 条件{cond}] crossing={dg['crossing_step']} t1/準安定開始={dg['metastable_start_step']} "
```

## 10. 時刻刻み・記録間隔（本実験は §9 で新規固定）

- 本実験の保存: step 0..1000 毎step, 1001.. 5step毎, 停止step 必ず。max_step: N5=2500/N40=4500/N300=10000。

## 11. 数値型・線形代数・丸め

- float64: eps=2.220446049250313e-16, tiny(smallest normal)=2.2250738585072014e-308, smallest_subnormal=5e-324, bits=64
- Q_Δ 丸め: half_to_even（numpy round はデフォルト banker's rounding = half to even）
- platform: macOS-26.3.1-arm64-arm-64bit

## polar retraction（採用・第7論文, 不変更 import）

```
run_transverse_stability_v1.py:210: def retract(W):
```

## 監査判定

第1予備実験コード・第7論文 retract・親基底・Gram を全て確認。欠落なし。read-only import で再利用。
