# 差異分析——complex_simplex_decompactification_N5_N16_20260826（原本 / 修正版 / baseline）

**作成日:** 2026-08-28　全表は `results/three_way_comparison.md`。

## 適用した修正（全パッケージ共通、判断 1〜6）

- A1 親の振幅正規化 `v/‖v‖` を除去（親の固有モード反復は位相のみ K のまま＝S3）
- A2/A3 初期状態は親そのもの `Z = v`（外部 seed δg も正規化も無し）
- A4 動力学の生成子を振幅込み `K_ij = Im(z̄ᵢ zⱼ)`（`set_state`）。baseline は `KMODE=phase` で位相のみ K
- A5 `zero_closure_generic` の正規化除去（未使用）／A6(b) σ は実際に回している K の厳密スペクトルで読む（冪反復は実行しない）
- R1 Cayley → `exp((2π/144)·K)`（`linear_rotation_step`、`cayley_step` は定義のみ残す）／R2 K/σ 枝廃止／R3 validate 書換え、γ は dphi=2π/n_den
- S1 `zero_closure_kernel_seed` を呼ばない／S2 `sigma_max_power`・`wp=rng.normal` を呼ばない（定義は残す）

**既知の問題（親探索）**：A1 により `_eigenmode_residual` が単位ノルム前提のため残差が意味を失い、収束判定が効かず、3 リスタートの選択も無効。N=6, 7, 10, 11 では**収束していない親**（真残差 0.47 / 0.39 / 0.081 / 1.4e-3）が選ばれている。summary.json の `parent_residual`（0.49〜73）は意味のない数。残差判定のスケール不変化は次回の判断項目。

## 3. decompactification N5 / N16（記事図 1 の源）

| N | 量 | 原本 | 修正版 | baseline |
|---|---|---|---|---|
| 5 | H⊥ start | 1.13437e-31 | 2.60001e-32 | 2.60001e-32 |
|   | H⊥ max | 0.718284 | 0.429825 | 0.54952 |
|   | H⊥ final | 0.403242 | 0.0632954 | 0.364591 |
|   | onset(H⊥≥1e-8) | 242 | 1 | 316 |
|   | R_perp_log_growth_rate_per_step | 0.0862598 | – | 0.0864059 |
|   | R_perp_fit_window | [89, 275] | [None, None] | [163, 349] |
|   | A_perp_final | 0.635013 | 0.251586 | 0.603814 |
|   | R_perp_takagi_final | 0.411405 | 0.203946 | 0.392728 |
|   | parent_residual | 2.13848e-13 | 0.4869 | 0.4869 |
| 16 | H⊥ start | 2.9563e-32 | 1.25474e-31 | 1.25474e-31 |
|   | H⊥ max | 0.127514 | 3.58977 | 0.463892 |
|   | H⊥ final | 0.127351 | 3.5305 | 0.459126 |
|   | onset(H⊥≥1e-8) | 94 | 1 | 103 |
|   | R_perp_log_growth_rate_per_step | 0.23608 | – | 0.251055 |
|   | R_perp_fit_window | [36, 103] | [None, None] | [50, 113] |
|   | A_perp_final | 0.356863 | 1.87896 | 0.677588 |
|   | R_perp_takagi_final | 0.236813 | 0.978747 | 0.479323 |
|   | parent_residual | 3.02778e-13 | 73.6618 | 73.6618 |


## 読み取り

- 記事図 1（10⁻³² から 31 桁の急拡大）の源。修正版では **潜伏相が無く**（onset = step 1）、指数成長の当てはめ窓（1e-10 < R⊥ < 1e-3 が 10 点以上）が存在しない（`R_perp_fit_window = [None, None]`）。H⊥ の出発値 10⁻³² は丸め床で、修正版ではそこから即座に立ち上がる。
- N=16 修正版：H⊥ max 3.59（H_total 3.62 の 99%）。親平面にほぼ何も残らない。原本は 0.13（12.7%）。
- baseline（位相のみ K・線形回転）は原本を再現：潜伏（onset 316 / 103）、成長率 0.0864 / 0.251（原本 0.0863 / 0.236）。


## 再生成した図

`results/N16_H_components.png`, `results/N16_decompactification_scale.png`, `results/N16_full_takagi_axes.png`, `results/N16_perp_takagi_axes.png`, `results/N5_H_components.png`, `results/N5_decompactification_scale.png`, `results/N5_full_takagi_axes.png`, `results/N5_perp_takagi_axes.png`

## このパッケージへの変更箇所

```
'GAMMA = math.tan(math.pi / 144.0)' -> 'GAMMA = math.tan(math.pi / 144.0)  # 旧 Cayley 刻み（記録用。力学には使わな' x1
'def set_theta(self, theta):\n        n = self.n\n        self.' -> 'def set_theta(self, theta):\n        """位相のみ（|z|=1）の生成子。make_' x1
'def cayley_step(self, z, sigma):' -> 'def dense_K(self):\n        """現在の状態（set_state）の生成子 K を密行列で返す' x1
'v = sys_lr.w(EV[:, idx].astype(complex))\n            v = v /' -> 'v = sys_lr.w(EV[:, idx].astype(complex))  # A1: 振幅正規化を除去' x1
'Z = X + 1j * Y\n    return Z / np.linalg.norm(Z)' -> 'Z = X + 1j * Y\n    return Z  # A5: 全体の振幅正規化を除去（|X|=|Y| の相対比は' x1
'theta0 = np.angle(Z)\n    sys_lr.set_theta(theta0)\n    Kd = A' -> 'sys_lr.set_state(Z)\n    zz = Z / np.abs(Z) if KMODE == "phas' x1
'sig_lr = sys_lr.sigma_spectrum()\n    sig_d = np.sort(np.lina' -> 'sig_lr = np.sort(np.linalg.eigvalsh(1j * sys_lr.dense_K()))\n' x1
'g = zero_closure_kernel_seed(sys_lr, rng)\n    Z = v + delta ' -> 'Z = v.copy()  # A2/A3/S1: 外部 seed も正規化も無し（zero_closure_kerne' x1
'wp = rng.normal(size=sys_lr.m)\n    t0 = time.time()\n    for ' -> 't0 = time.time()  # S2: 冪反復用乱数 wp は使わない\n    for t in range(c' x1
'sys_lr.set_theta(np.angle(Z))\n        sig_est, wp = sys_lr.s' -> 'sys_lr.set_state(Z)  # A4\n        Z = sys_lr.linear_rotation' x1
'wp = rng.normal(size=m)\n    plateau_tau = None' -> 'plateau_tau = None  # S2: wp は使わない' x1
'sys_lr.set_theta(np.angle(Z))\n        sig_est, wp = sys_lr.s' -> 'sys_lr.set_state(Z)  # A4\n        if t % sub == 0:\n         ' x1
'Z = sys_lr.cayley_step(Z, sig_est)\n    t_run = time.time() -' -> 'Z = sys_lr.linear_rotation_step(Z)  # R1\n    t_run = time.ti' x1
'out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA}' -> 'out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA, "angle"' x1
'(\\w+)\\.set_theta\\(np\\.angle\\((\\w+)\\)\\)' -> '\\1.set_state(\\2)  # A4' x1
'^\\s*se, ?wpn? ?= ?\\w+\\.sigma_max_power\\(wp\\)\\n' -> '' x1
'(\\w+)\\.cayley_step\\((\\w+), ?se\\)' -> '\\1.linear_rotation_step(\\2)  # R1' x1
'^\\s*wp ?= ?rng\\.normal\\(size=\\w+\\)\\n' -> '' x1
"'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMM" -> "'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMM" x1
```
