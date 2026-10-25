# 仕様書：初期化プログラム `floor_analytic_v1.py`（走行 1 の解析床）― コードと数式の対応（2026-10-25）

対象：`floor_analytic_v1.py`（SHA256 `2a2e6737de88384102db001dbae61854de142e43647d175ad7d1832f6016e1e4`、134 行）。
論文 3（2026-09-12）の `make_90deg_floor_analytic_v1.py` を写し、辺集合の選択と、辺リストからの接続行列の構成の 2 点だけを変えたもの。
乱数は使わない（決定論的）。種は置かない（Z0 ＝ 床そのもの）。

---

## 1. 入力

### 1.1 コマンドライン引数（111〜114 行）

| 位置 | 名前 | 型 | 意味 | 走行 1 での値 |
|---|---|---|---|---|
| 1 | `structure` | 文字列 | 辺集合の選択。`triu`（上三角 i＜j、M＝N(N−1)/2、論文 3 と同一）または `all`（対角込みの全順序対、M＝N²） | `all`（基準側は `triu`） |
| 2 | `outdir` | パス | 出力フォルダ | `floors_all`（基準側 `floors_triu`） |
| 3〜 | `N1 N2 …` | 整数 | 点の数 N の列 | `3 4 5 6 8 10 12 15 20 25` |

### 1.2 内部定数

| 名前 | 行 | 値 | 意味 |
|---|---|---|---|
| `CANON` | 22〜24 | 2026-09-05 正本 `run_N3_N40_stage123_v1.py` の絶対パス | 走行系の正本。ここから `one_step` と `adjacency` を取り出す |
| `PROG_SHA` | 25 | `1abf2353…` | 正本の SHA256。不一致なら 28 行の assert で停止 |
| `max_iter` | 69 | 500 | ラベル反復の上限 |

### 1.3 正本からの関数の取り出し（27〜32 行）

正本ファイルの先頭から `rows=[]; summaries=[]` の直前までを実行し（30 行）、`one_step`（31 行）と `adjacency`（32 行、triu 用の対照にのみ使用）を得る。`one_step` は本プログラムでは**残差の検算にだけ**使う（63〜66 行）。

正本の `one_step`（正本 21〜28 行）：
```
Ĥ_ef = A_ef · conj(ẑ_e) ẑ_f,  ẑ_e = e^{i arg z_e}          … H_of(exp(i·angle(z)), A)、対角 0
H₃  = i · Im Ĥ                                               … 1j*np.imag(H)
z ← V · exp(−i Δτ w) · V† z,  H₃ = V diag(w) V†,  Δτ = 2π/den … eigh と位相因子
```
すなわち K_ef ＝ A_ef sin(θ_f − θ_e) による実直交回転 z ← exp(Δτ K) z。

---

## 2. 記号と数式

- 点：0, 1, …, N−1。
- 辺（状態の添字）e ＝ (ea[e], eb[e])。配列の並びは
  - `triu`：`np.triu_indices(N, k=1)` の順（(0,1),(0,2),…,(0,N−1),(1,2),…）、M ＝ N(N−1)/2
  - `all`：`np.meshgrid(indexing='ij')` を ravel した順、e ＝ i·N ＋ j、(i, j) は対角込みの全順序対、M ＝ N²
- 接続行列（式 A）：
  ```
  A_ef = 1  if e ≠ f かつ {ea[e], eb[e]} ∩ {ea[f], eb[f]} ≠ ∅
       = 0  それ以外（三角の反対辺＝添字を共有しない対、および対角）
  ```
- ラベル（式 L）：k_e ＝ (ea[e] ＋ eb[e]) mod 4 ∈ {0, 1, 2, 3}。位相 θ_e ＝ 90°·k_e の離散化。
- 整数生成子（式 K）：
  ```
  K_ef = A_ef · s((k_f − k_e) mod 4),   s(0)=0, s(1)=+1, s(2)=0, s(3)=−1
       = A_ef · sin(90°·(k_f − k_e))
  ```
  A が対称、s が反対称（s(1)＝−s(3)）なので K は実反対称。iK はエルミート。
- 床（式 F）：iK の最大固有値 σ_max の固有ベクトル v、z ＝ v/‖v‖。
  ```
  iK z = σ_max z   ⟺   K z = −i σ_max z
  ```
  したがって連続流 dz/dτ ＝ Kz の上で z(τ) ＝ e^{−iσ_max τ} z(0)（剛体回転＝相対平衡）。
- 相対平衡の残差（式 R）：
  ```
  w = one_step(z, A, den=N),  φ = arg⟨z, w⟩,  r = ‖w − e^{iφ} z‖ / ‖z‖
  ```
  z が走行系の固定点（位相だけ進む）なら r は丸めの桁。

---

## 3. 処理の流れ（行番号 ⇔ 数式）

### S1 辺集合の生成（35〜42 行）

```python
  35  def edges_triu(N):
  36      a, b = np.triu_indices(N, k=1)
  37      return a.astype(np.int64), b.astype(np.int64)
  40  def edges_all(N):
  41      a, b = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
  42      return a.ravel().astype(np.int64), b.ravel().astype(np.int64)
```
出力 `ea, eb`（長さ M の整数配列）。`all` では ea[e] ＝ e // N、eb[e] ＝ e mod N。

### S2 接続行列（45〜52 行）＝ 式 A

```python
  45  def adjacency_from_edges(ea, eb):
  46      M = len(ea)
  47      A = np.zeros((M, M), dtype=np.float64)       # 全部 0 … 反対辺はここで 0 に決まる
  48      for e in range(M):
  49          share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])   # 行 ea[e]・行 eb[e]・列 ea[e]・列 eb[e]
  50          share[e] = False                          # 対角 0
  51          A[e, share] = 1.0                         # 添字を共有する対だけ 1
  52      return A
```
49 行は「f の行番号が e の行番号か列番号に等しい、または f の列番号が e の行番号か列番号に等しい」で、{ea[e],eb[e]} ∩ {ea[f],eb[f]} ≠ ∅ と同値（N＝5 で全組一致を確認済み）。
相手の数：`all` では非対角 (i, j) が 4N−5、対角 (i, i) が 2N−2。`triu` では 2(N−2)。

### S3 ラベルの種（72 行）＝ 式 L

```python
  72      labels = ((ea + eb) % 4).astype(int)          # 決定論的種
```

### S4 整数生成子（55〜60 行）＝ 式 K

```python
  55  def integer_K(labels, A):
  56      diff = (labels[None, :] - labels[:, None]) % 4   # diff[e,f] = (k_f − k_e) mod 4
  57      s = np.zeros_like(diff, dtype=float)
  58      s[diff == 1] = 1.0                                # sin(90°)  = +1
  59      s[diff == 3] = -1.0                               # sin(270°) = −1
  60      return A * s                                      # 反対辺（A=0）は K でも 0
```

### S5 床＝σ_max 固有ベクトル（76〜79 行）＝ 式 F

```python
  76          K = integer_K(labels, A)
  77          w, V = np.linalg.eigh(1j * K)                 # エルミート iK の固有値 w（昇順）と固有ベクトル V
  78          v = V[:, -1]                                  # 最大固有値 σ_max に属する固有ベクトル
  79          z = v / np.linalg.norm(v)                     # ‖z‖ = 1
```
`eigh` は昇順なので `V[:, -1]` が σ_max。固有値が縮退している場合の固有ベクトルの選択は LAPACK 任せ（決定論的だが環境依存）。

### S6 残差（63〜66 行、80 行）＝ 式 R

```python
  63  def releq_residual(z, A, N):
  64      w = one_step(z, A, N)                             # 走行系を 1 歩（den = N）
  65      ip = np.vdot(z, w)                                # ⟨z, w⟩
  66      return float(np.linalg.norm(w - np.exp(1j * np.angle(ip)) * z) / np.linalg.norm(z))
  80          res = releq_residual(z, A, N)
```

### S7 ラベル反復（73〜75、81〜93 行）

```python
  73      best = None                                    # (residual, z, labels)
  74      seen = {}
  75      for it in range(max_iter):
      …  S4, S5, S6 …
  81          if best is None or res < best[0]:
  82              best = (res, z.copy(), labels.copy())  # 最小残差の床を保持
  83          newlab = (np.round(np.degrees(np.angle(v)) % 360 / 90).astype(int)) % 4   # 固有ベクトルの位相を 4 象限に再量子化
  84          cand = None; agmax = -1
  85          for sh in range(4):                        # 大域位相 0/90/180/270 の 4 通り
  86              c = (newlab + sh) % 4; ag = int(np.sum(c == labels))
  87              if ag > agmax: agmax = ag; cand = c    # 前ラベルと最も一致するものを採る
  88          key = tuple(cand.tolist())
  89          if np.array_equal(cand, labels) or key in seen:
  90              break                                   # 固定点、またはゲージ等価サイクル
  91          seen[key] = it
  92          labels = cand
  93      res, z, labels = best
```
数式：k^{(t+1)} ＝ argmin_{s∈Z4} |{e : round(arg v_e^{(t)}/90°) ＋ s ≢ k_e^{(t)}}|、固定点 k^{(t+1)} ＝ k^{(t)} で停止。採用するのは反復中で残差最小の (z, k)。走行 1 では全 N で 2 回で停止（manifest の `iters`）。種ラベルがそのまま残った N と、一部の辺が k＋2 の符号ゲージに移った N がある（`seed_labels_kept`）。

### S8 検算量（94〜107 行）

```python
  94      a, b = z.real, z.imag
  95      ph = np.degrees(np.angle(z * np.exp(-1j * np.angle(z[np.argmax(np.abs(z))])))) % 360   # 最大振幅の要素を位相の基準に
  96      Kf = integer_K(labels, A)
  97      chk = dict(N=N, M=len(z), iters=it,
  98                 norm2=float(np.vdot(z, z).real),                      # ‖z‖²（=1）
  99                 centroid_abs=float(abs(z.sum())),                     # |Σ z_e|（重心）
  99                 re_closure=float(a @ a - b @ b), im_closure_cross=float(a @ b),   # Re Σz² = aᵀa − bᵀb、Im Σz²/2 = aᵀb
 100                 closure_abs=float(abs(np.sum(z * z))), releq_residual=res,         # |Σ z_e²|（二乗閉塞）、式 R
 101                 on_90deg_axes=bool(all(min(abs((p - x + 180) % 360 - 180) for x in (0, 90, 180, 270)) < 1e-8 for p in ph)),  # 全位相が軸上
 102                 n_quadrant=len(set((np.round(ph / 90).astype(int) % 4).tolist())),   # 使われている象限の数
 103                 n_distinct_amp=len(set(np.round(np.abs(z), 9).tolist())),             # 振幅の種類（9 桁丸め）
 104                 sigma=float(np.max(np.linalg.eigvalsh(1j * Kf))),                      # σ_max
 105                 rank_K=int(np.linalg.matrix_rank(Kf, tol=1e-9)),                       # rank K（平面数 × 2）
 106                 centroid_exact_zero=bool(N % 4 == 0),                                  # 論文 3 の分類（4|N で重心厳密 0）
 107                 seed_labels_kept=bool(np.array_equal(labels, ((ea + eb) % 4).astype(int))))
 108      return z, labels, chk, ea, eb
```

### S9 対照と出力（118〜134 行）

```python
 118      for N in Ns:
 119          z, labels, chk, ea, eb = build_analytic_floor(N, fn)
 120          if structure == 'triu':
 121              assert np.array_equal(adjacency_from_edges(ea, eb), canonical_adjacency(N))   # triu では正本の adjacency(N) と一致を要求
 122          name = f'floor_90deg_{structure}_N{N:05d}_analytic.npz'
 123          np.savez_compressed(os.path.join(outdir, name), v=z, g=np.zeros(len(z), dtype=np.complex128), Z0=z,
 124                              sigma=np.array([chk['sigma']]), residual=np.array(chk['releq_residual']),
 125                              n=np.array(N), seed=np.array(-1), delta=np.array(0.0), tol=np.array(0.0),
 126                              iters=np.array(chk['iters']), family=labels.astype(np.int64), theta=np.angle(z),
 127                              rule=np.array(f'theoretical_floor_90deg_integerK_sigmamax_analytic_{structure}'),
 128                              ea=ea, eb=eb, M=np.array(len(z)), structure=np.array(structure))
 129          manifest[name] = chk
 130〜132  print（1 行の要約）
 133      json.dump(manifest, …, f'manifest_floor_{structure}.json')
```

---

## 4. 出力データの構造

### 4.1 床ファイル `floors_<structure>/floor_90deg_<structure>_N<NNNNN>_analytic.npz`（N ごとに 1 本）

| キー | 形 | 型 | 意味 | 式 |
|---|---|---|---|---|
| `v` | (M,) | complex128 | 床 z（‖z‖＝1） | 式 F |
| `g` | (M,) | complex128 | 種。全て 0（種なし） | — |
| `Z0` | (M,) | complex128 | 走行の初期状態。z と同一 | Z0 ＝ z |
| `sigma` | (1,) | float64 | σ_max | 式 F |
| `residual` | () | float64 | 相対平衡の残差 r | 式 R |
| `n` | () | int64 | N | |
| `seed` | () | int64 | −1（乱数なし） | |
| `delta` | () | float64 | 0.0（種振幅なし） | |
| `tol` | () | float64 | 0.0（許容値なし、最小残差を採用） | |
| `iters` | () | int64 | ラベル反復の停止回数 | S7 |
| `family` | (M,) | int64 | 最終ラベル k_e ∈ {0,1,2,3} | 式 L＋S7 |
| `theta` | (M,) | float64 | arg z_e（rad） | |
| `rule` | () | str | `theoretical_floor_90deg_integerK_sigmamax_analytic_<structure>` | |
| `ea`, `eb` | (M,) | int64 | 辺の端点（状態の添字）。**走行本体はここから A を作る** | S1 |
| `M` | () | int64 | 状態数 | |
| `structure` | () | str | `triu` または `all` | |

走行本体 `run_sweep_from_parents_v1.py` が読むキーは `Z0`・`ea`・`eb` の 3 つ。A は走行本体の `adjacency(ea,eb)`（本プログラムの S2 と同じ規則）で再構成する。

### 4.2 manifest `floors_<structure>/manifest_floor_<structure>.json`

ファイル名をキーに、S8 の `chk` 辞書（N, M, iters, norm2, centroid_abs, re_closure, im_closure_cross, closure_abs, releq_residual, on_90deg_axes, n_quadrant, n_distinct_amp, sigma, rank_K, centroid_exact_zero, seed_labels_kept）。

### 4.3 標準出力（130〜132 行）

`N  M | 象限 軸上 振幅種 | |Σz| |Σz²| | 残差 σ rank | 反復 種ラベル保持` の 1 行／N。走行 1 の実測は `results_all/log_floors_all.txt`（all）、`results_control/log_floors_triu.txt`（triu）。

---

## 5. 走行 1 での実測値（出力の例）

| N | M | 象限 | 振幅種 | \|Σz\| | \|Σz²\| | 残差 | σ | rank K | 反復 | 種ラベル保持 |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 16 | 4 | 3 | 0.000 | 5.7e-16 | 5.8e-16 | 6.3246 | 6 | 2 | True |
| 6 | 36 | 4 | 3 | 0.333 | 7.4e-16 | 9.2e-16 | 10.3923 | 10 | 2 | False |
| 10 | 100 | 4 | 3 | 0.200 | 1.2e-15 | 1.1e-15 | 18.4391 | 18 | 2 | True |
| 25 | 625 | 4 | 5 | 0.040 | 7.3e-16 | 3.1e-15 | 48.4355 | 48 | 2 | False |

（全 N は `results_all/log_floors_all.txt`。σ² は偶 N で 4N²−6N、奇 N で 4N²−6N−4 に一致＝数値観察、rank K ＝ 2(N−1)。）

---

## 6. 床が走行系の固定点である根拠

論文 1 の恒等式：位相 θ_e で D ＝ diag(e^{iθ_e}) とすると D⁻¹ K D ＝ i W、W_ef ＝ A_ef sin²(θ_f − θ_e)。θ ＝ 90°·k なら sin² ∈ {0, 1} で W ＝ A ∘ 1[k_f − k_e が奇]。
W の Perron ベクトル r ＞ 0 に対し z ＝ D r は K z ＝ i σ z を満たし、走行系 one_step（位相のみ K、exp(ΔτK)）の下で位相だけが進む相対平衡になる。本プログラムは r を直接解かず iK の σ_max 固有ベクトルを取るが、両者は同じ床（論文 3 §2、符号ゲージを除く）。残差 r ≈ 1e-15 がその確認。

---

## 7. 注意

- 反対辺を 0 にしているのは S2 の 47 行（ゼロ初期化）と 51 行（共有対だけ 1）。他に A を書き換える箇所はない。
- 本プログラムは**初期値を作るだけ**で、走行系（one_step）は残差の検算にしか使わない。走行中の A は走行本体が同じ規則で別に作る。
- 乱数・種・許容値はない。同じ入力なら同じ出力（固有ベクトルの縮退時の選択は LAPACK の実装に依存）。
- 対角要素 (i, i)（`all` のみ）は、ラベル k ＝ 2i mod 4 ∈ {0, 2}、相手は行 i と列 i の 2N−2 本。
