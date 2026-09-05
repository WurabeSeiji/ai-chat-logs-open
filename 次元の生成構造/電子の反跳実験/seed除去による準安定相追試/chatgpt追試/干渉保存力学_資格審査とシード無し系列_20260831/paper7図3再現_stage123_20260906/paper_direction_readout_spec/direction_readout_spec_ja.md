# 三方向読出しの仕様 — 5色方向分解の定義・実装対応・段1+2+3 データでの読出し

（N=3..40 段1+2+3 スイープ再現論文の後続章・草稿 v0.2。DOI 未取得。
シリーズ正本: 総括＋3章 Version DOI 10.5281/zenodo.22317636 / Concept DOI
10.5281/zenodo.22317635。式番号は本稿内で独立に振る。第1章・第2章の式を参照する
場合は「第2章・式n」と明記する）

## 1. 目的

種なし系論文（第8論文系、Version DOI 10.5281/zenodo.21614402）とその後続の論文7系
図化は、図3「方向ごとの配分」——5色占有（P1 / direction 3 / direction 4 /
remaining other-rotation / kernel）——を用いて「立ち上がるのは決まった方向だけで、
残りの方向はほとんど空のまま」「第三の方向まで含めて、三つの方向の構造が、急拡大の
停止とともに定着していた」と結論した。この図化は単なる可視化ではなく、
**「方向」という状態空間の分解＝解釈装置**を測定に持ち込んでいる。装置の数学的定義・
実装行番号対応・再現ゲートを固定しなければ、後続の解釈（三方向構造の実在性・
空間読出しとの関係）の理論的背景が読めなくなる。

本稿の記述目的は次の3点の固定であり、**観測事実を超える物理的解釈は含まない**。

1. 5色方向分解の全数式（式1〜式14）と、原本プログラムの行番号レベルの対応
   （各式の直後に該当ソースコードを引用する）。
2. 同じ読出しを段1+2+3 の 10,000 歩データ（別力学・同一親）へ適用した再現実装
   `make_fig3_5color_stage123_v1.py` と、その対照ゲート。
3. この装置が持ち込む解釈上の前提の明示的列挙（§2.8）。

## 2. 理論的背景

記述順序は (2.1) 状態空間と頂点分解、(2.2) 固定親基底（P1/other/核）、
(2.3) グラム縮約、(2.4) 時間依存の支配平面、(2.5) 新方向 e₃e₄→f₃f₄、
(2.6) 占有と分裂量、(2.7) 成分分解と成分別位相軸、(2.8) 解釈上の前提、である。

### 2.1 状態空間と頂点分解（第1章・第2章と共通）

状態は完全グラフ K_N の辺上の複素ベクトル Z ∈ ℂ^M、M = N(N−1)/2、辺位相
θ_e = arg Z_e。親 v は make_parent（乱数種 40260722 + 1000N、iters=1200、tol=1e-12）
による固有状態、初期状態は Z₀ = normalize(v + δg)、δ = 10⁻¹⁵、
g = zero_closure_kernel_seed（第1章・式4〜式9）。

**(式1) 頂点分解と低ランク構造** — 位相差正弦生成子（第2章・式11）は頂点分解

    K = C Sᵀ − S Cᵀ = W J Wᵀ,   W = [C | S] ∈ ℝ^{M×2N},   rank K ≤ 2N

を持つ（第1章・式3）。ここで C, S は辺 e=(i,j) と頂点 i の接続に cos θ_e, sin θ_e を
置いた M×N 行列、J はシンプレクティック形

    J = [[0, I_N], [−I_N, 0]]

である。エンジン実装（`第5論文原本_自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py`。
段1+2+3 スイープ同梱コピーと bit 同一）:

```python
    64	        self.J = np.zeros((2 * n, 2 * n))
    65	        self.J[:n, n:] = np.eye(n)
    66	        self.J[n:, :n] = -np.eye(n)
```

```python
   110	    def kmatvec(self, z):
   111	        """K z = C(S^T z) - S(C^T z)"""
   112	        vs = self.vsum(self.s * z)
   113	        vc = self.vsum(self.c * z)
   114	        return self.c * (vs[self.ea] + vs[self.eb]) - self.s * (vc[self.ea] + vc[self.eb])
```

**(式2) グラム行列** — G = WᵀW ∈ ℝ^{2N×2N}（対称・半正定値）。`set_theta` が
θ 設定時に cos/sin のブロック（Gcc = CᵀC, Gcs = CᵀS, Gss = SᵀS）として構成する:

```python
    79	        Gcc = CT * CT
    80	        Gcs = CT * ST
    81	        Gss = ST * ST
    82	        np.fill_diagonal(Gcc, Gcc.sum(axis=1))
    83	        np.fill_diagonal(Gcs, Gcs.sum(axis=1))
    84	        np.fill_diagonal(Gss, Gss.sum(axis=1))
    85	        G = np.empty((2 * n, 2 * n))
    86	        G[:n, :n] = Gcc
    87	        G[:n, n:] = Gcs
    88	        G[n:, :n] = Gcs
    89	        G[n:, n:] = Gss
    90	        self.G = G
```

W の適用と転置適用は `w`（104-108行）・`wt`（100-102行）。

### 2.2 固定親基底 — P1・other 回転・核（式3〜式6）

親位相 θ⁰ = arg v で生成子 K を**密行列として**構成し、固有分解する。

**(式3) 親スペクトル** — K は実反対称なので固有値は純虚数対 ±iσ_k（σ_k ≥ 0）。
各固有ベクトル V_k（複素）の実部・虚部が σ_k の回転平面の実基底を張る。

**(式4) スペクトル群化（規約）** — 正の σ を **round(σ, 6)**（6桁丸め）で群化する。
丸めで 0 になる成分は核へ。連続的な閾値パラメータは持たず、明確なスペクトル
ギャップに依存する。

**(式5) P1 基底** — 最大 σ 群の {Re V, Im V} を QR 正規直交化（|R_jj| > 10⁻⁸ の列のみ
採用）したものを B_P1 とする。

**(式6) other 回転基底と核** — 残りの正 σ 群を合併し、P1 直交補へ射影
（R₀ − B_P1 B_P1ᵀ R₀）してから QR したものを B_rot。**核は基底を作らず**残差
E_ker = ‖Z‖² − E_P1 − E_other で扱う。

原本実装（`第5論文原本_自発的分裂予備実験_v1/run_plane_flow_exact_v1.py`）:

```python
    45	def parent_plane_split_exact(sys_lr, v, round_digits=6):
    46	    """密行列 eig(K) で (支配平面P₁, その他回転平面) の正規直交実基底 + σスペクトル。
    47	
    48	    分離は round(σ, round_digits) のグループ化のみ（明確なギャップに依存し、
    49	    連続的な閾値パラメータを持たない）。核基底は作らず残差で扱う。
    50	    """
    51	    M = sys_lr.m
    52	    sys_lr.set_theta(np.angle(v))
    53	    K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
    54	    w, V = np.linalg.eig(K)
    55	    sig_all = w.imag
    56	    groups = {}
    57	    for i in range(M):
    58	        key = round(float(sig_all[i]), round_digits)
    59	        if key > 0:                      # 正の虚部＝回転平面（round で 0 は核へ）
    60	            groups.setdefault(key, []).extend([np.real(V[:, i]), np.imag(V[:, i])])
    61	
    62	    def ortho(cols):
    63	        Q, R = np.linalg.qr(np.column_stack(cols))
    64	        return Q[:, np.abs(np.diag(R)) > 1e-8]
    65	
    66	    sig_sorted = sorted(groups, reverse=True)
    67	    p1_sigma = sig_sorted[0]
    68	    B_p1 = ortho(groups[p1_sigma])
    69	    rest = [c for k in sig_sorted[1:] for c in groups[k]]
    70	    B_rot = None
    71	    if rest:
    72	        R0 = np.column_stack(rest)
    73	        R0 = R0 - B_p1 @ (B_p1.T @ R0)
    74	        Q, R = np.linalg.qr(R0)
    75	        B_rot = Q[:, np.abs(np.diag(R)) > 1e-8]
    76	    spectrum = sorted((round(float(s), 6) for s in sig_all if s > 1e-10), reverse=True)
    77	    return p1_sigma, B_p1, B_rot, spectrum
```

**例外的挙動（式6の帰結）** — N=3 では正 σ 群が1つしかなく、70行の `B_rot = None` の
まま返る（M=3 = P1 の2次元＋核1次元）。**したがって N=3 に direction 3/4 は定義されず、
垂直成長はすべて核に計上される**。N=4 でも本稿の読出しでは E_new = 0 が実測される
（§6.2）。この場合分けは再現実装 91-104 行（§4.3.3）で明示的に扱う。

### 2.3 グラム縮約（式7）

現在状態 Z の位相 θ = arg Z における縮約生成子を、密行列を作らずグラム行列経由で
構成する。

**(式7) グラム縮約** — G（式2）を対称化・固有分解し、λ > τ_G λ_max（**τ_G = 10⁻¹²**）
の成分を保持する:

    G = V diag(λ) Vᵀ,   S = V_r diag(√λ_r),   K_r = ½(Sᵀ J S − (Sᵀ J S)ᵀ),
    (μ, U_r) = eigh(i K_r)

K = WJWᵀ かつ G = WᵀW より、K_r は K を W の主成分座標に縮約した実反対称行列で、
μ は K の非零固有値 ±σ の実数化に一致する（丸め誤差は診断量として記録）。

原本実装（`exact_lowN_eigenspectrum_v2/code/run_n300_dimension_saturation_v2.py`。
τ_G は34行 `TAU_G_PRIMARY = 1e-12`）:

```python
    44	def gram_reduce(sys_lr, Z, tau_G=TAU_G_PRIMARY):
    45	    sys_lr.set_theta(np.angle(Z))
    46	    G = 0.5 * (sys_lr.G + sys_lr.G.T)
    47	    Gnorm = np.linalg.norm(G) + 1e-300
    48	    lam, V = np.linalg.eigh(G)
    49	    order = np.argsort(-lam)
    50	    lam = lam[order]; V = V[:, order]
    51	    lmax = lam[0]
    52	    keep = lam > tau_G * lmax
    53	    lam_r = lam[keep]; V_r = V[:, keep]; sq = np.sqrt(lam_r)
    54	    S = V_r * sq[None, :]
    55	    Kr = S.T @ sys_lr.J @ S
    56	    Kr = 0.5 * (Kr - Kr.T)
    57	    Krnorm = np.linalg.norm(Kr) + 1e-300
    58	    mu, Ur = np.linalg.eigh(1j * Kr)
```

（59-67行は対称性誤差・再構成誤差・縮約反対称誤差・保持ランク r_G の診断記録。）

### 2.4 時間依存の支配平面（式8）

**(式8) 支配平面 B_dom(t)** — 最大 μ のモード u_r を辺空間へ持ち上げ

    u = W (V_r (u_r / √λ_r)) ∈ ℂ^M,   B_dom = QR[√2 Re u, √2 Im u] の先頭2列

とする。持ち上げ固有対残差 ‖Ku + iσu‖/|σ| と基底直交誤差を診断量として持つ。

```python
    70	def dominant_plane(sys_lr, gr):
    71	    mu, Ur, V_r, sq = gr["mu"], gr["Ur"], gr["V_r"], gr["sq"]
    72	    idx = int(np.argmax(mu))                 # 最大正 μ = σ_dom
    73	    sig = float(mu[idx]); ur = Ur[:, idx]
    74	    coef = V_r @ (ur / sq)                    # 2N 複素
    75	    u = sys_lr.w(coef.real) + 1j * sys_lr.w(coef.imag)   # 辺空間 M
    76	    b = np.column_stack([np.sqrt(2) * u.real, np.sqrt(2) * u.imag])
    77	    Q, R = np.linalg.qr(b)
    78	    B = Q[:, :2]
    79	    # 持ち上げ固有対残差（密行列を作らず Ku=WJ(Wᵀu)）
    80	    Wt_u = sys_lr.wt(u)                       # 2N 複素
    81	    Ku = sys_lr.w((sys_lr.J @ Wt_u.real)) + 1j * sys_lr.w((sys_lr.J @ Wt_u.imag))
    82	    lifted_res = float(np.linalg.norm(Ku + 1j * sig * u) / max(1.0, abs(sig)))
    83	    orthB = float(np.linalg.norm(B.T @ B - np.eye(2)))
    84	    return sig, B, u, lifted_res, orthB
```

B₀ は式8 を親 Z=v に適用したもの（B_P1 と同一平面。原本
`run_paper7_5color_timeseries.py` 56-57行 `gr0 = gram_reduce(sys_lr, v)` →
`dominant_plane`）。

### 2.5 新方向 e₃, e₄ → f₃, f₄（式9〜式11）

**(式9) 新方向（P1 直交補内）** — 時間依存の支配平面を親平面に直交化した2方向:

    R = B_dom − B₀ (B₀ᵀ B_dom),   e₃₄ = QR(R) の先頭2列

**(式10) other 空間への射影** — e₃₄ を固定 other 空間へ射影して正規直交化:

    f₃₄ = QR( B_rot (B_rotᵀ e₃₄) ) の先頭2列

**(式11) 連続整列（規約）** — 前時刻の f₃₄ との 2×2 直交プロクラステス整列:

    O = f_prevᵀ f_new,   O = UΣVᵀ（SVD）,   f_new ← f_new (UVᵀ)ᵀ

**この整列は d3/d4 の個別ラベルの割り振りにのみ影響し、平面 span(f₃, f₄) と
E_d3 + E_d4 には影響しない**（直交回転だから）。

原本実装（`paper7_longtime/code/run_paper7_5color_timeseries.py`）:

```python
    66	def s4_new_dirs(B0, Bdom):
    67	    """S4=orthonormalize[B0|Bdom] の B0 直交補2方向 e3,e4 を返す。"""
    68	    R = Bdom - B0 @ (B0.T @ Bdom)
    69	    Qr, _ = np.linalg.qr(R)
    70	    return Qr[:, :2]
```

```python
    73	def align_2d(f_prev, f_new):
    74	    """f_new(M×2) を前時刻 f_prev へ 2×2 回転で整列（連続基底固定・色反転防止）。"""
    75	    if f_prev is None:
    76	        return f_new
    77	    Ov = f_prev.T @ f_new                # 2×2
    78	    U, _, Vt = np.linalg.svd(Ov)
    79	    Rot = U @ Vt                          # 直交 2×2
    80	    return f_new @ Rot.T
```

### 2.6 占有・分裂量・5色（式12〜式14）

**(式12) 占有** — 実正規直交基底 B に対する複素状態 Z の占有を、実部・虚部への
同時適用で定義する:

    E_B(Z) = ‖Bᵀ Re Z‖² + ‖Bᵀ Im Z‖²

これは B の張る実部分空間の複素化 span_ℂ への直交射影のエネルギーに等しい。

```python
    42	def occ(B, Z):
    43	    if B is None or (hasattr(B, "shape") and B.shape[1] == 0):
    44	        return 0.0
    45	    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))
```

**(式13) 5色と分裂量** —

    E_P1 = E_{B_P1},   E_other = E_{B_rot},   E_ker = ‖Z‖² − E_P1 − E_other
    E_d3 = E_{f₃},   E_d4 = E_{f₄},   E_rem = max(0, E_other − E_d3 − E_d4)
    f = 1 − E_P1 / ‖Z‖²

f は第6論文と同一で、シリーズの休眠比 H⊥/H（第2章・式23）と同じ対象を固定親基底側
から読んだものである（p, q は Re v, Im v の正規直交化＝B_P1 と同一平面。実測でも
crossing = onset が一致する: §6.1）。5色の名称対応は direction 1, 2 = B_P1 の2列、
direction 3, 4 = f₃₄ の2列、remaining other-rotation = E_rem、kernel = E_ker。

原本の測定ループ（同プログラム）:

```python
   113	        if t % se_ev == 0 or t == XMAX:
   114	            totZ = float(np.real(np.conj(Zr) @ Zr))
   115	            E_P1 = occ(B_p1, Zr)
   116	            E_other = occ(B_rot, Zr)
   117	            E_ker = totZ - E_P1 - E_other
   118	            f = 1.0 - E_P1 / totZ
   119	            # 支配平面(gram) → 新方向 e3,e4 → other空間へ射影 f3,f4
   120	            gr = gram_reduce(sys_lr, Zr)
   121	            _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
   122	            e34 = s4_new_dirs(B0, Bdom)              # M×2, P1直交補
   123	            proj = B_rot @ (B_rot.T @ e34)           # other空間へ射影
   124	            fq, _ = np.linalg.qr(proj)
   125	            f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
   126	            f34 = align_2d(f_prev, f34); f_prev = f34
   127	            E_d3 = occ(f34[:, [0]], Zr)
   128	            E_d4 = occ(f34[:, [1]], Zr) if f34.shape[1] > 1 else 0.0
   129	            E_rem_other = max(0.0, E_other - E_d3 - E_d4)
```

### 2.7 成分分解と成分別位相軸（本稿の追加読出し・式15〜式16）

**(式15) 直交成分分解** — 実基底 B の複素適用 Z_B = B(BᵀZ) により

    Z = Z_P1 + Z_new + Z_rem + Z_ker,
    Z_P1 = B_P1(B_P1ᵀZ),  Z_new = f₃₄(f₃₄ᵀZ),  Z_rem = B_rot(B_rotᵀZ) − Z_new,
    Z_ker = Z − B_P1(B_P1ᵀZ) − B_rot(B_rotᵀZ)

**(式16) 成分別位相軸** — 各成分の波ごとの軸角 arg(z_e²)/2（mod 180°）を、
|z_e| > 10⁻⁶·max|z| の成分についてクラスタ化（2θ 空間の円環ギャップ > 0.15 rad で
分割）し、軸角・エネルギー重み・本数を記録する。実装は
`check_direction_axes_mapping_v1.py`:

```python
    45	def axes_of(z, amp_floor_rel=1e-6):
    46	    """位相軸（deg, mod 180）とその重み（Σ|z|²比）。微小成分は除外。"""
    47	    a = np.abs(z)
    48	    if a.max() == 0:
    49	        return []
    50	    keep = a > amp_floor_rel * a.max()
    51	    zz = z[keep]
    52	    ph2 = np.angle(zz ** 2)
    53	    order = np.argsort(ph2)
    54	    ps = ph2[order]
    55	    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    56	    cut = np.flatnonzero(gaps > 0.15)
```

（クラスタ集計は 57-73 行。）成分分解の実装は同プログラム 90-104 行（§4.3.3 に引用）。

### 2.8 この装置が持ち込む解釈上の前提（明示）

1. **「方向」は状態空間 ℝ^M（の複素化）の部分空間であり、複素平面の位相軸ではない。**
   図3の「方向ごとの配分」を空間方向として読むには、位相軸との対応を別途確立する
   必要がある。実測（§6.2）ではこの対応は自明でなく、60°3軸は成分間干渉が生成する。
2. **P1/other/核の3分割は親位相での K のスペクトル群化（式4、6桁丸め）に依存する。**
   丸め桁は規約であり、σ の縮退構造（第2章・2-9 の縮退整数対など）が桁境界に
   かかる系では分割が変わり得る。
3. **B_dom は時間依存で、τ_G = 10⁻¹²（式7）に依存する。** τ_G の掃引は原本の
   TAU_G_SWEEP（同ファイル35行）に系譜があるが、本読出しは PRIMARY のみを使う。
4. **連続整列（式11）は d3/d4 の個別ラベルの規約**であり、物理量は
   span(f₃,f₄) と E_d3+E_d4 のみが持つ。d3 と d4 を個別に解釈してはならない。
5. **シリーズ★★★留意（等振幅収束の循環導出の疑い）は本読出しにも波及する。**
   段2内在エンジンの軌道上の読出しであり、方向構造をノルム形式の独立傍証として
   提示してはならない。

## 3. 実装方法

再現実装 `make_fig3_5color_stage123_v1.py` の系譜は3種に区別される
（系列規約: 過去論文依拠はコピー→対照テスト→import）:

- **import（式1〜式8）**: 基底構成は原本モジュールをそのまま import する。
  `run_n_scaling_lowrank_v1`（LowRankSystem, make_parent, zero_closure_kernel_seed）、
  `run_plane_flow_exact_v1`（parent_plane_split_exact）、
  `run_n300_dimension_saturation_v2`（gram_reduce, dominant_plane）。
- **逐語コピー（式9〜式13）**: `occ` / `s4_new_dirs` / `align_2d` と測定ループの式は
  `run_paper7_5color_timeseries.py` 42-45・66-80・113-129 行のコピー
  （コピー先 60-75・113-127 行。差分は変数名 `se_ev`→`SAMPLE` 等の周辺のみで、
  式に関わる行は同一）。
- **置換（力学→読込み）**: 原本の時間発展（145行 `set_theta → sigma_max_power →
  cayley_step`）を、段1+2+3 の 10,000 歩状態 npz の読込み `Zr = Zs[t]`（105行）に
  置換した。**力学は一切走らせない**。

**式⇔実装対応表**:

| 式 | 内容 | 原本実装 | 再現実装（make_fig3_..._v1.py） |
|---|---|---|---|
| 式1 | K = WJWᵀ 頂点分解 | run_n_scaling_lowrank_v1.py 64-66, 110-114 | import |
| 式2 | G = WᵀW | 同 79-90（set_theta 内） | import |
| 式3〜6 | 固定親基底 P1/other/核 | run_plane_flow_exact_v1.py 45-77 | import（88行で呼出し） |
| 式7 | グラム縮約（τ_G=1e-12） | run_n300_dimension_saturation_v2.py 34, 44-67 | import（115行） |
| 式8 | 支配平面 B_dom | 同 70-84 | import（116行） |
| 式9 | 新方向 e₃₄ | run_paper7_5color_timeseries.py 66-70 | 65-69行（コピー） |
| 式10 | other 射影 f₃₄ | 同 123-125 | 118-120行（コピー） |
| 式11 | 連続整列 | 同 73-80 | 71-78行（コピー） |
| 式12 | 占有 E_B | 同 42-45 | 60-63行（コピー） |
| 式13 | 5色と f | 同 113-129 | 106-127行（コピー） |
| 式15〜16 | 成分分解・成分別軸 | —（本稿新設） | check_direction_axes_mapping_v1.py 45-73, 90-111 |

## 4. 詳細設計

### 4.1 全体フロー（make_fig3_5color_stage123_v1.py）

```
for N in {5, 40}:
  [ゲート]  入力 npz の SHA256 照合（式なし・§5）→ 基底構成（式3〜式8）
           → 親 v・Z0・states[0] の bit 一致検証
  [読出し]  t=0..10000: f を毎 step、5色（式13）を 25 step ごとに算出
  [出力]   bands CSV・meta JSON
[図化]     2段比較図（5色 log・床 1e-6・crossing 点線・共通横軸 0..10000）
```

### 4.2 全体データフロー

- **入力**:
  - `../N3_N40_long10000_20260905/results/hm_N{5,40}_den_{5,40}_states_10000.npz`
    （`Z`(10001×M)。SHA256 は同フォルダ SHA256SUMS.txt 正本と照合）
  - `../N3_N40_stage123_sweep_20260905/parents/parent_static_N{00005,00040}_makeparent_20260905.npz`
    （ゲート用 `v`, `Z0`）
- **パラメータ**: DELTA = 1e-15（原本と同一）、XMAX = 10000、SAMPLE = 25
  （原本の SAMPLE{5:25, 40:25} と同一）、NS = [5, 40]
- **出力**: `fig3_compare_stage123_N5_N40.png/.svg`、
  `bands5_stage123_N{00005,00040}.csv`（列: step, P1, dir3, dir4, remaining_other,
  kernel, f）、`make_fig3_5color_stage123_v1_meta.json`

### 4.3 個別処理

#### 4.3.1 ゲート（§5 の実装）

```python
    55	def gate_npz(rel):
    56	    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    57	    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
```

```python
    95	    par = np.load(PARENTS / f'parent_static_N{n:05d}_makeparent_20260905.npz')
    96	    assert np.array_equal(np.asarray(par['v']), v), f'PARENT GATE FAIL v N={n}'
    97	    assert np.array_equal(np.asarray(par['Z0']), Z0), f'PARENT GATE FAIL Z0 N={n}'
    98	    assert np.array_equal(Zs[0], Z0), f'STATE0 GATE FAIL N={n}'
```

#### 4.3.2 読出しループ（式13。力学の置換点は105行）

```python
   104	    for t in range(0, XMAX + 1):
   105	        Zr = Zs[t]
   106	        totZ = float(np.real(np.conj(Zr) @ Zr))
   107	        E_P1 = occ(B_p1, Zr)
   108	        f = 1.0 - E_P1 / totZ
   109	        if crossing is None and f > 0.05:
   110	            crossing = t
   111	        if t % SAMPLE != 0 and t != XMAX:
   112	            continue
   113	        E_other = occ(B_rot, Zr)
   114	        E_ker = totZ - E_P1 - E_other
   115	        gr = gram_reduce(sys_lr, Zr)
   116	        _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
   117	        e34 = s4_new_dirs(B0, Bdom)
   118	        proj = B_rot @ (B_rot.T @ e34)
   119	        fq, _ = np.linalg.qr(proj)
   120	        f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
   121	        f34 = align_2d(f_prev, f34); f_prev = f34
   122	        E_d3 = occ(f34[:, [0]], Zr)
   123	        E_d4 = occ(f34[:, [1]], Zr) if f34.shape[1] > 1 else 0.0
   124	        E_rem = max(0.0, E_other - E_d3 - E_d4)
```

原本との差: (i) crossing は原本では旧力学の事前走行（93-97行）で求めるが、本再現では
同じ定義（f > 0.05 の最初の step）を読出し列上で求める（109-110行）。定義は同一で、
実測値はスイープの onset（第2章・式24）と一致する（§6.1）。(ii) f は毎 step 算出
（クロッシング判定のため）、5色は 25 step ごと（原本と同じ標本間隔）。

#### 4.3.3 成分分解（式15。check_direction_axes_mapping_v1.py）

```python
    90	    Z_P1 = B_p1 @ (B_p1.T @ Z)
    91	    if B_rot is not None and B_rot.shape[1] > 0:
    92	        gr = gram_reduce(sys_lr, Z)
    93	        _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
    94	        e34 = s4_new_dirs(B0, Bdom)
    95	        proj = B_rot @ (B_rot.T @ e34)
    96	        fq, _ = np.linalg.qr(proj)
    97	        f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
    98	        Z_other = B_rot @ (B_rot.T @ Z)
    99	        Z_new = f34 @ (f34.T @ Z)
   100	        Z_rem = Z_other - Z_new
   101	    else:
   102	        # N=3: other 回転空間は空（M=3 = P1 2次元＋核1次元）。垂直成長は核に住む
   103	        Z_other = np.zeros_like(Z); Z_new = np.zeros_like(Z); Z_rem = np.zeros_like(Z)
   104	    Z_ker = Z - Z_P1 - Z_other
```

**停止条件・例外的挙動（数式化しない要素）**:
- 早期停止はない。全読出しは固定 10,001 step の走査（104行）。
- N=3 の B_rot = None 分岐（101-103行）は式6の帰結（§2.2）。
- 大 N で numpy の RuntimeWarning（matmul）が表示され得るが警告であり、結果への
  影響がないことは本系列の bit 一致ゲート群で確認済み（第2章 §4.3.2 と同じ注記）。

## 5. 対照ゲート（再現性）

1. **入力 SHA ゲート**: 10,000 歩スイープ npz を同フォルダ SHA256SUMS.txt と照合
   （N=5, 40、および成分分解では N=3, 4, 40。全一致）。
2. **親同一性ゲート（bit）**: 旧フレームワークの `make_parent`（乱数種
   40260722+1000N）で構成した v・Z₀ = normalize(v+10⁻¹⁵g) が、段1+2+3 スイープの
   静的親 npz の `v`・`Z0`、および状態列 `states[0]` と **bit 一致**
   （N=3, 4, 5, 40 で確認。両フレームワークの初期データは同一物である）。
3. **様式ゲート**: 図の様式（5色配色・log・床 10⁻⁶・crossing 点線・共通横軸）は
   正本図化 `make_paper7_figures_control_v1.py`（116行 suptitle、COLORS/LABELS
   28-29行）と同一。

## 6. 読出しの結果（解釈なし。数値は meta JSON / mapping JSON から転載）

### 6.1 5色時系列（fig3_compare_stage123_N5_N40.png）

| N | M | crossing | E_P1 | E_d3 | E_d4 | E_rem | E_ker | f 終端 |
|---|---|---|---|---|---|---|---|---|
| 5 | 10 | 64 | 0.5972 | 0.1142 | 0.0774 | 2.1×10⁻¹² | 0.2111 | 0.4028 |
| 40 | 780 | 358 | 0.8910 | 0.0218 | 0.0250 | 2.8×10⁻¹⁷ | 0.0622 | 0.1090 |

- crossing は 10,000 歩スイープの onset（N=5: 64、N=40: 358）と一致する。
- **remaining other-rotation は両系とも実質 0**（床 10⁻⁶ 未満）: 垂直成長は
  新平面 f₃₄ ＋核に閉じ込められ、other 空間の残り（N=40 では 776−2 次元）へは
  漏れない。
- 旧エンジン（σ時計＋Cayley、`figures_control/figure3_compare_N5_N40_N300.png`）と
  定性的に同一の構造（P1 高位・d3≈d4 の立ち上がりと定着・灰の消滅・核の定常帯）。
  **三方向構造はエンジン・読出し実装に依存しない**。

### 6.2 成分分解と成分別位相軸（check_direction_axes_mapping_v1.json、step 10000）

| N | 成分 | E 比 | 位相軸（deg mod 180 / 重み / 本数） |
|---|---|---|---|
| 3 | 総和 | 1 | 20.51/0.333/1・80.51/0.333/1・140.51/0.333/1（60°3軸） |
| 3 | P1 | **0.833333** | 80.51/0.4/1・170.51/0.6/2（90°2軸） |
| 3 | new | 0 | —（B_rot 空） |
| 3 | 核 | **0.166667** | 80.51/1.0/2（1軸） |
| 4 | 総和 | 1 | 29.26・89.26・149.26 各 0.333（60°3軸） |
| 4 | P1 | **0.833333** | 89.26/0.4/2・179.26/0.6/4（90°2軸） |
| 4 | new | 0 | — |
| 4 | 核 | **0.166667** | 89.26/1.0/4（1軸） |
| 40 | 総和 | 1 | 6.75/0.477/372・96.77/0.515/402（2軸支配）＋微小3クラスタ |
| 40 | P1 | 0.890951 | 4.03/0.494・94.51/0.506（90°2軸） |
| 40 | new | 0.046838 | 176.59/1.0/780（単一軸） |
| 40 | 核 | 0.062211 | 2.82/1.0/780（単一軸） |

- 結晶（N=3, 4）: E_P1 : E_ker = 5/6 : 1/6（表示桁で厳密）。**総和の 60°3軸は
  どの成分単独の軸集合にも含まれず**（P1 は 90°2軸、核は1軸）、成分の波ごとの和＝
  干渉が生成する。
- ガラス（N=40）: Z_new・Z_ker は各々**単一軸のコヒーレントベクトル**だが、その軸は
  親軸の近傍（176.59° は 4.03°+180° から 7.4°、2.82° は 1.2°）にあり、総和に第3の
  位相軸を作らない。

## 7. 主張範囲と反証可能性

本稿の主張は次に限定する。

1. 5色方向分解は式1〜式13 で完全に定義され、原本実装と行番号レベルで対応する（§2〜§3）。
2. 再現実装は §5 のゲートの下で正本と同一の装置であり、力学のみを段1+2+3 データの
   読込みに置換した（置換点は 4.3.2 の105行のみ）。
3. §6 の読出し事実。

「方向」の物理的意味づけ（空間次元との同一視・凝縮体との対応・第3方向の実在論）は
本稿の主張に含めない。反証条件: (i) 同一入力に対し引用行のコードが本稿の式と異なる
出力を与えること、(ii) §6 の数値が `run_all.sh` の再実行で再現しないこと、
(iii) §5 のゲートのいずれかが不一致になること。

## 補遺A プログラム系譜

```
論文6 固定親基底3分類（P1/other/核）
  └─ 論文7 5色分解: run_paper7_5color_timeseries.py（新方向 e3,e4 を追加。原本）
       └─ 図化: make_paper7_figures.py → paper7_f_projection_v1/ の5変種
            （control 版 = 本稿が様式を継承した figure3_compare の生成元）
                 └─ 本稿: make_fig3_5color_stage123_v1.py
                     （基底構成=import、5色式=逐語コピー、力学→10,000歩npz読込みに置換）
                     check_direction_axes_mapping_v1.py（式15〜16 を新設）
```

原本の所在: `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂
予備実験_v1/`（エンジン・plane_flow）、同 `exact_lowN_eigenspectrum_v2/code/`
（グラム縮約）、同 `.../paper7_longtime/code/`（5色時系列）。

## 補遺B 同梱物と台帳

本フォルダ `paper7図3再現_stage123_20260906/` の SHA256SUMS.txt を正とする
（プログラム2本・図 png/svg・bands CSV 2本・meta/mapping JSON・README・本稿）。
入力データの正本性は `../N3_N40_long10000_20260905/SHA256SUMS.txt`（状態 npz 38本を
含む55件、npz は Drive 正本）に依る。
