# 三方向読出しの仕様 — 5色方向分解の定義・実装対応・段1+2+3 データでの読出し

（N=3..40 段1+2+3 スイープ再現論文の後続章・草稿 v0.3。DOI 未取得。
シリーズ正本: 総括＋3章 Version DOI 10.5281/zenodo.22317636 / Concept DOI
10.5281/zenodo.22317635。式番号は本稿内で独立に振る。第1章・第2章の式を参照する
場合は「第2章・式n」と明記する）

## 1. 目的

種なし系論文（第8論文、§1.1 [P6]）とその後続の論文7系図化は、図3「方向ごとの配分」
——5色占有（P1 / direction 3 / direction 4 / remaining other-rotation / kernel）——を
用いて「立ち上がるのは決まった方向だけで、残りの方向はほとんど空のまま」「第三の方向
まで含めて、三つの方向の構造が、急拡大の停止とともに定着していた」と結論した。
この図化は単なる可視化ではなく、**「方向」という状態空間の分解＝解釈装置**を測定に
持ち込んでいる。装置の数学的定義・実装行番号対応・再現ゲートを固定しなければ、
後続の解釈（三方向構造の実在性・空間読出しとの関係）の理論的背景が読めなくなる。

本稿の記述目的は次の4点の固定であり、**観測事実を超える物理的解釈は含まない**。

1. 5色方向分解の全数式（式1〜式16）と、原本プログラムの行番号レベルの対応
   （各式の直後に該当ソースコードを引用する）。
2. 装置の系譜——由来する公開論文（自己引用、§1.1）と原本プログラムの所在
   （リポジトリ相対パス、§1.2）、および本稿で作成したプログラムの完全な一覧（§1.3）。
3. 同じ読出しを段1+2+3 の 10,000 歩データ（別力学・同一親）へ適用した再現実装と
   その対照ゲート・実行結果（§5）。
4. この装置が持ち込む解釈上の前提の明示的列挙（§2.8）と、本再現実装が論文2水準の
   計装に対して残す欠落の改修要件（§8・課題）。

### 1.1 系譜と自己引用（本稿が依拠する公開論文）

| # | 題名 | Version DOI / Concept DOI | 公開日 | 本稿との関係 |
|---|---|---|---|---|
| P1 | N体完全二体関係波における生成子ランクの線形上界と空間方向読出しの三方向飽和 | 10.5281/zenodo.21465899 / 21465898 | 2026-07-21 | 式1の低ランク構造 rank K ≤ 2N と「三方向飽和」概念の初出 |
| P2 | N体固定生成子系における平面分解読出し | 10.5281/zenodo.21468960 / 21468959 | 2026-07-21 | 固定基底による平面分解読出しの前身 |
| P3 | N体関係波閉鎖系における状態の自発的分裂の開始と帰結の三分類（7月正本） | 10.5281/zenodo.21486234 / 21486233 | 2026-07-22 | 親 v・種 g・分裂量 f の正本。エンジン（§1.2 [S1]）の由来 |
| P4 | 波の数は系の分解能である（分裂停止と直交平面創発。内部系譜の「論文6」） | 10.5281/zenodo.21543071 / 21543070 | 2026-07-25 | **固定親基底3分類（P1/other/核）と f = 1−E_P1 の定義元**（原本 [S4] 冒頭「論文6の固定親基底3分類…を維持」） |
| P5 | N体関係波閉鎖系における三方向空間の創発（内部系譜の「論文7」） | 10.5281/zenodo.21578402 / 21578401 | 2026-07-26 | **5色方向分解（direction 3/4 の新設）と図1/2/3 の初出**。本稿が仕様化する装置そのもの |
| P6 | N体関係波閉鎖系における三方向生成の時間構造——二段階seed除去による因果分離（第8論文） | 10.5281/zenodo.21614403 / 21614402 | 2026-07-27 | 種なし系での「三方向構造の定着」結論（本稿冒頭の引用文の出典）。再現パッケージに原本プログラム群を同梱 |
| P7 | 自己無撞着インフレーション機構の段階分解と N=3..40 全域再現（総括＋3章。シリーズ正本） | 10.5281/zenodo.22317636 / 22317635 | 2026-09-05 | 本稿の初期データ（静的親）・状態空間・第1章/第2章の式の定義元 |

### 1.2 原本プログラムの所在（リポジトリ相対パス。すべて git 管理下の正本）

以下、`ENGINE/` = `時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1` と略記する。

| # | 役割 | 相対パス | 本稿での扱い |
|---|---|---|---|
| S1 | エンジン（LowRankSystem・make_parent・zero_closure_kernel_seed。式1・式2の実装） | `ENGINE/run_n_scaling_lowrank_v1.py` | import |
| S2 | 固定親基底（parent_plane_split_exact。式3〜式6） | `ENGINE/run_plane_flow_exact_v1.py` | import |
| S2′ | 同・近似版（parent_plane_split_approx。N=300 用、σ_rel=1e-6） | `ENGINE/run_plane_flow_approx_v1.py` | **不使用**（本稿の対象は N ≤ 40。§4.2 参照） |
| S3 | グラム縮約・支配平面（gram_reduce・dominant_plane。式7・式8。τ_G は34行） | `ENGINE/exact_lowN_eigenspectrum_v2/code/run_n300_dimension_saturation_v2.py` | import |
| S4 | 5色時系列の原本（occ・s4_new_dirs・align_2d・測定ループ。式9〜式13） | `ENGINE/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py` | 逐語コピー |
| S5 | 論文7図化の原本 | `ENGINE/exact_lowN_eigenspectrum_v2/paper7_longtime/code/make_paper7_figures.py` | 様式の参照元 |
| S6 | 図3比較図の直接の見本（control 変種。COLORS/LABELS 28-29行・suptitle 116行） | `次元の生成構造/電子の反跳実験/paper7_f_projection_v1/make_paper7_figures_control_v1.py` | 様式ゲートの基準 |
| D1 | 原本の5色時系列データ | `ENGINE/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N{00005,00040,00300}/paper7_long_timeseries.csv` | §6.3 の対比・メタ転載元 |
| D2 | 原本のメタ（crossing・dims_P1/other/kernel） | 同 `summary/N*_5color_meta.json` | §6.3 に転載 |
| D3 | 第8論文再現パッケージ内の原本コピー群 | `次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/zenodo_package_v1/nbody_two_stage_seed_removal_reproduction_v1/programs/originals_paper7/` | 公開版の所在（照合用） |
| D4 | 本稿の入力状態（10,000歩スイープ） | `…/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_long10000_20260905/results/hm_N{N}_den_{N}_states_10000.npz` | 入力（SHA ゲート） |
| D5 | 静的親（ゲート用） | `…/N3_N40_stage123_sweep_20260905/parents/parent_static_N{N:05d}_makeparent_20260905.npz` | 入力（bit ゲート） |

### 1.3 本稿で作成したプログラム・成果物の完全一覧

所在はすべて本フォルダ `paper7図3再現_stage123_20260906/`。SHA256 全桁は補遺C。

| ファイル | 種別 | 内容 |
|---|---|---|
| `make_fig3_5color_stage123_v1.py` | プログラム（新規作成） | 5色読出しの段1+2+3 データへの適用（式9〜式13 は S4 の逐語コピー、式3〜式8 は S1〜S3 の import、力学は npz 読込みに置換）。SHA/親/初期状態の3重ゲート内蔵 |
| `check_direction_axes_mapping_v1.py` | プログラム（新規作成） | 成分分解（式15）と成分別位相軸（式16）。N=3,4,40 |
| `run_all.sh` | 実行スクリプト | 上記2本の一括実行（§5.1） |
| `fig3_compare_stage123_N5_N40.png` / `.svg` | 図（出力） | §5.6 に定義 |
| `bands5_stage123_N{00005,00040}.csv` | データ（出力） | §4.2 に列定義 |
| `make_fig3_5color_stage123_v1_meta.json` | メタ（出力） | §4.2 にスキーマ |
| `check_direction_axes_mapping_v1.json` | データ（出力） | §4.2 にスキーマ |
| `README.md` / 本稿 / `SHA256SUMS.txt` | 文書・台帳 | — |

## 2. 理論的背景

記述順序は (2.1) 状態空間と頂点分解、(2.2) 固定親基底（P1/other/核）、
(2.3) グラム縮約、(2.4) 時間依存の支配平面、(2.5) 新方向 e₃e₄→f₃f₄、
(2.6) 占有・分裂量と P1 平面恒等（証明）、(2.7) 成分分解と成分別位相軸、
(2.8) 解釈上の前提、である。

### 2.1 状態空間と頂点分解（第1章・第2章と共通）

状態は完全グラフ K_N の辺上の複素ベクトル Z ∈ ℂ^M、M = N(N−1)/2、辺位相
θ_e = arg Z_e。親 v は make_parent（乱数種 40260722 + 1000N、iters=1200、tol=1e-12）
による固有状態、初期状態は Z₀ = normalize(v + δg)、δ = 10⁻¹⁵、
g = zero_closure_kernel_seed（第1章・式4〜式9）。
実ベクトル空間への射影は、複素 Z の実部・虚部それぞれに実基底を適用する（式12）。
以下「方向」とは、断りのない限り**状態空間の部分空間**を指し、複素平面の位相軸
（シリーズ既報の 60°3軸）ではない。両者の関係は §6.2 の読出し事実。

**(式1) 頂点分解と低ランク構造** — 位相差正弦生成子（第2章・式11）は頂点分解

    K = C Sᵀ − S Cᵀ = W J Wᵀ,   W = [C | S] ∈ ℝ^{M×2N},   rank K ≤ 2N

を持つ（第1章・式3。三方向飽和との関係は [P1]）。ここで C, S は辺 e=(i,j) と頂点 i の
接続に cos θ_e, sin θ_e を置いた M×N 行列、J はシンプレクティック形

    J = [[0, I_N], [−I_N, 0]]

である。エンジン実装（[S1]。段1+2+3 スイープ同梱コピーと bit 同一）:

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

W の適用と転置適用は `w`（[S1] 104-108行）・`wt`（同 100-102行）。

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

原本実装（[S2]）:

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

各 N の3分割の次元（原本メタ [D2] より転載。crossing 列は**旧エンジンでの**値であり
本稿の値と混同しないこと——§6.3）:

| N | M | dims_P1 | dims_other | dims_kernel |
|---|---|---|---|---|
| 5 | 10 | **2** | 4 | 4 |
| 40 | 780 | **2** | 78 | 700 |

### 2.3 グラム縮約（式7）

現在状態 Z の位相 θ = arg Z における縮約生成子を、密行列を作らずグラム行列経由で
構成する。

**(式7) グラム縮約** — G（式2）を対称化・固有分解し、λ > τ_G λ_max（**τ_G = 10⁻¹²**）
の成分を保持する:

    G = V diag(λ) Vᵀ,   S = V_r diag(√λ_r),   K_r = ½(Sᵀ J S − (Sᵀ J S)ᵀ),
    (μ, U_r) = eigh(i K_r)

K = WJWᵀ かつ G = WᵀW より、K_r は K を W の主成分座標に縮約した実反対称行列で、
μ は K の非零固有値 ±σ の実数化に一致する（丸め誤差は診断量として記録）。

原本実装（[S3]。τ_G は34行 `TAU_G_PRIMARY = 1e-12`）:

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

B₀ は式8 を親 Z=v に適用したもの（B_P1 と同一平面。原本 [S4] 56-57行
`gr0 = gram_reduce(sys_lr, v)` → `dominant_plane`）。

### 2.5 新方向 e₃, e₄ → f₃, f₄（式9〜式11）

**(式9) 新方向（P1 直交補内）** — 時間依存の支配平面を親平面に直交化した2方向:

    R = B_dom − B₀ (B₀ᵀ B_dom),   e₃₄ = QR(R) の先頭2列

**(式10) other 空間への射影** — e₃₄ を固定 other 空間へ射影して正規直交化:

    f₃₄ = QR( B_rot (B_rotᵀ e₃₄) ) の先頭2列

**(式11) 連続整列（規約）** — 前時刻の f₃₄ との 2×2 直交プロクラステス整列:

    O = f_prevᵀ f_new,   O = UΣVᵀ（SVD）,   f_new ← f_new (UVᵀ)ᵀ

**この整列は d3/d4 の個別ラベルの割り振りにのみ影響し、平面 span(f₃, f₄) と
E_d3 + E_d4 には影響しない**（直交回転だから）。

原本実装（[S4]）:

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

### 2.6 占有・分裂量・5色と P1 平面恒等（式12〜式14）

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

5色の名称対応は direction 1, 2 = B_P1 の2列、direction 3, 4 = f₃₄ の2列、
remaining other-rotation = E_rem、kernel = E_ker（f = 1−E_P1 の定義元は [P4]）。

原本の測定ループ（[S4]）:

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

**(式14) P1 平面恒等 — f = H⊥/H の証明** — シリーズの休眠比（第2章・式22〜式23）は
p = Re Z₀/‖Re Z₀‖, q = GS(Im Z₀) の張る平面 Π への補直交エネルギー比である。
親 v は tol = 10⁻¹² で K(arg v) の複素固有ベクトル（第1章・式4〜式6）:

    K v = −i σ_max v   ⟹   K (Re v) = σ_max (Im v),   K (Im v) = −σ_max (Re v)

すなわち span_ℝ(Re v, Im v) は σ_max の回転平面そのものである。式4の群化で最大 σ 群が
**非縮退（dims_P1 = 2）**ならば、B_P1 はこの平面の正規直交基底であり、
δ = 10⁻¹⁵ より Z₀ ≈ v だから

    span_ℝ(p, q) = span_ℝ(Re v, Im v) = span_ℝ(B_P1)   ⟹   f = 1 − E_P1/‖Z‖² = H⊥/H

が**平面の一致として厳密に**成立する。非縮退の前提は実測で満たされる
（§2.2 の表: dims_P1 = 2、N=5, 40 とも）。縮退時（dims_P1 > 2）は B_P1 ⊋ Π となり
f ≤ H⊥/H に落ちるため、その場合は本恒等は使えない（§8・課題C8）。
実測上も crossing（f > 0.05 の最初の step）はスイープの onset（第2章・式24）と
一致する（§5.4・G4）。

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

1. **「方向」は状態空間の部分空間であり、複素平面の位相軸ではない。**
   図3の「方向ごとの配分」を空間方向として読むには、位相軸との対応を別途確立する
   必要がある。実測（§6.2）ではこの対応は自明でなく、60°3軸は成分間干渉が生成する。
2. **P1/other/核の3分割は親位相での K のスペクトル群化（式4、6桁丸め）に依存する。**
   丸め桁は規約であり、σ の縮退構造が桁境界にかかる系では分割が変わり得る。
3. **B_dom は時間依存で、τ_G = 10⁻¹²（式7）に依存する。** τ_G の掃引系譜は原本
   TAU_G_SWEEP（[S3] 35行）にあるが、本読出しは PRIMARY のみを使う。
4. **連続整列（式11）は d3/d4 の個別ラベルの規約**であり、物理量は
   span(f₃,f₄) と E_d3+E_d4 のみが持つ。d3 と d4 を個別に解釈してはならない。
5. **シリーズ★★★留意（等振幅収束の循環導出の疑い）は本読出しにも波及する。**
   段2内在エンジンの軌道上の読出しであり、方向構造をノルム形式の独立傍証として
   提示してはならない。

## 3. 実装方法

再現実装 `make_fig3_5color_stage123_v1.py` の系譜は3種に区別される
（系列規約: 過去論文依拠はコピー→対照テスト→import）:

- **import（式1〜式8）**: 基底構成は原本モジュール [S1][S2][S3] をそのまま import。
- **逐語コピー（式9〜式13）**: `occ` / `s4_new_dirs` / `align_2d` と測定ループの式は
  [S4] 42-45・66-80・113-129 行のコピー（コピー先 60-78・106-127 行。差分は変数名
  `se_ev`→`SAMPLE` 等の周辺のみで、式に関わる行は同一）。
- **置換（力学→読込み）**: 原本の時間発展（[S4] 145行 `set_theta → sigma_max_power →
  cayley_step`）を、段1+2+3 の 10,000 歩状態 npz の読込み `Zr = Zs[t]`（105行）に
  置換した。**力学は一切走らせない**。原本はさらに crossing を旧力学の事前走行で
  求める（[S4] 93-97行）:

```python
    93	    Zc = Z.copy(); wpc = wp.copy(); crossing = None; t = 0
    94	    while True:
    95	        if fval(Zc) > 0.05:
    96	            crossing = t; break
    97	        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc); Zc = sys_lr.cayley_step(Zc, se); t += 1
```

  本再現では同じ定義（f > 0.05 の最初の step）を読出し列上で求める（109-110行）。
  したがって crossing の**定義**は同一だが、**値**は力学が違うため原本と異なる
  （旧エンジン 1167/2011 ↔ 本再現 64/358。§6.3 の対比表）。

**式⇔実装対応表**:

| 式 | 内容 | 原本実装 | 再現実装（make_fig3_..._v1.py） |
|---|---|---|---|
| 式1 | K = WJWᵀ 頂点分解 | [S1] 64-66, 110-114 | import |
| 式2 | G = WᵀW | [S1] 79-90（set_theta 内） | import |
| 式3〜6 | 固定親基底 P1/other/核 | [S2] 45-77 | import（88行で呼出し） |
| 式7 | グラム縮約（τ_G=1e-12） | [S3] 34, 44-67 | import（115行） |
| 式8 | 支配平面 B_dom | [S3] 70-84 | import（116行） |
| 式9 | 新方向 e₃₄ | [S4] 66-70 | 65-69行（コピー） |
| 式10 | other 射影 f₃₄ | [S4] 123-125 | 118-120行（コピー） |
| 式11 | 連続整列 | [S4] 73-80 | 71-78行（コピー） |
| 式12 | 占有 E_B | [S4] 42-45 | 60-63行（コピー） |
| 式13 | 5色と f | [S4] 113-129 | 106-127行（コピー） |
| 式14 | P1 平面恒等 | —（本稿で証明） | crossing 判定 109-110行が帰結を使用 |
| 式15〜16 | 成分分解・成分別軸 | —（本稿新設） | check_direction_axes_mapping_v1.py 45-73, 90-111 |

## 4. 詳細設計

### 4.1 全体フロー（make_fig3_5color_stage123_v1.py）

```
for N in {5, 40}:
  [ゲート]  入力 npz の SHA256 照合（§5.4 G1）→ 基底構成（式3〜式8）
           → 親 v・Z0・states[0] の bit 一致検証（G2）
  [読出し]  t=0..10000: f（式13）を毎 step、5色を 25 step ごとに算出
  [出力]   bands CSV・meta JSON
[図化]     2段比較図（§5.6 の定義）
```

check_direction_axes_mapping_v1.py は N ∈ {3, 4, 40} の step 10000 に対し
式15・式16 を適用して JSON を出力する（同じ G1/G2 ゲート内蔵）。

### 4.2 全体データフロー

**入力**（詳細パスは §1.2 [D4][D5]）:

| ファイル | キーと形 | dtype |
|---|---|---|
| `hm_N{N}_den_{N}_states_10000.npz` | `Z` (10001, M)・`N`・`denominator`・`steps` | complex128 / int64 |
| `parent_static_N{N:05d}_makeparent_20260905.npz` | `v` (M)・`Z0` (M)（ゲート照合に使用） | complex128 |
| `../N3_N40_long10000_20260905/SHA256SUMS.txt` | 入力 npz の SHA256 正本（55件） | — |

**パラメータ**（すべてプログラム冒頭の定数。コマンドライン引数は持たない）:

| 定数 | 値 | 意味 | 原本との対応 |
|---|---|---|---|
| `DELTA` | 1e-15 | 種の振幅 δ | [S4] 36行と同一 |
| `XMAX` | 10000 | 読出し範囲 | 原本は 55000（データ長の違い） |
| `SAMPLE` | 25 | 5色の標本間隔 | [S4] 39行 SAMPLE{5:25, 40:25} と同一 |
| `NS` | [5, 40] | 対象 N | 原本は {5,40,300}。**N=300 は本スイープ（N≤40）に存在しないため対象外**。成分分解は {3,4,40}（結晶2系＋ガラス代表） |
| `COLORS`/`LABELS` | [S6] 28-29行と同一 | 5色の配色・凡例 | 様式ゲート |
| `amp_floor_rel` | 1e-6 | 式16 の微小成分除外 | 本稿新設 |
| 軸クラスタ閾値 | 0.15 rad | 式16 の分割 | 本稿新設 |

**出力**:

| ファイル | 形式 | 内容 |
|---|---|---|
| `bands5_stage123_N{00005,00040}.csv` | ヘッダ1行＋データ401行（t=0,25,…,10000）。列: `step, P1, dir3, dir4, remaining_other, kernel, f`。数値は `%.10e` | 式13 の時系列 |
| `make_fig3_5color_stage123_v1_meta.json` | キー: `{N: {M, crossing, final_bands_P1_d3_d4_rem_ker(5要素), final_f}}` | 終端値と crossing |
| `check_direction_axes_mapping_v1.json` | キー: `{N: {成分名: {energy_frac, axes_deg_weight_count[(軸角deg, 重み, 本数)]}}}`、成分名 ∈ {total, P1, new_plane_d3d4, remaining_other, kernel} | 式15・式16 の結果 |
| `fig3_compare_stage123_N5_N40.png/.svg` | §5.6 | 比較図 |

### 4.3 個別処理

#### 4.3.1 ゲート（§5.4 の実装）

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

## 5. 実行結果

### 5.1 再現コマンド

```bash
cd paper7図3再現_stage123_20260906
python3 make_fig3_5color_stage123_v1.py      # 5色読出し＋比較図（引数なし）
python3 check_direction_axes_mapping_v1.py   # 成分分解・成分別軸（引数なし）
# または ./run_all.sh（上記2本を順に実行）
```

### 5.2 実行環境

- Python 3.9.6、numpy 2.0.2（BLAS/LAPACK: macOS Accelerate）、matplotlib（Agg）
- macOS 26.3.1（arm64）

### 5.3 実行時間

run_all.sh 全体（fig3 N=5,40 ＋ 成分分解 N=3,4,40）で **実測 3.5 秒**（real、
2026-09-06）。支配項は状態 npz の展開読込みと、固定親基底の密 `eig(K)`（M=780、
式3。[S2] 53-54行）。5色ループ自体はグラム縮約が 2N×2N = 80×80 の `eigh` で済む
ため軽い（式7）。

### 5.4 検証ゲート

| ゲート | 合格条件 | 実測 | 合否 |
|---|---|---|---|
| G1: 入力 SHA | 入力 npz の SHA256 が 10,000歩スイープの台帳と一致 | fig3 2件＋成分分解 3件、全一致 | **PASS** |
| G2: 親同一性（bit） | make_parent の v・Z₀ = normalize(v+10⁻¹⁵g) が静的親 npz の `v`・`Z0` および `Z[0]` と bit 一致 | N ∈ {3,4,5,40} 全一致（`[gate]` 出力） | **PASS** |
| G3: 完走 | 両プログラムが `ALL DONE`・終了コード 0 | 確認（再実行 2026-09-06 でも同一） | **PASS** |
| G4: crossing = onset | f > 0.05 の最初の step が第2章スイープの onset と一致（式14 の帰結） | N=5: 64 = 64、N=40: 358 = 358 | **PASS** |

（注: G4 は現状、meta JSON とスイープ summary の目視突合である。機械突合は §8 課題C7。）

### 5.5 データ

| 項目 | 内容 |
|---|---|
| bands CSV | 402行（ヘッダ＋401）×2本（4.2 の列定義） |
| meta JSON | 530 bytes（4.2 のスキーマ） |
| mapping JSON | 2,910 bytes |
| 図 | png 66,850 bytes / svg 76,102 bytes |
| 本体プログラム SHA256 | `6e0054831ad2bf1ce06442df7598926af93c1afa3efbaedbbf10bc0aa08e0b66  make_fig3_5color_stage123_v1.py` |
| 全ファイルの正本値 | 同梱 `SHA256SUMS.txt`（全桁の転載は補遺C） |

### 5.6 図化の定義

`fig3_compare_stage123_N5_N40.png/.svg`: 2段グリッド（figsize 11×7.2、dpi 130、
sharex）。各段の**縦軸は5色占有比（式13）の log**（下限クリップ 10⁻⁶。旧論文の
control 図 [S6] と同じ床）、**黒線は f = 1−E_P1/‖Z‖²**、点線は crossing 位置。
横軸は絶対 step 0..10000（目盛 1000 刻み）。配色・凡例は [S6] 28-29行の
COLORS/LABELS と同一（P1=青 #4C78A8、d3=赤 #E45756、d4=橙 #F58518、
remaining=灰 #B0B0B0、kernel=緑 #54A24B）。

## 6. 実行分析（客観的報告と観察のみ）

### 6.1 5色時系列（数値は meta JSON より転載）

| N | M | crossing | E_P1 | E_d3 | E_d4 | E_rem | E_ker | f 終端 |
|---|---|---|---|---|---|---|---|---|
| 5 | 10 | 64 | 0.59724 | 0.11423 | 0.07744 | 2.14×10⁻¹² | 0.21109 | 0.40276 |
| 40 | 780 | 358 | 0.89095 | 0.02183 | 0.02500 | 2.78×10⁻¹⁷ | 0.06221 | 0.10905 |

- **remaining other-rotation は両系とも実質 0**: 垂直成長は新平面 f₃₄ ＋核に
  閉じ込められ、other 空間の残り（N=40 では 78−2 = 76 次元）へは漏れない。
- 旧エンジンの図（[S6] の figures_control/figure3_compare_N5_N40_N300.png）と
  定性的に同一の構造（P1 高位・d3≈d4 の立ち上がりと定着・灰の消滅・核の定常帯）。
  **三方向構造はエンジン・読出し実装に依存しない**。

### 6.2 成分分解と成分別位相軸（mapping JSON より転載、step 10000）

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
  どの成分単独の軸集合にも含まれず**、成分の波ごとの和＝干渉が生成する。
- ガラス（N=40）: Z_new・Z_ker は各々**単一軸のコヒーレントベクトル**だが、その軸は
  親軸の近傍（176.59° は 4.03°+180° から 7.4°、2.82° は 1.2°）にあり、総和に第3の
  位相軸を作らない。

### 6.3 原本（旧エンジン）との対比表（混同防止）

| 量 | 旧エンジン（[S4] 実測、[D2] より転載） | 本再現（段1+2+3、10,000歩） |
|---|---|---|
| 力学 | σ正規化 Cayley 時計（第2章・式12〜14） | 固定 Δτ=2π/N スペクトル写像（第2章・式16〜20） |
| crossing N=5 / N=40 | 1167 / 2011 | 64 / 358 |
| dims P1/other/kernel（N=5） | 2/4/4 | 同一（親が bit 一致・式は同一） |
| dims P1/other/kernel（N=40） | 2/78/700 | 同一 |
| 比較してよい対象 | —— | **構造（どの方向が立ち上がるか）のみ。crossing の値・時間スケールは時計が違うため比較不可** |

## 7. 主張範囲と反証可能性

本稿の主張は次に限定する。

1. 5色方向分解は式1〜式14 で完全に定義され、原本実装（§1.2）と行番号レベルで
   対応する（§2〜§3）。
2. 再現実装は §5.4 のゲートの下で正本と同一の装置であり、力学のみを段1+2+3 データの
   読込みに置換した（置換点は 4.3.2 の105行のみ）。
3. §6 の読出し事実。

「方向」の物理的意味づけ（空間次元との同一視・凝縮体との対応・第3方向の実在論）は
本稿の主張に含めない。反証条件: (i) 同一入力に対し引用行のコードが本稿の式と異なる
出力を与えること、(ii) §6 の数値が `run_all.sh` の再実行で再現しないこと、
(iii) §5.4 のゲートのいずれかが不一致になること。

## 8. 課題 — 再現実装の改修要件（論文2水準の計装との差分）

本稿の再現実装は式・ゲートは論文2水準だが、**計装（診断・メタデータ・集計）が
論文2の本体プログラム群に対して不足している**。以下を改修要件として固定する
（各項に合格条件を付す。改修実施時は v1→v2 の最小差分と再実行・SHA 更新を伴う）。

| # | 要件 | 合格条件 |
|---|---|---|
| C1 | dtype 断言の追加（論文2 の12行 `assert` 相当） | complex128/float64 の itemsize 断言がプログラム冒頭にあること |
| C2 | RUN_METADATA 出力（numpy/python/BLAS・全パラメータ・入力 npz の SHA256 を機械記録） | JSON が出力され §5.2 の記載と一致すること |
| C3 | 閉塞診断の記録: osum = (E_P1+E_d3+E_d4+E_rem+E_ker)/‖Z‖² と max\|osum−1\| を全標本で記録（原本 [S4] 132-134行相当。現状は E_rem の max(0,·) が誤差を吸収し検査不能） | bands CSV に osum 列が追加され、max 閉塞誤差が meta に記録されること |
| C4 | 原本 CSV 16列との列対応の復元: direction_1/2 個別・plane_1/plane_2・norm_error を追加し、原本 [S4] 101-106行のヘッダと対応表を文書化 | 列対応表が本稿 §4.2 に反映されること |
| C5 | ゲート結果の機械記録（現状 print のみ）: G1〜G4 の実測値を JSON に保存 | §5.4 の表が JSON からの転載になること |
| C6 | 毎 step の f 系列の成果物化（現状 crossing 判定に使うが 25 step 間引きでしか保存されず、CSV から crossing を再検証できない） | 毎 step f の CSV（10,001行）が出力されること |
| C7 | 集計プログラムの新設: crossing vs onset の機械突合・E_rem の全時間最大・d3/d4 の床離脱 step を JSON 出力（§6 の記載を全て転載制に） | §6 の全数値に出所 JSON があること |
| C8 | dims_P1/other/kernel の自前記録（式14 の非縮退前提の検査。現状は原本メタ [D2] からの転載） | 実行時に B_p1.shape[1] 等が meta に記録され、dims_P1 ≠ 2 の場合に警告すること |

## 補遺A プログラム系譜

```
[P2] 平面分解読出し → [P4] 固定親基底3分類（P1/other/核、f=1−E_P1）
  └─ [P5] 論文7 5色分解: run_paper7_5color_timeseries.py [S4]（新方向 e3,e4 を追加。原本）
       └─ 図化: make_paper7_figures.py [S5] → paper7_f_projection_v1/ の5変種
            （control 版 [S6] = 本稿が様式を継承した figure3_compare の生成元）
                 └─ 本稿: make_fig3_5color_stage123_v1.py
                     （基底構成=import [S1-S3]、5色式=逐語コピー [S4]、
                      力学→10,000歩npz読込みに置換）
                     check_direction_axes_mapping_v1.py（式15〜16 を新設）
```

## 補遺B 同梱物と入力の正本性

本フォルダの SHA256SUMS.txt を正とする。入力データの正本性は
`../N3_N40_long10000_20260905/SHA256SUMS.txt`（状態 npz 38本を含む55件、npz は
Drive 正本）、親の正本性は `../N3_N40_stage123_sweep_20260905/SHA256SUMS.txt` に依る。

## 補遺C SHA256（全桁、shasum -a 256 の出力を転載。本稿自身の行は台帳参照）

```
6e0054831ad2bf1ce06442df7598926af93c1afa3efbaedbbf10bc0aa08e0b66  make_fig3_5color_stage123_v1.py
33cf1e1ff016ee795158d0acb0bae9c32fa5e1296939fb7212e19a755c4b5875  check_direction_axes_mapping_v1.py
dfba3f33c3631d4c34d271e3ab2e3d3b3dd951fccad78f8f36032fea1909db24  run_all.sh
76f65f75e5c325b479ee0cdafad568c508b4eb7e7970249cf3c8adbc5b8f19e3  fig3_compare_stage123_N5_N40.png
c312925b6f3dc363321ffaf39dea0e0aa2b36a14b356d559fadb2a6dc67bc7da  fig3_compare_stage123_N5_N40.svg
cbf732a681b619b77da4b5df682534e2b0983cfe2d450dde4f422105218139c3  bands5_stage123_N00005.csv
08fc868b3b8562343edbb99984cda0c18a4f9323e574509569bb5d3fd728cb41  bands5_stage123_N00040.csv
f8c8555cdff168da45b799ce421bb374a858a863ede7c1d38c90c42560ce3e51  make_fig3_5color_stage123_v1_meta.json
fe396e3113a8a95023e74bb6f1f9961c9567e92a9df21f7c9baec04c6a66ec12  check_direction_axes_mapping_v1.json
```

---
（三方向読出しの仕様 v0.3 おわり）
