# 第2章 スイープ本体 — 段1+2+3 力学と N=3..40 走行（H⊥/H 分母コントロール図）

（N=3..40 段1+2+3 スイープ再現論文・第2章。総括論文と Concept DOI を共有する。
式番号は第1章（式1〜式9）からの連番）

## 1. 目的

本数値実験系列の中心目的は、**インフレーション的な発展——測定にかからないほど小さい種
（休眠比 H⊥/H ~ 10⁻³⁰）が、力学だけで何十桁も指数増幅して飽和に至る現象——が
起こる機構を明確にすること**である。7月正本（N=40, 300, 1000）でこの現象が観測されて
以降、系列は「どの力学構成要素がこの発展の必要条件か」を一因子ずつ切り分けてきた
（その確定結果が段1+2+3 の構成であり、監査記録は補遺A）。

本章はその到達点として、第1章で生成・検証した静的親38個（N=3..40）を初期データに、
段1+2+3 力学（位相のみ・虚部のみ生成子、実直交回転、固定 Δτ=2π/den）を
N=3..40 × 6分母 × 500 step で走行し、インフレーション的発展の指標である
**休眠比 H⊥/H**（定義と指標性の根拠は §2.5）・全エネルギー H_total・零閉塞 |zᵀz|/H を
全 step 記録して、この発展が N の全域と時計（分母）の掃引に対してどう現れるかを
目標図 `fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png` に固定する。

本論文の記述目的は完全な再現性と数式・プログラム・データの対応の固定であり、
観測事実を超える物理的解釈は含まない。段1+2+3 がこの形である根拠は、本文の流れを
妨げないよう**補遺A（削除対照の監査表）・補遺B（プログラム系譜）**に分離して記載する。

## 2. 理論的背景

本章の力学は、7月正本（旧プログラム）の力学を出発点とし、それを新アーキテクチャ上で
**段1・段2・段3** の3つの構成要素に分解して再構成したものである。したがって理論的背景は
(2.1) 共通の状態空間、(2.2) 旧プログラムの力学の数学、(2.3) 本プログラムの力学の数学
（段分解）、(2.4) 新旧の回転写像の数学的差異、(2.5) 測定量と保存則、の順で記述する。

### 2.1 状態空間と隣接構造（新旧共通）

**(式10) 線グラフ隣接行列** — 完全グラフ K_N の辺 e=(i,j), f=(k,l) が頂点を共有するとき
A_ef = 1、それ以外 0、対角 0（M×M、実対称）。辺順序は式1と同一の上三角順。
状態は z ∈ ℂ^M（M = N(N−1)/2）、辺位相は θ_e = arg z_e。

### 2.2 旧プログラム（7月正本）の力学の数学

旧プログラムは本パッケージ同梱のエンジン `run_n_scaling_lowrank_v1.py`（正本の bit 同一
コピー）で定義される。1ステップは次の3段の合成である。

**(式11) 旧生成子（位相差正弦・実反対称）** — 頂点共有辺対 (e,f) 上で

    K_ef(θ) = cos θ_e sin θ_f − sin θ_e cos θ_f = sin(θ_f − θ_e)

生成子は毎ステップ現在位相 θ = arg z から作り直される。振幅 |z_e| は入らない。
頂点分解 K = C Sᵀ − S Cᵀ = W J Wᵀ（rank ≤ 2N）は第1章の式3と同一。

**(式12) σ_max の冪反復推定（近似・履歴依存）** — 実ベクトル wp を持ち回り、

    wp ← −K(K wp)/‖−K(K wp)‖   （3回反復）
    σ̂_max = ‖K wp‖

−K² は半正定値対称でその最大固有値は σ_max² だから、σ̂_max は σ_max の**冪反復による
近似値**である。反復はわずか3回で、wp は前ステップから引き継がれる（warm start）。
したがって旧力学は厳密には z の写像ではなく、**隠れ状態 wp を持つ系**であり、σ̂_max(t) は
軌道履歴に依存する近似列である。

エンジン実装（`run_n_scaling_lowrank_v1.py`）:

```python
   122	    def sigma_max_power(self, wp, iters=3):
   123	        """warm-start 冪反復による σ_max 推定。wp は前ステップのベクトル（更新して返す）。"""
   124	        for _ in range(iters):
   125	            y = self.kmatvec(wp)
   126	            wp = -self.kmatvec(y)
   127	            nrm = np.linalg.norm(wp)
   128	            if nrm == 0.0:
   129	                return 0.0, wp
   130	            wp = wp / nrm
   131	        sig = np.linalg.norm(self.kmatvec(wp))
   132	        return float(sig), wp
```

**(式13) 旧の1ステップ写像（Cayley 変換・σ正規化時計）** — γ = tan(π/144) を定数として

    K̃ = K / σ̂_max
    z ← (I − γ K̃)⁻¹ (I + γ K̃) z

エンジン実装（Woodbury による 2N×2N 解）:

```python
   134	    def cayley_step(self, z, sigma):
   135	        """z ← (I-γK̃)^{-1}(I+γK̃) z, K̃ = K/σ。Woodbury で O(N^3)。"""
   136	        gn = GAMMA / sigma
   137	        r = z + gn * self.kmatvec(z)
   138	        A2 = (sigma / GAMMA) * self.J + self.G
   139	        rhs = self.wt(r)
   140	        y = np.linalg.solve(A2, rhs)
   141	        return r - self.w(y)
```

**(式14) 旧の固有平面回転角** — K の固有値は純虚数対 ±iσ_k（σ_k ≥ 0）。K̃ の固有値
±iσ_k/σ̂_max に対する Cayley 因子は

    (1 + iγσ_k/σ̂_max) / (1 − iγσ_k/σ̂_max)   （絶対値 1）

すなわち固有平面 k の1ステップ回転角は

    φ_k = 2 arctan( γ · σ_k / σ̂_max )

最速平面（σ_k = σ̂_max）の回転角は常に φ_max = 2 arctan(tan(π/144)) = **π/72** に固定される。
これは「系の最速モードで時間を刻む」**適応的（σ正規化）時計**であり、回転角の
スペクトル依存は arctan による**非線形圧縮**を受ける。

**(式15) 旧の保存量** — K は実反対称なので Cayley 変換 O = (I−γK̃)⁻¹(I+γK̃) は
**実直交行列**（Oᵀ O = I）。実直交行列は z = x + iy の実部・虚部に同一に作用するため

    ‖Oz‖² = ‖z‖²          （ノルム保存）
    (Oz)ᵀ(Oz) = zᵀ OᵀO z = zᵀz   （零閉塞 zᵀz の厳密保存）

**(式15') 相対平衡（親）** — 親 v は K(arg v) の固有平面上にあり（第1章・式4〜式6）、
流れ dZ/dt = KZ の上では KZ = −iσZ、すなわち Z(t) = e^{−iσt} Z(0)。全辺の位相が同量だけ
進むので位相差は不変、したがって K(θ(t)) = K(θ(0))。親は「位相差を変えずに剛体回転する」
相対平衡である。

### 2.3 本プログラムの力学の数学（段1・段2・段3 の定義）

段番号は N=40 一因子実験系列（補遺A）の定義に従う。**段1が基準アーキテクチャ**を与え、
**段2・段3がそれに施す2つの変更**である。

**(式16) 段1: 基準アーキテクチャ（明示行列・スペクトル写像・固定 Δτ 時計）** —
生成子をエルミート行列 H として明示的に構成し、1ステップを厳密スペクトル指数写像

    H = V diag(w) V†,   z ← V e^{−i Δτ w} V† z

で与える。時計は**固定の外部時計**

    Δτ = 2π/den,   den ∈ { N−2, N−1, N, N+1, N+2 } ∩ ℕ⁺ ∪ { 124 }

（小 N では den>0 の条件により系列が特殊になる。N=3 の系列は den=1,2,3,4,5,124。
den=1 は Δτ=2π を意味する。）段1のみの基準では H = A∘(z̄⊗z)（振幅込み）を用いる。
本写像は無記憶（1ステップは z と定数のみの純関数）であり、旧の隠れ状態 wp・
近似σ推定は存在しない。

#### 2.3.1 分母系列 den の設計背景（仮説）と実装対応

Δτ = 2π/den の分母系列 { N−2, N−1, N, N+1, N+2, 124 } は本章の発明ではなく、
本シリーズの**分母コントロール実験**（2026-09-03、設計指示書
`ChatGPT_denominator_controls_N3_N40_mixedseed_20260903/CLAUDE_CODE_RUN_INSTRUCTION_N3_N40_20260903.md`）
の設計をそのまま継承したものである。その背景は次の2系統に分かれる。

**(a) den=124（固定時計の対照）** — 干渉保存力学系列の刻み規約 Δτ = 2π/L, **L=124**
（例: `干渉保存力学_資格審査とシード無し系列_20260831/program/pass1_parents.py` 27行
`L=124`）の継承。時計を N に依らず固定した legacy 規約の対照である。

**(b) den ∈ {N−2..N+2}（系スケール時計の仮説系列）** — 仮説は式20 の回転角
ψ_k = (2π/den)·σ_k を通じて実装と結びつく。第1章で保存した親のσスペクトル
（npz の `sigma`）から σ_max を実測すると（集計 `analysis_sweep_summary_v1.json` の
`sigma_max_by_N`）:

| N | σ_max | ψ_max/2π（den=N） | ψ_max/2π（den=124） |
|---|---|---|---|
| 3 | 1.414 | 0.471 | 0.011 |
| 5 | 3.742 | 0.748 | 0.030 |
| 10 | 8.928 | 0.893 | 0.072 |
| 20 | 18.894 | 0.945 | 0.152 |
| 30 | 28.898 | 0.963 | 0.233 |
| 40 | 38.905 | 0.973 | 0.314 |

σ_max(N) は N とともにほぼ線形に増える（N−σ_max ≈ 1.1）。したがって

- **den ≈ N の系列**では最速固有平面の1ステップ回転が ψ_max ≈ 2π·(1 − O(1)/N)、
  すなわち**ほぼ一回転（2π との可換に近い・ストロボ的）領域**に置かれる。
  den を N±1, N±2 とずらす5点は、この可換性からの離調を掃引する対照である。
- **den = 124** では全固有平面が小角領域（本スイープの範囲で ψ_max/2π ≤ 0.314）に
  収まり、連続流 dz/dt = Kz の小刻み近似に近い**流れ的領域**の対照になる。

**実装対応**: den は 38行の系列生成（`pairs=[(N+o,…) for o in OFFSETS if N+o>0]+[(124,'124')]`）
で決まり、力学には 27行の `2.0*math.pi/den` を通じてのみ入る。den はそれ以外の箇所で
使われない（式20 の通り、den の効果は回転角の線形スケールとエイリアシングに尽きる）。

本章はこの設計仮説の当否を論じない。den 依存の観測結果は §6（onset の分母依存・
未交差の分布）に事実として記載する。

**(式17) 段2: 振幅正規化** — 生成子構成への入力を z から単位振幅化した ẑ に置換する:

    ẑ_e = e^{i θ_e} = e^{i arg z_e}
    Ĥ_ef = A_ef · conj(ẑ_e) ẑ_f = A_ef · e^{i(θ_f − θ_e)}

これにより生成子は旧（式11）と同じく**位相のみ**の関数になる。状態 z 自身の振幅は
保持され力学変数のまま進化する（正規化は生成子構成の中だけで毎ステップ行われる）。

**(式18) Ĥ の実虚分解** — Ĥ はエルミートであり、実対称部と実反対称部に一意に分解される:

    Ĥ = S + iK,
    S_ef = A_ef cos(θ_f − θ_e)   （実対称）
    K_ef = A_ef sin(θ_f − θ_e)   （実反対称）

この K は**旧生成子（式11）と同一の行列**である。

**(式19) 段3: 虚部抽出（cos対称部の除去）と実直交回転** — 生成子として Ĥ の代わりに

    H₃ = i K = i · Im(Ĥ)

を用いる（i×実反対称 = エルミート）。段1の写像（式16）に代入すると

    z ← exp(−i Δτ · iK) z = exp(Δτ K) z

exp(ΔτK) は実反対称行列の指数、すなわち**実直交行列**である。ゆえに旧（式15）と同じ
保存量 ‖z‖², zᵀz を持つ（証明は式15と同一）。段2により K が位相のみで作られるため、
式15' の相対平衡の議論もそのまま成立し、親軌道上で K は時間不変、親は剛体回転する。

**(式20) 新の固有平面回転角** — H₃ = iK の固有値 w_k は実数で、K の ±iσ_k に対応して
w = ±σ_k。固有平面 k の1ステップ回転角は

    ψ_k = Δτ · σ_k = (2π/den) · σ_k

回転角は固有値に**線形**で、上限がない。位相因子 e^{−iΔτw} は w について周期 2π/Δτ で
巻き戻る（**エイリアシング**）: Δτσ_k が 2π を超える固有平面では、実効回転角は
Δτσ_k mod 2π となる。

### 2.4 新旧の回転写像の数学的差異（本章の力学のキモ）

両者は「位相のみ実反対称生成子 K による実直交回転」という点で同一だが、
**回転角のスペクトル写像と時計が数学的に異なる**。

| 項目 | 旧（式12〜式14） | 本章＝段1+2+3（式16〜式20） |
|---|---|---|
| 固有平面回転角 | φ_k = 2 arctan(γ σ_k/σ̂_max) | ψ_k = (2π/den) σ_k |
| スペクトル依存 | arctan による非線形圧縮（\|φ\| < π） | 線形（mod 2π のエイリアシングあり） |
| 時計 | σ̂_max 正規化（最速平面 π/72 固定、系固有） | 固定 Δτ = 2π/den（外部、den は掃引パラメータ） |
| σ の取得 | 冪反復3回の近似 σ̂_max（warm start・履歴依存） | eigh による厳密固有値（推定なし） |
| 写像の実現 | Cayley 有理式を Woodbury 解で適用（O(N³)/step） | スペクトル分解で厳密指数を適用（O(M³)/step） |
| 記憶 | 隠れ状態 wp を持つ（純関数でない） | 無記憶（z の純関数） |
| 保存量 | ‖z‖², zᵀz（実直交） | 同左（実直交） |
| 相対平衡（親） | 成立（式15'） | 成立（式19の帰結） |

**(式21) 小角極限での対応** — Δτσ_k ≪ 1 かつ γσ_k/σ̂_max ≪ 1 の領域では
arctan x ≈ x より

    φ_k ≈ (2γ/σ̂_max) · σ_k = Δτ_eff · σ_k,   Δτ_eff = 2 tan(π/144)/σ̂_max ≈ (π/72)/σ̂_max

すなわち旧は「Δτ_eff が σ̂_max で毎ステップ再規格化される式20」と近似的に一致する。
両者の差は (i) Δτ の値（外部固定 2π/den か、系固有の (π/72)/σ̂_max か）、
(ii) arctan 圧縮とエイリアシングの有無、(iii) σ の厳密性、に整理される。
この対応の実測は補遺Aの最終行（σ正規化時計の参照実験: 後半勾配 65.8 steps/decade、
7月正本 64.0 steps/decade）に記録されている。本章の走行はあくまで**固定 Δτ（段1）**で
行われ、分母 den の掃引がこの時計差の影響を測る。

### 2.5 測定量と保存則 — 休眠比 H⊥/H の定義とインフレーション指標としての根拠

**(式22) 読出し平面** — 初期状態 z0 = v + δg（正規化後）から

    p = Re z0 / ‖Re z0‖,   q' = Im z0 − (Im z0·p) p,   q = q'/‖q'‖

p, q は走行を通じて**固定**される（初期状態の実部・虚部が張る実2次元平面
Π = span_ℝ(p, q) のグラム・シュミット正規直交基底）。δ = 10⁻¹⁵ なので Π は実質的に
**親 v の張る平面**である。

**(式23) 測定量** — 各 step の状態 z について

    H⊥/H = ‖ z − p(p·z) − q(q·z) ‖² / ‖z‖²   （読出し平面 Π の外にあるエネルギー比）
    H_total = z†z
    closure = |zᵀz| / z†z

**(式23') H⊥/H が「休眠比」であり、親の運動を拾わないこと** — 純粋な親（δ=0）の軌道は
相対平衡（式15'・式19）であり、Z(t) = e^{−iσt} Z(0)。展開すると

    e^{−iσt}(x + iy) = (x cos σt + y sin σt) + i (y cos σt − x sin σt),   x=Re Z(0), y=Im Z(0)

すなわち実部・虚部は常に span_ℝ(x, y) = Π の中に留まる。補直交射影
z − p(p·z) − q(q·z) は Π 内の任意の複素結合を厳密に消すので、**親軌道上で H⊥ は恒等的に
0** である。したがって H⊥/H は「親の剛体回転**以外**のすべて」——初期には種 g の
Π 外成分（H⊥/H(0) ≈ δ² 〜 10⁻³⁰）、以後は力学がそこから育てた成分——だけを測る。
種は親のスケールからは見えないエネルギーであることから、7月正本はこれを
**休眠フラクション（dormant fraction）**と呼ぶ。

**(式23'') 7月正本との定義の同一性** — 本章の H⊥/H は7月正本の f(τ) と同一の測定である。
正本 `自発的分裂予備実験_v1/run_spontaneous_splitting_largeN_v1.py`:

```python
    57	        Zp = Z - p * (p @ Z) - q * (q @ Z)
    58	        htot = float(np.real(np.conj(Z) @ Z))
    59	        f = float(np.real(np.conj(Zp) @ Zp)) / htot
```

（p, q の構成も同一。46-48行。）ゆえに本章の曲線は7月のインフレーション図
（dormant_growth_large_n_v1.png）と同じ量を同じ定義で描いたものであり、直接比較できる。

**インフレーション指標としての根拠** — 以上より H⊥/H の時間発展は
「(i) 初期値が種スケール 10⁻³⁰ 台にあること（親が平衡である限り第1歩で跳ばない）、
(ii) そこから半対数プロット上の直線＝一定レートの指数増幅が続くこと、
(iii) O(10⁻²〜1) で飽和すること」の3点でインフレーション的発展を定量化する。
(i) が破れる（第1歩で 10⁻⁸〜10⁻³ へ跳ぶ）場合はミスマッチ注入であり、種の増幅とは
区別される（補遺Aの削除対照はまさにこの (i) で判別される）。onset（式24）は (iii) への
到達の指標である。

保存則: 式19より写像は実直交だから、closure の分子 |zᵀz| と分母 z†z は理論上ともに
不変であり、closure の時間変化は数値丸めの蓄積のみを測る（実測は §6・G3）。
また ‖z‖ 保存により H⊥/H ∈ [0,1] が保証される。

**(式24) onset（記録用指標）** — H⊥/H > 0.05 となる最初の step 番号（なければ −1）。
力学には影響しない集計値である。

### 2.6 複素平面読出しとの対応 — インフレーション図と3種の複素平面図の設計意図

第3章の複素平面図（step0・終了時・凝縮中心の拡大）は、本章のインフレーション図
（H⊥/H 曲線）と対をなす読出しである。両者の関係を定める恒等式を先にここで与える。

**注意（混同しやすい2つの「平面」）** — 図の複素平面は「各辺の値 z_e ∈ ℂ を M 点
打った平面」であり、H⊥/H の基準である親平面 Π = span_ℝ(p,q) は「状態空間 ℂ^M の中の
実2次元部分空間」である。両者は別物だが、次の恒等式で結ばれる。

**(式25) Π 内の状態の per-edge 表示** — x = Re z0, y = Im z0, v ≈ z0（δ=10⁻¹⁵）とすると
x = (v+v̄)/2, y = (v−v̄)/(2i) より、

    z ∈ span_ℂ{x, y} = Π  ⟺  z_e = a·v_e + b·conj(v_e)   （a, b ∈ ℂ は辺に依らない定数）

すなわち **Π 内に留まる運動は、図の上では「親の星型の一斉回転・スケール（a·v_e）と
その鏡像の混合（b·v̄_e）」の2複素パラメータ族の中の変形しか起こせない**。特に純粋な
親軌道 Z(t)=e^{−iσt}Z(0)（式15'）は a=e^{−iσt}, b=0 で、**星型全体が形を変えずに原点まわりを
一斉回転するだけ**である。逆に、この族で表せない形（例: 全方位への位相分散＝リング）が
図に現れることは、状態が Π の外へ出たことの配置的な表現である。

この恒等式により、3種の複素平面図はインフレーション図（H⊥/H 曲線）の
始点・終点・終点内部を配置側から裏づける役割分担を持つ:

| 図 | H⊥/H 曲線との対応 | 読み取る内容 |
|---|---|---|
| step0 図 | 曲線の始点（H⊥/H ≈ 10⁻³¹） | Π の中身＝親の形（対蹠2ペア・4束の星型）。step1 の図が step0 と見分けがつかないことは f(1) が種スケールに留まること（§2.5 (i)）の可視化 |
| 終了時図 | 曲線の飽和域（H⊥/H ~ 0.05〜0.46） | 式25の族では表せない形（星→リング）への逸脱＝ Π 外成分が支配的になった配置。onset（式24）は形の崩れが目視可能になる時期に概ね対応 |
| 凝縮中心の拡大図 | 曲線が言わない情報の補完 | H⊥/H は「どれだけ Π を出たか」しか測らない。拡大図は「出た先の組織化」——角クラスターが厳密な同一複素数への凝縮（厳密縮退）か有限幅の束か——を分解する |

図の実装（グリッド描画・クラスター抽出の算法・行番号対応）は第3章で記述する。

## 3. 実装方法

- スイープ本体 `run_N3_N40_stage123_v1.py` は自己完結（エンジン import なし）。
  同梱エンジン `run_n_scaling_lowrank_v1.py` は §2.2（旧の数学）の定義体・第1章の
  親生成にのみ使われ、本章の力学ループには関与しない。
- 本プログラムは ChatGPT 作成の正本 `run_and_plot_N3_N40_mixedseed_20260903.py` から
  **文書化された最小差分の鎖**（補遺B）だけを経て導出されており、力学関数 `one_step` は
  N=40 一因子実験で3ゲート合格済みの段3版と同一である。
- 初期データは第1章の静的親 npz の `Z0` を per-N で読み込む（再生成・変更なし）。
- 入力ゲート `check_sweep_inputs_v1.py` が、全228走行の保存済み `Z[0]` と静的親 `Z0` の
  bit 一致を事後検証する。
- §6 の数値は集計プログラム `analyze_sweep_summary_v1.py` の出力
  `results/analysis_sweep_summary_v1.json` から転載する（手計算・手集計なし）。

**段⇔数式⇔実装の対応表**（実装は `run_N3_N40_stage123_v1.py`）:

| 段 | 内容 | 数式 | 実装箇所 |
|---|---|---|---|
| 段1 | 明示行列＋厳密スペクトル写像 | 式16 | 21-22行（`H_of`）・27-28行（`eigh` と位相因子適用） |
| 段1 | 固定 Δτ=2π/den 時計と分母系列 | 式16 | 27行の `2.0*math.pi/den`・38行（系列生成） |
| 段2 | 振幅正規化（ẑ=e^{iθ}） | 式17 | 26行の `np.exp(1j*np.angle(z))` |
| 段3 | 虚部抽出 H₃=i·Im(Ĥ) → 実直交回転 | 式18・式19 | 26行の `(1j*np.imag(H))` |
| — | 読出し平面・測定量・onset | 式22・式23・式24 | 29-30行・31-32行・45-46行 |

（対比: 旧の式12〜式14 はエンジン 122-132 行・134-141 行にのみ存在し、本章の
ループでは呼ばれない。）

## 4. 詳細設計

### 4.1 全体フロー

```
for N in 3..40:
  [初期化]   静的親 Z0 読込 → 隣接行列 A（式10）→ 読出し平面 p,q（式22）→ 分母系列（式16）
  for den in 系列:
    [ループ処理] t=0..500: 状態と測定量（式23）を記録し、t<500 なら one_step（式16〜式19）
    [終了処理(den)] 状態 npz 保存、onset（式24）等を summary 行に追加
[終了処理(全体)] timeseries/summary CSV 出力 → 8×5 グリッド図 → RUN_METADATA 出力
```

### 4.2 全体データフロー

- **入力**: `parents/parent_static_N{N:05d}_makeparent_20260905.npz` の `Z0`（38個、第1章の
  成果物。SHA256 は SHA256SUMS.txt 正本）
- **パラメータ**:
  - `STEPS = 500` … 1走行のステップ数（早期停止なし、固定）
  - `OFFSETS = (-2,-1,0,1,2)` … 分母系列の N からのオフセット（式13）
  - 分母 124 … 全 N 共通の対照分母
  - dtype: 状態 complex128／実数 float64（12行目の assert で固定）
- **出力**:
  - `results/hm_N{N}_den_{den}_states_500.npz` × 228（`Z`(501×M), `N`, `denominator`, `steps`）
  - `results/timeseries_64bit_with124_N3_N40.csv`（列: N, series, denominator, step,
    Hperp_frac, H_total, global_closure。228×501 = 114,228 行）
  - `results/summary_64bit_with124_N3_N40.csv`（列: N, series, denominator, onset_gt_0.05,
    initial, step1, final, max。228行）
  - `results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`（目標図）
  - `results/RUN_METADATA_N3_N40_stage123.json`

### 4.3 個別処理

#### 4.3.1 初期化（`run_N3_N40_stage123_v1.py`）

定数と dtype 固定:

```python
    11	STEPS=500; OFFSETS=(-2,-1,0,1,2)
    12	assert np.dtype(np.float64).itemsize==8 and np.dtype(np.complex128).itemsize==16
```

辺・隣接行列（式1・式10）:

```python
    14	def edges(N):
    15	    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
    16	def adjacency(N):
    17	    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    18	    for e in range(M):
    19	        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    20	    return A
```

親読込・読出し平面（式22）・分母系列（式16）:

```python
    34	for N in range(3,41):
    35	    # 初期データ: 各 N の静的親ファイルの Z0 を使用
    36	    z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)
    37	    A=adjacency(N); p,q=plane(z0)
    38	    pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]
```

```python
    29	def plane(v):
    30	    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p); q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q
```

#### 4.3.2 ループ処理（1ステップの力学 — 式16〜式20）

```python
    21	def H_of(z,A):
    22	    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
    23	def one_step(z,A,den):
    24	    # 段3の最小変更（唯一の力学変更点）: 位相のみ生成子 Ĥ の虚部だけを取る H=i·K（K=sin(Δθ) 実反対称）。
    25	    # exp(-iΔτ·iK)=exp(Δτ·K) の実直交回転となり、Z^T Z（零閉塞）と ‖Z‖ を厳密保存する。
    26	    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    27	    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    28	    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
```

- 26行前半 `np.exp(1j*np.angle(z))`: **段2**＝式17（振幅正規化。生成子入力を ẑ=e^{iθ} に置換）。
- 26行後半 `(1j*np.imag(H))`: **段3**＝式18の虚部抽出→式19の生成子 H₃=i·K
  （`np.imag` は要素ごとの虚部。エルミート行列 Ĥ の虚部は実反対称行列 K になる）。
- 27-28行: **段1**＝式16（`np.linalg.eigh` による厳密スペクトル分解、固定 Δτ=2π/den の
  位相因子 e^{−i(2π/den)w} の適用）。固有平面回転角は式20 ψ_k=(2π/den)σ_k であり、
  旧の式14（arctan 圧縮・σ̂_max 正規化）とは異なる（§2.4）。

測定と記録（式23）:

```python
    31	def metrics(z,p,q):
    32	    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real; return float(hp/h),float(h),float(abs(z@z)/h)
```

```python
    40	        z=z0.copy(); vals=np.empty(STEPS+1,np.float64); states=np.empty((STEPS+1,z.size),np.complex128); closures=np.empty(STEPS+1,np.float64); htot=np.empty(STEPS+1,np.float64)
    41	        for t in range(STEPS+1):
    42	            states[t]=z; vals[t],htot[t],closures[t]=metrics(z,p,q)
    43	            if t<STEPS: z=one_step(z,A,den)
```

**停止条件・例外的挙動（数式化しない要素）**:
- 早期停止はない。全走行が**固定 501 回**（t=0..500）の記録と 500 回の写像適用で構成される
  （41-43行）。onset（式24）は集計にのみ使われ、力学を止めない（45-46行の
  `np.flatnonzero(vals>0.05)`）。
- 大きい N の走行では `V@(phase*(V†z))` の行列積に対して numpy の RuntimeWarning
  （divide by zero / overflow / invalid value encountered in matmul）が表示されることがある。
  これは警告であって例外ではなく、実行は停止しない。本系列では同一入力に対する再走行の
  bit 一致（第1章 G1、および本章の入力ゲート）で結果への影響がないことを確認している。
- N=3,4 では式16の `N+o>0` フィルタにより den=1（Δτ=2π）等の小さな分母が系列に入る。
  プログラム上の特別扱いはない（38行の内包表記の帰結）。式20の観点では小さな den は
  大きな Δτ を意味し、エイリアシング（Δτσ_k mod 2π）が強く働く領域である。

#### 4.3.3 終了処理（保存・図化・メタデータ）

```python
    44	        np.savez_compressed(os.path.join(OUT,f'hm_N{N}_den_{den}_states_500.npz'),Z=states,N=np.int64(N),denominator=np.int64(den),steps=np.int64(STEPS))
    45	        rows.extend((N,label,den,t,vals[t],htot[t],closures[t]) for t in range(STEPS+1)); ix=np.flatnonzero(vals>0.05)
    46	        summaries.append((N,label,den,int(ix[0]) if ix.size else -1,float(vals[0]),float(vals[1]),float(vals[-1]),float(vals.max())))
```

CSV・図・メタデータ（48-66行）は測定値の書き出しと描画のみで力学に影響しない。
図は 8×5 グリッド（55-64行）、38パネル使用・2パネル off（64行）。

#### 4.3.4 入力ゲート（`check_sweep_inputs_v1.py`）

```python
    18	for N in range(3, 41):
    19	    Z0 = np.load(os.path.join(PARENT_DIR, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0']
    20	    dens = [N + o for o in (-2, -1, 0, 1, 2) if N + o > 0] + [124]
    21	    for den in dens:
    22	        p = os.path.join(RESULT_DIR, f'hm_N{N}_den_{den}_states_500.npz')
    23	        same = bool(np.array_equal(np.load(p)['Z'][0], Z0))
    24	        n_checked += 1
    25	        if not same:
    26	            print(f'MISMATCH: N={N} den={den}')
    27	            ok = False
...
    35	sys.exit(0 if ok else 1)
```

保存済み全 npz の `Z[0]` が対応する静的親 `Z0` と `np.array_equal`（bit 一致）であることを
検証し、1件でも不一致なら終了コード 1。

## 5. 実行結果

### 5.1 再現コマンド

```bash
cd N3_N40_stage123_sweep_20260905
python3 run_N3_N40_stage123_v1.py        # スイープ本体
python3 check_sweep_inputs_v1.py         # 入力ゲート
python3 analyze_sweep_summary_v1.py      # §6 の数値の集計
# または ./run_all.sh（親生成→スイープ→ゲート→集計→図化）
```

### 5.2 実行環境

- Python 3.9.6（`.venv/bin/python3`）、numpy 2.0.2（BLAS/LAPACK: macOS Accelerate）、
  matplotlib（グリッド図描画）
- macOS 26.3.1（arm64）

### 5.3 実行時間

スイープ本体（38N×6分母×500 step、状態全保存込み）で**約45〜50分**（実測 2026-09-05、
支配項は M×M の `eigh` × 500 × 6。N=40 単独で約3分）。入力ゲートと集計は各1分未満。

### 5.4 検証ゲート

| ゲート | 合格条件 | 実測 | 合否 |
|---|---|---|---|
| G1: 入力同一性 | 全228 npz の `Z[0]` が対応する静的親 `Z0` と bit 一致 | checked 228 / MISMATCH 0 | **PASS** |
| G2: 完走 | `ALL DONE` 出力・終了コード 0 | 確認 | **PASS** |
| G3: 保存量（観察） | closure = \|zᵀz\|/H が全走行・全 step で小さいまま | step0 最大 5.19e-15、全域最大 3.72e-13 | **PASS** |

（G3 の数値は `analysis_sweep_summary_v1.json` の `global_closure_step0_max` /
`global_closure_all_max`）

### 5.5 データ

| 項目 | 内容 |
|---|---|
| フォルダ | `N3_N40_stage123_sweep_20260905/results/` |
| 状態 npz | `hm_N{N}_den_{den}_states_500.npz` × 228（数十KB〜約5.9MB/個、合計約900MB 相当のブロック、実サイズは SHA256SUMS.txt と同梱物参照） |
| 時系列 CSV | `timeseries_64bit_with124_N3_N40.csv`（約8.7MB、114,228行） |
| サマリ CSV | `summary_64bit_with124_N3_N40.csv`（228行） |
| メタデータ | `RUN_METADATA_N3_N40_stage123.json` |
| 集計 | `analysis_sweep_summary_v1.json`（§6 の数値の出所） |
| SHA256 | 全ファイルの正本値は同梱 `SHA256SUMS.txt`。本体プログラム: `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567  run_N3_N40_stage123_v1.py` |

### 5.6 図化

- 目標図: `results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`
  （8×5 グリッド。各パネルの**縦軸は休眠比 H⊥/H**〔式23・意味は §2.5〕の半対数、
  横軸は step。6分母の曲線を重ね描き。半対数上の直線区間がインフレーション的
  指数増幅、水平区間が飽和に対応する）
- 複素平面読出し図（step0・終了時・拡大）は第3章で扱う。

## 6. 実行分析（客観的報告と観察のみ。数値は analysis_sweep_summary_v1.json より）

1. 228走行のうち **212走行が 500 step 内に H⊥/H > 0.05 を交差**した（onset 範囲 45〜481）。
2. **全228走行で step1 は種スケールに留まった**: H⊥/H(0) ∈ [1.92e-33, 1.60e-31]、
   H⊥/H(1) ∈ [3.64e-30, 7.32e-25]。ミスマッチ注入（10⁻⁸〜10⁻³ 級の跳び）は 0 件である。
3. 交差した走行の終値は 0.0501〜0.458。
4. 未交差16走行の内訳: 15走行は final 0.025〜0.048 で**閾値 0.05 のわずかに下で飽和**して
   おり（N≥26 の 124 系列と N+2 系列、および N=39 の3系列）、成長が遅かったのは
   N=3・den=124 の1走行のみ（final 1.75e-8）。
5. 124 系列の onset は N=4 の 318 から N の増加とともに短くなり、N≥15 では 89〜160 の帯に
   入る（交差した N に限る。表は `analysis_sweep_summary_v1.json` の
   `onset_by_series['124']`）。
6. 零閉塞 |zᵀz|/H は全 228×501 記録点で最大 3.72e-13（初期最大 5.19e-15）に留まった。
   式12の保存則の数値的な確認である。
7. 図の目視観察: 全パネルで「種スケールからの直線的指数増幅 → 10⁻³〜10⁻¹ 台での飽和」
   という同型の曲線族が現れる。

## 補遺A 段構成の由来 — 一因子実験と削除対照の監査表

本章の力学（段1+2+3）は、N=40・静的親（第1章と同一の Z0）に対する一因子実験系列で
確定した構成である。実験一式（プログラム・データ・図・README・SHA256SUMS）は
`ChatGPT_denominator_controls_N40_selfcontrol_20260904/` に完備されており、
**本論文のアップロードに同梱する**。数値の出所は各 results ディレクトリの
summary CSV および README である。

| 構成 | 力学（生成子／時計） | f(1) の水準 | 0.05 交差（500step, 6分母） | 終値・観察 | 出力ディレクトリ |
|---|---|---|---|---|---|
| 段1のみ（基準） | H=A∘(z̄⊗z)（振幅込み・cos入り）／固定Δτ | 5.3〜6.5e-8（注入） | 0/6 | 天井 ~1.16e-3 で安定 | results_staticparent/ |
| 段1+2（段3削除） | 位相のみ Ĥ（cos入り）／固定Δτ | 1.4e-9〜8.9e-3（注入・分母依存） | 6/6（τ=4〜198） | 0.94〜0.9999 まで離脱 | results_staticparent_phaseonly/ |
| 段1+3（段2削除） | 振幅込み i·Im(H)／固定Δτ | 1.8〜2.3e-8（注入） | 0/6 | ~3.5e-3 へ緩慢に漂う | results_staticparent_ampimK/ |
| **段1+2+3（本章）** | 位相のみ i·K／固定Δτ | **4.0〜9.1e-29（種スケール）** | 4/6（τ=276〜456） | 0.038〜0.10、緩和曲線あり | results_staticparent_imK/ |
| 段2初期化のみ（参考） | z0 を1回等振幅化＋振幅込み i·Im(H)／固定Δτ | 2.9〜3.4e-3（注入） | 6/6（τ=3〜6） | 0.96〜0.99、即時離脱 | results_staticparent_stage2init/ |
| 段2+3＋σ時計（参考） | 位相のみ i·K／σ正規化 Cayley 時計 | 5.84e-30（種スケール） | 0/6（500step内） | 後半勾配 65.8 steps/decade（7月正本 64.0） | results_staticparent_sigmaclock/ |

観察（事実のみ）: f(1) が種スケール（~10⁻²⁹）に留まるのは段2と段3が**同時に**入っている
構成のみであり、いずれか一方を削除した構成では f(1) が 10⁻⁸〜10⁻³ に跳ぶ。段2を初期化時
1回に移した構成では跳びが最大（10⁻³）になる。

## 補遺B プログラム系譜（最小差分の鎖と SHA256）

各段の差分は該当パッケージの README に全文記録されている。

| # | プログラム | 直前からの差分 |
|---|---|---|
| 0 | `ChatGPT_denominator_controls_N3_N40_mixedseed_20260903/run_and_plot_N3_N40_mixedseed_20260903.py`（ChatGPT 作成正本） | — |
| 1 | `ChatGPT_denominator_controls_N40_selfcontrol_20260904/run_and_plot_N40_only_selfcontrol.py` | OUT 先変更・`range(40,41)` の2行（N=40 出力が正本と bit 一致を確認） |
| 2 | `…/run_N40_staticparent_v1.py` | 初期データを第1章静的親 Z0 に差替（入力のみ） |
| 3 | `…/run_N40_staticparent_phaseonly_v1.py` | one_step 1行: H を exp(i·arg z) から構成（段2） |
| 4 | `…/run_N40_staticparent_imK_v1.py` | one_step 1行: H=1j·Im(Ĥ)（段3） |
| 5 | `N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py`（本章） | ループ range(3,41)・per-N 親読込・出力先・図名・メタデータ名のみ（力学無変更） |

SHA256（全桁、shasum -a 256 の出力を転載）:

```
c709c56335d4c67373bff9a3ef6414ea17564d5fbf9c10b7bc9c3724ff091b92  run_and_plot_N3_N40_mixedseed_20260903.py
5a07f354e19985dca0f5de89217e2aa22ac511afaa4b2bb4aa5d93e0c7f9706f  run_and_plot_N40_only_selfcontrol.py
cb5a0ab6db9ae5719eac7aee3b539eb04ca09f5ebfcab6dcce36e8a6e727719e  run_N40_staticparent_v1.py
c1d6d2e60e101f0a99585eefdc22990a087c46246aba112042ac97a7fcbd1a71  run_N40_staticparent_phaseonly_v1.py
a67912b77f7f112731c1eac7612f21b464aed7e2c42431f7b3a7afabb2cd051d  run_N40_staticparent_imK_v1.py
1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567  run_N3_N40_stage123_v1.py
```

---
（第2章おわり。第3章「複素平面読出し図」に続く）
