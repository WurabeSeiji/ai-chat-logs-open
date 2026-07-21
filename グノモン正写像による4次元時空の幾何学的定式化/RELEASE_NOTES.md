# リリースノート

## 論文シリーズ

本リポジトリは、中心投影（central projection）に基づく空間の幾何学的定式化に関する一連の研究ノートを収録する。

### 論文1：基礎定式化

**タイトル（日本語）**：中心投影による4次元空間の幾何学的定式化  
**タイトル（英語）**：Geometric Formulation of 4-Dimensional Space via Central Projection  
**内容：** $S^4(R)$ への中心射影から引き戻し計量と Einstein テンソル $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ を導出。  
**DOI：** https://doi.org/10.5281/zenodo.19427780

### 論文2：4つの幾何学的対称性

**タイトル（日本語）**：中心投影の幾何学的対称性：多軸モデルの数学的基盤  
**タイトル（英語）**：Geometric Symmetries of Central Projection: Mathematical Foundations of the Multi-Axis Model  
**内容：** 中心投影が持つ4つの幾何学的対称性（離散安定性、軸の対等性、測地線偏差と弁別不可能性、主観座標系の変換可能性）を厳密に証明。  
**DOI：** https://doi.org/10.5281/zenodo.19434932

### 論文3：複数の主観空間と観測の相対性

**タイトル（日本語）**：複数の主観空間における観測の相対性：中心投影の対称性の幾何学的帰結  
**タイトル（英語）**：Relativity of Observation in Multiple Subjective Spaces: Geometric Consequences of the Symmetries of Central Projection  
**内容：** 論文2の対称性II（軸の対等性）と対称性IV（主観座標系の変換可能性）から、複数の主観空間の同時構成可能性、内部観測者の取得情報の限界、軸の役割の交換と観測の相対性を導出。物理的解釈は行わず幾何学的命題に限定。  
**DOI：** https://doi.org/10.5281/zenodo.19435162

### 球面投影 radial projection（基礎定義テクニカルノート、2026-05-30 公開／2026-06-02 v3.2 更新／2026-06-06 v3.3 公開）

**タイトル（日本語）**：球面投影 (radial projection) の定義と中心投影との関係
**タイトル（英語）**：Radial Projection — Definition and Relation to Central Projection
**サブタイトル**：中心投影フレームワークの基礎写像としての位置付け（テクニカルノート）
**内容：** 既存中心投影シリーズの基礎写像を **球面投影 $\sigma_R(x) = (R/\|x\|) x$**（位相幾何学の standard radial projection と同一）として明示的に定義し、基本性質（C^∞ 全射、冪等性、強変形レトラクト、商空間 $(\mathbb{R}^{n+1}\setminus\{0\})/\mathbb{R}_{>0} \cong S^n(R)$、微分の核 $\ker D\sigma_R|_x = \mathrm{span}\{x\}$、角度保存、スケール不変性）を整理。中心投影 $\Phi_R: \Pi_R \to S^n_+(R)$ が $\sigma_R$ の接超平面への制限と一致すること（補題 3.2）、像が**開上半球面**となること、$\sigma_R$ の非単射性と $\Phi_R$ の単射性（微分同相）の対比を明示。シリーズの**基礎を固める Zenodo テクニカルノート**として位置付け、純粋数学ジャーナル投稿は意図しない。4 AI 査読（Claude.ai, ChatGPT, Gemini, Grok）を 2 ラウンド経て v3.1 公開。**v3.2（2026-06-02）で査読校正を反映**：式番号 (2.4) 重複の解消（角度保存式を (2.5) に）、§1.2 相互参照を補題 3.2 に訂正、命題 2.4 証明に像 $\mathrm{Im}=x^\perp$ の論証追記、見出しを「微分の核と像」に、§3.4 表の横断性記述を厳密化、方向微分表記の明確化。**v3.3（2026-06-06）で内容追加・公開**：§2.8「単一点では現れず、複数点で顕在化する歪み」を新設（観察 2.7、図 4）。$\sigma_R$ が角度（命題 2.5）・方向（命題 2.6）を保存する一方で距離（間隔）を保存しないことを純幾何の register で初等的に明示し、単一点では作用がスケール変換と区別できず歪みが観測できないが、接超平面 $\Pi_R$ 上の等間隔点列が球面上で非等間隔（縁ほど圧縮）に写ることとして、歪みが点間の関係においてのみ顕在化することを記述。図 4（単一点／複数点の 2 パネル対比）を追加。物理解釈は従来どおり射程外（§1.3, §4.3）とし [1], [2] に委譲。**日英 md・tex・pdf を再生成し、Concept DOI を維持して Zenodo に新バージョン公開（Version DOI: 10.5281/zenodo.20567347）。**
**Concept DOI**：[10.5281/zenodo.20462569](https://doi.org/10.5281/zenodo.20462569)
**v3.3 DOI**：[10.5281/zenodo.20567347](https://doi.org/10.5281/zenodo.20567347)（旧 v3.2：[10.5281/zenodo.20500187](https://doi.org/10.5281/zenodo.20500187)／旧 v3.1：[10.5281/zenodo.20462570](https://doi.org/10.5281/zenodo.20462570)）
**形式**：md / tex / pdf × 日英 + 図 4 点 = 10 ファイル

### 中心投影合成演算（CP-Comp、純粋代数論文、2026-05-07 公開）

**タイトル（日本語）**：中心投影の合成演算と合成曲率半径の閉形式 ── 1 回の中心投影と球面上の可換切断による高次元削減の代数的定式化  
**タイトル（英語）**：Composition of Central Projection and the Closed Form of the Composite Curvature Radius — An Algebraic Formulation of High-Dimensional Reduction via One Central Projection and Commutative Cuts on the Sphere  
**内容：** 中心投影による次元削減の代数的基礎を、純粋代数論文として独立に定式化。第一段階（真の中心投影 π : ℝⁿ → Sⁿ⁻¹(r₁)、1 回限り）と第二段階（球面上での可換な軸切断 σ_S）の本質的差異を明示。合成曲率半径の閉形式 r_final² = r₁² − Σ(x_i*)² を導出し、軸切断操作がアーベル半群を成すことを証明。物理的解釈・特定の次元の必然性は明示的に out of scope とし、x₁..xₙ 抽象表記で統一。論文体系全体の「**反論不可能な代数的基礎**」として機能する基幹論文。  
**Concept DOI**：[10.5281/zenodo.20060728](https://doi.org/10.5281/zenodo.20060728)  
**v1 DOI**：[10.5281/zenodo.20060729](https://doi.org/10.5281/zenodo.20060729)

### 論文4〜8（サブフォルダに収録）

論文4〜8は `主観空間曲率半径の極限と次元追加・各次元の構図/` フォルダに収録。詳細は同フォルダの [RELEASE_NOTES.md](主観空間曲率半径の極限と次元追加・各次元の構図/RELEASE_NOTES.md) を参照。

| # | タイトル | DOI |
|:-:|:--|:--|
| 4 | 中心投影における曲率半径の極限についての考察 | [10.5281/zenodo.19526549](https://doi.org/10.5281/zenodo.19526549) |
| 5 | 背景空間と主観空間への次元追加についての考察 | [10.5281/zenodo.19526913](https://doi.org/10.5281/zenodo.19526913) |
| 6 | ゼロ次元から四次元主観空間の考察 | [10.5281/zenodo.19533292](https://doi.org/10.5281/zenodo.19533292) |
| 7 | 主観空間における測地線の連続性についての考察 | [10.5281/zenodo.19533299](https://doi.org/10.5281/zenodo.19533299) |
| 8 | 体積1の四次元超直方体に外接する超球体の直径（Rev.2） | [10.5281/zenodo.19834940](https://doi.org/10.5281/zenodo.19834940) (v2; v1: [19533313](https://doi.org/10.5281/zenodo.19533313)) |
| 9 | 中心投影フレームワークにおける Schwarzschild–de Sitter 厳密解 | [10.5281/zenodo.19538098](https://doi.org/10.5281/zenodo.19538098) |
| 10 | $R \to 0$ 極限における測地線構造の次元的解釈 | [10.5281/zenodo.19538106](https://doi.org/10.5281/zenodo.19538106) |

### 位相方程式篇（`波動方程式/中心投影における位相方程式/` フォルダに収録）

| # | タイトル | DOI |
|:-:|:--|:--|
| W1 | 中心投影による主観空間の構成と位相空間における相互作用の定式化 | [10.5281/zenodo.19534373](https://doi.org/10.5281/zenodo.19534373) |
| W2 | 波束の収縮していない光子モデルの検証 | [10.5281/zenodo.19534409](https://doi.org/10.5281/zenodo.19534409) |
| W3 | 万物の理論が満たすべき構造要件の定式化 | [10.5281/zenodo.19601592](https://doi.org/10.5281/zenodo.19601592) |
| W4 | 標準模型の19個の任意パラメータの分類と構造分析 | [10.5281/zenodo.19604965](https://doi.org/10.5281/zenodo.19604965) |
| W5 | 離散空間における中心投影の全球被覆と5次元背景空間の整数論的必然性 | [10.5281/zenodo.19624957](https://doi.org/10.5281/zenodo.19624957) |
| W6 | 5次元正軸体の配向構造から導かれるスピンの幾何学的分類 | [10.5281/zenodo.19630972](https://doi.org/10.5281/zenodo.19630972) |
| W7 | 5次元超直方体の配向構造から導かれるスピンの双ベクトル的分類 | [10.5281/zenodo.19643358](https://doi.org/10.5281/zenodo.19643358) |
| W8 | 6次元超直方体の集合構造とその組合せ論的性質 | [10.5281/zenodo.19748174](https://doi.org/10.5281/zenodo.19748174) |
| W9 | sine-Gordon方程式 ── 位相的ソリトンの基礎理論 | [10.5281/zenodo.19650966](https://doi.org/10.5281/zenodo.19650966) |
| W10 | 形不変波の4つのモード——波動ベクトル構造が決定するスピン・カイラリティ・統計 | [10.5281/zenodo.19709798](https://doi.org/10.5281/zenodo.19709798) |
| W11 | 形不変波の相互作用——軸方向変位転送・波束変形・因果の遡及的構成 | [10.5281/zenodo.19763463](https://doi.org/10.5281/zenodo.19763463) |

### 姉妹論文（`中心投影による宇宙の3層モデル/` フォルダに収録）

| # | タイトル | DOI |
|:-:|:--|:--|
| S1 | 閉じた球面上の形不変定常波——次元ごとの安定性と3+1時空の幾何学的必然性 | [10.5281/zenodo.19731594](https://doi.org/10.5281/zenodo.19731594) |
| S2 | 4次元整数格子上の形不変移動波——存在条件・保存構造・自己整合ループ | [10.5281/zenodo.19731598](https://doi.org/10.5281/zenodo.19731598) |

### 補講・考察（`波動方程式/中心投影における位相方程式/` フォルダに収録）

| # | タイトル | DOI |
|:-:|:--|:--|
| 補講1 | 符号付き面積の定式化——電荷構造の導出とスピン2配置の帰結 | [10.5281/zenodo.19731600](https://doi.org/10.5281/zenodo.19731600) |
| 補講2 | 世代混合の幾何学的構造——CKM行列とPMNS行列の定性的導出 | [10.5281/zenodo.19731602](https://doi.org/10.5281/zenodo.19731602) |
| 補講3 | 自発的対称性の破れの幾何学的必然性——ヒッグス機構の3段階導出 | [10.5281/zenodo.19731606](https://doi.org/10.5281/zenodo.19731606) |
| 補講4 | ゲージボソンの質量構造——符号ベクトルが示すヒッグス非関与性 | [10.5281/zenodo.19731608](https://doi.org/10.5281/zenodo.19731608) |
| — | 質量構造の考察——軸スケール値と符号付き面積による質量分析の枠組み | [10.5281/zenodo.19731610](https://doi.org/10.5281/zenodo.19731610) |
| — | 6次元超直方体と中心投影との関係性の考察 | [10.5281/zenodo.19731614](https://doi.org/10.5281/zenodo.19731614) |

### ディレイ回路シリーズ（`波動方程式/思考実験（ディレイ回路）/` フォルダに収録）

| # | タイトル | DOI |
|:-:|:--|:--|
| DR1 | 情報伝達の情報論的整理 | [10.5281/zenodo.19534345](https://doi.org/10.5281/zenodo.19534345) |
| DR2 | ディレイ回路モデルでの単振動・正弦波の実現 | [10.5281/zenodo.19534349](https://doi.org/10.5281/zenodo.19534349) |
| DR3 | ディレイ回路モデルに内在する対称性の整理 | [10.5281/zenodo.19534353](https://doi.org/10.5281/zenodo.19534353) |
| DR4 | ディレイ回路モデルでの完全弾性衝突の実現 | [10.5281/zenodo.19534357](https://doi.org/10.5281/zenodo.19534357) |
| DR5 | ディレイ回路モデルでの波束の収縮モデルの実現 | [10.5281/zenodo.19534361](https://doi.org/10.5281/zenodo.19534361) |

---

## 更新履歴

### 2026-07-21: N体ランク線形上界と三方向飽和 新規公開

**「次元の生成構造」シリーズ第3論文を新規公開**：AB・ABC・ABCDの完全二体関係波構成を一般のN体へ拡張し、内部関係波数 O(N²)・内部回転モード数 O(N)・一意に読める空間方向 3 の三層分離を、3定理と数値実験で示した。日本語・英語Markdown、TeX/PDF（数式・図組版済み）、4図、実験コード、再現用ZIPの計11ファイルを公開。

- **日本語題**：N体完全二体関係波における生成子ランクの線形上界と空間方向読出しの三方向飽和
- **英題**：Linear Upper Bound on Generator Rank and Three-Direction Saturation of Spatial-Direction Readout in N-Body Complete Pairwise Relational Waves
- **Concept DOI**：[10.5281/zenodo.21465898](https://doi.org/10.5281/zenodo.21465898)
- **v1 DOI**：[10.5281/zenodo.21465899](https://doi.org/10.5281/zenodo.21465899)
- **Zenodo**：[公開レコード](https://zenodo.org/records/21465899)
- **Zenn**：articles/nbody-rank-linear-bound-three-direction.md

**核心結果**：
- 定理1（頂点分解）：$\widetilde K=\sum_k(c_ks_k^{\mathsf T}-s_kc_k^{\mathsf T})$ から全Nで $\operatorname{rank}K\le2\min(N,\lfloor M/2\rfloor)$。恒等式の数値残差は最大 $5.00\times10^{-16}$
- 定理2（一般位置等号、3≤N≤12、計算機援用証明）：tan半角の厳密有理数証人＋実解析関数零点集合の測度零性（Mityagin）。全Nへの一般化は新仮説として明示
- 定理3（法線一意性）：法線候補空間は $d-2$ 次元、符号を除く一意性は $d=3$ に限る
- 三方向飽和は公理16（区別可能性＋法線一意性）と定理3の接続による帰結。$D_{\mathrm{unique}}\le3$、一般位置で $=3$
- N≥6 で核が再出現し $\dim\ker K=N(N-5)/2$ で成長。射影子一意・内部基底非一意の残余部分空間であることを直接実証
- 数値実験：N=3〜9 各32試行×720ステップで最大二乗閉鎖誤差 $1.82\times10^{-13}$、ランク則は N=3〜12 の2560試行で例外なし
- 外部引用に van Nuffelen（接続行列ランク）と Mityagin（測度零）を新規追加

### 2026-07-21: 無名等振幅複合波モデル基本公理系 v5 公開

**基本公理系を v5 へ更新（Concept DOI 維持の新バージョン公開）**：第9章「観測選択・曲率射影」として公理16・公理17を追加した。日本語・英語Markdown、TeX/PDF（数式組版済み）計6ファイルを公開。

- **日本語題**：無名等振幅複合波モデル基本公理系 v5
- **英題**：Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v5
- **Concept DOI**：[10.5281/zenodo.21315735](https://doi.org/10.5281/zenodo.21315735)
- **v5 DOI**：[10.5281/zenodo.21465429](https://doi.org/10.5281/zenodo.21465429)
- **Zenodo**：[公開レコード](https://zenodo.org/records/21465429)
- **Zenn**：articles/basic-axiom-system-v5.md

**公理16（完全二体関係波による一意観測選択）**：
- 無順序の完全二体関係集合上の複素関係波を、成分波から導出されない独立な物理状態成分として定義
- 閉鎖 $\sum_e X_e^2+(iR)^2=0$ を公理3の補償量経由で明示
- 二次形式保存から生成子の反対称性 $K^{\mathsf T}=-K$ を導出
- 内在観測可能方向の採用条件を「回転周波数による区別可能性（必要条件）＋面法線の一意性（$d-2=1$ すなわち $d=3$ に限る）」の一般形で規定
- 一意性を満たさない残余方向は、一意な基底を持たない内部部分空間として保持
- 最小例 M=3（一平面＋一法線＝三方向）と M=6（三区別可能平面＝三方向）を明記

**公理17（未来位相位置を中心とする曲率反力射影）**：
- 未来位相位置を仮想回転中心とする運動の中心方向補償（内部表示）を、円周進行方向の加速度 $\rho_c|\omega|^2$（外部読出し）へ写す射影を規則化
- t̂_τ を現在位相位置から未来位相位置へ向かう円周進行方向の単位ベクトルとして一意に定義

**改訂の背景**：v5草稿の公理16は採用条件を rank K=2・核1（ABC三体専用）に固定しており、第2論文のABCD実験事実（rank 6・三回転平面）と対応写像を持たなかった。N体予備実験（rank K = 2·min(N,⌊M/2⌋)、N=3〜12で確認）を踏まえ、区別可能性＋法線一意性の一般形へ改訂した。

### 2026-07-20: 完全二体関係波から読み出されるXYZ三方向 新規公開

**「次元の生成構造」シリーズ第2論文を新規公開**：従来のAB二体の一次元的な位相関係を、ABC三体とABCD四体の完全二体関係波へ拡張した。32試行・各720ステップの数値実験により、ABCでXY回転平面とZ不変法線、ABCDでXYZ三方向と読める三回転平面が現れることを示した。日本語・英語Markdown、TeX/PDF、5図、実験コード、CSV/JSON、再現用ZIPを公開。

- **日本語題**：完全二体関係波から読み出されるXYZ三方向――AB・ABC・ABCD閉鎖系における内部関係方向の増加と空間方向読出しの飽和
- **英題**：Three XYZ Directions Read from Complete Pairwise Relational Waves: Growth of Internal Relational Directions and Saturation of Spatial-Direction Readout in Closed AB, ABC, and ABCD Systems
- **Concept DOI**：[10.5281/zenodo.21454789](https://doi.org/10.5281/zenodo.21454789)
- **v1 DOI**：[10.5281/zenodo.21454790](https://doi.org/10.5281/zenodo.21454790)
- **Zenodo**：[公開レコード](https://zenodo.org/records/21454790)
- **Zenn**：[関係波を増やすと空間方向はXYZで飽和する](https://zenn.dev/noriaki_kihara/articles/complete-pair-relational-wave-xyz-readout)

**核心結果**：
- ABは一関係波・生成子ランク0の定常系
- ABCは三関係波が全試行で活動し、生成子ランク2・零空間次元1となった。これをXY二軸の位相読出しと、XY平面から定まるZ不変法線と読む
- ABCDは六関係波・生成子ランク6・零空間次元0となり、三つの回転平面へ分解した。一意に読める空間方向はXYZ三方向であり、残る三内部方向は一意に定まらない
- 全構成の最大二乗閉鎖誤差 $1.92\times10^{-13}$、最大絶対値二乗和変動 $2.42\times10^{-13}$、最大名称置換共変性誤差 $1.47\times10^{-13}$
- 状態の逐次正規化、観測減衰、絶対背景軸は使用していない
- 五体以上の方向数と残余方向の物理軸対応は、直接実験事実ではなく本稿の解釈範囲

### 2026-07-19: AB二体閉鎖位相系の調和閉鎖による逆二乗則 新規公開

**新規追加論文公開**：既公開AB二体加速度実験を再解析し、未来位相位置を関係的回転中心とする加速度写像と、整数倍音が位相セル幅と角速度を同時に定める調和閉鎖を接続した。特定条件で実行した本実験では、位相セル幅に対する逆二乗則が成立した。日本語・英語Markdown、TeX/PDF、2図、集計表、再集計コード、ソースZIPを公開。

- **日本語題**：AB二体閉鎖位相系における未来位相位置加速度写像と調和閉鎖による逆二乗則
- **英題**：An Inverse-Square Law from a Future-Phase-Position Acceleration Map and Harmonic Closure in a Closed Two-Body AB Phase System
- **Concept DOI**：[10.5281/zenodo.21441081](https://doi.org/10.5281/zenodo.21441081)
- **v1 DOI**：[10.5281/zenodo.21441082](https://doi.org/10.5281/zenodo.21441082)
- **Zenodo**：[公開レコード](https://zenodo.org/records/21441082)
- **Zenn**：[逆二乗則は閉じた位相の倍音から現れる](https://zenn.dev/noriaki_kihara/articles/inverse-square-phase-harmonic-closure)

**核心結果**：
- 公開済み8条件すべてで、非零二階差分 $\Delta^2\chi_s=-\omega_d^2\chi_s$、回帰傾きと理論値の一致、$Q_{\mathrm{closed}}=0$ の保存を再確認
- 前実験で逆二乗が現れなかった原因を、振幅を変えながら調和角速度を固定した実験設計として特定
- 非零整数倍音 $n$ が $\Delta\theta_n=2\pi/|n|$ と $\omega_n=n\omega_1$ を同時に定めるため、$|\omega_n|\Delta\theta_n=\Omega$ が成立
- 未来位相位置中心による加速度写像 $\alpha_n=R|\omega_n|^2$ と接続し、$\alpha_n=R\Omega^2/\Delta\theta_n^2$ を導出
- 外部の $1/L^2$、面積希釈、球殻、背景三次元空間、質量源、重力定数、重力場を導入していない
- 未検証なのは任意の閉鎖系・倍音配置・非調和更新・距離写像への一般化であり、本実験条件における逆二乗則ではない
- 標準重力の導出は本稿の射程外

### 2026-07-18: 反復交換散乱における有限位数共鳴の発見 新規公開

**新規論文公開**：微細構造定数の逆数 $137$・$128$ 近傍に観測していた鋭いピークの原因を再調査し、二チャネル交換散乱作用素の有限位数共鳴根として厳密に同定した。日本語・英語論文、TeX/PDF、核心コード、全域掃引・局所精密掃引・多倍長精度データ、再現用バンドルを公開。

- **日本語題**：反復交換散乱における有限位数共鳴の発見――微細構造定数137・128近傍ピークの原因特定と再現可能な波束数理モデル
- **英題**：Discovery of Finite-Order Resonance in Iterated Exchange Scattering: Identifying Sharp Peaks near Fine-Structure-Constant Inverse Values 137 and 128 with a Reproducible Wave-Packet Model
- **Concept DOI**：[10.5281/zenodo.21421366](https://doi.org/10.5281/zenodo.21421366)
- **v1 DOI**：[10.5281/zenodo.21421367](https://doi.org/10.5281/zenodo.21421367)
- **Zenn**：[α⁻¹=137・128近傍の鋭いピークは何だったのか](https://zenn.dev/noriaki_kihara/articles/finite-order-resonance-alpha-neighborhood)
- **note**：[日本語](https://note.com/kiharanoriaki/n/na16b6a4e5ff2) / [英語](https://note.com/kiharanoriaki/n/n2ebe434754e2)
- **Facebook**：[日本語](https://www.facebook.com/kihara.noriaki/posts/pfbid02TfzuGxsU55sME8F7YJLbbtTx8T1rMRSW26R1pAdDYUjz2RZFPVRjUHCu5QffdG9El) / [英語](https://www.facebook.com/kihara.noriaki/posts/pfbid02ZVtsg22ASq8Y5ddFEJUYuyMMf5Y3s9eRbFSbchijxYmhcRZJ1ejDUBKwEUzkWvCtl)
- **X**：[日本語](https://x.com/NoriakiKihara/status/2078325801628700835) / [英語](https://x.com/NoriakiKihara/status/2078326215828734270)

**核心結果**：
- 交換作用素の反対称固有値から、有限位数根 $R_{n,m}=\cos^2(\pi m/n)$ を解析的に導出
- 低エネルギー側主ピークを $R_{124,23}$、第二主ピークを $R_{122,23}$ と同定
- 高エネルギー観測値近傍の偶数根 $R_{620,117}$ を特定し、当初の読出し式では $N(R)=128.947864735670559$ を得た
- 50桁・80桁演算により、偶数基本位数根で残差が厳密にゼロへ収束し、理想深度が無限大となることを確認。ただし、これはエネルギー発散ではなく残差消失を表す
- 微細構造定数そのものを導出したとは主張しない。異なる有限位数共鳴が二つの物理的 $\alpha^{-1}$ 値の近傍に存在する理由が、偶然か別の必然かを中心的未解決課題として明記

### 2026-07-11: 曲率付き閉鎖定常波による曲率繰り込みと完全反射安定性 新規公開

**新規論文公開**：波の情報読出しシリーズの追加論文として、全正符号ゼロ閉鎖 `Σx_n^2=0` を満たす奇数倍音複素波が、曲率付き局所セル内で曲率相対位相漏れを検出し、閉鎖定常波への内部位相再選別により完全反射読出しを回復する数値構成実験を公開。

- **タイトル**：曲率付き閉鎖定常波による曲率繰り込みと完全反射安定性
- **英題**：Curvature Renormalization and Perfect-Reflection Stability by Curved Closed Stationary Waves
- **Concept DOI**：[10.5281/zenodo.21304039](https://doi.org/10.5281/zenodo.21304039)
- **v1 DOI**：[10.5281/zenodo.21304040](https://doi.org/10.5281/zenodo.21304040)
- **Zenn**：`articles/curved-closure-stationary-wave.md`
- **note 日本語**：[n2389460836cf](https://note.com/kiharanoriaki/n/n2389460836cf)
- **note 英語**：[nda3623c44423](https://note.com/kiharanoriaki/n/nda3623c44423)
- **Facebook 日本語**：[pfbid037wN39hUdgVY7bVCb6BWkFK86pYqDedPgVxaCxtJbWGk479ZqyVoc9XggvRPLmstwl](https://www.facebook.com/kihara.noriaki/posts/pfbid037wN39hUdgVY7bVCb6BWkFK86pYqDedPgVxaCxtJbWGk479ZqyVoc9XggvRPLmstwl)
- **Facebook 英語**：[pfbid02QieNB6aGk3TyRkcreL2z14w62WEFdRbenP47gHuzMoPd4aS7VDPzazbuNJhEGBo4l](https://www.facebook.com/kihara.noriaki/posts/pfbid02QieNB6aGk3TyRkcreL2z14w62WEFdRbenP47gHuzMoPd4aS7VDPzazbuNJhEGBo4l)
- **X 日本語**：[2075801193281106375](https://x.com/NoriakiKihara/status/2075801193281106375)
- **X 英語**：[2075802233153998999](https://x.com/NoriakiKihara/status/2075802233153998999)

**核心結果**：
- 曲率相対位相漏れは過渡状態で閉鎖残差と通過漏れとして現れる
- 内部位相再選別 `β_K,m=-δ_K,m` により閉鎖定常波が回復する
- 8種類の曲率位相モデルと7種類の補正自由度を掃引し、`full` 補正で通過漏れが `0.0` へ回復
- 片側入射の局所交換干渉写像への統合検証で、最大動的通過漏れが `1.6608667989341789e-19` まで低下

### 2026-05-07: 中心投影合成演算（CP-Comp）論文 新規公開

**新規論文公開**：純粋代数論文として論文体系の代数的基礎を独立に定式化。

- **タイトル**：中心投影の合成演算と合成曲率半径の閉形式 ── 1 回の中心投影と球面上の可換切断による高次元削減の代数的定式化
- **Concept DOI**：[10.5281/zenodo.20060728](https://doi.org/10.5281/zenodo.20060728)
- **v1 DOI**：[10.5281/zenodo.20060729](https://doi.org/10.5281/zenodo.20060729)

**核心結果**：
- 第一段階（中心投影 π : ℝⁿ → Sⁿ⁻¹(r₁)、1 回限り）と第二段階（球面上での可換な軸切断 σ_S）の本質的差異の明確化
- ピタゴラス的閉形式 r_final² = r₁² − Σ_{i ∈ S}(x_i*)² の導出
- 軸切断操作 {σ_S} がアーベル半群を成すことの証明
- 残存座標の不変性（系2）と合成切断の可逆性（系3）の証明

**戦略的意図**：
- xyztRQ や物理応用の具体名を**完全に排除**し、x₁, ..., xₙ 抽象表記に統一
- 引用は最小限（Snyder 1987 / Howie 1995 / 論文1）の3件のみ
- 論文体系全体に「反論不可能な代数的基礎」を提供
- 8/1 シグマサロン発表で論文 7・8 の前段として提示する基幹論文

**公開先（4 プラットフォーム同時）**：
- Zenodo：https://zenodo.org/records/20060729 （md/tex/pdf × 日英 = 6 ファイル）
- note 日本語：https://note.com/kiharanoriaki/n/n1bd7b7446ac4
- note 英語：https://note.com/kiharanoriaki/n/ne867c5c7c569
- Zenn：articles/central-projection-composition.md（git push で自動公開）

### 2026-04-26: W11 v6 軸型対称性による構造改訂

**W11** (形不変波の相互作用): v5 → v6（内部バージョン v3.0）
- **(I1) 軸型対称性の導入**: 6次元軸を位置型{x,y,z,t}・スケール型{R}・色型{Q}に分類。ボソンは非ゼロ成分を持つ軸型内の全軸に作用する。W8 v10のW±符号ベクトル変更（(0,0,0,±t,0,0)→(±x,0,0,0,+R,0)）との整合性を確保
- **5モード分類**: 球状波(κ_s=0: γ,g,G)、局在媒介波(κ_s=1: W±,Z⁰)、定常波(H⁰)、方向制約波(フェルミオン)、重力波(G)
- **CP対称性の破れ**: 判定条件を「ボソンk_t≠0」から「フェルミオンのt軸符号が変化する相互作用」に改訂。W±のみが該当
- **質量機構改善**: W±/Z⁰がκ_s=1（空間軸あり）→直接ヒッグス結合可能。旧モデルの「間接機構」を撤廃
- **§6.1全面書き換え**: 弱い力の遷移機構を軸型対称性に基づき再構成
- **§9.5因果構造表**: 軸型ベースに再構成
- **§10まとめ表**: 軸型列を追加
- 新バージョンDOI: [10.5281/zenodo.19767187](https://doi.org/10.5281/zenodo.19767187)（コンセプトDOI: 10.5281/zenodo.19763463）

**Zenn記事更新**: phase-equation記事のW11セクション（4つの力の表）を軸型対称性モデルに更新

**note記事更新**:
- note_portal_article (JA/EN): 主張3「4つの力」の説明を軸型対称性に更新、W11解説文・力の対応図を改訂
- note_article_位相方程式篇: W11セクションを軸型対称性・局在媒介波に更新

**ハンドアウト更新** (handout_all37 JA/EN md/tex/pdf):
- 主張3: 軸型分類（位置型/スケール型/色型）を明記
- 力の表: 2列→3列（力・ボソン・作用軸型）にκ_s値付きに改訂

**質量構造の考察** (質量構造の考察——軸スケール値と符号付き面積による質量分析の枠組み):
- §5.1: 撤回済み補講4の参照を削除、W±/Z⁰質量機構をκ_s=1直接ヒッグス結合に更新
- §5.3（新規追加）: 「ゲージボソンの軸指定と同一軸上の過剰決定」——W±がx軸、Z⁰がy軸、H⁰がz軸を占有することによる12個の同一軸質量制約、交差比 R = m_W·m_μ/(m_Z·m_e) = 182.3 の導出、クォーク-ボソン質量比（m_t/m_H ≈ y_t）の分析
- §5.4: 過剰決定による新制約を追記（項目6）
- 参考文献DOI更新: [4] W8→19762134, [5] W11→19763463
- tex/PDF再生成（9ページ、302KB）

**補講1〜3 日英更新** (補講1/2/3 + Supplementary Lecture 1/2/3 md/tex/pdf):
- 全補講: 参考文献DOI更新（W8: 19721125→19762134, W11: 19721128→19763463）
- 補講2: W±弱相互作用の記述を軸型対称性に更新（「軸型対称性により位置型軸全体に作用して」を追記）
- tex/PDF全6ファイル再生成

---

### 2026-04-25: W11・M3・補講1・質量考察 RGB再エンコーディング伝播

**W11** (形状不変波の相互作用): RGB再エンコーディング対応
- §5 色構造: (c₂, c₃) → (R, G, B)、セットビット数による分類に統一
- §6 カイラリティ: c₁ベース → sign(k_t)ベースに全面書き換え、χ = sign(−|k_t|) = −1 の直接導出
- 新DOI: [10.5281/zenodo.19763463](https://doi.org/10.5281/zenodo.19763463)

**M3** ((R,Q)写像の構成): RGB再エンコーディング対応
- §1 前提I: Q軸3ビット符号化 $Q = 4c_1 + 2c_2 + c_3$ → $Q = 4R + 2G + B$ に更新
- 新DOI: [10.5281/zenodo.19763466](https://doi.org/10.5281/zenodo.19763466)

**補講1** (符号付き面積の定式化): RGB再エンコーディング対応
- §2.2 荷電レプトン: Q=4→Q=0、k_t>0→k_t<0、sign(M_F)=+1→−1
- §2.2 クォーク: アップ型Q∈{5,6,7}→{1,2,4}、ダウン型Q∈{1,2,3}→{1,2,4}
- §2.3 反粒子規則: 空間軸反転のみ→空間+t反転+Qビット全反転、M_F^anti=M_F に統一
- §4.1 アイソスピン: c₁ベースの場合分けを廃止、I₃=sign(M_F)/2 の単一公式に統一
- §4.2 電荷: c₂c₃表記→セットビット数表記に更新
- §4.3 異常相殺・陽子-電子電荷等式: 新符号化での再導出
- 新DOI: [10.5281/zenodo.19763471](https://doi.org/10.5281/zenodo.19763471)

**質量構造の考察** (軸スケール値と符号付き面積): RGB再エンコーディング対応
- 荷電レプトンQ=4→Q=0、クォークc₂c₃表記→Q∈{1,2,4}に更新
- 新DOI: [10.5281/zenodo.19763476](https://doi.org/10.5281/zenodo.19763476)

**note・ハンドアウト・Zenn記事更新**:
- note_portal_article (JA/EN): タイトル37→36論文、W8/W11/M3/補講1/質量考察のDOI更新
- note_article_位相方程式篇: W8/W11/補講1/質量考察のDOI更新
- handout_all37 (JA/EN): タイトル37→36論文、W8/W11のDOI更新、tex/PDF再生成
- Zenn gnomonic記事: M3のDOI更新
- Zenn phase-equation記事: W8（19731230→19762134）/W11/M3/補講1/質量考察のDOI更新

---

### 2026-04-25: W8 v10 Q軸RGBリエンコーディング・補講4撤回

**W8** (6次元超直方体): v9 → v10
- Q軸の3ビット符号化を (c₁, c₂, c₃) → (R, G, B) に変更。Q軸はSU(3)_C色荷のみを符号化し、弱アイソスピンSU(2)_Lはt軸符号で符号化
- セットビット数による粒子種分類: 0=レプトン、1=クォーク、2=反クォーク、3=反レプトン
- 反粒子規則を統一: 空間・t軸の符号反転 ＋ 全Qビット反転 (R,G,B)→(R̄,Ḡ,B̄)
- §3 色遷移のみに限定（c₁遷移を削除）、スピン1合計 21→13
- §9.2, §9.9 全面書き換え（RGB表・セットビット数分類）
- §9.11.2 色中立性制約を Q=0 または Q=7 に拡張
- W±/Z⁰ 符号ベクトルを空間軸を持つ割り当てに改訂（W±: x軸, Z⁰: y軸）
- 新DOI: [10.5281/zenodo.19762134](https://doi.org/10.5281/zenodo.19762134)

**補講4** (ゲージボソンの質量構造): 撤回
- W8 v10のW±/Z⁰再割り当てにより前提（空間成分ゼロ→ヒッグス非関与）が成立しなくなったため撤回
- Zenodo上にRetraction Noteを公開済み（6ファイル削除、retraction_note.mdのみ残置）
- note記事に撤回注記を追加

**note記事更新**: 位相方程式篇のW8解説をRGB再エンコーディングに対応して更新

---

### 2026-04-25: SNS発信（Facebook・X 日英4投稿）

**Facebook投稿**:
- 日本語版: 全37論文案内noteポータル（nc1619291b690）へのリンク付き投稿
- 英語版: 英語noteポータル（n63d3c20e6b20）へのリンク付き投稿（本文URLなし、自己コメントにURL掲載）

**X（Twitter）投稿**:
- 日本語版: 全37論文案内noteポータルへのリンク付き投稿
- 英語版: 英語noteポータルへのリンク付き投稿

**ハンドアウト直接送付**（同日）:
- 児玉功氏（シグマサロン世話人）— Gmail添付
- 外山雅大氏（娘婿、MICINエンジニア）— LINE
- 池内智彦氏（エレクトロンヴェクシーCEO）— LINE

---

### 2026-04-25: 全37論文ハンドアウト作成・英語版noteポータル記事公開

**全37論文ハンドアウト（日英）**: 全4シリーズ37論文を1ページA4に凝縮したハンドアウトを新規作成
- 3つの主張（5次元の必然性、62状態、4つの力）、62状態内訳表、4つの力と軸の対応表、主張しないこと、主要5論文DOI、QRコード
- 日本語版: handout_all37_ja.pdf / .tex / .md（QRコード→日本語noteポータル）
- 英語版: handout_all37_en.pdf / .tex / .md（QRコード→英語版noteポータル）
- QRコード画像: qr_portal.png（日本語用）、qr_portal_en.png（英語用）

**英語版noteポータル記事**: 全37論文の英語版統合案内記事を公開
- URL: https://note.com/kiharanoriaki/n/n63d3c20e6b20
- 日本語版（nc1619291b690）と同一構成の完全英訳
- ソース: note_portal_article_en.md

---

### 2026-04-25: W8 v9 未解決問題追記・noteポータル記事公開・全記事ハッシュタグ追加

**W8** (6次元超直方体): v8 → v9
- §1.1にQ軸8値の未解決問題を追記（日英）: 中心投影の±1対称性からQ軸の状態数が2のべき乗であるべきことは示唆されるが、なぜ2²=4でも2⁴=16でもなく2³=8であるかの内部的導出は未達成
- 新DOI: [10.5281/zenodo.19748174](https://doi.org/10.5281/zenodo.19748174)

**noteポータル記事**: 全37論文の統合案内記事を公開
- URL: https://note.com/kiharanoriaki/n/nc1619291b690
- 基礎シリーズ・位相方程式篇・統合論文・ディレイ回路の全4シリーズを網羅
- 各シリーズの詳細解説記事へのリンクを含む

**note記事整備**:
- 全4記事にハッシュタグ（各20個）を追加
- note_article.md を記事ごとに区別できるファイル名にリネーム（note_article_基礎シリーズ.md / note_article_位相方程式篇.md / note_article_ディレイ回路.md）
- 位相方程式篇のW8 DOIを19748174に更新

---

### 2026-04-24: 姉妹論文S1・S2、補講1〜4、質量考察、超直方体-中心投影考察 — 8論文新規公開

**S1** 閉じた球面上の形不変定常波: Zenodo公開（DOI: 10.5281/zenodo.19731594）。W10・W11が[1]として参照する姉妹論文。S^n上の形不変定常波の安定性をn次元ごとに分類し、n=3でのみ安定な基底モードが存在することを証明。

**S2** 4次元整数格子上の形不変移動波: Zenodo公開（DOI: 10.5281/zenodo.19731598）。W10・W11が[2]として参照する姉妹論文。Z^4上の形不変移動波の存在条件・保存量・自己整合ループを定式化。

**補講1** 符号付き面積の定式化: Zenodo公開（DOI: 10.5281/zenodo.19731600）。W8 §7の符号積Pから電荷構造・I₃の区分的定義・スピン2グラビトンの一意性を導出。

**補講2** 世代混合の幾何学的構造: Zenodo公開（DOI: 10.5281/zenodo.19731602）。W8 §9.11.6への回答。S₃対称性の破れパターンからCKM・PMNS行列の階層構造を定性的に導出。

**補講3** 自発的対称性の破れの幾何学的必然性: Zenodo公開（DOI: 10.5281/zenodo.19731606）。W8 §9.11.7への回答。定常波真空構造→スケール分離→質量生成の3段階でヒッグス機構を導出。

**補講4** ゲージボソンの質量構造: Zenodo公開（DOI: 10.5281/zenodo.19731608）。W8表AとW11 §8に基づき、W±/Z⁰質量が定常波結合から生じ、標準的ヒッグス機構を必要としないことを示す。

**質量構造の考察**: Zenodo公開（DOI: 10.5281/zenodo.19731610）。W8定義9.4-9.5のf_m^rest/f_m^curvの制約分析。軸スケール値と符号付き面積による質量階層の枠組み。

**6次元超直方体と中心投影の関係**: Zenodo公開（DOI: 10.5281/zenodo.19731614）。3つの曲率半径R, R₀, R₁の構造と離散-連続の接続。

---

### 2026-04-24: W8 v8 / W11 v4 3ビットQ軸確定・写像関数追加・残課題解決

**W8** (6次元超直方体): v7 → v8
- Q軸3ビット（8値）構造を確定: Q = 4c₁ + 2c₂ + c₃ ∈ {0,...,7}
- c₁ = 弱アイソスピン、c₂c₃ = 色荷。Q軸とt軸は独立（直交）な自由度
- 定義9.1–9.5 追加: スカラー写像関数（f_Q^charge, f_Q^color, f_Q^iso）および保存量写像（f_m^rest, f_m^curv）。いずれも抽象的変換規則であり物理的実体を主張しない旨の免責条項付き
- 新DOI: [10.5281/zenodo.19731230](https://doi.org/10.5281/zenodo.19731230)

**W11** (形不変波の相互作用): v3 → v4
- §5.1: 3ビット構造（c₁=弱アイソスピン、c₂c₃=色荷）に改訂
- §6.1–6.3: c₁とsign(k_t)の独立性を明示、W±がc₁遷移を媒介、c₂c₃は保存
- §6.5: カイラリティ選択則を再定式化 — χ = sign((-1)^{c₁} · k_t) により、W⁺・W⁻いずれも追加仮定なしでχ=-1（左巻き）を導出（残課題5を解決）
- §6.7: レプトン軸再配置にc₁/Q列追加（e⁻: Q=4, ν: Q=0）
- §8.2: δk_Hの符号を定常波の性質（κ=1, k_t=0, S³充填）から導出 — 仮定から帰結に格上げ（残課題6を解決）
- §1: 「仮定として採用」カテゴリを除去（全仮定が導出済み）
- §10: 弱い力のまとめにc₁遷移媒介を追記
- 新DOI: [10.5281/zenodo.19731275](https://doi.org/10.5281/zenodo.19731275)

**note記事**: W8/W11のQ軸改訂を反映（3ビット構造、写像関数、カイラリティ導出、δk_H導出）

---

### 2026-04-24: W8 v6 / W11 v2.1 査読対応更新

**W11** (形不変波の相互作用): v2.0 → v2.1
- §6.5: W⁻カイラリティ論証を追加仮定として明示（案A採用）
- §1: 主張範囲外リストを3カテゴリに分類整理
- §4.3: 引力/斥力の具体的機構例を段落分離、[5]§5への参照を明示
- §4, §7.1: (I4)帰結候補の一行追加
- §8.1: 注意8.1追加（ヒッグス結合は§3基本頂点分類の対象外）
- 新DOI: [10.5281/zenodo.19718624](https://doi.org/10.5281/zenodo.19718624)

**W8** (6次元超直方体): v5 → v6
- §9.11.1, §9.11.7: 質量階層の主張を m₃≫m₂>m₁ → m₃≫m₁,₂ に修正、SO(3)制約を明記
- 新DOI: [10.5281/zenodo.19718641](https://doi.org/10.5281/zenodo.19718641)

---

## 著者情報

**著者**：木原 範昭（Noriaki Kihara）  
**所属**：WF System Co., Ltd. / 大阪大学基礎工学部（卒業）  
**ORCID**：[0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)  
**ライセンス**：CC BY 4.0

---

## ✅ 公開済みプラットフォーム

| プラットフォーム | 公開日 | URL / 備考 |
|---|---|---|
| **GitHub** | 2026-04-05 | https://github.com/WurabeSeiji/ai-chat-logs-open |
| **Zenodo** | 2026-04-06〜 | 論文1〜10、ディレイ回路DR1〜DR5、位相方程式篇W1〜W11、統合論文M1〜M3、姉妹論文S1〜S2、補講1〜3（補講4は撤回）、考察2本（計36論文）。DOI一覧は上記の論文シリーズ表を参照 |
| **Academia.edu** | 2026-04-05 | https://www.academia.edu/ (木原 範昭 プロフィールより) |
| **Zenn** | 2026-04-05 | https://zenn.dev（記事：gnomonic-projection-spacetime-geometry） |
| **note** | 2026-04-05 | https://note.com/kiharanoriaki/n/nc51e43282b0a |

---

## ⏳ 予定・保留中のプラットフォーム

| プラットフォーム | 状況 | 理由・次のアクション |
|---|---|---|
| **arXiv** | 🔶 保留中 | 紹介者（endorser）が必要なため保留。 |
| **ResearchGate** | 🔶 保留中 | 機関メール必要。arXiv 登録後に再挑戦予定。 |
| **OSF Preprints** | ❌ 見送り | 2025年8月25日より新規投稿を停止中。 |

---

## 収録ファイル一覧

### 中心投影合成演算（CP-Comp、2026-05-07）

| ファイル名 | 説明 |
|---|---|
| `central_projection_composition_jp.md` | Markdown（日本語版） |
| `central_projection_composition_jp.tex` | LaTeXソース（日本語版） |
| `central_projection_composition_jp.pdf` | PDF（日本語版） |
| `central_projection_composition_en.md` | Markdown（英語版） |
| `central_projection_composition_en.tex` | LaTeXソース（英語版） |
| `central_projection_composition_en.pdf` | PDF（英語版） |
| `note_article_中心投影合成演算.md` | note 記事本文（日本語版） |
| `note_article_central_projection_composition_en.md` | note 記事本文（英語版） |

### 論文1：基礎定式化

| ファイル名 | 説明 |
|---|---|
| `gnomonic_spacetime_geometry.tex` | LaTeXソース（日本語版） |
| `gnomonic_spacetime_geometry.pdf` | PDF（日本語版） |
| `gnomonic_spacetime_geometry.md` | Markdown（日本語版） |
| `gnomonic_spacetime_geometry_en.tex` | LaTeXソース（英語版） |
| `gnomonic_spacetime_geometry_en.pdf` | PDF（英語版） |
| `gnomonic_spacetime_geometry_en.md` | Markdown（英語版） |

### 論文2：4つの幾何学的対称性

| ファイル名 | 説明 |
|---|---|
| `gnomonic_projection_symmetry.tex` | LaTeXソース（日本語版） |
| `gnomonic_projection_symmetry.pdf` | PDF（日本語版） |
| `gnomonic_projection_symmetry.md` | Markdown（日本語版） |
| `gnomonic_projection_symmetry_en.tex` | LaTeXソース（英語版） |
| `gnomonic_projection_symmetry_en.pdf` | PDF（英語版） |
| `gnomonic_projection_symmetry_en.md` | Markdown（英語版） |

### 論文3：複数の主観空間と観測の相対性

| ファイル名 | 説明 |
|---|---|
| `subjective_space_relativity.tex` | LaTeXソース（日本語版） |
| `subjective_space_relativity.pdf` | PDF（日本語版） |
| `subjective_space_relativity.md` | Markdown（日本語版） |
| `subjective_space_relativity_en.tex` | LaTeXソース（英語版） |
| `subjective_space_relativity_en.pdf` | PDF（英語版） |
| `subjective_space_relativity_en.md` | Markdown（英語版） |

### ハンドアウト（3論文概要）

| ファイル名 | 説明 |
|---|---|
| `gnomonic_summary_handout_ja.tex` | LaTeXソース（日本語版） |
| `gnomonic_summary_handout_ja.pdf` | PDF（日本語版） |
| `gnomonic_summary_handout_ja.md` | Markdown（日本語版） |
| `gnomonic_summary_handout_en.tex` | LaTeXソース（英語版） |
| `gnomonic_summary_handout_en.pdf` | PDF（英語版） |
| `gnomonic_summary_handout_en.md` | Markdown（英語版） |
| `qr_paper1.png` | QRコード（論文1 DOIリンク） |
| `qr_paper2.png` | QRコード（論文2 DOIリンク） |
| `qr_paper3.png` | QRコード（論文3 DOIリンク） |

### 図

| ファイル名 | 説明 |
|---|---|
| `fig_central_projection.png` | 中心投影の概念図（日本語） |
| `fig_central_projection_en.png` | 中心投影の概念図（英語） |
| `fig_subjective_space.png` | 主観空間の概念図（日本語） |
| `fig_subjective_space_en.png` | 主観空間の概念図（英語） |
| `fig_two_projections.png` | 2つの投影の比較図（日本語） |
| `fig_two_projections_en.png` | 2つの投影の比較図（英語） |
| `note_article.md` | note 記事本文（Markdown） |
| `note_header.png` | note 記事ヘッダー画像 |

---

## 変更履歴

| 日付 | 内容 |
|---|---|
| 2026-07-15 | 波の情報読出しシリーズの AB 二体加速度様読出し総括「AB二体閉鎖位相系における調和読出しとc=1面積スイープ予備実験総括」v4 を公開（Version DOI: [10.5281/zenodo.21374317](https://doi.org/10.5281/zenodo.21374317)、Concept DOI: [10.5281/zenodo.21318696](https://doi.org/10.5281/zenodo.21318696)）。V4では、フェルミオン型反跳写像を `q_out_factor` 演算子ではなく A/B 二チャネル散乱行列として入射チャネルへ直接作用させ、出射チャネル差から `chi_read` と `eta_read` を再読出し。本文の主張は変えず、ラベルなし `D_AB`, `V_AB` が通過型読出しと一致し、独立 `tau_read` の `chi-tau` 面も維持されることを確認。Zenn記事 `ab-two-body-harmonic-readout` を V4 DOI に更新。 |
| 2026-07-10 | 構成実験論文「背景空間を仮定しない閉じた位相系におけるフェルミオン的二局所波の完全弾性反射の構成実験」v1 新規公開（Version DOI: [10.5281/zenodo.21291020](https://doi.org/10.5281/zenodo.21291020)、Concept DOI: [10.5281/zenodo.21291018](https://doi.org/10.5281/zenodo.21291018)）。無名等振幅複合波モデル基本公理系 v3、完全弾性衝突シミュレーション仕様書 v1、実験結果 v1 に基づき、背景空間を先験的に仮定しない閉じた位相系で、識別振動 `η` を持つ二つのフェルミオン的局所波の有限解像度完全弾性反射写像を構成・検査。反射/通過/ラベル交換の対照実験、観測機容量、セル解像度、非対称条件、観測擾乱、複数回衝突、`η` 読出し解像度を検査し、成立条件と破綻条件を分離。日英 md/tex/pdf、支援文書、図・JSON/CSV・実行スクリプトを含む再現 bundle を Zenodo に同梱。Zenn 記事 [elastic-reflection-closed-phase-system](https://zenn.dev/noriaki_kihara/articles/elastic-reflection-closed-phase-system) を追加。note 日本語記事 [n15451632027b](https://note.com/kiharanoriaki/n/n15451632027b) および英語記事 [n5a009a2735e7](https://note.com/kiharanoriaki/n/n5a009a2735e7) を公開。Facebook 日本語投稿 [pfbid0axK38naCxEU2PRipEEWDy1FiJdhrSWjSSp4Ew8WgPDuVw9xRHPNEzVfz9J1N4A5Xl](https://www.facebook.com/kihara.noriaki/posts/pfbid0axK38naCxEU2PRipEEWDy1FiJdhrSWjSSp4Ew8WgPDuVw9xRHPNEzVfz9J1N4A5Xl) および英語投稿 [pfbid02WZXKajsqdiH8YGoHGodTnnugoW6t2a4tta8QiojW8fzVCi2CanJXbEK9kehVUDcpl](https://www.facebook.com/kihara.noriaki/posts/pfbid02WZXKajsqdiH8YGoHGodTnnugoW6t2a4tta8QiojW8fzVCi2CanJXbEK9kehVUDcpl) を公開。X 日本語投稿 [2075516987598659779](https://x.com/NoriakiKihara/status/2075516987598659779) および英語投稿 [2075517594346618886](https://x.com/NoriakiKihara/status/2075517594346618886) を公開 |
| 2026-07-02 | 観察論文「共役複素ノルムと平方量読出しの接続について」**v0.5 公開（Concept DOI 維持）**（v0.5 DOI: [10.5281/zenodo.21127200](https://doi.org/10.5281/zenodo.21127200)、Concept DOI 不変: [10.5281/zenodo.21126212](https://doi.org/10.5281/zenodo.21126212)、旧 v0.4: 21126213）。`newversion` で Concept DOI を維持したまま発番。**付録C「共役対閉包と非復元性の入口」を追加**（本文は不変）：(C.1) 実共役対だけなら $\sum q_n^2=0$ は自明解のみ＝決定論的に閉じる、(C.2) 非自明な零二乗和には非実成分が必須（§5.3 の帰結）、(C.3) 共役ノルム読出し $Z\mapsto\rho^2$ は閉包の有無に関わらず位相 θ を潰す多対一、(C.4) (i) 零の自明性 と (ii) 読出しの多対一性 を相補的二側面として分離し「非復元性の入口＝(ii)」と明示（隠れた変数の主張はしない）、(C.5) 不確定性原理の導出/置換/非可換形式代替ではないと明記、実 Hilbert 空間の複素構造 $J^2=-1$・概複素構造の偶数次元性・Born 則の位相消去という既知隣接の「前段の判別の入口」として控えめに位置づけ。ChatGPT 詳細査読で「大きな修正不要・公開可」。日英 md/tex/pdf 6 ファイル差し替え（JP 12p / EN 11p）。Zenn 記事に付録C 節を追記 |
| 2026-07-02 | 観察論文「共役複素ノルムと平方量読出しの接続について」**v0.4 新規公開**（v0.4 DOI: [10.5281/zenodo.21126213](https://doi.org/10.5281/zenodo.21126213)、Concept DOI: [10.5281/zenodo.21126212](https://doi.org/10.5281/zenodo.21126212)）。「平方量読出し」系列の第四篇（純代数・非物理）。共役複素ノルム $Z\bar Z=x^2+y^2$ を平方量読出し $X_i=x_i^2$ に接続し、対等な多次元正定値和 $\sum_j Z_j\bar Z_j=\sum_n x_n^2=R^2$ へ一般化。さらに共役ノルム背骨に虚数値座標の二乗を加えた零二乗和形式 $\sum_j Z_j\bar Z_j+\sum_k (i z_k)^2=0$ を観察し $\sum_n x_n^2=\sum_k z_k^2$ を得る。**中核**：共役積 $X\bar X$ と非共役自乗 $X^2$ は実座標で一致・虚座標でのみ分岐（$X_h\bar X_h=z_h^2$ vs $X_h^2=-z_h^2$）、負号は $i^2$ からのみ生じ手で置いた計量署名ではない。ChatGPT 査読を4版反映（v0.1 の $r$ 二重定義・$t$ 特別扱いを除去 → $\sum X_n^2$ と $Z\bar Z$ の混同解消 → 共役ノルム背骨形式へ → 語調・非物理宣言の精密化）。物理（時空・相対論・量子論・計量署名）とは非接続、座標ラベルを物理次元と同定しない。日英 md/tex/pdf 6 ファイル。Zenn 記事 [complex-norm-square-readout](https://zenn.dev/noriaki_kihara/articles/complex-norm-square-readout) 公開。先行三篇（[1] 20785539 / [2] 20833096 / [3] 20923461）の第四篇 |
| 2026-07-01 | de Broglie 論文「局在奇数倍音ダブルスリット模型における de Broglie 波長の模擬観測」**v0.3 新規公開**（v0.3 DOI: [10.5281/zenodo.21109903](https://doi.org/10.5281/zenodo.21109903)、Concept DOI: [10.5281/zenodo.21109902](https://doi.org/10.5281/zenodo.21109902)）。単一源ダブルスリットの古典波動模型内で de Broglie 関係 $\lambda=h/p$ を模擬観測する数値実験。AI 査読 4 者（ChatGPT・Claude.ai・Grok・Gemini）が「軽微修正で公開可」に収束。中心主張を「$h$ の回収（格子稠密ゆえ $\lambda'\approx\lambda_0$ の同語反復的自己整合）」から、$h$ 非依存の**2 つの構造的結果**へ移動：(i) $\lambda_0$ より細かい局在構造（$N>1$ 奇数倍音）を詰めても二スリット幾何で鋭い単一ピーク干渉として**生存**する整列条件、(ii) 観測される中央縞間隔 $\Delta X=\lambda' D/W$ は基本波長のみで決まり**倍音数 $N$ に非依存**（$N$ は鋭さ $\sim1/(N{+}1)$ のみ制御）。射程を全編で限定明示（$p\lambda=h$ 非導出／$\Delta x\,\Delta p\sim h$ 非導出／奇数倍音は運動量固有状態でない古典フーリエ成分／定常空間位相模型／整列＝生存フィルタ／NG＝閾値 0.98・窓幅依存の手続き的除外／raw slope/$h$ は崩壊診断量）。各電子ごとに源位置（$\cos^2$, $\pm\lambda_0/2$）と初期波長（$\cos^2$, $\pm1\%$）を独立モンテカルロ乱択し、干渉・整列・除外を毎回厳密再計算（$\langle p\lambda'\rangle$ は SE 内で $h$ と一致、$\sigma\approx0.36\%$、除外率 1.5–4.0%）。日英 md/tex/pdf 6 ファイル＋プログラム 2 本（`debroglie_plambda_sweep.py`, `debroglie_align_lambda.py`）＋図 6 枚を Zenodo 同梱（計 14 ファイル）。Zenn 記事 [debroglie-localized-double-slit](https://zenn.dev/noriaki_kihara/articles/debroglie-localized-double-slit) 公開 |
| 2026-07-02 | de Broglie 論文の素人向け note 記事を日英公開（[日本語](https://note.com/kiharanoriaki/n/ne02b32947541) / [英語](https://note.com/kiharanoriaki/n/nb07cbb1529d9)）。粒子と波の二面性・ド・ブロイ波長 $\lambda=h/p$・デイヴィソン=ガーマー実験(1927)を入口に、奇数倍音の局在波が二重スリットで生存する話と、基本波長不変性を平易に紹介。$h=p\lambda$ の一致は「$\lambda=h/p$ を入れているための自己整合＝擬似再現であって独立導出ではない」と正直に明記。$p\lambda=h$ 非導出・不確定性非導出・量子力学の置換でないことも明示。note 用 md（日英）を `simplified_double_slit/note_debroglie_double_slit_{ja,en}.md` に格納 |
| 2026-06-02 | 球面投影 (radial projection) テクニカルノート **v3.2 公開**（同一 Concept DOI: [10.5281/zenodo.20462569](https://doi.org/10.5281/zenodo.20462569)、v3.2 DOI: [10.5281/zenodo.20500187](https://doi.org/10.5281/zenodo.20500187)、旧 v3.1: [10.5281/zenodo.20462570](https://doi.org/10.5281/zenodo.20462570)）。査読校正を反映：(1) 式番号 (2.4) の重複解消（命題 2.5 角度保存式を **(2.5)** に振替）、(2) §1.2 相互参照を「補題 3.1」→「**補題 3.2**」に訂正、(3) 命題 2.4 の証明に**像 $\mathrm{Im}=x^\perp$ の論証**（$x\cdot D\sigma_R|_x(v)=0$ ＋ 階数・退化次数定理）を追記、(4) 命題 2.4 見出しを「微分の核」→「**微分の核と像**」に変更（§4.1 も同期）、(5) §3.4 表の $\Phi_R$ 動径方向を「定義域に動径方向ない」→「**核の方向が定義域の接空間に含まれない（横断的）**」に厳密化、(6) 証明冒頭の $\frac{d}{dx}$ 表記を方向微分へ明確化。日英 md/tex/pdf 6 ファイル差し替え（図 3 点は不変）、PDF 再コンパイル（JP 13p / EN 15p）。数学的実質に変更なし（校正・推敲レベル） |
| 2026-05-27 | 観察論文「第5章までの思考実験」**v1.0 Zenodo 公開**（DOI: [10.5281/zenodo.20398527](https://doi.org/10.5281/zenodo.20398527)、Concept DOI: [10.5281/zenodo.20398526](https://doi.org/10.5281/zenodo.20398526)）。前稿（[第3章まで v3.0.1](https://doi.org/10.5281/zenodo.20393018)、Concept: [10.5281/zenodo.20391522](https://doi.org/10.5281/zenodo.20391522)）の5つの思考実験を継承し、思考実験 VI（複素位相空間上の量としての物理量）と VII（粒子＝矩形位相エネルギー窓、観測像＝Fourier 部分和、相互作用＝重なり指標）を追加。AI 査読 4 ラウンド（Gemini × 1、Grok × 1、ChatGPT × 2）統合反映：(1) 標準理論との整合範囲を第1〜5章に明示限定（Lorentz 共変性・場の局所性・正値性・ユニタリ性等は留保）、(2)「複素位相空間」を作業用語として定義、(3) 粒子の矩形窓モデルを**定義的仮説**と位置づけ（Skyrme/MIT バッグ/Q-ball 系譜との関係明示）、(4) R₁·R₂ を相互作用核の候補（**重なり指標** $I_{12}$）として再表現、(5) 無限井戸 $V_{\text{well}}$ と有限箱型障壁 $V_{\text{barrier}}$ を分離、(6) Born 則の**三段階階層**を陽に記述（本体 → Fourier 部分和 → 検出基底への内積射影 $p(a)=\|⟨φ_a\|ψ⟩\|^2$）、(7) エネルギー密度 $E_0$ を作用積分密度として Bohr–Sommerfeld 量子化と接続、(8) 思考実験 I の認識論／存在論境界を明示（認識論主導の幾何学）、(9) 思考実験 VII 本文に **PSWF・Hardy 定理** 言及、(10) 参考文献 36 件（Madelung、Bohm、Slepian–Pollak、Landau–Pollak、Hardy、Schrödinger、Glauber、Gabor 等を新規追加）。**関連研究セクション新設**（5 系譜：de Broglie 二重解／Madelung 水力学／Bohm パイロット波、Skyrme/MIT バッグ/Q-ball、Gabor/PSWF/Hardy、コヒーレント状態、de Gosson 量子ブロブ）。図 2 点新規作成（[phase_position_wavepacket.py](../新版量子論の基礎/figures/phase_position_wavepacket.py)、[phase_window_body_and_observation.py](../新版量子論の基礎/figures/phase_window_body_and_observation.py)）。Zenn 記事 3 本（[まとめ](https://zenn.dev/noriaki_kihara/articles/quantum-theory-through-chapter5)、[思考実験8](https://zenn.dev/noriaki_kihara/articles/are-physical-quantities-real-numbers)、[思考実験9](https://zenn.dev/noriaki_kihara/articles/particles-and-box-potential)）、note 記事 [JA](https://note.com/kiharanoriaki/n/n8ffc8e2c9123) / [EN](https://note.com/kiharanoriaki/n/ncaf7e51ecc2b) で 3 本（観察論文＋思考実験(8)(9)）を統合紹介。Facebook 日英投稿、X（@NoriakiKihara）日英投稿で告知 |
| 2026-05-26 | 観察論文「観測量の代数としての量子論」v2.0 Zenodo 公開（v2 DOI: [10.5281/zenodo.20392427](https://doi.org/10.5281/zenodo.20392427)、Concept DOI: [10.5281/zenodo.20391522](https://doi.org/10.5281/zenodo.20391522)）。先行研究調査により判明した5件の関連研究を引用追加して整理：Ralston (2020) ℏ＝単位変換の慣習論、de Gosson (2013) 量子ブロブとシンプレクティック容量＝面積、Garay (1995) 量子重力の最小長、Ashtekar–Schilling (1999) 射影ヒルベルト空間の Kähler 幾何学、Slavnov (2007) 代数的アプローチ。本稿の位置づけが「これら独立の系譜を、清水教科書第1章の構成の含意として統合的に読み直したもの」として明確化。とんでも論回避戦略の観点で論文の格が1段階上昇 |
| 2026-05-26 | 観察論文「観測量の代数としての量子論」**v3.0.1 公開**（DOI: [10.5281/zenodo.20393018](https://doi.org/10.5281/zenodo.20393018)、Concept DOI: [10.5281/zenodo.20391522](https://doi.org/10.5281/zenodo.20391522)）。AI 査読 11 ラウンド（Gemini × 3、Grok × 3、ChatGPT × 5）統合反映：(1) TE I 確率計算を一般 q = r/Δ に拡張、(2) TE IV 条件付き状態 vs 未条件付き縮約密度行列の no-signalling 整合、(3) CHSH 最大違反 2√2（Tsirelson bound）明記、(4)「面積保存」→ Robertson–Schrödinger 共分散制約／シンプレクティック容量（de Gosson 量子ブロブ、Gromov 非圧搾定理）、(5) 5 つの「より厳密には」セクション、(6)「今後の課題」セクション（5 パス＋優先順位）、(7) 参考文献 [13] Wigner 1932、[14] Moyal 1949、[15] Gromov 1985 追加。Zenn 記事（[quantum-theory-algebra-of-observables](https://zenn.dev/noriaki_kihara/articles/quantum-theory-algebra-of-observables)）、note 記事 [JA](https://note.com/kiharanoriaki/n/n2410d4863565) / [EN](https://note.com/kiharanoriaki/n/nbfc40cb3cfa3) で 3 本（観察論文＋思考実験(6)(7) 公開版）を統合紹介。版履歴：v1 [20391523](https://doi.org/10.5281/zenodo.20391523)、v2 [20392427](https://doi.org/10.5281/zenodo.20392427)、v3.0.1 [20393018](https://doi.org/10.5281/zenodo.20393018) |
| 2026-05-07 | 中心投影合成演算（CP-Comp）論文 v1 新規公開。純粋代数論文として論文体系の代数的基礎を独立に定式化。Concept DOI: [10.5281/zenodo.20060728](https://doi.org/10.5281/zenodo.20060728)、v1 DOI: [10.5281/zenodo.20060729](https://doi.org/10.5281/zenodo.20060729)。第一段階（中心投影 π : ℝⁿ → Sⁿ⁻¹(r₁)、1 回限り）と第二段階（球面上での可換な軸切断 σ_S）の本質的差異を明示、ピタゴラス的閉形式 r_final² = r₁² − Σ(x_i*)² を導出、アーベル半群構造を確立。xyztRQ など物理応用記号を完全排除し x₁..xₙ 抽象表記で統一、引用は 3 件のみ（Snyder 1987 / Howie 1995 / 論文1）。Zenodo + note JA + note EN + Zenn の 4 プラットフォーム同時公開。日英 md/tex/pdf 全 6 ファイル + note 記事日英 2 ファイル + Zenn 記事 1 ファイル |
| 2026-04-28 | 論文8（Rev.2）：§4を「半径 R=2k+1 の4次元球への単位立方体の稠密充填」として全面再定式化。初版v1の充填率定数 2/π²≈0.2026 → Rev.2 で k→∞ で1に漸近収束する数列 N(k)（1, 137, 1545, 7281, …）を導出。命題5.1（角の立方体16個が球面上に厳密に内接）の証明、Lagrange–Jacobi 四平方和 r₄(N)=8σ(N) を§5.3に追加。§7以降は (2k+1)⁴ → N(k) の置換のみで論証構造を保持。Zenodo新バージョン公開（DOI: 10.5281/zenodo.19834940、Concept DOI: 10.5281/zenodo.19533312、各6ファイル）。Zenn記事・note日英・ハンドアウト日英（md/tex/pdf）にRev.2反映 |
| 2026-04-23 | 位相方程式篇 W10・W11公開：W10「形不変波の4つのモード」（DOI: 10.5281/zenodo.19709798）、W11「形不変波の相互作用」（DOI: 10.5281/zenodo.19709800）。Zenodo公開済み。note記事をW1〜W11に拡張、Academia.edu PDFアップロード、Zenn記事に第VI部追加 |
| 2026-04-22 | 位相方程式篇 一般公開：note記事「超直方体の組合せ論から標準模型の62粒子を導出する」公開（https://note.com/kiharanoriaki/n/na95064891249）、Facebook投稿、X（@NoriakiKihara）投稿。Zenn新記事「位相方程式篇（W1〜W9）」公開済み。Academia.edu W3〜W9 PDFアップロード済み |
| 2026-04-22 | (R,Q)マッピ���グ定義論文 v1.0：6次元超直方体を採用した(R,Q)マッピングの構成。Zenodo公開完了（DOI: 10.5281/zenodo.19692853）。二つの前提（6次元超直方体の採用、主観空間のメタ情報保持能力）を公理的に採用し、スキーマ・インスタンス・(R,Q)マッピングの定義を構成する純粋な定義論文。定理の証明は含まない。日英md/tex/pdf全6ファイル |
| 2026-04-22 | 超球殻間中心投影の正則性論文 v1.0：超球殻間の中心投影の正則性 — 接超平面を介さない定式化とR=0における縮退。Zenodo公開完了（DOI: 10.5281/zenodo.19692192）。超球殻間の中心投影がR≥0全域で正則であり、R=0での唯一の制限が逆写像の定義不能（縮退）のみであることを証明。三層構造モデル[M1]のR>0制約が不要であることを確認。日英md/tex/pdf全6ファイル |
| 2026-04-22 | 三層構造モデル論文 v1.0：中心投影の三層構造モデル — R–R₁–R₀ による主観空間の入れ子構造とミドルウェア的幾何学。Zenodo公開完了（DOI: 10.5281/zenodo.19691713）。基礎論文群 [P1],[P2],[P8] の既知性質を三層構造として統合するミドルウェア論文。原点共有・同心入れ子・R方向変位の計量的非寄与・相互写像の正則性を主張。実装例（6次元超直方体、ソリトン波、マイクロブラックホール解釈）との四重の分離装置。日英md/tex/pdf全6ファイル |
| 2026-04-22 | 位相方程式篇 論文8（W8）v5：§1–§8全面リストラクチャリング（純粋組合せ論の本体）と§9（解釈例）の厳密分離。Q軸を4値→8値（3ビット符号化）に拡張。c₂c₃遷移9通り＋c₁遷移8通りの明示的数え上げ。符号積Pの完全分類表と証明（命題7.1）。62状態＝SM61＋グラビトン1のマッピング。タイトル変更：「集合構造とその組合せ論的性質」。Zenodo v5公開（DOI: 10.5281/zenodo.19688521） |
| 2026-04-21 | 位相方程式篇 論文9（W9）v5：§11–§12全面リストラクチャリング。§11を変換不変量（W1–W11）、§12を抽象的状態変換モデル（W12–W30）に再構成。物理的メタファー（ばね-質量、流体、縦波）を全面削除。用語統一：「一次形式の保存」「二次形式の保存」。Zenodo公開（DOI: 10.5281/zenodo.19666249） |
| 2026-04-20 | 位相方程式篇 論文8（W8）v4：§8符号付き面積M(σ)を完全削除（SM粒子に対し自明な値しか取らず情報量ゼロのため）。§10.9・表AのM列も削除、セクション番号繰り上げ。Zenodo公開（DOI: 10.5281/zenodo.19665041） |
| 2026-04-20 | 位相方程式篇 論文8（W8）v3全面リライト：全導出を自己完結的に再構成、定義と物理的解釈を厳密分離（本体§1–§9＝組合せ論、§10＝解釈例）、スピンを8種に拡張（n=0,1,2,3）、「配向」→「集合(configuration)」用語統一、表Bに具体的状態数（tQ:18, RQ:18, tRQ:36）、63状態＝SM62＋色シングレット1の明示。Zenodo新バージョン公開（DOI: 10.5281/zenodo.19657042、コンセプトDOI: 10.5281/zenodo.19646651） |
| 2026-04-19 | 位相方程式篇 論文9（W9）v3：§11.6式参照修正（W19）→（W13）。Zenodo v3公開（DOI: 10.5281/zenodo.19651284） |
| 2026-04-19 | 位相方程式篇 論文9（W9）v2：査読指摘対応（§11.2「衝突面」の幾何学的概念を「保存則を介した瞬間的状態変換」に抽象化、§11.6 φ=0ケースの自己整合的説明を追加、§11.8衝突面参照を修正）。Zenodo v2公開（DOI: 10.5281/zenodo.19651253） |
| 2026-04-19 | 位相方程式篇 論文8（W8）v3：査読指摘5点対応（§2.1空間軸パリティ文言修正、§10結論#7の結合強度主張をトーンダウン、§2.1色軸の改段落、結論#6「連続的」→「非有界な」、§9.6電荷値の幾何学的導出を追加）。Zenodo v3公開（DOI: 10.5281/zenodo.19651124） |
| 2026-04-19 | 位相方程式篇 論文9（W9）：sine-Gordon方程式 ── 位相的ソリトンの基礎理論 — Zenodo公開完了（DOI: 10.5281/zenodo.19650966）。sine-Gordon方程式の6つの暗黙の前提を明示（命題9.6a）、矩形波極限での構造保存性を証明（定理9.1）、§11で完全流体の任意位相弾性衝突モデルを独自構成（定理11.1–11.3）。日英md/tex/pdf全6ファイル |
| 2026-04-19 | 位相方程式篇 論文8（W8）v2：R軸を非有界に一般化（W5の離散化はオプション）、電荷・色荷の自由度の明確化（分数電荷1/3=色数、7D拡張の棄却理由）、C行列のQ̄記法導入、§6.5（質量-周波数節）を削除、§7.1結合定数を定性的示唆に修正（α_W≈1/30）、62状態を「SM61+G1」と明記。Zenodo v2公開（DOI: 10.5281/zenodo.19648782） |
| 2026-04-19 | 位相方程式篇 論文8（W8）：6次元超直方体の配向構造から導かれる標準模型粒子の幾何学的分類 — Zenodo公開完了（DOI: 10.5281/zenodo.19646652）。5次元幾何学的超直方体＋離散ラベル軸Q（={0,r,g,b}）による6次元構成。標準模型全62状態（ボソン14＋フェルミオン48）の完全な幾何学的導出。グルーオン8状態=3×3−1、グラビトン=tR型（n=2）の自然な帰結。質量=R含有面の符号付き面積、力の強さ=Qチャネル数として導出。日英md/tex/pdf全6ファイル |
| 2026-04-19 | 位相方程式篇 論文7（W7）：5次元超直方体の配向構造から導かれるスピンの双ベクトル的分類 — Zenodo公開完了（DOI: 10.5281/zenodo.19643358）。W6（正軸体版）の双対多胞体による再定式化。3-cellの法線双ベクトルとしてスピンを直接的に表現。双対定理（定理8.1）により正軸体・超直方体のスピン分類の同型性を証明。4-facet中心=正軸体頂点の一致、32頂点=Cl(5)次元の対応を明示。日英md/tex/pdf全6ファイル |
| 2026-04-18 | 位相方程式篇 論文6（W6）v2：前提論文への依存を排除し自己完結型に改訂。仮定A（軸の3+2物理的分類）を§1.3に新設、中心投影を前提としない構成に変更。§1.2（計量の整数性）を削除。命題4.1をスケール非依存の一般的表現に修正（a₅≫1も許容）。§9.1から中心投影の言及を除去。参考文献を外部文献5件のみに整理。Zenodo v2公開（DOI: 10.5281/zenodo.19643261） |
| 2026-04-17 | 位相方程式篇 論文6（W6）：5次元正軸体の配向構造から導かれるスピンの幾何学的分類 — Zenodo公開完了（DOI: 10.5281/zenodo.19630972）。スピン値6種の完全分類、フェルミオン/ボソン区別の幾何学的導出、力の方向性の符号偶奇則、S₃対称性による3+3+3+1世代構造、物質/反物質非対称性10:9を導出。日英md/tex/pdf全6ファイル |
| 2026-04-17 | 位相方程式篇 論文5（W5）：離散空間における中心投影の全球被覆と5次元背景空間の整数論的必然性 — Zenodo公開完了（DOI: 10.5281/zenodo.19624957）。ラグランジュ四平方定理→4次元必然性、射影幾何→5次元必然性、5軸10面リレーモデル、稠密充填R²整数条件を導出。日英md/tex/pdf全6ファイル |
| 2026-04-16 | 位相方程式篇 論文4（W4）v2：β（時空の最小面積単位, 1.641×10⁻⁶⁹ [L²]）をグループIIに追加、P2をℏ→hに変更、付録A（A.1〜A.11）を全削除。Zenodo v2公開（DOI: 10.5281/zenodo.19609691） |
| 2026-04-16 | 位相方程式篇 論文3（W3）v2：要件8（次元の単一性）・要件9（スケール対称性）を§3.4に追加、§4.3を削除（要件9に昇格）。Zenodo v2公開（DOI: 10.5281/zenodo.19607631） |
| 2026-04-16 | 位相方程式篇 論文4（W4）：標準模型の19個の任意パラメータの分類と構造分析 — Zenodo公開完了（DOI: 10.5281/zenodo.19604965）。日英md/tex/pdf全6ファイル |
| 2026-04-16 | 位相方程式篇 論文3（W3）：万物の理論が満たすべき構造要件の定式化 — Zenodo公開完了（DOI: 10.5281/zenodo.19601592） |
| 2026-04-14 | 位相方程式篇 論文1（W1）・論文2（W2）：Zenodo公開完了。ディレイ回路シリーズ DR1〜DR5：Zenodo公開完了 |
| 2026-04-12 | 論文9・10：Zenodo公開完了（DOI確定）、日英 Markdown・LaTeX・PDF を作成。論文4〜8に読み順ガイド追加。ハンドアウト・RELEASE_NOTES.md を論文10本体制に更新 |
| 2026-04-12 | 論文4〜8：Zenodo公開完了（DOI確定）、日英 Markdown・LaTeX・PDF を作成、RELEASE_NOTES.md・Zenn記事を論文8本体制に更新 |
| 2026-04-09 | ピアレビュー対応：論文2から対称性V（向心加速度）を削除、5つ→4つの対称性に修正。全論文・ハンドアウト・Zenn・note・Zenodoに波及修正を反映 |
| 2026-04-09 | 軽微修正：論文1 命題5.1 添字清書、論文2 定理5.2(ii) 証明簡略化、論文3 命題4.1 記号注記追加 |
| 2026-04-09 | 全8件の tex→PDF を再生成、Zenodo 全3論文のファイル差し替え・メタデータ更新完了 |
| 2026-04-09 | Academia.edu 全3論文：PDF差し替え、論文2 Abstract更新、旧論文3（Charge Quantization）を新論文3（Relativity of Observation）に置換 |
| 2026-04-09 | Zenodo 全3論文：メタデータ修正（タイトル・Description を最新版に統一、日本語タイトル・日本語 Abstract を追加） |
| 2026-04-09 | GitHub リポジトリ：Description・Topics 更新、README.md を論文一覧付きに改訂 |
| 2026-04-09 | Zenn 記事：論文3リライトに対応して全面更新 |
| 2026-04-09 | 論文3：全面リライト（旧版 gnomonic_charge_dimension を削除、新版 subjective_space_relativity として再構成。論文2の対称性II・IVから観測の相対性を導出する構成に変更） |
| 2026-04-09 | ハンドアウト（日英 md/tex/pdf）を論文3リライトに対応して更新 |
| 2026-04-09 | 全プラットフォーム上書き更新（Zenodo・Academia.edu・Zenn・note） |
| 2026-04-08 | 論文2：全面改訂（タイトル・内容変更、ファイル名 gnomonic_lorentzian_origin → gnomonic_projection_symmetry） |
| 2026-04-08 | 論文3：全面改訂（タイトル・内容変更、物理的同定を排除し幾何学的命題に限定） |
| 2026-04-08 | ハンドアウト（日英 md/tex/pdf）新規作成、QRコード追加 |
| 2026-04-08 | 英語版概念図（fig_*_en.png）3点追加 |
| 2026-04-08 | note ヘッダー画像更新、Zenn 記事更新 |
| 2026-04-06 | 論文2・3：Zenodo 公開完了（DOI確定）|
| 2026-04-06 | 論文2・3：日英 LaTeX/PDF を追加 |
| 2026-04-06 | 論文3：テスト粒子応答を「仮定」→「幾何学的帰結」に修正（注意 8.3）|
| 2026-04-06 | 論文3：Maxwell 方程式の完全形を導出済みとして格上げ（命題 8.4）|
| 2026-04-06 | 論文3：$e$ → $e_0$ 統一、概要・結果表の更新 |
| 2026-04-06 | 論文2・論文3の英訳版（`_en.md`）を追加 |
| 2026-04-06 | 論文2：Lorentzian 符号の起源を完成（全10章）|
| 2026-04-06 | 論文3：電荷次元・Coulomb 法則・偏向構造を完成（全8章）|
| 2026-04-05 | Zenodo 既存レコード（19427780）を編集・全6ファイルを最新版に更新・MD5検証完了 |
| 2026-04-05 | GitHub main ブランチに全ファイルをプッシュ |
| 2026-04-05 | Academia.edu にプレプリント（英語版PDF）を登録 |
| 2026-04-05 | Zenn に GitHub連携で記事を公開（articles/ ディレクトリ追加） |
| 2026-04-05 | note に解説記事を公開（バナー画像付き） |
