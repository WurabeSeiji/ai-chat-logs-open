# 第五思考実験：隠れた状態はなかったか
## ― 第四思考実験の二状態軌道生成を自己監査し、五状態から第六永続状態 $\mathcal N=\nu$ まで状態閉包を追跡する ―

**著者:** 木原範昭 (Noriaki Kihara)  
**ORCID:** 0009-0004-6753-4020  
**Version DOI:** 10.5281/zenodo.22973087  
**Concept DOI:** 10.5281/zenodo.22973086  
**版:** v2.1  
**日付:** 2026-09-26  

---

## 要旨

第四思考実験 [4] では、二成分内部状態

$$
X_n=\begin{pmatrix}a_n\\b_n\end{pmatrix}
$$

に状態依存作用 $S_n$ を逐次作用させ、二次読み出し

$$
\Phi(X_n)=
\begin{pmatrix}
a_n^2-b_n^2\\
2a_nb_n
\end{pmatrix}
$$

から、近点移動と leading 2.5PN 放射反作用による縮退を含む重力軌道族を再構成した。しかし、その実装を情報の流れとして自己監査すると、$X_n$ の外側で処理位相 $\chi$、半直弦 $p$、離心率 $e$、方位角 $\phi$、時計 $t$ が保持・更新されており、$X_n$ はそれらから構成される派生表現であった。したがって「二状態だけで軌道を生成した」という読み方は厳密ではない。

第一段階の監査では、マクロstepを越えて保持される五つの永続論理状態

$$
(\chi,p,e,\phi,t)
$$

を明示し、

$$
Z_5=(U,P,E,H,Q)
=(e^{i\chi},p,e,e^{i\phi},e^{t/\tau_0})
$$

へ再符号化した。Paper 4 の RK4 逐次処理は、評価点、$k_1,\ldots,k_4$、中点、合成増分、および11相 one-hot 位相を明示 machine/work state とする同期状態機械へ再配置できた。12条件の長時間照合では、基準更新との差は $H$ で最大 $3.38\times10^{-16}$、$Q$ で最大 $4.00\times10^{-15}$ であった。

しかし追加監査により、放射反作用核

$$
C_{\rm RR}=\frac85\nu r^{-3}
$$

の対称質量比 $\nu$ が状態外の外部引数として残っていることを発見した。したがって第一監査後の系は実際には

$$
Z_5'=F(Z_5;\nu)
$$

であり、完全には状態閉包していなかった。本稿では

$$
\boxed{
Z_6=(U,P,E,H,Q,\mathcal N),\qquad \mathcal N=\nu
}
$$

として第六永続状態へ内部化し、11相すべてで

$$
\mathcal N'=
\sum_{j=0}^{10}q_j(1\cdot\mathcal N)
=\mathcal N
$$

という恒等行として同じ同期写像内に保持した。新遷移関数から外部 `nu` 引数を除去した対照実験では、12条件 $\times$ 4000 macrostep で旧5状態と新6状態の最初の5成分、ならびに $r,x,y,t$ の保存値が bitwise identical となり、$\mathcal N$ の bitwise mismatch も0であった。代表3条件 $\times$ 12000 macrostep でも比較量の最大差は0、$\mathcal N$ の microstep 恒等保持違反は0件であった。保存済み Paper 4 C1 データに対する誤差プロファイルも旧系と完全に一致した。

また、角度級数を $N=12,16,20,24$ と高次化すると、最も厳しい $p=24.3,e=0.55$ の3軌道後の方位角誤差は $1.56\times10^{-7}$、$9.03\times10^{-10}$、$5.54\times10^{-12}$、$3.55\times10^{-14}$ と減少した。時間積分についても同じ物理状態を保ったまま4、6、8次への高次化を検証した。したがって、少なくとも本実験では、状態閉包の修正と数値高精度化は別の操作である。

本稿は、六状態が数学的に最小であること、11相機械が自然界で必然であること、既存B方式全体が後から定めた最も厳しい「完全乗法内部演算」条件を満たすこと、あるいは一般相対論をこの状態系から導出したことを主張しない。示したのは、第四思考実験を自己監査すると、見かけの二状態から五永続状態、さらに外部 $\nu$ を含む六永続状態へと状態閉包が修正され、既知の重力モデルと比較可能な軌道生成を、明示状態・同期更新・非帰還読み出しとして再構成できることである。

**キーワード:** 隠れた状態、状態閉包、同期更新、乗法状態、状態依存作用、対称質量比、ポストニュートン近似、離散力学、数値監査、思考実験

---
# 0. 本稿の位置づけ

本稿は研究プロジェクト設計 [0] のもとでの第五思考実験である。

第一思考実験 [1] では、背景時空や計量を先に置かず、意味を定めない二自由度と面積読み出しの保存から、慣性・一様加速・回転を区別する最小構成を調べた。

第二思考実験 [2] では、その二つの値に固定作用 $S$ を反復し、積の読み出しと面積時計から逆二乗型運動を構成した。

第三思考実験 [3] では向きを逆転し、既知の Newtonian / PN 二体軌道から質量比・相対距離・相対時間を読み出せるかを調べた。

第四思考実験 [4] では、第三思考実験で外部計算していた軌道生成を再び内部相互作用へ戻すため、

$$
X_{n+1}=S_nX_n,
\qquad
S_n=\mathcal S(\mathcal R_n),
\qquad
\det S_n=1
$$

という状態依存作用を構成した。

第五思考実験の目的は、この結果をさらに拡張することではない。問いは

$$
\boxed{
\text{第四思考実験は、本当に明示した状態だけで閉じていたのか。}
}
$$

である。

本稿では札を次の意味で使う。

| 札 | 意味 |
|---|---|
| 【原則】 | 本稿で自己監査のために維持する計算原則 |
| 【設計仮定】 | 現在の表現・実装のために置く仮定。将来外してよい |
| 【構成】 | 既知式または明示状態から作用を構成する手順 |
| 【数値】 | プログラムによる数値照合で確認した結果 |
| 【照合】 | 既知物理・既存数学との比較。基礎原理からの導出ではない |
| 【未証明】 | 現在の実験では決められない事項 |

本稿の最も重要な方針は、**第四思考実験の成功を前提にしない**ことである。隠れた状態が見つかれば、結果に都合よく無視せず状態へ移す。実際、本稿では第一監査で五永続状態を明示した後、さらに外部 $\nu$ を発見し、第六永続状態へ内部化した。この修正過程そのものを結果として記録する。

# 1. 動機：二状態だけで本当に軌道を生成していたのか

第四思考実験の表面上の状態は

$$
X_n=(a_n,b_n)^{\mathsf T}
$$

であり、軌道は

$$
Y_n=\Phi(X_n)
$$

から読み出されていた。

しかし実際のプログラムでは、各 step で

$$
\chi_n,
\quad p_n,
\quad e_n,
\quad \phi_n,
\quad t_n
$$

を別に保持し、そこから

$$
r_n=\frac{p_n}{1+e_n\cos\chi_n},
\qquad
\theta_n=\frac{\phi_n}{2}
$$

を求め、さらに $S_n$ と $X_n$ を構成していた [4]。したがって実際の情報流は

$$
(\chi,p,e,\phi,t)
\longrightarrow
(r,\theta)
\longrightarrow
S_n
\longrightarrow
X
\longrightarrow
Y
$$

であり、$X$ は上流状態から導かれた派生表現であった。

第一監査ではこの上流情報を五永続状態として明示した。しかしさらに放射反作用核まで状態依存性を追跡すると、対称質量比 $\nu$ が状態外から渡されていた。したがって状態監査の経路は

$$
\boxed{
X=(a,b)
\;\longrightarrow\;
Z_5=(U,P,E,H,Q)
\;\longrightarrow\;
Z_6=(U,P,E,H,Q,\mathcal N)
}
$$

となる。

![図1　状態監査の経路：2→5→6](paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure01.png)

**図1.** Paper 4 の見かけの二状態を上流へ追跡すると五つの永続状態が現れ、追加監査で外部 $\nu$ が見つかった。改訂後は $\mathcal N=\nu$ を第六永続状態として同一状態写像に含める。RK4 work state と11相 one-hot 位相 $q$ は別の machine/work state であり、「全実スカラー自由度が6」という主張ではない。

監査対象そのものを見失わないため、Paper 4 で逐次相互作用から再生成した代表条件 C1 の軌道を図2に示す。本稿はこの軌道を新たな重力法則から導出する論文ではなく、この既存の軌道生成がどの永続情報に依存していたかを自己監査する論文である。図2では、保存済み Paper 4 C1 軌道 12001点と、外部 $\nu$ を第六状態 $\mathcal N$ へ内部化した生成器で独立に再生成した 12000 macrostep（3 radial periods）の軌道を重ねた。生成時には保存済み C1 データを状態更新へ帰還させていない。

![図2　Paper 4 C1軌道と第六状態内部化後の再生成軌道](08_ν第6状態内部化_対照実験_20260925/figures/figure_orbit_c1_paper4_vs_6state.png)

**図2.** Paper 4 の保存済み C1 軌道と、$Z_6=(U,P,E,H,Q,\mathcal N)$ による再生成軌道。12000 macrostep は3 radial periodsに対応し、近点移動と放射反作用による軌道縮退を含む。両曲線は描画分解能では重なる。保存済み C1 に対する最大絶対差は $|\Delta x|=3.85\times10^{-12}$、$|\Delta y|=5.11\times10^{-12}$、$|\Delta r|=6.89\times10^{-13}$、$|\Delta t|=1.94\times10^{-10}$ である。これは新しい物理モデルとの独立検証ではなく、Paper 4 の軌道生成を第六状態内部化後も同じ出力として保持できることを可視化したものである。

本稿の問いを四つに分ける。

$$
\boxed{
\text{(1) 次状態を決める情報をすべて状態として数えると何が残るか。}
}
$$

$$
\boxed{
\text{(2) ケース依存の物理パラメータが状態外から流入していないか。}
}
$$

$$
\boxed{
\text{(3) 更新途中値と処理段階を隠さず同期写像として書けるか。}
}
$$

$$
\boxed{
\text{(4) 近似次数の高次化と物理状態の追加を分離できるか。}
}
$$

# 2. 本稿の目的と主張範囲

本稿の目的は一般相対論の新しい厳密解を導出することではない。また、ここで構成する状態遷移則が自然界の基本法則であると主張するものでもない。

目的は、通常の数値計算では暗黙に利用しやすい外部時計、過去状態、処理step番号、更新途中値、ケース依存パラメータ、暗黙のレジスタを監査し、

$$
\boxed{
\text{現在の明示状態}
\rightarrow
\text{相互作用}
\rightarrow
\text{次の明示状態}
\rightarrow
\text{非帰還の読み出し}
}
$$

という形へどこまで閉じられるかを調べることである。

### 本稿で確認したこと

1. Paper 4 の二状態 $X=(a,b)$ だけでは実装は閉じていなかった。
2. 第一監査では、$\chi,p,e,\phi,t$ の五つを永続状態として明示すると Paper 4 の軌道計算を再構成できた。
3. しかし放射反作用核で $\nu$ が外部引数として残っており、五状態系は完全には閉じていなかった。
4. $\mathcal N=\nu$ を第六永続状態として内部化し、$\mathcal N'=\mathcal N$ を同じ one-hot 積和写像の恒等行として実装すると、旧5状態・軌道・時計を bitwise に変更せず外部 `nu` 引数を除去できた。
5. Paper 4 の RK4 逐次計算は、処理段階と work register を明示した11相同期状態機械として再構成できた。
6. 角度級数と数値積分の高次化は、第六状態を追加したこととは独立であり、精度向上のために新たな永続物理状態を導入する必要はなかった。

### 本稿では証明していないこと

1. 六つの永続論理状態が数学的に最小であること。
2. 11相 machine state が自然界で必然であること。
3. 一つの固定作用 $S$ では原理的に不可能であり、作用列が必ず必要であること。
4. 既存B方式全体が、後から定めた最も厳しい「完全乗法内部演算」条件を満たすこと。
5. 近似次数を無限に上げた極限が一般相対論の厳密解へ一致すること。
6. 本稿の相互作用則が自然界の基本法則であること。

# 3. 自己監査規則

## 3.1 現在状態だけから次状態を作る

状態集合を

$$
\Psi_n=(\psi_{1,n},\ldots,\psi_{N,n})
$$

とする。

一回の相互作用では、現在状態を一度に読み、次状態を一度に書く。

$$
\boxed{
\psi_{i,n+1}=F_i(\Psi_n)
}
$$

でなければならない。

同じ一回の相互作用中に先に計算した $\psi_{j,n+1}$ を別成分の計算に使うことは禁止する。

---

## 3.2 増分代入と積和を区別する

プログラム上の

```text
psi += delta
```

は、旧値と新値を暗黙に保持する逐次更新として読めるため、基本相互作用としては採用しない。

一方、現在状態だけから

$$
\psi_i'
=
\sum_{\alpha}
c_{i\alpha}
\prod_j \psi_j^{m_{i\alpha j}}
$$

を一度に計算して書き込む**積和**は許容する。

例えば

$$
P'=P+D_p
$$

は、逐次代入ではなく

$$
P'=1\cdot P+1\cdot D_p
$$

という二つの現在状態の積和として一度に計算する限り、本稿の規則に反しない。

---

## 3.3 演算回数を外部から参照しない

外部の loop counter によって作用を選んではならない。

処理段階の区別が必要なら、one-hot 位相状態

$$
q_n=(q_{0,n},\ldots,q_{m-1,n})
$$

を系内部に置き、巡回置換

$$
\boxed{
q_{n+1}=C_m q_n,
\qquad C_m^m=I
}
$$

で更新する。

作用の選択も数式上は

$$
\boxed{
\Psi_{n+1}
=
\sum_{j=0}^{m-1} q_{j,n}F_j(\Psi_n)
}
$$

と積和で書く。

実際の監査コードでは $q$ が厳密な one-hot であるため高速化として `argmax(q)` と分岐を使用するが、これは物理的な外部情報ではなく、上式の one-hot 選択を実装したものとして扱う。したがって、**数式上の基本定義は $q$ による積和選択**であり、`argmax` そのものを基本相互作用とはみなさない。

---

## 3.4 読み出しは状態へ書き戻さない

観測のための写像

$$
O=\mathcal O(\Psi)
$$

は許される。

しかし、特に

$$
\log
$$

を用いて加法座標へ戻し、旧来の加法計算を行って再び指数化する

$$
\Psi
\rightarrow
\log\Psi
\rightarrow
\text{additive update}
\rightarrow
\exp
\rightarrow
\Psi'
$$

という手順は基礎モデルとして認めない。

読み出しは

$$
\boxed{
\text{state}\rightarrow\text{observation}
}
$$

の一方向とする。

---

## 3.5 定数

固定定数との積和は、由来が明確なら許容する。

例えば RK4 の

$$
\frac16(1,2,2,1)
$$

や、級数

$$
c_{m+1}=\frac{2m+1}{m+1}c_m
$$

の係数は導出根拠が明確である。

一方、既知軌道に近づけるためだけの自由な fitting coefficient は用いない。

![図3　改訂版の状態閉包と同期更新](paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure02.png)

**図3.** 改訂後の監査範囲。旧 full state を一度に読み、同一 microstep 写像で次 full state を生成する。読み出し $x,y,r,t$ は状態更新へ戻さない。$\mathcal N$ の恒等保持も同じ写像内に含める。一方、六状態の最小性、11相機械の自然界での必然性、既存B方式全体の完全乗法性は主張しない。

---

# 4. 既存理論での位置づけ

本稿の個々の計算概念そのものを新規数学として主張するものではない。

少数の変数へ射影した有効記述では、消去した自由度の影響が memory term として現れ得ることは Zwanzig の射影形式に代表される [5]。本稿ではこの結果を六状態の導出に使うのではなく、**現在状態だけで閉じないときに、欠落した情報の所在を監査する**という背景としてのみ参照する。

また、発展を複数作用の積へ分解する表現には Trotter の積公式 [6] が古典的な先例を与える。本稿の11相状態機械は Trotter 分解そのものではないが、複数の作用を明示列として扱うこと自体は既存数学の範囲にある。

数値積分については、高次 Runge--Kutta 過程の系統的構成は Butcher により古くから研究されている [7]。本稿で4、6、8次への高次化を調べる目的は新しい積分法を提案することではなく、**数値近似の深度と永続物理状態数を区別すること**にある。

重力側では、Schwarzschild 束縛軌道の角度進行に Darwin のパラメータ表示 [8] を用い、時間発展する軌道要素という構成の近縁な標準的背景として osculating-orbit 形式 [9] を参照する。ポストニュートン近似と放射反作用の一般的背景は Blanchet のレビュー [10] に、leading 2.5PN radiation-reaction 項を含む放射反作用式の系譜は Blanchet, Faye, Trestini [11] に置く。

したがって本稿固有の点は、これら既存要素の新規発明ではない。第四思考実験の実装に対して、状態漏れ・更新途中値・処理位相・読み出し帰還を同時に監査し、途中で見つかった外部 $\nu$ も含めて状態閉包を修正したことにある。

# 5. Paper 4 の状態監査

## 5.1 表面状態 $X=(a,b)$ は基本状態ではなかった

Paper 4 では

$$
X_{n+1}=S_nX_n
$$

を主要な状態遷移として表示した [4]。

しかし実装では先に

$$
(\chi_n,p_n,e_n,\phi_n,t_n)
$$

を更新し、

$$
r_n=\frac{p_n}{1+e_n\cos\chi_n},
\qquad
\theta_n=\frac{\phi_n}{2}
$$

から

$$
X_n
=
\sqrt{r_n}
\begin{pmatrix}
\cos\theta_n\\
\sin\theta_n
\end{pmatrix}
$$

を構成していた。したがって $X_n$ は上流状態の圧縮された派生表現であり、それだけでは次状態を決定できない。

---

## 5.2 第一監査：五つの永続状態

第一監査では、現在の表現でマクロstepをまたいで必要な状態として

$$
(\chi,p,e,\phi,t)
$$

を明示した。

- $\chi$：処理位相、
- $p$：半直弦、
- $e$：離心率、
- $\phi$：軌道方位角、
- $t$：時計読み出しの累積座標。

これを

$$
Z_5=(U,P,E,H,Q)
=(e^{i\chi},p,e,e^{i\phi},e^{t/\tau_0})
$$

へ再符号化した。五状態同期写像を Paper 4 の元計算と12条件で照合したところ、軌道読み出しまで含めて既存計算を再現できた。この段階で一度「五状態で閉じた」と判断した。

しかし、この結論は追加監査で修正される。

### 表1　状態監査の最終整理

| 量 | 役割 | 分類 | 監査理由 | 監査段階 |
|---|---|---|---|---|
| $X=(a,b)$ | Paper 4で表示された二成分 | 派生表現 | $r,\phi$ から構成され、単独では次状態を閉じない | Paper 4→第1監査 |
| $U=e^{i\chi}$ | 処理位相 | 永続論理状態 | $\chi$ の周期情報を保持 | 5状態監査 |
| $P=p$ | 軌道要素 | 永続論理状態 | 放射反作用で更新 | 5状態監査 |
| $E=e$ | 軌道要素 | 永続論理状態 | 放射反作用で更新 | 5状態監査 |
| $H=e^{i\phi}$ | 軌道方位位相 | 永続論理状態 | 角度進行を乗法的に保持 | 5状態監査 |
| $Q=e^{t/\tau_0}$ | 時計状態 | 永続論理状態 | $t$ は読み出し時のみ $\log$ で取得 | 5状態監査 |
| $\mathcal N=\nu$ | 対称質量比 | 永続論理状態 | 旧実装では外部引数。追加監査で状態漏れと判定 | 追加 $\nu$ 監査 |
| $q_0,\ldots,q_{10}$ | 11相 one-hot 処理位相 | machine state | 外部 loop counter の代替 | 11相同期監査 |
| RK work registers | $k_1,\ldots,k_4$ 等 | work/machine state | 更新途中値を暗黙保持しないため明示 | 11相同期監査 |

---

## 5.3 追加監査：外部 $\nu$ は第六状態だった

放射反作用核では

$$
C_{\rm RR}=\frac85\nu r^{-3}
$$

を使うため、第一監査後の実装は実際には

$$
\boxed{Z_5'=F(Z_5;\nu)}
$$

であった。$\nu$ は固定質量比のケースごとに与えられる物理情報であり、状態写像の外から供給されていた。したがって、状態閉包を厳密に読む限り五状態系は閉じていない。

そこで

$$
\boxed{
Z_6=(U,P,E,H,Q,\mathcal N),
\qquad
\mathcal N=\nu
}
$$

とし、放射反作用核を

$$
C_{\rm RR}=\frac85\mathcal N r^{-3}
$$

へ置き換えた。$\mathcal N$ は11相すべてで係数1を持つ恒等行として

$$
\boxed{
\mathcal N_{m+1}
=
\sum_{j=0}^{10}q_{j,m}(1\cdot\mathcal N_m)
=
\mathcal N_m
}
$$

と保持する。これは状態機械の外で `N_next=N_old` とコピーするのではなく、one-hot 積和写像の第六行として実装した。

旧遷移

```text
micro_fast(s, nu)
macro_fast(s, nu)
```

は

```text
micro_fast6(s)
macro_fast6(s)
```

となり、新遷移関数は外部 `nu` 引数を持たない。

### $\nu$ 内部化の対照結果

12条件 $\times$ 4000 macrostep の48,012保存行で、旧5状態と新6状態の最初の5成分は

$$
\boxed{\Delta U=\Delta P=\Delta E=\Delta H=\Delta Q=0}
$$

であり、読み出しも

$$
\boxed{\Delta r=\Delta x=\Delta y=\Delta t=0}
$$

であった。complex128 / float64 の保存ビット列比較でも、旧5状態と新系最初の5状態の bitwise mismatch は0、$\mathcal N$ と初期 $\nu_0$ の bitwise mismatch も0であった。

さらに代表3条件 $\times$ 12000 macrostep でも

$$
\boxed{
\max|\Delta(U,P,E,H,Q,r,x,y,t)|=0
}
$$

かつ $\mathcal N$ の microstep 恒等保持違反は0件であった。保存済み Paper 4 C1 CSV に対する $p,e,t,x,y$ の誤差プロファイルも旧5状態系と新6状態系で完全に同一だった。【数値】

この対照実験が示すのは、**物理式を変更せず、情報の置き場所だけを外部引数から第六状態へ移せた**ことである。荷電相互作用や別の物理自由度を追加した結果ではない。

---

## 5.4 「六状態」の意味

改訂後の正確な主張は

$$
\boxed{
\text{現在の表現では、マクロstepを越えて保持される永続論理状態が六つ}
}
$$

である。

これは全状態機械の実スカラー成分が六個という意味ではない。RK4 を隠れた逐次レジスタなしに実装するための $k_1,\ldots,k_4$、評価点、中点、合成増分、11相 one-hot 位相は別の明示 machine/work state として存在する。

また、六が数学的最小であることも未証明である。したがって

$$
\boxed{
\text{現在の Paper 4 表現を自己監査すると、少なくとも六つの永続論理状態を明示する必要があった}
}
$$

というのが本稿の結論である。

# 6. 乗法状態への再符号化

第一監査で得た加法座標

$$
(\chi,p,e,\phi,t)
$$

に、追加監査で見つかった対称質量比 $\nu$ を加え、改訂後の永続論理状態を

$$
\boxed{
Z_6=(U,P,E,H,Q,\mathcal N)
}
$$

とする。定義は

$$
U=e^{i\chi},
\qquad P=p,
\qquad E=e,
\qquad H=e^{i\phi},
\qquad Q=e^{t/\tau_0},
\qquad \mathcal N=\nu.
$$

位相と時計は

$$
U'=Ue^{i\Delta\chi},
$$

$$
H'=He^{i\Delta\phi},
$$

$$
Q'=Qe^{\Delta t/\tau_0}
$$

という乗法更新になる。$\mathcal N$ は固定質量比モデルでは

$$
\mathcal N'=\mathcal N
$$

であり、同一状態写像の恒等行として保持される。

一方、$P,E$ は現在の物理式で正の軌道要素として直接現れるので、本稿では無理に $\log P,\log E$ を永続状態として持たない。

【設計仮定】現在の検証プログラムは $U,H$ を複素数で簡潔に実装している。しかしこれは基本状態が本質的に複素数でなければならないという結論ではない。例えば $U=(\cos\chi,\sin\chi)$ を実二成分とし、回転行列で更新することもできる。したがって「六状態」は**論理状態量の数**についての記述であり、最小実スカラー成分数を確定するものではない。

# 7. B方式：11相の同期内部状態機械

## 7.1 なぜ内部相が必要か

Paper 4 の放射反作用更新は RK4 であり、通常のコードでは

$$
k_1\rightarrow y_2\rightarrow k_2\rightarrow y_3\rightarrow k_3\rightarrow y_4\rightarrow k_4
$$

と逐次処理する。

これを一つの物理相互作用だとみなすと、処理途中の値が隠れたレジスタになる。

そこで各中間値を明示状態とし、11 回の microstep に分解した。

内部状態を模式的に

$$
\mathcal Z
=(Z,W,q)
$$

と書く。$Z_6=(U,P,E,H,Q,\mathcal N)$ が永続状態、$W$ が $k_1,\ldots,k_4$、評価点、合成増分、中点などの work state、$q$ が 11 相 one-hot 状態である。

一 microstep は

$$
\boxed{
\mathcal Z_{m+1}
=
\sum_{j=0}^{10}q_{j,m}F_j(\mathcal Z_m)
}
$$

であり、同時に

$$
\boxed{
q_{m+1}=C_{11}q_m
}
$$

とする。

11 microstep 後に

$$
q_{m+11}=q_m
$$

となり、一つのマクロ状態更新が完了する。

---

## 7.2 実際の監査コード

以下は $\nu$ 内部化後の対照コードの中核である。`o` は更新前 full state、`n` は次状態の書込先であり、各分岐は OLD state だけを読む。放射反作用は外部 `nu` ではなく `o['N']` を読む。

```python
def micro_fast6(s):
    j = int(np.argmax(s['q']))
    o = clone_state(s)          # OLD full state
    n = clone_state(s)          # next-state buffer

    U, P, E, N = o['U'], o['P'], o['E'], o['N']

    if j == 0:
        n['k1'] = rr(U, P, E, N)
    elif j == 1:
        n['Us'] = U * PH2
        n['Ps'] = P + 0.5 * H * o['k1'][0]
        n['Es'] = E + 0.5 * H * o['k1'][1]
    elif j == 2:
        n['k2'] = rr(o['Us'], o['Ps'], o['Es'], N)
    elif j == 3:
        n['Us'] = U * PH2
        n['Ps'] = P + 0.5 * H * o['k2'][0]
        n['Es'] = E + 0.5 * H * o['k2'][1]
    elif j == 4:
        n['k3'] = rr(o['Us'], o['Ps'], o['Es'], N)
    elif j == 5:
        n['Us'] = U * PH
        n['Ps'] = P + H * o['k3'][0]
        n['Es'] = E + H * o['k3'][1]
    elif j == 6:
        n['k4'] = rr(o['Us'], o['Ps'], o['Es'], N)
    elif j == 7:
        n['d'] = (H / 6) * (
            o['k1'] + 2*o['k2'] + 2*o['k3'] + o['k4']
        )
    elif j == 8:
        n['Um'] = U * PH2
        n['Pm'] = P + 0.5 * o['d'][0]
        n['Em'] = E + 0.5 * o['d'][1]
    elif j == 9:
        n['dphi'] = H * gN(o['Um'], o['Pm'], o['Em'])
    elif j == 10:
        n['U'] = U * PH
        n['P'] = P + o['d'][0]
        n['E'] = E + o['d'][1]
        n['H'] = o['H'] * complex(
            math.cos(o['dphi']), math.sin(o['dphi'])
        )
        n['Q'] = o['Q'] * math.exp(o['d'][2] / TAU0)

    # sixth persistent state: same one-hot interaction, identity row
    Nc = [1.0 * N for _ in range(11)]
    n['N'] = float(blend(o['q'], Nc))

    n['q'] = np.roll(o['q'], 1)
    return n
```

コード中の

```python
j = int(np.argmax(s['q']))
```

は one-hot 状態の非零成分を高速に選ぶ実装上の短縮である。数式上の相互作用定義は外部 loop counter ではなく

$$
\sum_{j=0}^{10}q_jF_j
$$

である。

一マクロstepは

```python
def macro_fast6(s):
    for _ in range(11):
        s = micro_fast6(s)
    return s
```

で検証した。ここでも `for` の回数を状態遷移の入力には使わない。作用選択に使われるのは $q$ だけである。

---

## 7.3 各内部相の数式

$h=\Delta\chi$、

$$
\Phi_{1/2}=e^{ih/2},
\qquad
\Phi_1=e^{ih}
$$

とする。

第1評価は

$$
K_1=R(U,P,E,\mathcal N),
$$

第2評価点は

$$
U_2=U\Phi_{1/2},
$$

$$
P_2=P+\frac h2 K_{1p},
\qquad
E_2=E+\frac h2 K_{1e},
$$

第2評価は

$$
K_2=R(U_2,P_2,E_2,\mathcal N).
$$

同様に

$$
U_3=U\Phi_{1/2},
\quad
P_3=P+\frac h2K_{2p},
\quad
E_3=E+\frac h2K_{2e},
$$

$$
K_3=R(U_3,P_3,E_3,\mathcal N),
$$

$$
U_4=U\Phi_1,
\quad
P_4=P+hK_{3p},
\quad
E_4=E+hK_{3e},
$$

$$
K_4=R(U_4,P_4,E_4,\mathcal N).
$$

合成は

$$
\boxed{
D
=\frac h6
(K_1+2K_2+2K_3+K_4)
}
$$

である。

角度評価点を

$$
U_m=U\Phi_{1/2},
\qquad
P_m=P+\frac12D_p,
\qquad
E_m=E+\frac12D_e
$$

とし、

$$
\Delta\phi=h\,g_N(U_m,P_m,E_m)
$$

から最終的に

$$
\boxed{
U'=U\Phi_1
}
$$

$$
\boxed{
P'=P+D_p,
\qquad
E'=E+D_e
}
$$

$$
\boxed{
H'=He^{i\Delta\phi}
}
$$

$$
\boxed{
Q'=Qe^{D_t/\tau_0}
}
$$

$$
\boxed{
\mathcal N'=\mathcal N
}
$$

とする。

これは `P += Dp` のような逐次代入ではない。$P'$ と $E'$ は microstep 開始時の旧状態だけから一度に書かれる積和である。

---

# 8. 相互作用核の完全展開

## 8.1 $U$ から三角量を読む

$$
U=e^{i\chi}
$$

より

$$
C_\chi:=\cos\chi
=\frac{U+U^{-1}}2,
$$

$$
S_\chi:=\sin\chi
=\frac{U-U^{-1}}{2i}.
$$

したがって $\chi=\arg U$ を状態更新内で復元する必要はない。

---

## 8.2 2.5PN 放射反作用核

Paper 4 と同じ leading 2.5PN 放射反作用を用いる。PN理論の背景は [10]、放射反作用式の系譜は [11] を参照する。

$$
d=1+EC_\chi,
\qquad
r=\frac Pd,
\qquad
\sqrt P=P^{1/2},
$$

$$
\dot r=\frac{ES_\chi}{\sqrt P},
\qquad
v_t=\frac d{\sqrt P},
$$

$$
v^2=\dot r^2+v_t^2,
\qquad
r^{-1}=\frac dP.
$$

放射反作用係数を

$$
C_{\rm RR}=\frac85\mathcal N r^{-3},
$$

$$
A
=\dot r
\left(
18v^2+\frac23r^{-1}-25\dot r^2
\right),
$$

$$
B
=-\left(
6v^2-2r^{-1}-15\dot r^2
\right)
$$

とすると、

$$
\dot E_{\rm orb}
=C_{\rm RR}(A\dot r+Bv^2),
$$

$$
\dot p=2C_{\rm RR}Bp,
$$

$$
E_{\rm orb}=\frac{e^2-1}{2p},
$$

$$
\dot e
=\frac{p\dot E_{\rm orb}+E_{\rm orb}\dot p}{e}.
$$

また

$$
\frac{dt}{d\chi}
=\frac{r^2}{\sqrt p}.
$$

実際の次数検証コードの核は次である。

```python
def rr_numeric(U, P, E, N):
    Cc = 0.5 * (U + 1.0 / U)
    Ss = (U - 1.0 / U) / (2j)

    den = 1.0 + E * Cc
    r = P / den
    sqP = P**0.5

    rdot = E * Ss / sqP
    vt = den / sqP
    v2 = rdot * rdot + vt * vt
    invr = den / P

    Crr = 1.6 * N * invr**3

    an = rdot * (
        18.0 * v2
        + (2.0 / 3.0) * invr
        - 25.0 * rdot * rdot
    )

    bv = -(
        6.0 * v2
        - 2.0 * invr
        - 15.0 * rdot * rdot
    )

    Edot = Crr * (an * rdot + bv * v2)
    pdot = 2.0 * Crr * bv * P

    Eorb = (E * E - 1.0) / (2.0 * P)
    edot = (P * Edot + Eorb * pdot) / E

    dt = r * r / sqP

    return pdot * dt, edot * dt, dt
```

この関数は旧状態の $U,P,E,\mathcal N$ を入力に取り、$dp/d\chi,de/d\chi,dt/d\chi$ を返す。

---

## 8.3 角度進行核

Schwarzschild 束縛軌道の Darwin パラメータ表示 [8] から

$$
\frac{d\phi}{d\chi}
=
\left[
1-\frac{2(3+e\cos\chi)}p
\right]^{-1/2}
$$

を用いる。

$$
u
:=
\frac{3+EC_\chi}{P}
$$

と置けば

$$
(1-2u)^{-1/2}
=
\sum_{m=0}^{\infty}c_mu^m,
$$

$$
c_0=1,
\qquad
c_{m+1}=\frac{2m+1}{m+1}c_m.
$$

$N$ 次で

$$
\boxed{
g_N(U,P,E)
=\sum_{m=0}^{N}c_m
\left[
\frac{3+\frac E2(U+U^{-1})}{P}
\right]^m
}
$$

である。

コードでは Horner 法で

```python
def gN_numeric(U, P, E, N):
    C = 0.5 * (U + 1.0 / U)
    u = (3.0 + E * C) / P

    c = COEFF[N]
    y = complex(c[-1])
    for a in c[-2::-1]:
        y = complex(a) + u * y
    return y
```

と評価する。

---

## 8.4 補助関数を消した完全展開の照合

途中の検討では `rr()` や `gN()` を残した段階を「完全展開」と誤って呼んだ。これは完全展開ではなかった。

最終監査では $dp/d\chi$ と $de/d\chi$ を $U^{\pm1}$ の有限 Laurent 多項式へ展開し、$g_{12}$ も係数を明示した上で、元式と 100,000 ランダム点で照合した。

最大誤差は

$$
|\Delta(dp/d\chi)|_{\max}
\approx1.25\times10^{-16},
$$

$$
|\Delta(de/d\chi)|_{\max}
\approx2.08\times10^{-17},
$$

$$
\left|\frac{\Delta(dt/d\chi)}{dt/d\chi}\right|_{\max}
\approx8.61\times10^{-16},
$$

$$
|\Delta g_{12}|_{\max}
\approx8.88\times10^{-16}
$$

であり、浮動小数点丸めの範囲で一致した。【数値】

---

# 9. 読み出し

## 9.1 軌道位置

状態 $U,P,E,H$ から

$$
C_\chi=\frac{U+U^{-1}}2
$$

を計算し、

$$
\boxed{
r=\frac{P}{1+EC_\chi}
}
$$

とする。

$$
H=e^{i\phi}
$$

なので

$$
\boxed{
x=r\operatorname{Re}H,
\qquad
y=r\operatorname{Im}H
}
$$

である。

---

## 9.2 時計

$$
Q=e^{t/\tau_0}
$$

なので観測時だけ

$$
\boxed{
t_{\rm obs}
=\tau_0\log |Q|
}
$$

と読む。

ここで得た $t_{\rm obs}$ を状態更新へ戻さない。

監査コード中の読み出し演算を関数として整理すると次である。これは監査で実際に用いた演算と同値である。

```python
def readout(U, P, E, H, Q, tau0):
    cos_chi = 0.5 * (U + 1.0 / U)
    r = (P / (1.0 + E * cos_chi)).real

    x = r * H.real
    y = r * H.imag

    # log is used only here, never in the state update
    t = tau0 * math.log(abs(Q))

    return x, y, r, t
```

したがって本稿の基本的な情報流は

$$
\boxed{
(U,P,E,H,Q)
\xrightarrow{\text{interaction}}
(U',P',E',H',Q')
}
$$

と

$$
\boxed{
(U,P,E,H,Q)
\xrightarrow{\text{readout}}
(x,y,r,t)
}
$$

に分離される。

---

# 10. B方式と第六状態のプログラム監査結果

B方式について、単に最終軌道が似ているかではなく、状態機械としての規則を直接監査した。

第一監査の11相機械では、

1. 同一 microstep で更新後状態を読む箇所がないか、
2. one-hot 位相が壊れないか、
3. commit 相以外で永続状態が不正に変化しないか、
4. 11 microstep 後に位相が元へ戻るか、
5. 一マクロstepが元の Paper 4 RK4 step と一致するか、

を検査した。

500個のランダム状態による一マクロstep比較では、第一監査の五成分について

$$
\boxed{
\max|\Delta(U,P,E,H,Q)|=(0,0,0,0,0)
}
$$

であった。12条件の長時間比較では

$$
\boxed{
\max|\Delta(U,P,E,H,Q)|
=
(0,0,0,3.38\times10^{-16},4.00\times10^{-15})
}
$$

であった。

追加 $\nu$ 監査では、12条件 $\times$ 4000 macrostep で旧5状態と新6状態の最初の5成分、および $r,x,y,t$ が保存値の bit 列まで一致し、$\mathcal N$ の bitwise mismatch も0であった。代表3条件 $\times$ 12000 macrostep でも全比較量の最大差は0で、$\mathcal N$ の microstep 恒等保持違反は0件だった。

保存済み C1 CSV に対する最大絶対差は旧系・新系で同一であり、

$$
|\Delta p|_{\max}=4.97\times10^{-14},
\qquad
|\Delta e|_{\max}=7.77\times10^{-16},
$$

$$
|\Delta t|_{\max}=1.94\times10^{-10},
$$

$$
|\Delta x|_{\max}=3.85\times10^{-12},
\qquad
|\Delta y|_{\max}=5.11\times10^{-12}
$$

であった。この微小差は今回の $\nu$ 内部化で生じたものではなく、旧B監査系が保存済みC1に対してもともと持っていた丸め・実装経路上の差である。【数値】

![図4　12条件におけるB方式の厳密監査](paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure03.png)

**図4.** 12条件での B方式と基準更新の差。非零成分もほぼ浮動小数点丸めの範囲にある。$\nu$ 内部化の対照では、旧5状態と新6状態の最初の5成分は bitwise identical であり、第六状態 $\mathcal N$ も恒等保持された。

以上から、現在の限定されたモデルについては

$$
\boxed{
\text{Paper 4 の逐次処理を、外部 }\nu\text{ を含まない明示状態機械へ再配置できる}
}
$$

ことを確認した。ただし、これは11相状態機械が自然界の基本作用であること、あるいは既存B方式全体が後から定めた最厳格の完全乗法内部演算条件を満たすことを意味しない。

# 11. 近似次数を分離する

Paper 4 と本稿には複数の「次数」が存在する。これらを混同してはならない。

### 表2　近似次数と状態閉包の区別

| 層 | 現在/検証次数 | 誤差・未保持項 | 意味 | 状態数との関係 |
|---|---|---|---|---|
| 散逸物理モデル | leading 2.5PN radiation reaction | 高次PN項は未導入 | 物理モデル近似次数 | 六永続状態とは別問題 |
| 角度級数 | $N=12$ baseline; tested $12,16,20,24$ | $O(u^{N+1})$ | 級数打切り | 状態数を増やさず高次化 |
| Paper 4時間積分 | RK4 | local $O(h^5)$, global $O(h^4)$ | 数値積分次数 | 11相展開は処理順の明示化 |
| B 11相展開 | RK4の代数的再配置 | 新たな打切り誤差なし | machine-state深度 | 物理的11段階性の主張ではない |
| Taylor次数検証 | 4,6,8 | 実測収束を別途確認 | 数値近似深度 | 永続状態数を増やさず高次化 |
| $\nu$ 内部化 | $\mathcal N=\nu$ identity state | 近似ではない | 状態閉包の修正 | 5→6は精度向上ではなく情報漏れ除去 |

特に

$$
\boxed{
\text{数値積分次数を上げること}
\neq
\text{PN物理次数を上げること}
}
$$

である。

また

$$
\boxed{
\text{近似次数を上げること}
\neq
\text{物理状態を追加すること}
}
$$

であるかを次に検証する。

---

# 12. 角度級数の高次化

角度因子の厳密式は

$$
g_\infty(u)=(1-2u)^{-1/2}
$$

であり、有限 $N$ では

$$
g_N(u)=\sum_{m=0}^{N}c_mu^m.
$$

永続状態 $Z_6=(U,P,E,H,Q,\mathcal N)$ を固定し、角度級数の打切り次数だけを

$$
N=12,16,20,24
$$

へ増加させた。

最も厳しい試験条件の一つである

$$
p_0=24.3,
\qquad e_0=0.55
$$

では、3軌道後の方位角誤差が

| $N$ | $|\Delta\phi|$ |
|---:|---:|
| 12 | $1.56\times10^{-7}$ |
| 16 | $9.03\times10^{-10}$ |
| 20 | $5.54\times10^{-12}$ |
| 24 | $3.55\times10^{-14}$ |

と減少した。

![図5　角度級数の次数と誤差](paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure04.png)

**図5.** 9条件における $N=12,16,20,24$ の収束。状態数は変更せず、同じ級数生成則の打切り位置だけを変えている。

この結果から、少なくとも有限範囲では

$$
\boxed{
N\uparrow
\quad\text{while}\quad
\dim Z_6=6
}
$$

が成立する。

---

# 13. 六永続状態のまま数値積分次数を上げられるか

## 13.1 直接 Taylor 漸化

RK4 の段数を増やすと、内部 stage state も増える。

そこで「高次化のために新しい永続物理状態が必要か」を分離して調べるため、同じ六永続状態のうち $\mathcal N$ を恒等保持したまま、進化する五成分

$$
Z=(U,P,E,H,Q)
$$

に対して Taylor 係数を一 step 内だけで生成した。

連続ベクトル場を

$$
\frac{dZ}{d\chi}=F_N(Z)
$$

と書くと、

$$
Z(\chi+h)
=\sum_{k=0}^{p}A_k h^k+O(h^{p+1}),
$$

$$
A_0=Z_n,
$$

$$
\boxed{
(k+1)A_{k+1}
=[\xi^k]
F_N\!\left(\sum_{j=0}^{k}A_j\xi^j\right)
}
$$

とする。

コードの中核は次である。

```python
def taylor_step(z, h, order, N, tau0):
    # Same five dynamically evolving components.
    # Revised closure appends N as an identity persistent state.
    # No RK stage is promoted to a new persistent state.
    a = np.zeros((5, order + 1), complex)
    a[:, 0] = z

    for k in range(order):
        series = [TPS(a[j], order) for j in range(5)]
        U, P, E, H, Q = series

        dp, de, dt = rr_tps(U, P, E)
        g = gN_tps(U, P, E, N)

        F = [
            1j * U,
            dp,
            de,
            1j * g * H,
            (dt / tau0) * Q,
        ]

        for j in range(5):
            a[j, k + 1] = F[j].c[k] / (k + 1)

    powers = np.array(
        [h**k for k in range(order + 1)],
        complex
    )
    return a @ powers
```

ここで Taylor 係数 $A_k$ は一 step の計算係数であり、次のマクロstepへ保持しない。

この数値収束実験では進化する五成分だけを配列化しているが、改訂後の閉包では $\mathcal N$ を恒等状態として併置する。したがって物理式を変えず、六永続状態のまま $p=4,6,8$ を比較できる。

---

## 13.2 実測収束次数

9条件で step/orbit を 8、16、32、64 と変えた。

丸め誤差の影響が比較的小さい 16→32 step/orbit の収束率を 9条件平均で見ると、

### 表3　16→32 step/orbit の平均実測次数

| 設計次数 | $P$ | $E$ | $\phi$ |
|---:|---:|---:|---:|
| 4 | 3.998 | 3.993 | 3.928 |
| 6 | 5.997 | 6.010 | 5.921 |
| 8 | 7.918 | 8.072 | 7.916 |

となった。

4次・6次ではさらに 32→64 でもほぼ設計次数を保つ。8次では $P,E$ が機械精度へ近づき始めるため、より細かい step では参照解誤差と丸め誤差の影響が無視できなくなる。

![図6　同じ永続状態での4・6・8次収束](paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure05.png)

**図6.** 同じベクトル場・同じ進化する五成分に対して近似深度だけを変更した収束試験。改訂後の第六状態 $\mathcal N$ は恒等保持され、この高次化とは独立である。

9条件の全体像を図7に示す。

![図7　9条件における高次化の比較](paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure06.png)

**図7.** 強場・中間・弱場および離心率の異なる条件でも、高次化に伴う誤差低下が確認できる。

したがって本稿で確認されたのは

$$
\boxed{
\text{近似次数の増加}
\not\Rightarrow
\text{新しい永続物理状態の追加}
}
$$

である。

ただし計算量や一時 Taylor 係数の個数は次数とともに増える。

---

# 14. 固定作用 $S$ と状態依存作用列について

本稿の監査では、現在の Paper 4 の表現をそのまま明示状態へ展開した場合、一つの固定作用だけでは処理を表せず、処理位相に応じた複数の内部作用を使用した。

数式上は

$$
\boxed{
\mathcal Z_{m+1}
=F_{\sigma_m}(\mathcal Z_m),
\qquad
\sigma_m\ \text{is encoded in }q_m
}
$$

である。

しかし、この結果から

$$
\boxed{
\text{自然界では必ず }S_n\text{ が必要}
}
$$

とは結論できない。

理由は三つある。

第一に、別の状態座標を選べば、一つの固定作用へ統合できる可能性がある。

第二に、内部11相は RK4 という数値アルゴリズムの処理構造を明示化したものであり、物理そのものの11段階性を証明するものではない。

第三に、Taylor 漸化では同じ永続状態集合と一つのベクトル場 $F_N$ のまま高次化できた。

したがって本稿の結論は

$$
\boxed{
\text{現在の B 構成では作用列が必要だったが、その原理的必然性は未証明}
}
$$

である。

---

# 15. 本稿から分かったこと

## 15.1 「二状態モデル」への自己修正

Paper 4 の $X=(a,b)$ は有用な読み出し表現だったが、次状態を決める全情報を含む基本状態ではなかった。これは Paper 4 の数値結果を無効にするものではない。Paper 5 により、その結果がどの上流情報に依存していたかが明示された。

## 15.2 「五状態で閉じた」という第一監査も再修正された

第一監査では $Z_5=(U,P,E,H,Q)$ まで明示して閉じたと判断した。しかし放射反作用核を追跡すると $\nu$ が外部から供給されていた。

したがって本稿自身が、監査途中の結論を

$$
\boxed{
Z_5'=F(Z_5;\nu)
\quad\longrightarrow\quad
Z_6'=\widetilde F(Z_6)
}
$$

へ修正した。ここで重要なのは「6」という数字より、**次状態を決める情報を計算補助値として隠さず、見つかれば状態候補へ移す**という監査手順である。

## 15.3 隠れた更新途中値を明示すると状態機械になる

通常の逐次数値コードでは、評価点や積分途中の値は「計算の都合」として扱われる。しかし本稿の規則では、次の演算に使うなら情報を保持している以上、明示 work/machine state として扱う。その結果、Paper 4 の更新は one-hot 位相を持つ11相状態機械として書き直せた。

## 15.4 読み出しと状態遷移を分離できる

特に時計について

$$
Q=e^{t/\tau_0}
$$

を状態とすれば、状態更新中に $t$ や $\log Q$ を必要としない。

$$
t=\tau_0\log|Q|
$$

は読み出し時だけに使う。同様に $x,y,r$ も読み出しであり、状態更新へ戻さない。

$$
\boxed{
\text{内部状態の演算形式}
\neq
\text{観測者が読む座標形式}
}
$$

を具体的に実装できた。

## 15.5 高精度化と状態閉包の修正を区別できる

角度級数の次数と数値積分次数を上げても、永続状態を追加する必要はなかった。一方、$\nu$ の内部化は精度向上とは無関係に、外部情報漏れを除去するため5→6へ状態数を増やした。

したがって少なくとも本モデルでは

$$
\boxed{
\text{高精度化}
\neq
\text{状態閉包の修正}
}
$$

である。

# 16. 限界

本稿には以下の限界がある。

1. **六永続状態の最小性は未証明である。**  
   別座標、保存量、あるいはより根本的な状態表現を用いれば、より少ない永続状態で閉じる可能性がある。

2. **11相 machine state の最小性・自然界での必然性は未証明である。**  
   11相は RK4 を隠れた途中値なしに展開した現在実装の構造であり、自然界が11段階で作用するという主張ではない。

3. **既存B方式全体の「完全乗法内部演算」適合は未証明である。**  
   本稿で追加確認したのは、外部 $\nu$ を第六状態 $\mathcal N$ へ内部化しても旧出力が厳密に保たれることまでである。

4. **Paper 4 の物理近似を引き継ぐ。**  
   保存部分は Schwarzschild/test-particle 型、散逸は leading 2.5PN であり [4,8-11]、完全な有限質量比一般相対論ではない。

5. **2.5PN より高い物理次数は検証していない。**  
   数値積分を8次へ上げても、物理モデル自身が2.5PNから高次化するわけではない。

6. **$U,H$ の複素表現は実装上の選択である。**  
   複素数を基本的実在と主張していない。

7. **全物理モデルの厳密解への収束は証明していない。**  
   角度級数は既知解析関数に対して明確な収束を示したが、モデル全体の無限次数極限が厳密GRと一致することは別問題である。

したがって本稿の到達点は、自然界の最小レジスタや万能相互作用の発見ではなく、**現在の二体重力モデルにおける状態閉包監査を2→5→6まで実行し、外部 $\nu$ を除去した**ところに限定される。

# 17. 再現性と研究記録

本稿ではプログラム監査そのものが実験であるため、最終コードだけでなく途中経路も保存した。

Google Drive の研究フォルダ

```text
匿名頂点状態生成幾何-匿名内部観測者から不変な関係量の体系/
第五思考実験_乗法状態内部展開と近似次数検証_20260924/
```

に、少なくとも以下を保存している。

- `00_reference_from_paper4/` — Paper 4 の基準コードと C1 生データ
- `01_五状態同期写像の再構成/` — 第一監査：五永続状態
- `02_乗法状態表現と同期インライン化の経緯/` — 途中案・却下案を含む
- `03_A局所exp-log_vs_B内部状態展開/` — A/B比較
- `04_B方式_11相同期状態機械_厳密監査/` — B方式の静的・数値監査
- `05_演算ブロック完全展開_Laurent_N12/` — Laurent 展開と10万点照合
- `06_近似次数と状態_S不増加検証/` — $N=12,16,20,24$ および4,6,8次
- `08_ν第6状態内部化_対照実験_20260925/` — 外部 $\nu$ の第六状態内部化、12条件・長時間対照、bitwise監査
- `paper5_revision_v2_figures_tables_20260926/` — 改訂版図表、元データ、図表再生成プログラム、再現性ZIP

改訂版図表は保存済み CSV / JSON から再生成できる。図表再生成プログラムは

```text
paper5_revision_v2_figures_tables_20260926/scripts/
generate_paper5_revision_v2.py
```

に保存した。$\nu$ 内部化の主要再現コードは

```text
08_ν第6状態内部化_対照実験_20260925/
run_nu_internal_state_control_v1.py
```

である。図2の軌道図は、同じフォルダの

```text
plot_c1_orbit_6state_audit_v1.py
```

で再生成できる。このプログラムは保存済み C1 CSV を参照比較用に読み込み、第六状態生成器の12000 macrostepを別途生成してから重ね描画する。生成中に参照軌道を状態更新へ帰還させない。

本稿の改訂で、荷電8状態実験（09以降）の結果は使用しない。第五思考実験は重力系の状態閉包監査だけで閉じる。

# 18. 結論

第四思考実験 [4] は、表面上は

$$
X_{n+1}=S_nX_n
$$

という二成分状態の相互作用から相対論的重力軌道を再構成する実験であった。

しかし本稿の第一監査により、実際の軌道生成には

$$
\chi,p,e,\phi,t
$$

という五つの永続状態が使われ、$X=(a,b)$ はそれらから作られる派生表現であることが明らかになった。そこで

$$
Z_5=(U,P,E,H,Q)
$$

へ再符号化し、RK4 の途中値と処理位相も明示 state として展開することで、Paper 4 の逐次処理を11相同期状態機械へ再構成した。

しかし、監査をそこで止めると不十分だった。放射反作用核には対称質量比 $\nu$ が外部引数として残っていたため、実際には

$$
Z_5'=F(Z_5;\nu)
$$

であった。

そこで

$$
\boxed{
Z_6=(U,P,E,H,Q,\mathcal N),
\qquad
\mathcal N=\nu
}
$$

とし、$\mathcal N$ を11相 one-hot 相互作用の恒等行として内部化した。新遷移関数は外部 `nu` 引数を持たない。

12条件 $\times$ 4000 macrostep では旧5状態と新6状態の最初の5成分、および $r,x,y,t$ が bitwise identical であり、$\mathcal N$ の bitwise mismatch も0であった。代表3条件 $\times$ 12000 macrostep でも比較量の最大差は0、$\mathcal N$ の microstep 恒等保持違反は0件だった。したがって今回の5→6修正は物理式や軌道を変える改変ではなく、**外部に残っていた物理情報の置き場所を状態内部へ移した状態閉包の修正**である。

また、角度級数の打切り次数や数値積分次数を高めても、新しい永続物理状態を追加する必要はなかった。したがって本稿では、

$$
\boxed{
\text{状態閉包の修正}
\neq
\text{数値高精度化}
}
$$

を具体的に分離できた。

本稿の到達点を最も限定して書けば、

$$
\boxed{
\begin{array}{c}
\text{Paper 4 の見かけの二状態}\
\downarrow\\
\text{第一監査：五永続状態}\
\downarrow\\
\text{追加監査：外部 }\nu\text{ を発見}\
\downarrow\\
\text{六永続状態 }Z_6\text{ へ内部化}
\end{array}
}
$$

という自己修正の経路を、数式・コード・bitwise 対照実験で追跡したことである。

これは六状態が自然界の最小レジスタであること、11相が基本相互作用の実在構造であること、あるいは一般相対論を六状態から導出したことを意味しない。現在確定したのは、既知の二体重力モデルを逆構成したこの限定系について、**次状態を決める情報を隠れた外部パラメータとして残さず、六つの永続論理状態と明示 machine/work state を持つ同期状態写像へ再配置できた**ということまでである。

# 参考文献

## 本研究系列

[0] 木原範昭, **状態と関係性から時空・粒子・相互作用の創発を探索する研究プロジェクト設計 ― 不確実性下の基礎物理モデル探索に対する予備設計と実装可能性検証の枠組み ―**, v1.0 (2026-09-20). Concept DOI: 10.5281/zenodo.22851944; Version DOI: 10.5281/zenodo.22851945.

[1] 木原範昭, **第一思考実験：等価原理を無名な二自由度まで遡る ― 面積読み出しの保存から、二乗和の符号・曲率半径・三つの系（慣性／一様加速／回転）を導く ―**, v1.0 (2026-09-20). Concept DOI: 10.5281/zenodo.22857952; Version DOI: 10.5281/zenodo.22857953.

[2] 木原範昭, **第二思考実験：クーロン型の逆二乗の力を、積の読み出しと面積の時計から導く ― 二つの値と線形の法則のまま、焦点・引力と斥力・中性を出し、二つの値では書けないものを確定する ―**, v1.2 (2026-09-21). Concept DOI: 10.5281/zenodo.22867336; Version DOI: 10.5281/zenodo.22876596.

[3] 木原範昭, **第三思考実験：a b 2体の天体の軌道から a b 2体の状態を読み出す ― 意味を定めない二つの値と $G=c=M=1$ の知識だけで、軌道と重力波の読み出しから質量比・相対距離・相対時間を得る ―**, v1.5 (2026-09-23). Concept DOI: 10.5281/zenodo.22909660; Version DOI: 10.5281/zenodo.22909661.

[4] 木原範昭, **第四思考実験：固定作用 $S$ から状態依存作用 $S_n$ へ ― 第三思考実験で外部計算していた相対論的重力軌道を、$X_{n+1}=S_nX_n$ の逐次相互作用から再生成できるか ―**, v1.0 (2026-09-23). Concept DOI: 10.5281/zenodo.22919787; Version DOI: 10.5281/zenodo.22919788.

## 外部文献

[5] R. Zwanzig, “Memory Effects in Irreversible Thermodynamics,” *Physical Review* **124**, 983–992 (1961). DOI: 10.1103/PhysRev.124.983.

[6] H. F. Trotter, “On the Product of Semi-Groups of Operators,” *Proceedings of the American Mathematical Society* **10**(4), 545–551 (1959). DOI: 10.1090/S0002-9939-1959-0108732-6.

[7] J. C. Butcher, “On Runge-Kutta Processes of High Order,” *Journal of the Australian Mathematical Society* **4**(2), 179–194 (1964). DOI: 10.1017/S1446788700023387.

[8] C. G. Darwin, “The Gravity Field of a Particle,” *Proceedings of the Royal Society of London A* **249**, 180–194 (1959). DOI: 10.1098/rspa.1959.0015.

[9] A. Pound and E. Poisson, “Osculating Orbits in Schwarzschild Spacetime, with an Application to Extreme Mass-Ratio Inspirals,” *Physical Review D* **77**, 044013 (2008). DOI: 10.1103/PhysRevD.77.044013.

[10] L. Blanchet, “Post-Newtonian Theory for Gravitational Waves,” *Living Reviews in Relativity* **27**, 4 (2024). DOI: 10.1007/s41114-024-00050-z.

[11] L. Blanchet, G. Faye, and D. Trestini, “Gravitational Radiation Reaction for Compact Binary Systems at the Fourth-and-a-Half Post-Newtonian Order,” *Classical and Quantum Gravity* **42**, 065015 (2025). DOI: 10.1088/1361-6382/adac9d.

---

# 付録A：論文本文で使用する図

改訂版の本文図は

```text
paper5_revision_v2_figures_tables_20260926/
figures_png_preview/
```

に保存している。

1. `paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure01.png` — 状態監査の経路 2→5→6
2. `08_ν第6状態内部化_対照実験_20260925/figures/figure_orbit_c1_paper4_vs_6state.png` — Paper 4 C1 と第六状態内部化後の再生成軌道
3. `paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure02.png` — 状態閉包・同期更新・読み出し非帰還の範囲
4. `paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure03.png` — 12条件 B方式 strict audit
5. `paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure04.png` — 角度級数 $N=12,16,20,24$ の収束
6. `paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure05.png` — 4・6・8次数値積分の収束
7. `paper5_revision_v2_figures_tables_20260926/figures_png_preview/figure06.png` — 9条件における高次化比較

同じ改訂パッケージには SVG 原図、表CSV、元データ、図表再生成プログラム、および再現性ZIPを含める。

---

# 付録B：主要な再現用ファイル

- `08_ν第6状態内部化_対照実験_20260925/run_nu_internal_state_control_v1.py`
- `08_ν第6状態内部化_対照実験_20260925/plot_c1_orbit_6state_audit_v1.py`
- `08_ν第6状態内部化_対照実験_20260925/figures/figure_orbit_c1_paper4_vs_6state.svg`
- `08_ν第6状態内部化_対照実験_20260925/figures/figure_orbit_c1_paper4_vs_6state_summary.txt`
- `08_ν第6状態内部化_対照実験_20260925/nu_internal_state_control_case_summary.csv`
- `08_ν第6状態内部化_対照実験_20260925/nu_internal_state_control_summary.json`
- `08_ν第6状態内部化_対照実験_20260925/long_run_checks.json`
- `08_ν第6状態内部化_対照実験_20260925/saved_C1_check.json`
- `paper5_revision_v2_figures_tables_20260926/EXPERIMENT_INDEX_v2.csv`
- `paper5_revision_v2_figures_tables_20260926/FIGURE_INDEX_v2_ja.md`
- `paper5_revision_v2_figures_tables_20260926/tables/table01_state_audit_summary_v2.csv`
- `paper5_revision_v2_figures_tables_20260926/tables/table02_approximation_and_state_scope_v2.csv`
- `paper5_revision_v2_figures_tables_20260926/tables/table03_key_numerical_results_v2.csv`
- `paper5_revision_v2_figures_tables_20260926/scripts/generate_paper5_revision_v2.py`
- `paper5_revision_v2_figures_tables_20260926/paper5_revision_v2_figures_tables_20260926.zip`

荷電8状態実験（09以降）は本稿の証拠として使用しない。
